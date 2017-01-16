c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa4(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p
     &   3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa4
     &   )
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBoxCross2FhlFunctions/F
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
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
      Complex*16 Fa(392),F(1853)
      Real*8 P(1355) 
        COMMON/HexBoxCross2PFunctions/P
       Complex*16 Fa4(58:76)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa4(58) = p1mup6*F(166)+p5mup6*F(167)+p6mup6*F(168)+p4mup6*F(16
     &   9)
       Fa4(59) = p6mup6*F(185)+p5mup6*F(186)+p1mup6*F(187)+p4mup6*F(18
     &   8)
       Fa4(60) = p5mup6*F(189)+p6mup6*F(190)+p1mup6*F(191)+p4mup6*F(19
     &   2)
       Fa4(61) = p6mup6*F(194)-p5mup6*F(195)+p1mup6*F(196)+p4mup6*F(19
     &   7)
       Fa4(62) = p5mup6*F(198)+p6mup6*F(199)-p1mup6*F(200)+p4mup6*F(20
     &   1)
       Fa4(63) = p6mup6*F(202)+p5mup6*F(203)+p1mup6*F(204)+p4mup6*F(20
     &   5)
       Fa4(64) = p5mup6*F(206)+p6mup6*F(207)+p1mup6*F(208)+p4mup6*F(20
     &   9)
       Fa4(65) = p1mup6*F(179)-p5mup6*F(210)-p6mup6*F(211)+p4mup6*F(21
     &   2)
       Fa4(66) = p1mup6*F(179)+p4mup6*F(213)-p5mup6*F(214)-p6mup6*F(21
     &   5)
       Fa4(67) = p4mup6*F(216)-p5mup6*F(217)
       Fa4(68) = -(p6mup6*F(218))+p4mup6*F(219)+p5mup6*F(220)+p1mup6*F
     &   (221)
       Fa4(69) = p4mup6*F(225)+p5mup6*F(226)+p6mup6*F(227)+p1mup6*F(22
     &   8)
       Fa4(70) = p1mup6*F(229)+p5mup6*F(230)+p6mup6*F(231)+p4mup6*F(23
     &   2)
       Fa4(71) = p1mup6*F(233)+p5mup6*F(234)+p6mup6*F(235)-2*p4mup6*F(
     &   236)
       Fa4(72) = p6mup6*F(240)+p4mup6*F(241)+p1mup6*F(242)+p5mup6*F(24
     &   3)
       Fa4(73) = p4mup6*F(244)+p6mup6*F(245)-p1mup6*F(246)+p5mup6*F(24
     &   7)
       Fa4(74) = p5mup6*F(248)+p6mup6*F(249)+p1mup6*F(250)+p4mup6*F(25
     &   1)
       Fa4(75) = p5mup6*F(252)+p6mup6*F(253)+p1mup6*F(254)+p4mup6*F(25
     &   5)
       Fa4(76) = -(p5mup6*F(220))+p6mup6*F(256)+p1mup6*F(257)+p4mup6*F
     &   (258)
       Return
       End
