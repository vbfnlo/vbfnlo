c************** Calling the Fa functions*************************
       subroutine HexBox2Fa4(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa4)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBox2FhlFunctions/F
       COMMON/HexBox2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,
     &   s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(385),F(1646)
      Real*8 P(1185) 
        COMMON/HexBox2PFunctions/P
       Complex*16 Fa4(58:76)
       COMMON/HexBox2FaFunctions/Fa
       Fa4(58) = p4mup6*F(148)+p5mup6*F(149)+p6mup6*F(150)+p1mup6*F(15
     &   1)
       Fa4(59) = 8*(p4mup6*F(158)+p5mup6*F(159)-p6mup6*F(160))+p1mup6*
     &   F(161)
       Fa4(60) = p5mup6*F(162)+p4mup6*F(163)+p6mup6*F(164)+p1mup6*F(16
     &   5)
       Fa4(61) = -(p1mup6*F(168))+2*(p4mup6*F(166)-p5mup6*F(167)+p6mup
     &   6*F(169))
       Fa4(62) = p4mup6*F(170)-p5mup6*F(171)+p6mup6*F(172)-p1mup6*F(17
     &   3)
       Fa4(63) = p4mup6*F(174)+p5mup6*F(175)+p6mup6*F(176)-p1mup6*F(17
     &   7)
       Fa4(64) = p1mup6*F(178)+p4mup6*F(179)-p6mup6*F(180)+p5mup6*F(18
     &   1)
       Fa4(65) = p4mup6*F(187)+p5mup6*F(188)+p6mup6*F(189)+p1mup6*F(19
     &   0)
       Fa4(66) = p4mup6*F(191)+p5mup6*F(192)+p6mup6*F(193)+p1mup6*F(19
     &   4)
       Fa4(67) = p4mup6*F(197)+p5mup6*F(198)+p6mup6*F(199)+p1mup6*F(20
     &   0)
       Fa4(68) = p4mup6*F(215)+p5mup6*F(216)+p6mup6*F(217)+p1mup6*F(21
     &   8)
       Fa4(69) = p5mup6*F(222)-p4mup6*F(223)+p6mup6*F(224)-p1mup6*F(22
     &   5)
       Fa4(70) = p5mup6*F(226)+p4mup6*F(227)+p6mup6*F(228)+p1mup6*F(22
     &   9)
       Fa4(71) = p5mup6*F(231)+p4mup6*F(232)+p6mup6*F(233)+p1mup6*F(23
     &   4)
       Fa4(72) = p5mup6*F(235)+p4mup6*F(236)+p6mup6*F(237)+p1mup6*F(23
     &   8)
       Fa4(73) = p4mup6*F(239)+p5mup6*F(240)+p6mup6*F(241)+p1mup6*F(24
     &   2)
       Fa4(74) = p4mup6*F(243)+p5mup6*F(244)+p6mup6*F(245)+p1mup6*F(24
     &   6)
       Fa4(75) = p4mup6*F(248)+p5mup6*F(249)+p6mup6*F(250)+p1mup6*F(25
     &   1)
       Fa4(76) = p4mup6*F(252)+p5mup6*F(253)+p6mup6*F(254)+p1mup6*F(25
     &   5)
       Return
       End
