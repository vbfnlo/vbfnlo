c************** Calling the Fa functions*************************
       subroutine HexBoxFa6(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,
     &   p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa6)
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
       Complex*16 Fa6(61:72)
       COMMON/HexBoxFaFunctions/Fa
       Fa6(61) = -(p1mup4*F(227))+p4mup4*F(228)+p5mup4*F(229)+p6mup4*F
     &   (1310)
       Fa6(62) = p1mup4*F(230)+p6mup4*F(231)+p4mup4*F(232)+p5mup4*F(23
     &   3)
       Fa6(63) = p6mup4*F(234)+p4mup4*F(235)+p5mup4*F(236)+p1mup4*F(23
     &   7)
       Fa6(64) = mup3mup4*F(241)+p1mup3*((p4mup4+p5mup4)*F(242)+p6mup4
     &   *F(244)+8*(p2mup4*F(243)-p1mup4*F(246)))+p4mup3*p4mup4*F(250
     &   )+p4mup4*p5mup3*F(250)+p4mup3*p5mup4*F(250)+p5mup3*p5mup4*F(
     &   250)+p4mup3*p6mup4*F(251)+p5mup3*p6mup4*F(251)+p6mup3*p6mup4
     &   *F(252)+p2mup3*p6mup4*F(255)+8*(p2mup3*p4mup4*F(248)+p2mup3*
     &   p5mup4*F(248)-p2mup4*(p4mup3+p5mup3)*F(249)-p2mup3*p2mup4*F(
     &   253)-p1mup4*(p2mup3*F(245)-(p4mup3+p5mup3)*F(247)+p6mup3*F(2
     &   54))+p2mup4*p6mup3*F(256))+(p4mup4+p5mup4)*p6mup3*F(257)
       Fa6(65) = -(p1mup4*F(260))+8*((p4mup4+p5mup4)*F(259)+p2mup4*F(2
     &   61)+p6mup4*F(262))
       Fa6(66) = p1mup3*F(263)+8*(p2mup3*F(264)+p6mup3*F(265)-(p4mup3+
     &   p5mup3)*F(266))
       Fa6(67) = 2*(p4mup4+p5mup4)*(-(p2mup3*F(268))+(p4mup3+p5mup3)*F
     &   (269)-p6mup3*F(270)+p1mup3*F(271))+mup3mup4*F(267)*P(116)
       Fa6(68) = mup3mup4*F(272)-2*(p1mup3*(p2mup4*F(273)+p6mup4*F(275
     &   )-p1mup4*F(276))-p2mup3*p6mup4*F(278)-p2mup4*(p4mup3+p5mup3)
     &   *F(279)+p4mup3*p6mup4*F(280)+p5mup3*p6mup4*F(280)+p6mup3*p6m
     &   up4*F(281)+p2mup3*p2mup4*F(282)+p1mup4*(-(p2mup3*F(274))+(p4
     &   mup3+p5mup3)*F(277)-p6mup3*F(283))+p2mup4*p6mup3*F(284))
       Fa6(69) = 4*mup3mup4*F(261)+8*(p4mup4+p5mup4)*(p1mup3*F(273)-(p
     &   4mup3+p5mup3)*F(279)+p2mup3*F(282)+p6mup3*F(284))
       Fa6(70) = 4*mup3mup4*F(287)-8*(-(p4mup3*p6mup4*F(280))-p5mup3*p
     &   6mup4*F(280)-p6mup3*p6mup4*F(281)-p1mup3*p6mup4*F(290)+p1mup
     &   4*p2mup3*F(291)+p1mup3*p1mup4*F(292)-p1mup4*p4mup3*F(293)-p1
     &   mup4*p5mup3*F(293)-p1mup4*p6mup3*F(294)-p2mup3*p6mup4*F(295)
     &   )
       Fa6(71) = p6mup4*F(296)+8*(p1mup4*p4sq*F(267)+(p4mup4+p5mup4)*F
     &   (297))
       Fa6(72) = p1mup3*F(298)+p2mup3*F(299)-(p4mup3+p5mup3)*F(300)+p6
     &   mup3*F(301)
       Return
       End
