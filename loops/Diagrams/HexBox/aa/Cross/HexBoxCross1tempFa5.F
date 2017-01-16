c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa5(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p
     &   3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa5
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
       Complex*16 Fa5(53:65)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa5(53) = p4mup4*F(178)+p1mup4*F(179)-p5mup4*F(180)+p6mup4*F(18
     &   1)
       Fa5(54) = p4mup4*F(197)+p5mup4*F(198)+p6mup4*F(199)+p1mup4*F(20
     &   0)
       Fa5(55) = -(p1mup4*F(201))+p4mup4*F(202)+2*p5mup4*F(203)+p6mup4
     &   *F(204)
       Fa5(56) = p4mup4*F(205)+p5mup4*F(206)+p6mup4*F(207)-p1mup4*F(20
     &   8)
       Fa5(57) = p4mup4*F(213)+p5mup4*F(214)-p1mup4*F(215)-p6mup4*F(21
     &   6)
       Fa5(58) = p4mup4*F(220)+p5mup4*F(221)+p6mup4*F(222)-p1mup4*F(22
     &   3)
       Fa5(59) = -(p1mup4*F(224))+p4mup4*F(225)+p5mup4*F(226)+p6mup4*F
     &   (1457)
       Fa5(60) = p4mup4*F(227)+p1mup4*F(228)+p5mup4*F(229)+p6mup4*F(23
     &   0)
       Fa5(61) = p4mup4*F(231)+p5mup4*F(232)+p6mup4*F(233)+p1mup4*F(23
     &   4)
       Fa5(62) = -(p6mup4*F(159))-p1mup4*F(235)+p5mup4*F(236)+p4mup4*F
     &   (237)
       Fa5(63) = p4mup4*F(239)+p5mup4*F(240)+p6mup4*F(241)+p1mup4*F(24
     &   2)
       Fa5(64) = p4mup4*F(243)+p5mup4*F(244)+p6mup4*F(245)+p1mup4*F(24
     &   6)
       Fa5(65) = p1mup4*F(247)+p5mup4*F(248)+p6mup4*F(249)+p4mup4*F(25
     &   0)
       Return
       End
