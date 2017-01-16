c************** Calling the Fa functions*************************
       subroutine HexBoxFa12(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa12)
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
       Complex*16 Fa12(133:144)
       COMMON/HexBoxFaFunctions/Fa
       Fa12(133) = p2mup3*F(311)+p1mup3*F(312)+(p4mup3+p5mup3)*F(315)+
     &   p6mup3*F(316)
       Fa12(134) = (p4mup4+p5mup4)*F(335)-2*p2mup4*F(445)
       Fa12(135) = -(p2mup3*F(309))+(p4mup3+p5mup3)*F(329)+p1mup3*F(33
     &   0)-p6mup3*F(332)
       Fa12(136) = -(p2mup3*F(334))+p1mup3*F(335)
       Fa12(137) = (p4mup4+p5mup4)*F(443)-p2mup4*F(445)
       Fa12(138) = (p4mup3+p5mup3)*F(329)-p6mup3*F(332)+p2mup3*F(377)-
     &   p1mup3*F(449)
       Fa12(139) = p2mup3*F(273)+p1mup3*F(327)-(p4mup3+p5mup3)*F(328)+
     &   p6mup3*F(331)
       Fa12(140) = (p4mup4+p5mup4)*F(288)+p6mup4*F(445)
       Fa12(141) = (p4mup4+p5mup4)*F(334)-p1mup4*F(445)
       Fa12(142) = p6mup4*F(450)+p1mup4*F(453)
       Fa12(143) = -8*(p1mup3*F(337)-(p4mup3+p5mup3)*F(338)+p6mup3*F(3
     &   39))+p2mup3*F(454)
       Fa12(144) = -2*(p1mup3*F(337)-(p4mup3+p5mup3)*F(338)+p6mup3*F(3
     &   39))+p2mup3*F(455)
       Return
       End
