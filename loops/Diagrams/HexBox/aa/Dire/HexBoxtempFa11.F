c************** Calling the Fa functions*************************
       subroutine HexBoxFa11(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa11)
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
       Complex*16 Fa11(121:132)
       COMMON/HexBoxFaFunctions/Fa
       Fa11(121) = p1mup3*F(275)-p2mup3*F(278)+(p4mup3+p5mup3)*F(280)+
     &   p6mup3*F(281)
       Fa11(122) = (p4mup3+p5mup3)*F(280)+p6mup3*F(281)+p1mup3*F(290)+
     &   p2mup3*F(295)
       Fa11(123) = -(p1mup3*F(273))+(p4mup3+p5mup3)*F(279)-p2mup3*F(28
     &   2)+p6mup3*F(444)
       Fa11(124) = p6mup3*F(281)+p1mup3*F(290)+p2mup3*F(295)-(p4mup3+p
     &   5mup3)*F(379)
       Fa11(125) = (p4mup4+p5mup4)*F(343)-p1mup4*F(445)
       Fa11(126) = -4*(p1mup3*F(337)-(p4mup3+p5mup3)*F(338)+p6mup3*F(3
     &   39))+p2mup3*F(407)
       Fa11(127) = (p4mup4+p5mup4)*F(267)+p6mup4*F(446)
       Fa11(128) = (p4mup4+p5mup4)*F(343)-p1mup4*F(446)
       Fa11(129) = -(p1mup3*F(309))+(p4mup3+p5mup3)*F(310)+p2mup3*F(31
     &   3)-p6mup3*F(314)
       Fa11(130) = p2mup4*F(267)+p6mup4*F(335)
       Fa11(131) = p1mup4*F(335)-p2mup4*F(343)
       Fa11(132) = (p4mup4+p5mup4)*F(267)+2*p6mup4*F(445)
       Return
       End
