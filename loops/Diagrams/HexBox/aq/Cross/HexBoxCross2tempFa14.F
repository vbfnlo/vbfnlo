c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa14(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   14)
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
       Complex*16 Fa14(248:266)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa14(248) = p5mup6*(p1mup3*F(802)+p2mup3*F(803)-p4mup3*F(804)+(
     &   p5mup3+p6mup3)*F(805))
       Fa14(249) = (p6mup3*F(806))/s56-2*(p5mup3*F(763)+p2mup3*F(782)+
     &   p4mup3*F(797)-p1mup3*F(807))
       Fa14(250) = (p5mup3+p6mup3)*F(741)+p4mup3*F(742)
       Fa14(251) = p1mup3*F(745)-p4mup3*F(746)+p2mup3*F(747)+(p5mup3+p
     &   6mup3)*F(748)
       Fa14(252) = p2mup3*F(736)+(p5mup3+p6mup3)*F(737)+p4mup3*F(738)+
     &   p1mup3*F(739)
       Fa14(253) = p2mup3*F(743)-p1mup3*F(744)
       Fa14(254) = -(p2mup3*F(767))+p4mup3*F(768)+(p5mup3+p6mup3)*F(77
     &   0)+p1mup3*F(818)
       Fa14(255) = p5mup6*(p2mup3*F(820)+p1mup3*F(821)+(p5mup3+p6mup3)
     &   *F(822)+p4mup3*F(823))
       Fa14(256) = p5mup6*(p1mup3*F(824)-p4mup3*F(825)+(p5mup3+p6mup3)
     &   *F(826)+p2mup3*F(827))
       Fa14(257) = p5mup6*((p5mup3+p6mup3)*F(797)-p2mup3*F(829)-p1mup3
     &   *F(830)+p4mup3*F(831))
       Fa14(258) = p5mup6*(p2mup3*F(775)+(p5mup3+p6mup3)*F(799)+p1mup3
     &   *F(832)-p4mup3*F(833))
       Fa14(259) = -(p1mup3*F(802))-p2mup3*F(803)+p4mup3*F(804)+(p5mup
     &   3+p6mup3)*F(835)
       Fa14(260) = p5mup6*(-(p1mup3*F(836))-p2mup3*F(837)+p4mup3*F(838
     &   )+(p5mup3+p6mup3)*F(1850))
       Fa14(261) = p2mup3*F(820)+p1mup3*F(821)+(p5mup3+p6mup3)*F(822)+
     &   p4mup3*F(823)
       Fa14(262) = p1mup3*F(824)-p4mup3*F(825)+(p5mup3+p6mup3)*F(826)+
     &   p2mup3*F(827)
       Fa14(263) = -(p1mup3*F(836))-p2mup3*F(837)+p4mup3*F(838)+(p5mup
     &   3+p6mup3)*F(1850)
       Fa14(264) = (p5mup3+p6mup3)*F(797)-p2mup3*F(829)-p1mup3*F(830)+
     &   p4mup3*F(831)
       Fa14(265) = p2mup3*F(775)+(p5mup3+p6mup3)*F(799)+p1mup3*F(832)-
     &   p4mup3*F(833)
       Fa14(266) = p1mup6*F(840)+p6mup6*F(841)+p2mup6*F(842)+p4mup6*F(
     &   843)+p5mup6*F(844)
       Return
       End
