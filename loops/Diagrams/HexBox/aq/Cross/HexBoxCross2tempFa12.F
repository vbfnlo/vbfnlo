c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa12(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   12)
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
       Complex*16 Fa12(210:228)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa12(210) = p2mup6*F(717)-p4mup6*F(718)
       Fa12(211) = p5mup6*F(562)+(p1mup6+p6mup6)*F(716)
       Fa12(212) = p2mup6*F(554)-p5mup6*F(680)
       Fa12(213) = (p1mup6+p6mup6)*F(680)+p2mup6*F(689)
       Fa12(214) = p5mup3*F(681)+(p1mup3+p6mup3)*F(682)+p4mup3*F(684)+
     &   p2mup3*F(686)
       Fa12(215) = p5mup3*F(567)+p4mup3*F(570)+(p1mup3+p6mup3)*F(577)-
     &   p2mup3*F(579)
       Fa12(216) = (p1mup3+p6mup3)*F(722)+p4mup3*F(723)+p5mup3*F(724)-
     &   p2mup3*F(725)
       Fa12(217) = 2*(p5mup3*F(567)+p4mup3*F(570))+p2mup3*F(727)+(p1mu
     &   p3+p6mup3)*F(728)
       Fa12(218) = (p1mup3+p4mup3+p6mup3)*F(530)+p5mup3*F(531)-p2mup3*
     &   F(532)
       Fa12(219) = -2*p5mup6*F(555)+p4mup6*F(716)
       Fa12(220) = 8*(p1mup6+p6mup6)*F(555)+4*p4mup6*F(562)
       Fa12(221) = -4*p2mup6*F(562)+8*(p1mup6+p6mup6)*F(680)
       Fa12(222) = -(p5mup6*F(555))+p4mup6*F(716)
       Fa12(223) = p5mup6*F(680)+p2mup6*F(716)
       Fa12(224) = -(p2mup3*F(729))+p4mup3*F(730)+p5mup3*F(731)+(p1mup
     &   3+p6mup3)*F(1849)
       Fa12(225) = p1mup3*(p2mup6*F(536)+p4mup6*F(682))+p2mup6*(p6mup3
     &   *F(536)+p2mup3*F(540)+p5mup3*F(683)+p4mup3*F(685))+p4mup6*(p
     &   5mup3*F(681)+p6mup3*F(682)+p4mup3*F(684)+p2mup3*F(686))-mup3
     &   mup6*F(709)
       Fa12(226) = (p1mup3+p6mup3)*F(533)+p4mup3*F(534)+p5mup3*F(535)+
     &   p2mup3*F(644)
       Fa12(227) = (p1mup3+p6mup3)*F(533)+p4mup3*F(534)+p5mup3*F(535)-
     &   p2mup3*F(732)
       Fa12(228) = p4mup6*(p5mup3*F(567)+p4mup3*F(570)+(p1mup3+p6mup3)
     &   *F(577)-p2mup3*F(579))+p2mup6*((p1mup3+p6mup3)*F(592)+p4mup3
     &   *F(593)+p5mup3*F(594)-p2mup3*F(595))-mup3mup6*F(710)
       Return
       End
