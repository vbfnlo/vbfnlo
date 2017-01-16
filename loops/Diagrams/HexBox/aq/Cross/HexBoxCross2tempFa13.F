c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa13(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   13)
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
       Complex*16 Fa13(229:247)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa13(229) = (p1mup6+p6mup6)*((p1mup3+p6mup3)*F(536)+p2mup3*F(54
     &   0)+p5mup3*F(683)+p4mup3*F(685))
       Fa13(230) = (p1mup6+p6mup6)*((p1mup3+p6mup3)*F(592)+p4mup3*F(59
     &   3)+p5mup3*F(594)-p2mup3*F(595))
       Fa13(231) = (p1mup6+p6mup6)*(p5mup3*F(681)+(p1mup3+p6mup3)*F(68
     &   2)+p4mup3*F(684)+p2mup3*F(686))
       Fa13(232) = (p1mup6+p6mup6)*(p5mup3*F(567)+p4mup3*F(570)+(p1mup
     &   3+p6mup3)*F(577)-p2mup3*F(579))
       Fa13(233) = (p1mup6+p6mup6)*(-(p2mup3*F(729))+p4mup3*F(730)+p5m
     &   up3*F(731)+(p1mup3+p6mup3)*F(1849))
       Fa13(234) = p5mup6*(p2mup3*F(736)+(p5mup3+p6mup3)*F(737)+p4mup3
     &   *F(738)+p1mup3*F(739))
       Fa13(235) = p5mup6*(p5mup3+p6mup3)*F(741)+p4mup3*p5mup6*F(742)
       Fa13(236) = p5mup6*(p2mup3*F(743)-p1mup3*F(744))
       Fa13(237) = p5mup6*(p1mup3*F(745)-p4mup3*F(746)+p2mup3*F(747)+(
     &   p5mup3+p6mup3)*F(748))
       Fa13(238) = p5mup6*(p1mup3*F(750)+p2mup3*F(751)+(p5mup3+p6mup3)
     &   *F(752)+p4mup3*F(753))
       Fa13(239) = (p6mup3*F(755))/s56-8*(p2mup3*F(756)-p4mup3*F(757)-
     &   p5mup3*F(758))
       Fa13(240) = (p6mup3*F(759))/s56+8*(p2mup3*F(760)-p4mup3*F(761)-
     &   p5mup3*F(762)+p1mup3*F(763))
       Fa13(241) = -4*mup3mup6*s56*F(758)+8*p5mup6*(-(p2mup3*F(767))+p
     &   4mup3*F(768)-p1mup3*F(769)+(p5mup3+p6mup3)*F(770))
       Fa13(242) = p5mup3*F(761)-p1mup3*F(772)+p4mup3*F(773)+(p6mup3*F
     &   (774))/s56-p2mup3*F(775)
       Fa13(243) = 4*p6mup6*F(758)+(p5mup6*F(776))/s56
       Fa13(244) = (p6mup3*F(781))/s56+4*(p2mup3*F(747)+p5mup3*F(748)+
     &   p4mup3*F(775)+p1mup3*F(782))
       Fa13(245) = p1mup3*F(783)+p2mup3*F(784)+p4mup3*F(785)+(p5mup3+p
     &   6mup3)*F(786)
       Fa13(246) = p5mup3*F(763)+p2mup3*F(782)+(p6mup3*F(795))/s56+p1m
     &   up3*F(796)+p4mup3*F(797)
       Fa13(247) = p2mup3*F(748)-p5mup3*F(762)+p1mup3*F(798)+p4mup3*F(
     &   799)-(p6mup3*F(800))/s56
       Return
       End
