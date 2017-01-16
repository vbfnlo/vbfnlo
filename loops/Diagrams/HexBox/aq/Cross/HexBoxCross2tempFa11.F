c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa11(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   11)
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
       Complex*16 Fa11(191:209)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa11(191) = (p1mup3+p6mup3)*F(690)+p4mup3*F(691)+p5mup3*F(692)+
     &   p2mup3*F(693)
       Fa11(192) = (p1mup3+p6mup3)*F(694)+p4mup3*F(695)+p5mup3*F(696)+
     &   p2mup3*F(697)
       Fa11(193) = (p1mup3+p6mup3)*F(536)+p2mup3*F(540)+p5mup3*F(683)+
     &   p4mup3*F(685)
       Fa11(194) = (p1mup3+p6mup3)*F(592)+p4mup3*F(593)+p5mup3*F(594)-
     &   p2mup3*F(595)
       Fa11(195) = (p1mup6+p6mup6)*((p1mup3+p6mup3)*F(574)+p5mup3*F(57
     &   5)+p4mup3*F(576)-p2mup3*F(578))
       Fa11(196) = 8*(p1mup6*p6mup3*F(536)+p6mup3*p6mup6*F(536)+p1mup6
     &   *p2mup3*F(540)+p2mup3*p6mup6*F(540)+p4mup6*p6mup3*F(574)+p4m
     &   up6*p5mup3*F(575)+p4mup3*p4mup6*F(576)-p2mup3*p4mup6*F(578)+
     &   p5mup6*p6mup3*F(592)+p1mup3*((p1mup6+p6mup6)*F(536)+p4mup6*F
     &   (574)+p5mup6*F(592))+p4mup3*p5mup6*F(593)+p5mup3*p5mup6*F(59
     &   4)-p2mup3*p5mup6*F(595)+p1mup6*p5mup3*F(683)+p5mup3*p6mup6*F
     &   (683)+p4mup3*(p1mup6+p6mup6)*F(685))-4*mup3mup6*F(698)
       Fa11(197) = 8*(-(p2mup6*(p1mup3+p6mup3)*F(574))+p1mup3*p5mup6*F
     &   (577)+p5mup6*p6mup3*F(577)+p2mup3*p2mup6*F(578)-p2mup3*p5mup
     &   6*F(579)+p5mup3*(p5mup6*F(567)-p2mup6*F(575)+(p1mup6+p6mup6)
     &   *F(681))+p1mup6*p6mup3*F(682)+p6mup3*p6mup6*F(682)+p1mup3*(p
     &   1mup6+p6mup6)*F(682)+p4mup3*(p5mup6*F(570)-p2mup6*F(576)+(p1
     &   mup6+p6mup6)*F(684))+p2mup3*(p1mup6+p6mup6)*F(686))-mup3mup6
     &   *F(701)
       Fa11(198) = -((p1mup3+p6mup3)*F(703))+p4mup3*F(704)+p5mup3*F(70
     &   5)+p2mup3*F(706)
       Fa11(199) = -(p2mup6*F(688))-p4mup6*F(708)+(p1mup6+p6mup6)*F(70
     &   9)+p5mup6*F(710)
       Fa11(200) = -2*p2mup6*F(699)+p4mup6*F(712)+12*(p1mup6+p5mup6+p6
     &   mup6)*F(1847)
       Fa11(201) = p5mup3*F(537)+p4mup3*F(538)+2*((p1mup3+p6mup3)*F(53
     &   6)+p2mup3*F(540))
       Fa11(202) = (p1mup3+p6mup3)*F(568)+p5mup3*F(569)+p4mup3*F(571)+
     &   p2mup3*F(572)
       Fa11(203) = (p1mup3+p6mup3)*F(574)+p5mup3*F(575)+p4mup3*F(576)-
     &   p2mup3*F(578)
       Fa11(204) = (p1mup3+p6mup3)*F(536)+p2mup3*F(540)+p5mup3*F(581)+
     &   p4mup3*F(582)
       Fa11(205) = 16*p5mup6*F(555)-8*p4mup6*F(716)
       Fa11(206) = p4mup6*F(562)+(p1mup6+p6mup6)*F(717)
       Fa11(207) = 2*p5mup6*F(680)+p2mup6*F(716)
       Fa11(208) = p2mup6*F(562)+(p1mup6+p6mup6)*F(719)
       Fa11(209) = p4mup6*F(719)-p2mup6*F(720)
       Return
       End
