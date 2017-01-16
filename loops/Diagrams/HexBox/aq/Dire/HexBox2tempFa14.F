c************** Calling the Fa functions*************************
       subroutine HexBox2Fa14(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa14)
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
       Complex*16 Fa14(248:266)
       COMMON/HexBox2FaFunctions/Fa
       Fa14(248) = 8*(p5mup6+p6mup6)*F(670)+4*p4mup6*F(790)
       Fa14(249) = -(p1mup6*F(670))+p4mup6*F(671)
       Fa14(250) = p4mup6*F(834)+(p5mup6+p6mup6)*F(835)+p2mup6*F(836)
       Fa14(251) = 8*(p5mup6+p6mup6)*F(837)+p2mup6*F(838)+p1mup6*F(839
     &   )+p4mup6*F(840)
       Fa14(252) = p2mup6*F(841)+p4mup6*F(842)+8*(p5mup6+p6mup6)*F(843
     &   )+p1mup6*F(844)
       Fa14(253) = p1mup3*F(645)-p2mup3*F(646)+(p5mup3+p6mup3)*F(650)-
     &   p4mup3*F(728)
       Fa14(254) = p1mup3*F(643)-p4mup3*F(766)-p2mup3*F(767)+(p5mup3+p
     &   6mup3)*F(768)
       Fa14(255) = (p5mup6+p6mup6)*F(681)+p4mup6*F(790)
       Fa14(256) = p1mup3*F(854)+(p5mup3+p6mup3)*F(855)-2*p4mup3*F(856
     &   )+p2mup3*F(857)
       Fa14(257) = -8*(p1mup3*F(725)-p4mup3*F(727)+(p5mup3+p6mup3)*F(7
     &   30))+p2mup3*F(826)
       Fa14(258) = -(p1mup3*F(743))-p2mup3*F(744)+p4mup3*F(746)+(p5mup
     &   3+p6mup3)*F(861)
       Fa14(259) = (p5mup6+p6mup6)*F(670)+p4mup6*F(862)
       Fa14(260) = p4mup6*F(639)-p1mup6*F(670)
       Fa14(261) = p1mup3*F(743)+p2mup3*F(744)-p4mup3*F(746)+(p5mup3+p
     &   6mup3)*F(748)
       Fa14(262) = (p5mup6+p6mup6)*F(676)+p2mup6*F(862)
       Fa14(263) = p2mup6*F(639)-p1mup6*F(676)
       Fa14(264) = p1mup3*F(725)-p4mup3*F(727)+(p5mup3+p6mup3)*F(730)-
     &   p2mup3*F(866)
       Fa14(265) = (p5mup6+p6mup6)*(p1mup3*F(867)+p2mup3*F(868)-p4mup3
     &   *F(869)-(p5mup3+p6mup3)*F(1641))
       Fa14(266) = (p5mup6+p6mup6)*(-4*(p1mup3*F(725)-p4mup3*F(727)+(p
     &   5mup3+p6mup3)*F(730))+p2mup3*F(870))
       Return
       End
