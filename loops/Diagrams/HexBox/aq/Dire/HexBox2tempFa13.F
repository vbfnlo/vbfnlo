c************** Calling the Fa functions*************************
       subroutine HexBox2Fa13(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa13)
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
       Complex*16 Fa13(229:247)
       COMMON/HexBox2FaFunctions/Fa
       Fa13(229) = 8*(p5mup6*p6mup3*F(650)+p6mup3*p6mup6*F(650)+p2mup3
     &   *p5mup6*F(667)+p2mup3*p6mup6*F(667)+p1mup3*p5mup6*F(726)+p1m
     &   up3*p6mup6*F(726)+p1mup3*p2mup6*F(743)+p2mup3*p2mup6*F(744)-
     &   p1mup6*p2mup3*F(745)-p1mup3*p1mup6*F(747)+p2mup6*p6mup3*F(74
     &   8)+p1mup6*p6mup3*F(749)+p5mup3*((p5mup6+p6mup6)*F(650)+p2mup
     &   6*F(748)+p1mup6*F(749))+p4mup3*((p5mup6+p6mup6)*F(649)-p2mup
     &   6*F(746)+p1mup6*F(750)))+mup3mup6*F(791)
       Fa13(230) = p2mup3*F(792)-p4mup3*F(793)+p1mup3*F(794)+(p5mup3+p
     &   6mup3)*F(795)
       Fa13(231) = p2mup3*F(796)+p1mup3*F(797)+(p5mup3+p6mup3)*F(798)+
     &   p4mup3*F(799)
       Fa13(232) = -(p2mup3*F(745))-p1mup3*F(747)+(p5mup3+p6mup3)*F(74
     &   9)+p4mup3*F(750)
       Fa13(233) = -((p5mup6+p6mup6)*(p2mup3*F(745)+p1mup3*F(747)-(p5m
     &   up3+p6mup3)*F(749)-p4mup3*F(750)))
       Fa13(234) = p4mup6*F(801)-(p5mup6+p6mup6)*F(803)+p2mup6*F(804)-
     &   p1mup6*F(805)
       Fa13(235) = p2mup6*F(806)+p4mup6*F(807)+12*(p1mup6+p5mup6+p6mup
     &   6)*F(1635)
       Fa13(236) = -(p2mup3*p4mup6*F(745))-p1mup3*p4mup6*F(747)+p4mup6
     &   *p5mup3*F(749)+p4mup6*p6mup3*F(749)+p4mup3*p4mup6*F(750)+p2m
     &   up3*p2mup6*F(760)+p1mup3*p2mup6*F(761)+p2mup6*p4mup3*F(763)+
     &   p2mup6*(p5mup3+p6mup3)*F(764)+mup3mup6*F(805)
       Fa13(237) = p2mup6*F(670)-p4mup6*F(676)
       Fa13(238) = 2*(p5mup6+p6mup6)*F(676)+p2mup6*F(790)
       Fa13(239) = p2mup6*F(671)-p1mup6*F(676)
       Fa13(240) = p1mup3*F(811)+p2mup3*F(812)+p4mup3*F(813)+(p5mup3+p
     &   6mup3)*F(814)
       Fa13(241) = 16*(p4mup3*p4mup6*F(746)-p4mup6*p5mup3*F(748)-p4mup
     &   6*p6mup3*F(748)-p1mup3*((p5mup6+p6mup6)*F(643)+p4mup6*F(743)
     &   -2*p1mup6*F(761))+2*(p1mup6*p4mup3*F(763)+p1mup6*p5mup3*F(76
     &   4)+p1mup6*p6mup3*F(764))+p4mup3*p5mup6*F(766)+p4mup3*p6mup6*
     &   F(766)+p2mup3*(-(p4mup6*F(744))+2*p1mup6*F(760)+(p5mup6+p6mu
     &   p6)*F(767))-(p5mup3+p6mup3)*(p5mup6+p6mup6)*F(768))+mup3mup6
     &   *F(816)
       Fa13(242) = (p5mup6+p6mup6)*F(671)+p1mup6*F(790)
       Fa13(243) = p6mup3*F(820)+p1mup3*F(821)+p2mup3*F(822)+p4mup3*F(
     &   823)+p5mup3*F(824)
       Fa13(244) = (p5mup6+p6mup6)*(p1mup3*F(743)+p2mup3*F(744)-p4mup3
     &   *F(746)+(p5mup3+p6mup3)*F(748))
       Fa13(245) = -8*(p2mup3*p4mup6*F(744)-p1mup6*p2mup3*F(760)+p1mup
     &   3*((p5mup6+p6mup6)*F(725)+p4mup6*F(743)-p1mup6*F(761))-p4mup
     &   3*((p5mup6+p6mup6)*F(727)+p4mup6*F(746)+p1mup6*F(763))+(p5mu
     &   p3+p6mup3)*((p5mup6+p6mup6)*F(730)+p4mup6*F(748)-p1mup6*F(76
     &   4)))+mup3mup6*F(825)+p2mup3*(p5mup6+p6mup6)*F(826)
       Fa13(246) = p2mup3*F(760)+p1mup3*F(761)+p4mup3*F(763)+(p5mup3+p
     &   6mup3)*F(764)
       Fa13(247) = (p5mup6+p6mup6)*(p2mup3*F(760)+p1mup3*F(761)+p4mup3
     &   *F(763)+(p5mup3+p6mup3)*F(764))
       Return
       End
