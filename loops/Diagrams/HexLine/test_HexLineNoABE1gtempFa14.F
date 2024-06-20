c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g14(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa14)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/test_FhlFunctionsNoAbe1g/F
       COMMON/Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45
     -   ,s56,s16,s123,s234,s345
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
      Complex*16 Fa(587),F(6395)
      Real*8 P(7129) 
        COMMON/test_PFunctionsNoAbe1g/P
       Complex*16 Fa14(183:196)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa14(183) = p1mup2*(p1mup4*p3mup3*p5mup6*F(245)+p2mup4*p3mup3*p
     -   5mup6*F(245)+p3mup3*p3mup4*p5mup6*F(245)+p1mup6*p2mup3*p3mup
     -   4*F(246)+p1mup3*p2mup6*p3mup4*F(246)+p2mup3*p2mup6*p3mup4*F(
     -   246)+p1mup6*p3mup3*p3mup4*F(246)+p1mup3*p3mup4*p3mup6*F(246)
     -   +p3mup3*p3mup4*p3mup6*F(246)+p1mup4*((p1mup3+p2mup3)*p2mup6+
     -   p1mup6*(p2mup3+p3mup3)+(p1mup3+p3mup3)*p3mup6)*F(246)+p2mup4
     -   *((p1mup3+p2mup3)*p2mup6+p1mup6*(p2mup3+p3mup3)+(p1mup3+p3mu
     -   p3)*p3mup6)*F(246)-p3mup3*p4mup4*p5mup6*F(248)+((p1mup3+p2mu
     -   p3)*p2mup6+p1mup6*(p2mup3+p3mup3)+(p1mup3+p3mup3)*p3mup6)*p4
     -   mup4*F(249))+p2mup2*(p1mup4*p3mup3*p5mup6*F(245)+p2mup4*p3mu
     -   p3*p5mup6*F(245)+p3mup3*p3mup4*p5mup6*F(245)+p1mup6*p2mup3*p
     -   3mup4*F(246)+p1mup3*p2mup6*p3mup4*F(246)+p2mup3*p2mup6*p3mup
     -   4*F(246)+p1mup6*p3mup3*p3mup4*F(246)+p1mup3*p3mup4*p3mup6*F(
     -   246)+p3mup3*p3mup4*p3mup6*F(246)+p1mup4*((p1mup3+p2mup3)*p2m
     -   up6+p1mup6*(p2mup3+p3mup3)+(p1mup3+p3mup3)*p3mup6)*F(246)+p2
     -   mup4*((p1mup3+p2mup3)*p2mup6+p1mup6*(p2mup3+p3mup3)+(p1mup3+
     -   p3mup3)*p3mup6)*F(246)+p3mup3*p3mup4*p4mup6*F(247)-p3mup3*p4
     -   mup4*p5mup6*F(248)+p1mup6*p2mup3*p4mup4*F(249)+p1mup3*p2mup6
     -   *p4mup4*F(249)+p2mup3*p2mup6*p4mup4*F(249)+p1mup3*p3mup6*p4m
     -   up4*F(249)+p3mup3*p3mup6*p4mup4*F(249)+p3mup3*p4mup4*p4mup6*
     -   F(250))
       Fa14(184) = p1mup2*p1mup3*(p1mup4*p5mup6*F(245)+p2mup4*p5mup6*F
     -   (245)+p3mup4*p5mup6*F(245)+p1mup4*p4mup6*F(247)+p2mup4*p4mup
     -   6*F(247)+p3mup4*p4mup6*F(247)-p4mup4*p5mup6*F(248)+p4mup4*p4
     -   mup6*F(250))+p2mup2*(p1mup6*p3mup3*p4mup4*F(249)+p1mup3*(p1m
     -   up4*p5mup6*F(245)+p2mup4*p5mup6*F(245)+p3mup4*p5mup6*F(245)+
     -   p1mup4*p4mup6*F(247)+p2mup4*p4mup6*F(247)+p3mup4*p4mup6*F(24
     -   7)-p4mup4*p5mup6*F(248)+p4mup4*p4mup6*F(250)))
       Fa14(185) = p5mup4*(16*((p1mup3*p2mup2+p1mup2*(p1mup3+p2mup3+p3
     -   mup3))*p4mup6*F(253)+(p1mup3*p2mup2+p1mup2*(p1mup3+p2mup3))*
     -   p5mup6*F(256))+(p1mup2*(p1mup3+p2mup3+p3mup3)*(p1mup6+p2mup6
     -   +p3mup6)+p2mup2*(p1mup6*(p2mup3+p3mup3)+p1mup3*(p1mup6+p2mup
     -   6+p3mup6)))*F(269))
       Fa14(186) = p5mup4*(16*p1mup2*p3mup3*p5mup6*F(256)+p2mup2*(p2mu
     -   p3+p3mup3)*(16*(p4mup6*F(253)+p5mup6*F(256))+(p2mup6+p3mup6)
     -   *F(269)))
       Fa14(187) = (p1mup2+p2mup2)*(p1mup4*p2mup3*F(264)+p2mup3*p2mup4
     -   *F(264)+p1mup4*p3mup3*F(264)+p2mup4*p3mup3*F(264)+p2mup3*p3m
     -   up4*F(264)+p3mup3*p3mup4*F(264)-p2mup3*p5mup4*F(271)-p3mup3*
     -   p5mup4*F(271)+p2mup3*p4mup4*F(272)+p3mup3*p4mup4*F(272)+p1mu
     -   p3*((p1mup4+p2mup4+p3mup4)*F(264)-p5mup4*F(271)+p4mup4*F(272
     -   ))-p4mup3*p5mup4*s123*F(273)-p4mup3*p4mup4*s123*F(274)-(p1mu
     -   p4+p2mup4+p3mup4)*p4mup3*s123*F(275))
       Fa14(188) = 8*p4mup3*(p3mup4*p5mup6*F(245)+p1mup6*p3mup4*F(246)
     -   +p2mup6*p3mup4*F(246)+p3mup4*p3mup6*F(246)+p3mup4*p4mup6*F(2
     -   47)+p1mup4*(p5mup6*F(245)+(p1mup6+p2mup6+p3mup6)*F(246)+p4mu
     -   p6*F(247))+p2mup4*(p5mup6*F(245)+(p1mup6+p2mup6+p3mup6)*F(24
     -   6)+p4mup6*F(247))-p4mup4*p5mup6*F(248)+p1mup6*p4mup4*F(249)+
     -   p2mup6*p4mup4*F(249)+p3mup6*p4mup4*F(249)+p4mup4*p4mup6*F(25
     -   0)+p4mup6*p5mup4*F(253)+p5mup4*p5mup6*F(256)-(p1mup6+p2mup6+
     -   p3mup6)*p5mup4*F(257))-mup4mup6*p4mup3*F(276)
       Fa14(189) = (p1mup3+p2mup3+p3mup3)*(8*(p3mup4*p5mup6*F(245)+p1m
     -   up6*p3mup4*F(246)+p2mup6*p3mup4*F(246)+p3mup4*p3mup6*F(246)+
     -   p3mup4*p4mup6*F(247)+p1mup4*(p5mup6*F(245)+(p1mup6+p2mup6+p3
     -   mup6)*F(246)+p4mup6*F(247))+p2mup4*(p5mup6*F(245)+(p1mup6+p2
     -   mup6+p3mup6)*F(246)+p4mup6*F(247))-p4mup4*p5mup6*F(248)+p1mu
     -   p6*p4mup4*F(249)+p2mup6*p4mup4*F(249)+p3mup6*p4mup4*F(249)+p
     -   4mup4*p4mup6*F(250)+p4mup6*p5mup4*F(253)+p5mup4*p5mup6*F(256
     -   )+(p1mup6+p2mup6+p3mup6)*p5mup4*F(260))+mup4mup6*F(277))
       Fa14(190) = p1mup4*p2mup3*F(264)+p2mup3*p2mup4*F(264)+p1mup4*p3
     -   mup3*F(264)+p2mup4*p3mup3*F(264)+p2mup3*p3mup4*F(264)+p3mup3
     -   *p3mup4*F(264)-p2mup3*p5mup4*F(271)-p3mup3*p5mup4*F(271)+p2m
     -   up3*p4mup4*F(272)+p3mup3*p4mup4*F(272)+p1mup3*((p1mup4+p2mup
     -   4+p3mup4)*F(264)-p5mup4*F(271)+p4mup4*F(272))-p4mup3*p5mup4*
     -   s123*F(273)-p4mup3*p4mup4*s123*F(274)-(p1mup4+p2mup4+p3mup4)
     -   *p4mup3*s123*F(275)
       Fa14(191) = 8*p3mup2*(p3mup4*p5mup6*F(245)+p1mup6*p3mup4*F(246)
     -   +p2mup6*p3mup4*F(246)+p3mup4*p3mup6*F(246)+p3mup4*p4mup6*F(2
     -   47)+p1mup4*(p5mup6*F(245)+(p1mup6+p2mup6+p3mup6)*F(246)+p4mu
     -   p6*F(247))+p2mup4*(p5mup6*F(245)+(p1mup6+p2mup6+p3mup6)*F(24
     -   6)+p4mup6*F(247))-p4mup4*p5mup6*F(248)+p1mup6*p4mup4*F(249)+
     -   p2mup6*p4mup4*F(249)+p3mup6*p4mup4*F(249)+p4mup4*p4mup6*F(25
     -   0)+p4mup6*p5mup4*F(253)+p5mup4*p5mup6*F(256)-(p1mup6+p2mup6+
     -   p3mup6)*p5mup4*F(257))-mup4mup6*p3mup2*F(276)
       Fa14(192) = (p1mup2+p2mup2)*(8*(p3mup4*p5mup6*F(245)+p1mup6*p3m
     -   up4*F(246)+p2mup6*p3mup4*F(246)+p3mup4*p3mup6*F(246)+p3mup4*
     -   p4mup6*F(247)+p1mup4*(p5mup6*F(245)+(p1mup6+p2mup6+p3mup6)*F
     -   (246)+p4mup6*F(247))+p2mup4*(p5mup6*F(245)+(p1mup6+p2mup6+p3
     -   mup6)*F(246)+p4mup6*F(247))-p4mup4*p5mup6*F(248)+p1mup6*p4mu
     -   p4*F(249)+p2mup6*p4mup4*F(249)+p3mup6*p4mup4*F(249)+p4mup4*p
     -   4mup6*F(250)+p4mup6*p5mup4*F(253)+p5mup4*p5mup6*F(256)+(p1mu
     -   p6+p2mup6+p3mup6)*p5mup4*F(260))+mup4mup6*F(277))
       Fa14(193) = 4*(p1mup4*p3mup2*F(264)+p2mup4*p3mup2*F(264)+p3mup2
     -   *p3mup4*F(264)-p3mup2*p5mup4*F(271)-p4mup2*p5mup4*s123*F(273
     -   )-p4mup2*p4mup4*s123*F(274)-p1mup4*p4mup2*s123*F(275)-p2mup4
     -   *p4mup2*s123*F(275)-p3mup4*p4mup2*s123*F(275))+p3mup2*p4mup4
     -   *F(278)
       Fa14(194) = (p1mup2+p2mup2)*((p1mup4+p2mup4+p3mup4)*F(264)-p5mu
     -   p4*F(271)+p4mup4*F(272))
       Fa14(195) = 8*p4mup2*(p3mup4*p5mup6*F(245)+p1mup6*p3mup4*F(246)
     -   +p2mup6*p3mup4*F(246)+p3mup4*p3mup6*F(246)+p3mup4*p4mup6*F(2
     -   47)+p1mup4*(p5mup6*F(245)+(p1mup6+p2mup6+p3mup6)*F(246)+p4mu
     -   p6*F(247))+p2mup4*(p5mup6*F(245)+(p1mup6+p2mup6+p3mup6)*F(24
     -   6)+p4mup6*F(247))-p4mup4*p5mup6*F(248)+p1mup6*p4mup4*F(249)+
     -   p2mup6*p4mup4*F(249)+p3mup6*p4mup4*F(249)+p4mup4*p4mup6*F(25
     -   0)+p4mup6*p5mup4*F(253)+p5mup4*p5mup6*F(256)+(p1mup6+p2mup6+
     -   p3mup6)*p5mup4*F(260))+mup4mup6*p4mup2*F(277)
       Fa14(196) = -(p5mup4*p5mup6*F(251))-p1mup6*p5mup4*F(252)-p2mup6
     -   *p5mup4*F(252)-p3mup6*p5mup4*F(252)+4*(p3mup4*p5mup6*F(245)+
     -   p1mup6*p3mup4*F(246)+p2mup6*p3mup4*F(246)+p3mup4*p3mup6*F(24
     -   6)+p3mup4*p4mup6*F(247)+p1mup4*(p5mup6*F(245)+(p1mup6+p2mup6
     -   +p3mup6)*F(246)+p4mup6*F(247))+p2mup4*(p5mup6*F(245)+(p1mup6
     -   +p2mup6+p3mup6)*F(246)+p4mup6*F(247))-p4mup4*p5mup6*F(248)+p
     -   1mup6*p4mup4*F(249)+p2mup6*p4mup4*F(249)+p3mup6*p4mup4*F(249
     -   )+p4mup4*p4mup6*F(250)+p4mup6*p5mup4*F(253))-mup4mup6*F(279)
       Return
       End
