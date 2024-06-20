c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g12(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa12)
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
       Complex*16 Fa12(155:168)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa12(155) = 4*(p1mup6*p2mup3*F(187)+p2mup3*p2mup6*F(187)+p2mup3
     -   *p3mup6*F(188)+p2mup3*p4mup6*F(189)+p2mup3*p5mup6*F(189)+p1m
     -   up3*((p1mup6+p2mup6)*F(187)+p3mup6*F(188)+(p4mup6+p5mup6)*F(
     -   189))+p1mup6*p3mup3*F(192)+p2mup6*p3mup3*F(192)+p3mup3*p3mup
     -   6*F(193)-p3mup3*p4mup6*F(194)-p3mup3*p5mup6*F(194)+p3mup6*p4
     -   mup3*F(198)+p3mup6*p5mup3*F(198))+mup3mup6*F(217)+p4mup3*p4m
     -   up6*F(218)+p4mup6*p5mup3*F(218)+p4mup3*p5mup6*F(218)+p5mup3*
     -   p5mup6*F(218)+(p1mup6+p2mup6)*(p4mup3+p5mup3)*F(219)
       Fa12(156) = 2*(p4mup6+p5mup6)*F(119)+(p1mup6+p2mup6)*F(220)+p3m
     -   up6*F(221)
       Fa12(157) = (p4mup6+p5mup6)*F(149)+(p1mup6+p2mup6)*F(150)+p3mup
     -   6*F(151)
       Fa12(158) = mup4mup6*F(222)+2*(-(p3mup4*p5mup6*F(223))-p4mup4*p
     -   5mup6*F(224)+p1mup6*p3mup4*F(225)+p2mup6*p3mup4*F(225)+p3mup
     -   4*p3mup6*F(225)+p3mup4*p4mup6*F(226)+p1mup4*(-(p5mup6*F(223)
     -   )+(p1mup6+p2mup6+p3mup6)*F(225)+p4mup6*F(226))+p2mup4*(-(p5m
     -   up6*F(223))+(p1mup6+p2mup6+p3mup6)*F(225)+p4mup6*F(226))+p1m
     -   up6*p4mup4*F(227)+p2mup6*p4mup4*F(227)+p3mup6*p4mup4*F(227)+
     -   p4mup4*p4mup6*F(228)-p5mup4*p5mup6*F(229)+p1mup6*p5mup4*F(23
     -   0)+p2mup6*p5mup4*F(230)+p3mup6*p5mup4*F(230)+p4mup6*p5mup4*F
     -   (231))
       Fa12(159) = p2mup3*p2mup6*F(232)+p2mup6*p3mup3*F(232)+p1mup6*(p
     -   2mup3+p3mup3)*F(232)+p2mup3*p3mup6*F(232)+p3mup3*p3mup6*F(23
     -   2)+p1mup3*((p1mup6+p2mup6+p3mup6)*F(232)+4*(p4mup6*F(233)-p5
     -   mup6*F(234)))+4*(p2mup3*p4mup6*F(233)+p3mup3*p4mup6*F(233)-p
     -   2mup3*p5mup6*F(234)-p3mup3*p5mup6*F(234)-p4mup3*p4mup6*s123*
     -   F(235)-p1mup6*p4mup3*s123*F(236)-p2mup6*p4mup3*s123*F(236)-p
     -   3mup6*p4mup3*s123*F(236)-p4mup3*p5mup6*s123*F(237))
       Fa12(160) = (p1mup6+p2mup6+p3mup6)*F(238)+p4mup6*F(239)+p5mup6*
     -   F(240)
       Fa12(161) = 4*p4mup2*s123*(p4mup6*F(235)+(p1mup6+p2mup6+p3mup6)
     -   *F(236)+p5mup6*F(237))+p3mup2*(4*(p5mup6*F(234)-p4mup6*F(241
     -   ))+(p1mup6+p2mup6+p3mup6)*F(242))
       Fa12(162) = (p1mup2+p2mup2)*((p1mup6+p2mup6+p3mup6)*F(232)+4*(p
     -   4mup6*F(233)-p5mup6*F(234)))
       Fa12(163) = -(mup4mup6*F(244))-p5mup4*p5mup6*F(251)-p1mup6*p5mu
     -   p4*F(252)-p2mup6*p5mup4*F(252)-p3mup6*p5mup4*F(252)+4*(p3mup
     -   4*p5mup6*F(245)+p1mup6*p3mup4*F(246)+p2mup6*p3mup4*F(246)+p3
     -   mup4*p3mup6*F(246)+p3mup4*p4mup6*F(247)+p1mup4*(p5mup6*F(245
     -   )+(p1mup6+p2mup6+p3mup6)*F(246)+p4mup6*F(247))+p2mup4*(p5mup
     -   6*F(245)+(p1mup6+p2mup6+p3mup6)*F(246)+p4mup6*F(247))-p4mup4
     -   *p5mup6*F(248)+p1mup6*p4mup4*F(249)+p2mup6*p4mup4*F(249)+p3m
     -   up6*p4mup4*F(249)+p4mup4*p4mup6*F(250)+p4mup6*p5mup4*F(253))
       Fa12(164) = (p1mup2+p2mup2)*(p4mup6*F(235)+(p1mup6+p2mup6+p3mup
     -   6)*F(236)+p5mup6*F(237))
       Fa12(165) = p4mup3*(-(mup4mup6*F(255))+8*(p1mup4*(p5mup6*F(245)
     -   +(p1mup6+p2mup6+p3mup6)*F(246)+p4mup6*F(247))+p2mup4*(p5mup6
     -   *F(245)+(p1mup6+p2mup6+p3mup6)*F(246)+p4mup6*F(247))+p3mup4*
     -   (p5mup6*F(245)+(p1mup6+p2mup6+p3mup6)*F(246)+p4mup6*F(247))+
     -   p4mup4*(-(p5mup6*F(248))+(p1mup6+p2mup6+p3mup6)*F(249)+p4mup
     -   6*F(250))+p5mup4*(p4mup6*F(253)+p5mup6*F(256)-(p1mup6+p2mup6
     -   +p3mup6)*F(257))))
       Fa12(166) = 4*(mup4mup6*p2mup2*s123*F(258)+p1mup2*s123*(2*(p1mu
     -   p4*(-(p5mup6*F(223))+(p1mup6+p2mup6+p3mup6)*F(225)+p4mup6*F(
     -   226))+p2mup4*(-(p5mup6*F(223))+(p1mup6+p2mup6+p3mup6)*F(225)
     -   +p4mup6*F(226))+p3mup4*(-(p5mup6*F(223))+(p1mup6+p2mup6+p3mu
     -   p6)*F(225)+p4mup6*F(226))+p4mup4*(-(p5mup6*F(224))+(p1mup6+p
     -   2mup6+p3mup6)*F(227)+p4mup6*F(228))+p5mup4*(-(p5mup6*F(229))
     -   +(p1mup6+p2mup6+p3mup6)*F(230)+p4mup6*F(231)))+mup4mup6*F(25
     -   8)))+mup4mup6*p4mup2*s12*F(259)+8*(-(p2mup2*p2mup4*p5mup6*s1
     -   23*F(223))-p2mup2*p3mup4*p5mup6*s123*F(223)-p2mup2*p4mup4*p5
     -   mup6*s123*F(224)+p1mup6*p2mup2*p2mup4*s123*F(225)+p2mup2*p2m
     -   up4*p2mup6*s123*F(225)+p1mup6*p2mup2*p3mup4*s123*F(225)+p2mu
     -   p2*p2mup6*p3mup4*s123*F(225)+p2mup2*p2mup4*p3mup6*s123*F(225
     -   )+p2mup2*p3mup4*p3mup6*s123*F(225)+p2mup2*p2mup4*p4mup6*s123
     -   *F(226)+p2mup2*p3mup4*p4mup6*s123*F(226)+p1mup6*p2mup2*p4mup
     -   4*s123*F(227)+p2mup2*p2mup6*p4mup4*s123*F(227)+p2mup2*p3mup6
     -   *p4mup4*s123*F(227)+p2mup2*p4mup4*p4mup6*s123*F(228)-p2mup2*
     -   p5mup4*p5mup6*s123*F(229)+p1mup6*p2mup2*p5mup4*s123*F(230)+p
     -   2mup2*p2mup6*p5mup4*s123*F(230)+p2mup2*p3mup6*p5mup4*s123*F(
     -   230)+p2mup2*p4mup6*p5mup4*s123*F(231)+p2mup4*p4mup2*p5mup6*s
     -   12*F(245)+p3mup4*p4mup2*p5mup6*s12*F(245)+p1mup6*p2mup4*p4mu
     -   p2*s12*F(246)+p2mup4*p2mup6*p4mup2*s12*F(246)+p1mup6*p3mup4*
     -   p4mup2*s12*F(246)+p2mup6*p3mup4*p4mup2*s12*F(246)+p2mup4*p3m
     -   up6*p4mup2*s12*F(246)+p3mup4*p3mup6*p4mup2*s12*F(246)+p2mup4
     -   *p4mup2*p4mup6*s12*F(247)+p3mup4*p4mup2*p4mup6*s12*F(247)+p1
     -   mup4*(p2mup2*s123*(-(p5mup6*F(223))+(p1mup6+p2mup6+p3mup6)*F
     -   (225)+p4mup6*F(226))+p4mup2*s12*(p5mup6*F(245)+(p1mup6+p2mup
     -   6+p3mup6)*F(246)+p4mup6*F(247)))-p4mup2*p4mup4*p5mup6*s12*F(
     -   248)+p1mup6*p4mup2*p4mup4*s12*F(249)+p2mup6*p4mup2*p4mup4*s1
     -   2*F(249)+p3mup6*p4mup2*p4mup4*s12*F(249)+p4mup2*p4mup4*p4mup
     -   6*s12*F(250)+p4mup2*p4mup6*p5mup4*s12*F(253)+p4mup2*p5mup4*p
     -   5mup6*s12*F(256)+(p1mup6+p2mup6+p3mup6)*p4mup2*p5mup4*s12*F(
     -   260))
       Fa12(167) = (p1mup2+p2mup2)*(p1mup3*p5mup6*F(234)+p2mup3*p5mup6
     -   *F(234)+p3mup3*p5mup6*F(234)+p4mup3*s123*(p4mup6*F(235)+(p1m
     -   up6+p2mup6+p3mup6)*F(236))+p4mup3*p5mup6*s123*F(237)-p1mup3*
     -   p4mup6*F(241)-p2mup3*p4mup6*F(241)-p3mup3*p4mup6*F(241)-(p1m
     -   up3+p2mup3+p3mup3)*(p1mup6+p2mup6+p3mup6)*F(261))
       Fa12(168) = p5mup4*F(262)+p4mup4*F(263)-2*(p1mup4+p2mup4+p3mup4
     -   )*F(264)
       Return
       End
