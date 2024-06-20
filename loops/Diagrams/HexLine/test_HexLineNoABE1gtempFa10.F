c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g10(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa10)
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
       Complex*16 Fa10(127:140)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa10(127) = 2*mup3mup6*s45*F(140)+4*(p4mup3*(p5mup6*F(149)+(p1m
     -   up6+p2mup6)*F(150)+p3mup6*F(151)+p4mup6*F(152))-p4mup6*(p3mu
     -   p3*F(153)+p5mup3*F(154)+(p1mup3+p2mup3)*F(155)))
       Fa10(128) = -4*(p1mup6*p2mup2*s45*F(144)+p2mup2*p2mup6*s45*F(14
     -   4)+p2mup2*p5mup6*s45*F(145)+p2mup2*p3mup6*s45*F(146)+p1mup2*
     -   s45*((p1mup6+p2mup6)*F(144)+p5mup6*F(145)+p3mup6*F(146))+p4m
     -   up2*p4mup6*s12*F(149)+p4mup2*p5mup6*s12*F(149)+p1mup6*p4mup2
     -   *s12*F(150)+p2mup6*p4mup2*s12*F(150)+p3mup6*p4mup2*s12*F(151
     -   ))+(p1mup2+p2mup2)*p4mup6*F(156)
       Fa10(129) = mup3mup6*F(158)+4*(-(p4mup3*p4mup6*F(135))-p4mup6*p
     -   5mup3*F(135)-p4mup3*p5mup6*F(135)-p5mup3*p5mup6*F(135)-p3mup
     -   3*p4mup6*F(159)-p3mup3*p5mup6*F(159)+p1mup3*(-((p4mup6+p5mup
     -   6)*F(130))+(p1mup6+p2mup6)*F(160)+p3mup6*F(161))+p2mup3*(-((
     -   p4mup6+p5mup6)*F(130))+(p1mup6+p2mup6)*F(160)+p3mup6*F(161))
     -   +p3mup3*p3mup6*F(162)+p1mup6*p3mup3*F(163)+p2mup6*p3mup3*F(1
     -   63)+p3mup6*p4mup3*F(164)+p3mup6*p5mup3*F(164)+(p1mup6+p2mup6
     -   )*(p4mup3+p5mup3)*F(165))
       Fa10(130) = p1mup2*((p1mup6+p2mup6)*F(120)+4*((p4mup6+p5mup6)*F
     -   (119)-p3mup6*F(121)))+p2mup2*((p1mup6+p2mup6)*F(120)+4*((p4m
     -   up6+p5mup6)*F(119)-p3mup6*F(121)))+4*p3mup2*s12*((p4mup6+p5m
     -   up6)*F(149)+(p1mup6+p2mup6)*F(150)+p3mup6*F(151))
       Fa10(131) = -8*(p4mup3*p4mup6*p5mup4*F(135)+p4mup6*p5mup3*p5mup
     -   4*F(135)+p4mup3*p5mup4*p5mup6*F(135)+p5mup3*p5mup4*p5mup6*F(
     -   135)+p3mup3*p4mup6*p5mup4*F(159)+p3mup3*p5mup4*p5mup6*F(159)
     -   -p3mup3*p3mup6*p5mup4*F(162)-p1mup6*p3mup3*p5mup4*F(163)-p2m
     -   up6*p3mup3*p5mup4*F(163)-p3mup6*p4mup3*p5mup4*F(164)-p3mup6*
     -   p5mup3*p5mup4*F(164)-p1mup6*p4mup3*p5mup4*F(165)-p2mup6*p4mu
     -   p3*p5mup4*F(165)-p1mup6*p5mup3*p5mup4*F(165)-p2mup6*p5mup3*p
     -   5mup4*F(165))-mup3mup6*p5mup4*F(167)+4*mup4mup6*p4mup3*F(168
     -   )+p1mup3*(8*p5mup4*(-((p4mup6+p5mup6)*F(130))+(p1mup6+p2mup6
     -   )*F(160)+p3mup6*F(161))+12*mup4mup6*s45*F(169))+p2mup3*(8*p5
     -   mup4*(-((p4mup6+p5mup6)*F(130))+(p1mup6+p2mup6)*F(160)+p3mup
     -   6*F(161))+12*mup4mup6*s45*F(169))+12*(mup4mup6*p3mup3*s45*F(
     -   170)+mup4mup6*p5mup3*s45*F(171))
       Fa10(132) = 4*p3mup2*s12*(3*mup4mup6*s45*F(148)-2*p5mup4*((p4mu
     -   p6+p5mup6)*F(149)+(p1mup6+p2mup6)*F(150)+p3mup6*F(151)))+mup
     -   4mup6*p4mup2*s12*F(173)+p1mup2*(mup4mup6*s45*F(172)+8*p5mup4
     -   *(-((p4mup6+p5mup6)*F(119))+p3mup6*F(143)+(p1mup6+p2mup6)*F(
     -   174)))+p2mup2*(mup4mup6*s45*F(172)+8*p5mup4*(-((p4mup6+p5mup
     -   6)*F(119))+p3mup6*F(143)+(p1mup6+p2mup6)*F(174)))
       Fa10(133) = -4*mup3mup6*p3mup2*s12*s45*F(140)+mup3mup6*p4mup2*s
     -   12*F(167)+mup3mup6*p2mup2*s45*F(175)+8*(-(p2mup2*p2mup3*p4mu
     -   p6*F(127))-p2mup2*p4mup3*p5mup6*F(129)+p2mup3*p4mup2*p4mup6*
     -   s12*F(130)+p2mup3*p4mup2*p5mup6*s12*F(130)+p1mup6*p2mup2*p2m
     -   up3*s45*F(130)+p2mup2*p2mup3*p2mup6*s45*F(130)+p1mup6*p2mup2
     -   *p3mup3*s45*F(131)+p2mup2*p2mup6*p3mup3*s45*F(131)+p2mup2*p3
     -   mup3*p3mup6*s45*F(132)+p2mup2*p2mup3*p3mup6*s45*F(133)+p2mup
     -   2*p3mup6*p5mup3*s45*F(134)+p4mup2*p4mup3*p4mup6*s12*F(135)+p
     -   4mup2*p4mup6*p5mup3*s12*F(135)+p4mup2*p4mup3*p5mup6*s12*F(13
     -   5)+p4mup2*p5mup3*p5mup6*s12*F(135)+p1mup6*p2mup2*p5mup3*s45*
     -   F(135)+p2mup2*p2mup6*p5mup3*s45*F(135)-p2mup2*p2mup3*p5mup6*
     -   s45*F(136)-p2mup2*p5mup3*p5mup6*s45*F(137)-p2mup2*p3mup3*p5m
     -   up6*s45*F(138)-p3mup2*p4mup3*p5mup6*s12*F(149)-p1mup6*p3mup2
     -   *p4mup3*s12*F(150)-p2mup6*p3mup2*p4mup3*s12*F(150)-p3mup2*p3
     -   mup6*p4mup3*s12*F(151)-p3mup2*p4mup3*p4mup6*s12*F(152)+p3mup
     -   2*p3mup3*p4mup6*s12*F(153)+p3mup2*p4mup6*p5mup3*s12*F(154)+p
     -   2mup3*p3mup2*p4mup6*s12*F(155)+p3mup3*p4mup2*p4mup6*s12*F(15
     -   9)+p3mup3*p4mup2*p5mup6*s12*F(159)-p1mup6*p2mup3*p4mup2*s12*
     -   F(160)-p2mup3*p2mup6*p4mup2*s12*F(160)-p2mup3*p3mup6*p4mup2*
     -   s12*F(161)-p1mup3*(p2mup2*(p4mup6*F(127)-s45*((p1mup6+p2mup6
     -   )*F(130)+p3mup6*F(133)-p5mup6*F(136)))+s12*(-(p3mup2*p4mup6*
     -   F(155))+p4mup2*(-((p4mup6+p5mup6)*F(130))+(p1mup6+p2mup6)*F(
     -   160)+p3mup6*F(161))))-p3mup3*p3mup6*p4mup2*s12*F(162)-p1mup6
     -   *p3mup3*p4mup2*s12*F(163)-p2mup6*p3mup3*p4mup2*s12*F(163)-p3
     -   mup6*p4mup2*p4mup3*s12*F(164)-p3mup6*p4mup2*p5mup3*s12*F(164
     -   )-p1mup6*p4mup2*p4mup3*s12*F(165)-p2mup6*p4mup2*p4mup3*s12*F
     -   (165)-p1mup6*p4mup2*p5mup3*s12*F(165)-p2mup6*p4mup2*p5mup3*s
     -   12*F(165)-p2mup2*p4mup3*p4mup6*F(176)-p2mup2*p3mup3*p4mup6*F
     -   (177)+p2mup2*p3mup6*p4mup3*F(178)+p2mup2*p4mup6*p5mup3*F(179
     -   )+p2mup2*(p1mup6+p2mup6)*p4mup3*F(180))+p1mup2*(mup3mup6*s45
     -   *F(175)+8*(-(p4mup3*p5mup6*F(129))+s45*(p1mup6*p3mup3*F(131)
     -   +p2mup6*p3mup3*F(131)+p3mup3*p3mup6*F(132)+p3mup6*p5mup3*F(1
     -   34)+(p1mup6+p2mup6)*p5mup3*F(135))+p1mup3*(-(p4mup6*F(127))+
     -   s45*((p1mup6+p2mup6)*F(130)+p3mup6*F(133)-p5mup6*F(136)))+p2
     -   mup3*(-(p4mup6*F(127))+s45*((p1mup6+p2mup6)*F(130)+p3mup6*F(
     -   133)-p5mup6*F(136)))-p5mup3*p5mup6*s45*F(137)-p3mup3*p5mup6*
     -   s45*F(138)-p4mup3*p4mup6*F(176)-p3mup3*p4mup6*F(177)+p3mup6*
     -   p4mup3*F(178)+p4mup6*p5mup3*F(179)+(p1mup6+p2mup6)*p4mup3*F(
     -   180)))
       Fa10(134) = (6*s45*((p1mup3+p2mup3)*F(169)+p3mup3*F(170)+p5mup3
     -   *F(171))+p4mup3*F(181))/s45
       Fa10(135) = ((p1mup2+p2mup2)*F(118))/s12-6*p3mup2*F(148)+(p4mup
     -   2*F(157))/s45
       Fa10(136) = 8*mup4mup6*(p4mup2*((p4mup3+p5mup3)*F(182)-p3mup3*F
     -   (184)-(p1mup3+p2mup3)*F(185))+p3mup2*(p4mup3*F(183)+s45*((p1
     -   mup3+p2mup3)*F(190)+p3mup3*F(191)+p5mup3*F(195))))+p3mup2*p5
     -   mup4*(8*mup3mup6*F(186)+p2mup6*p4mup3*F(197)+p2mup6*p5mup3*F
     -   (197)+p1mup6*(p4mup3+p5mup3)*F(197)+16*(p2mup3*p2mup6*F(187)
     -   +p2mup3*p3mup6*F(188)+p2mup3*p4mup6*F(189)+p2mup3*p5mup6*F(1
     -   89)+p1mup3*((p1mup6+p2mup6)*F(187)+p3mup6*F(188)+(p4mup6+p5m
     -   up6)*F(189))+p2mup6*p3mup3*F(192)+p1mup6*(p2mup3*F(187)+p3mu
     -   p3*F(192))+p3mup3*p3mup6*F(193)-p3mup3*p4mup6*F(194)-p3mup3*
     -   p5mup6*F(194)+p4mup3*p4mup6*F(196)+p4mup6*p5mup3*F(196)+p4mu
     -   p3*p5mup6*F(196)+p5mup3*p5mup6*F(196)+p3mup6*(p4mup3+p5mup3)
     -   *F(198)))
       Fa10(137) = 4*p3mup2*(p4mup3*F(183)+s45*((p1mup3+p2mup3)*F(190)
     -   +p3mup3*F(191)+p5mup3*F(195)))-p4mup2*(4*((p1mup3+p2mup3)*F(
     -   185)-(p4mup3+p5mup3)*F(199))+p3mup3*F(200))
       Fa10(138) = 8*mup4mup6*(-(p3mup3*p4mup2*s12*F(153))+p4mup2*s12*
     -   (-((p4mup3+p5mup3)*F(154))-(p1mup3+p2mup3)*F(155))+(p1mup2+p
     -   2mup2)*p3mup3*s45*F(191)+(p1mup2+p2mup2)*(p4mup3*F(183)+(p1m
     -   up3+p2mup3)*s45*F(190)+p5mup3*s45*F(195)))+(p1mup2+p2mup2)*p
     -   5mup4*(-(mup3mup6*F(201))+16*(p1mup6*p2mup3*F(187)+p2mup3*p2
     -   mup6*F(187)+p2mup3*p3mup6*F(188)+p2mup3*p4mup6*F(189)+p2mup3
     -   *p5mup6*F(189)+p1mup3*((p1mup6+p2mup6)*F(187)+p3mup6*F(188)+
     -   (p4mup6+p5mup6)*F(189))+p1mup6*p3mup3*F(192)+p2mup6*p3mup3*F
     -   (192)+p3mup3*p3mup6*F(193)-p3mup3*p4mup6*F(194)-p3mup3*p5mup
     -   6*F(194)+p4mup3*p4mup6*F(196)+p4mup6*p5mup3*F(196)+p4mup3*p5
     -   mup6*F(196)+p5mup3*p5mup6*F(196)+p3mup6*p4mup3*F(198)+p3mup6
     -   *p5mup3*F(198)-(p1mup6+p2mup6)*(p4mup3+p5mup3)*F(202)))
       Fa10(139) = -4*(p3mup3*p4mup2*s12*F(153)+p4mup2*(p4mup3+p5mup3)
     -   *s12*F(154)+(p1mup3+p2mup3)*p4mup2*s12*F(155)-(p1mup2+p2mup2
     -   )*p3mup3*s45*F(191))+(p1mup2+p2mup2)*(4*s45*((p1mup3+p2mup3)
     -   *F(190)+p5mup3*F(195))+p4mup3*F(203))
       Fa10(140) = mup4mup6*(p3mup2*s12*(p3mup3*F(153)+(p4mup3+p5mup3)
     -   *F(154)+(p1mup3+p2mup3)*F(155))+(p1mup2+p2mup2)*(-(p3mup3*F(
     -   184))-(p1mup3+p2mup3)*F(185)+(p4mup3+p5mup3)*F(199)))
       Return
       End
