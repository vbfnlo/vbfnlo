c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa3(p1mup2,p1mup3,p1mup5,p1mup6,p2mup
     -   2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4
     -   mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3
     -   ,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,
     -   mup5mup6,Fa3)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe2gaHexLineFhlFunctions/F
       COMMON/NoAbe2gaHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,
     -   s12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(872),F(9706)
      Real*8 P(9687) 
        COMMON/NoAbe2gaHexLinePFunctions/P
       Complex*16 Fa3(49:72)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa3(49) = p4mup3*p6mup5*(p1mup6*(p2mup2+p3mup2+p4mup2+p5mup2)*F
     -   (154)-4*(p2mup6*(p2mup2+p3mup2)+p2mup2*p3mup6)*F(9570)+p1mup
     -   2*(p1mup6*F(154)+p5mup6*F(155)-4*(p2mup6+p3mup6+p4mup6)*F(95
     -   70)))
       Fa3(50) = p4mup3*p6mup5*(-((p2mup2+p3mup2+p4mup2+p5mup2)*p5mup6
     -   *F(155))+4*(p2mup2*p4mup6+p3mup2*(p3mup6+p4mup6)+(p2mup6+p3m
     -   up6+p4mup6)*(p4mup2+p5mup2))*F(9570))
       Fa3(51) = p4mup3*p6mup2*p6mup5*(-(p1mup6*F(154))-p5mup6*F(155)+
     -   4*(p2mup6+p3mup6+p4mup6)*F(9570))
       Fa3(52) = p4mup3*p6mup6*(p1mup5*(p2mup2+p3mup2+p4mup2+p5mup2)*F
     -   (156)+2*p2mup2*p2mup5*F(9570)+p1mup2*(p1mup5*F(156)-p5mup5*F
     -   (157)+2*(p2mup5+p3mup5+p4mup5)*F(9570)))
       Fa3(53) = p4mup3*p6mup6*(p1mup5*p6mup2*F(156)-(p2mup2+p3mup2+p4
     -   mup2+p5mup2)*p5mup5*F(157)+2*((p3mup5+p4mup5)*(p2mup2+p3mup2
     -   +p4mup2+p5mup2)+p4mup5*p6mup2+p2mup5*(p3mup2+p4mup2+p5mup2+p
     -   6mup2))*F(9570))
       Fa3(54) = p4mup3*p6mup6*(p5mup5*p6mup2*F(157)+2*((p1mup2+p2mup2
     -   +p3mup2+p4mup2+p5mup2+p6mup2)*p6mup5*F(9568)-p3mup5*p6mup2*F
     -   (9570)))
       Fa3(55) = p4mup3*(p1mup2+p2mup2+p3mup2+p4mup2+p5mup2+p6mup2)
       Fa3(56) = mup5mup6*F(160)+p5mup5*p6mup6*F(166)+p5mup6*p6mup5*F(
     -   167)-(p2mup5+p3mup5+p4mup5)*p5mup6*F(169)+p2mup5*p6mup6*F(17
     -   0)+p3mup5*p6mup6*F(170)+p4mup5*p6mup6*F(170)+p1mup5*(p1mup6*
     -   F(161)-p5mup6*F(162)+2*(p6mup6*F(163)+(p2mup6+p3mup6+p4mup6)
     -   *F(171)))+p1mup6*(p5mup5*F(164)+2*p6mup5*F(165)+(p2mup5+p3mu
     -   p5+p4mup5)*F(172))+p2mup6*p6mup5*F(174)+p3mup6*p6mup5*F(174)
     -   +p4mup6*p6mup5*F(174)+2*(p6mup5*p6mup6*F(168)+p2mup6*p5mup5*
     -   F(173)+p3mup6*p5mup5*F(173)+p4mup6*p5mup5*F(173)-p5mup5*p5mu
     -   p6*F(175))+8*(p2mup5+p3mup5+p4mup5)*(p2mup6+p3mup6+p4mup6)*F
     -   (176)
       Fa3(57) = mup3mup6*(p1mup5*F(178)+p6mup5*F(182)+2*(p5mup5*F(180
     -   )-(p2mup5+p3mup5+p4mup5)*F(184)))+mup3mup5*(p6mup6*F(183)+2*
     -   (p1mup6*F(179)-p5mup6*F(181)-(p2mup6+p3mup6+p4mup6)*F(185)))
       Fa3(58) = mup3mup6*(p6mup5*F(190)+2*(p1mup5*F(186)+p5mup5*F(188
     -   )-(p2mup5+p3mup5+p4mup5)*F(192)))+mup3mup5*(p1mup6*F(187)-p5
     -   mup6*F(189)+p6mup6*F(191)+(p2mup6+p3mup6+p4mup6)*F(193))
       Fa3(59) = mup5mup6*F(194)+p2mup5*(p1mup6*F(196)+p5mup6*F(197)+p
     -   6mup6*F(198))+p2mup6*(p5mup5*F(199)+p6mup5*F(200)+p2mup5*F(2
     -   01)+2*(p3mup5+p4mup5)*F(202))+p2mup5*p3mup6*F(203)+p2mup5*p4
     -   mup6*F(203)+p1mup6*p5mup5*F(205)+p5mup5*p5mup6*F(209)+p5mup5
     -   *p6mup6*F(210)+2*(p1mup6*p6mup5*F(206)+p6mup5*p6mup6*F(211)+
     -   p1mup6*p3mup5*F(212)+p1mup6*p4mup5*F(212))+p1mup5*(2*p2mup6*
     -   F(195)+p6mup6*F(204)+p1mup6*F(207)-4*p5mup6*F(208)+(p3mup6+p
     -   4mup6)*F(213))+(p3mup5+p4mup5)*p6mup6*F(215)+p3mup6*p5mup5*F
     -   (216)+p4mup6*p5mup5*F(216)+p3mup6*p6mup5*F(217)+p4mup6*p6mup
     -   5*F(217)+p5mup6*p6mup5*F(218)+8*(p3mup5*p5mup6*F(214)+p4mup5
     -   *p5mup6*F(214)+(p3mup5+p4mup5)*(p3mup6+p4mup6)*F(219))
       Fa3(60) = mup5mup6*F(220)+p2mup5*(p1mup6*F(222)+p5mup6*F(223)+p
     -   6mup6*F(224))+p2mup6*(p6mup5*F(226)+2*(p5mup5*F(225)+p2mup5*
     -   F(227))+(p3mup5+p4mup5)*F(228))+p1mup6*p5mup5*F(231)+p5mup5*
     -   p5mup6*F(235)+2*(p1mup6*p6mup5*F(232)-p5mup5*p6mup6*F(236)+p
     -   6mup5*p6mup6*F(237)+p1mup6*p3mup5*F(238)+p1mup6*p4mup5*F(238
     -   ))+p1mup5*(2*(p2mup6*F(221)+p6mup6*F(230))+p1mup6*F(233)+p5m
     -   up6*F(234)+(p3mup6+p4mup6)*F(239))+p5mup6*p6mup5*F(240)+p3mu
     -   p5*p5mup6*F(241)+p4mup5*p5mup6*F(241)+4*(p3mup5+p4mup5)*p6mu
     -   p6*F(242)+(p3mup6+p4mup6)*(p2mup5*F(229)+p5mup5*F(243)+8*(p6
     -   mup5*F(244)+(p3mup5+p4mup5)*F(245)))
       Fa3(61) = p1mup6*F(58)+p5mup6*F(246)+p6mup6*F(247)+p2mup6*F(248
     -   )+(p3mup6+p4mup6)*F(249)
       Fa3(62) = -2*mup2mup6*F(84)-(p1mup2+p2mup2+p3mup2+p4mup2+p5mup2
     -   +p6mup2)*(2*(p1mup6*F(11)-p5mup6*F(12)-(p2mup6+p3mup6+p4mup6
     -   )*F(14))+p6mup6*F(250))
       Fa3(63) = 4*(p1mup2*p2mup6*F(97)+p2mup2*p2mup6*F(97)+p2mup6*p5m
     -   up2*F(97)+p2mup6*p6mup2*F(97))+p1mup6*p3mup2*F(251)+p1mup6*p
     -   4mup2*F(251)+p3mup2*p6mup6*F(253)+p4mup2*p6mup6*F(253)-2*(p2
     -   mup6*p3mup2*F(98)+p2mup6*p4mup2*F(98)-mup2mup6*F(105)+p3mup2
     -   *p5mup6*F(252)+p4mup2*p5mup6*F(252)+p3mup2*p3mup6*F(254)+p3m
     -   up6*p4mup2*F(254)+p3mup2*p4mup6*F(254)+p4mup2*p4mup6*F(254)+
     -   (p1mup2+p2mup2+p5mup2+p6mup2)*(-(p6mup6*F(255))+2*(p5mup6*F(
     -   100)-(p3mup6+p4mup6)*F(256)-p1mup6*F(257)))*P(59))
       Fa3(64) = p4mup2*p6mup6*F(258)-2*(p1mup6*(p3mup2+p4mup2)*F(127)
     -   +p4mup2*p5mup6*F(128)+p2mup6*p4mup2*F(134)-mup2mup6*F(141)+p
     -   3mup6*p4mup2*F(259)+p4mup2*p4mup6*F(259))+p3mup2*(p6mup6*F(2
     -   58)-2*(p5mup6*F(128)+p2mup6*F(134)+(p3mup6+p4mup6)*F(259)))+
     -   p1mup2*p6mup6*F(261)+p2mup2*p6mup6*F(261)+p5mup2*p6mup6*F(26
     -   1)+p6mup2*p6mup6*F(261)+4*(p2mup6*(p1mup2+p2mup2+p5mup2+p6mu
     -   p2)*F(133)+p1mup2*p5mup6*F(260)+p2mup2*p5mup6*F(260)+p5mup2*
     -   p5mup6*F(260)+p5mup6*p6mup2*F(260)+(p3mup6+p4mup6)*(p1mup2+p
     -   2mup2+p5mup2+p6mup2)*F(262)+p1mup6*(p1mup2+p2mup2+p5mup2+p6m
     -   up2)*F(263))
       Fa3(65) = p1mup5*F(264)+p5mup5*F(265)+p6mup5*F(266)+p2mup5*F(26
     -   7)+(p3mup5+p4mup5)*F(268)
       Fa3(66) = mup5mup6*(p2mup2+p3mup2+p4mup2+p5mup2)
       Fa3(67) = -4*(p2mup5*p3mup2*p5mup6*F(113)+p2mup5*p4mup2*p5mup6*
     -   F(113)+p2mup5*p5mup2*p5mup6*F(113))+p1mup5*p1mup6*p3mup2*F(1
     -   14)+p1mup5*p1mup6*p4mup2*F(114)-p1mup6*p2mup5*p3mup2*F(116)-
     -   p1mup6*p2mup5*p4mup2*F(116)+p1mup6*p3mup2*p5mup5*F(117)+p1mu
     -   p6*p4mup2*p5mup5*F(117)-p1mup6*p3mup2*p6mup5*F(119)-p1mup6*p
     -   4mup2*p6mup5*F(119)-p1mup5*p3mup2*p6mup6*F(120)-p1mup5*p4mup
     -   2*p6mup6*F(120)-p3mup2*p5mup5*p6mup6*F(121)-p4mup2*p5mup5*p6
     -   mup6*F(121)+p2mup5*p3mup2*p6mup6*F(122)+p2mup5*p4mup2*p6mup6
     -   *F(122)+p2mup2*(p1mup5*(p1mup6*F(114)+2*p5mup6*F(118))+2*p2m
     -   up6*(p5mup5*F(110)-2*p2mup5*F(112)+p1mup5*F(115)-p6mup5*F(95
     -   72)))-2*(p3mup2*p5mup5*p5mup6*F(111)+p4mup2*p5mup5*p5mup6*F(
     -   111)+p5mup2*p5mup5*p5mup6*F(111)-p1mup5*p3mup2*p5mup6*F(118)
     -   -p1mup5*p4mup2*p5mup6*F(118)-p1mup5*p5mup2*p5mup6*F(118)-p3m
     -   up2*p6mup5*p6mup6*F(9569)-p4mup2*p6mup5*p6mup6*F(9569)+(p3mu
     -   p2+p4mup2+p5mup2)*p5mup6*p6mup5*F(9571)-p2mup6*(p3mup2+p4mup
     -   2)*(p5mup5*F(110)-2*p2mup5*F(112)+p1mup5*F(115)-p6mup5*F(957
     -   2)))
       Fa3(68) = -(mup5mup6*p6mup2*F(271))+mup2mup6*p6mup5*F(274)+2*(m
     -   up2mup6*p1mup5*F(272)+mup2mup6*p5mup5*F(273)+mup2mup6*p2mup5
     -   *F(275))+4*(p1mup5*p1mup6*p5mup2*F(114)-p1mup6*p2mup5*p5mup2
     -   *F(116)+p1mup6*p5mup2*p5mup5*F(117)-p1mup6*p5mup2*p6mup5*F(1
     -   19)-p1mup5*p5mup2*p6mup6*F(120)-p5mup2*p5mup5*p6mup6*F(121)+
     -   p2mup5*p5mup2*p6mup6*F(122)-p2mup2*(4*p2mup5*p5mup6*F(113)+p
     -   1mup6*p2mup5*F(116)-p1mup6*p5mup5*F(117)+p1mup6*p6mup5*F(119
     -   )+p1mup5*p6mup6*F(120)+p5mup5*p6mup6*F(121)-p2mup5*p6mup6*F(
     -   122)+2*(p5mup5*p5mup6*F(111)-p6mup5*p6mup6*F(9569)+p5mup6*p6
     -   mup5*F(9571))))+8*(p5mup2*p6mup5*p6mup6*F(9569)+p2mup6*p5mup
     -   2*(p5mup5*F(110)-2*p2mup5*F(112)+p1mup5*F(115)-p6mup5*F(9572
     -   )))
       Fa3(69) = p4mup5*(p2mup2*(-8*(p4mup6*F(112)+p5mup6*F(113))+2*p6
     -   mup6*F(122))+mup2mup6*F(275))-2*p1mup2*(4*(p4mup5*p4mup6*F(1
     -   12)+p2mup5*p5mup6*F(113)+p4mup5*p5mup6*F(113))-p1mup5*p1mup6
     -   *F(114)+p1mup6*p2mup5*F(116)-p1mup6*p5mup5*F(117)+p1mup6*p6m
     -   up5*F(119)+p1mup5*p6mup6*F(120)+p5mup5*p6mup6*F(121)-p2mup5*
     -   p6mup6*F(122)-p4mup5*p6mup6*F(122)+2*(p5mup5*p5mup6*F(111)-p
     -   1mup5*p5mup6*F(118)-p6mup5*p6mup6*F(9569)+p5mup6*p6mup5*F(95
     -   71)-p2mup6*(p5mup5*F(110)-2*p2mup5*F(112)+p1mup5*F(115)-p6mu
     -   p5*F(9572))))
       Fa3(70) = 16*(p2mup6*p3mup2*p3mup5*F(112)+p2mup6*p3mup5*p4mup2*
     -   F(112)+p2mup6*p3mup2*p4mup5*F(112)+p2mup6*p4mup2*p4mup5*F(11
     -   2)+p2mup6*p3mup5*p5mup2*F(112)+p2mup6*p4mup5*p5mup2*F(112)+p
     -   3mup2*p3mup5*p5mup6*F(113)+p3mup5*p4mup2*p5mup6*F(113)+p3mup
     -   2*p4mup5*p5mup6*F(113)+p4mup2*p4mup5*p5mup6*F(113)+p3mup5*p5
     -   mup2*p5mup6*F(113)+p4mup5*p5mup2*p5mup6*F(113)+p2mup5*p5mup6
     -   *p6mup2*F(113))+4*(p1mup6*p3mup2*p3mup5*F(116)+p1mup6*p3mup5
     -   *p4mup2*F(116)+p1mup6*p3mup2*p4mup5*F(116)+p1mup6*p4mup2*p4m
     -   up5*F(116)+p1mup6*p3mup5*p5mup2*F(116)+p1mup6*p4mup5*p5mup2*
     -   F(116)+p1mup6*p2mup5*p6mup2*F(116)-p3mup2*p3mup5*p6mup6*F(12
     -   2)-p3mup5*p4mup2*p6mup6*F(122)-p3mup2*p4mup5*p6mup6*F(122)-p
     -   4mup2*p4mup5*p6mup6*F(122)-p3mup5*p5mup2*p6mup6*F(122)-p4mup
     -   5*p5mup2*p6mup6*F(122)-p2mup5*p6mup2*p6mup6*F(122)+p1mup2*(4
     -   *(p2mup6*(p3mup5+p4mup5)*F(112)+p3mup5*p5mup6*F(113))+p1mup6
     -   *p3mup5*F(116)+p1mup6*p4mup5*F(116)-p3mup5*p6mup6*F(122))+p2
     -   mup2*(4*(p2mup6*(p3mup5+p4mup5)*F(112)+p3mup5*p5mup6*F(113))
     -   +p1mup6*p3mup5*F(116)+p1mup6*p4mup5*F(116)-p3mup5*p6mup6*F(1
     -   22)))-2*mup2mup6*p3mup5*F(275)+mup2mup5*p1mup6*F(276)-mup2mu
     -   p5*p5mup6*F(277)-mup2mup5*p6mup6*F(278)-mup2mup5*p2mup6*F(27
     -   9)
       Fa3(71) = (p1mup2+p2mup2)*p4mup6*(p5mup5*F(110)-2*(p2mup5+p3mup
     -   5)*F(112)+p1mup5*F(115))
       Fa3(72) = p6mup2*(-2*(p2mup6*p5mup5*F(110)-p5mup5*p5mup6*F(111)
     -   )+4*p2mup6*(p2mup5+p3mup5+p4mup5)*F(112)+(p3mup5+p4mup5)*(4*
     -   p5mup6*F(113)+p1mup6*F(116)-p6mup6*F(122)))
       Return
       End
