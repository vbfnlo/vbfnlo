c************** Calling the Fa functions*************************
       subroutine NoAbe2gaPenLineFa2(p1mup2,p1mup4,p1mup5,p2mup2,p2mup
     &   4,p2mup5,p3mup2,p3mup4,p3mup5,p4mup2,p4mup4,p4mup5,p5mup2,p5
     &   mup4,p5mup5,mup2mup4,mup2mup5,mup4mup5,Fa2)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup4, p1mup5, p2mup2, p2mup4, p2mup5, 
     -          p3mup2, p3mup4, p3mup5, p4mup2, p4mup4, p4mup5, 
     -          p5mup2, p5mup4, p5mup5
       Complex*16   mup2mup4, mup2mup5, mup4mup5
        common/NoAbe2gaPenLineFhlFunctions/F
       COMMON/NoAbe2gaPenLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s
     &   23,s34,s45,s15
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12, s15 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45
      Complex*16 Fa(84),F(426)
      Real*8 P(767) 
        COMMON/NoAbe2gaPenLinePFunctions/P
       Complex*16 Fa2(29:56)
       COMMON/NoAbe2gaPenLineFaFunctions/Fa
       Fa2(29) = p3mup5*F(138)+p5mup5*F(143)+F(141)*Fa(16)
       Fa2(30) = p3mup5*p5mup4*F(136)+p5mup4*p5mup5*F(144)+mup4mup5*F(
     &   145)+(p3mup5*F(137)+p5mup5*F(142))*Fa(17)+Fa(16)*Fa(28)+p3mu
     &   p4*Fa2(29)
       Fa2(31) = p1mup4*p1mup5*F(147)+p3mup5*p5mup4*F(148)+p1mup5*p3mu
     &   p4*F(149)+p1mup4*p3mup5*F(150)+p2mup4*p3mup5*F(151)+p3mup4*p
     &   3mup5*F(152)+p1mup5*p2mup4*F(153)+p2mup5*p3mup4*F(154)+p2mup
     &   4*p2mup5*F(155)+p1mup4*p2mup5*F(156)+p5mup4*p5mup5*F(157)+p2
     &   mup4*p5mup5*F(158)+p2mup5*p5mup4*F(159)+p3mup4*p5mup5*F(160)
     &   +p1mup5*p5mup4*F(161)+p1mup4*p5mup5*F(162)+mup4mup5*F(163)
       Fa2(32) = p2mup2*p3mup5*F(164)+p3mup2*p3mup5*F(165)+p1mup2*p3mu
     &   p5*F(166)+p1mup2*p2mup5*F(167)+p1mup5*p2mup2*F(168)+p2mup2*p
     &   2mup5*F(169)+p1mup5*p3mup2*F(170)+p2mup5*p3mup2*F(171)+p1mup
     &   2*p1mup5*F(172)+p1mup5*p5mup2*F(173)+p2mup5*p5mup2*F(174)+p3
     &   mup5*p5mup2*F(175)+p3mup2*p5mup5*F(176)+p2mup2*p5mup5*F(177)
     &   +p1mup2*p5mup5*F(178)+p5mup2*p5mup5*F(179)+mup2mup5*F(180)
       Fa2(33) = p2mup2*p5mup4*F(181)+p1mup4*p2mup2*F(182)+p2mup2*p2mu
     &   p4*F(183)+p3mup2*p5mup4*F(184)+p1mup4*p3mup2*F(185)+p2mup4*p
     &   3mup2*F(186)+p1mup2*p2mup4*F(187)+p1mup2*p3mup4*F(188)+p2mup
     &   2*p3mup4*F(189)+p3mup2*p3mup4*F(190)+p1mup2*p1mup4*F(191)+p1
     &   mup2*p5mup4*F(192)+p2mup4*p5mup2*F(193)+p3mup4*p5mup2*F(194)
     &   +p1mup4*p5mup2*F(195)+p5mup2*p5mup4*F(196)+mup2mup4*F(197)
       Fa2(34) = p1mup2*F(201)+p5mup2*F(202)
       Fa2(35) = p1mup2*F(203)+p5mup2*F(204)
       Fa2(36) = p3mup5*F(206)+p5mup5*F(207)
       Fa2(37) = p5mup2*F(208)+p1mup2*F(209)
       Fa2(38) = p2mup5*F(210)+p5mup5*F(211)+p3mup5*F(212)
       Fa2(39) = p5mup2*F(213)+p1mup2*F(214)+p3mup2*F(215)
       Fa2(40) = p5mup4*F(216)+p3mup4*F(217)
       Fa2(41) = p3mup5*F(218)+p5mup5*F(219)
       Fa2(42) = p3mup5*F(221)+p5mup5*F(222)
       Fa2(43) = p3mup4*F(223)+p5mup4*F(224)
       Fa2(44) = p1mup2*F(225)+p5mup2*F(226)+p3mup2*F(227)
       Fa2(45) = p1mup2*F(228)+p5mup2*F(229)+p3mup2*F(230)
       Fa2(46) = p2mup5*F(231)+p1mup5*F(232)+p5mup5*F(233)+p3mup5*F(23
     &   4)
       Fa2(47) = p2mup4*F(235)+p5mup4*F(236)+p3mup4*F(237)
       Fa2(48) = p5mup2*F(238)+p1mup2*F(239)
       Fa2(49) = p5mup2*F(240)+p1mup2*F(241)
       Fa2(50) = p3mup4*F(245)+p5mup4*F(246)
       Fa2(51) = p3mup5*F(247)+p1mup5*F(248)+p5mup5*F(249)
       Fa2(52) = p1mup2*F(252)+p5mup2*F(253)
       Fa2(53) = p1mup2*F(256)+p5mup2*F(257)
       Fa2(54) = p2mup2*p5mup4*p5mup5*F(198)+p2mup5*p5mup2*p5mup4*F(19
     &   9)+p1mup2*p2mup5*p5mup4*F(200)+p3mup2*p3mup5*p5mup4*F(205)+p
     &   2mup2*p3mup5*p5mup4*F(220)+p3mup5*p5mup2*p5mup4*F(242)+p3mup
     &   2*p5mup4*p5mup5*F(243)+p1mup2*p1mup4*p3mup5*F(244)+p1mup4*p3
     &   mup5*p5mup2*F(250)+p1mup2*p3mup5*p5mup4*F(251)+p1mup2*p2mup4
     &   *p5mup5*F(254)+p1mup2*p3mup4*p5mup5*F(255)+p1mup2*p5mup4*p5m
     &   up5*F(258)+p2mup4*p5mup2*p5mup5*F(259)+p3mup4*p5mup2*p5mup5*
     &   F(260)+p5mup2*p5mup4*p5mup5*F(261)+mup2mup4*p3mup5*F(262)+mu
     &   p2mup4*p2mup5*F(263)+mup2mup4*p1mup5*F(264)+mup2mup4*p5mup5*
     &   F(265)+mup2mup5*p3mup4*F(266)+mup2mup5*p2mup4*F(267)+mup2mup
     &   5*p1mup4*F(268)+mup2mup5*p5mup4*F(269)+mup4mup5*p1mup2*F(270
     &   )+mup4mup5*p2mup2*F(271)+mup4mup5*p3mup2*F(272)+mup4mup5*p5m
     &   up2*F(273)+p1mup4*p1mup5*Fa2(34)+p1mup5*p5mup4*Fa2(35)+p3mup
     &   2*p3mup4*Fa2(36)+p3mup4*p3mup5*Fa2(37)+p2mup2*p2mup4*Fa2(38)
     &   +p2mup4*p2mup5*Fa2(39)+p2mup2*p2mup5*Fa2(40)+p2mup4*p3mup2*F
     &   a2(41)+p2mup2*p3mup4*Fa2(42)+p2mup5*p3mup2*Fa2(43)+p1mup4*p2
     &   mup5*Fa2(44)+p1mup5*p2mup4*Fa2(45)+p1mup4*p2mup2*Fa2(46)+p1m
     &   up5*p2mup2*Fa2(47)+p2mup5*p3mup4*Fa2(48)+p2mup4*p3mup5*Fa2(4
     &   9)+p1mup5*p3mup2*Fa2(50)+p1mup4*p3mup2*Fa2(51)+p1mup5*p3mup4
     &   *Fa2(52)+p1mup4*p5mup5*Fa2(53)
       Fa2(55) = p3mup5*F(275)+p5mup5*F(276)
       Fa2(56) = p5mup2*F(278)+p1mup2*F(279)
       Return
       End
