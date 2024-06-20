c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim1(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,
     -   p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mu
     -   p3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup
     -   3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa1)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenOddGSim/F
      Complex*16 EpsF(126)
       common/EPSFPenOddGSim/EPSF
       COMMON/InvariantsPenOddGSim/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s3
     -   4,s45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       Real*8   p1sq 
        Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s234,m1,m1sq,m3,m3sq,m4,m4sq
      Complex*16 Fa(3196),F(4476)
      Real*8 P(1383) 
        COMMON/PPenOddGSim/P
       Complex*16 Fa1(1:99)
       COMMON/FaPenOddGSim/Fa
       Fa1(1) = p1mup5*F(1)+p3mup5*F(2)+p2mup5*F(3)+p4mup5*F(4)
       Fa1(2) = mup4mup5*(2*(p1mup3*F(5)+p2mup3*F(6))+p3mup3*F(7))+2*(
     -   mup3mup5*p1mup4*F(8)+mup3mup5*p2mup4*F(9)+mup3mup4*p1mup5*F(
     -   10)+mup3mup4*p2mup5*F(11))-4*(p1mup4*p3mup3*p4mup5*F(23)-p1m
     -   up4*p3mup5*p4mup3*F(24)+p2mup3*p3mup4*p4mup5*F(40)+p2mup3*p3
     -   mup5*p4mup4*F(41)+p2mup4*p3mup3*p4mup5*F(42)-p2mup4*p3mup5*p
     -   4mup3*F(43)+p2mup5*p3mup3*p4mup4*F(44)+p2mup5*p3mup4*p4mup3*
     -   F(45)+p3mup3*p3mup4*p4mup5*F(49)+p3mup3*p3mup5*p4mup4*F(50)+
     -   p1mup3*(4*p1mup4*p1mup5*F(12)+p3mup4*(2*p3mup5*F(14)+p4mup5*
     -   F(21)+4*p2mup5*F(28))+2*p2mup4*(p3mup5*F(16)+2*p2mup5*F(27)+
     -   p4mup5*F(31))+p4mup4*(p3mup5*F(22)+2*(p2mup5*F(32)+p4mup5*F(
     -   51))))+p3mup4*p3mup5*p4mup3*F(52)+p3mup3*p4mup4*p4mup5*F(54)
     -   +p3mup4*p4mup3*p4mup5*F(55)+p3mup5*p4mup3*p4mup4*F(56)+p1mup
     -   5*(2*p2mup3*(2*(p2mup4*F(13)+p3mup4*F(17))+p4mup4*F(20))+p3m
     -   up3*(2*(p3mup4*F(15)+p2mup4*F(18))+p4mup4*F(25))+p4mup3*(p3m
     -   up4*F(26)+2*p4mup4*F(60)+4*p2mup4*F(61)))+4*(p1mup4*p2mup3*p
     -   2mup5*F(27)+p2mup3*p2mup5*p3mup4*F(34)+p1mup4*p2mup5*p4mup3*
     -   F(62)+p2mup4*p2mup5*p4mup3*F(63))+2*(p1mup4*p2mup3*p4mup5*F(
     -   19)+p1mup4*p2mup3*p3mup5*F(29)+p1mup4*p2mup5*p3mup3*F(30)+p2
     -   mup3*p2mup4*p3mup5*F(33)+p2mup4*p2mup5*p3mup3*F(35)+p2mup3*p
     -   2mup4*p4mup5*F(36)+p2mup3*p2mup5*p4mup4*F(37)+p2mup3*p3mup4*
     -   p3mup5*F(38)+p2mup5*p3mup3*p3mup4*F(39)+p3mup3*p3mup4*p3mup5
     -   *F(46)+p1mup4*p3mup3*p3mup5*F(47)+p2mup4*p3mup3*p3mup5*F(48)
     -   +p2mup3*p4mup4*p4mup5*F(53)+p4mup3*p4mup4*p4mup5*F(57)+p1mup
     -   4*p4mup3*p4mup5*F(58)+p2mup4*p4mup3*p4mup5*F(59)+p2mup5*p4mu
     -   p3*p4mup4*F(64)))
       Fa1(3) = mup3mup4*p1mup5*F(67)+mup3mup5*(p1mup4*F(66)+p2mup4*F(
     -   68))+2*(mup4mup5*p2mup3*F(69)+mup3mup4*p2mup5*F(70)-mup4mup5
     -   *p3mup3*F(71))+p1mup3*(mup4mup5*F(65)+8*p1mup4*p1mup5*F(79)+
     -   4*(p2mup5*p4mup4*F(44)+2*(p2mup5*p3mup4*F(39)+p3mup4*p3mup5*
     -   F(80))+p2mup4*(p4mup5*F(74)+2*(p2mup5*F(35)+p3mup5*F(81)))+p
     -   3mup4*p4mup5*F(82)+p3mup5*p4mup4*F(83)+p4mup4*p4mup5*F(86)))
     -   +4*(p1mup5*(p2mup3*(-2*(p2mup4*F(72)-p3mup4*F(73))+p4mup4*F(
     -   77))+p2mup4*p4mup3*F(78)+p3mup4*p4mup3*F(85)+p4mup3*p4mup4*F
     -   (88)+2*(p3mup3*p3mup4*F(101)+p2mup4*p3mup3*F(103))+p3mup3*p4
     -   mup4*F(104))+p1mup4*(-(p2mup3*p4mup5*F(75))-p4mup3*(p2mup5*F
     -   (76)+p3mup5*F(84))+p4mup3*p4mup5*F(87)+p3mup3*p4mup5*F(89)+4
     -   *(p3mup3*p3mup5*F(100)+p2mup3*p3mup5*F(102)+p2mup3*p2mup5*F(
     -   109)+p2mup5*p3mup3*F(110)))+2*(p2mup3*(p2mup4*p4mup5*F(90)+p
     -   4mup4*p4mup5*F(93)+p3mup4*p4mup5*F(96)+p3mup5*p4mup4*F(108)+
     -   2*(p2mup4*p3mup5*F(105)+p3mup4*p3mup5*F(106)+p2mup5*p3mup4*F
     -   (111))+p2mup5*p4mup4*F(112))+p3mup3*(p3mup4*p4mup5*F(91)+p4m
     -   up4*p4mup5*F(94)+p2mup4*p4mup5*F(97)+p3mup5*p4mup4*F(107)+2*
     -   (p3mup4*p3mup5*F(99)+p2mup4*p3mup5*F(106)+p2mup4*p2mup5*F(11
     -   1)+p2mup5*p3mup4*F(113))+p2mup5*p4mup4*F(114))+p4mup3*(p3mup
     -   4*p3mup5*F(92)+p3mup4*p4mup5*F(95)+p2mup4*p3mup5*F(98)+p2mup
     -   4*p2mup5*F(115)+p2mup4*p4mup5*F(116)+p2mup5*p3mup4*F(117)+p3
     -   mup5*p4mup4*F(118)+p2mup5*p4mup4*F(119)+2*p4mup4*p4mup5*F(12
     -   0))))
       Fa1(4) = mup3mup4*p1mup5*F(123)+mup4mup5*p2mup3*F(124)+mup4mup5
     -   *p3mup3*F(125)+2*(mup3mup5*p2mup4*F(126)-mup3mup4*p2mup5*F(1
     -   27))+p1mup4*(mup3mup5*F(122)+8*p1mup3*p1mup5*F(136)+4*(p2mup
     -   5*p3mup3*F(44)+p3mup3*p3mup5*F(134)+p2mup3*(p3mup5*F(131)+2*
     -   (p2mup5*F(37)+p4mup5*F(138)))+p3mup3*p4mup5*F(140)+p3mup5*p4
     -   mup3*F(143)+2*(p2mup5*p4mup3*F(64)+p4mup3*p4mup5*F(152))))-4
     -   *(p1mup5*p2mup3*p3mup4*F(132)-p1mup5*p2mup4*p3mup3*F(133)-p1
     -   mup5*p3mup3*p3mup4*F(135)-p1mup5*p3mup4*p4mup3*F(144)-p1mup5
     -   *p3mup3*p4mup4*F(156))+p1mup3*(mup4mup5*F(121)+4*(-(p2mup4*p
     -   3mup5*F(129))+p2mup5*p3mup4*F(130)+p3mup4*p3mup5*F(137)+p3mu
     -   p4*p4mup5*F(139)+p3mup5*p4mup4*F(142)+4*(p2mup4*p2mup5*F(63)
     -   +p4mup4*p4mup5*F(151)+p2mup4*p4mup5*F(153)+p2mup5*p4mup4*F(1
     -   66))))+8*(p3mup3*p3mup5*p4mup4*F(95)+p2mup4*p2mup5*p3mup3*F(
     -   115)-p1mup5*p2mup3*p2mup4*F(128)+p2mup3*p2mup4*p3mup5*F(145)
     -   +p2mup3*p3mup4*p3mup5*F(146)+p2mup4*p3mup3*p3mup5*F(147)+p3m
     -   up3*p3mup4*p4mup5*F(148)+p2mup3*p3mup4*p4mup5*F(149)+p2mup3*
     -   p3mup5*p4mup4*F(150)+p1mup5*p2mup3*p4mup4*F(154)+p1mup5*p2mu
     -   p4*p4mup3*F(155)+p3mup4*p3mup5*p4mup3*F(158)+p3mup3*p4mup4*p
     -   4mup5*F(160)+p3mup4*p4mup3*p4mup5*F(161)+p2mup4*p3mup3*p4mup
     -   5*F(162)+p2mup4*p3mup5*p4mup3*F(163)+p2mup3*p2mup5*p3mup4*F(
     -   164)+p2mup5*p3mup3*p3mup4*F(165)+p2mup5*p3mup3*p4mup4*F(168)
     -   +p2mup5*p3mup4*p4mup3*F(169)+p1mup5*p4mup3*p4mup4*F(170)+p3m
     -   up5*p4mup3*p4mup4*F(173))+16*(p3mup3*p3mup4*p3mup5*F(141)+p2
     -   mup3*p2mup4*p4mup5*F(157)+p2mup3*p4mup4*p4mup5*F(159)+p2mup3
     -   *p2mup5*p4mup4*F(167)+p2mup4*p2mup5*p4mup3*F(171)+p2mup5*p4m
     -   up3*p4mup4*F(172)+p4mup3*p4mup4*p4mup5*F(174)+p2mup4*p4mup3*
     -   p4mup5*F(175))
       Fa1(5) = (mup4mup5*F(176))/2.+p2mup5*p3mup4*F(178)+p2mup5*p4mup
     -   4*F(180)+p3mup4*p4mup5*F(183)+p3mup5*p4mup4*F(184)+p3mup4*p3
     -   mup5*F(186)+p2mup4*(p3mup5*F(177)+p4mup5*F(179)+p2mup5*F(182
     -   )+2*p1mup5*F(187))+4*p1mup4*p1mup5*F(188)+p1mup5*p3mup4*F(18
     -   9)+2*(p1mup4*p2mup5*F(181)-p1mup4*p3mup5*F(185)+p1mup4*p4mup
     -   5*F(190))+p1mup5*p4mup4*F(191)+p4mup4*p4mup5*F(192)
       Fa1(6) = (mup4mup5*F(193))/2.+p1mup5*p3mup4*F(195)+p1mup5*p4mup
     -   4*F(199)+p3mup4*p4mup5*F(200)+p3mup5*p4mup4*F(201)+p2mup4*p3
     -   mup5*F(202)+p3mup4*p3mup5*F(203)+4*p2mup4*p2mup5*F(207)+p1mu
     -   p4*(p3mup5*F(194)+p4mup5*F(198)-2*(p1mup5*F(196)+p2mup5*F(20
     -   8)))-2*(p1mup5*p2mup4*F(197)-p2mup5*p3mup4*F(204)-p2mup4*p4m
     -   up5*F(205)-p2mup5*p4mup4*F(206)-p4mup4*p4mup5*F(209))
       Fa1(7) = (mup4mup5*F(210))/2.+p1mup5*p2mup4*F(212)+p1mup5*p4mup
     -   4*F(216)+p2mup4*p4mup5*F(217)+p2mup5*p4mup4*F(218)+p2mup4*p2
     -   mup5*F(219)+p2mup5*p3mup4*F(220)+p1mup4*(p2mup5*F(211)+p4mup
     -   5*F(215)+2*(p1mup5*F(213)+p3mup5*F(223)))-4*p3mup4*p3mup5*F(
     -   225)+2*(p1mup5*p3mup4*F(214)+p3mup4*p4mup5*F(221)+p3mup5*p4m
     -   up4*F(222)+p2mup4*p3mup5*F(224)+p4mup4*p4mup5*F(226))
       Fa1(8) = (mup4mup5*F(227))/2.+p1mup5*p2mup4*F(229)+p1mup5*p3mup
     -   4*F(231)+p2mup4*p3mup5*F(232)+p2mup5*p3mup4*F(233)+p2mup4*p2
     -   mup5*F(236)+2*(p2mup4*p4mup5*F(234)+p2mup5*p4mup4*F(235)+p3m
     -   up4*p4mup5*F(238)-p3mup5*p4mup4*F(239))+p1mup4*(p2mup5*F(228
     -   )+p3mup5*F(230)+p1mup5*F(237)+p4mup5*F(240))+p3mup4*p3mup5*F
     -   (241)+p1mup5*p4mup4*F(242)+p4mup4*p4mup5*F(243)
       Fa1(9) = p1mup4*F(244)+p3mup4*F(245)+p2mup4*F(246)+p4mup4*F(247
     -   )
       Fa1(10) = mup4mup5*p2mup3*F(250)+mup3mup5*p2mup4*F(251)+mup3mup
     -   4*p2mup5*F(252)+p2mup5*p4mup3*p4mup4*F(278)+8*(p1mup3*p2mup4
     -   *p2mup5*F(254)-p1mup5*p2mup3*p2mup4*F(255)+p1mup3*p2mup5*p3m
     -   up4*F(256)-p1mup5*p2mup3*p3mup4*F(257)-p1mup4*p2mup3*p2mup5*
     -   F(265)-p2mup3*p2mup4*p3mup5*F(265)-p1mup4*p2mup3*p3mup5*F(26
     -   6)-p2mup3*p2mup5*p3mup4*F(267)-p2mup3*p3mup4*p3mup5*F(270)-p
     -   1mup4*p2mup5*p3mup3*F(279)-p2mup4*p2mup5*p3mup3*F(280)-p2mup
     -   5*p3mup3*p3mup4*F(281))+4*(p1mup3*p2mup4*p4mup5*F(258)+p1mup
     -   3*p2mup5*p4mup4*F(259)-p1mup4*p2mup3*p4mup5*F(260)-p1mup4*p2
     -   mup5*p4mup3*F(261)-p1mup5*p2mup3*p4mup4*F(262)+p1mup5*p2mup4
     -   *p4mup3*F(263)-p2mup3*p2mup4*p4mup5*F(268)-p2mup3*p2mup5*p4m
     -   up4*F(269)-p2mup3*p3mup4*p4mup5*F(271)-p2mup3*p3mup5*p4mup4*
     -   F(272)+p2mup4*p3mup5*p4mup3*F(273)-p2mup5*p3mup4*p4mup3*F(27
     -   4)+p2mup4*p2mup5*p4mup3*F(275)-p2mup3*p4mup4*p4mup5*F(276)+p
     -   2mup4*p4mup3*p4mup5*F(277)-p2mup4*p3mup3*p4mup5*F(282))+p2mu
     -   p5*p3mup3*p4mup4*F(283)-2*(mup3mup5*p1mup4*F(253)-(-(mup3mup
     -   4*p1mup5*F(263))+mup4mup5*p1mup3*F(264)+mup4mup5*p3mup3*F(28
     -   2))*P(7))
       Fa1(11) = mup3mup4*p1mup5*F(286)+mup3mup5*p2mup4*F(288)+mup3mup
     -   4*p2mup5*F(289)+mup4mup5*(p1mup3*F(284)+p2mup3*F(287)-2*p3mu
     -   p3*F(290))+p1mup3*p2mup5*p4mup4*F(300)+p1mup4*(mup3mup5*F(28
     -   5)+4*(p2mup3*(p3mup5*F(272)+p2mup5*F(291)+p4mup5*F(298))+p2m
     -   up5*p3mup3*F(299)+p2mup5*p4mup3*F(301)))-4*(p1mup5*p2mup3*p2
     -   mup4*F(292)+p1mup3*p2mup4*p3mup5*F(293)-p1mup3*p2mup5*p3mup4
     -   *F(294)-p1mup5*p2mup3*p3mup4*F(295)+p1mup3*p2mup4*p2mup5*F(2
     -   96)+p1mup3*p2mup4*p4mup5*F(297)+p1mup5*p2mup4*p4mup3*F(302)+
     -   p1mup5*p2mup4*p3mup3*F(310)-p1mup5*p2mup3*p4mup4*F(313))+8*(
     -   p2mup3*p2mup4*p3mup5*F(303)-p2mup4*p2mup5*p3mup3*F(304)+p2mu
     -   p3*p2mup5*p4mup4*F(304)-p2mup3*p2mup4*p4mup5*F(305)+p2mup3*p
     -   3mup4*p3mup5*F(306)-p2mup4*p3mup5*p4mup3*F(306)-p2mup4*p3mup
     -   3*p4mup5*F(306)-p2mup5*p3mup3*p3mup4*F(307)+p2mup3*p3mup5*p4
     -   mup4*F(308)+p2mup5*p3mup4*p4mup3*F(309)+p2mup3*p2mup5*p3mup4
     -   *F(311)-p2mup4*p3mup3*p3mup5*F(312)-p2mup4*p2mup5*p4mup3*F(3
     -   14)+p2mup3*p4mup4*p4mup5*F(315)+p2mup3*p3mup4*p4mup5*F(316)-
     -   p2mup4*p4mup3*p4mup5*F(316)+p2mup5*p3mup3*p4mup4*F(317)+p2mu
     -   p5*p4mup3*p4mup4*F(318))
       Fa1(12) = p1mup2*F(319)+p3mup2*F(320)+p4mup2*F(321)
       Fa1(13) = p1mup3*F(324)+p2mup3*F(325)
       Fa1(14) = mup4mup5*p2mup3*F(336)+mup3mup4*p2mup5*F(338)+mup4mup
     -   5*p1mup3*F(339)+mup4mup5*p3mup3*F(340)+mup3mup5*p1mup4*F(341
     -   )-p2mup4*(-(mup3mup5*F(337))+4*(p3mup5*(p4mup3*F(302)+p1mup3
     -   *F(344)+p2mup3*F(349))+p3mup3*(p4mup5*F(297)+p3mup5*F(351)+p
     -   2mup5*F(355)))+8*(p2mup5*p4mup3*F(314)+p1mup3*p2mup5*F(353)+
     -   p2mup3*p4mup5*F(353)+p1mup3*p4mup5*F(354)+p4mup3*p4mup5*F(36
     -   1)))+4*(p1mup3*p2mup5*p3mup4*F(345)+p1mup4*p2mup3*p3mup5*F(3
     -   46)+p1mup4*p2mup5*p3mup3*F(347)+p2mup3*p2mup5*p3mup4*F(350)+
     -   p2mup5*p3mup3*p3mup4*F(352)-p2mup3*p3mup4*p3mup5*F(356)-p2mu
     -   p3*p3mup4*p4mup5*F(357)+p2mup5*p3mup4*p4mup3*F(358)-p2mup3*p
     -   3mup5*p4mup4*F(362))+p2mup5*p3mup3*p4mup4*F(363)-p1mup5*(-(m
     -   up3mup4*F(342))+4*(p2mup3*p3mup4*F(263)+p2mup4*p3mup3*F(348)
     -   )+8*(p2mup3*p2mup4*F(255)+p2mup4*p4mup3*F(364)))-8*(-(p1mup4
     -   *p2mup3*p2mup5*F(343))+p1mup4*p2mup5*p4mup3*F(359)+p1mup3*p2
     -   mup5*p4mup4*F(359)+p2mup3*p2mup5*p4mup4*F(360)-p2mup5*p4mup3
     -   *p4mup4*F(365))
       Fa1(15) = p2mup1*F(334)+p3mup1*F(335)+p4mup1*F(370)
       Fa1(16) = p2mup5*F(372)+p1mup5*F(373)+p3mup5*F(374)+p4mup5*F(37
     -   5)
       Fa1(17) = mup4mup5*(2*(p1mup3*F(376)+p2mup3*F(377))+p3mup3*F(37
     -   8))-2*((mup3mup5*p1mup4-mup3mup4*p1mup5)*F(263)+(mup3mup5*p2
     -   mup4-mup3mup4*p2mup5)*F(379))
       Fa1(18) = -((mup3mup5*p1mup4-mup3mup4*p1mup5)*F(380))-2*(mup3mu
     -   p5*p2mup4-mup3mup4*p2mup5)*F(381)+mup4mup5*(p1mup3*F(380)+2*
     -   (p2mup3*F(381)+p3mup3*F(382)))
       Fa1(19) = mup4mup5*p1mup3*F(383)-mup3mup5*p1mup4*F(383)+mup3mup
     -   4*p1mup5*F(383)+2*(mup4mup5*p3mup3*F(4465)+(mup4mup5*p2mup3-
     -   mup3mup5*p2mup4+mup3mup4*p2mup5)*F(4469))
       Fa1(20) = -4*(p1mup5*p2mup4*F(263)-p1mup4*p2mup5*F(263))+mup4mu
     -   p5*F(384)+2*(p2mup5*p4mup4*F(383)-p2mup4*p4mup5*F(383)+p2mup
     -   4*p3mup5*F(385)+p2mup5*p3mup4*F(386))
       Fa1(21) = mup4mup5*F(387)-2*(p1mup5*(2*p2mup4*F(379)-p3mup4*F(3
     -   88)+p4mup4*F(4464))-p1mup4*(2*p2mup5*F(379)-p3mup5*F(388)+p4
     -   mup5*F(4464)))+4*(p2mup5*p3mup4*F(381)-p2mup4*p3mup5*F(381)-
     -   p3mup5*p4mup4*F(4464)+p3mup4*p4mup5*F(4464)+(p2mup5*p4mup4-p
     -   2mup4*p4mup5)*F(4469))
       Fa1(22) = 2*(p2mup5*(p1mup4*F(380)+2*(p4mup4*F(356)+p3mup4*F(38
     -   2)))-p2mup4*(p1mup5*F(380)+2*(p4mup5*F(356)+p3mup5*F(382))))
     -   +mup4mup5*F(389)
       Fa1(23) = -2*(p1mup5*p2mup4*F(383)-p1mup4*p2mup5*F(383))+mup4mu
     -   p5*F(390)+4*((p2mup5*p3mup4-p2mup4*p3mup5)*F(4465)+(p2mup5*p
     -   4mup4-p2mup4*p4mup5)*F(4467))
       Fa1(24) = p2mup4*F(391)+p1mup4*F(392)+p3mup4*F(393)+p4mup4*F(39
     -   4)
       Fa1(25) = mup4mup5*p2mup3-mup3mup5*p2mup4+mup3mup4*p2mup5
       Fa1(26) = p1mup2*F(397)+p3mup2*F(398)-p4mup2*F(399)
       Fa1(27) = -(p2mup3*F(402))+p1mup3*F(403)
       Fa1(28) = p3mup1*F(389)+p4mup1*F(390)+p2mup1*F(395)
       Fa1(29) = p1mup5*F(414)+p2mup5*F(415)+p3mup5*F(416)+p4mup5*F(41
     -   7)
       Fa1(30) = 2*((mup4mup5*p1mup3-mup3mup5*p1mup4+mup3mup4*p1mup5)*
     -   F(263)+(mup4mup5*p2mup3-mup3mup5*p2mup4+mup3mup4*p2mup5)*F(3
     -   79))+mup4mup5*p3mup3*F(418)
       Fa1(31) = mup4mup5*p1mup3*F(418)-mup3mup5*p1mup4*F(418)+mup3mup
     -   4*p1mup5*F(418)+2*(mup4mup5*p3mup3*F(4466)+(mup4mup5*p2mup3-
     -   mup3mup5*p2mup4+mup3mup4*p2mup5)*F(4468))
       Fa1(32) = -4*(p1mup5*p2mup4*F(263)-p1mup4*p2mup5*F(263))+2*((p2
     -   mup5*p4mup4-p2mup4*p4mup5)*F(383)+(p2mup5*p3mup4-p2mup4*p3mu
     -   p5)*F(418))+mup4mup5*F(419)
       Fa1(33) = mup4mup5*F(420)-2*(p1mup5*(2*p2mup4*F(379)+p3mup4*F(4
     -   463)+p4mup4*F(4464))-p1mup4*(2*p2mup5*F(379)+p3mup5*F(4463)+
     -   p4mup5*F(4464)))+4*((p2mup5*p3mup4-p2mup4*p3mup5)*F(4468)+(p
     -   2mup5*p4mup4-p2mup4*p4mup5)*F(4469))
       Fa1(34) = -2*(p1mup5*p2mup4*F(418)-p1mup4*p2mup5*F(418))+mup4mu
     -   p5*F(421)+4*((p2mup5*p4mup4-p2mup4*p4mup5)*F(4465)+(p2mup5*p
     -   3mup4-p2mup4*p3mup5)*F(4466))
       Fa1(35) = p1mup4*F(422)+p2mup4*F(423)+p3mup4*F(424)+p4mup4*F(42
     -   5)
       Fa1(36) = p1mup2*F(426)+p3mup2*F(427)+p4mup2*F(428)
       Fa1(37) = p1mup3*F(414)+p2mup3*F(415)
       Fa1(38) = p4mup1*F(390)+p3mup1*F(421)+p2mup1*F(435)
       Fa1(39) = p2mup5*F(439)+p1mup5*F(440)+p3mup5*F(441)+p4mup5*F(44
     -   2)
       Fa1(40) = mup4mup5*(2*(p1mup3*F(376)+p2mup3*F(377))+p3mup3*F(37
     -   8))-2*(mup3mup5*(p1mup4*F(376)+p2mup4*F(377))-mup3mup4*(p1mu
     -   p5*F(263)+p2mup5*F(379)))
       Fa1(41) = -((mup3mup5*p1mup4-mup3mup4*p1mup5)*F(443))-2*(mup3mu
     -   p5*p2mup4-mup3mup4*p2mup5)*F(444)+mup4mup5*(p1mup3*F(443)+2*
     -   (p2mup3*F(444)+p3mup3*F(445)))
       Fa1(42) = 4*p1mup5*p2mup4*F(263)+mup4mup5*F(446)+2*(p2mup4*(p3m
     -   up5*F(386)+p4mup5*F(447))-p2mup5*(p3mup4*F(378)+p4mup4*F(448
     -   )+2*(p1mup4*F(376)+p2mup4*F(449))))
       Fa1(43) = mup4mup5*F(450)+2*(p1mup5*p3mup4*F(452)+p1mup5*p4mup4
     -   *F(454)+p1mup4*(2*(p2mup5*F(379)+p1mup5*F(451))+p3mup5*F(452
     -   )+p4mup5*F(454)))+4*(p1mup5*p2mup4*F(453)+p2mup4*p4mup5*F(45
     -   5)+p2mup5*(p3mup4*F(381)+p4mup4*F(444)+2*p2mup4*F(456))+p2mu
     -   p4*p3mup5*F(457)+p3mup5*p4mup4*F(458)+p3mup4*p4mup5*F(458))+
     -   8*(p3mup4*p3mup5*F(382)+p4mup4*p4mup5*F(459))
       Fa1(44) = -2*(p2mup5*(p1mup4*F(380)+2*(p4mup4*F(356)+p3mup4*F(3
     -   82)))-p2mup4*(p1mup5*F(380)+2*(p4mup5*F(356)+p3mup5*F(382)))
     -   )+mup4mup5*F(412)
       Fa1(45) = 2*(p1mup5*p2mup4*F(443)-p1mup4*p2mup5*F(443))-4*((p2m
     -   up5*p3mup4-p2mup4*p3mup5)*F(445)+(p2mup5*p4mup4-p2mup4*p4mup
     -   5)*F(459))+mup4mup5*F(460)
       Fa1(46) = -(p2mup4*F(461))+p1mup4*F(462)+p3mup4*F(463)+p4mup4*F
     -   (464)
       Fa1(47) = p1mup2*F(466)+p3mup2*F(467)+p4mup2*F(468)
       Fa1(48) = p2mup3*F(470)+p1mup3*F(471)
       Fa1(49) = p3mup1*F(412)+p2mup1*F(450)+p4mup1*F(460)
       Fa1(50) = p2mup5*F(481)+p1mup5*F(482)+p3mup5*F(483)+p4mup5*F(48
     -   4)
       Fa1(51) = mup4mup5*F(478)+2*(p2mup5*(2*p1mup4*F(376)+p3mup4*F(3
     -   78)+p4mup4*F(448))+p2mup4*(p3mup5*F(385)-2*(p1mup5*F(263)-p2
     -   mup5*F(449))+p4mup5*F(485)))
       Fa1(52) = 2*(p1mup5*p3mup4*F(452)+p1mup5*p4mup4*F(454)+p1mup4*(
     -   2*(p2mup5*F(379)+p1mup5*F(451))+p3mup5*F(452)+p4mup5*F(454))
     -   )+4*(p1mup5*p2mup4*F(453)+p2mup4*p4mup5*F(455)+p2mup5*(p3mup
     -   4*F(381)+p4mup4*F(444)+2*p2mup4*F(456))+p2mup4*p3mup5*F(457)
     -   +p3mup5*p4mup4*F(458)+p3mup4*p4mup5*F(458))+8*(p3mup4*p3mup5
     -   *F(382)+p4mup4*p4mup5*F(459))-mup4mup5*F(479)
       Fa1(53) = -2*(p1mup5*p2mup4*F(443)-p1mup4*p2mup5*F(443))+4*((p2
     -   mup5*p3mup4-p2mup4*p3mup5)*F(445)+(p2mup5*p4mup4-p2mup4*p4mu
     -   p5)*F(459))+mup4mup5*F(480)
       Fa1(54) = p2mup4*F(470)+p1mup4*F(471)+p3mup4*F(486)+p4mup4*F(48
     -   7)
       Fa1(55) = p1mup2*F(489)+p3mup2*F(490)+p4mup2*F(491)
       Fa1(56) = -(p2mup3*F(461))+p1mup3*F(462)
       Fa1(57) = p3mup1*F(389)+p2mup1*F(479)+p4mup1*F(480)
       Fa1(58) = p1mup5*F(422)+p2mup5*F(423)+p3mup5*F(424)+p4mup5*F(42
     -   5)
       Fa1(59) = 4*p1mup5*p2mup4*F(263)-2*(2*p1mup4*p2mup5*F(263)+(p2m
     -   up5*p4mup4-p2mup4*p4mup5)*F(383)+(p2mup5*p3mup4-p2mup4*p3mup
     -   5)*F(418))+mup4mup5*F(436)
       Fa1(60) = mup4mup5*F(492)+2*(p1mup5*(2*p2mup4*F(379)+p3mup4*F(4
     -   463)+p4mup4*F(4464))-p1mup4*(2*p2mup5*F(379)+p3mup5*F(4463)+
     -   p4mup5*F(4464)))+4*((-(p2mup5*p3mup4)+p2mup4*p3mup5)*F(4468)
     -   +(-(p2mup5*p4mup4)+p2mup4*p4mup5)*F(4469))
       Fa1(61) = 2*(p1mup5*p2mup4*F(418)-p1mup4*p2mup5*F(418))+mup4mup
     -   5*F(438)-4*((p2mup5*p4mup4-p2mup4*p4mup5)*F(4465)+(p2mup5*p3
     -   mup4-p2mup4*p3mup5)*F(4466))
       Fa1(62) = 2*(p1mup5*p2mup4*F(383)-p1mup4*p2mup5*F(383))+mup4mup
     -   5*F(413)-4*((p2mup5*p3mup4-p2mup4*p3mup5)*F(4465)+(p2mup5*p4
     -   mup4-p2mup4*p4mup5)*F(4467))
       Fa1(63) = p1mup4*F(414)+p2mup4*F(415)+p3mup4*F(416)+p4mup4*F(41
     -   7)
       Fa1(64) = p1mup2*F(493)+p3mup2*F(494)+p4mup2*F(495)
       Fa1(65) = p1mup3*F(422)+p2mup3*F(423)
       Fa1(66) = p4mup1*F(413)+p2mup1*F(437)+p3mup1*F(438)
       Fa1(67) = p2mup5*F(391)+p1mup5*F(392)+p3mup5*F(393)+p4mup5*F(39
     -   4)
       Fa1(68) = 4*(p1mup5*p2mup4*F(263)-p1mup4*p2mup5*F(263))-2*(p2mu
     -   p5*p4mup4*F(383)-p2mup4*p4mup5*F(383)-p2mup5*p3mup4*F(385)-p
     -   2mup4*p3mup5*F(386))+mup4mup5*F(410)
       Fa1(69) = mup4mup5*F(496)+2*(p1mup5*(2*p2mup4*F(379)-p3mup4*F(3
     -   88)+p4mup4*F(4464))-p1mup4*(2*p2mup5*F(379)-p3mup5*F(388)+p4
     -   mup5*F(4464)))+4*(-(p2mup5*p3mup4*F(381))+p2mup4*p3mup5*F(38
     -   1)+p3mup5*p4mup4*F(4464)-p3mup4*p4mup5*F(4464)+(-(p2mup5*p4m
     -   up4)+p2mup4*p4mup5)*F(4469))
       Fa1(70) = p2mup4*F(372)+p1mup4*F(373)+p3mup4*F(374)+p4mup4*F(37
     -   5)
       Fa1(71) = p4mup2*F(399)+p1mup2*F(497)+p3mup2*F(498)
       Fa1(72) = p2mup3*F(500)+p1mup3*F(501)
       Fa1(73) = p2mup1*F(411)+p3mup1*F(412)+p4mup1*F(413)
       Fa1(74) = p1mup5*F(505)+p3mup5*F(506)+p2mup5*F(507)+p4mup5*F(50
     -   8)
       Fa1(75) = mup4mup5*(-2*(p1mup3*F(5)+p2mup3*F(6))+p3mup3*F(509))
     -   +2*(8*(p1mup4*p2mup3*p2mup5*F(27)+p2mup3*p2mup5*p3mup4*F(34)
     -   +p1mup4*p2mup5*p4mup3*F(62)+p2mup4*p2mup5*p4mup3*F(63))+4*(p
     -   1mup4*p2mup3*p4mup5*F(19)+p1mup4*p2mup3*p3mup5*F(29)+p1mup4*
     -   p2mup5*p3mup3*F(30)+p2mup3*p2mup4*p3mup5*F(33)+p2mup4*p2mup5
     -   *p3mup3*F(35)+p2mup3*p2mup4*p4mup5*F(36)+p2mup3*p2mup5*p4mup
     -   4*F(37)+p2mup3*p3mup4*p3mup5*F(38)+p2mup5*p3mup3*p3mup4*F(39
     -   )+p3mup3*p3mup4*p3mup5*F(46)+p1mup4*p3mup3*p3mup5*F(47)+p2mu
     -   p4*p3mup3*p3mup5*F(48)+p2mup3*p4mup4*p4mup5*F(53)+p4mup3*p4m
     -   up4*p4mup5*F(57)+p1mup4*p4mup3*p4mup5*F(58)+p2mup4*p4mup3*p4
     -   mup5*F(59)+p2mup5*p4mup3*p4mup4*F(64))+mup3mup5*p1mup4*F(510
     -   )+mup3mup5*p2mup4*F(511)+p1mup5*(2*(p3mup3*(2*(p3mup4*F(15)+
     -   p2mup4*F(18))+p4mup4*F(25))+p3mup4*p4mup3*F(26))+4*(p2mup3*(
     -   2*(p2mup4*F(13)+p3mup4*F(17))+p4mup4*F(20))+p4mup3*p4mup4*F(
     -   60))+8*p2mup4*p4mup3*F(61)+mup3mup4*F(512))+mup3mup4*p2mup5*
     -   F(513)+2*(p1mup4*p3mup3*p4mup5*F(23)+p2mup3*p3mup4*p4mup5*F(
     -   40)+p2mup3*p3mup5*p4mup4*F(41)+p2mup4*p3mup3*p4mup5*F(42)+p2
     -   mup5*p3mup3*p4mup4*F(44)+p2mup5*p3mup4*p4mup3*F(45)+p3mup3*p
     -   3mup4*p4mup5*F(49)+p3mup3*p3mup5*p4mup4*F(50)+p1mup3*(4*p1mu
     -   p4*p1mup5*F(12)+p3mup4*(2*p3mup5*F(14)+p4mup5*F(21)+4*p2mup5
     -   *F(28))+2*p2mup4*(p3mup5*F(16)+2*p2mup5*F(27)+p4mup5*F(31))+
     -   p4mup4*(p3mup5*F(22)+2*(p2mup5*F(32)+p4mup5*F(51))))+p3mup4*
     -   p3mup5*p4mup3*F(52)+p3mup3*p4mup4*p4mup5*F(54)+p3mup4*p4mup3
     -   *p4mup5*F(55)+p3mup5*p4mup3*p4mup4*F(56)+p1mup4*p3mup5*p4mup
     -   3*F(514)+p2mup4*p3mup5*p4mup3*F(515)))
       Fa1(76) = -8*(p2mup3*p2mup4*p4mup5*F(90)+p3mup3*p3mup4*p4mup5*F
     -   (91)+p3mup4*p3mup5*p4mup3*F(92)+p2mup3*p4mup4*p4mup5*F(93)+p
     -   3mup3*p4mup4*p4mup5*F(94)+p3mup4*p4mup3*p4mup5*F(95)+p2mup3*
     -   p3mup4*p4mup5*F(96)+p2mup4*p3mup3*p4mup5*F(97)+p2mup4*p3mup5
     -   *p4mup3*F(98)+p3mup3*p3mup5*p4mup4*F(107)+p2mup3*p3mup5*p4mu
     -   p4*F(108)+p2mup3*p2mup5*p4mup4*F(112)+p2mup5*p3mup3*p4mup4*F
     -   (114)+p2mup4*p2mup5*p4mup3*F(115)+p2mup4*p4mup3*p4mup5*F(116
     -   )+p2mup5*p3mup4*p4mup3*F(117)+p3mup5*p4mup3*p4mup4*F(118)+p2
     -   mup5*p4mup3*p4mup4*F(119))-16*(p3mup3*p3mup4*p3mup5*F(99)+p1
     -   mup4*p3mup3*p3mup5*F(100)+p1mup4*p2mup3*p3mup5*F(102)+p2mup3
     -   *p2mup4*p3mup5*F(105)+p2mup4*p3mup3*p3mup5*F(106)+p2mup3*p3m
     -   up4*p3mup5*F(106)+p1mup4*p2mup3*p2mup5*F(109)+p1mup4*p2mup5*
     -   p3mup3*F(110)+p2mup4*p2mup5*p3mup3*F(111)+p2mup3*p2mup5*p3mu
     -   p4*F(111)+p2mup5*p3mup3*p3mup4*F(113)+p4mup3*p4mup4*p4mup5*F
     -   (120))+mup4mup5*p1mup3*F(516)+mup3mup5*p1mup4*F(517)+mup3mup
     -   5*p2mup4*F(519)+2*(mup4mup5*p3mup3*F(71)+mup4mup5*p2mup3*F(5
     -   20)+mup3mup4*p2mup5*F(521))+p1mup5*(mup3mup4*F(518)-4*(p2mup
     -   4*p4mup3*F(78)+p3mup4*p4mup3*F(85)+p4mup3*p4mup4*F(88)+2*(p3
     -   mup3*p3mup4*F(101)+p2mup4*p3mup3*F(103))+p3mup3*p4mup4*F(104
     -   )-p2mup3*(2*(p2mup4*F(72)+p3mup4*F(522))+p4mup4*F(523))))+4*
     -   (p1mup4*p2mup3*p4mup5*F(75)+p1mup4*p2mup5*p4mup3*F(76)+p1mup
     -   4*p3mup5*p4mup3*F(84)-p1mup3*(p2mup5*p4mup4*F(44)+2*(p2mup5*
     -   p3mup4*F(39)+p1mup4*p1mup5*F(79)+p3mup4*p3mup5*F(80))+p2mup4
     -   *(p4mup5*F(74)+2*(p2mup5*F(35)+p3mup5*F(81)))+p3mup4*p4mup5*
     -   F(82)+p3mup5*p4mup4*F(83)+p4mup4*p4mup5*F(86))-p1mup4*p3mup3
     -   *p4mup5*F(89)+p1mup4*p4mup3*p4mup5*F(524))
       Fa1(77) = -8*(p3mup3*p3mup5*p4mup4*F(95)+p2mup4*p2mup5*p3mup3*F
     -   (115)-p1mup5*p2mup3*p2mup4*F(128)+p2mup3*p2mup4*p3mup5*F(145
     -   )+p2mup3*p3mup4*p3mup5*F(146)+p2mup4*p3mup3*p3mup5*F(147)+p3
     -   mup3*p3mup4*p4mup5*F(148)+p2mup3*p3mup4*p4mup5*F(149)+p2mup3
     -   *p3mup5*p4mup4*F(150)+p1mup5*p2mup3*p4mup4*F(154)+p1mup5*p2m
     -   up4*p4mup3*F(155)+p3mup4*p3mup5*p4mup3*F(158)+p3mup3*p4mup4*
     -   p4mup5*F(160)+p3mup4*p4mup3*p4mup5*F(161)+p2mup4*p3mup3*p4mu
     -   p5*F(162)+p2mup4*p3mup5*p4mup3*F(163)+p2mup3*p2mup5*p3mup4*F
     -   (164)+p2mup5*p3mup3*p3mup4*F(165)+p2mup5*p3mup3*p4mup4*F(168
     -   )+p2mup5*p3mup4*p4mup3*F(169)+p1mup5*p4mup3*p4mup4*F(170)+p3
     -   mup5*p4mup3*p4mup4*F(173))-16*(p3mup3*p3mup4*p3mup5*F(141)+p
     -   2mup3*p2mup4*p4mup5*F(157)+p2mup3*p4mup4*p4mup5*F(159)+p2mup
     -   3*p2mup5*p4mup4*F(167)+p2mup4*p2mup5*p4mup3*F(171)+p2mup5*p4
     -   mup3*p4mup4*F(172)+p4mup3*p4mup4*p4mup5*F(174)+p2mup4*p4mup3
     -   *p4mup5*F(175))+mup4mup5*p1mup3*F(525)+mup3mup5*p1mup4*F(526
     -   )+mup3mup4*p1mup5*F(527)+mup4mup5*p2mup3*F(528)+mup4mup5*p3m
     -   up3*F(529)-2*(mup3mup5*p2mup4*F(126)-mup3mup4*p2mup5*F(530))
     -   +4*(p1mup5*p2mup3*p3mup4*F(132)-p1mup4*(p2mup5*p3mup3*F(44)+
     -   p3mup3*p3mup5*F(134)+p2mup3*(p3mup5*F(131)+2*(p2mup5*F(37)+p
     -   4mup5*F(138)))+p3mup3*p4mup5*F(140)+p3mup5*p4mup3*F(143)+2*(
     -   p2mup5*p4mup3*F(64)+p1mup3*p1mup5*F(136)+p4mup3*p4mup5*F(152
     -   )))-p1mup5*p3mup3*p4mup4*F(156)-p1mup3*(p2mup5*p3mup4*F(130)
     -   +p3mup4*p3mup5*F(137)+p3mup4*p4mup5*F(139)+p3mup5*p4mup4*F(1
     -   42)+p2mup4*(-(p3mup5*F(129))+4*(p2mup5*F(63)+p4mup5*F(153)))
     -   +4*(p4mup4*p4mup5*F(151)+p2mup5*p4mup4*F(166)))+p1mup5*p2mup
     -   4*p3mup3*F(531)+p1mup5*p3mup3*p3mup4*F(532)+p1mup5*p3mup4*p4
     -   mup3*F(533))
       Fa1(78) = (mup4mup5*F(534))/2.+p2mup4*p3mup5*F(535)+p2mup5*p3mu
     -   p4*F(536)+p2mup4*p4mup5*F(537)+p2mup5*p4mup4*F(538)+p3mup4*p
     -   4mup5*F(541)+p3mup5*p4mup4*F(542)-2*(p2mup4*p2mup5*F(540)-p3
     -   mup4*p3mup5*F(543))+p1mup4*(2*(p3mup5*F(185)-p2mup5*F(539))+
     -   p1mup5*F(544)+p4mup5*F(546))+p1mup5*(-2*p2mup4*F(187)+p3mup4
     -   *F(545)+p4mup4*F(547))+p4mup4*p4mup5*F(548)
       Fa1(79) = (mup4mup5*F(549))/2.+p1mup4*p3mup5*F(550)+p1mup5*p3mu
     -   p4*F(551)+p1mup4*p4mup5*F(554)+p1mup5*p4mup4*F(555)+p3mup4*p
     -   4mup5*F(556)+p3mup5*p4mup4*F(557)+p2mup4*(-2*(p4mup5*F(205)-
     -   p1mup5*F(553))+p3mup5*F(558))+p3mup4*p3mup5*F(559)-2*(p2mup5
     -   *(p3mup4*F(204)+p4mup4*F(206)+2*p2mup4*F(207)-p1mup4*F(208))
     -   -p1mup4*p1mup5*F(552)-p4mup4*p4mup5*F(560))
       Fa1(80) = 4*p3mup4*p3mup5*F(225)+(mup4mup5*F(561))/2.+p1mup5*(-
     -   2*p3mup4*F(214)+p2mup4*F(563)+p4mup4*F(565))+p2mup4*p4mup5*F
     -   (566)+p2mup5*p4mup4*F(567)+p2mup4*p2mup5*F(568)+p2mup5*p3mup
     -   4*F(569)+p1mup4*(p2mup5*F(562)+p4mup5*F(564)-2*(p1mup5*F(213
     -   )-p3mup5*F(571)))-2*(p3mup4*p4mup5*F(221)-p3mup5*p4mup4*F(57
     -   0)-p2mup4*p3mup5*F(572)-p4mup4*p4mup5*F(573))
       Fa1(81) = (mup4mup5*F(574))/2.+p1mup4*p2mup5*F(575)+p1mup5*p3mu
     -   p4*F(578)+p2mup5*p3mup4*F(580)+p2mup4*(p1mup5*F(576)+p3mup5*
     -   F(579)-2*(p4mup5*F(234)-p2mup5*F(582)))+p1mup4*p1mup5*F(583)
     -   +2*(p2mup5*p4mup4*F(581)+p3mup4*p4mup5*F(584))+p1mup4*p4mup5
     -   *F(585)+p3mup5*(p1mup4*F(577)+2*(p4mup4*F(239)+p3mup4*F(586)
     -   ))+p1mup5*p4mup4*F(587)+p4mup4*p4mup5*F(588)
       Fa1(82) = p1mup4*F(589)+p3mup4*F(590)+p2mup4*F(591)+p4mup4*F(59
     -   2)
       Fa1(83) = 8*(p2mup3*p2mup4*p3mup5*F(265)+p2mup3*p2mup5*p3mup4*F
     -   (267)+p2mup3*p3mup4*p3mup5*F(270)+p2mup4*p2mup5*p3mup3*F(280
     -   )+p2mup5*p3mup3*p3mup4*F(281))+mup4mup5*p2mup3*F(595)+mup3mu
     -   p5*(2*p1mup4*F(253)+p2mup4*F(596))+mup3mup4*p2mup5*F(597)-4*
     -   (-2*p1mup5*p2mup3*(p2mup4*F(255)+p3mup4*F(257))+p1mup3*(2*(p
     -   2mup4*p2mup5*F(254)+p2mup5*p3mup4*F(256))+p2mup4*p4mup5*F(25
     -   8)+p2mup5*p4mup4*F(259))-p1mup5*p2mup3*p4mup4*F(262)+p1mup5*
     -   p2mup4*p4mup3*F(263)-p2mup3*p2mup4*p4mup5*F(268)-p2mup3*p2mu
     -   p5*p4mup4*F(269)-p2mup3*p3mup4*p4mup5*F(271)-p2mup3*p3mup5*p
     -   4mup4*F(272)+p2mup4*p3mup5*p4mup3*F(273)-p2mup5*p3mup4*p4mup
     -   3*F(274)+p2mup4*p2mup5*p4mup3*F(275)-p2mup3*p4mup4*p4mup5*F(
     -   276)+p2mup4*p4mup3*p4mup5*F(277)-p1mup4*(p2mup3*p4mup5*F(260
     -   )+p2mup5*p4mup3*F(261)+2*(p2mup3*p2mup5*F(265)+p2mup3*p3mup5
     -   *F(266)+p2mup5*p3mup3*F(279)))-p2mup4*p3mup3*p4mup5*F(282)-p
     -   2mup5*p4mup3*p4mup4*F(598)-p2mup5*p3mup3*p4mup4*F(599))+2*(m
     -   up4mup5*p1mup3*F(258)+mup3mup4*p1mup5*F(263)+mup4mup5*p3mup3
     -   *F(600))*P(7)
       Fa1(84) = 2*mup4mup5*p3mup3*F(290)-8*(p2mup3*p2mup4*p3mup5*F(30
     -   3)-p2mup4*p2mup5*p3mup3*F(304)+p2mup3*p2mup5*p4mup4*F(304)-p
     -   2mup3*p2mup4*p4mup5*F(305)+p2mup3*p3mup4*p3mup5*F(306)-p2mup
     -   4*p3mup5*p4mup3*F(306)-p2mup4*p3mup3*p4mup5*F(306)-p2mup5*p3
     -   mup3*p3mup4*F(307)+p2mup3*p3mup5*p4mup4*F(308)+p2mup5*p3mup4
     -   *p4mup3*F(309)+p2mup3*p2mup5*p3mup4*F(311)-p2mup4*p3mup3*p3m
     -   up5*F(312)-p2mup4*p2mup5*p4mup3*F(314)+p2mup3*p4mup4*p4mup5*
     -   F(315)+p2mup3*p3mup4*p4mup5*F(316)-p2mup4*p4mup3*p4mup5*F(31
     -   6)+p2mup5*p3mup3*p4mup4*F(317)+p2mup5*p4mup3*p4mup4*F(318))+
     -   mup4mup5*p1mup3*F(601)+mup4mup5*p2mup3*F(604)+mup3mup5*p2mup
     -   4*F(605)+mup3mup4*p2mup5*F(606)+p1mup5*(mup3mup4*F(603)+4*(p
     -   2mup4*(p4mup3*F(302)+p3mup3*F(310))+p2mup3*(p2mup4*F(292)-p4
     -   mup4*F(313)+p3mup4*F(607))))+p1mup4*(-4*(p2mup3*(p3mup5*F(27
     -   2)+p2mup5*F(291)+p4mup5*F(298))+p2mup5*p4mup3*F(301))+mup3mu
     -   p5*F(602)+p2mup5*p3mup3*F(608))+4*(p1mup3*p2mup4*p3mup5*F(29
     -   3)-p1mup3*p2mup5*p3mup4*F(294)+p1mup3*p2mup4*p2mup5*F(296)+p
     -   1mup3*p2mup4*p4mup5*F(297)+p1mup3*p2mup5*p4mup4*F(609))
       Fa1(85) = p1mup2*F(610)+p3mup2*F(611)+p4mup2*F(612)
       Fa1(86) = p1mup3*F(615)+p2mup3*F(616)
       Fa1(87) = 8*(-(p1mup4*p2mup3*p2mup5*F(343))+p1mup4*p2mup5*p4mup
     -   3*F(359)+p1mup3*p2mup5*p4mup4*F(359)+p2mup3*p2mup5*p4mup4*F(
     -   360)-p2mup5*p4mup3*p4mup4*F(365))+mup4mup5*p2mup3*F(625)+mup
     -   3mup5*p2mup4*F(626)+mup3mup4*p2mup5*F(627)+mup4mup5*p3mup3*F
     -   (629)+2*(mup4mup5*p1mup3*F(628)+mup3mup5*p1mup4*F(630))+p1mu
     -   p5*(4*(p2mup3*p3mup4*F(263)+p2mup4*p3mup3*F(348))+8*(p2mup3*
     -   p2mup4*F(255)+p2mup4*p4mup3*F(364))+mup3mup4*F(631))-4*(p1mu
     -   p3*p2mup5*p3mup4*F(345)+p1mup4*p2mup3*p3mup5*F(346)+p1mup4*p
     -   2mup5*p3mup3*F(347)+p2mup3*p2mup5*p3mup4*F(350)+p2mup5*p3mup
     -   3*p3mup4*F(352)-p2mup3*p3mup4*p3mup5*F(356)-p2mup3*p3mup4*p4
     -   mup5*F(357)+p2mup5*p3mup4*p4mup3*F(358)-p2mup4*(p3mup5*(p4mu
     -   p3*F(302)+p1mup3*F(344)+p2mup3*F(349))+p3mup3*(p4mup5*F(297)
     -   +p3mup5*F(351)+p2mup5*F(355))+2*(p2mup5*(p4mup3*F(314)+p1mup
     -   3*F(353))+p4mup5*(p2mup3*F(353)+p1mup3*F(354)+p4mup3*F(361))
     -   ))-p2mup3*p3mup5*p4mup4*F(362)-p2mup5*p3mup3*p4mup4*F(632))
       Fa1(88) = p2mup1*F(368)+p3mup1*F(369)+p4mup1*F(371)
       Fa1(89) = p2mup5*p3mup4*p4mup2*F(379)+2*(p1mup5*p3mup2*(p2mup4*
     -   F(255)+p3mup4*F(257)+p1mup4*F(634))+p1mup2*p3mup5*(p1mup4*F(
     -   633)+p2mup4*F(635)+p3mup4*F(636)))+p3mup2*(p2mup4*p4mup5*F(6
     -   39)+2*(p2mup5*p3mup4*F(637)+p1mup4*p2mup5*F(638)+p3mup4*p3mu
     -   p5*F(642)+p2mup4*p3mup5*F(643)+p1mup4*p3mup5*F(644))+p3mup4*
     -   p4mup5*F(645)+p4mup4*(p4mup5*F(276)+p2mup5*F(641)+p3mup5*F(6
     -   46)))+p3mup5*p4mup2*(-(p2mup4*F(640))+p3mup4*F(647))
       Fa1(90) = 4*(p3mup2*p3mup4*p3mup5*F(106)+p2mup5*p3mup2*p3mup4*F
     -   (111)+p4mup2*p4mup4*p4mup5*F(172))+p1mup2*p1mup5*p4mup4*F(65
     -   2)+p1mup2*p2mup5*p4mup4*F(655)+p1mup2*p4mup4*p4mup5*F(662)-p
     -   1mup5*p4mup2*p4mup4*F(669)+p1mup4*(p1mup2*(p4mup5*F(651)+2*(
     -   p2mup5*F(648)+p3mup5*F(658)))-p4mup2*(p2mup5*F(656)+p1mup5*F
     -   (661)+p4mup5*F(663))+2*p3mup2*(2*(p3mup5*F(102)+p2mup5*F(109
     -   ))+p1mup5*F(670)))+p2mup4*(p1mup2*p4mup5*F(654)+p1mup5*p4mup
     -   2*F(657)+2*(p1mup2*p1mup5*F(649)+p1mup2*p3mup5*F(659)+p3mup2
     -   *(2*p3mup5*F(105)+p4mup5*F(664)+p1mup5*F(671))+p4mup2*(p3mup
     -   5*F(665)+p2mup5*F(673)+p4mup5*F(674))))+2*(p3mup2*p3mup5*p4m
     -   up4*F(108)+p2mup5*p3mup2*p4mup4*F(112)+p1mup2*p1mup5*p3mup4*
     -   F(650)+p1mup2*p2mup5*p3mup4*F(653)+p1mup2*p3mup4*p3mup5*F(66
     -   0)+p3mup2*p3mup4*p4mup5*F(666)+p3mup4*p3mup5*p4mup2*F(667)+p
     -   3mup2*p4mup4*p4mup5*F(668)+p1mup5*p3mup2*p3mup4*F(672)+p2mup
     -   5*p3mup4*p4mup2*F(675)+p2mup5*p4mup2*p4mup4*F(676))
       Fa1(91) = -(p3mup5*(p1mup2*p1mup4*F(677)+p1mup2*p2mup4*F(679)-2
     -   *p2mup4*p4mup2*F(689)))+p3mup2*(p1mup4*(p2mup5*F(641)+p3mup5
     -   *F(682)+p1mup5*F(684))+p2mup4*(-(p1mup5*F(681))+2*(p4mup5*F(
     -   314)-p3mup5*F(686)))+p3mup4*(4*p3mup5*F(306)-p1mup5*F(683)+2
     -   *(p4mup5*F(687)+p2mup5*F(690)))+2*p4mup4*(p4mup5*F(315)+p3mu
     -   p5*F(688)+p2mup5*F(691)))+p3mup4*(p1mup2*p1mup5*F(678)+p1mup
     -   2*p2mup5*F(680)+p1mup2*p3mup5*F(685)+2*(p3mup5*p4mup2*F(316)
     -   +p2mup5*p4mup2*F(692)))
       Fa1(92) = (mup2mup5*F(693))/2.+p1mup2*p3mup5*F(695)+p3mup2*p4mu
     -   p5*F(698)+p3mup5*p4mup2*F(699)+4*p1mup2*p1mup5*F(700)+p2mup5
     -   *(p3mup2*F(694)+p4mup2*F(696)+2*p1mup2*F(701))+p1mup5*p4mup2
     -   *F(702)-2*(p1mup2*p4mup5*F(697)-p3mup2*p3mup5*F(703)+p4mup2*
     -   p4mup5*F(704))+p1mup5*p3mup2*F(705)
       Fa1(93) = (mup2mup5*F(706))/2.+p1mup5*p3mup2*F(708)+p1mup5*p4mu
     -   p2*F(710)+p3mup2*p4mup5*F(711)+p3mup5*p4mup2*F(712)+p1mup2*(
     -   p3mup5*F(707)+p4mup5*F(709)-2*p1mup5*F(714)+p2mup5*F(716))+p
     -   3mup2*p3mup5*F(717)+2*(p2mup5*p3mup2*F(713)+p2mup5*p4mup2*F(
     -   715)+p4mup2*p4mup5*F(718))
       Fa1(94) = (mup2mup5*F(719))/2.+p1mup5*p4mup2*F(723)+p2mup5*p4mu
     -   p2*F(724)+p3mup2*p4mup5*F(725)+p1mup2*(p2mup5*F(721)+p4mup5*
     -   F(722)-2*(p1mup5*F(720)-p3mup5*F(728)))-2*(p3mup5*p4mup2*F(7
     -   26)-p2mup5*p3mup2*F(727)-p1mup5*p3mup2*F(729)+p4mup2*p4mup5*
     -   F(730))+4*p3mup2*p3mup5*F(731)
       Fa1(95) = (mup2mup5*F(732))/2.+p1mup5*p3mup2*F(736)+p2mup5*p3mu
     -   p2*F(737)+p1mup2*(p2mup5*F(734)+p3mup5*F(735)+2*(p1mup5*F(73
     -   3)-p4mup5*F(740)))+2*(p3mup2*p3mup5*F(738)+p2mup5*p4mup2*F(7
     -   39)+p1mup5*p4mup2*F(741)+p3mup2*p4mup5*F(742))+p3mup5*p4mup2
     -   *F(743)-4*p4mup2*p4mup5*F(744)
       Fa1(96) = mup4mup5*p4mup3*F(747)+mup3mup5*(p3mup4*F(745)+p4mup4
     -   *F(748))+mup3mup4*(p3mup5*F(746)+p4mup5*F(749))-16*p1mup5*(p
     -   1mup4*p2mup3+p1mup3*p2mup4)*F(750)-8*(p1mup4*p1mup5*p3mup3*F
     -   (753)+p1mup3*(p1mup4*p3mup5*F(751)+p1mup4*p4mup5*F(754)+p1mu
     -   p5*p4mup4*F(755)+2*(p1mup5*p3mup4*F(752)+p1mup4*p2mup5*F(756
     -   )))+2*(p2mup3*p2mup4*p2mup5*F(757)+p1mup4*p1mup5*p4mup3*F(75
     -   8)))
       Fa1(97) = mup4mup5*p4mup3*F(760)+mup3mup5*(p3mup4*F(759)+p4mup4
     -   *F(761))+mup3mup4*p4mup5*F(762)+2*mup3mup4*p3mup5*F(763)+4*(
     -   p1mup3*(p1mup4*(2*(p2mup5*F(30)+p3mup5*F(766))+p4mup5*F(767)
     -   )+p1mup5*(2*(p3mup4*F(47)+p2mup4*F(765))+p4mup4*F(768)))-p1m
     -   up4*p1mup5*(p4mup3*F(769)+2*(p2mup3*F(764)-p3mup3*F(770))))+
     -   16*p2mup3*p2mup4*p2mup5*F(771)
       Fa1(98) = mup3mup5*p3mup4*F(772)+mup3mup4*p3mup5*F(773)+2*(mup4
     -   mup5*p4mup3*F(774)+mup3mup5*p4mup4*F(775)-mup3mup4*p4mup5*F(
     -   776)+2*(p1mup3*(p1mup5*(2*(p4mup4*F(58)-p2mup4*F(777))-p3mup
     -   4*F(779))+p1mup4*(p3mup5*F(778)+2*(p2mup5*F(32)+p4mup5*F(782
     -   ))))+p1mup4*p1mup5*(p3mup3*F(781)+2*(p2mup3*F(780)+p4mup3*F(
     -   783))))+8*p2mup3*p2mup4*p2mup5*F(784))
       Fa1(99) = (mup3mup5*F(785))/2.+p2mup5*p3mup3*F(789)+p2mup5*p4mu
     -   p3*F(791)+p1mup3*p4mup5*F(792)+p3mup3*p4mup5*F(793)+p3mup5*p
     -   4mup3*F(794)+p2mup3*(p2mup5*F(786)+p3mup5*F(788)+p4mup5*F(79
     -   0)+p1mup5*F(795))+p1mup5*p4mup3*F(796)+2*(p1mup3*p2mup5*F(78
     -   7)+p4mup3*p4mup5*F(797)-(p1mup3*p3mup5*F(344)+p1mup5*p3mup3*
     -   F(348)+p3mup3*p3mup5*F(351)+2*p1mup3*p1mup5*F(634))*P(7))
       Return
       End
