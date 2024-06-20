c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa35(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa35)
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
       Complex*16 Fa35(817:872)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa35(817) = mup2mup3*(-(p5mup5*p5mup6*F(5688))-p2mup6*p5mup5*F(
     -   9403)-p1mup6*p5mup5*F(9487)+p3mup6*(-(p5mup5*F(5684))+2*(p1m
     -   up5*F(9632)+p2mup5*F(9639)+p3mup5*F(9644)+p4mup5*F(9646)))+p
     -   4mup6*(-(p5mup5*F(5687))+2*(p1mup5*F(9633)+p2mup5*F(9640)+p3
     -   mup5*F(9646)+p4mup5*F(9648)))+2*(p1mup6*p2mup5*F(9629)+p1mup
     -   5*p2mup6*F(9629)+p1mup6*p3mup5*F(9632)+p1mup6*p4mup5*F(9633)
     -   +p1mup5*p5mup6*F(9634)+p2mup5*p2mup6*F(9636)+p2mup6*p3mup5*F
     -   (9639)+p2mup6*p4mup5*F(9640)+p2mup5*p5mup6*F(9641)+p3mup5*p5
     -   mup6*F(9647)+p4mup5*p5mup6*F(9649)+p1mup5*p1mup6*F(9657)))
       Fa35(818) = mup2mup6*(4*(p4mup3*p4mup5*F(8733)+p4mup5*p5mup3*F(
     -   8734)+p2mup3*p4mup5*F(8735)+p3mup3*p4mup5*F(8746)+p1mup5*p4m
     -   up3*F(8747)+p1mup5*p2mup3*F(8748)+p2mup5*(p4mup3*F(7521)+p5m
     -   up3*F(7554)+p2mup3*F(8744)+p3mup3*F(8745)+p1mup3*F(8748))+p1
     -   mup5*p5mup3*F(8753)+p1mup5*p3mup3*F(8754)+p3mup5*(-(p4mup3*F
     -   (7564))-p5mup3*F(7569)+p3mup3*F(8737)+p2mup3*F(8745)+p1mup3*
     -   F(8754))+p1mup3*p4mup5*F(8755)+p1mup3*p1mup5*F(8756))-p4mup3
     -   *p5mup5*F(9483)+p1mup3*p5mup5*F(9485)-p5mup3*p5mup5*F(9489)+
     -   p2mup3*p5mup5*F(9515)+p3mup3*p5mup5*F(9516))+4*mup2mup5*(p3m
     -   up6*p5mup3*F(4061)+p4mup3*p5mup6*F(6932)+p3mup6*p4mup3*F(710
     -   9)+p5mup3*p5mup6*F(9517)+p1mup6*p4mup3*F(9518)+p1mup6*p5mup3
     -   *F(9519)+p1mup6*p2mup3*F(9628)+p1mup6*p3mup3*F(9630)+p1mup3*
     -   p3mup6*F(9630)+p1mup3*p5mup6*F(9632)+p2mup3*p3mup6*F(9637)+p
     -   2mup6*(p5mup3*F(4045)+p4mup3*F(9514)+p1mup3*F(9628)+p2mup3*F
     -   (9635)+p3mup3*F(9637))+p2mup3*p5mup6*F(9639)+p4mup6*(p5mup3*
     -   F(4057)+p4mup3*F(6971)+p1mup3*F(9631)+p2mup3*F(9638)+p3mup3*
     -   F(9642))+p3mup3*p3mup6*F(9643)+p3mup3*p5mup6*F(9644)+p1mup3*
     -   p1mup6*F(9656))
       Fa35(819) = mup2mup3*p5mup5*p6mup6*F(9520)+4*mup2mup6*p6mup3*(-
     -   (p5mup5*F(5723))+4*(p2mup5*F(9618)+p4mup5*F(9620)+p3mup5*F(9
     -   624)+p1mup5*F(9626)))+8*(2*(mup5mup6*p4mup2*p6mup3*F(8631)-m
     -   up2mup5*p6mup3*(p2mup6*F(9618)+p4mup6*F(9620)+p5mup6*F(9621)
     -   +p3mup6*F(9624)+p1mup6*F(9626)))+mup2mup3*(p6mup6*(p2mup5*F(
     -   9619)+p3mup5*F(9621)+p4mup5*F(9622)+p1mup5*F(9627))+2*(p2mup
     -   6*p6mup5*F(9619)+p6mup5*(p3mup6*F(9621)+p4mup6*F(9622)+p5mup
     -   6*F(9625)+p1mup6*F(9627)))))
       Fa35(820) = mup2mup6*p6mup5*(p2mup3*F(9619)+p3mup3*F(9621)+p4mu
     -   p3*F(9622)+p5mup3*F(9625)+p1mup3*F(9627))
       Fa35(821) = 2*mup5mup6*p3mup2*p6mup3*F(8682)+mup2mup5*p5mup3*p6
     -   mup6*F(9521)-mup2mup3*p6mup5*p6mup6*F(9617)
       Fa35(822) = mup2mup5*p6mup6*(-(p4mup3*F(9522))+p2mup3*F(9618)+p
     -   3mup3*F(9624)+p1mup3*F(9626))
       Fa35(823) = p1mup6*F(9523)+p5mup6*F(9524)+p6mup6*F(9525)+p2mup6
     -   *F(9526)+p3mup6*F(9527)+p4mup6*F(9528)
       Fa35(824) = p6mup2*p6mup6*F(9539)-2*(mup2mup6*F(9529)-p2mup2*(8
     -   *(p3mup6*F(5664)+p2mup6*F(8343)+p1mup6*F(8389)+p4mup6*F(8674
     -   ))-4*p5mup6*F(9530)+p6mup6*F(9531))-p1mup2*p6mup6*F(9533)-p4
     -   mup2*p6mup6*F(9536)-p5mup2*p6mup6*F(9537)-p3mup2*p6mup6*F(95
     -   41))-8*(p1mup2*p5mup6*F(9532)-p5mup2*p5mup6*F(9534)-p4mup2*p
     -   5mup6*F(9535)-p5mup6*p6mup2*F(9538)-p3mup2*p5mup6*F(9540)-p4
     -   mup6*p5mup2*F(9543)-p3mup6*p5mup2*F(9544)-p1mup6*p5mup2*F(95
     -   45)-p2mup6*(p5mup2*F(9542)+2*(p4mup2*F(8650)+p1mup2*F(8667)-
     -   p6mup2*F(9618)+p3mup2*F(9635))))+16*(p1mup6*p4mup2*F(8656)+p
     -   1mup2*p1mup6*F(8668)+p1mup2*p4mup6*F(8673)+p1mup2*p3mup6*F(8
     -   676)+p3mup6*p4mup2*F(8678)+p4mup2*p4mup6*F(8679)-p4mup6*p6mu
     -   p2*F(9620)-p3mup6*p6mup2*F(9624)-p1mup6*p6mup2*F(9626)+p1mup
     -   6*p3mup2*F(9628)+p3mup2*p3mup6*F(9637)+p3mup2*p4mup6*F(9638)
     -   )
       Fa35(825) = p1mup5*F(9546)+p5mup5*F(9547)+p6mup5*F(9548)+p2mup5
     -   *F(9549)+p3mup5*F(9550)+p4mup5*F(9551)
       Fa35(826) = 16*(p2mup5*p3mup2*F(4323)+p3mup2*p4mup5*F(4326)+p3m
     -   up2*p3mup5*F(4498)-p3mup5*p5mup2*F(5681)-p3mup5*p4mup2*F(573
     -   3)+p2mup5*p5mup2*F(8194)-p1mup2*p4mup5*F(8723)-p1mup2*p3mup5
     -   *F(8729)+p1mup2*p2mup5*F(8738)+p2mup5*p4mup2*F(8757)+p4mup2*
     -   p4mup5*F(8758)+p4mup5*p5mup2*F(8759)+p1mup5*(p3mup2*F(4320)-
     -   p1mup2*F(8742)+p4mup2*F(8760)+p5mup2*F(8761)))+p5mup5*p6mup2
     -   *F(9560)+2*(p2mup2*(8*(p2mup5*F(4109)-p3mup5*F(4111)-p4mup5*
     -   F(4114)-p1mup5*F(4122))+p5mup5*F(9552)+4*p6mup5*F(9553))+p1m
     -   up2*p5mup5*F(9554)+p4mup2*p5mup5*F(9556)+p5mup2*p5mup5*F(955
     -   8)+p3mup2*p5mup5*F(9561))+p2mup5*p6mup2*F(9563)+p4mup5*p6mup
     -   2*F(9564)-8*(p6mup2*p6mup5*F(8331)+p3mup5*p6mup2*F(9457)-p1m
     -   up2*p6mup5*F(9555)-p4mup2*p6mup5*F(9557)-p5mup2*p6mup5*F(955
     -   9)-p3mup2*p6mup5*F(9562)+mup2mup5*F(9565))+p1mup5*p6mup2*F(9
     -   566)
       Fa35(827) = 2*mup5mup6*p1mup2*F(17)+8*p6mup2*p6mup5*p6mup6*F(54
     -   )+Fa(3)+Fa(4)+Fa(5)+Fa(6)+Fa(7)-4*Fa(8)+Fa(9)-Fa(10)+Fa(11)-
     -   16*Fa(12)
       Fa35(828) = mup5mup6*p4mup3*p6mup2*F(87)+32*p3mup5*p3mup6*p4mup
     -   3*p6mup2*F(112)+F(87)*(Fa(16)+Fa(17))+Fa(18)+Fa(19)+16*(Fa(2
     -   3)+Fa(29)+Fa(32))+8*(Fa(24)-Fa(25)+Fa(26)+Fa(27)+Fa(28)+Fa(3
     -   0)+Fa(31)+Fa(33)+Fa(34)-Fa(35)+Fa(36))+4*(Fa(20)+Fa(22)+F(10
     -   7)*Fa(21)*P(59))
       Fa35(829) = -(mup5mup6*p4mup3*p6mup2*F(123))+32*p3mup5*p3mup6*p
     -   4mup3*p6mup2*F(9574)-F(123)*(Fa(16)+Fa(17))+Fa(37)+Fa(38)+Fa
     -   (39)+4*(F(143)*Fa(21)+Fa(40))+8*(Fa(41)+Fa(42)-Fa(43)+Fa(44)
     -   +Fa(45)+Fa(46)+Fa(47)+Fa(48)+Fa(49)-Fa(50)-Fa(51)-Fa(52)-Fa(
     -   53)+Fa(54))
       Fa35(830) = -(mup5mup6*p1mup2*F(271))-F(271)*Fa(66)+Fa(68)+2*Fa
     -   (69)-Fa(70)+8*(p6mup2*p6mup5*p6mup6*F(9569)+Fa(71))+Fa(73)+4
     -   *(Fa(67)-Fa(72)-Fa(74))-16*F(112)*Fa(75)
       Fa35(831) = mup5mup6*p1mup2*F(280)+F(280)*Fa(66)-16*F(9574)*Fa(
     -   75)+Fa(77)-Fa(78)+Fa(79)-8*(p6mup2*p6mup5*p6mup6*F(9568)+Fa(
     -   80))+Fa(82)-4*(Fa(76)-Fa(81)+Fa(83))
       Fa35(832) = F(305)*Fa(90)+Fa(92)+Fa(93)+Fa(94)-4*(Fa(91)-2*F(10
     -   9)*Fa(95)-Fa(96)-Fa(97)-Fa(98)+4*F(9575)*Fa(99))
       Fa35(833) = -2*(mup5mup6*p1mup2*F(375)-Fa(116))+Fa(117)+Fa(118)
     -   +Fa(119)-4*(p6mup2*p6mup5*p6mup6*s12*F(425)-Fa(120)-Fa(121)-
     -   Fa(123)-s12*(Fa(122)+Fa(124)))+16*Fa(125)
       Fa35(834) = mup5mup6*(p1mup2+p2mup2)*F(607)-16*p6mup2*p6mup5*p6
     -   mup6*s12*F(9567)+Fa(162)+4*(Fa(160)+Fa(161)+Fa(163)+Fa(164)+
     -   Fa(166)+s12*(Fa(165)+Fa(167))+Fa(168))
       Fa35(835) = 8*(Fa(169)+Fa(170)+F(586)*Fa(171))-4*(mup5mup6*(-(p
     -   5mup3*p6mup2)+p5mup2*p6mup3)*F(630)-2*(p1mup3*p1mup5*p1mup6*
     -   p5mup2*F(586)-Fa(172)+Fa(173)+Fa(174)-Fa(175)+Fa(176)-Fa(177
     -   )-Fa(178)-Fa(179)-Fa(180)+Fa(181)-Fa(182)))
       Fa35(836) = 4*F(630)*Fa(185)+8*(p1mup2*p1mup3*p1mup6*p6mup5*F(5
     -   71)+p1mup2*p1mup3*p1mup5*p6mup6*F(576)-p2mup2*p2mup3*p5mup5*
     -   p5mup6*F(585)-F(586)*Fa(186)-Fa(187)-Fa(188)-Fa(189)-Fa(190)
     -   -Fa(191)-Fa(192)+Fa(193)+Fa(194)+Fa(195)+Fa(196))
       Fa35(837) = -4*F(630)*Fa(197)-8*(Fa(198)+Fa(199)-Fa(200)-Fa(201
     -   )-Fa(202)-Fa(203)-Fa(204)-Fa(205)-Fa(206)+Fa(207)+Fa(208)+2*
     -   Fa(209))
       Fa35(838) = mup5mup6*p6mup3*F(684)+4*(Fa(214)+Fa(215)+Fa(216)-F
     -   a(217)-Fa(218)+Fa(219)+Fa(220)+Fa(221)+Fa(222)+Fa(223))+8*(p
     -   2mup3*p6mup5*p6mup6*F(581)+Fa(224))
       Fa35(839) = -2*mup5mup6*p5mup3*F(685)+4*(Fa(225)+Fa(226)+Fa(227
     -   )+Fa(228)+Fa(229)+Fa(230)+Fa(231)-Fa(232)-Fa(233))+8*Fa(234)
       Fa35(840) = mup5mup6*p1mup2*F(751)-4*p6mup2*p6mup5*p6mup6*F(887
     -   )+Fa(256)+Fa(257)+Fa(258)+Fa(259)+Fa(260)+Fa(261)+Fa(262)+Fa
     -   (263)+Fa(264)+16*Fa(265)
       Fa35(841) = Fa(283)+2*(mup5mup6*p1mup2*F(1235)+Fa(281)-Fa(284)-
     -   Fa(285))-4*(Fa(282)-Fa(288))+8*(p6mup2*p6mup5*p6mup6*F(1361)
     -   +Fa(286)-Fa(287)+Fa(289))-16*Fa(290)
       Fa35(842) = Fa(292)+Fa(293)+2*(-(mup5mup6*p1mup2*F(1365))+Fa(29
     -   1)+Fa(294)+Fa(295)+4*Fa(296)+2*(Fa(297)-Fa(298))+8*(Fa(299)+
     -   Fa(300)))
       Fa35(843) = Fa(312)+Fa(313)+Fa(314)+4*(-(mup5mup6*p1mup2*F(1563
     -   ))+Fa(310)+Fa(311)+2*Fa(315)+Fa(316)+Fa(317)-4*(Fa(318)-Fa(3
     -   19)))
       Fa35(844) = mup5mup6*p1mup3*F(1799)+mup3mup6*p3mup5*F(1810)+Fa(
     -   334)+Fa(335)+Fa(336)+Fa(338)+4*(mup3mup5*p3mup6*F(1828)+Fa(3
     -   37)+Fa(339)+Fa(340)+Fa(341)-Fa(342)-Fa(343)+Fa(347))-8*(p5mu
     -   p5*p6mup3*p6mup6*F(1911)-Fa(344)-Fa(345)-Fa(346)-Fa(348))
       Fa35(845) = Fa(349)+Fa(350)-Fa(351)+4*(s12*Fa(352)-Fa(353)+Fa(3
     -   54)+s12*(4*p5mup2*p6mup3*p6mup6*F(1911)+mup3mup6*p6mup2*F(19
     -   40)+p1mup6*p2mup3*p6mup2*F(9603)-Fa(355)+Fa(356)+Fa(357)))
       Fa35(846) = -8*p6mup2*p6mup5*p6mup6*s12*F(1911)+mup5mup6*p1mup2
     -   *F(1984)+Fa(360)+Fa(361)+Fa(362)+Fa(363)+Fa(364)+Fa(366)+s12
     -   *(Fa(359)-4*(Fa(365)-Fa(367)))+Fa(368)
       Fa35(847) = 4*(Fa(370)+Fa(371)+Fa(372))+Fa(373)+Fa(374)+Fa(375)
     -   +Fa(376)+Fa(377)+16*Fa(379)+32*Fa(381)+Fa(384)-Fa(387)+8*(Fa
     -   (378)+Fa(380)+Fa(382)+Fa(383)+Fa(385)+Fa(386)-4*(Fa(388)-Fa(
     -   389)+Fa(390))-Fa(391)-Fa(392)-Fa(393))
       Fa35(848) = Fa(395)+Fa(396)+Fa(397)-8*F(2314)*Fa(398)+Fa(400)+F
     -   a(401)+4*(Fa(399)+Fa(402))+Fa(403)+s12*(Fa(404)+Fa(405)-8*(p
     -   1mup5*p4mup3*p6mup2*F(2342)-p1mup5*p2mup3*p6mup2*F(2349)+Fa(
     -   406)))
       Fa35(849) = -Fa(407)-Fa(408)+Fa(410)+Fa(411)+4*Fa(412)+Fa(414)+
     -   8*(mup5mup6*p6mup2*p6mup3*s12*F(2376)+s12*(Fa(409)-2*Fa(413)
     -   )-Fa(416)-Fa(417)-Fa(418)+Fa(419)+Fa(420)+Fa(421)-Fa(422)+4*
     -   (Fa(423)-Fa(424)+F(9579)*Fa(425))+Fa(426)+Fa(427)+2*(Fa(415)
     -   +Fa(428)))
       Fa35(850) = 4*(Fa(429)+Fa(430))+Fa(431)+Fa(432)+Fa(433)-16*(Fa(
     -   436)+Fa(438))+Fa(439)+Fa(441)+Fa(442)+8*(Fa(434)-Fa(435)-Fa(
     -   437)-Fa(440)+4*(Fa(443)+Fa(444)-Fa(445))+Fa(446)+Fa(447)-Fa(
     -   448))
       Fa35(851) = -2*mup5mup6*p1mup3*F(2058)+mup3mup6*p3mup5*F(2468)+
     -   mup3mup5*p3mup6*F(2470)-Fa(453)+Fa(454)+Fa(455)+Fa(456)-16*(
     -   Fa(463)-F(9579)*Fa(464))-4*(Fa(458)+Fa(459)+Fa(460)-Fa(462)-
     -   Fa(466))+8*(Fa(457)-Fa(461)-Fa(465)+Fa(467)+Fa(468))
       Fa35(852) = mup3mup6*p3mup5*F(2478)+2*(mup5mup6*p1mup3*F(2063)+
     -   Fa(469))+Fa(470)+Fa(471)+4*(mup3mup5*p3mup6*F(2408)+Fa(472)-
     -   Fa(474)-Fa(475)+Fa(476)+Fa(477)+Fa(478))-16*(Fa(479)+Fa(480)
     -   )+Fa(482)+Fa(483)-8*(Fa(473)+Fa(481)-Fa(484))
       Fa35(853) = mup5mup6*p1mup3*F(2706)+mup3mup6*p3mup5*F(2717)+Fa(
     -   521)+Fa(522)+Fa(523)+Fa(525)+Fa(526)+8*(Fa(532)+Fa(533))+Fa(
     -   534)+4*(mup3mup5*p3mup6*F(2744)+Fa(524)+Fa(527)+Fa(528)-Fa(5
     -   29)-Fa(530)-Fa(531)+Fa(535))+Fa(536)
       Fa35(854) = Fa(537)+Fa(538)+Fa(539)+Fa(540)+Fa(541)+4*(mup3mup6
     -   *p6mup2*F(3046)+p1mup6*p2mup3*p6mup2*F(3129)+Fa(542)+Fa(543)
     -   )+Fa(544)+Fa(545)+Fa(546)
       Fa35(855) = mup5mup6*p1mup2*F(3173)+8*p6mup2*p6mup5*p6mup6*F(34
     -   02)+Fa(548)+Fa(549)+Fa(550)+Fa(551)+Fa(552)+Fa(553)+Fa(555)+
     -   Fa(556)+4*(Fa(554)-Fa(557))
       Fa35(856) = Fa(560)+Fa(562)+Fa(563)+Fa(564)+Fa(565)+Fa(566)+Fa(
     -   569)+16*Fa(571)-8*(Fa(567)-Fa(568)-Fa(570)-Fa(572)-Fa(573)-F
     -   a(574))+Fa(575)-4*(Fa(559)-Fa(561)-Fa(577)-8*Fa(578)-Fa(579)
     -   -Fa(580)-Fa(582)-2*(Fa(576)-Fa(581)+Fa(583)))
       Fa35(857) = p1mup5*p4mup3*p6mup2*F(4699)+Fa(585)+Fa(586)+Fa(587
     -   )+Fa(589)+Fa(590)+Fa(591)+Fa(593)+4*(Fa(592)+Fa(594))+Fa(595
     -   )+8*(p1mup5*p2mup3*p6mup2*F(4700)-Fa(588)-Fa(596))
       Fa35(858) = Fa(597)+Fa(598)+8*Fa(599)+Fa(600)+Fa(601)+Fa(602)+4
     -   *(4*Fa(603)-Fa(607)+Fa(615)+2*(Fa(604)+Fa(605)-Fa(606)+Fa(60
     -   8)+Fa(609)+Fa(610)+Fa(611)-Fa(612)+Fa(613)+Fa(614)+4*Fa(616)
     -   +2*Fa(617)-Fa(618)+Fa(619)-Fa(620)+Fa(621)))
       Fa35(859) = Fa(622)+Fa(623)+Fa(625)+Fa(626)+Fa(627)+Fa(629)+Fa(
     -   632)+4*Fa(640)+32*Fa(641)+16*(Fa(628)+Fa(634)+Fa(642))+Fa(64
     -   5)-8*(Fa(624)-Fa(630)+Fa(631)-Fa(633)-Fa(635)+Fa(636)-Fa(637
     -   )-Fa(638)-Fa(639)+Fa(643)-Fa(644)+Fa(646))
       Fa35(860) = 2*mup5mup6*p1mup3*F(7375)+mup3mup6*p3mup5*F(7386)+m
     -   up3mup5*p3mup6*F(7412)+Fa(651)+Fa(652)+Fa(653)+Fa(654)+4*(Fa
     -   (657)+Fa(660))+16*(Fa(661)+Fa(662))+Fa(664)-8*(Fa(655)+Fa(65
     -   6)+Fa(658)+Fa(659)-Fa(663)-Fa(665)-Fa(666))
       Fa35(861) = 2*mup5mup6*p1mup3*F(7586)+mup3mup6*p3mup5*F(7597)+F
     -   a(667)+Fa(668)+Fa(669)+16*(Fa(677)+Fa(678))+4*(mup3mup5*p3mu
     -   p6*F(7623)+Fa(670)-Fa(672)+Fa(673)+Fa(674)-Fa(675)-Fa(676)+F
     -   a(680))+Fa(681)+8*(Fa(671)+Fa(679)-Fa(682))
       Fa35(862) = -2*(mup5mup6*p1mup2*F(8143)-Fa(697))+Fa(698)+Fa(699
     -   )+Fa(700)+Fa(701)+Fa(704)+8*(Fa(702)+Fa(703)+Fa(705))+16*(p6
     -   mup2*p6mup5*p6mup6*F(8331)+Fa(706))
       Fa35(863) = mup5mup6*p1mup2*F(8336)+Fa(707)+Fa(709)+Fa(710)+Fa(
     -   711)+4*(Fa(708)+Fa(713)+Fa(714))+8*(p6mup2*p6mup5*p6mup6*F(9
     -   616)+Fa(712)+Fa(715))+16*Fa(716)
       Fa35(864) = Fa(717)+4*(mup3mup6*p6mup2*F(7788)+Fa(718)+4*(p1mup
     -   6*p5mup3*p6mup2*F(8630)-Fa(719)+Fa(720)+Fa(721)+Fa(722))+Fa(
     -   723)+2*(p5mup3*p6mup2*p6mup6*F(8331)-Fa(724)))
       Fa35(865) = Fa(725)+Fa(728)+4*(mup3mup6*p6mup2*F(8672)+2*p5mup3
     -   *p6mup2*p6mup6*F(9616)+Fa(726)-Fa(727)-4*(p1mup6*p5mup3*p6mu
     -   p2*F(9626)-Fa(729)+Fa(730))+Fa(731)-Fa(732))
       Fa35(866) = Fa(736)-Fa(737)+Fa(738)+4*(Fa(734)+Fa(735)+Fa(739))
     -   +Fa(740)+2*(Fa(741)+2*(p1mup5*p4mup3*p6mup2*F(8747)+p1mup5*p
     -   2mup3*p6mup2*F(8748)+mup2mup3*p4mup5*F(9653)+Fa(742)+Fa(743)
     -   ))
       Fa35(867) = 8*Fa(745)+Fa(748)+Fa(749)+Fa(750)+Fa(752)+2*Fa(753)
     -   +4*(-(p1mup5*p4mup3*p6mup2*F(7704))-p1mup5*p2mup3*p6mup2*F(9
     -   628)+Fa(744)+Fa(746)-Fa(747)+Fa(751)-Fa(754)+Fa(755))
       Fa35(868) = Fa(760)+Fa(761)+Fa(762)+Fa(764)+4*(-(mup5mup6*p1mup
     -   2*F(8958))+Fa(763)+2*Fa(765)+Fa(766)+Fa(767)-4*(Fa(768)+Fa(7
     -   69)))
       Fa35(869) = 4*(Fa(770)+Fa(771)+Fa(772)+4*(Fa(773)+Fa(774)-Fa(77
     -   5)))+Fa(776)
       Fa35(870) = mup2mup5*p6mup3*F(9369)+32*(Fa(777)+Fa(778))+4*(Fa(
     -   781)-Fa(782)+Fa(783))+Fa(785)+Fa(786)-16*(Fa(779)-Fa(780)-Fa
     -   (784)-Fa(787))
       Fa35(871) = Fa(796)+Fa(797)+Fa(798)+4*(mup5mup6*p1mup3*F(9469)+
     -   Fa(801)-Fa(802)+Fa(803)-Fa(805)-4*(Fa(800)-Fa(804)-Fa(806))-
     -   2*(mup3mup6*p3mup5*F(5712)-mup3mup5*p3mup6*F(9482)+p2mup5*p6
     -   mup3*p6mup6*F(9616)-Fa(799)-Fa(807)-Fa(808)))
       Fa35(872) = -16*(Fa(809)-Fa(810)-Fa(811)+Fa(815)+Fa(816))+4*(Fa
     -   (812)-Fa(813)+Fa(818))+Fa(819)+8*(mup2mup5*p6mup3*p6mup6*F(9
     -   616)+Fa(814)-Fa(817)-2*Fa(820)+Fa(821)-Fa(822))
       Return
       End
