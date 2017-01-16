c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim6(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,
     -   p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mu
     -   p3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup
     -   3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa6)
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
       Complex*16 Fa6(496:594)
       COMMON/FaPenOddGSim/Fa
       Fa6(496) = mup1mup4*(-2*(p1mup5*F(263)+p2mup5*F(379))+p3mup5*F(
     -   385)+p4mup5*F(485))
       Fa6(497) = -4*p1mup5*p2mup1*F(263)-2*p2mup5*p3mup1*F(396)+mup1m
     -   up5*F(446)+p2mup1*(-4*p2mup5*F(379)+2*(p3mup5*F(385)+p4mup5*
     -   F(485)))
       Fa6(498) = 2*p3mup1*(p1mup5*F(380)+2*(p4mup5*F(356)+p2mup5*F(38
     -   1)+p3mup5*F(382)))+mup1mup5*F(450)
       Fa6(499) = 2*p1mup5*p2mup1*F(443)+4*p2mup1*(p2mup5*F(444)+p3mup
     -   5*F(445)+p4mup5*F(459))-mup1mup5*F(460)
       Fa6(500) = mup1mup5*F(419)+2*(p2mup1*(2*(p1mup5*F(263)+p2mup5*F
     -   (379))+p4mup5*F(383)+p3mup5*F(418))+p2mup5*p3mup1*F(4463))
       Fa6(501) = 2*p1mup5*p3mup1*F(418)-mup1mup5*F(435)+4*p3mup1*(p4m
     -   up5*F(4465)+p3mup5*F(4466)+p2mup5*F(4468))
       Fa6(502) = 2*p1mup5*p2mup1*F(418)+mup1mup5*F(421)+4*p2mup1*(p4m
     -   up5*F(4465)+p3mup5*F(4466)+p2mup5*F(4468))
       Fa6(503) = 2*p1mup5*p2mup1*F(383)+mup1mup5*F(390)+4*p2mup1*(p3m
     -   up5*F(4465)+p4mup5*F(4467)+p2mup5*F(4469))
       Fa6(504) = mup1mup4*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F(3
     -   83)+p3mup5*F(386))
       Fa6(505) = mup1mup5*F(384)+2*(p2mup1*(2*(p1mup5*F(263)+p2mup5*F
     -   (379))+p4mup5*F(383)+p3mup5*F(386))+p2mup5*p3mup1*F(396))
       Fa6(506) = 2*p3mup1*(p1mup5*F(380)+2*(p4mup5*F(356)+p2mup5*F(38
     -   1)+p3mup5*F(382)))-mup1mup5*F(395)
       Fa6(507) = 2*(p1mup4*p2mup1*F(263)+p2mup1*p2mup4*F(379))+p2mup1
     -   *p4mup4*F(383)-p2mup1*p3mup4*F(385)+p2mup4*p3mup1*F(396)+p2m
     -   up4*p4mup1*F(4464)
       Fa6(508) = -(mup1mup4*(2*(p1mup3*F(376)+p2mup3*F(377))+p3mup3*F
     -   (378)))+mup1mup3*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(3
     -   83)+p3mup4*F(386))
       Fa6(509) = mup1mup2*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(3
     -   83)-p3mup4*F(385))
       Fa6(510) = mup1mup4*(p4mup5*F(1363)+p1mup5*F(1365)+p2mup5*F(136
     -   6)+p3mup5*F(1642))
       Fa6(511) = mup1mup4*(p1mup5*F(1649)+2*(p2mup5*F(1654)+p3mup5*F(
     -   1656))+p4mup5*F(1827))
       Fa6(512) = mup1mup4*(p1mup5*F(1390)+2*(p2mup5*F(1395)+p4mup5*F(
     -   1398))+p3mup5*F(1828))
       Fa6(513) = (mup1mup5*F(366))/2.+p2mup1*(p4mup5*F(1363)+p1mup5*F
     -   (1365)+p2mup5*F(1366)+p3mup5*F(1642))+p3mup1*(p4mup5*F(1217)
     -   +p2mup5*F(1772)+2*(p1mup5*F(263)+p3mup5*F(273))*P(7))
       Fa6(514) = (mup1mup5*F(334))/2.+p3mup1*(p1mup5*F(1649)+2*(p2mup
     -   5*F(1654)+p3mup5*F(1656))+p4mup5*F(1827))
       Fa6(515) = (mup1mup5*F(335))/2.+p2mup1*(p1mup5*F(1280)+2*(p2mup
     -   5*F(1285)+p3mup5*F(1288))+p4mup5*F(1825))
       Fa6(516) = p1mup5*p3mup1*F(286)+(mup1mup5*F(370))/2.+p1mup5*p2m
     -   up1*F(1699)+2*(p3mup1*p3mup5*F(885)+p3mup1*p4mup5*F(888)-p2m
     -   up1*p2mup5*F(1703)-p2mup1*p4mup5*F(1706))+p2mup1*p3mup5*F(18
     -   26)+p2mup5*p3mup1*F(2190)
       Fa6(517) = 4*mup2mup3*(p3mup4*p4mup5*F(271)+p3mup5*p4mup4*F(272
     -   )+p4mup4*p4mup5*F(276)+p1mup5*(p4mup4*F(262)+2*(p2mup4*F(255
     -   )+p3mup4*F(257)+p1mup4*F(634)))+p1mup4*(p4mup5*F(260)+2*(p3m
     -   up5*F(266)+p2mup5*F(638)))+p2mup4*p4mup5*F(639)+p2mup5*p4mup
     -   4*F(641)+2*(p2mup4*p3mup5*F(265)+p3mup4*p3mup5*F(270)+p2mup5
     -   *p3mup4*F(637)+p2mup4*p2mup5*F(1264)))+mup2mup3*mup4mup5*F(2
     -   199)
       Fa6(518) = mup2mup3*(4*(p1mup4*(p3mup5*F(272)+p4mup5*F(298)+p2m
     -   up5*F(641)+p1mup5*F(684))+p1mup5*(p3mup4*F(295)+p4mup4*F(313
     -   )+p2mup4*F(1069)))+8*(p2mup4*p3mup5*F(304)+p3mup4*p3mup5*F(3
     -   06)+p3mup5*p4mup4*F(308)+p2mup4*p4mup5*F(314)+p4mup4*p4mup5*
     -   F(315)+p3mup4*p4mup5*F(316)+p2mup5*p4mup4*F(691)+p2mup4*p2mu
     -   p5*F(692)+p2mup5*p3mup4*F(1831))+mup4mup5*F(2329))
       Fa6(519) = p1mup4*p4mup1*F(1389)-2*(p2mup4*p4mup1*F(1394)-p4mup
     -   1*p4mup4*F(1397))+p1mup4*p3mup1*F(2312)+p2mup4*p3mup1*F(2313
     -   )+p3mup1*p3mup4*F(2314)+p3mup1*p4mup4*F(2315)+p3mup4*p4mup1*
     -   F(2322)
       Fa6(520) = p1mup4*p4mup1*F(602)-2*(p4mup1*p4mup4*F(887)-p3mup4*
     -   p4mup1*F(1218))+p1mup4*p2mup1*F(2307)+p2mup1*p2mup4*F(2308)+
     -   p2mup1*p3mup4*F(2309)+p2mup1*p4mup4*F(2310)+p2mup4*p4mup1*F(
     -   2321)
       Fa6(521) = -8*(p1mup5*p2mup1*(p2mup4*F(255)+p1mup4*F(634))+p2mu
     -   p1*(p2mup4*p3mup5*F(265)+p2mup5*p3mup4*F(637)+p1mup4*p2mup5*
     -   F(638)))+mup4mup5*p2mup1*F(1829)+2*mup1mup4*p1mup5*F(263)*P(
     -   7)+mup1mup4*(p4mup5*F(1217)+p2mup5*F(1772)+2*p3mup5*F(273)*P
     -   (7))
       Fa6(522) = p2mup1*(4*(p1mup4*p2mup5*F(641)+p1mup4*p1mup5*F(684)
     -   +p1mup5*p2mup4*F(1069))+mup4mup5*F(1830)+8*(p2mup4*p3mup5*F(
     -   304)+p2mup5*p3mup4*F(1831)))+mup1mup4*(p1mup5*F(603)-2*(p3mu
     -   p5*F(885)+p4mup5*F(888))+p2mup5*F(2189))
       Fa6(523) = mup1mup4*(p1mup3*F(1648)-2*(p3mup3*F(1287)-p2mup3*F(
     -   1653)))+mup1mup3*(p1mup4*F(2307)+p2mup4*F(2308)+p3mup4*F(230
     -   9)+p4mup4*F(2310))
       Fa6(524) = -4*mup1mup3*(p3mup4*p4mup5*F(271)+p3mup5*p4mup4*F(27
     -   2)+p4mup4*p4mup5*F(276)+p1mup5*(p4mup4*F(262)+2*(p2mup4*F(25
     -   5)+p3mup4*F(257)+p1mup4*F(634)))+p1mup4*(p4mup5*F(260)+2*(p3
     -   mup5*F(266)+p2mup5*F(638)))+p2mup4*p4mup5*F(639)+p2mup5*p4mu
     -   p4*F(641)+2*(p2mup4*p3mup5*F(265)+p3mup4*p3mup5*F(270)+p2mup
     -   5*p3mup4*F(637)+p2mup4*p2mup5*F(1264)))+mup1mup3*mup4mup5*F(
     -   1829)
       Fa6(525) = -(mup1mup3*(4*(p1mup5*p4mup4*F(313)-p1mup5*p3mup4*F(
     -   607)+p1mup4*(p3mup5*F(272)+p4mup5*F(298)+p2mup5*F(641)+p1mup
     -   5*F(684)))+p2mup4*(-4*p1mup5*F(681)+8*(p3mup5*F(304)+p4mup5*
     -   F(314)+p2mup5*F(692)))+8*(p3mup4*p3mup5*F(306)+p3mup5*p4mup4
     -   *F(308)+p4mup4*p4mup5*F(315)+p3mup4*p4mup5*F(316)+p2mup5*p4m
     -   up4*F(691)+p2mup5*p3mup4*F(1831))-mup4mup5*F(2200)))
       Fa6(526) = mup1mup2*(p1mup4*F(2312)+p2mup4*F(2313)+p3mup4*F(231
     -   4)+p4mup4*F(2315))
       Fa6(527) = 4*mup1mup2*(p3mup4*p4mup5*F(271)+p3mup5*p4mup4*F(272
     -   )+p4mup4*p4mup5*F(276)+p1mup5*(p4mup4*F(262)+2*(p2mup4*F(255
     -   )+p3mup4*F(257)+p1mup4*F(634)))+p1mup4*(p4mup5*F(260)+2*(p3m
     -   up5*F(266)+p2mup5*F(638)))+p2mup4*p4mup5*F(639)+p2mup5*p4mup
     -   4*F(641)+2*(p2mup4*p3mup5*F(265)+p3mup4*p3mup5*F(270)+p2mup5
     -   *p3mup4*F(637)+p2mup4*p2mup5*F(1264)))+mup1mup2*mup4mup5*F(2
     -   199)
       Fa6(528) = mup1mup2*(4*(p1mup4*(p3mup5*F(272)+p4mup5*F(298)+p2m
     -   up5*F(641)+p1mup5*F(684))+p1mup5*(p3mup4*F(295)+p4mup4*F(313
     -   )+p2mup4*F(1069)))+mup4mup5*F(1830)+8*(p2mup4*p3mup5*F(304)+
     -   p3mup4*p3mup5*F(306)+p3mup5*p4mup4*F(308)+p2mup4*p4mup5*F(31
     -   4)+p4mup4*p4mup5*F(315)+p3mup4*p4mup5*F(316)+p2mup5*p4mup4*F
     -   (691)+p2mup4*p2mup5*F(692)+p2mup5*p3mup4*F(1831)))
       Fa6(529) = 4*mup2mup3*(p1mup5*p3mup4*F(263)-p1mup4*p3mup5*F(346
     -   )+p3mup4*p3mup5*F(356)+p3mup4*p4mup5*F(357)+p3mup5*p4mup4*F(
     -   362)+p2mup5*p3mup4*F(379)-p2mup4*p3mup5*F(2318))+mup2mup3*mu
     -   p4mup5*F(2330)
       Fa6(530) = mup1mup4*(p1mup5*F(631)+p3mup5*F(1247)+p4mup5*F(1248
     -   )+p2mup5*F(1464))+p2mup1*(-4*(p2mup5*p3mup4*F(379)-p2mup4*p3
     -   mup5*F(2318))+mup4mup5*F(2320))
       Fa6(531) = -4*mup1mup3*(p1mup5*p3mup4*F(263)-p1mup4*p3mup5*F(34
     -   6)+p3mup4*p3mup5*F(356)+p3mup4*p4mup5*F(357)+p3mup5*p4mup4*F
     -   (362)+p2mup5*p3mup4*F(379)-p2mup4*p3mup5*F(2318))+mup1mup3*m
     -   up4mup5*F(2320)
       Fa6(532) = 4*mup1mup2*(p1mup5*p3mup4*F(263)-p1mup4*p3mup5*F(346
     -   )+p3mup4*p3mup5*F(356)+p3mup4*p4mup5*F(357)+p3mup5*p4mup4*F(
     -   362)+p2mup5*p3mup4*F(379)-p2mup4*p3mup5*F(2318))+mup1mup2*mu
     -   p4mup5*F(2319)
       Fa6(533) = p3mup1*p4mup4*F(1216)+p3mup4*p4mup1*F(1246)+p4mup1*p
     -   4mup4*F(1249)+p2mup1*p3mup4*F(1357)+2*(p3mup1*p3mup4*F(1214)
     -   +p1mup4*(p3mup1*F(253)+p4mup1*F(630)+p2mup1*F(1360)))+p2mup1
     -   *p4mup4*F(1362)+p2mup4*(p3mup1*F(596)+2*p2mup1*F(1359)+p4mup
     -   1*F(1463))
       Fa6(534) = p2mup1*(-2*(p2mup3*F(1250)+p1mup3*F(1251))+p3mup3*F(
     -   1253))+p2mup3*p3mup1*F(1492)+2*p3mup1*(p1mup3*F(258)+p3mup3*
     -   F(600))*P(7)
       Fa6(535) = p3mup1*(p1mup3*F(1279)+2*(p2mup3*F(1284)+p3mup3*F(12
     -   87)))
       Fa6(536) = p2mup1*(p1mup3*F(1648)-2*(p3mup3*F(1287)-p2mup3*F(16
     -   53)))
       Fa6(537) = p1mup4*p3mup1*F(285)+p1mup4*p2mup1*F(1698)+2*(p3mup1
     -   *p3mup4*F(884)+p3mup1*p4mup4*F(887)+p2mup1*p2mup4*F(1394)+p2
     -   mup1*p4mup4*F(1705))+p2mup1*p3mup4*F(2311)+p2mup4*p3mup1*F(2
     -   316)
       Fa6(538) = p1mup3*p3mup1*F(284)-2*p3mup1*p3mup3*F(290)+p1mup3*p
     -   2mup1*F(2327)+p2mup1*p2mup3*F(2331)+p2mup3*p3mup1*F(2332)+p2
     -   mup1*p3mup3*F(2333)
       Fa6(539) = mup1mup3*(p1mup4*F(1698)+2*(p2mup4*F(1394)+p4mup4*F(
     -   1705))+p3mup4*F(2311))+mup1mup4*(p1mup3*F(2327)+p2mup3*F(233
     -   1)+p3mup3*F(2333))
       Fa6(540) = mup2mup3*(4*(p1mup4*p2mup5*F(346)-p2mup5*p4mup4*F(36
     -   2)-p2mup4*(p1mup5*F(263)+p3mup5*F(273)+p4mup5*F(357)+p2mup5*
     -   F(1470))+p2mup5*p3mup4*F(2188))+mup4mup5*F(2334))
       Fa6(541) = mup1mup3*(p1mup4*F(285)+2*(p3mup4*F(884)+p4mup4*F(88
     -   7))+p2mup4*F(2316))+mup1mup4*(p1mup3*F(284)-2*p3mup3*F(290)+
     -   p2mup3*F(2332))
       Fa6(542) = mup1mup3*(-4*(p1mup4*p2mup5*F(346)-p2mup5*p4mup4*F(3
     -   62)-p2mup4*(p1mup5*F(263)+p3mup5*F(273)+p4mup5*F(357)+p2mup5
     -   *F(1470))+p2mup5*p3mup4*F(2188))+mup4mup5*F(2305))
       Fa6(543) = mup1mup2*(mup4mup5*F(2187)+4*(p1mup4*p2mup5*F(346)-p
     -   2mup5*p4mup4*F(362)-p2mup4*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(1470))+p2mup5*p3mup4*F(2188)))
       Fa6(544) = mup1mup2*(p1mup4*F(602)-2*(p4mup4*F(887)-p3mup4*F(12
     -   18))+p2mup4*F(2321))
       Fa6(545) = mup1mup2*(p1mup4*F(1389)-2*(p2mup4*F(1394)-p4mup4*F(
     -   1397))+p3mup4*F(2322))
       Fa6(546) = mup1mup4*(2*(p2mup3*F(1250)+p1mup3*F(1251))+p3mup3*F
     -   (1641))+mup1mup3*(p3mup4*F(1678)+2*(p2mup4*F(1680)+p1mup4*F(
     -   1681))+p4mup4*F(1683))
       Fa6(547) = mup1mup3*(p2mup4*F(251)-2*(p1mup4*F(253)+p3mup4*F(87
     -   8))+p4mup4*F(880))+mup1mup4*p2mup3*F(1771)+2*mup1mup4*(p1mup
     -   3*F(264)+p3mup3*F(282))*P(7)
       Fa6(548) = mup1mup2*(p2mup4*F(596)+2*(p1mup4*F(253)+p3mup4*F(12
     -   14))+p4mup4*F(1216))
       Fa6(549) = mup1mup2*(p3mup4*F(1357)+2*(p2mup4*F(1359)+p1mup4*F(
     -   1360))+p4mup4*F(1362))
       Fa6(550) = mup1mup4*(p2mup3*F(336)+p1mup3*F(339)+p3mup3*F(340))
     -   +mup1mup3*(p1mup4*F(341)+p3mup4*F(984)+p4mup4*F(988)+p2mup4*
     -   F(1743))
       Fa6(551) = mup1mup2*(2*p1mup4*F(630)+p3mup4*F(1246)+p4mup4*F(12
     -   49)+p2mup4*F(1463))
       Fa6(552) = -4*p4mup1*(p3mup3*p4mup5*F(297)+p3mup5*p4mup3*F(302)
     -   +p1mup3*p3mup5*F(344)+p3mup3*p3mup5*F(351)+p1mup5*(p3mup3*F(
     -   348)+2*(p2mup3*F(255)+p4mup3*F(364)+p1mup3*F(634)))+p2mup3*p
     -   3mup5*F(1266)+p2mup5*p3mup3*F(1267)+2*(p2mup3*p4mup5*F(353)+
     -   p1mup3*p4mup5*F(354)+p4mup3*p4mup5*F(361)+p1mup3*p2mup5*F(63
     -   8)+p2mup3*p2mup5*F(1264)+p2mup5*p4mup3*F(1277)))+mup3mup5*p4
     -   mup1*F(2335)
       Fa6(553) = p4mup1*(4*(p1mup3*(p3mup5*F(293)+p4mup5*F(297)+p1mup
     -   5*F(828)+p2mup5*F(1267))+p1mup5*(p4mup3*F(302)+p3mup3*F(310)
     -   +p2mup3*F(1657)))+8*(p2mup3*p4mup5*F(304)+p3mup5*p4mup3*F(30
     -   6)+p3mup3*p4mup5*F(306)+p2mup3*p3mup5*F(311)+p3mup3*p3mup5*F
     -   (312)+p4mup3*p4mup5*F(316)+p2mup3*p2mup5*F(1303)+p2mup5*p3mu
     -   p3*F(1304)+p2mup5*p4mup3*F(1831))+mup3mup5*F(2336))
       Fa6(554) = p4mup1*(p1mup5*F(342)+p3mup5*F(985)+p4mup5*F(987)+p2
     -   mup5*F(1744))
       Fa6(555) = p4mup1*(p1mup5*F(1699)-2*(p2mup5*F(1703)+p4mup5*F(17
     -   06))+p3mup5*F(1826))
       Fa6(556) = p4mup1*(p1mup5*F(286)+2*(p3mup5*F(885)+p4mup5*F(888)
     -   )+p2mup5*F(2190))
       Fa6(557) = mup2mup4*(mup3mup5*F(2337)+4*(p1mup3*p4mup5*F(258)+p
     -   1mup5*p4mup3*F(263)+p3mup5*p4mup3*F(273)+p4mup3*p4mup5*F(277
     -   )-p3mup3*p4mup5*F(282)+p2mup5*p4mup3*F(379)+p2mup3*p4mup5*F(
     -   2338)))
       Fa6(558) = mup2mup4*(4*(p1mup3*(p3mup5*F(293)+p4mup5*F(297)+p1m
     -   up5*F(828)+p2mup5*F(1267))+p1mup5*(p4mup3*F(302)+p3mup3*F(31
     -   0)+p2mup3*F(1657)))+8*(p2mup3*p4mup5*F(304)+p3mup5*p4mup3*F(
     -   306)+p3mup3*p4mup5*F(306)+p2mup3*p3mup5*F(311)+p3mup3*p3mup5
     -   *F(312)+p4mup3*p4mup5*F(316)+p2mup3*p2mup5*F(1303)+p2mup5*p3
     -   mup3*F(1304)+p2mup5*p4mup3*F(1831))+mup3mup5*F(2339))
       Fa6(559) = -4*mup1mup4*(p1mup3*p4mup5*F(258)+p1mup5*p4mup3*F(26
     -   3)+p3mup5*p4mup3*F(273)+p4mup3*p4mup5*F(277)-p3mup3*p4mup5*F
     -   (282)+p2mup5*p4mup3*F(379)+p2mup3*p4mup5*F(2338))+mup1mup4*m
     -   up3mup5*F(2341)
       Fa6(560) = -4*mup1mup4*(p1mup3*(p3mup5*F(293)+p4mup5*F(297)+p1m
     -   up5*F(828)+p2mup5*F(1267))+p1mup5*(p4mup3*F(302)+p3mup3*F(31
     -   0)-p2mup3*F(1291))+2*(p2mup3*p4mup5*F(304)+p3mup5*p4mup3*F(3
     -   06)+p3mup3*p4mup5*F(306)+p2mup3*p3mup5*F(311)+p3mup3*p3mup5*
     -   F(312)+p4mup3*p4mup5*F(316)+p2mup3*p2mup5*F(1303)+p2mup5*p3m
     -   up3*F(1304)+p2mup5*p4mup3*F(1831)))+mup1mup4*mup3mup5*F(2342
     -   )
       Fa6(561) = mup1mup4*(p1mup2*F(1278)+2*p3mup2*F(1286))
       Fa6(562) = -4*mup2mup4*(p3mup3*p4mup5*F(297)+p3mup5*p4mup3*F(30
     -   2)+p1mup3*p3mup5*F(344)+p3mup3*p3mup5*F(351)+p1mup5*(p3mup3*
     -   F(348)+2*(p2mup3*F(255)+p4mup3*F(364)+p1mup3*F(634)))+p2mup3
     -   *p3mup5*F(1266)+p2mup5*p3mup3*F(1267)+2*(p2mup3*p4mup5*F(353
     -   )+p1mup3*p4mup5*F(354)+p4mup3*p4mup5*F(361)+p1mup3*p2mup5*F(
     -   638)+p2mup3*p2mup5*F(1264)+p2mup5*p4mup3*F(1277)))+mup2mup4*
     -   mup3mup5*F(2335)
       Fa6(563) = 4*mup1mup4*(p3mup3*p4mup5*F(297)+p3mup5*p4mup3*F(302
     -   )+p1mup3*p3mup5*F(344)+p3mup3*p3mup5*F(351)+p1mup5*(p3mup3*F
     -   (348)+2*(p2mup3*F(255)+p4mup3*F(364)+p1mup3*F(634)))+p2mup3*
     -   p3mup5*F(1266)+p2mup5*p3mup3*F(1267)+2*(p2mup3*p4mup5*F(353)
     -   +p1mup3*p4mup5*F(354)+p4mup3*p4mup5*F(361)+p1mup3*p2mup5*F(6
     -   38)+p2mup3*p2mup5*F(1264)+p2mup5*p4mup3*F(1277)))+mup1mup4*m
     -   up3mup5*F(2343)
       Fa6(564) = p4mup1*(p2mup3*F(336)+p1mup3*F(339)+p3mup3*F(340)-2*
     -   p4mup3*F(986))+p2mup1*p4mup3*F(1644)
       Fa6(565) = p4mup1*(p1mup3*F(2327)+p2mup3*F(2331)+p3mup3*F(2333)
     -   +2*p4mup3*F(2344))
       Fa6(566) = p4mup1*(p1mup3*F(284)-2*(p3mup3*F(290)+p4mup3*F(886)
     -   )+p2mup3*F(2332))+p2mup1*p4mup3*F(2340)
       Fa6(567) = p2mup1*p3mup2*F(1640)+p1mup2*p3mup1*F(1773)+2*(p1mup
     -   2*p2mup1*F(1364)+p3mup1*p3mup2*F(1777))
       Fa6(568) = p3mup1*(p1mup2*F(1647)+2*p3mup2*F(1655))
       Fa6(569) = p2mup1*(p1mup2*F(1278)+2*p3mup2*F(1286))
       Fa6(570) = mup3mup5*p4mup1*F(2187)+mup3mup5*p2mup1*F(2341)-4*(-
     -   (p1mup3*(p2mup5*p4mup1-p2mup1*p4mup5)*F(258))+p1mup5*(p2mup3
     -   *p4mup1+p2mup1*p4mup3)*F(263)+p2mup3*p3mup5*p4mup1*F(273)+p2
     -   mup1*p3mup5*p4mup3*F(273)+p2mup1*p4mup3*p4mup5*F(277)+p2mup5
     -   *p3mup3*p4mup1*F(282)-p2mup1*p3mup3*p4mup5*F(282)+p2mup3*p4m
     -   up1*p4mup5*F(357)+p2mup1*p2mup5*p4mup3*F(379)+p2mup3*p2mup5*
     -   p4mup1*F(1500)+p2mup1*p2mup3*p4mup5*F(2338)+p2mup5*p4mup1*p4
     -   mup3*F(2345))
       Fa6(571) = 4*p2mup1*(p3mup3*p4mup5*F(297)+p3mup5*p4mup3*F(302)+
     -   p1mup3*p3mup5*F(344)+p3mup3*p3mup5*F(351)+p1mup5*(p3mup3*F(3
     -   48)+2*(p2mup3*F(255)+p4mup3*F(364)+p1mup3*F(634)))+p2mup3*p3
     -   mup5*F(1266)+p2mup5*p3mup3*F(1267)+2*(p2mup3*p4mup5*F(353)+p
     -   1mup3*p4mup5*F(354)+p4mup3*p4mup5*F(361)+p1mup3*p2mup5*F(638
     -   )+p2mup3*p2mup5*F(1264)+p2mup5*p4mup3*F(1277)))+mup3mup5*p2m
     -   up1*F(2343)
       Fa6(572) = p2mup1*(4*(p1mup3*(p3mup5*F(293)+p4mup5*F(297)+p1mup
     -   5*F(828)+p2mup5*F(1267))+p1mup5*(p4mup3*F(302)+p3mup3*F(310)
     -   +p2mup3*F(1657)))+8*(p2mup3*p4mup5*F(304)+p3mup5*p4mup3*F(30
     -   6)+p3mup3*p4mup5*F(306)+p2mup3*p3mup5*F(311)+p3mup3*p3mup5*F
     -   (312)+p4mup3*p4mup5*F(316)+p2mup3*p2mup5*F(1303)+p2mup5*p3mu
     -   p3*F(1304)+p2mup5*p4mup3*F(1831))+mup3mup5*F(2336))
       Fa6(573) = p1mup2*p2mup1*F(1697)+p1mup2*p3mup1*F(2346)+p2mup1*p
     -   3mup2*F(2347)+p3mup1*p3mup2*F(2348)
       Fa6(574) = -2*mup1mup4*p4mup3*F(2344)+mup3mup4*p3mup1*F(2349)
       Fa6(575) = mup1mup3*p1mup2*F(1388)+mup1mup2*p2mup3*F(2331)+mup1
     -   mup2*p3mup3*F(2333)+2*mup1mup2*p4mup3*F(2344)+mup1mup3*p3mup
     -   2*F(2350)+mup2mup3*p3mup1*F(2351)
       Fa6(576) = mup1mup2*(p3mup3*p3mup5*F(351)+2*(p4mup3*p4mup5*F(36
     -   1)+p1mup3*p1mup5*F(634))+p2mup3*p3mup5*F(1266)+p2mup5*p3mup3
     -   *F(1267))
       Fa6(577) = 4*mup1mup2*(p1mup3*p1mup5*F(828)+2*(p2mup3*p3mup5*F(
     -   311)+p4mup3*p4mup5*F(316)+p3mup3*(p3mup5*F(312)+p2mup5*F(130
     -   4))))
       Fa6(578) = mup1mup3*p1mup2*F(1647)-2*(mup1mup2*p3mup3*F(1287)-m
     -   up1mup2*p2mup3*F(1653)-mup1mup3*p3mup2*F(1655))+mup1mup2*p4m
     -   up3*F(2352)
       Fa6(579) = mup2mup4*(4*(p1mup3*p2mup5*F(258)-p2mup5*p3mup3*F(28
     -   2)-p2mup3*(p1mup5*F(263)+p3mup5*F(273)+p4mup5*F(357)+p2mup5*
     -   F(1500))-p2mup5*p4mup3*F(2345))+mup3mup5*F(2353))
       Fa6(580) = 2*mup1mup4*p4mup3*F(886)+mup3mup4*p2mup1*F(2354)
       Fa6(581) = mup1mup4*(mup3mup5*F(2305)+4*(p2mup3*(p1mup5*F(263)+
     -   p3mup5*F(273)+p4mup5*F(357))+p2mup5*(-(p1mup3*F(258))+p3mup3
     -   *F(282)+p2mup3*F(1500)+p4mup3*F(2345))))
       Fa6(582) = mup1mup2*(-2*(p3mup3*F(290)+p4mup3*F(886))+p2mup3*F(
     -   2332))+mup1mup3*p1mup2*F(2355)+mup2mup3*p2mup1*F(2356)+mup1m
     -   up3*p3mup2*F(2357)
       Fa6(583) = mup1mup2*(p2mup3*p3mup5*F(273)+p2mup5*p3mup3*F(282))
       Fa6(584) = mup1mup4*p1mup2*F(2346)+mup1mup4*p3mup2*F(2348)+mup2
     -   mup4*p2mup1*F(2358)
       Fa6(585) = mup1mup4*p1mup2*F(1697)+mup1mup4*p3mup2*F(2347)+mup2
     -   mup4*p3mup1*F(2359)
       Fa6(586) = mup1mup4*p4mup3*F(1256)+mup3mup4*p3mup1*F(2360)
       Fa6(587) = 2*mup1mup3*p1mup2*F(1364)+mup1mup3*p3mup2*F(1640)+mu
     -   p1mup2*(2*p2mup3*F(1250)+p3mup3*F(1641)+p4mup3*F(1644))+mup2
     -   mup3*p3mup1*F(2341)
       Fa6(588) = mup1mup4*p4mup3*F(1215)+mup3mup4*p2mup1*F(2361)
       Fa6(589) = mup1mup3*p1mup2*F(1773)+2*mup1mup3*p3mup2*F(1777)+mu
     -   p2mup3*p2mup1*F(2362)+mup1mup2*(p4mup3*F(879)+p2mup3*F(1771)
     -   +2*p3mup3*F(282)*P(7))
       Fa6(590) = 2*mup1mup4*(p1mup2*F(1494)+p3mup2*F(1498))+mup2mup4*
     -   p2mup1*F(2363)
       Fa6(591) = mup1mup4*p3mup2*F(1252)+2*mup1mup4*p1mup2*F(1258)+mu
     -   p2mup4*p3mup1*F(2337)
       Fa6(592) = 2*mup1mup4*p4mup3*F(986)+mup3mup4*p3mup1*F(2187)+mup
     -   3mup4*p2mup1*F(2364)
       Fa6(593) = 2*mup1mup4*p1mup2*F(1745)+mup1mup4*p3mup2*F(1746)+mu
     -   p2mup4*p2mup1*F(2320)+mup2mup4*p3mup1*F(2365)
       Fa6(594) = mup1mup2*(p2mup3*F(336)+p3mup3*F(340)-2*p4mup3*F(986
     -   ))+2*mup1mup3*p1mup2*F(1465)+mup1mup3*p3mup2*F(1466)+mup2mup
     -   3*p2mup1*F(2330)+mup2mup3*p3mup1*F(2334)
       Return
       End
