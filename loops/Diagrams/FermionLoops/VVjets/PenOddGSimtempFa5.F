c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim5(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,
     -   p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mu
     -   p3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup
     -   3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa5)
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
       Complex*16 Fa5(397:495)
       COMMON/FaPenOddGSim/Fa
       Fa5(397) = mup1mup3*(p1mup5*F(286)+2*(p3mup5*F(885)+p4mup5*F(88
     -   8))+p2mup5*F(2190))
       Fa5(398) = mup1mup2*(p1mup5*F(631)+p3mup5*F(1247)+p4mup5*F(1248
     -   )+p2mup5*F(1464))
       Fa5(399) = mup1mup2*(p1mup5*F(603)-2*(p3mup5*F(885)+p4mup5*F(88
     -   8))+p2mup5*F(2189))
       Fa5(400) = mup1mup4*(p3mup5*F(1254)+2*(p1mup5*F(1259)+p2mup5*F(
     -   1260))+p4mup5*F(1684))
       Fa5(401) = mup1mup4*(p1mup5*F(1280)+2*(p2mup5*F(1285)+p3mup5*F(
     -   1288))+p4mup5*F(1825))
       Fa5(402) = mup1mup4*(p1mup5*F(1699)-2*(p2mup5*F(1703)+p4mup5*F(
     -   1706))+p3mup5*F(1826))
       Fa5(403) = (mup1mup5*F(367))/2.+p2mup1*(p3mup5*F(1254)+2*(p1mup
     -   5*F(1259)+p2mup5*F(1260))+p4mup5*F(1684))+p3mup1*(p4mup5*F(8
     -   81)+p2mup5*F(1493)-2*(p1mup5*F(263)+p3mup5*F(273))*P(7))
       Fa5(404) = (mup1mup5*F(368))/2.+p3mup1*(p1mup5*F(1280)+2*(p2mup
     -   5*F(1285)+p3mup5*F(1288))+p4mup5*F(1825))
       Fa5(405) = (mup1mup5*F(369))/2.+p2mup1*(p1mup5*F(1649)+2*(p2mup
     -   5*F(1654)+p3mup5*F(1656))+p4mup5*F(1827))
       Fa5(406) = (mup1mup5*F(371))/2.+p1mup5*p3mup1*F(603)+p1mup5*p2m
     -   up1*F(1390)-2*(p3mup1*p3mup5*F(885)+p3mup1*p4mup5*F(888)-p2m
     -   up1*p2mup5*F(1395)-p2mup1*p4mup5*F(1398))+p2mup1*p3mup5*F(18
     -   28)+p2mup5*p3mup1*F(2189)
       Fa5(407) = -4*mup2mup3*(p3mup4*p4mup5*F(271)+p3mup5*p4mup4*F(27
     -   2)+p4mup4*p4mup5*F(276)+p1mup5*(p4mup4*F(262)+2*(p2mup4*F(25
     -   5)+p3mup4*F(257)+p1mup4*F(634)))+p1mup4*(p4mup5*F(260)+2*(p3
     -   mup5*F(266)+p2mup5*F(638)))+p2mup4*p4mup5*F(639)+p2mup5*p4mu
     -   p4*F(641)+2*(p2mup4*p3mup5*F(265)+p3mup4*p3mup5*F(270)+p2mup
     -   5*p3mup4*F(637)+p2mup4*p2mup5*F(1264)))+mup2mup3*mup4mup5*F(
     -   1829)
       Fa5(408) = -(mup2mup3*(4*(p1mup5*p4mup4*F(313)-p1mup5*p3mup4*F(
     -   607)+p1mup4*(p3mup5*F(272)+p4mup5*F(298)+p2mup5*F(641)+p1mup
     -   5*F(684)))+p2mup4*(-4*p1mup5*F(681)+8*(p3mup5*F(304)+p4mup5*
     -   F(314)+p2mup5*F(692)))+8*(p3mup4*p3mup5*F(306)+p3mup5*p4mup4
     -   *F(308)+p4mup4*p4mup5*F(315)+p3mup4*p4mup5*F(316)+p2mup5*p4m
     -   up4*F(691)+p2mup5*p3mup4*F(1831))-mup4mup5*F(2306)))
       Fa5(409) = p1mup4*p4mup1*F(1698)+2*(p2mup4*p4mup1*F(1394)+p4mup
     -   1*p4mup4*F(1705))+p1mup4*p3mup1*F(2307)+p2mup4*p3mup1*F(2308
     -   )+p3mup1*p3mup4*F(2309)+p3mup1*p4mup4*F(2310)+p3mup4*p4mup1*
     -   F(2311)
       Fa5(410) = p1mup4*p4mup1*F(285)+2*(p3mup4*p4mup1*F(884)+p4mup1*
     -   p4mup4*F(887))+p1mup4*p2mup1*F(2312)+p2mup1*p2mup4*F(2313)+p
     -   2mup1*p3mup4*F(2314)+p2mup1*p4mup4*F(2315)+p2mup4*p4mup1*F(2
     -   316)
       Fa5(411) = 8*(p1mup5*p2mup1*p2mup4*F(255)+p2mup1*p2mup4*p3mup5*
     -   F(265)+p1mup4*p1mup5*p2mup1*F(634)+p2mup1*p2mup5*p3mup4*F(63
     -   7)+p1mup4*p2mup1*p2mup5*F(638))+mup1mup4*p4mup5*F(881)+mup1m
     -   up4*p2mup5*F(1493)+mup4mup5*p2mup1*F(2199)-2*mup1mup4*(p1mup
     -   5*F(263)+p3mup5*F(273))*P(7)
       Fa5(412) = -4*p2mup1*(p1mup4*p2mup5*F(641)-p1mup5*p2mup4*F(681)
     -   +p1mup4*p1mup5*F(684)+2*(p2mup4*p3mup5*F(304)+p2mup5*p3mup4*
     -   F(1831)))+mup1mup4*(p1mup5*F(286)+2*(p3mup5*F(885)+p4mup5*F(
     -   888))+p2mup5*F(2190))+mup4mup5*p2mup1*F(2200)
       Fa5(413) = mup1mup4*(p1mup3*F(1279)+2*(p2mup3*F(1284)+p3mup3*F(
     -   1287)))+mup1mup3*(p1mup4*F(2312)+p2mup4*F(2313)+p3mup4*F(231
     -   4)+p4mup4*F(2315))
       Fa5(414) = 4*mup1mup3*(p3mup4*p4mup5*F(271)+p3mup5*p4mup4*F(272
     -   )+p4mup4*p4mup5*F(276)+p1mup5*(p4mup4*F(262)+2*(p2mup4*F(255
     -   )+p3mup4*F(257)+p1mup4*F(634)))+p1mup4*(p4mup5*F(260)+2*(p3m
     -   up5*F(266)+p2mup5*F(638)))+p2mup4*p4mup5*F(639)+p2mup5*p4mup
     -   4*F(641)+2*(p2mup4*p3mup5*F(265)+p3mup4*p3mup5*F(270)+p2mup5
     -   *p3mup4*F(637)+p2mup4*p2mup5*F(1264)))+mup1mup3*mup4mup5*F(2
     -   199)
       Fa5(415) = mup1mup3*(4*(p1mup4*(p3mup5*F(272)+p4mup5*F(298)+p2m
     -   up5*F(641)+p1mup5*F(684))+p1mup5*(p3mup4*F(295)+p4mup4*F(313
     -   )+p2mup4*F(1069)))+mup4mup5*F(1830)+8*(p2mup4*p3mup5*F(304)+
     -   p3mup4*p3mup5*F(306)+p3mup5*p4mup4*F(308)+p2mup4*p4mup5*F(31
     -   4)+p4mup4*p4mup5*F(315)+p3mup4*p4mup5*F(316)+p2mup5*p4mup4*F
     -   (691)+p2mup4*p2mup5*F(692)+p2mup5*p3mup4*F(1831)))
       Fa5(416) = mup1mup2*(p1mup4*F(2307)+p2mup4*F(2308)+p3mup4*F(230
     -   9)+p4mup4*F(2310))
       Fa5(417) = -4*mup1mup2*(p3mup4*p4mup5*F(271)+p3mup5*p4mup4*F(27
     -   2)+p4mup4*p4mup5*F(276)+p1mup5*(p4mup4*F(262)+2*(p2mup4*F(25
     -   5)+p3mup4*F(257)+p1mup4*F(634)))+p1mup4*(p4mup5*F(260)+2*(p3
     -   mup5*F(266)+p2mup5*F(638)))+p2mup4*p4mup5*F(639)+p2mup5*p4mu
     -   p4*F(641)+2*(p2mup4*p3mup5*F(265)+p3mup4*p3mup5*F(270)+p2mup
     -   5*p3mup4*F(637)+p2mup4*p2mup5*F(1264)))+mup1mup2*mup4mup5*F(
     -   1829)
       Fa5(418) = -(mup1mup2*(4*(p1mup5*p4mup4*F(313)-p1mup5*p3mup4*F(
     -   607)+p1mup4*(p3mup5*F(272)+p4mup5*F(298)+p2mup5*F(641)+p1mup
     -   5*F(684)))+p2mup4*(-4*p1mup5*F(681)+8*(p3mup5*F(304)+p4mup5*
     -   F(314)+p2mup5*F(692)))+8*(p3mup4*p3mup5*F(306)+p3mup5*p4mup4
     -   *F(308)+p4mup4*p4mup5*F(315)+p3mup4*p4mup5*F(316)+p2mup5*p4m
     -   up4*F(691)+p2mup5*p3mup4*F(1831))-mup4mup5*F(2200)))
       Fa5(419) = mup2mup3*(mup4mup5*F(2317)+4*(-(p1mup5*p3mup4*F(263)
     -   )+p1mup4*p3mup5*F(346)-p3mup4*p3mup5*F(356)-p3mup4*p4mup5*F(
     -   357)-p3mup5*p4mup4*F(362)-p2mup5*p3mup4*F(379)+p2mup4*p3mup5
     -   *F(2318)))
       Fa5(420) = mup1mup4*(p1mup5*F(342)+p3mup5*F(985)+p4mup5*F(987)+
     -   p2mup5*F(1744))+p2mup1*(4*(p2mup5*p3mup4*F(379)-p2mup4*p3mup
     -   5*F(2318))+mup4mup5*F(2319))
       Fa5(421) = 4*mup1mup3*(p1mup5*p3mup4*F(263)-p1mup4*p3mup5*F(346
     -   )+p3mup4*p3mup5*F(356)+p3mup4*p4mup5*F(357)+p3mup5*p4mup4*F(
     -   362)+p2mup5*p3mup4*F(379)-p2mup4*p3mup5*F(2318))+mup1mup3*mu
     -   p4mup5*F(2319)
       Fa5(422) = -4*mup1mup2*(p1mup5*p3mup4*F(263)-p1mup4*p3mup5*F(34
     -   6)+p3mup4*p3mup5*F(356)+p3mup4*p4mup5*F(357)+p3mup5*p4mup4*F
     -   (362)+p2mup5*p3mup4*F(379)-p2mup4*p3mup5*F(2318))+mup1mup2*m
     -   up4mup5*F(2320)
       Fa5(423) = -2*p3mup1*p3mup4*F(878)+p3mup1*p4mup4*F(880)+p3mup4*
     -   p4mup1*F(984)+p4mup1*p4mup4*F(988)+p2mup1*p3mup4*F(1678)+p1m
     -   up4*(p4mup1*F(341)-2*(p3mup1*F(253)-p2mup1*F(1681)))+p2mup1*
     -   p4mup4*F(1683)+p2mup4*(p3mup1*F(251)+2*p2mup1*F(1680)+p4mup1
     -   *F(1743))
       Fa5(424) = p2mup1*(2*(p2mup3*F(1250)+p1mup3*F(1251))+p3mup3*F(1
     -   641))+p2mup3*p3mup1*F(1771)+2*p3mup1*(p1mup3*F(264)+p3mup3*F
     -   (282))*P(7)
       Fa5(425) = p3mup1*(p1mup3*F(1648)-2*(p3mup3*F(1287)-p2mup3*F(16
     -   53)))
       Fa5(426) = p2mup1*(p1mup3*F(1279)+2*(p2mup3*F(1284)+p3mup3*F(12
     -   87)))
       Fa5(427) = p1mup4*p3mup1*F(602)+p1mup4*p2mup1*F(1389)-2*(p3mup1
     -   *p4mup4*F(887)-p3mup1*p3mup4*F(1218)+p2mup1*p2mup4*F(1394)-p
     -   2mup1*p4mup4*F(1397))+p2mup4*p3mup1*F(2321)+p2mup1*p3mup4*F(
     -   2322)
       Fa5(428) = p3mup1*(2*p3mup3*F(290)+p1mup3*F(601)+p2mup3*F(2325)
     -   )+p2mup1*(p1mup3*F(2323)+p2mup3*F(2324)+p3mup3*F(2326))
       Fa5(429) = mup1mup3*(p1mup4*F(1389)-2*(p2mup4*F(1394)-p4mup4*F(
     -   1397))+p3mup4*F(2322))+mup1mup4*(p1mup3*F(2323)+p2mup3*F(232
     -   4)+p3mup3*F(2326))
       Fa5(430) = mup2mup3*(-4*(p1mup4*p2mup5*F(346)-p2mup5*p4mup4*F(3
     -   62)-p2mup4*(p1mup5*F(263)+p3mup5*F(273)+p4mup5*F(357)+p2mup5
     -   *F(1470))+p2mup5*p3mup4*F(2188))+mup4mup5*F(2328))
       Fa5(431) = mup1mup3*(p1mup4*F(602)-2*(p4mup4*F(887)-p3mup4*F(12
     -   18))+p2mup4*F(2321))+mup1mup4*(2*p3mup3*F(290)+p1mup3*F(601)
     -   +p2mup3*F(2325))
       Fa5(432) = mup1mup3*(mup4mup5*F(2187)+4*(p1mup4*p2mup5*F(346)-p
     -   2mup5*p4mup4*F(362)-p2mup4*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(1470))+p2mup5*p3mup4*F(2188)))
       Fa5(433) = mup1mup2*(-4*(p1mup4*p2mup5*F(346)-p2mup5*p4mup4*F(3
     -   62)-p2mup4*(p1mup5*F(263)+p3mup5*F(273)+p4mup5*F(357)+p2mup5
     -   *F(1470))+p2mup5*p3mup4*F(2188))+mup4mup5*F(2305))
       Fa5(434) = mup1mup2*(p1mup4*F(285)+2*(p3mup4*F(884)+p4mup4*F(88
     -   7))+p2mup4*F(2316))
       Fa5(435) = mup1mup2*(p1mup4*F(1698)+2*(p2mup4*F(1394)+p4mup4*F(
     -   1705))+p3mup4*F(2311))
       Fa5(436) = mup1mup4*(-2*(p2mup3*F(1250)+p1mup3*F(1251))+p3mup3*
     -   F(1253))+mup1mup3*(p3mup4*F(1357)+2*(p2mup4*F(1359)+p1mup4*F
     -   (1360))+p4mup4*F(1362))
       Fa5(437) = mup1mup3*(p2mup4*F(596)+2*(p1mup4*F(253)+p3mup4*F(12
     -   14))+p4mup4*F(1216))+mup1mup4*p2mup3*F(1492)+2*mup1mup4*(p1m
     -   up3*F(258)+p3mup3*F(600))*P(7)
       Fa5(438) = mup1mup2*(p2mup4*F(251)-2*(p1mup4*F(253)+p3mup4*F(87
     -   8))+p4mup4*F(880))
       Fa5(439) = mup1mup2*(p3mup4*F(1678)+2*(p2mup4*F(1680)+p1mup4*F(
     -   1681))+p4mup4*F(1683))
       Fa5(440) = mup1mup4*(p2mup3*F(625)+2*p1mup3*F(628)+p3mup3*F(629
     -   ))+mup1mup3*(2*p1mup4*F(630)+p3mup4*F(1246)+p4mup4*F(1249)+p
     -   2mup4*F(1463))
       Fa5(441) = mup1mup2*(p1mup4*F(341)+p3mup4*F(984)+p4mup4*F(988)+
     -   p2mup4*F(1743))
       Fa5(442) = mup1mup4*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F(3
     -   83)-p3mup5*F(385))
       Fa5(443) = mup1mup4*(p1mup5*F(380)+2*(p4mup5*F(356)+p2mup5*F(38
     -   1)+p3mup5*F(382)))
       Fa5(444) = mup1mup4*(p1mup5*F(383)+2*(p3mup5*F(4465)+p4mup5*F(4
     -   467)+p2mup5*F(4469)))
       Fa5(445) = -2*(p2mup1*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F
     -   (383)-p3mup5*F(385))+p2mup5*p3mup1*F(396))+mup1mup5*F(410)
       Fa5(446) = 2*p3mup1*(p1mup5*F(380)+2*(p4mup5*F(356)+p2mup5*F(38
     -   1)+p3mup5*F(382)))+mup1mup5*F(411)
       Fa5(447) = 2*p2mup1*(p1mup5*F(380)+2*(p4mup5*F(356)+p2mup5*F(38
     -   1)+p3mup5*F(382)))-mup1mup5*F(412)
       Fa5(448) = 2*p1mup5*p2mup1*F(383)-mup1mup5*F(413)+4*p2mup1*(p3m
     -   up5*F(4465)+p4mup5*F(4467)+p2mup5*F(4469))
       Fa5(449) = p3mup1*(p1mup4*F(380)+2*(p4mup4*F(356)+p2mup4*F(381)
     -   +p3mup4*F(382)))+p4mup1*(p1mup4*F(383)+2*(p3mup4*F(4465)+p4m
     -   up4*F(4467)+p2mup4*F(4469)))
       Fa5(450) = p2mup1*(p1mup4*F(380)+2*(p4mup4*F(356)+p2mup4*F(381)
     -   +p3mup4*F(382)))
       Fa5(451) = mup4mup5*p2mup1-mup1mup4*p2mup5
       Fa5(452) = -(mup1mup4*(p1mup3*F(380)+2*(p2mup3*F(381)+p3mup3*F(
     -   382))))+mup1mup3*(p1mup4*F(380)+2*(p4mup4*F(356)+p2mup4*F(38
     -   1)+p3mup4*F(382)))
       Fa5(453) = mup1mup2*(p1mup4*F(380)+2*(p4mup4*F(356)+p2mup4*F(38
     -   1)+p3mup4*F(382)))
       Fa5(454) = 2*(p1mup4*p2mup1*F(263)+p2mup1*p2mup4*F(379))+p2mup1
     -   *p4mup4*F(383)+p2mup1*p3mup4*F(386)+p2mup4*p3mup1*F(396)+p2m
     -   up4*p4mup1*F(4464)
       Fa5(455) = 2*(p1mup3*p2mup1*F(376)+p2mup1*p2mup3*F(377))+p2mup1
     -   *p3mup3*F(378)+p2mup3*p3mup1*F(396)
       Fa5(456) = p3mup1*(p1mup3*F(380)+2*(p2mup3*F(381)+p3mup3*F(382)
     -   ))
       Fa5(457) = p2mup1*(p1mup3*F(380)+2*(p2mup3*F(381)+p3mup3*F(382)
     -   ))
       Fa5(458) = p2mup1*(p1mup4*F(383)+2*(p3mup4*F(4465)+p4mup4*F(446
     -   7)+p2mup4*F(4469)))
       Fa5(459) = p2mup1*(p1mup3*F(383)+2*(p3mup3*F(4465)+p2mup3*F(446
     -   9)))
       Fa5(460) = -(mup1mup4*(p1mup3*F(383)+2*(p3mup3*F(4465)+p2mup3*F
     -   (4469))))+mup1mup3*(p1mup4*F(383)+2*(p3mup4*F(4465)+p4mup4*F
     -   (4467)+p2mup4*F(4469)))
       Fa5(461) = mup1mup2*(p1mup4*F(383)+2*(p3mup4*F(4465)+p4mup4*F(4
     -   467)+p2mup4*F(4469)))
       Fa5(462) = mup1mup4*(2*(p1mup3*F(376)+p2mup3*F(377))+p3mup3*F(3
     -   78))-mup1mup3*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(383)
     -   -p3mup4*F(385))
       Fa5(463) = mup1mup4*p2mup3-mup1mup3*p2mup4
       Fa5(464) = mup1mup2*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(3
     -   83)+p3mup4*F(386))
       Fa5(465) = mup1mup4*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F(3
     -   83)+p3mup5*F(418))
       Fa5(466) = mup1mup4*(p1mup5*F(418)+2*(p4mup5*F(4465)+p3mup5*F(4
     -   466)+p2mup5*F(4468)))
       Fa5(467) = -(mup1mup5*F(436))+2*(p2mup1*(2*(p1mup5*F(263)+p2mup
     -   5*F(379))+p4mup5*F(383)+p3mup5*F(418))+p2mup5*p3mup1*F(4463)
     -   )
       Fa5(468) = 2*p1mup5*p3mup1*F(418)+mup1mup5*F(437)+4*p3mup1*(p4m
     -   up5*F(4465)+p3mup5*F(4466)+p2mup5*F(4468))
       Fa5(469) = 2*p1mup5*p2mup1*F(418)-mup1mup5*F(438)+4*p2mup1*(p4m
     -   up5*F(4465)+p3mup5*F(4466)+p2mup5*F(4468))
       Fa5(470) = p1mup4*(p4mup1*F(383)+p3mup1*F(418))+2*(p3mup4*p4mup
     -   1*F(4465)+p3mup1*p4mup4*F(4465)+p3mup1*p3mup4*F(4466)+p4mup1
     -   *p4mup4*F(4467)+p2mup4*p3mup1*F(4468)+p2mup4*p4mup1*F(4469))
       Fa5(471) = p2mup1*(p1mup4*F(418)+2*(p4mup4*F(4465)+p3mup4*F(446
     -   6)+p2mup4*F(4468)))
       Fa5(472) = -(mup1mup4*(p1mup3*F(418)+2*(p3mup3*F(4466)+p2mup3*F
     -   (4468))))+mup1mup3*(p1mup4*F(418)+2*(p4mup4*F(4465)+p3mup4*F
     -   (4466)+p2mup4*F(4468)))
       Fa5(473) = mup1mup2*(p1mup4*F(418)+2*(p4mup4*F(4465)+p3mup4*F(4
     -   466)+p2mup4*F(4468)))
       Fa5(474) = 2*(p1mup4*p2mup1*F(263)+p2mup1*p2mup4*F(379))+p2mup1
     -   *p4mup4*F(383)+p2mup1*p3mup4*F(418)+p2mup4*p3mup1*F(4463)+p2
     -   mup4*p4mup1*F(4464)
       Fa5(475) = 2*(p1mup3*p2mup1*F(263)+p2mup1*p2mup3*F(379))+p2mup1
     -   *p3mup3*F(418)+p2mup3*p3mup1*F(4463)
       Fa5(476) = p3mup1*(p1mup3*F(418)+2*(p3mup3*F(4466)+p2mup3*F(446
     -   8)))
       Fa5(477) = p2mup1*(p1mup3*F(418)+2*(p3mup3*F(4466)+p2mup3*F(446
     -   8)))
       Fa5(478) = -(mup1mup4*(2*(p1mup3*F(263)+p2mup3*F(379))+p3mup3*F
     -   (418)))+mup1mup3*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(3
     -   83)+p3mup4*F(418))
       Fa5(479) = mup1mup2*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(3
     -   83)+p3mup4*F(418))
       Fa5(480) = mup1mup4*(2*(p1mup5*F(263)+p2mup5*F(379))+p3mup5*F(3
     -   86)+p4mup5*F(447))
       Fa5(481) = mup1mup4*(p1mup5*F(443)+2*(p2mup5*F(444)+p3mup5*F(44
     -   5)+p4mup5*F(459)))
       Fa5(482) = 4*p1mup5*p2mup1*F(263)+2*(p2mup5*p3mup1*F(396)+p2mup
     -   1*(2*p2mup5*F(379)+p3mup5*F(386)+p4mup5*F(447)))+mup1mup5*F(
     -   478)
       Fa5(483) = 2*p3mup1*(p1mup5*F(380)+2*(p4mup5*F(356)+p2mup5*F(38
     -   1)+p3mup5*F(382)))-mup1mup5*F(479)
       Fa5(484) = 2*p2mup1*(p1mup5*F(380)+2*(p4mup5*F(356)+p2mup5*F(38
     -   1)+p3mup5*F(382)))+mup1mup5*F(389)
       Fa5(485) = 2*p1mup5*p2mup1*F(443)+4*p2mup1*(p2mup5*F(444)+p3mup
     -   5*F(445)+p4mup5*F(459))+mup1mup5*F(480)
       Fa5(486) = p3mup1*(p1mup4*F(380)+2*(p4mup4*F(356)+p2mup4*F(381)
     -   +p3mup4*F(382)))+p4mup1*(p1mup4*F(443)+2*(p2mup4*F(444)+p3mu
     -   p4*F(445)+p4mup4*F(459)))
       Fa5(487) = -(mup4mup5*p2mup1)+mup1mup4*p2mup5
       Fa5(488) = 2*(p1mup4*p2mup1*F(376)+p2mup1*p2mup4*F(377))+p2mup1
     -   *p3mup4*F(378)+p2mup4*p3mup1*F(396)+p2mup1*p4mup4*F(448)+p2m
     -   up4*p4mup1*F(477)
       Fa5(489) = p2mup1*(p1mup4*F(443)+2*(p2mup4*F(444)+p3mup4*F(445)
     -   +p4mup4*F(459)))
       Fa5(490) = p2mup1*(p1mup3*F(443)+2*(p2mup3*F(444)+p3mup3*F(445)
     -   ))
       Fa5(491) = -(mup1mup4*(p1mup3*F(443)+2*(p2mup3*F(444)+p3mup3*F(
     -   445))))+mup1mup3*(p1mup4*F(443)+2*(p2mup4*F(444)+p3mup4*F(44
     -   5)+p4mup4*F(459)))
       Fa5(492) = mup1mup2*(p1mup4*F(443)+2*(p2mup4*F(444)+p3mup4*F(44
     -   5)+p4mup4*F(459)))
       Fa5(493) = -(mup1mup4*(2*(p1mup3*F(376)+p2mup3*F(377))+p3mup3*F
     -   (378)))+mup1mup3*(2*(p1mup4*F(376)+p2mup4*F(377))+p3mup4*F(3
     -   78)+p4mup4*F(448))
       Fa5(494) = -(mup1mup4*p2mup3)+mup1mup3*p2mup4
       Fa5(495) = mup1mup2*(2*(p1mup4*F(376)+p2mup4*F(377))+p3mup4*F(3
     -   78)+p4mup4*F(448))
       Return
       End
