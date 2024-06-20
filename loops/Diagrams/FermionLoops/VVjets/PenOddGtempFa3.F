c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddG3(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1
     -   ,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3m
     -   up5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,
     -   p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,m
     -   up2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa3)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenOddG/F
      Complex*16 EpsF(126)
       common/EPSFPenOddG/EPSF
       COMMON/InvariantsPenOddG/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s
     -   45,s234,m1,m1sq,m3,m3sq,m4,m4sq
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
      Complex*16 Fa(4390),F(7459)
      Real*8 P(1989) 
        COMMON/PPenOddG/P
       Complex*16 Fa3(585:876)
       COMMON/FaPenOddG/Fa
       Fa3(585) = p1mup1*(p2mup4*F(878)-2*(p1mup4*F(880)+p3mup4*F(881)
     -   )+p4mup4*F(883))+p4mup1*(p1mup4*F(1237)-2*(p4mup4*F(906)-p3m
     -   up4*F(1242))+p2mup4*F(2410))
       Fa3(586) = p1mup1*(p3mup4*F(1005)+p1mup4*F(1007)+p4mup4*F(1011)
     -   +p2mup4*F(1770))
       Fa3(587) = mup1mup4*(p3mup5*F(1273)+p1mup5*F(1275)+p4mup5*F(127
     -   6)+p2mup5*F(1512))
       Fa3(588) = mup1mup4*(-2*(p3mup5*F(904)+p4mup5*F(907))+p1mup5*F(
     -   1238)+p2mup5*F(2401))
       Fa3(589) = p2mup1*(-2*(p2mup3*F(1830)+p1mup3*F(1831))+p3mup3*F(
     -   1833))+p2mup3*p3mup1*F(1929)+2*p3mup1*(p1mup3*F(271)+p3mup3*
     -   F(1235))*P(7)
       Fa3(590) = p3mup1*(p1mup3*F(1852)+2*(p2mup3*F(1858)+p3mup3*F(18
     -   61)))+p1mup1*(2*(p2mup3*F(1830)+p1mup3*F(1831))+p3mup3*F(228
     -   0))
       Fa3(591) = p2mup1*(p1mup3*F(2288)-2*(p3mup3*F(1861)-p2mup3*F(22
     -   94)))+p1mup1*p2mup3*F(2339)+2*p1mup1*(p1mup3*F(885)+p3mup3*F
     -   (894))*P(7)
       Fa3(592) = -2*p3mup1*p3mup3*F(902)+p1mup1*p2mup3*F(1000)+p1mup1
     -   *p3mup3*F(1004)+p1mup3*(p3mup1*F(896)+p1mup1*F(1003)+p2mup1*
     -   F(2680))+p2mup1*p3mup3*F(2682)+p2mup1*p2mup3*F(2683)+p2mup3*
     -   p3mup1*F(2687)
       Fa3(593) = mup1mup4*(p2mup3*F(1000)+p1mup3*F(1003)+p3mup3*F(100
     -   4))+mup1mup3*(p3mup4*F(1005)+p1mup4*F(1007)+p4mup4*F(1011)+p
     -   2mup4*F(1770))+mup3mup4*(p3mup1*F(2590)+p2mup1*F(2669))
       Fa3(594) = mup1mup3*(p1mup4*F(1697)+2*(p2mup4*F(1380)+p4mup4*F(
     -   1705))+p3mup4*F(2588))+mup3mup4*p1mup1*F(2668)+mup3mup4*p3mu
     -   p1*F(2671)+mup1mup4*p1mup3*F(2680)+mup1mup4*p3mup3*F(2682)+m
     -   up1mup4*p2mup3*F(2683)
       Fa3(595) = mup1mup3*(p1mup4*F(897)+2*(p3mup4*F(903)+p4mup4*F(90
     -   6))+p2mup4*F(2587))+mup3mup4*(p1mup1*F(2413)+p2mup1*F(2670))
     -   +mup1mup4*(p1mup3*F(896)-2*p3mup3*F(902)+p2mup3*F(2687))
       Fa3(596) = p2mup1*p2mup2*F(1679)+p2mup1*p3mup2*F(1832)+p2mup2*p
     -   3mup1*F(1928)+2*(p1mup2*(p2mup1*F(1676)+p3mup1*F(1931))+p3mu
     -   p1*p3mup2*F(1935))
       Fa3(597) = p3mup1*(p1mup2*F(1851)+2*(p2mup2*F(1857)+p3mup2*F(18
     -   60)))+p1mup1*(2*(p1mup2*F(1320)+p2mup2*F(1323))+p3mup2*F(227
     -   9))
       Fa3(598) = p1mup2*p2mup1*F(2287)+p1mup1*p2mup2*F(2338)+p1mup1*p
     -   1mup2*F(2341)+2*(p2mup1*p2mup2*F(2293)+p2mup1*p3mup2*F(2296)
     -   +p1mup1*p3mup2*F(2345))
       Fa3(599) = 2*p2mup1*p2mup2*F(1379)+p1mup1*p2mup2*F(1510)+p1mup1
     -   *p3mup2*F(1514)+p2mup1*p3mup2*F(2684)+p1mup2*(p2mup1*F(1373)
     -   +2*p1mup1*F(1513)+p3mup1*F(2685))+p2mup2*p3mup1*F(2688)+p3mu
     -   p1*p3mup2*F(2689)
       Fa3(600) = mup1mup4*p2mup2*F(1510)+mup1mup2*(p3mup4*F(1272)+2*p
     -   1mup4*F(1274)+p4mup4*F(1277)+p2mup4*F(1511))+2*mup1mup4*p1mu
     -   p2*F(1513)+mup1mup4*p3mup2*F(1514)+mup2mup4*p2mup1*F(2412)+m
     -   up2mup4*p3mup1*F(2678)
       Fa3(601) = mup1mup2*(p1mup4*F(1374)-2*(p2mup4*F(1380)-p4mup4*F(
     -   1383))+p3mup4*F(2411))+mup2mup4*(p1mup1*F(2589)+p3mup1*F(267
     -   9))+mup1mup4*(p1mup2*F(1373)+2*p2mup2*F(1379)+p3mup2*F(2684)
     -   )
       Fa3(602) = mup1mup2*(p1mup4*F(1237)-2*(p4mup4*F(906)-p3mup4*F(1
     -   242))+p2mup4*F(2410))+mup2mup4*p1mup1*F(2676)+mup2mup4*p2mup
     -   1*F(2677)+mup1mup4*p1mup2*F(2685)+mup1mup4*p2mup2*F(2688)+mu
     -   p1mup4*p3mup2*F(2689)
       Fa3(603) = mup1mup2*(p2mup3*F(1267)+2*p1mup3*F(1270)+p3mup3*F(1
     -   271))+mup1mup3*(p2mup2*F(1769)+2*p1mup2*F(1772)+p3mup2*F(177
     -   3))
       Fa3(604) = mup1mup2*(p1mup3*F(2664)+p2mup3*F(2665)+p3mup3*F(266
     -   7))+mup1mup3*(p1mup2*F(1696)+2*p2mup2*F(1702)+p3mup2*F(2674)
     -   )+mup2mup3*p1mup1*F(2696)
       Fa3(605) = mup1mup2*(2*p3mup3*F(902)+p1mup3*F(1236)+p2mup3*F(26
     -   66))+mup1mup3*(p1mup2*F(2672)+p2mup2*F(2673)+p3mup2*F(2675))
     -   +mup2mup3*p1mup1*F(2697)
       Fa3(606) = mup1mup4*(p4mup5*F(1232)+p2mup5*F(2340)+2*(p1mup5*F(
     -   276)+p3mup5*F(888))*P(7))
       Fa3(607) = mup1mup3*(p3mup4*F(1669)+2*(p2mup4*F(1671)+p1mup4*F(
     -   1672))+p4mup4*F(1674))+mup1mup4*(2*(p2mup3*F(1830)+p1mup3*F(
     -   1831))+p3mup3*F(2280))+mup3mup4*p3mup1*F(2692)
       Fa3(608) = mup1mup2*(p2mup4*F(1227)+2*(p1mup4*F(880)+p3mup4*F(1
     -   229))+p4mup4*F(1231))+mup1mup4*(p2mup2*F(2338)+p1mup2*F(2341
     -   )+2*p3mup2*F(2345))+mup2mup4*p2mup1*F(2693)
       Fa3(609) = 4*mup1mup2*(p1mup5*(2*(p2mup4*F(268)+p3mup4*F(270))+
     -   p4mup4*F(275))+p3mup5*p4mup4*F(290)+p1mup4*(p4mup5*F(273)+2*
     -   (p1mup5*F(130)+p3mup5*F(278)+p2mup5*F(637)))+p2mup5*p4mup4*F
     -   (671)+p3mup4*p4mup5*F(887)+p4mup4*p4mup5*F(891)+2*(p2mup4*p3
     -   mup5*F(277)+p3mup4*p3mup5*F(282)+p2mup5*p3mup4*F(636)+p2mup4
     -   *p2mup5*F(1279))+p2mup4*p4mup5*F(1328))+mup1mup2*mup4mup5*F(
     -   2586)
       Fa3(610) = mup1mup2*(mup4mup5*F(2413)+4*(p1mup4*p2mup5*F(327)-p
     -   2mup5*p4mup4*F(1016)-p2mup4*(p1mup5*F(276)+p3mup5*F(888)+p4m
     -   up5*F(1012)+p2mup5*F(1519))+p2mup5*p3mup4*F(2415)))
       Fa3(611) = mup1mup3*(2*p1mup2*F(1676)+p2mup2*F(1679)+p3mup2*F(1
     -   832))+mup1mup2*(-2*(p2mup3*F(1830)+p1mup3*F(1831))+p3mup3*F(
     -   1833))
       Fa3(612) = mup1mup2*p2mup3*F(1929)+mup1mup3*(p2mup2*F(1928)+2*(
     -   p1mup2*F(1931)+p3mup2*F(1935)))+2*mup1mup2*(p1mup3*F(271)+p3
     -   mup3*F(1235))*P(7)
       Fa3(613) = mup1mup3*(p2mup4*F(878)-2*(p1mup4*F(880)+p3mup4*F(88
     -   1))+p4mup4*F(883))+mup1mup4*p2mup3*F(2339)+mup3mup4*p2mup1*F
     -   (2690)+2*mup1mup4*(p1mup3*F(885)+p3mup3*F(894))*P(7)
       Fa3(614) = mup1mup2*(p3mup4*F(1313)+2*(p2mup4*F(1315)+p1mup4*F(
     -   1316))+p4mup4*F(1318))+mup1mup4*(2*(p1mup2*F(1320)+p2mup2*F(
     -   1323))+p3mup2*F(2279))+mup2mup4*p3mup1*F(2691)
       Fa3(615) = 4*mup2mup3*(p1mup5*(2*(p2mup4*F(268)+p3mup4*F(270))+
     -   p4mup4*F(275))+p3mup5*p4mup4*F(290)+p1mup4*(p4mup5*F(273)+2*
     -   (p1mup5*F(130)+p3mup5*F(278)+p2mup5*F(637)))+p2mup5*p4mup4*F
     -   (671)+p3mup4*p4mup5*F(887)+p4mup4*p4mup5*F(891)+2*(p2mup4*p3
     -   mup5*F(277)+p3mup4*p3mup5*F(282)+p2mup5*p3mup4*F(636)+p2mup4
     -   *p2mup5*F(1279))+p2mup4*p4mup5*F(1328))+mup2mup3*mup4mup5*F(
     -   2586)
       Fa3(616) = mup2mup3*(4*(p1mup4*(p3mup5*F(290)+p4mup5*F(294)+p2m
     -   up5*F(671)+p1mup5*F(673))+p1mup5*(p3mup4*F(291)+p4mup4*F(305
     -   )+p2mup4*F(1102)))+8*(p2mup4*p3mup5*F(300)+p3mup4*p3mup5*F(3
     -   02)+p2mup4*p4mup5*F(306)+p2mup5*p3mup4*F(678)+p3mup5*p4mup4*
     -   F(909)+p4mup4*p4mup5*F(911)+p3mup4*p4mup5*F(912)+p2mup4*p2mu
     -   p5*F(1049)+p2mup5*p4mup4*F(1395))+mup4mup5*F(2698))
       Fa3(617) = mup1mup4*(p1mup3*F(2288)-2*(p3mup3*F(1861)-p2mup3*F(
     -   2294)))+mup1mup3*(p1mup4*F(2402)+p2mup4*F(2403)+p3mup4*F(240
     -   4)+p4mup4*F(2405))+mup3mup4*p1mup1*F(2692)
       Fa3(618) = mup1mup3*(p1mup4*p4mup5*F(294)+p1mup5*p4mup4*F(305)-
     -   p1mup5*p3mup4*F(613)+2*p2mup4*p2mup5*F(1049))
       Fa3(619) = mup1mup4*(p1mup2*F(2287)+2*(p2mup2*F(2293)+p3mup2*F(
     -   2296)))+mup1mup2*(p1mup4*F(2406)+p2mup4*F(2407)+p3mup4*F(240
     -   8)+p4mup4*F(2409))+mup2mup4*p1mup1*F(2691)
       Fa3(620) = mup1mup2*(4*(p1mup4*(p3mup5*F(290)+p4mup5*F(294)+p2m
     -   up5*F(671)+p1mup5*F(673))+p1mup5*(p3mup4*F(291)+p4mup4*F(305
     -   )+p2mup4*F(1102)))+8*(p2mup4*p3mup5*F(300)+p3mup4*p3mup5*F(3
     -   02)+p2mup4*p4mup5*F(306)+p2mup5*p3mup4*F(678)+p3mup5*p4mup4*
     -   F(909)+p4mup4*p4mup5*F(911)+p3mup4*p4mup5*F(912)+p2mup4*p2mu
     -   p5*F(1049)+p2mup5*p4mup4*F(1395))+mup4mup5*F(2400))
       Fa3(621) = mup1mup2*(mup4mup5*F(2412)-4*(-(p1mup5*p3mup4*F(276)
     -   )+p1mup4*p3mup5*F(327)-p3mup4*p3mup5*F(335)-p2mup5*p3mup4*F(
     -   359)-p3mup4*p4mup5*F(1012)-p3mup5*p4mup4*F(1016)+p2mup4*p3mu
     -   p5*F(2414)))
       Fa3(622) = 4*mup2mup3*(p1mup5*p3mup4*F(276)-p1mup4*p3mup5*F(327
     -   )+p3mup4*p3mup5*F(335)+p2mup5*p3mup4*F(359)+p3mup4*p4mup5*F(
     -   1012)+p3mup5*p4mup4*F(1016)-p2mup4*p3mup5*F(2414))+mup2mup3*
     -   mup4mup5*F(2696)
       Fa3(623) = mup1mup3*(p1mup2*F(1851)+2*(p2mup2*F(1857)+p3mup2*F(
     -   1860)))+mup1mup2*(p1mup3*F(1852)+2*(p2mup3*F(1858)+p3mup3*F(
     -   1861)))+mup2mup3*p1mup1*F(2699)
       Fa3(624) = -4*mup2mup3*(-(p1mup4*p2mup5*F(327))+p2mup5*p4mup4*F
     -   (1016)+p2mup4*(p1mup5*F(276)+p3mup5*F(888)+p4mup5*F(1012)+p2
     -   mup5*F(1519))-p2mup5*p3mup4*F(2415))+mup2mup3*mup4mup5*F(269
     -   7)
       Fa3(625) = -4*p4mup1*(p3mup3*p3mup5*F(133)+p3mup3*p4mup5*F(293)
     -   +p3mup5*p4mup3*F(298)+p1mup5*(p3mup3*F(132)+2*(p2mup3*F(268)
     -   +p4mup3*F(338)))+p1mup3*(p3mup5*F(131)+2*(p1mup5*F(130)+p4mu
     -   p5*F(333)+p2mup5*F(637)))+2*(p2mup3*p4mup5*F(332)+p4mup3*p4m
     -   up5*F(1015)+p2mup3*p2mup5*F(1279)+p2mup5*p4mup3*F(1283))+p2m
     -   up3*p3mup5*F(1841)+p2mup5*p3mup3*F(1842))+mup3mup5*p4mup1*F(
     -   2700)
       Fa3(626) = p4mup1*(4*(p1mup5*p4mup3*F(298)+2*(p2mup3*p4mup5*F(3
     -   00)+p3mup5*p4mup3*F(302)+p2mup3*p3mup5*F(304)+p2mup5*p4mup3*
     -   F(678)+p4mup3*p4mup5*F(912)+p2mup3*p2mup5*F(1060))+p1mup3*(p
     -   1mup5*F(164)+p3mup5*F(165)+p4mup5*F(293)+p2mup5*F(1842))+p3m
     -   up3*(p1mup5*F(167)+2*(p3mup5*F(166)+p4mup5*F(302)+p2mup5*F(1
     -   874)))+p1mup5*p2mup3*F(2298))+mup3mup5*F(2701))
       Fa3(627) = p3mup1*p4mup3*F(882)+p4mup1*(p2mup3*F(1000)+p1mup3*F
     -   (1003)+p3mup3*F(1004)-2*p4mup3*F(1009))+p2mup1*p4mup3*F(2283
     -   )
       Fa3(628) = p1mup1*p4mup3*F(1836)+p3mup1*p4mup3*F(2702)+p4mup1*(
     -   p1mup3*F(2680)+p3mup3*F(2682)+p2mup3*F(2683)+2*p4mup3*F(2703
     -   ))
       Fa3(629) = p1mup3*p4mup1*F(896)-2*(p3mup3*p4mup1*F(902)+p4mup1*
     -   p4mup3*F(905))+p1mup1*p4mup3*F(1230)+p2mup3*p4mup1*F(2687)+p
     -   2mup1*p4mup3*F(2704)
       Fa3(630) = p4mup3*(p3mup1*F(905)+p1mup1*F(1009)-p2mup1*F(2703))
       Fa3(631) = 4*p2mup1*(p3mup3*p3mup5*F(133)+p3mup3*p4mup5*F(293)+
     -   p3mup5*p4mup3*F(298)+p1mup5*(p3mup3*F(132)+2*(p2mup3*F(268)+
     -   p4mup3*F(338)))+p1mup3*(p3mup5*F(131)+2*(p1mup5*F(130)+p4mup
     -   5*F(333)+p2mup5*F(637)))+2*(p2mup3*p4mup5*F(332)+p4mup3*p4mu
     -   p5*F(1015)+p2mup3*p2mup5*F(1279)+p2mup5*p4mup3*F(1283))+p2mu
     -   p3*p3mup5*F(1841)+p2mup5*p3mup3*F(1842))+mup3mup5*p2mup1*F(2
     -   705)
       Fa3(632) = p2mup1*(4*(p1mup5*p4mup3*F(298)+2*(p2mup3*p4mup5*F(3
     -   00)+p3mup5*p4mup3*F(302)+p2mup3*p3mup5*F(304)+p2mup5*p4mup3*
     -   F(678)+p4mup3*p4mup5*F(912)+p2mup3*p2mup5*F(1060))+p1mup3*(p
     -   1mup5*F(164)+p3mup5*F(165)+p4mup5*F(293)+p2mup5*F(1842))+p3m
     -   up3*(p1mup5*F(167)+2*(p3mup5*F(166)+p4mup5*F(302)+p2mup5*F(1
     -   874)))+p1mup5*p2mup3*F(2298))+mup3mup5*F(2701))
       Fa3(633) = p2mup2*p4mup1*F(1510)+2*p1mup2*p4mup1*F(1513)+p3mup2
     -   *p4mup1*F(1514)+p3mup1*p4mup2*F(2342)
       Fa3(634) = p1mup2*p4mup1*F(1373)+2*p2mup2*p4mup1*F(1379)+p3mup2
     -   *p4mup1*F(2684)+p3mup1*p4mup2*F(2706)
       Fa3(635) = p4mup1*(p1mup2*F(2685)+p2mup2*F(2688)+p3mup2*F(2689)
     -   )
       Fa3(636) = -2*mup1mup2*p4mup3*F(1009)+mup2mup3*p2mup1*F(2696)+m
     -   up2mup3*p3mup1*F(2697)
       Fa3(637) = 2*mup1mup2*p4mup3*F(2703)+mup2mup3*p3mup1*F(2708)
       Fa3(638) = -2*mup1mup2*p4mup3*F(905)+mup2mup3*p2mup1*F(2709)
       Fa3(639) = -(mup1mup4*(-(mup3mup5*F(2695))+4*(p4mup3*(p1mup5*F(
     -   276)+p2mup5*F(359)+p3mup5*F(888))+p4mup5*(p1mup3*F(271)+p4mu
     -   p3*F(892)-p3mup3*F(894)+p2mup3*F(2710)))))
       Fa3(640) = 4*mup1mup4*(p3mup3*p3mup5*F(133)+p3mup3*p4mup5*F(293
     -   )+p3mup5*p4mup3*F(298)+p1mup5*(p3mup3*F(132)+2*(p2mup3*F(268
     -   )+p4mup3*F(338)))+p1mup3*(p3mup5*F(131)+2*(p1mup5*F(130)+p4m
     -   up5*F(333)+p2mup5*F(637)))+2*(p2mup3*p4mup5*F(332)+p4mup3*p4
     -   mup5*F(1015)+p2mup3*p2mup5*F(1279)+p2mup5*p4mup3*F(1283))+p2
     -   mup3*p3mup5*F(1841)+p2mup5*p3mup3*F(1842))+mup1mup4*mup3mup5
     -   *F(2705)
       Fa3(641) = mup1mup2*p4mup3*F(2283)+mup2mup3*p3mup1*F(2695)
       Fa3(642) = mup1mup2*p4mup3*F(882)+mup2mup3*p2mup1*F(2699)
       Fa3(643) = -4*mup1mup2*(p3mup3*p3mup5*F(133)+p3mup3*p4mup5*F(29
     -   3)+p3mup5*p4mup3*F(298)+p1mup5*(p3mup3*F(132)+2*(p2mup3*F(26
     -   8)+p4mup3*F(338)))+p1mup3*(p3mup5*F(131)+2*(p1mup5*F(130)+p4
     -   mup5*F(333)+p2mup5*F(637)))+2*(p2mup3*p4mup5*F(332)+p4mup3*p
     -   4mup5*F(1015)+p2mup3*p2mup5*F(1279)+p2mup5*p4mup3*F(1283))+p
     -   2mup3*p3mup5*F(1841)+p2mup5*p3mup3*F(1842))+mup1mup2*mup3mup
     -   5*F(2700)
       Fa3(644) = mup1mup2*(mup3mup5*F(2413)+4*(p1mup3*p2mup5*F(271)-p
     -   2mup5*p3mup3*F(894)-p2mup3*(p1mup5*F(276)+p3mup5*F(888)+p4mu
     -   p5*F(1012)+p2mup5*F(1554))-p2mup5*p4mup3*F(2711)))
       Fa3(645) = mup2mup4*(mup3mup5*F(2693)+4*(p4mup3*(p1mup5*F(276)+
     -   p2mup5*F(359)+p3mup5*F(888))+p4mup5*(p1mup3*F(271)+p4mup3*F(
     -   892)-p3mup3*F(894)+p2mup3*F(2710))))
       Fa3(646) = mup2mup4*(4*(p1mup5*p4mup3*F(298)+2*(p2mup3*p4mup5*F
     -   (300)+p3mup5*p4mup3*F(302)+p2mup3*p3mup5*F(304)+p2mup5*p4mup
     -   3*F(678)+p4mup3*p4mup5*F(912)+p2mup3*p2mup5*F(1060))+p1mup3*
     -   (p1mup5*F(164)+p3mup5*F(165)+p4mup5*F(293)+p2mup5*F(1842))+p
     -   3mup3*(p1mup5*F(167)+2*(p3mup5*F(166)+p4mup5*F(302)+p2mup5*F
     -   (1874)))+p1mup5*p2mup3*F(2298))+mup3mup5*F(2712))
       Fa3(647) = -4*mup1mup4*(p1mup5*p4mup3*F(298)+2*(p2mup3*p4mup5*F
     -   (300)+p3mup5*p4mup3*F(302)+p2mup3*p3mup5*F(304)+p2mup5*p4mup
     -   3*F(678)+p4mup3*p4mup5*F(912)+p2mup3*p2mup5*F(1060))+p1mup3*
     -   (p1mup5*F(164)+p3mup5*F(165)+p4mup5*F(293)+p2mup5*F(1842))-p
     -   1mup5*p2mup3*F(1863)+p3mup3*(p1mup5*F(167)+2*(p3mup5*F(166)+
     -   p4mup5*F(302)+p2mup5*F(1874))))+mup1mup4*mup3mup5*F(2713)
       Fa3(648) = mup3mup5*p4mup1*F(2413)+mup3mup5*p2mup1*F(2695)-4*(-
     -   (p1mup3*(p2mup5*p4mup1-p2mup1*p4mup5)*F(271))+p1mup5*(p2mup3
     -   *p4mup1+p2mup1*p4mup3)*F(276)+p2mup1*p2mup5*p4mup3*F(359)+p2
     -   mup3*p3mup5*p4mup1*F(888)+p2mup1*p3mup5*p4mup3*F(888)+p2mup1
     -   *p4mup3*p4mup5*F(892)+p2mup5*p3mup3*p4mup1*F(894)-p2mup1*p3m
     -   up3*p4mup5*F(894)+p2mup3*p4mup1*p4mup5*F(1012)+p2mup3*p2mup5
     -   *p4mup1*F(1554)+p2mup1*p2mup3*p4mup5*F(2710)+p2mup5*p4mup1*p
     -   4mup3*F(2711))
       Fa3(649) = mup1mup4*(mup3mup5*F(2590)+4*(p2mup3*(p1mup5*F(276)+
     -   p3mup5*F(888)+p4mup5*F(1012))+p2mup5*(-(p1mup3*F(271))+p3mup
     -   3*F(894)+p2mup3*F(1554)+p4mup3*F(2711))))
       Fa3(650) = mup1mup2*(mup3mup5*F(2699)+4*(p4mup3*(p1mup5*F(276)+
     -   p2mup5*F(359)+p3mup5*F(888))+p4mup5*(p1mup3*F(271)+p4mup3*F(
     -   892)-p3mup3*F(894)+p2mup3*F(2710))))
       Fa3(651) = -4*mup2mup4*(p3mup3*p3mup5*F(133)+p3mup3*p4mup5*F(29
     -   3)+p3mup5*p4mup3*F(298)+p1mup5*(p3mup3*F(132)+2*(p2mup3*F(26
     -   8)+p4mup3*F(338)))+p1mup3*(p3mup5*F(131)+2*(p1mup5*F(130)+p4
     -   mup5*F(333)+p2mup5*F(637)))+2*(p2mup3*p4mup5*F(332)+p4mup3*p
     -   4mup5*F(1015)+p2mup3*p2mup5*F(1279)+p2mup5*p4mup3*F(1283))+p
     -   2mup3*p3mup5*F(1841)+p2mup5*p3mup3*F(1842))+mup2mup4*mup3mup
     -   5*F(2700)
       Fa3(652) = mup1mup2*(4*(p1mup5*p4mup3*F(298)+2*(p2mup3*p4mup5*F
     -   (300)+p3mup5*p4mup3*F(302)+p2mup3*p3mup5*F(304)+p2mup5*p4mup
     -   3*F(678)+p4mup3*p4mup5*F(912)+p2mup3*p2mup5*F(1060))+p1mup3*
     -   (p1mup5*F(164)+p3mup5*F(165)+p4mup5*F(293)+p2mup5*F(1842))+p
     -   3mup3*(p1mup5*F(167)+2*(p3mup5*F(166)+p4mup5*F(302)+p2mup5*F
     -   (1874)))+p1mup5*p2mup3*F(2298))+mup3mup5*F(2701))
       Fa3(653) = mup2mup4*(mup3mup5*F(2678)+4*(p1mup3*p2mup5*F(271)-p
     -   2mup5*p3mup3*F(894)-p2mup3*(p1mup5*F(276)+p3mup5*F(888)+p4mu
     -   p5*F(1012)+p2mup5*F(1554))-p2mup5*p4mup3*F(2711)))
       Fa3(654) = p2mup1*p4mup3*F(353)+p2mup3*p4mup1*F(7441)
       Fa3(655) = 2*p3mup1*p4mup3*F(335)+p1mup1*p4mup3*F(353)+p4mup1*(
     -   p1mup3*F(352)+2*(p3mup3*F(7442)+p4mup3*F(7444)+p2mup3*F(7448
     -   )))
       Fa3(656) = 2*p3mup1*p4mup2*F(335)+p4mup1*(p2mup2*F(1626)+2*(p1m
     -   up2*F(1012)+p3mup2*F(7442)))
       Fa3(657) = 2*p3mup1*p4mup2*F(335)+p4mup1*(p1mup2*F(352)+2*(p3mu
     -   p2*F(7442)+p2mup2*F(7448)))
       Fa3(658) = mup1mup2*p4mup3*F(353)+mup2mup3*p3mup1*F(1051)
       Fa3(659) = p2mup1*p4mup3*F(352)+p2mup3*p4mup1*F(7441)
       Fa3(660) = p1mup3*p4mup1*F(352)+p1mup1*p4mup3*F(352)+2*(p3mup3*
     -   p4mup1*F(7442)+p3mup1*p4mup3*F(7442)+p4mup1*p4mup3*F(7444)+p
     -   2mup3*p4mup1*F(7448))
       Fa3(661) = p2mup2*p4mup1*F(1626)+2*(p1mup2*p4mup1*F(1012)+(p3mu
     -   p2*p4mup1+p3mup1*p4mup2)*F(7442))
       Fa3(662) = p1mup2*p4mup1*F(352)+2*((p3mup2*p4mup1+p3mup1*p4mup2
     -   )*F(7442)+p2mup2*p4mup1*F(7448))
       Fa3(663) = mup1mup2*p4mup3*F(352)+mup2mup3*p3mup1*F(7440)
       Fa3(664) = p2mup1*p4mup3*F(430)+p2mup3*p4mup1*F(1095)
       Fa3(665) = 2*p3mup1*p4mup3*F(335)+p1mup1*p4mup3*F(430)+p4mup1*(
     -   p1mup3*F(438)+2*(p2mup3*F(434)+p4mup3*F(436)+p3mup3*F(439)))
       Fa3(666) = 2*p3mup1*p4mup2*F(335)+p4mup1*(p2mup2*F(1653)+2*(p3m
     -   up2*F(439)+p1mup2*F(1654)))
       Fa3(667) = 2*p3mup1*p4mup2*F(335)+p4mup1*(p1mup2*F(438)+2*(p2mu
     -   p2*F(434)+p3mup2*F(439)))
       Fa3(668) = mup1mup2*p4mup3*F(430)+mup2mup3*p3mup1*F(1051)
       Fa3(669) = 4*p4mup1*(p3mup3*p3mup5*F(133)+p3mup3*p4mup5*F(293)+
     -   p3mup5*p4mup3*F(298)+p1mup5*(p3mup3*F(132)+2*(p2mup3*F(268)+
     -   p4mup3*F(338)))+p1mup3*(p3mup5*F(131)+2*(p1mup5*F(130)+p4mup
     -   5*F(333)+p2mup5*F(637)))+2*(p2mup3*p4mup5*F(332)+p4mup3*p4mu
     -   p5*F(1015)+p2mup3*p2mup5*F(1279)+p2mup5*p4mup3*F(1283))+p2mu
     -   p3*p3mup5*F(1841)+p2mup5*p3mup3*F(1842))+mup3mup5*p4mup1*F(2
     -   705)
       Fa3(670) = -4*p4mup1*(p1mup5*p4mup3*F(298)+2*(p2mup3*p4mup5*F(3
     -   00)+p3mup5*p4mup3*F(302)+p2mup3*p3mup5*F(304)+p2mup5*p4mup3*
     -   F(678)+p4mup3*p4mup5*F(912)+p2mup3*p2mup5*F(1060))+p1mup3*(p
     -   1mup5*F(164)+p3mup5*F(165)+p4mup5*F(293)+p2mup5*F(1842))-p1m
     -   up5*p2mup3*F(1863)+p3mup3*(p1mup5*F(167)+2*(p3mup5*F(166)+p4
     -   mup5*F(302)+p2mup5*F(1874))))+mup3mup5*p4mup1*F(2713)
       Fa3(671) = p4mup1*(p2mup3*F(1267)+2*(p4mup3*F(1009)+p1mup3*F(12
     -   70))+p3mup3*F(1271))+p4mup3*(p3mup1*F(1230)+p2mup1*F(1836))
       Fa3(672) = p1mup1*p4mup3*F(2283)+p4mup1*(p1mup3*F(2664)+p2mup3*
     -   F(2665)+p3mup3*F(2667)-2*p4mup3*F(2703))+p3mup1*p4mup3*F(270
     -   4)
       Fa3(673) = p1mup1*p4mup3*F(882)+p4mup1*(2*(p3mup3*F(902)+p4mup3
     -   *F(905))+p1mup3*F(1236)+p2mup3*F(2666))+p2mup1*p4mup3*F(2702
     -   )
       Fa3(674) = -4*p2mup1*(p3mup3*p3mup5*F(133)+p3mup3*p4mup5*F(293)
     -   +p3mup5*p4mup3*F(298)+p1mup5*(p3mup3*F(132)+2*(p2mup3*F(268)
     -   +p4mup3*F(338)))+p1mup3*(p3mup5*F(131)+2*(p1mup5*F(130)+p4mu
     -   p5*F(333)+p2mup5*F(637)))+2*(p2mup3*p4mup5*F(332)+p4mup3*p4m
     -   up5*F(1015)+p2mup3*p2mup5*F(1279)+p2mup5*p4mup3*F(1283))+p2m
     -   up3*p3mup5*F(1841)+p2mup5*p3mup3*F(1842))+mup3mup5*p2mup1*F(
     -   2700)
       Fa3(675) = -4*p2mup1*(p1mup5*p4mup3*F(298)+2*(p2mup3*p4mup5*F(3
     -   00)+p3mup5*p4mup3*F(302)+p2mup3*p3mup5*F(304)+p2mup5*p4mup3*
     -   F(678)+p4mup3*p4mup5*F(912)+p2mup3*p2mup5*F(1060))+p1mup3*(p
     -   1mup5*F(164)+p3mup5*F(165)+p4mup5*F(293)+p2mup5*F(1842))-p1m
     -   up5*p2mup3*F(1863)+p3mup3*(p1mup5*F(167)+2*(p3mup5*F(166)+p4
     -   mup5*F(302)+p2mup5*F(1874))))+mup3mup5*p2mup1*F(2713)
       Fa3(676) = p2mup2*p4mup1*F(1769)+2*p1mup2*p4mup1*F(1772)+p3mup2
     -   *p4mup1*F(1773)+p3mup1*p4mup2*F(1932)
       Fa3(677) = p1mup2*p4mup1*F(1696)+2*p2mup2*p4mup1*F(1702)+p3mup2
     -   *p4mup1*F(2674)+p3mup1*p4mup2*F(2714)
       Fa3(678) = p4mup1*(p1mup2*F(2672)+p2mup2*F(2673)+p3mup2*F(2675)
     -   )
       Fa3(679) = 2*mup1mup2*p4mup3*F(1009)+mup2mup3*p2mup1*F(2681)+mu
     -   p2mup3*p3mup1*F(2686)
       Fa3(680) = -2*mup1mup2*p4mup3*F(2703)+mup2mup3*p3mup1*F(2709)
       Fa3(681) = 2*mup1mup2*p4mup3*F(905)+mup2mup3*p2mup1*F(2708)
       Fa3(682) = mup1mup4*(mup3mup5*F(2699)+4*(p4mup3*(p1mup5*F(276)+
     -   p2mup5*F(359)+p3mup5*F(888))+p4mup5*(p1mup3*F(271)+p4mup3*F(
     -   892)-p3mup3*F(894)+p2mup3*F(2710))))
       Fa3(683) = -4*mup1mup4*(p3mup3*p3mup5*F(133)+p3mup3*p4mup5*F(29
     -   3)+p3mup5*p4mup3*F(298)+p1mup5*(p3mup3*F(132)+2*(p2mup3*F(26
     -   8)+p4mup3*F(338)))+p1mup3*(p3mup5*F(131)+2*(p1mup5*F(130)+p4
     -   mup5*F(333)+p2mup5*F(637)))+2*(p2mup3*p4mup5*F(332)+p4mup3*p
     -   4mup5*F(1015)+p2mup3*p2mup5*F(1279)+p2mup5*p4mup3*F(1283))+p
     -   2mup3*p3mup5*F(1841)+p2mup5*p3mup3*F(1842))+mup1mup4*mup3mup
     -   5*F(2700)
       Fa3(684) = mup1mup2*p4mup3*F(1836)+mup2mup3*p3mup1*F(2699)
       Fa3(685) = mup1mup2*p4mup3*F(1230)+mup2mup3*p2mup1*F(2695)
       Fa3(686) = 4*mup1mup2*(p3mup3*p3mup5*F(133)+p3mup3*p4mup5*F(293
     -   )+p3mup5*p4mup3*F(298)+p1mup5*(p3mup3*F(132)+2*(p2mup3*F(268
     -   )+p4mup3*F(338)))+p1mup3*(p3mup5*F(131)+2*(p1mup5*F(130)+p4m
     -   up5*F(333)+p2mup5*F(637)))+2*(p2mup3*p4mup5*F(332)+p4mup3*p4
     -   mup5*F(1015)+p2mup3*p2mup5*F(1279)+p2mup5*p4mup3*F(1283))+p2
     -   mup3*p3mup5*F(1841)+p2mup5*p3mup3*F(1842))+mup1mup2*mup3mup5
     -   *F(2705)
       Fa3(687) = mup1mup2*(mup3mup5*F(2590)+4*(p2mup3*(p1mup5*F(276)+
     -   p3mup5*F(888)+p4mup5*F(1012))+p2mup5*(-(p1mup3*F(271))+p3mup
     -   3*F(894)+p2mup3*F(1554)+p4mup3*F(2711))))
       Fa3(688) = -(mup2mup4*(-(mup3mup5*F(2691))+4*(p4mup3*(p1mup5*F(
     -   276)+p2mup5*F(359)+p3mup5*F(888))+p4mup5*(p1mup3*F(271)+p4mu
     -   p3*F(892)-p3mup3*F(894)+p2mup3*F(2710)))))
       Fa3(689) = -4*mup2mup4*(p1mup5*p4mup3*F(298)+2*(p2mup3*p4mup5*F
     -   (300)+p3mup5*p4mup3*F(302)+p2mup3*p3mup5*F(304)+p2mup5*p4mup
     -   3*F(678)+p4mup3*p4mup5*F(912)+p2mup3*p2mup5*F(1060))+p1mup3*
     -   (p1mup5*F(164)+p3mup5*F(165)+p4mup5*F(293)+p2mup5*F(1842))-p
     -   1mup5*p2mup3*F(1863)+p3mup3*(p1mup5*F(167)+2*(p3mup5*F(166)+
     -   p4mup5*F(302)+p2mup5*F(1874))))+mup2mup4*mup3mup5*F(2716)
       Fa3(690) = mup1mup4*(4*(p1mup5*p4mup3*F(298)+2*(p2mup3*p4mup5*F
     -   (300)+p3mup5*p4mup3*F(302)+p2mup3*p3mup5*F(304)+p2mup5*p4mup
     -   3*F(678)+p4mup3*p4mup5*F(912)+p2mup3*p2mup5*F(1060))+p1mup3*
     -   (p1mup5*F(164)+p3mup5*F(165)+p4mup5*F(293)+p2mup5*F(1842))+p
     -   3mup3*(p1mup5*F(167)+2*(p3mup5*F(166)+p4mup5*F(302)+p2mup5*F
     -   (1874)))+p1mup5*p2mup3*F(2298))+mup3mup5*F(2701))
       Fa3(691) = mup3mup5*p4mup1*F(2590)+mup3mup5*p2mup1*F(2699)+4*(p
     -   1mup3*(-(p2mup5*p4mup1)+p2mup1*p4mup5)*F(271)+p1mup5*(p2mup3
     -   *p4mup1+p2mup1*p4mup3)*F(276)+p2mup1*p2mup5*p4mup3*F(359)+p2
     -   mup3*p3mup5*p4mup1*F(888)+p2mup1*p3mup5*p4mup3*F(888)+p2mup1
     -   *p4mup3*p4mup5*F(892)+p2mup5*p3mup3*p4mup1*F(894)-p2mup1*p3m
     -   up3*p4mup5*F(894)+p2mup3*p4mup1*p4mup5*F(1012)+p2mup3*p2mup5
     -   *p4mup1*F(1554)+p2mup1*p2mup3*p4mup5*F(2710)+p2mup5*p4mup1*p
     -   4mup3*F(2711))
       Fa3(692) = mup1mup4*(mup3mup5*F(2413)+4*(p1mup3*p2mup5*F(271)-p
     -   2mup5*p3mup3*F(894)-p2mup3*(p1mup5*F(276)+p3mup5*F(888)+p4mu
     -   p5*F(1012)+p2mup5*F(1554))-p2mup5*p4mup3*F(2711)))
       Fa3(693) = -(mup1mup2*(-(mup3mup5*F(2695))+4*(p4mup3*(p1mup5*F(
     -   276)+p2mup5*F(359)+p3mup5*F(888))+p4mup5*(p1mup3*F(271)+p4mu
     -   p3*F(892)-p3mup3*F(894)+p2mup3*F(2710)))))
       Fa3(694) = 4*mup2mup4*(p3mup3*p3mup5*F(133)+p3mup3*p4mup5*F(293
     -   )+p3mup5*p4mup3*F(298)+p1mup5*(p3mup3*F(132)+2*(p2mup3*F(268
     -   )+p4mup3*F(338)))+p1mup3*(p3mup5*F(131)+2*(p1mup5*F(130)+p4m
     -   up5*F(333)+p2mup5*F(637)))+2*(p2mup3*p4mup5*F(332)+p4mup3*p4
     -   mup5*F(1015)+p2mup3*p2mup5*F(1279)+p2mup5*p4mup3*F(1283))+p2
     -   mup3*p3mup5*F(1841)+p2mup5*p3mup3*F(1842))+mup2mup4*mup3mup5
     -   *F(2705)
       Fa3(695) = -4*mup1mup2*(p1mup5*p4mup3*F(298)+2*(p2mup3*p4mup5*F
     -   (300)+p3mup5*p4mup3*F(302)+p2mup3*p3mup5*F(304)+p2mup5*p4mup
     -   3*F(678)+p4mup3*p4mup5*F(912)+p2mup3*p2mup5*F(1060))+p1mup3*
     -   (p1mup5*F(164)+p3mup5*F(165)+p4mup5*F(293)+p2mup5*F(1842))-p
     -   1mup5*p2mup3*F(1863)+p3mup3*(p1mup5*F(167)+2*(p3mup5*F(166)+
     -   p4mup5*F(302)+p2mup5*F(1874))))+mup1mup2*mup3mup5*F(2713)
       Fa3(696) = mup2mup4*(mup3mup5*F(2676)+4*(p2mup3*(p1mup5*F(276)+
     -   p3mup5*F(888)+p4mup5*F(1012))+p2mup5*(-(p1mup3*F(271))+p3mup
     -   3*F(894)+p2mup3*F(1554)+p4mup3*F(2711))))
       Fa3(697) = mup2mup5*p4mup1*F(2589)+mup2mup5*p3mup1*F(2695)+4*(p
     -   1mup5*(p3mup2*p4mup1-p3mup1*p4mup2)*F(276)+p2mup5*(p3mup2*p4
     -   mup1-p3mup1*p4mup2)*F(359)-p3mup1*p3mup5*p4mup2*F(888)+p3mup
     -   2*p4mup1*p4mup5*F(1012)+p1mup2*p3mup5*p4mup1*F(1325)-p1mup2*
     -   p3mup1*p4mup5*F(1325)-p3mup1*p4mup2*p4mup5*F(1334)-p2mup2*p3
     -   mup5*p4mup1*F(1335)+p2mup2*p3mup1*p4mup5*F(1335)+p3mup2*p3mu
     -   p5*p4mup1*F(1843)-p3mup1*p3mup2*p4mup5*F(2717)+p3mup5*p4mup1
     -   *p4mup2*F(2718))
       Fa3(698) = -4*p4mup1*(p1mup2*p1mup5*F(668)+p1mup2*p2mup5*F(670)
     -   +p1mup5*p2mup2*F(676)+2*(p4mup2*p4mup5*F(306)+p2mup2*p2mup5*
     -   F(1036)+p2mup5*p4mup2*F(1049)+p2mup2*p4mup5*F(1049))+p1mup2*
     -   p4mup5*F(1387)+p1mup5*p4mup2*F(1394)+p3mup5*(2*(p4mup2*F(300
     -   )+p3mup2*F(304)+p2mup2*F(1060))+p1mup2*F(1864))+p3mup2*(2*(p
     -   4mup5*F(300)+p2mup5*F(1060))+p1mup5*F(1873)))+mup2mup5*p4mup
     -   1*F(2719)
       Fa3(699) = 4*p3mup1*(p1mup2*p1mup5*F(668)+p1mup2*p2mup5*F(670)+
     -   p1mup5*p2mup2*F(676)+2*(p4mup2*p4mup5*F(306)+p2mup2*p2mup5*F
     -   (1036)+p2mup5*p4mup2*F(1049)+p2mup2*p4mup5*F(1049))+p1mup2*p
     -   4mup5*F(1387)+p1mup5*p4mup2*F(1394)+p3mup5*(2*(p4mup2*F(300)
     -   +p3mup2*F(304)+p2mup2*F(1060))+p1mup2*F(1864))+p3mup2*(2*(p4
     -   mup5*F(300)+p2mup5*F(1060))+p1mup5*F(1873)))+mup2mup5*p3mup1
     -   *F(2720)
       Fa3(700) = p4mup2*(p2mup1*F(1317)+2*p4mup1*F(1517))
       Fa3(701) = p4mup2*(2*p4mup1*F(1382)+p1mup1*F(1673))
       Fa3(702) = p4mup2*(p1mup1*F(1932)+p2mup1*F(2714)+p4mup1*F(2715)
     -   )
       Fa3(703) = p4mup2*(p2mup1*F(1704)+p1mup1*F(1776))
       Fa3(704) = 4*(p3mup1*(p1mup2*(2*(p1mup5*F(130)+p3mup5*F(278)+p4
     -   mup5*F(333))+p2mup5*F(962))+p2mup2*p2mup5*F(1022)+p1mup5*(p2
     -   mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028)))+2*(p4mup2*p4mu
     -   p5*F(1015)+p3mup5*p4mup2*F(1047)+p3mup2*p4mup5*F(1047)+p3mup
     -   2*p3mup5*F(1059))+p2mup2*p4mup5*F(1387)+p2mup5*p4mup2*F(1394
     -   )+p2mup2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))-p2mup1*p3mup5
     -   *p4mup2*F(2718))+mup2mup5*p3mup1*F(2721)
       Fa3(705) = p2mup1*p4mup2*(2*(p3mup5*F(300)+p4mup5*F(306)+p2mup5
     -   *F(1049))+p1mup5*F(1394))
       Fa3(706) = -4*(-(p2mup1*p2mup5*p4mup2*F(359))-p2mup1*p3mup5*p4m
     -   up2*F(888)+p1mup2*p4mup1*(2*(p1mup5*F(130)+p3mup5*F(278)+p4m
     -   up5*F(333))+p2mup5*F(962))+p2mup2*p2mup5*p4mup1*F(1022)+p1mu
     -   p5*(-(p2mup1*p4mup2*F(276))+p2mup2*p4mup1*F(964)+2*(p4mup1*p
     -   4mup2*F(338)+p3mup2*p4mup1*F(1028)))+2*(p4mup1*p4mup2*p4mup5
     -   *F(1015)+p3mup5*p4mup1*p4mup2*F(1047)+p3mup2*p4mup1*p4mup5*F
     -   (1047)+p3mup2*p3mup5*p4mup1*F(1059))-p2mup1*p4mup2*p4mup5*F(
     -   1334)+p2mup2*p4mup1*p4mup5*F(1387)+p2mup5*p4mup1*p4mup2*F(13
     -   94)+p2mup2*p3mup5*p4mup1*F(1864)+p2mup5*p3mup2*p4mup1*F(1873
     -   ))+mup2mup5*p4mup1*F(2722)
       Fa3(707) = mup1mup4*p4mup2*F(1932)+mup2mup4*p4mup1*F(2678)
       Fa3(708) = mup1mup4*(mup2mup5*F(2699)+4*(p1mup5*p4mup2*F(276)+p
     -   2mup5*p4mup2*F(359)+p3mup5*p4mup2*F(888)+p1mup2*p4mup5*F(132
     -   5)+p4mup2*p4mup5*F(1334)-p2mup2*p4mup5*F(1335)+p3mup2*p4mup5
     -   *F(2717)))
       Fa3(709) = 4*mup1mup4*(p1mup2*(2*(p1mup5*F(130)+p3mup5*F(278)+p
     -   4mup5*F(333))+p2mup5*F(962))+p2mup2*p2mup5*F(1022)+p1mup5*(p
     -   2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028)))+2*(p4mup2*p4m
     -   up5*F(1015)+p3mup5*p4mup2*F(1047)+p3mup2*p4mup5*F(1047)+p3mu
     -   p2*p3mup5*F(1059))+p2mup2*p4mup5*F(1387)+p2mup5*p4mup2*F(139
     -   4)+p2mup2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))+mup1mup4*mup
     -   2mup5*F(2721)
       Fa3(710) = mup1mup3*p4mup2*F(1317)+mup2mup3*p4mup1*F(2681)
       Fa3(711) = mup1mup3*p4mup2*F(2342)+mup2mup3*p4mup1*F(2686)
       Fa3(712) = mup1mup3*(mup2mup5*F(2589)+4*(p1mup2*p3mup5*F(1325)-
     -   p2mup2*p3mup5*F(1335)+p3mup2*(p1mup5*F(276)+p2mup5*F(359)+p4
     -   mup5*F(1012)+p3mup5*F(1843))+p3mup5*p4mup2*F(2718)))
       Fa3(713) = -4*mup1mup3*(p1mup2*(2*(p1mup5*F(130)+p3mup5*F(278)+
     -   p4mup5*F(333))+p2mup5*F(962))+p2mup2*p2mup5*F(1022)+p1mup5*(
     -   p2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028)))+2*(p4mup2*p4
     -   mup5*F(1015)+p3mup5*p4mup2*F(1047)+p3mup2*p4mup5*F(1047)+p3m
     -   up2*p3mup5*F(1059))+p2mup2*p4mup5*F(1387)+p2mup5*p4mup2*F(13
     -   94)+p2mup2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))+mup1mup3*mu
     -   p2mup5*F(2722)
       Fa3(714) = mup1mup4*p4mup2*F(1673)+mup2mup4*p4mup1*F(2412)
       Fa3(715) = mup1mup2*mup3mup4*F(2690)+mup1mup3*mup2mup4*F(2693)+
     -   mup1mup4*mup2mup3*F(2695)
       Fa3(716) = mup1mup3*mup2mup4*F(2412)+mup1mup2*mup3mup4*F(2669)+
     -   mup1mup4*mup2mup3*F(2681)
       Fa3(717) = mup1mup2*mup3mup4*F(2590)+mup1mup3*mup2mup4*F(2678)+
     -   mup1mup4*mup2mup3*F(2686)
       Fa3(718) = -(mup3mup4*(-(mup2mup5*F(2690))+4*(p1mup5*p4mup2*F(2
     -   76)+p2mup5*p4mup2*F(359)+p3mup5*p4mup2*F(888)+p1mup2*p4mup5*
     -   F(1325)+p4mup2*p4mup5*F(1334)-p2mup2*p4mup5*F(1335)+p3mup2*p
     -   4mup5*F(2717))))
       Fa3(719) = -4*mup3mup4*(p1mup2*p1mup5*F(668)+p1mup2*p2mup5*F(67
     -   0)+p1mup5*p2mup2*F(676)+2*(p4mup2*p4mup5*F(306)+p2mup2*p2mup
     -   5*F(1036)+p2mup5*p4mup2*F(1049)+p2mup2*p4mup5*F(1049))+p1mup
     -   2*p4mup5*F(1387)+p1mup5*p4mup2*F(1394)+p3mup5*(2*(p4mup2*F(3
     -   00)+p3mup2*F(304)+p2mup2*F(1060))+p1mup2*F(1864))+p3mup2*(2*
     -   (p4mup5*F(300)+p2mup5*F(1060))+p1mup5*F(1873)))+mup2mup5*mup
     -   3mup4*F(2716)
       Fa3(720) = mup2mup4*p4mup1*F(2679)+mup1mup4*p4mup2*F(2714)
       Fa3(721) = 4*mup1mup4*(p1mup2*p1mup5*F(668)+p1mup2*p2mup5*F(670
     -   )+p1mup5*p2mup2*F(676)+2*(p4mup2*p4mup5*F(306)+p2mup2*p2mup5
     -   *F(1036)+p2mup5*p4mup2*F(1049)+p2mup2*p4mup5*F(1049))+p1mup2
     -   *p4mup5*F(1387)+p1mup5*p4mup2*F(1394)+p3mup5*(2*(p4mup2*F(30
     -   0)+p3mup2*F(304)+p2mup2*F(1060))+p1mup2*F(1864))+p3mup2*(2*(
     -   p4mup5*F(300)+p2mup5*F(1060))+p1mup5*F(1873)))+mup1mup4*mup2
     -   mup5*F(2720)
       Fa3(722) = -(mup1mup4*(-(mup2mup5*F(2412))+4*(p1mup2*p3mup5*F(1
     -   325)-p2mup2*p3mup5*F(1335)+p3mup2*(p1mup5*F(276)+p2mup5*F(35
     -   9)+p4mup5*F(1012)+p3mup5*F(1843))+p3mup5*p4mup2*F(2718))))
       Fa3(723) = -(mup1mup3*(-(mup2mup5*F(2695))+4*(p1mup5*p4mup2*F(2
     -   76)+p2mup5*p4mup2*F(359)+p3mup5*p4mup2*F(888)+p1mup2*p4mup5*
     -   F(1325)+p4mup2*p4mup5*F(1334)-p2mup2*p4mup5*F(1335)+p3mup2*p
     -   4mup5*F(2717))))
       Fa3(724) = mup3mup4*(mup2mup5*F(2669)+4*(p1mup2*p3mup5*F(1325)-
     -   p2mup2*p3mup5*F(1335)+p3mup2*(p1mup5*F(276)+p2mup5*F(359)+p4
     -   mup5*F(1012)+p3mup5*F(1843))+p3mup5*p4mup2*F(2718)))
       Fa3(725) = -4*mup1mup3*(p1mup2*p1mup5*F(668)+p1mup2*p2mup5*F(67
     -   0)+p1mup5*p2mup2*F(676)+2*(p4mup2*p4mup5*F(306)+p2mup2*p2mup
     -   5*F(1036)+p2mup5*p4mup2*F(1049)+p2mup2*p4mup5*F(1049))+p1mup
     -   2*p4mup5*F(1387)+p1mup5*p4mup2*F(1394)+p3mup5*(2*(p4mup2*F(3
     -   00)+p3mup2*F(304)+p2mup2*F(1060))+p1mup2*F(1864))+p3mup2*(2*
     -   (p4mup5*F(300)+p2mup5*F(1060))+p1mup5*F(1873)))+mup1mup3*mup
     -   2mup5*F(2719)
       Fa3(726) = mup1mup3*p4mup2*F(2706)+mup2mup3*p4mup1*F(2709)
       Fa3(727) = -4*mup3mup4*(p1mup2*(2*(p1mup5*F(130)+p3mup5*F(278)+
     -   p4mup5*F(333))+p2mup5*F(962))+p2mup2*p2mup5*F(1022)+p1mup5*(
     -   p2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028)))+2*(p4mup2*p4
     -   mup5*F(1015)+p3mup5*p4mup2*F(1047)+p3mup2*p4mup5*F(1047)+p3m
     -   up2*p3mup5*F(1059))+p2mup2*p4mup5*F(1387)+p2mup5*p4mup2*F(13
     -   94)+p2mup2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))+mup2mup5*mu
     -   p3mup4*F(2722)
       Fa3(728) = mup1mup2*mup3mup4*F(2671)+mup1mup3*mup2mup4*F(2679)+
     -   mup1mup4*mup2mup3*F(2709)
       Fa3(729) = mup2mup5*(p3mup1*F(1051)+p4mup1*F(7441))
       Fa3(730) = p4mup2*(p2mup1*F(1622)-2*p4mup1*F(7444))
       Fa3(731) = p4mup2*(p1mup1*F(1622)+2*p4mup1*F(7444))
       Fa3(732) = (p1mup1+p2mup1)*p4mup2
       Fa3(733) = mup1mup3*p4mup2*F(1622)-mup2mup3*p4mup1*F(7441)
       Fa3(734) = -(mup1mup4*p4mup2*F(1622))+mup2mup4*p4mup1*F(7441)
       Fa3(735) = mup1mup4*mup2mup3-mup1mup3*mup2mup4+mup1mup2*mup3mup
     -   4
       Fa3(736) = mup2mup5*(p3mup1*F(7440)+p4mup1*F(7441))
       Fa3(737) = mup2mup5*(p3mup1*F(1051)+p4mup1*F(1095))
       Fa3(738) = p4mup2*(2*p4mup1*F(436)-p2mup1*F(1645))
       Fa3(739) = p4mup2*(2*p4mup1*F(436)+p1mup1*F(1645))
       Fa3(740) = mup2mup3*p4mup1*F(1095)-mup1mup3*p4mup2*F(1645)
       Fa3(741) = -(mup2mup4*p4mup1*F(1095))+mup1mup4*p4mup2*F(1645)
       Fa3(742) = -2*p4mup1*p4mup2*F(436)+p2mup1*p4mup2*F(1645)
       Fa3(743) = -(mup2mup3*p4mup1*F(1095))+mup1mup3*p4mup2*F(1645)
       Fa3(744) = mup2mup4*p4mup1*F(1095)-mup1mup4*p4mup2*F(1645)
       Fa3(745) = -(p2mup1*p4mup2*F(1622))+2*p4mup1*p4mup2*F(7444)
       Fa3(746) = -(mup1mup3*p4mup2*F(1622))+mup2mup3*p4mup1*F(7441)
       Fa3(747) = mup1mup4*p4mup2*F(1622)-mup2mup4*p4mup1*F(7441)
       Fa3(748) = mup2mup5*p4mup1*F(2412)+mup2mup5*p3mup1*F(2699)+4*(p
     -   1mup5*(-(p3mup2*p4mup1)+p3mup1*p4mup2)*F(276)+p2mup5*(-(p3mu
     -   p2*p4mup1)+p3mup1*p4mup2)*F(359)+p3mup1*p3mup5*p4mup2*F(888)
     -   -p3mup2*p4mup1*p4mup5*F(1012)-p1mup2*p3mup5*p4mup1*F(1325)+p
     -   1mup2*p3mup1*p4mup5*F(1325)+p3mup1*p4mup2*p4mup5*F(1334)+p2m
     -   up2*p3mup5*p4mup1*F(1335)-p2mup2*p3mup1*p4mup5*F(1335)-p3mup
     -   2*p3mup5*p4mup1*F(1843)+p3mup1*p3mup2*p4mup5*F(2717)-p3mup5*
     -   p4mup1*p4mup2*F(2718))
       Fa3(749) = 4*p4mup1*(p1mup2*p1mup5*F(668)+p1mup2*p2mup5*F(670)+
     -   p1mup5*p2mup2*F(676)+2*(p4mup2*p4mup5*F(306)+p2mup2*p2mup5*F
     -   (1036)+p2mup5*p4mup2*F(1049)+p2mup2*p4mup5*F(1049))+p1mup2*p
     -   4mup5*F(1387)+p1mup5*p4mup2*F(1394)+p3mup5*(2*(p4mup2*F(300)
     -   +p3mup2*F(304)+p2mup2*F(1060))+p1mup2*F(1864))+p3mup2*(2*(p4
     -   mup5*F(300)+p2mup5*F(1060))+p1mup5*F(1873)))+mup2mup5*p4mup1
     -   *F(2720)
       Fa3(750) = -4*p3mup1*(p1mup2*p1mup5*F(668)+p1mup2*p2mup5*F(670)
     -   +p1mup5*p2mup2*F(676)+2*(p4mup2*p4mup5*F(306)+p2mup2*p2mup5*
     -   F(1036)+p2mup5*p4mup2*F(1049)+p2mup2*p4mup5*F(1049))+p1mup2*
     -   p4mup5*F(1387)+p1mup5*p4mup2*F(1394)+p3mup5*(2*(p4mup2*F(300
     -   )+p3mup2*F(304)+p2mup2*F(1060))+p1mup2*F(1864))+p3mup2*(2*(p
     -   4mup5*F(300)+p2mup5*F(1060))+p1mup5*F(1873)))+mup2mup5*p3mup
     -   1*F(2719)
       Fa3(751) = p4mup2*(p2mup1*F(1673)+2*p4mup1*F(1776))
       Fa3(752) = p4mup2*(p1mup1*F(1317)+2*p4mup1*F(1704))
       Fa3(753) = p4mup2*(p1mup1*F(2342)+p2mup1*F(2706)+p4mup1*F(2707)
     -   )
       Fa3(754) = p4mup2*(p2mup1*F(1382)+p1mup1*F(1517))
       Fa3(755) = -4*(p3mup1*(p1mup2*(2*(p1mup5*F(130)+p3mup5*F(278)+p
     -   4mup5*F(333))+p2mup5*F(962))+p2mup2*p2mup5*F(1022)+p1mup5*(p
     -   2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028)))+2*(p4mup2*p4m
     -   up5*F(1015)+p3mup5*p4mup2*F(1047)+p3mup2*p4mup5*F(1047)+p3mu
     -   p2*p3mup5*F(1059))+p2mup2*p4mup5*F(1387)+p2mup5*p4mup2*F(139
     -   4)+p2mup2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))-p2mup1*p3mup
     -   5*p4mup2*F(2718))+mup2mup5*p3mup1*F(2722)
       Fa3(756) = 4*p2mup1*p4mup2*(2*(p3mup5*F(300)+p4mup5*F(306)+p2mu
     -   p5*F(1049))+p1mup5*F(1394))
       Fa3(757) = 4*(-(p2mup1*p2mup5*p4mup2*F(359))-p2mup1*p3mup5*p4mu
     -   p2*F(888)+p1mup2*p4mup1*(2*(p1mup5*F(130)+p3mup5*F(278)+p4mu
     -   p5*F(333))+p2mup5*F(962))+p2mup2*p2mup5*p4mup1*F(1022)+p1mup
     -   5*(-(p2mup1*p4mup2*F(276))+p2mup2*p4mup1*F(964)+2*(p4mup1*p4
     -   mup2*F(338)+p3mup2*p4mup1*F(1028)))+2*(p4mup1*p4mup2*p4mup5*
     -   F(1015)+p3mup5*p4mup1*p4mup2*F(1047)+p3mup2*p4mup1*p4mup5*F(
     -   1047)+p3mup2*p3mup5*p4mup1*F(1059))-p2mup1*p4mup2*p4mup5*F(1
     -   334)+p2mup2*p4mup1*p4mup5*F(1387)+p2mup5*p4mup1*p4mup2*F(139
     -   4)+p2mup2*p3mup5*p4mup1*F(1864)+p2mup5*p3mup2*p4mup1*F(1873)
     -   )+mup2mup5*p4mup1*F(2721)
       Fa3(758) = mup1mup4*p4mup2*F(2342)+mup2mup4*p4mup1*F(2676)
       Fa3(759) = -(mup1mup4*(-(mup2mup5*F(2695))+4*(p1mup5*p4mup2*F(2
     -   76)+p2mup5*p4mup2*F(359)+p3mup5*p4mup2*F(888)+p1mup2*p4mup5*
     -   F(1325)+p4mup2*p4mup5*F(1334)-p2mup2*p4mup5*F(1335)+p3mup2*p
     -   4mup5*F(2717))))
       Fa3(760) = -4*mup1mup4*(p1mup2*(2*(p1mup5*F(130)+p3mup5*F(278)+
     -   p4mup5*F(333))+p2mup5*F(962))+p2mup2*p2mup5*F(1022)+p1mup5*(
     -   p2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028)))+2*(p4mup2*p4
     -   mup5*F(1015)+p3mup5*p4mup2*F(1047)+p3mup2*p4mup5*F(1047)+p3m
     -   up2*p3mup5*F(1059))+p2mup2*p4mup5*F(1387)+p2mup5*p4mup2*F(13
     -   94)+p2mup2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))+mup1mup4*mu
     -   p2mup5*F(2722)
       Fa3(761) = mup1mup3*p4mup2*F(1673)+mup2mup3*p4mup1*F(2696)
       Fa3(762) = mup1mup3*p4mup2*F(1932)+mup2mup3*p4mup1*F(2697)
       Fa3(763) = -(mup1mup3*(-(mup2mup5*F(2412))+4*(p1mup2*p3mup5*F(1
     -   325)-p2mup2*p3mup5*F(1335)+p3mup2*(p1mup5*F(276)+p2mup5*F(35
     -   9)+p4mup5*F(1012)+p3mup5*F(1843))+p3mup5*p4mup2*F(2718))))
       Fa3(764) = 4*mup1mup3*(p1mup2*(2*(p1mup5*F(130)+p3mup5*F(278)+p
     -   4mup5*F(333))+p2mup5*F(962))+p2mup2*p2mup5*F(1022)+p1mup5*(p
     -   2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028)))+2*(p4mup2*p4m
     -   up5*F(1015)+p3mup5*p4mup2*F(1047)+p3mup2*p4mup5*F(1047)+p3mu
     -   p2*p3mup5*F(1059))+p2mup2*p4mup5*F(1387)+p2mup5*p4mup2*F(139
     -   4)+p2mup2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))+mup1mup3*mup
     -   2mup5*F(2721)
       Fa3(765) = mup1mup4*p4mup2*F(1317)+mup2mup4*p4mup1*F(2589)
       Fa3(766) = mup1mup3*mup2mup4*F(2691)+mup1mup2*mup3mup4*F(2692)+
     -   mup1mup4*mup2mup3*F(2699)
       Fa3(767) = mup1mup3*mup2mup4*F(2589)+mup1mup2*mup3mup4*F(2668)+
     -   mup1mup4*mup2mup3*F(2696)
       Fa3(768) = mup1mup2*mup3mup4*F(2413)+mup1mup3*mup2mup4*F(2676)+
     -   mup1mup4*mup2mup3*F(2697)
       Fa3(769) = mup3mup4*(mup2mup5*F(2692)+4*(p1mup5*p4mup2*F(276)+p
     -   2mup5*p4mup2*F(359)+p3mup5*p4mup2*F(888)+p1mup2*p4mup5*F(132
     -   5)+p4mup2*p4mup5*F(1334)-p2mup2*p4mup5*F(1335)+p3mup2*p4mup5
     -   *F(2717)))
       Fa3(770) = 4*mup3mup4*(p1mup2*p1mup5*F(668)+p1mup2*p2mup5*F(670
     -   )+p1mup5*p2mup2*F(676)+2*(p4mup2*p4mup5*F(306)+p2mup2*p2mup5
     -   *F(1036)+p2mup5*p4mup2*F(1049)+p2mup2*p4mup5*F(1049))+p1mup2
     -   *p4mup5*F(1387)+p1mup5*p4mup2*F(1394)+p3mup5*(2*(p4mup2*F(30
     -   0)+p3mup2*F(304)+p2mup2*F(1060))+p1mup2*F(1864))+p3mup2*(2*(
     -   p4mup5*F(300)+p2mup5*F(1060))+p1mup5*F(1873)))+mup2mup5*mup3
     -   mup4*F(2712)
       Fa3(771) = mup2mup4*p4mup1*F(2677)+mup1mup4*p4mup2*F(2706)
       Fa3(772) = -4*mup1mup4*(p1mup2*p1mup5*F(668)+p1mup2*p2mup5*F(67
     -   0)+p1mup5*p2mup2*F(676)+2*(p4mup2*p4mup5*F(306)+p2mup2*p2mup
     -   5*F(1036)+p2mup5*p4mup2*F(1049)+p2mup2*p4mup5*F(1049))+p1mup
     -   2*p4mup5*F(1387)+p1mup5*p4mup2*F(1394)+p3mup5*(2*(p4mup2*F(3
     -   00)+p3mup2*F(304)+p2mup2*F(1060))+p1mup2*F(1864))+p3mup2*(2*
     -   (p4mup5*F(300)+p2mup5*F(1060))+p1mup5*F(1873)))+mup1mup4*mup
     -   2mup5*F(2719)
       Fa3(773) = mup1mup4*(mup2mup5*F(2589)+4*(p1mup2*p3mup5*F(1325)-
     -   p2mup2*p3mup5*F(1335)+p3mup2*(p1mup5*F(276)+p2mup5*F(359)+p4
     -   mup5*F(1012)+p3mup5*F(1843))+p3mup5*p4mup2*F(2718)))
       Fa3(774) = mup1mup3*(mup2mup5*F(2699)+4*(p1mup5*p4mup2*F(276)+p
     -   2mup5*p4mup2*F(359)+p3mup5*p4mup2*F(888)+p1mup2*p4mup5*F(132
     -   5)+p4mup2*p4mup5*F(1334)-p2mup2*p4mup5*F(1335)+p3mup2*p4mup5
     -   *F(2717)))
       Fa3(775) = -(mup3mup4*(-(mup2mup5*F(2668))+4*(p1mup2*p3mup5*F(1
     -   325)-p2mup2*p3mup5*F(1335)+p3mup2*(p1mup5*F(276)+p2mup5*F(35
     -   9)+p4mup5*F(1012)+p3mup5*F(1843))+p3mup5*p4mup2*F(2718))))
       Fa3(776) = 4*mup1mup3*(p1mup2*p1mup5*F(668)+p1mup2*p2mup5*F(670
     -   )+p1mup5*p2mup2*F(676)+2*(p4mup2*p4mup5*F(306)+p2mup2*p2mup5
     -   *F(1036)+p2mup5*p4mup2*F(1049)+p2mup2*p4mup5*F(1049))+p1mup2
     -   *p4mup5*F(1387)+p1mup5*p4mup2*F(1394)+p3mup5*(2*(p4mup2*F(30
     -   0)+p3mup2*F(304)+p2mup2*F(1060))+p1mup2*F(1864))+p3mup2*(2*(
     -   p4mup5*F(300)+p2mup5*F(1060))+p1mup5*F(1873)))+mup1mup3*mup2
     -   mup5*F(2720)
       Fa3(777) = mup2mup3*p4mup1*F(2708)+mup1mup3*p4mup2*F(2714)
       Fa3(778) = 4*mup3mup4*(p1mup2*(2*(p1mup5*F(130)+p3mup5*F(278)+p
     -   4mup5*F(333))+p2mup5*F(962))+p2mup2*p2mup5*F(1022)+p1mup5*(p
     -   2mup2*F(964)+2*(p4mup2*F(338)+p3mup2*F(1028)))+2*(p4mup2*p4m
     -   up5*F(1015)+p3mup5*p4mup2*F(1047)+p3mup2*p4mup5*F(1047)+p3mu
     -   p2*p3mup5*F(1059))+p2mup2*p4mup5*F(1387)+p2mup5*p4mup2*F(139
     -   4)+p2mup2*p3mup5*F(1864)+p2mup5*p3mup2*F(1873))+mup2mup5*mup
     -   3mup4*F(2721)
       Fa3(779) = mup1mup2*mup3mup4*F(2670)+mup1mup3*mup2mup4*F(2677)+
     -   mup1mup4*mup2mup3*F(2708)
       Fa3(780) = mup1mup5*(2*(p1mup2*F(1320)+p2mup2*F(1323))+p3mup2*F
     -   (2279))
       Fa3(781) = mup1mup5*(p1mup2*F(2287)+2*(p2mup2*F(2293)+p3mup2*F(
     -   2296)))+mup2mup5*p1mup1*F(2695)
       Fa3(782) = mup2mup5*p1mup1*F(2589)+mup1mup5*(p1mup2*F(1373)+2*p
     -   2mup2*F(1379)+p3mup2*F(2684))
       Fa3(783) = 2*(mup1mup5*p2mup3*F(1830)+mup1mup5*p1mup3*F(1831))+
     -   mup1mup5*p3mup3*F(2280)+mup1mup5*p4mup3*F(2283)+mup3mup5*p3m
     -   up1*F(2699)+4*(p1mup5*p3mup1*p4mup3*F(276)+p2mup5*p3mup1*p4m
     -   up3*F(359)+p3mup1*p3mup5*p4mup3*F(888)+p3mup1*p4mup3*p4mup5*
     -   F(892)+p2mup3*p3mup1*p4mup5*F(2710))
       Fa3(784) = mup3mup5*p1mup1*F(2699)+mup1mup5*(p1mup3*F(2288)-2*(
     -   p3mup3*F(1861)-p2mup3*F(2294))+p4mup3*F(2702))
       Fa3(785) = -8*(p3mup1*p4mup3*(p3mup5*F(302)+p2mup5*F(678)+p4mup
     -   5*F(912))+p2mup3*p3mup1*(p4mup5*F(300)+p3mup5*F(304)+p2mup5*
     -   F(1060)))-4*p1mup5*p3mup1*(p4mup3*F(298)-p2mup3*F(1863))+mup
     -   1mup5*p1mup3*F(2680)+mup1mup5*p3mup3*F(2682)+mup1mup5*p2mup3
     -   *F(2683)+mup3mup5*p1mup1*F(2700)+2*mup1mup5*p4mup3*F(2703)+m
     -   up3mup5*p3mup1*F(2713)
       Fa3(786) = mup1mup5*p3mup4*F(1669)+2*(mup1mup5*p2mup4*F(1671)+m
     -   up1mup5*p1mup4*F(1672))+mup1mup5*p4mup4*F(1674)-4*(p1mup5*p3
     -   mup4*p4mup1*F(276)+p3mup4*p4mup1*(p3mup5*F(335)+p2mup5*F(359
     -   )+p4mup5*F(1012))-p2mup4*p3mup5*p4mup1*F(2414))
       Fa3(787) = -4*(p3mup4*p4mup1*(-(p1mup5*F(613))+2*(p3mup5*F(302)
     -   +p2mup5*F(678)+p4mup5*F(912)))+p2mup4*p4mup1*(-(p1mup5*F(672
     -   ))+2*(p3mup5*F(300)+p4mup5*F(306)+p2mup5*F(1049))))+mup1mup5
     -   *p1mup4*F(2402)+mup1mup5*p2mup4*F(2403)+mup1mup5*p3mup4*F(24
     -   04)+mup1mup5*p4mup4*F(2405)
       Fa3(788) = mup1mup5*(p1mup4*F(1697)+2*(p2mup4*F(1380)+p4mup4*F(
     -   1705))+p3mup4*F(2588))
       Fa3(789) = p1mup1*(p2mup5*F(1322)+p3mup5*F(2281))
       Fa3(790) = p1mup1*(p2mup5*F(2340)+2*p3mup5*F(888)*P(7))
       Fa3(791) = p1mup1*(p3mup5*F(1273)+p2mup5*F(1512))
       Fa3(792) = mup1mup5*(p3mup4*F(1005)+p1mup4*F(1007)+p4mup4*F(101
     -   1)+p2mup4*F(1770))
       Fa3(793) = mup1mup5*(p1mup4*F(897)+2*(p3mup4*F(903)+p4mup4*F(90
     -   6))+p2mup4*F(2587))
       Fa3(794) = mup1mup5*p1mup3*F(1003)+mup1mup5*p3mup3*F(1004)-2*mu
     -   p1mup5*p4mup3*F(1009)+p2mup3*(mup1mup5*F(1000)+4*p3mup1*(p3m
     -   up5*F(888)+p4mup5*F(1012)+p2mup5*F(1554)))+mup3mup5*p3mup1*F
     -   (2590)+4*(p1mup5*p2mup3*p3mup1*F(276)+p2mup5*p3mup1*p4mup3*F
     -   (2711))
       Fa3(795) = mup3mup5*p1mup1*F(2413)+mup1mup5*(p1mup3*F(896)-2*(p
     -   3mup3*F(902)+p4mup3*F(905))+p2mup3*F(2687))
       Fa3(796) = p2mup1*(p3mup4*p4mup3*F(1014)+p2mup3*p3mup4*F(2724)+
     -   p3mup3*(p3mup4*F(331)-p4mup4*F(1278)+p2mup4*F(2725))-2*(p1mu
     -   p4*p4mup3*F(336)+p1mup3*p4mup4*F(336)+p2mup3*p4mup4*F(337)-p
     -   1mup3*p1mup4*F(960)-p4mup3*p4mup4*F(1018)-p2mup3*p2mup4*F(27
     -   23)+p2mup4*p4mup3*F(2726)))
       Fa3(797) = p2mup1*(-(p1mup3*p4mup4*F(615))+p1mup4*(p4mup3*F(297
     -   )+p1mup3*F(908))+2*(p3mup4*p4mup3*F(910)+p4mup3*p4mup4*F(914
     -   )+p2mup3*p4mup4*F(2727)-p3mup3*(p3mup4*F(303)-p4mup4*F(913)+
     -   p2mup4*F(2728))+p2mup3*p3mup4*F(2729)-p2mup4*p4mup3*F(2730)+
     -   p2mup3*p2mup4*F(2731)))
       Fa3(798) = mup1mup5*(p2mup2*F(1510)+2*p1mup2*F(1513)+p3mup2*F(1
     -   514))-4*(p1mup5*p2mup1*p3mup2*F(276)+p2mup1*p3mup2*(p2mup5*F
     -   (359)+p4mup5*F(1012)+p3mup5*F(1843)))+mup2mup5*p2mup1*F(2412
     -   )
       Fa3(799) = -4*p2mup1*p3mup2*(2*(p4mup5*F(300)+p3mup5*F(304)+p2m
     -   up5*F(1060))+p1mup5*F(1873))+mup1mup5*(p1mup2*F(2685)+p2mup2
     -   *F(2688)+p3mup2*F(2689))+mup2mup5*p2mup1*F(2719)+mup2mup5*p1
     -   mup1*F(2722)
       Fa3(800) = mup1mup5*(p2mup4*F(878)-2*(p1mup4*F(880)+p3mup4*F(88
     -   1))+p4mup4*F(883))+4*(p1mup5*p2mup4*p4mup1*F(276)+p4mup1*(p2
     -   mup4*p3mup5*F(888)+p2mup4*p4mup5*F(1012)+p2mup4*p2mup5*F(151
     -   9)-p2mup5*p3mup4*F(2415)))
       Fa3(801) = 4*mup3mup4*p4mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mup
     -   5*F(888)+p4mup5*F(1012))+mup1mup5*mup3mup4*F(2692)
       Fa3(802) = -4*mup3mup4*p3mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mu
     -   p5*F(888)+p4mup5*F(1012))+mup1mup5*mup3mup4*F(2668)
       Fa3(803) = mup1mup5*p2mup2*F(2338)+mup1mup5*p1mup2*F(2341)+2*mu
     -   p1mup5*p3mup2*F(2345)+mup2mup5*p2mup1*F(2699)+4*p2mup1*p3mup
     -   2*p4mup5*F(2717)
       Fa3(804) = -4*mup2mup4*p4mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mu
     -   p5*F(888)+p4mup5*F(1012))+mup1mup5*mup2mup4*F(2691)
       Fa3(805) = -4*mup2mup4*p2mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mu
     -   p5*F(888)+p4mup5*F(1012))+mup1mup5*mup2mup4*F(2676)
       Fa3(806) = -4*mup2mup3*p3mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mu
     -   p5*F(888)+p4mup5*F(1012))+mup1mup5*mup2mup3*F(2696)
       Fa3(807) = 4*mup2mup3*p2mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mup
     -   5*F(888)+p4mup5*F(1012))+mup1mup5*mup2mup3*F(2697)
       Fa3(808) = p2mup1*(p1mup4*p4mup3*F(274)-p1mup3*(p4mup4*F(272)+2
     -   *p1mup4*F(886))+p3mup4*p4mup3*F(889)+p4mup3*p4mup4*F(1233)+p
     -   3mup3*p4mup4*F(1234)-p2mup3*p4mup4*F(2733)+p2mup4*p4mup3*F(2
     -   734)+2*(p2mup4*p3mup3*F(285)+p3mup3*p3mup4*F(286)-p2mup3*p2m
     -   up4*F(2732)-p2mup3*p3mup4*F(2735)))
       Fa3(809) = mup1mup5*(p4mup3*F(882)+p2mup3*F(2339)+2*(p1mup3*F(8
     -   85)+p3mup3*F(894))*P(7))
       Fa3(810) = 4*mup3mup4*p2mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mup
     -   5*F(888)+p4mup5*F(1012))+mup1mup5*mup3mup4*F(2413)
       Fa3(811) = 4*mup2mup4*p3mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mup
     -   5*F(888)+p4mup5*F(1012))+mup1mup5*mup2mup4*F(2589)
       Fa3(812) = 4*mup2mup3*p4mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mup
     -   5*F(888)+p4mup5*F(1012))+mup1mup5*mup2mup3*F(2699)
       Fa3(813) = (mup1mup4*mup2mup3-mup1mup3*mup2mup4+mup1mup2*mup3mu
     -   p4)*(p1mup5*F(276)+p2mup5*F(359)+p3mup5*F(888)+p4mup5*F(1012
     -   ))
       Fa3(814) = mup1mup5*(p3mup2*F(1621)+2*(p1mup2*F(1325)+p2mup2*F(
     -   1623)))
       Fa3(815) = mup1mup5*(2*(p2mup2*F(361)+p3mup2*F(362))+p1mup2*F(3
     -   67))+mup2mup5*p1mup1*F(1051)
       Fa3(816) = mup2mup5*p1mup1*F(7441)+mup1mup5*(p1mup2*F(352)+2*(p
     -   3mup2*F(7442)+p2mup2*F(7448)))
       Fa3(817) = mup1mup5*(p3mup3*F(351)+p4mup3*F(353)+2*(p1mup3*F(34
     -   9)+p2mup3*F(1050)))-mup3mup5*p3mup1*F(1051)
       Fa3(818) = mup1mup5*(2*(p4mup3*F(335)+p2mup3*F(361)+p3mup3*F(36
     -   2))+p1mup3*F(367))+mup3mup5*p1mup1*F(1051)
       Fa3(819) = mup3mup5*p1mup1*F(7441)+mup1mup5*(p1mup3*F(352)+2*(p
     -   3mup3*F(7442)+p4mup3*F(7444)+p2mup3*F(7448)))
       Fa3(820) = mup1mup5*(p3mup4*F(350)+p4mup4*F(352)+2*(p1mup4*F(27
     -   6)+p2mup4*F(359)))
       Fa3(821) = mup1mup5*(2*(p4mup4*F(335)+p2mup4*F(361)+p3mup4*F(36
     -   2))+p1mup4*F(367))
       Fa3(822) = mup1mup5*(p1mup4*F(352)+2*(p3mup4*F(7442)+p4mup4*F(7
     -   444)+p2mup4*F(7448)))
       Fa3(823) = p1mup1*(p3mup5*F(350)+2*p2mup5*F(359))
       Fa3(824) = mup2mup5*p2mup1*F(7441)-mup1mup5*(p2mup2*F(1626)+2*(
     -   p1mup2*F(1012)+p3mup2*F(7442)))
       Fa3(825) = -(mup2mup5*p2mup1*F(1051))+mup1mup5*(p2mup2*F(1624)+
     -   2*(p3mup2*F(362)+p1mup2*F(1625)))
       Fa3(826) = mup1mup5*(2*(p1mup2*F(1325)+p2mup2*F(1623))+p3mup2*F
     -   (1643))
       Fa3(827) = mup2mup5*p1mup1*F(7440)+mup1mup5*(p1mup2*F(402)+2*(p
     -   3mup2*F(7443)+p2mup2*F(7447)))
       Fa3(828) = mup1mup5*(p4mup3*F(352)+2*(p1mup3*F(276)+p2mup3*F(35
     -   9))+p3mup3*F(402))-mup3mup5*p3mup1*F(7440)
       Fa3(829) = mup3mup5*p1mup1*F(7440)+mup1mup5*(p1mup3*F(402)+2*(p
     -   4mup3*F(7442)+p3mup3*F(7443)+p2mup3*F(7447)))
       Fa3(830) = mup1mup5*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(35
     -   9))+p3mup4*F(402))
       Fa3(831) = mup1mup5*(p1mup4*F(402)+2*(p4mup4*F(7442)+p3mup4*F(7
     -   443)+p2mup4*F(7447)))
       Fa3(832) = p1mup1*(2*p2mup5*F(359)+p3mup5*F(402))
       Fa3(833) = mup2mup5*p2mup1*F(7440)-mup1mup5*(p2mup2*F(1644)+2*(
     -   p1mup2*F(888)+p3mup2*F(7443)))
       Fa3(834) = mup1mup5*(-(p3mup2*F(1621))+2*(p2mup2*F(1335)+p1mup2
     -   *F(1646)))
       Fa3(835) = mup1mup5*(p1mup2*F(438)+2*(p2mup2*F(434)+p3mup2*F(43
     -   9)))+mup2mup5*p1mup1*F(1095)
       Fa3(836) = mup1mup5*(p3mup3*F(351)+p4mup3*F(430)+2*(p1mup3*F(34
     -   9)+p2mup3*F(1050)))-mup3mup5*p3mup1*F(1051)
       Fa3(837) = mup1mup5*(p1mup3*F(438)+2*(p2mup3*F(434)+p4mup3*F(43
     -   6)+p3mup3*F(439)))+mup3mup5*p1mup1*F(1095)
       Fa3(838) = mup1mup5*(p3mup4*F(351)+p4mup4*F(430)+2*(p1mup4*F(34
     -   9)+p2mup4*F(1050)))
       Fa3(839) = mup1mup5*(p1mup4*F(438)+2*(p2mup4*F(434)+p4mup4*F(43
     -   6)+p3mup4*F(439)))
       Fa3(840) = p1mup1*(-2*p2mup5*F(359)+p3mup5*F(370))
       Fa3(841) = -(mup2mup5*p2mup1*F(1095))+mup1mup5*(p2mup2*F(1653)+
     -   2*(p3mup2*F(439)+p1mup2*F(1654)))
       Fa3(842) = mup1mup5*(-(p3mup2*F(1643))+2*(p2mup2*F(1335)+p1mup2
     -   *F(1646)))
       Fa3(843) = -(mup2mup5*p2mup1*F(7441))+mup1mup5*(p2mup2*F(1626)+
     -   2*(p1mup2*F(1012)+p3mup2*F(7442)))
       Fa3(844) = -(mup2mup5*p2mup1*F(7440))+mup1mup5*(p2mup2*F(1644)+
     -   2*(p1mup2*F(888)+p3mup2*F(7443)))
       Fa3(845) = mup1mup5*(p4mup4*F(352)+2*(p1mup4*F(276)+p2mup4*F(35
     -   9))-p3mup4*F(370))
       Fa3(846) = mup1mup5*(2*p1mup2*F(1676)+p2mup2*F(1679)+p3mup2*F(1
     -   832))
       Fa3(847) = mup1mup5*(p1mup2*F(1851)+2*(p2mup2*F(1857)+p3mup2*F(
     -   1860)))+mup2mup5*p1mup1*F(2699)
       Fa3(848) = mup2mup5*p1mup1*F(2412)+mup1mup5*(p1mup2*F(1696)+2*p
     -   2mup2*F(1702)+p3mup2*F(2674))
       Fa3(849) = -2*(mup1mup5*p2mup3*F(1830)+mup1mup5*p1mup3*F(1831))
     -   +mup1mup5*p3mup3*F(1833)+mup1mup5*p4mup3*F(1836)+mup3mup5*p3
     -   mup1*F(2695)-4*(p1mup5*p3mup1*p4mup3*F(276)+p2mup5*p3mup1*p4
     -   mup3*F(359)+p3mup1*p3mup5*p4mup3*F(888)+p3mup1*p4mup3*p4mup5
     -   *F(892)+p2mup3*p3mup1*p4mup5*F(2710))
       Fa3(850) = mup3mup5*p1mup1*F(2695)+mup1mup5*(p1mup3*F(1852)+2*(
     -   p2mup3*F(1858)+p3mup3*F(1861))+p4mup3*F(2704))
       Fa3(851) = 8*(p3mup1*p3mup5*p4mup3*F(302)+p2mup5*p3mup1*p4mup3*
     -   F(678)+p3mup1*p4mup3*p4mup5*F(912)+p2mup3*p3mup1*(p4mup5*F(3
     -   00)+p3mup5*F(304)+p2mup5*F(1060)))+4*p1mup5*p3mup1*(p4mup3*F
     -   (298)+p2mup3*F(2298))+mup1mup5*p1mup3*F(2664)+mup1mup5*p2mup
     -   3*F(2665)+mup1mup5*p3mup3*F(2667)+mup3mup5*p3mup1*F(2701)-2*
     -   mup1mup5*p4mup3*F(2703)+mup3mup5*p1mup1*F(2705)
       Fa3(852) = mup1mup5*p3mup4*F(1313)+2*(mup1mup5*p2mup4*F(1315)+m
     -   up1mup5*p1mup4*F(1316))+mup1mup5*p4mup4*F(1318)+4*(p1mup5*p3
     -   mup4*p4mup1*F(276)+p3mup4*p4mup1*(p3mup5*F(335)+p2mup5*F(359
     -   )+p4mup5*F(1012))-p2mup4*p3mup5*p4mup1*F(2414))
       Fa3(853) = 8*(p3mup4*p4mup1*(p3mup5*F(302)+p2mup5*F(678)+p4mup5
     -   *F(912))+p2mup4*p4mup1*(p3mup5*F(300)+p4mup5*F(306)+p2mup5*F
     -   (1049)))+4*p1mup5*p4mup1*(p3mup4*F(291)+p2mup4*F(1102))+mup1
     -   mup5*p1mup4*F(2406)+mup1mup5*p2mup4*F(2407)+mup1mup5*p3mup4*
     -   F(2408)+mup1mup5*p4mup4*F(2409)
       Fa3(854) = mup1mup5*(p1mup4*F(1374)-2*(p2mup4*F(1380)-p4mup4*F(
     -   1383))+p3mup4*F(2411))
       Fa3(855) = p1mup1*(2*p2mup5*F(1678)+p3mup5*F(1834))
       Fa3(856) = p1mup1*(p2mup5*F(1930)-2*p3mup5*F(888)*P(7))
       Fa3(857) = p1mup1*(p3mup5*F(1006)+p2mup5*F(1771))
       Fa3(858) = mup1mup5*(p3mup4*F(1272)+2*p1mup4*F(1274)+p4mup4*F(1
     -   277)+p2mup4*F(1511))
       Fa3(859) = mup1mup5*(p1mup4*F(1237)-2*(p4mup4*F(906)-p3mup4*F(1
     -   242))+p2mup4*F(2410))
       Fa3(860) = 2*(mup1mup5*p4mup3*F(1009)+mup1mup5*p1mup3*F(1270))+
     -   mup1mup5*p3mup3*F(1271)+p2mup3*(mup1mup5*F(1267)-4*p3mup1*(p
     -   3mup5*F(888)+p4mup5*F(1012)+p2mup5*F(1554)))+mup3mup5*p3mup1
     -   *F(2413)-4*(p1mup5*p2mup3*p3mup1*F(276)+p2mup5*p3mup1*p4mup3
     -   *F(2711))
       Fa3(861) = mup3mup5*p1mup1*F(2590)+mup1mup5*(2*(p3mup3*F(902)+p
     -   4mup3*F(905))+p1mup3*F(1236)+p2mup3*F(2666))
       Fa3(862) = p2mup1*(p3mup3*p4mup4*F(1017)+4*(p3mup3*p3mup4*F(331
     -   )+p3mup4*p4mup3*F(1014)+p2mup3*p3mup4*F(2724)+p2mup4*p3mup3*
     -   F(2725))-8*(p1mup4*p4mup3*F(336)+p1mup3*p4mup4*F(336)+p2mup3
     -   *p4mup4*F(337)-p1mup3*p1mup4*F(960)-p4mup3*p4mup4*F(1018)-p2
     -   mup3*p2mup4*F(2723)+p2mup4*p4mup3*F(2726)))
       Fa3(863) = p2mup1*(4*p1mup4*p4mup3*F(297)+p1mup3*(p4mup4*F(296)
     -   +4*p1mup4*F(908))+8*(-(p3mup3*p3mup4*F(303))+p3mup4*p4mup3*F
     -   (910)+p3mup3*p4mup4*F(913)+p4mup3*p4mup4*F(914)+p2mup3*p4mup
     -   4*F(2727)-p2mup4*p3mup3*F(2728)+p2mup3*p3mup4*F(2729)-p2mup4
     -   *p4mup3*F(2730)+p2mup3*p2mup4*F(2731)))
       Fa3(864) = mup1mup5*(p2mup2*F(1769)+2*p1mup2*F(1772))+p3mup2*(m
     -   up1mup5*F(1773)+4*p2mup1*(p1mup5*F(276)+p2mup5*F(359)+p4mup5
     -   *F(1012)+p3mup5*F(1843)))+mup2mup5*p2mup1*F(2589)
       Fa3(865) = 4*p2mup1*p3mup2*(2*(p4mup5*F(300)+p3mup5*F(304)+p2mu
     -   p5*F(1060))+p1mup5*F(1873))+mup1mup5*(p1mup2*F(2672)+p2mup2*
     -   F(2673)+p3mup2*F(2675))+mup2mup5*p2mup1*F(2720)+mup2mup5*p1m
     -   up1*F(2721)
       Fa3(866) = mup1mup5*(p2mup4*F(1227)+2*(p1mup4*F(880)+p3mup4*F(1
     -   229))+p4mup4*F(1231))-4*(p1mup5*p2mup4*p4mup1*F(276)+p2mup4*
     -   p4mup1*(p3mup5*F(888)+p4mup5*F(1012)+p2mup5*F(1519))-p2mup5*
     -   p3mup4*p4mup1*F(2415))
       Fa3(867) = -4*mup3mup4*p4mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mu
     -   p5*F(888)+p4mup5*F(1012))+mup1mup5*mup3mup4*F(2690)
       Fa3(868) = 4*mup3mup4*p3mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mup
     -   5*F(888)+p4mup5*F(1012))+mup1mup5*mup3mup4*F(2669)
       Fa3(869) = mup1mup5*p2mup2*F(1928)+2*(mup1mup5*p1mup2*F(1931)+m
     -   up1mup5*p3mup2*F(1935))+mup2mup5*p2mup1*F(2695)-4*p2mup1*p3m
     -   up2*p4mup5*F(2717)
       Fa3(870) = 4*mup2mup4*p4mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mup
     -   5*F(888)+p4mup5*F(1012))+mup1mup5*mup2mup4*F(2693)
       Fa3(871) = 4*mup2mup4*p2mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mup
     -   5*F(888)+p4mup5*F(1012))+mup1mup5*mup2mup4*F(2678)
       Fa3(872) = 4*mup2mup3*p3mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mup
     -   5*F(888)+p4mup5*F(1012))+mup1mup5*mup2mup3*F(2681)
       Fa3(873) = -4*mup2mup3*p2mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mu
     -   p5*F(888)+p4mup5*F(1012))+mup1mup5*mup2mup3*F(2686)
       Fa3(874) = p2mup1*(p4mup3*p4mup4*F(893)+p3mup3*p4mup4*F(895)+4*
     -   (p1mup3*p4mup4*F(272)-p1mup4*p4mup3*F(274)-p3mup4*p4mup3*F(8
     -   89)-2*(p2mup4*p3mup3*F(285)+p3mup3*p3mup4*F(286)-p2mup3*p2mu
     -   p4*F(2732))+p2mup3*p4mup4*F(2733)-p2mup4*p4mup3*F(2734))+8*(
     -   p1mup3*p1mup4*F(886)+p2mup3*p3mup4*F(2735)))
       Fa3(875) = mup1mup5*(p4mup3*F(1230)+p2mup3*F(1929)+2*(p1mup3*F(
     -   271)+p3mup3*F(1235))*P(7))
       Fa3(876) = -4*mup3mup4*p2mup1*(p1mup5*F(276)+p2mup5*F(359)+p3mu
     -   p5*F(888)+p4mup5*F(1012))+mup1mup5*mup3mup4*F(2590)
       Return
       End
