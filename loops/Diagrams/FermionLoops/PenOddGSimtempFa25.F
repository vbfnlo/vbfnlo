c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim25(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa25)
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
       Complex*16 Fa25(2377:2475)
       COMMON/FaPenOddGSim/Fa
       Fa25(2377) = mup1mup2*F(450)+2*(Fa(600)+Fa(737))
       Fa25(2378) = -4*p2mup1*p4mup2*F(356)+mup1mup2*F(412)-2*Fa(601)
       Fa25(2379) = -4*p2mup1*p4mup2*F(459)+mup1mup2*F(460)-2*Fa(628)
       Fa25(2380) = -(mup1mup2*p1mup3*F(443))+2*mup1mup3*p4mup2*F(459)
     -   +Fa(629)
       Fa25(2381) = -(mup1mup5*mup2mup3)+Fa(738)
       Fa25(2382) = 2*mup1mup3*p4mup2*F(356)-mup1mup2*p1mup3*F(380)+Fa
     -   (630)
       Fa25(2383) = 2*mup1mup2*p1mup3*F(376)+Fa(633)+Fa(739)
       Fa25(2384) = 2*mup1mup3*p4mup2*F(459)+Fa(637)
       Fa25(2385) = -Fa(40)+Fa(160)
       Fa25(2386) = Fa(41)+2*Fa(161)
       Fa25(2387) = Fa(56)+Fa(166)
       Fa25(2388) = mup3mup4*F(488)+2*Fa(165)
       Fa25(2389) = mup3mup4*F(478)+2*Fa(167)
       Fa25(2390) = mup3mup4*F(480)+2*Fa(168)
       Fa25(2391) = mup2mup4*F(473)+2*Fa(277)
       Fa25(2392) = mup2mup3*F(475)+2*(Fa(269)+Fa(360))
       Fa25(2393) = mup2mup3*F(469)+2*(Fa(244)+Fa(368))
       Fa25(2394) = mup2mup3*F(476)+2*Fa(271)+4*Fa(362)
       Fa25(2395) = mup2mup4*F(472)-2*Fa(273)
       Fa25(2396) = mup2mup4*F(474)-2*Fa(274)
       Fa25(2397) = mup1mup5*p4mup2*F(1623)+Fa(366)-Fa(735)+Fa(814)
       Fa25(2398) = -2*mup1mup5*p4mup2*F(459)+Fa(357)-Fa(833)
       Fa25(2399) = mup3mup5*p4mup1*F(477)+Fa(367)-Fa(834)
       Fa25(2400) = -Fa(359)+Fa(835)
       Fa25(2401) = -(mup4mup5*(p3mup1*F(396)+p4mup1*F(477)))+Fa(496)+
     -   Fa(836)
       Fa25(2402) = Fa(481)-Fa(837)
       Fa25(2403) = -2*mup1mup5*p4mup2*F(459)-mup1mup2*p2mup5*F(477)+F
     -   a(840)
       Fa25(2404) = -2*p2mup5*p4mup1*F(477)+Fa(497)
       Fa25(2405) = Fa(498)+2*Fa(624)
       Fa25(2406) = mup1mup4*F(479)-2*Fa(486)
       Fa25(2407) = mup1mup4*F(478)+2*Fa(488)
       Fa25(2408) = mup1mup3*F(446)-2*(Fa(455)+Fa(625))
       Fa25(2409) = 4*p3mup1*p4mup3*F(356)+mup1mup3*F(450)+2*(Fa(456)+
     -   Fa(626))
       Fa25(2410) = mup1mup4*F(480)+2*Fa(489)
       Fa25(2411) = -4*p2mup1*p4mup3*F(459)+mup1mup3*F(460)-2*Fa(490)
       Fa25(2412) = 2*mup1mup4*p4mup2*F(459)-Fa(492)+Fa(631)
       Fa25(2413) = mup3mup4*p4mup1*F(477)+Fa(493)-Fa(639)
       Fa25(2414) = Fa(495)+Fa(610)+Fa(742)
       Fa25(2415) = 2*mup1mup4*p4mup2*F(459)+mup1mup2*p2mup4*F(477)+Fa
     -   (641)
       Fa25(2416) = mup1mup2*F(478)+2*(Fa(599)-Fa(736))
       Fa25(2417) = mup1mup2*F(479)-2*(Fa(600)+Fa(737))
       Fa25(2418) = 4*p2mup1*p4mup2*F(459)+mup1mup2*F(480)+2*Fa(628)
       Fa25(2419) = mup1mup2*p1mup3*F(443)-2*mup1mup3*p4mup2*F(459)+Fa
     -   (638)
       Fa25(2420) = 2*mup1mup2*p1mup3*F(376)+Fa(640)-Fa(741)
       Fa25(2421) = Fa(19)-2*Fa(173)
       Fa25(2422) = Fa(65)+Fa(180)
       Fa25(2423) = mup3mup4*F(492)+2*Fa(178)
       Fa25(2424) = mup3mup4*F(438)+2*Fa(179)
       Fa25(2425) = mup3mup4*F(436)+2*Fa(181)
       Fa25(2426) = mup3mup4*F(413)+2*Fa(182)
       Fa25(2427) = mup2mup4*F(433)+2*Fa(255)
       Fa25(2428) = mup2mup4*F(432)+2*Fa(257)
       Fa25(2429) = mup2mup3*F(434)-2*(Fa(258)+Fa(350))
       Fa25(2430) = mup2mup3*F(430)-2*(Fa(259)+Fa(351))
       Fa25(2431) = mup2mup3*F(429)-2*Fa(260)-4*Fa(352)
       Fa25(2432) = mup2mup3*F(409)-2*Fa(247)-4*Fa(340)
       Fa25(2433) = mup2mup4*F(431)+2*Fa(262)
       Fa25(2434) = mup2mup4*F(407)+2*Fa(250)
       Fa25(2435) = -(mup1mup5*p4mup2*F(1609))+Fa(346)+Fa(732)+Fa(841)
       Fa25(2436) = Fa(349)-Fa(828)
       Fa25(2437) = Fa(336)-Fa(819)
       Fa25(2438) = -Fa(444)+Fa(822)
       Fa25(2439) = mup1mup2*p2mup5*F(4463)+2*mup1mup5*p4mup2*F(4465)+
     -   Fa(842)
       Fa25(2440) = mup1mup2*p2mup5*F(4464)+2*mup1mup5*p4mup2*F(4467)+
     -   Fa(843)
       Fa25(2441) = -Fa(466)+Fa(830)
       Fa25(2442) = 2*p2mup5*p4mup1*F(4464)+Fa(500)
       Fa25(2443) = Fa(501)+2*Fa(595)
       Fa25(2444) = mup1mup4*F(437)+2*Fa(470)
       Fa25(2445) = mup1mup4*F(438)-2*Fa(471)
       Fa25(2446) = 2*mup1mup4*p4mup3*F(4465)-Fa(472)
       Fa25(2447) = mup1mup4*F(436)-2*Fa(474)
       Fa25(2448) = mup1mup3*F(419)+2*(Fa(475)+Fa(614))
       Fa25(2449) = mup1mup3*F(435)-4*p3mup1*p4mup3*F(4465)-2*(Fa(476)
     -   +Fa(598))
       Fa25(2450) = mup1mup3*F(421)+4*p2mup1*p4mup3*F(4465)+2*Fa(477)
       Fa25(2451) = mup1mup4*F(413)-2*Fa(458)
       Fa25(2452) = mup1mup3*F(390)+4*p2mup1*p4mup3*F(4467)+2*Fa(459)
       Fa25(2453) = 2*mup1mup4*p4mup3*F(4467)-Fa(460)
       Fa25(2454) = mup3mup4*p4mup1*F(4464)+Fa(478)+Fa(645)
       Fa25(2455) = -(mup1mup2*p2mup4*F(4463))-2*mup1mup4*p4mup2*F(446
     -   5)+Fa(647)
       Fa25(2456) = -Fa(479)+Fa(648)+Fa(744)
       Fa25(2457) = -(mup1mup2*p2mup4*F(4464))-2*mup1mup4*p4mup2*F(446
     -   7)+Fa(649)
       Fa25(2458) = mup1mup2*F(436)+2*(Fa(642)+Fa(733))
       Fa25(2459) = mup1mup2*F(437)+2*(Fa(616)+Fa(734))
       Fa25(2460) = mup1mup2*F(438)-4*p2mup1*p4mup2*F(4465)-2*Fa(617)
       Fa25(2461) = mup1mup2*F(413)-4*p2mup1*p4mup2*F(4467)-2*Fa(602)
       Fa25(2462) = -(mup1mup2*p1mup3*F(383))+2*mup1mup3*p4mup2*F(4467
     -   )+Fa(643)
       Fa25(2463) = -(mup1mup2*p1mup3*F(418))+2*mup1mup3*p4mup2*F(4465
     -   )+Fa(644)
       Fa25(2464) = 2*mup1mup2*p1mup3*F(263)+Fa(646)+Fa(743)
       Fa25(2465) = Fa(17)+Fa(185)
       Fa25(2466) = Fa(72)+Fa(188)
       Fa25(2467) = mup3mup4*F(411)+2*Fa(125)
       Fa25(2468) = mup3mup4*F(410)+2*Fa(189)
       Fa25(2469) = mup2mup4*F(499)+2*Fa(282)
       Fa25(2470) = mup2mup3*F(502)-2*(Fa(243)+Fa(337))
       Fa25(2471) = mup2mup3*F(503)+2*(Fa(270)+Fa(338))
       Fa25(2472) = mup2mup4*F(504)+2*Fa(283)
       Fa25(2473) = -(mup1mup5*p4mup2*F(1609))+Fa(371)+Fa(725)-Fa(832)
       Fa25(2474) = mup3mup5*p4mup1*F(4464)+Fa(372)-Fa(817)
       Fa25(2475) = mup4mup5*p3mup1*F(396)+mup4mup5*p4mup1*F(4464)+Fa(
     -   504)-Fa(844)
       Return
       End
