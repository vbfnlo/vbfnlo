c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim26(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa26)
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
       Complex*16 Fa26(2476:2574)
       COMMON/FaPenOddGSim/Fa
       Fa26(2476) = 2*p2mup5*p4mup1*F(4464)+Fa(505)
       Fa26(2477) = Fa(506)+2*Fa(595)
       Fa26(2478) = mup1mup4*F(411)+2*Fa(449)
       Fa26(2479) = mup1mup4*F(410)-2*Fa(507)
       Fa26(2480) = mup1mup3*F(384)+2*(Fa(455)+Fa(597))
       Fa26(2481) = -4*p3mup1*p4mup3*F(356)+mup1mup3*F(395)-2*(Fa(456)
     -   +Fa(598))
       Fa26(2482) = mup3mup4*p4mup1*F(4464)+Fa(508)+Fa(650)
       Fa26(2483) = -Fa(509)+Fa(635)+Fa(744)
       Fa26(2484) = mup1mup2*F(410)+2*(Fa(627)+Fa(726))
       Fa26(2485) = mup1mup2*F(411)+2*(Fa(600)+Fa(727))
       Fa26(2486) = 2*mup1mup2*p1mup3*F(376)+Fa(651)+Fa(743)
       Fa26(2487) = Fa(75)+Fa(196)
       Fa26(2488) = Fa(76)+Fa(197)
       Fa26(2489) = Fa(77)+Fa(198)
       Fa26(2490) = Fa(83)+Fa(205)
       Fa26(2491) = Fa(84)-2*Fa(206)
       Fa26(2492) = Fa(86)+Fa(209)
       Fa26(2493) = Fa(87)+Fa(212)
       Fa26(2494) = -4*Fa(89)+Fa(287)
       Fa26(2495) = 4*Fa(190)+Fa(288)
       Fa26(2496) = 4*Fa(191)+Fa(289)
       Fa26(2497) = (mup3mup4*F(593))/2.+Fa(203)
       Fa26(2498) = (mup3mup4*F(594))/2.+Fa(204)
       Fa26(2499) = 4*Fa(207)+Fa(291)
       Fa26(2500) = -4*Fa(108)+Fa(292)
       Fa26(2501) = Fa(208)+Fa(294)
       Fa26(2502) = (mup3mup4*F(621))/2.+Fa(210)
       Fa26(2503) = (mup3mup4*F(622))/2.+Fa(211)
       Fa26(2504) = (mup2mup4*F(613))/2.+Fa(290)
       Fa26(2505) = (mup2mup4*F(614))/2.+Fa(293)
       Fa26(2506) = (mup2mup3*F(617))/2.+Fa(295)+Fa(382)
       Fa26(2507) = (mup2mup3*F(618))/2.+Fa(296)+Fa(383)
       Fa26(2508) = (mup2mup3*F(619))/2.+Fa(297)+Fa(384)
       Fa26(2509) = (mup2mup3*F(620))/2.+Fa(299)+Fa(385)
       Fa26(2510) = (mup2mup4*F(623))/2.+Fa(302)
       Fa26(2511) = (mup2mup4*F(624))/2.+Fa(303)
       Fa26(2512) = mup1mup5*p4mup2*F(1682)+Fa(374)+Fa(745)+Fa(845)
       Fa26(2513) = mup1mup5*p4mup2*F(2372)+Fa(375)+Fa(746)+Fa(846)
       Fa26(2514) = 2*mup1mup5*p4mup2*F(1704)+Fa(376)+Fa(747)+Fa(847)
       Fa26(2515) = Fa(377)+Fa(652)+Fa(848)-4*(p1mup5*p3mup1*p4mup3*F(
     -   263)+Fa(976))
       Fa26(2516) = Fa(378)+Fa(653)+Fa(849)
       Fa26(2517) = Fa(379)+Fa(850)+4*(p1mup5*p3mup1*p4mup3*F(302)+Fa(
     -   977)+Fa(1007))
       Fa26(2518) = Fa(380)+Fa(510)+Fa(851)+Fa(997)
       Fa26(2519) = Fa(381)+Fa(512)+Fa(853)
       Fa26(2520) = -4*Fa(322)+Fa(521)+Fa(856)+Fa(999)
       Fa26(2521) = 4*Fa(392)+Fa(522)+Fa(857)
       Fa26(2522) = Fa(393)+Fa(530)+Fa(862)
       Fa26(2523) = Fa(394)+Fa(670)+Fa(865)
       Fa26(2524) = -4*p1mup2*p2mup1*p4mup5*F(1262)+mup1mup5*p4mup2*F(
     -   1495)+Fa(395)+Fa(758)+Fa(866)
       Fa26(2525) = Fa(396)+Fa(671)+Fa(868)+4*Fa(986)
       Fa26(2526) = Fa(397)+Fa(672)+Fa(869)
       Fa26(2527) = 4*p1mup2*p2mup1*p3mup5*F(1262)+2*mup1mup5*p4mup2*F
     -   (1749)+Fa(398)+Fa(761)+Fa(870)
       Fa26(2528) = mup1mup5*p4mup2*F(2376)+Fa(399)+Fa(871)+4*(Fa(710)
     -   +Fa(1012))
       Fa26(2529) = Fa(511)+Fa(852)+Fa(998)+4*Fa(1008)
       Fa26(2530) = Fa(513)+Fa(654)
       Fa26(2531) = Fa(514)+Fa(655)
       Fa26(2532) = Fa(515)+Fa(656)
       Fa26(2533) = Fa(517)+Fa(657)+Fa(748)+4*(Fa(854)+Fa(901)+Fa(950)
     -   )
       Fa26(2534) = Fa(518)+Fa(658)+Fa(749)+Fa(855)+Fa(902)+4*Fa(951)
       Fa26(2535) = (mup1mup4*F(368))/2.+Fa(519)
       Fa26(2536) = (mup1mup4*F(369))/2.+Fa(520)
       Fa26(2537) = mup3mup4*p4mup1*F(2349)+mup1mup4*p4mup3*F(2352)+Fa
     -   (523)+4*Fa(858)+Fa(903)
       Fa26(2538) = Fa(524)+Fa(659)+Fa(859)+Fa(904)+4*Fa(981)
       Fa26(2539) = mup1mup5*mup3mup4*F(2354)+Fa(525)+Fa(660)+Fa(905)
       Fa26(2540) = Fa(526)+Fa(661)+Fa(750)+8*Fa(882)+Fa(952)
       Fa26(2541) = Fa(527)+Fa(751)+Fa(860)+Fa(953)-4*(Fa(982)-Fa(1009
     -   ))
       Fa26(2542) = mup1mup5*mup2mup4*F(2358)+Fa(528)+Fa(752)+4*Fa(954
     -   )
       Fa26(2543) = Fa(529)+Fa(662)+Fa(753)+4*(-Fa(861)+Fa(906)+Fa(955
     -   ))
       Fa26(2544) = Fa(531)+Fa(663)+Fa(863)+Fa(907)+4*Fa(1000)
       Fa26(2545) = Fa(532)+Fa(754)+Fa(864)-4*(Fa(956)-Fa(984)+Fa(1010
     -   ))
       Fa26(2546) = (mup1mup4*F(367))/2.+Fa(533)
       Fa26(2547) = (mup1mup3*F(366))/2.+p3mup1*p4mup3*F(1215)+Fa(534)
     -   +Fa(664)
       Fa26(2548) = (mup1mup3*F(334))/2.+p3mup1*p4mup3*F(2340)+Fa(535)
     -   +Fa(665)
       Fa26(2549) = (mup1mup3*F(335))/2.+Fa(536)+Fa(666)
       Fa26(2550) = (mup1mup4*F(371))/2.+Fa(537)
       Fa26(2551) = (mup1mup3*F(370))/2.-2*(p3mup1*p4mup3*F(886)-p2mup
     -   1*p4mup3*F(2344))+Fa(538)
       Fa26(2552) = Fa(539)+Fa(674)+Fa(872)-4*Fa(908)
       Fa26(2553) = Fa(540)+Fa(678)+Fa(765)+Fa(909)-4*(Fa(805)-Fa(961)
     -   )
       Fa26(2554) = Fa(541)+Fa(679)+8*Fa(806)+Fa(910)
       Fa26(2555) = Fa(542)+Fa(680)+Fa(874)+4*(Fa(888)+Fa(989))
       Fa26(2556) = Fa(543)+Fa(767)+Fa(876)+4*(Fa(964)-Fa(991)+Fa(1015
     -   ))
       Fa26(2557) = mup1mup4*p4mup2*F(2374)+Fa(544)+Fa(682)+Fa(911)+Fa
     -   (965)
       Fa26(2558) = 2*mup1mup4*p4mup2*F(1396)+Fa(545)+Fa(683)+4*(Fa(91
     -   2)+Fa(966))
       Fa26(2559) = mup3mup4*p4mup1*F(2364)+Fa(546)+Fa(684)-4*(Fa(811)
     -   +Fa(914))
       Fa26(2560) = mup3mup4*p4mup1*F(2187)+Fa(547)+Fa(686)+4*Fa(812)+
     -   Fa(915)
       Fa26(2561) = Fa(548)+Fa(688)+Fa(771)+Fa(916)+Fa(970)
       Fa26(2562) = Fa(549)+Fa(689)+Fa(772)+4*(Fa(917)+Fa(944))
       Fa26(2563) = Fa(550)+Fa(690)-8*Fa(813)+Fa(919)
       Fa26(2564) = 2*mup1mup4*p4mup2*F(1469)+Fa(551)+Fa(691)+4*(Fa(92
     -   0)-Fa(972))
       Fa26(2565) = (mup1mup2*F(367))/2.+Fa(667)+Fa(755)
       Fa26(2566) = (mup1mup2*F(368))/2.+Fa(668)+Fa(756)
       Fa26(2567) = (mup1mup2*F(369))/2.+Fa(669)+Fa(757)
       Fa26(2568) = (mup1mup2*F(371))/2.+Fa(673)+Fa(760)
       Fa26(2569) = 2*mup1mup3*p4mup2*F(1704)+mup1mup2*p1mup3*F(2323)+
     -   Fa(675)+Fa(957)
       Fa26(2570) = Fa(762)+Fa(873)+Fa(958)+Fa(1001)+4*(Fa(576)+Fa(101
     -   3))
       Fa26(2571) = mup1mup5*mup2mup3*F(2356)-4*Fa(676)+Fa(763)+Fa(959
     -   )+Fa(1002)
       Fa26(2572) = mup1mup2*p1mup3*F(1279)+Fa(677)+Fa(764)+4*Fa(960)
       Fa26(2573) = mup1mup2*p1mup3*F(601)+mup1mup3*p4mup2*F(2376)+Fa(
     -   681)+4*Fa(962)
       Fa26(2574) = Fa(766)+Fa(875)+Fa(963)+Fa(1003)+4*(Fa(583)-Fa(101
     -   4))
       Return
       End
