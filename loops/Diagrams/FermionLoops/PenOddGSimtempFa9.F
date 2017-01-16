c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim9(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,
     -   p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mu
     -   p3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup
     -   3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa9)
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
       Complex*16 Fa9(793:891)
       COMMON/FaPenOddGSim/Fa
       Fa9(793) = mup1mup5*(p1mup4*F(341)+p3mup4*F(984)+p4mup4*F(988)+
     -   p2mup4*F(1743))
       Fa9(794) = -4*mup3mup4*p3mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(379))+mup1mup5*mup3mup4*F(2364)
       Fa9(795) = 4*mup2mup4*p3mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mup
     -   5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup4*F(2320)
       Fa9(796) = mup1mup5*(p4mup3*F(879)+p2mup3*F(1771)+2*(p1mup3*F(2
     -   64)+p3mup3*F(282))*P(7))
       Fa9(797) = mup1mup5*p1mup2*F(1773)+2*mup1mup5*p3mup2*F(1777)+mu
     -   p2mup5*p2mup1*F(2362)
       Fa9(798) = 4*mup2mup3*p4mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mup
     -   5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup3*F(2362)
       Fa9(799) = mup1mup5*p1mup3*F(339)+mup1mup5*p3mup3*F(340)-2*mup1
     -   mup5*p4mup3*F(986)+p2mup3*(mup1mup5*F(336)+4*p3mup1*(p3mup5*
     -   F(273)+p4mup5*F(357)+p2mup5*F(1500)))+mup3mup5*p3mup1*F(2305
     -   )+4*(p1mup5*p2mup3*p3mup1*F(263)+p2mup5*p3mup1*p4mup3*F(2345
     -   ))
       Fa9(800) = mup1mup5*(p1mup3*F(284)-2*(p3mup3*F(290)+p4mup3*F(88
     -   6))+p2mup3*F(2332))
       Fa9(801) = 2*mup1mup5*p1mup2*F(1465)+mup1mup5*p3mup2*F(1466)+mu
     -   p2mup5*p2mup1*F(2319)
       Fa9(802) = mup1mup5*p1mup2*F(2355)+mup1mup5*p3mup2*F(2357)+mup2
     -   mup5*p2mup1*F(2370)
       Fa9(803) = p3mup1*(p1mup3*p3mup4*F(294)+p1mup4*p3mup3*F(299))
       Fa9(804) = -4*mup2mup3*p3mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup3*F(2330)
       Fa9(805) = mup2mup5*(p2mup3*(p1mup4*F(346)-p4mup4*F(362)+p3mup4
     -   *F(2188))+p2mup4*(-(p1mup3*F(258))+p3mup3*F(282)+p4mup3*F(23
     -   45)+p2mup3*F(2395)))
       Fa9(806) = p2mup1*(p2mup3*p4mup4*F(2386)-p2mup4*p4mup3*F(2388))
       Fa9(807) = 4*mup3mup4*p2mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mup
     -   5*F(357)+p2mup5*F(379))+mup1mup5*mup3mup4*F(2187)
       Fa9(808) = 4*mup2mup3*p2mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mup
     -   5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup3*F(2334)
       Fa9(809) = -4*mup2mup4*p2mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup4*F(2365)
       Fa9(810) = (mup1mup4*mup2mup3-mup1mup3*mup2mup4+mup1mup2*mup3mu
     -   p4)*(p1mup5*F(263)+p3mup5*F(273)+p4mup5*F(357)+p2mup5*F(379)
     -   )
       Fa9(811) = p1mup3*p3mup4*(2*p3mup1*F(256)+p4mup1*F(345))+p1mup4
     -   *p3mup3*(-2*p3mup1*F(279)+p4mup1*F(347))
       Fa9(812) = p2mup1*(p2mup3*p4mup4*F(2380)-p2mup4*p4mup3*F(2381))
       Fa9(813) = p2mup1*(p2mup3*p4mup4*F(360)+p2mup4*p4mup3*F(2394))
       Fa9(814) = mup1mup5*(2*p1mup2*F(1262)+p3mup2*F(1607))
       Fa9(815) = mup1mup5*(p1mup2*F(380)+2*p3mup2*F(382))
       Fa9(816) = mup1mup5*(p1mup2*F(383)+2*p3mup2*F(4465))
       Fa9(817) = -(mup3mup5*p3mup1*F(396))+mup1mup5*(2*(p1mup3*F(376)
     -   +p2mup3*F(377))+p3mup3*F(378)+p4mup3*F(1012))
       Fa9(818) = mup1mup5*(p1mup3*F(380)+2*(p4mup3*F(356)+p2mup3*F(38
     -   1)+p3mup3*F(382)))
       Fa9(819) = mup1mup5*(p1mup3*F(383)+2*(p3mup3*F(4465)+p4mup3*F(4
     -   467)+p2mup3*F(4469)))
       Fa9(820) = mup1mup5*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(3
     -   83)+p3mup4*F(386))
       Fa9(821) = mup1mup5*(p1mup4*F(380)+2*(p4mup4*F(356)+p2mup4*F(38
     -   1)+p3mup4*F(382)))
       Fa9(822) = mup1mup5*(p1mup4*F(383)+2*(p3mup4*F(4465)+p4mup4*F(4
     -   467)+p2mup4*F(4469)))
       Fa9(823) = mup2mup5*p2mup1*F(396)-2*mup1mup5*(p3mup2*F(382)+p1m
     -   up2*F(1610))
       Fa9(824) = mup2mup5*p2mup1*F(4464)-2*mup1mup5*(p1mup2*F(357)+p3
     -   mup2*F(4465))
       Fa9(825) = mup1mup5*(2*p1mup2*F(1262)+p3mup2*F(1622))
       Fa9(826) = mup1mup5*(p1mup2*F(418)+2*p3mup2*F(4466))
       Fa9(827) = mup1mup5*(2*(p1mup3*F(263)+p2mup3*F(379))+p4mup3*F(3
     -   83)+p3mup3*F(418))-mup3mup5*p3mup1*F(4463)
       Fa9(828) = mup1mup5*(p1mup3*F(418)+2*(p4mup3*F(4465)+p3mup3*F(4
     -   466)+p2mup3*F(4468)))
       Fa9(829) = mup1mup5*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(3
     -   83)+p3mup4*F(418))
       Fa9(830) = mup1mup5*(p1mup4*F(418)+2*(p4mup4*F(4465)+p3mup4*F(4
     -   466)+p2mup4*F(4468)))
       Fa9(831) = mup2mup5*p2mup1*F(4463)-2*mup1mup5*(p1mup2*F(273)+p3
     -   mup2*F(4466))
       Fa9(832) = mup1mup5*(p3mup2*F(1607)-2*p1mup2*F(1608))
       Fa9(833) = mup1mup5*(p1mup2*F(443)+2*p3mup2*F(445))
       Fa9(834) = -(mup3mup5*p3mup1*F(396))+mup1mup5*(2*(p1mup3*F(376)
     -   +p2mup3*F(377))+p3mup3*F(378)+p4mup3*F(448))
       Fa9(835) = mup1mup5*(p1mup3*F(443)+2*(p2mup3*F(444)+p3mup3*F(44
     -   5)+p4mup3*F(459)))
       Fa9(836) = mup1mup5*(2*(p1mup4*F(376)+p2mup4*F(377))+p3mup4*F(3
     -   78)+p4mup4*F(448))
       Fa9(837) = mup1mup5*(p1mup4*F(443)+2*(p2mup4*F(444)+p3mup4*F(44
     -   5)+p4mup4*F(459)))
       Fa9(838) = -(mup2mup5*p2mup1*F(396))+2*mup1mup5*(p3mup2*F(382)+
     -   p1mup2*F(1610))
       Fa9(839) = -(mup2mup5*p2mup1*F(477))+2*mup1mup5*(p3mup2*F(445)+
     -   p1mup2*F(1629))
       Fa9(840) = mup2mup5*p2mup1*F(477)-2*mup1mup5*(p3mup2*F(445)+p1m
     -   up2*F(1629))
       Fa9(841) = mup1mup5*(2*p1mup2*F(1608)-p3mup2*F(1622))
       Fa9(842) = -(mup2mup5*p2mup1*F(4463))+2*mup1mup5*(p1mup2*F(273)
     -   +p3mup2*F(4466))
       Fa9(843) = -(mup2mup5*p2mup1*F(4464))+2*mup1mup5*(p1mup2*F(357)
     -   +p3mup2*F(4465))
       Fa9(844) = mup1mup5*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(3
     -   83)-p3mup4*F(385))
       Fa9(845) = mup1mup5*(p3mup2*F(1252)+2*p1mup2*F(1258))
       Fa9(846) = mup1mup5*(p1mup2*F(1278)+2*p3mup2*F(1286))
       Fa9(847) = mup1mup5*(p1mup2*F(1697)+p3mup2*F(2347))
       Fa9(848) = mup1mup5*(-2*(p2mup3*F(1250)+p1mup3*F(1251))+p3mup3*
     -   F(1253)+p4mup3*F(1256))-4*p3mup1*(p3mup5*p4mup3*F(273)+p4mup
     -   3*p4mup5*F(277)+p2mup5*p4mup3*F(379)+p2mup3*p4mup5*F(2338))+
     -   mup3mup5*p3mup1*F(2341)
       Fa9(849) = mup1mup5*(p1mup3*F(1279)+2*(p2mup3*F(1284)+p3mup3*F(
     -   1287))+p4mup3*F(2340))
       Fa9(850) = 4*p2mup3*p3mup1*(2*(p4mup5*F(304)+p3mup5*F(311)+p2mu
     -   p5*F(1303))+p1mup5*F(1657))+8*p3mup1*p4mup3*(p3mup5*F(306)+p
     -   4mup5*F(316)+p2mup5*F(1831))+mup1mup5*p2mup3*F(2324)+mup3mup
     -   5*p3mup1*F(2336)+mup1mup5*(p1mup3*F(2323)+p3mup3*F(2326)-2*p
     -   4mup3*F(2344))
       Fa9(851) = mup1mup5*p3mup4*F(1357)+2*(mup1mup5*p2mup4*F(1359)+m
     -   up1mup5*p1mup4*F(1360))+mup1mup5*p4mup4*F(1362)+4*(p1mup5*p3
     -   mup4*p4mup1*F(263)+p3mup4*p4mup1*(p3mup5*F(356)+p4mup5*F(357
     -   )+p2mup5*F(379))-p2mup4*p3mup5*p4mup1*F(2318))
       Fa9(852) = 4*p1mup5*p4mup1*(p3mup4*F(295)+p2mup4*F(1069))+8*(p2
     -   mup4*p4mup1*(p3mup5*F(304)+p4mup5*F(314)+p2mup5*F(692))+p3mu
     -   p4*p4mup1*(p3mup5*F(306)+p4mup5*F(316)+p2mup5*F(1831)))+mup1
     -   mup5*p1mup4*F(2312)+mup1mup5*p2mup4*F(2313)+mup1mup5*p3mup4*
     -   F(2314)+mup1mup5*p4mup4*F(2315)
       Fa9(853) = mup1mup5*(p1mup4*F(1389)-2*(p2mup4*F(1394)-p4mup4*F(
     -   1397))+p3mup4*F(2322))
       Fa9(854) = mup2mup5*(p3mup4*p4mup3*F(274)+p4mup3*p4mup4*F(598)+
     -   p3mup3*p4mup4*F(599)-p1mup3*(p4mup4*F(259)+2*(p2mup4*F(254)+
     -   p3mup4*F(256)+p1mup4*F(882)))+p1mup4*(p4mup3*F(261)+2*(p3mup
     -   3*F(279)-p2mup3*F(2379)))-p2mup3*p4mup4*F(2380)+p2mup4*p4mup
     -   3*F(2381)+2*(p2mup4*p3mup3*F(280)+p3mup3*p3mup4*F(281)-p2mup
     -   3*p3mup4*F(2378)-p2mup3*p2mup4*F(2382)))
       Fa9(855) = mup2mup5*(p1mup3*(p4mup4*F(300)+4*(p3mup4*F(294)+p1m
     -   up4*F(889)+p2mup4*F(2383)))+4*(p1mup4*(p3mup3*F(299)+p4mup3*
     -   F(301))-p1mup4*p2mup3*F(2384))-8*(-(p3mup4*p4mup3*F(309))-p4
     -   mup3*p4mup4*F(318)+p3mup3*(p3mup4*F(307)-p4mup4*F(317)+p2mup
     -   4*F(2385))-p2mup3*p4mup4*F(2386)-p2mup3*p3mup4*F(2387)+p2mup
     -   4*p4mup3*F(2388)-p2mup3*p2mup4*F(2389)))
       Fa9(856) = mup1mup5*(p2mup4*F(596)+2*(p1mup4*F(253)+p3mup4*F(12
     -   14))+p4mup4*F(1216))-4*(p2mup4*p4mup1*(p1mup5*F(263)+p3mup5*
     -   F(273)+p4mup5*F(357)+p2mup5*F(1470))-p2mup5*p3mup4*p4mup1*F(
     -   2188))
       Fa9(857) = mup1mup5*(p1mup4*F(602)-2*(p4mup4*F(887)-p3mup4*F(12
     -   18))+p2mup4*F(2321))
       Fa9(858) = p4mup1*(p1mup3*p3mup4*F(294)+p1mup4*p3mup3*F(299))
       Fa9(859) = -4*mup3mup4*p4mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(379))+mup1mup5*mup3mup4*F(2360)
       Fa9(860) = 4*mup2mup4*p4mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mup
     -   5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup4*F(2337)
       Fa9(861) = mup2mup5*(p1mup3*p3mup4*F(345)+p3mup4*p4mup3*F(358)+
     -   p1mup4*(p3mup3*F(347)+2*(p2mup3*F(343)-p4mup3*F(359)+p1mup3*
     -   F(938)))+p2mup3*p3mup4*F(2390)+p3mup3*(p3mup4*F(352)-p4mup4*
     -   F(632)+p2mup4*F(2391))-2*(p4mup4*(p1mup3*F(359)+p2mup3*F(360
     -   ))-p4mup3*p4mup4*F(365)-p1mup3*p2mup4*F(2392)-p2mup3*p2mup4*
     -   F(2393)+p2mup4*p4mup3*F(2394)))
       Fa9(862) = mup1mup5*(2*p1mup4*F(630)+p3mup4*F(1246)+p4mup4*F(12
     -   49)+p2mup4*F(1463))
       Fa9(863) = 4*mup3mup4*p3mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mup
     -   5*F(357)+p2mup5*F(379))+mup1mup5*mup3mup4*F(2367)
       Fa9(864) = -4*mup2mup4*p3mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup4*F(2319)
       Fa9(865) = mup1mup5*(p4mup3*F(1215)+p2mup3*F(1492)+2*(p1mup3*F(
     -   258)+p3mup3*F(600))*P(7))
       Fa9(866) = 2*mup1mup5*(p1mup2*F(1494)+p3mup2*F(1498))+mup2mup5*
     -   p2mup1*F(2341)
       Fa9(867) = -4*mup2mup3*p4mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup3*F(2341)
       Fa9(868) = mup1mup5*p3mup3*F(629)+2*(mup1mup5*p1mup3*F(628)+mup
     -   1mup5*p4mup3*F(986))+p2mup3*(mup1mup5*F(625)-4*p3mup1*(p3mup
     -   5*F(273)+p4mup5*F(357)+p2mup5*F(1500)))+mup3mup5*p3mup1*F(21
     -   87)-4*(p1mup5*p2mup3*p3mup1*F(263)+p2mup5*p3mup1*p4mup3*F(23
     -   45))
       Fa9(869) = mup1mup5*(p1mup3*F(601)+2*(p3mup3*F(290)+p4mup3*F(88
     -   6))+p2mup3*F(2325))
       Fa9(870) = 2*mup1mup5*p1mup2*F(1745)+mup1mup5*p3mup2*F(1746)+mu
     -   p2mup5*p2mup1*F(2320)
       Fa9(871) = mup1mup5*p1mup2*F(2346)+mup1mup5*p3mup2*F(2348)+mup2
     -   mup5*p2mup1*F(2371)
       Fa9(872) = -4*p1mup3*p3mup1*p3mup4*F(294)+p1mup4*p3mup1*p3mup3*
     -   F(608)
       Fa9(873) = 4*mup2mup3*p3mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mup
     -   5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup3*F(2317)
       Fa9(874) = -4*mup3mup4*p2mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(379))+mup1mup5*mup3mup4*F(2305)
       Fa9(875) = -4*mup2mup3*p2mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mu
     -   p5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup3*F(2328)
       Fa9(876) = 4*mup2mup4*p2mup1*(p1mup5*F(263)+p3mup5*F(273)+p4mup
     -   5*F(357)+p2mup5*F(379))+mup1mup5*mup2mup4*F(2353)
       Fa9(877) = mup3mup5*(p4mup2*p4mup4*F(1685)+4*(p2mup4*p4mup2*F(6
     -   40)-p1mup2*(2*(p1mup4*F(633)+p2mup4*F(635)+p3mup4*F(636))+p4
     -   mup4*F(1367))+p1mup4*p4mup2*F(1368)+p3mup4*p4mup2*F(2396))-8
     -   *(p3mup2*p3mup4*F(2397)+p1mup4*p3mup2*F(2398)+p2mup4*p3mup2*
     -   F(2399))+p3mup2*p4mup4*F(2400))
       Fa9(878) = mup3mup5*(p1mup4*p4mup2*F(1401)+p1mup2*(p1mup4*F(677
     -   )+p2mup4*F(679)-p4mup4*F(1400)+p3mup4*F(2401))+p1mup4*p3mup2
     -   *F(2403)+2*(p4mup2*p4mup4*F(1402)-p3mup2*p3mup4*F(2402)-p2mu
     -   p4*(p4mup2*F(689)+p3mup2*F(2404))+p3mup4*p4mup2*F(2405)-p3mu
     -   p2*p4mup4*F(2406)))
       Fa9(879) = p4mup1*(-(p1mup3*p4mup4*F(609))+p1mup3*p2mup4*F(2383
     -   )+p1mup4*(p4mup3*F(301)+p1mup3*F(889)-p2mup3*F(2384))+2*(p3m
     -   up4*p4mup3*F(309)+p4mup3*p4mup4*F(318)-p3mup3*(p3mup4*F(307)
     -   -p4mup4*F(317)+p2mup4*F(2385))+p2mup3*p4mup4*F(2386)+p2mup3*
     -   p3mup4*F(2387)-p2mup4*p4mup3*F(2388)+p2mup3*p2mup4*F(2389)))
       Fa9(880) = mup1mup5*(p3mup4*p4mup3*F(274)+p4mup3*p4mup4*F(598)+
     -   p3mup3*p4mup4*F(599)-p1mup3*(p4mup4*F(259)+2*(p2mup4*F(254)+
     -   p3mup4*F(256)+p1mup4*F(882)))+p1mup4*(p4mup3*F(261)+2*(p3mup
     -   3*F(279)-p2mup3*F(2379)))-p2mup3*p4mup4*F(2380)+p2mup4*p4mup
     -   3*F(2381)+2*(p2mup4*p3mup3*F(280)+p3mup3*p3mup4*F(281)-p2mup
     -   3*p3mup4*F(2378)-p2mup3*p2mup4*F(2382)))
       Fa9(881) = mup1mup5*(p1mup3*(p4mup4*F(300)+4*(p3mup4*F(294)+p1m
     -   up4*F(889)+p2mup4*F(2383)))+4*(p1mup4*(p3mup3*F(299)+p4mup3*
     -   F(301))-p1mup4*p2mup3*F(2384))-8*(-(p3mup4*p4mup3*F(309))-p4
     -   mup3*p4mup4*F(318)+p3mup3*(p3mup4*F(307)-p4mup4*F(317)+p2mup
     -   4*F(2385))-p2mup3*p4mup4*F(2386)-p2mup3*p3mup4*F(2387)+p2mup
     -   4*p4mup3*F(2388)-p2mup3*p2mup4*F(2389)))
       Fa9(882) = p4mup1*(p2mup4*p4mup2*F(689)-p4mup2*p4mup4*F(1402)+p
     -   3mup2*p3mup4*F(2402))
       Fa9(883) = mup3mup5*(p1mup4*p3mup2*F(346)-p3mup2*p4mup4*F(362)+
     -   p1mup2*p3mup4*F(1262)+p2mup4*p3mup2*F(2318)+p3mup4*p4mup2*F(
     -   2369)-p3mup2*p3mup4*F(2407))
       Fa9(884) = mup1mup5*(p1mup3*p3mup4*F(345)+p3mup4*p4mup3*F(358)+
     -   p1mup4*(p3mup3*F(347)+2*(p2mup3*F(343)-p4mup3*F(359)+p1mup3*
     -   F(938)))+p2mup3*p3mup4*F(2390)+p3mup3*(p3mup4*F(352)-p4mup4*
     -   F(632)+p2mup4*F(2391))-2*(p4mup4*(p1mup3*F(359)+p2mup3*F(360
     -   ))-p4mup3*p4mup4*F(365)-p1mup3*p2mup4*F(2392)-p2mup3*p2mup4*
     -   F(2393)+p2mup4*p4mup3*F(2394)))
       Fa9(885) = p3mup1*(p1mup3*(p4mup4*F(300)+4*(p1mup4*F(889)+p2mup
     -   4*F(2383)))+4*p1mup4*(p4mup3*F(301)-p2mup3*F(2384))+8*(-(p3m
     -   up3*p3mup4*F(307))+p3mup4*p4mup3*F(309)+p3mup3*p4mup4*F(317)
     -   +p4mup3*p4mup4*F(318)-p2mup4*p3mup3*F(2385)+p2mup3*p4mup4*F(
     -   2386)+p2mup3*p3mup4*F(2387)-p2mup4*p4mup3*F(2388)+p2mup3*p2m
     -   up4*F(2389)))
       Fa9(886) = mup3mup5*(p2mup4*p4mup2*F(946)+p1mup2*(p2mup4*F(939)
     -   -2*(p4mup4*F(359)-p3mup4*F(2408)))+p2mup4*p3mup2*F(2410)+2*(
     -   p4mup2*p4mup4*F(365)-p1mup4*(p4mup2*F(359)-p1mup2*F(938)+p3m
     -   up2*F(944))-p3mup2*p4mup4*F(949)-p3mup2*p3mup4*F(2409)+p3mup
     -   4*p4mup2*F(2411)))
       Fa9(887) = -(p2mup1*(4*(p1mup4*p4mup3*F(301)+p1mup3*(p3mup4*F(2
     -   94)-p4mup4*F(609)+p1mup4*F(889)+p2mup4*F(2383))-p1mup4*p2mup
     -   3*F(2384))-p3mup3*(p1mup4*F(608)+8*(p3mup4*F(307)-p4mup4*F(3
     -   17)+p2mup4*F(2385)))+8*(p3mup4*p4mup3*F(309)+p4mup3*p4mup4*F
     -   (318)+p2mup3*p3mup4*F(2387)+p2mup3*p2mup4*F(2389))))
       Fa9(888) = mup1mup5*(p2mup3*(p1mup4*F(346)-p4mup4*F(362)+p3mup4
     -   *F(2188))+p2mup4*(-(p1mup3*F(258))+p3mup3*F(282)+p4mup3*F(23
     -   45)+p2mup3*F(2395)))
       Fa9(889) = 4*p2mup1*(p1mup2*p3mup4*F(2401)-2*(p2mup4*p4mup2*F(6
     -   89)-p4mup2*p4mup4*F(1402)+p3mup2*p3mup4*F(2402)))
       Fa9(890) = p3mup1*(-4*p1mup2*p3mup4*F(2401)-8*p4mup2*(p4mup4*F(
     -   1402)+p3mup4*F(2405))+p3mup2*(8*(p2mup4*F(2404)+p4mup4*F(240
     -   6))+p1mup4*F(2412)))
       Fa9(891) = -(mup1mup2*(4*(p1mup4*p4mup3*F(301)+p1mup3*(p3mup4*F
     -   (294)-p4mup4*F(609)+p1mup4*F(889)+p2mup4*F(2383))-p1mup4*p2m
     -   up3*F(2384))-p3mup3*(p1mup4*F(608)+8*(p3mup4*F(307)-p4mup4*F
     -   (317)+p2mup4*F(2385)))+8*(p3mup4*p4mup3*F(309)+p4mup3*p4mup4
     -   *F(318)+p2mup3*p4mup4*F(2386)+p2mup3*p3mup4*F(2387)-p2mup4*p
     -   4mup3*F(2388)+p2mup3*p2mup4*F(2389))))
       Return
       End
