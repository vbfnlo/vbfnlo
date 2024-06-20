c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim24(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa24)
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
       Complex*16 Fa24(2278:2376)
       COMMON/FaPenOddGSim/Fa
       Fa24(2278) = mup3mup4*p2mup1+Fa(463)
       Fa24(2279) = 2*mup1mup4*p4mup2*F(356)+mup1mup2*p2mup4*F(396)+Fa
     -   (609)
       Fa24(2280) = Fa(464)+Fa(610)+Fa(730)
       Fa24(2281) = mup1mup2*p2mup4*F(4464)+2*mup1mup4*p4mup2*F(4467)+
     -   Fa(611)
       Fa24(2282) = mup1mup2*F(384)+2*(Fa(599)-Fa(726))
       Fa24(2283) = mup1mup2*F(395)-2*(Fa(600)+Fa(727))
       Fa24(2284) = 4*p2mup1*p4mup2*F(356)+mup1mup2*F(389)+2*Fa(601)
       Fa24(2285) = mup1mup2*F(390)+4*p2mup1*p4mup2*F(4467)+2*Fa(602)
       Fa24(2286) = mup1mup2*p1mup3*F(383)-2*mup1mup3*p4mup2*F(4467)+F
     -   a(603)
       Fa24(2287) = mup1mup5*mup2mup3+Fa(728)
       Fa24(2288) = -2*mup1mup3*p4mup2*F(356)+mup1mup2*p1mup3*F(380)+F
     -   a(604)
       Fa24(2289) = 2*mup1mup2*p1mup3*F(376)+Fa(607)-Fa(729)
       Fa24(2290) = 2*mup1mup3*p4mup2*F(356)+Fa(608)
       Fa24(2291) = 2*mup1mup3*p4mup2*F(4467)+Fa(612)
       Fa24(2292) = Fa(30)+Fa(133)
       Fa24(2293) = Fa(31)+2*Fa(134)
       Fa24(2294) = Fa(37)+Fa(140)
       Fa24(2295) = mup3mup4*F(420)+2*Fa(138)
       Fa24(2296) = mup3mup4*F(421)+2*Fa(139)
       Fa24(2297) = mup3mup4*F(419)+2*Fa(141)
       Fa24(2298) = mup2mup4*F(429)-2*Fa(255)
       Fa24(2299) = mup2mup4*F(430)-2*Fa(257)
       Fa24(2300) = mup2mup3*F(431)+2*(Fa(258)+Fa(350))
       Fa24(2301) = mup2mup3*F(432)+2*(Fa(259)+Fa(351))
       Fa24(2302) = mup2mup3*F(433)+2*Fa(260)+4*Fa(352)
       Fa24(2303) = mup2mup4*F(434)-2*Fa(262)
       Fa24(2304) = -(mup1mup5*p4mup2*F(1609))+Fa(346)+Fa(732)-Fa(825)
       Fa24(2305) = 2*mup1mup5*p4mup2*F(4465)-Fa(347)+Fa(826)
       Fa24(2306) = mup3mup5*p4mup1*F(4464)+Fa(348)-Fa(827)
       Fa24(2307) = -Fa(349)+Fa(828)
       Fa24(2308) = mup4mup5*p3mup1*F(4463)+mup4mup5*p4mup1*F(4464)+Fa
     -   (465)-Fa(829)
       Fa24(2309) = -(mup1mup2*p2mup5*F(4463))-2*mup1mup5*p4mup2*F(446
     -   5)+Fa(831)
       Fa24(2310) = Fa(466)-Fa(830)
       Fa24(2311) = 2*p2mup5*p4mup1*F(4464)+Fa(467)
       Fa24(2312) = Fa(468)+2*Fa(595)
       Fa24(2313) = mup1mup4*F(435)-2*Fa(470)
       Fa24(2314) = mup1mup4*F(421)+2*Fa(471)
       Fa24(2315) = -2*mup1mup4*p4mup3*F(4465)+Fa(472)
       Fa24(2316) = 2*mup1mup4*p4mup2*F(4465)-Fa(473)+Fa(613)
       Fa24(2317) = mup1mup4*F(419)+2*Fa(474)
       Fa24(2318) = mup1mup3*F(436)-2*(Fa(475)+Fa(614))
       Fa24(2319) = mup1mup3*F(437)+4*p3mup1*p4mup3*F(4465)+2*(Fa(476)
     -   +Fa(598))
       Fa24(2320) = mup1mup3*F(438)-4*p2mup1*p4mup3*F(4465)-2*Fa(477)
       Fa24(2321) = mup3mup4*p4mup1*F(4464)+Fa(478)-Fa(619)
       Fa24(2322) = mup1mup2*p2mup4*F(4463)+2*mup1mup4*p4mup2*F(4465)+
     -   Fa(622)
       Fa24(2323) = Fa(479)+Fa(623)+Fa(730)
       Fa24(2324) = mup1mup2*F(419)+2*(Fa(615)-Fa(733))
       Fa24(2325) = mup1mup2*F(435)-2*(Fa(616)+Fa(734))
       Fa24(2326) = mup1mup2*F(421)+4*p2mup1*p4mup2*F(4465)+2*Fa(617)
       Fa24(2327) = mup1mup2*p1mup3*F(418)-2*mup1mup3*p4mup2*F(4465)+F
     -   a(618)
       Fa24(2328) = 2*mup1mup2*p1mup3*F(263)+Fa(620)-Fa(729)
       Fa24(2329) = 2*mup1mup3*p4mup2*F(4465)+Fa(621)
       Fa24(2330) = Fa(40)+Fa(146)
       Fa24(2331) = Fa(41)-2*Fa(147)
       Fa24(2332) = Fa(48)+Fa(154)
       Fa24(2333) = mup3mup4*F(465)+2*Fa(152)
       Fa24(2334) = mup3mup4*F(412)+2*Fa(153)
       Fa24(2335) = mup3mup4*F(446)+2*Fa(155)
       Fa24(2336) = mup3mup4*F(460)+2*Fa(156)
       Fa24(2337) = mup2mup4*F(406)+2*Fa(239)
       Fa24(2338) = mup2mup4*F(469)+2*Fa(267)
       Fa24(2339) = mup2mup3*F(472)-2*(Fa(269)+Fa(360))
       Fa24(2340) = mup2mup3*F(473)+2*(Fa(270)-Fa(361))
       Fa24(2341) = mup2mup3*F(400)-2*Fa(245)-4*Fa(339)
       Fa24(2342) = mup2mup3*F(474)-2*Fa(271)-4*Fa(362)
       Fa24(2343) = mup2mup4*F(475)+2*Fa(273)
       Fa24(2344) = mup2mup4*F(476)+2*Fa(274)
       Fa24(2345) = -(mup1mup5*p4mup2*F(1623))+Fa(356)+Fa(735)-Fa(832)
       Fa24(2346) = 2*mup1mup5*p4mup2*F(356)-Fa(332)+Fa(815)
       Fa24(2347) = 2*mup1mup5*p4mup2*F(459)-Fa(357)+Fa(833)
       Fa24(2348) = -(mup3mup5*p4mup1*F(477))+Fa(358)+Fa(834)
       Fa24(2349) = Fa(335)-Fa(818)
       Fa24(2350) = Fa(359)-Fa(835)
       Fa24(2351) = mup4mup5*p3mup1*F(396)+mup4mup5*p4mup1*F(477)+Fa(4
     -   80)-Fa(836)
       Fa24(2352) = -Fa(481)+Fa(837)
       Fa24(2353) = -(mup1mup5*p2mup4)+Fa(487)
       Fa24(2354) = mup3mup5*p2mup1+mup1mup5*p2mup3-mup1mup3*p2mup5
       Fa24(2355) = 2*mup1mup5*p4mup2*F(356)+mup1mup2*p2mup5*F(396)+Fa
     -   (838)
       Fa24(2356) = 2*mup1mup5*p4mup2*F(459)+mup1mup2*p2mup5*F(477)+Fa
     -   (839)
       Fa24(2357) = -Fa(443)+Fa(821)
       Fa24(2358) = 2*p2mup5*p4mup1*F(477)+Fa(482)
       Fa24(2359) = Fa(483)+2*Fa(624)
       Fa24(2360) = mup1mup4*F(450)+2*Fa(486)
       Fa24(2361) = mup1mup4*F(412)-2*Fa(450)
       Fa24(2362) = -2*mup1mup4*p4mup2*F(356)+Fa(453)-Fa(596)
       Fa24(2363) = mup1mup4*F(446)-2*Fa(488)
       Fa24(2364) = mup1mup3*F(478)+2*(Fa(455)+Fa(625))
       Fa24(2365) = -4*p3mup1*p4mup3*F(356)+mup1mup3*F(479)-2*(Fa(456)
     -   +Fa(626))
       Fa24(2366) = 4*p2mup1*p4mup3*F(356)+mup1mup3*F(389)+2*Fa(457)
       Fa24(2367) = mup1mup4*F(460)-2*Fa(489)
       Fa24(2368) = 4*p2mup1*p4mup3*F(459)+mup1mup3*F(480)+2*Fa(490)
       Fa24(2369) = -2*mup1mup4*p4mup3*F(459)+Fa(491)
       Fa24(2370) = -2*mup1mup4*p4mup2*F(459)+Fa(492)-Fa(631)
       Fa24(2371) = mup3mup4*p4mup1*F(477)+Fa(493)+Fa(632)
       Fa24(2372) = -(mup3mup4*p2mup1)+Fa(494)
       Fa24(2373) = -2*mup1mup4*p4mup2*F(356)-mup1mup2*p2mup4*F(396)+F
     -   a(634)
       Fa24(2374) = -Fa(495)+Fa(635)+Fa(740)
       Fa24(2375) = -2*mup1mup4*p4mup2*F(459)-mup1mup2*p2mup4*F(477)+F
     -   a(636)
       Fa24(2376) = mup1mup2*F(446)+2*(Fa(627)+Fa(736))
       Return
       End
