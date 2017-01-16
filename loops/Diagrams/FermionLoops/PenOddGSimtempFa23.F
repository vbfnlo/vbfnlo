c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim23(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa23)
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
       Complex*16 Fa23(2179:2277)
       COMMON/FaPenOddGSim/Fa
       Fa23(2179) = Fa(405)+Fa(556)
       Fa23(2180) = Fa(407)+Fa(557)+Fa(696)+Fa(785)+Fa(877)+4*Fa(922)
       Fa23(2181) = Fa(408)+Fa(558)+Fa(697)+Fa(786)+4*Fa(878)+Fa(923)
       Fa23(2182) = (mup1mup4*F(334))/2.+Fa(409)
       Fa23(2183) = (mup1mup4*F(335))/2.+Fa(410)
       Fa23(2184) = mup1mup4*p4mup3*F(2340)+mup3mup4*p4mup1*F(2354)+Fa
     -   (413)+Fa(789)-4*Fa(879)
       Fa23(2185) = Fa(414)+Fa(559)+Fa(790)+4*(Fa(880)-Fa(981))
       Fa23(2186) = mup1mup5*mup3mup4*F(2349)+Fa(415)+Fa(560)+Fa(881)
       Fa23(2187) = Fa(416)+Fa(561)+Fa(698)-8*Fa(882)+4*Fa(924)
       Fa23(2188) = Fa(417)+Fa(699)+Fa(791)+4*(Fa(925)+Fa(982)-Fa(1009
     -   ))
       Fa23(2189) = mup1mup5*mup2mup4*F(2359)+Fa(418)+Fa(700)+Fa(926)
       Fa23(2190) = Fa(419)+Fa(562)+Fa(701)+Fa(792)-4*Fa(883)+Fa(927)
       Fa23(2191) = Fa(421)+Fa(563)+Fa(794)-4*(Fa(884)-Fa(983))
       Fa23(2192) = Fa(422)+Fa(702)+Fa(795)+4*(Fa(928)-Fa(984)+Fa(1010
     -   ))
       Fa23(2193) = (mup1mup4*F(366))/2.+Fa(423)
       Fa23(2194) = (mup1mup3*F(367))/2.+p3mup1*p4mup3*F(879)+Fa(424)+
     -   Fa(564)
       Fa23(2195) = (mup1mup3*F(368))/2.+p3mup1*p4mup3*F(2352)+Fa(425)
     -   +Fa(565)
       Fa23(2196) = (mup1mup3*F(369))/2.+Fa(426)+Fa(566)
       Fa23(2197) = (mup1mup4*F(370))/2.+Fa(427)
       Fa23(2198) = (mup1mup3*F(371))/2.+2*(p3mup1*p4mup3*F(886)-p2mup
     -   1*p4mup3*F(2344))+Fa(428)
       Fa23(2199) = Fa(429)+Fa(574)+4*Fa(803)+Fa(885)
       Fa23(2200) = Fa(430)+Fa(579)+Fa(714)+4*(Fa(805)+Fa(886))+Fa(934
     -   )
       Fa23(2201) = Fa(431)+Fa(580)-8*Fa(806)+Fa(887)
       Fa23(2202) = Fa(432)+Fa(581)+Fa(807)-4*(Fa(888)+Fa(989))
       Fa23(2203) = Fa(433)+Fa(716)+Fa(809)+Fa(937)+4*(Fa(991)-Fa(1015
     -   ))
       Fa23(2204) = mup1mup4*p4mup2*F(2376)+Fa(434)+Fa(584)+Fa(889)+4*
     -   Fa(938)
       Fa23(2205) = 2*mup1mup4*p4mup2*F(1704)+Fa(435)+Fa(585)+Fa(890)+
     -   Fa(939)
       Fa23(2206) = mup3mup4*p4mup1*F(2367)+Fa(436)+Fa(586)+4*Fa(811)+
     -   Fa(892)
       Fa23(2207) = mup3mup4*p4mup1*F(2305)+Fa(437)+Fa(588)-4*(Fa(812)
     -   -Fa(893))
       Fa23(2208) = Fa(438)+Fa(590)+Fa(720)+4*(Fa(894)+Fa(943))
       Fa23(2209) = Fa(439)+Fa(591)+Fa(721)+Fa(895)-4*Fa(944)
       Fa23(2210) = Fa(440)+Fa(592)+8*Fa(813)-4*Fa(897)
       Fa23(2211) = 2*mup1mup4*p4mup2*F(1749)+Fa(441)+Fa(593)+Fa(898)+
     -   Fa(946)
       Fa23(2212) = (mup1mup2*F(366))/2.+Fa(567)+Fa(703)
       Fa23(2213) = (mup1mup2*F(334))/2.+Fa(568)+Fa(704)
       Fa23(2214) = (mup1mup2*F(335))/2.+Fa(569)+Fa(705)
       Fa23(2215) = (mup1mup2*F(370))/2.+Fa(573)+Fa(708)
       Fa23(2216) = 2*mup1mup3*p4mup2*F(1396)+mup1mup2*p1mup3*F(2327)+
     -   Fa(575)+4*Fa(930)
       Fa23(2217) = Fa(711)+Fa(804)-4*(Fa(576)-Fa(931)-Fa(987)+Fa(1013
     -   ))
       Fa23(2218) = mup1mup5*mup2mup3*F(2351)+Fa(577)+Fa(712)+Fa(932)+
     -   4*Fa(988)
       Fa23(2219) = mup1mup2*p1mup3*F(1648)+Fa(578)+Fa(713)+Fa(933)
       Fa23(2220) = mup1mup2*p1mup3*F(284)+mup1mup3*p4mup2*F(2374)+Fa(
     -   582)+Fa(935)
       Fa23(2221) = Fa(715)+Fa(808)+Fa(936)-4*(Fa(583)+Fa(990)-Fa(1014
     -   ))
       Fa23(2222) = 2*mup1mup2*p1mup3*F(1251)+Fa(587)+Fa(718)+Fa(941)
       Fa23(2223) = Fa(589)+Fa(719)+Fa(942)+2*mup1mup2*p1mup3*F(264)*P
     -   (7)
       Fa23(2224) = mup1mup2*p1mup3*F(339)+2*mup1mup3*p4mup2*F(1469)+F
     -   a(594)+4*Fa(949)
       Fa23(2225) = Fa(707)+Fa(798)+4*(mup1mup2*p4mup3*p4mup5*F(277)+F
     -   a(929)-Fa(985)+Fa(1011))
       Fa23(2226) = Fa(717)+Fa(891)+4*Fa(940)+Fa(993)
       Fa23(2227) = Fa(722)+Fa(896)+Fa(945)-4*(Fa(994)+Fa(1017))
       Fa23(2228) = Fa(723)+Fa(899)+Fa(995)-4*(Fa(947)-Fa(1018))
       Fa23(2229) = Fa(724)+Fa(996)+4*(Fa(900)+Fa(948)+Fa(1019))
       Fa23(2230) = Fa(810)-Fa(992)+Fa(1016)
       Fa23(2231) = -Fa(17)+Fa(118)
       Fa23(2232) = Fa(18)+2*Fa(119)
       Fa23(2233) = Fa(19)+2*Fa(120)
       Fa23(2234) = Fa(27)+Fa(127)
       Fa23(2235) = mup3mup4*F(395)-2*Fa(125)
       Fa23(2236) = mup3mup4*F(389)+2*Fa(126)
       Fa23(2237) = mup3mup4*F(384)+2*Fa(128)
       Fa23(2238) = mup3mup4*F(390)+2*Fa(129)
       Fa23(2239) = mup2mup4*F(400)-2*Fa(239)
       Fa23(2240) = mup2mup4*F(401)+2*Fa(241)
       Fa23(2241) = mup2mup3*F(404)+2*(Fa(243)+Fa(337))
       Fa23(2242) = mup2mup3*F(405)+2*(Fa(244)+Fa(338))
       Fa23(2243) = mup2mup3*F(406)+2*Fa(245)+4*Fa(339)
       Fa23(2244) = mup2mup3*F(407)+2*Fa(247)+4*Fa(340)
       Fa23(2245) = mup2mup4*F(408)-2*Fa(249)
       Fa23(2246) = mup2mup4*F(409)-2*Fa(250)
       Fa23(2247) = -(mup1mup5*p4mup2*F(1609))+Fa(331)+Fa(725)-Fa(814)
       Fa23(2248) = -2*mup1mup5*p4mup2*F(356)+Fa(332)-Fa(815)
       Fa23(2249) = 2*mup1mup5*p4mup2*F(4467)-Fa(333)+Fa(816)
       Fa23(2250) = mup3mup5*p4mup1*F(4464)+Fa(334)-Fa(817)
       Fa23(2251) = -Fa(335)+Fa(818)
       Fa23(2252) = -Fa(336)+Fa(819)
       Fa23(2253) = mup4mup5*p3mup1*F(396)+mup4mup5*p4mup1*F(4464)+Fa(
     -   442)-Fa(820)
       Fa23(2254) = Fa(444)-Fa(822)
       Fa23(2255) = mup1mup5*p2mup4+Fa(451)
       Fa23(2256) = -(mup3mup5*p2mup1)-mup1mup5*p2mup3+mup1mup3*p2mup5
       Fa23(2257) = -2*mup1mup5*p4mup2*F(356)-mup1mup2*p2mup5*F(396)+F
     -   a(823)
       Fa23(2258) = -(mup1mup2*p2mup5*F(4464))-2*mup1mup5*p4mup2*F(446
     -   7)+Fa(824)
       Fa23(2259) = Fa(443)-Fa(821)
       Fa23(2260) = -2*p2mup5*p4mup1*F(4464)+Fa(445)
       Fa23(2261) = Fa(446)+2*Fa(595)
       Fa23(2262) = mup2mup5*mup3mup4-mup2mup4*mup3mup5+mup2mup3*mup4m
     -   up5
       Fa23(2263) = mup1mup4*F(395)-2*Fa(449)
       Fa23(2264) = mup1mup4*F(389)+2*Fa(450)
       Fa23(2265) = -2*mup1mup4*p4mup3*F(356)+Fa(452)
       Fa23(2266) = mup1mup5*mup3mup4-mup1mup4*mup3mup5+mup1mup3*mup4m
     -   up5
       Fa23(2267) = 2*mup1mup4*p4mup2*F(356)-Fa(453)+Fa(596)
       Fa23(2268) = mup1mup5*mup2mup4-mup1mup4*mup2mup5+mup1mup2*mup4m
     -   up5
       Fa23(2269) = mup1mup4*F(384)+2*Fa(454)
       Fa23(2270) = mup1mup3*F(410)-2*(Fa(455)+Fa(597))
       Fa23(2271) = 4*p3mup1*p4mup3*F(356)+mup1mup3*F(411)+2*(Fa(456)+
     -   Fa(598))
       Fa23(2272) = -4*p2mup1*p4mup3*F(356)+mup1mup3*F(412)-2*Fa(457)
       Fa23(2273) = mup1mup4*F(390)+2*Fa(458)
       Fa23(2274) = mup1mup3*F(413)-4*p2mup1*p4mup3*F(4467)-2*Fa(459)
       Fa23(2275) = -2*mup1mup4*p4mup3*F(4467)+Fa(460)
       Fa23(2276) = 2*mup1mup4*p4mup2*F(4467)-Fa(461)+Fa(605)
       Fa23(2277) = -(mup3mup4*p4mup1*F(4464))+Fa(462)+Fa(606)
       Return
       End
