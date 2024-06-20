c************** Calling the Fa functions*************************
       subroutine HexBoxFa16(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa16)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxFhlFunctions/F
       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(259),F(1314)
      Real*8 P(1362) 
        COMMON/HexBoxPFunctions/P
       Complex*16 Fa16(181:192)
       COMMON/HexBoxFaFunctions/Fa
       Fa16(181) = p1mup4*F(842)+p6mup4*F(843)+p4mup4*F(844)+p5mup4*F(
     &   845)+p2mup4*F(846)
       Fa16(182) = p1mup3*F(847)+2*(p6mup3*F(848)+p4mup3*F(849)+p5mup3
     &   *F(850)+p2mup3*F(851))
       Fa16(183) = mup3mup4*F(852)-4*(p1mup3*(p5mup4*F(500)+p4mup4*F(8
     &   53)+p6mup4*F(856)+p2mup4*F(857)-p1mup4*F(858))+p2mup3*(p5mup
     &   4*F(501)+p4mup4*F(854)-p6mup4*F(855)+p1mup4*F(857)+p2mup4*F(
     &   865))+p4mup3*(p5mup4*F(511)+p4mup4*F(859)+p6mup4*F(860)-p1mu
     &   p4*F(866)-p2mup4*F(868))+p5mup3*(p5mup4*F(512)-p1mup4*F(519)
     &   +p4mup4*F(861)+p6mup4*F(863)-p2mup4*F(869))+p6mup3*(p5mup4*F
     &   (513)+p4mup4*F(862)+p6mup4*F(864)-p1mup4*F(867)-p2mup4*F(870
     &   )))
       Fa16(184) = 4*mup3mup4*F(871)+8*(-(p5mup4*p6mup3*F(591))-p5mup3
     &   *(p5mup4*F(589)+p6mup4*F(590)-p1mup4*F(592)-p2mup4*F(593))-p
     &   6mup3*p6mup4*F(618)+p1mup4*p6mup3*F(625)+p2mup4*p6mup3*F(628
     &   )-p4mup4*(p5mup3*F(588)+p6mup3*F(617)-p2mup3*F(873))+p2mup3*
     &   p5mup4*F(874)+p2mup3*p6mup4*F(875)-p1mup4*p2mup3*F(877)+p1mu
     &   p3*(p5mup4*F(583)+p4mup4*F(872)+p6mup4*F(876)-p1mup4*F(878)-
     &   p2mup4*F(879))-p2mup3*p2mup4*F(880)-p4mup3*(p5mup4*F(546)+p4
     &   mup4*F(881)+p6mup4*F(882)-p1mup4*F(883)-p2mup4*F(884)))
       Fa16(185) = mup3mup4*F(885)-4*(-(p5mup4*p6mup3*F(591))-p6mup3*p
     &   6mup4*F(618)+p1mup4*p6mup3*F(625)+p1mup3*(p5mup4*F(583)+p2mu
     &   p4*F(854)+p4mup4*F(872)+p6mup4*F(876)-p1mup4*F(878))-p4mup4*
     &   (p5mup3*F(588)+p6mup3*F(617)-p2mup3*F(884))-p1mup4*p2mup3*F(
     &   886)-p4mup3*(p5mup4*F(546)+p4mup4*F(881)+p6mup4*F(882)-p1mup
     &   4*F(883)-p2mup4*F(887))-p5mup3*(p5mup4*F(589)+p6mup4*F(590)-
     &   p1mup4*F(592)-p2mup4*F(888))+p2mup4*p6mup3*F(889)-p2mup3*p5m
     &   up4*F(890)-p2mup3*p6mup4*F(891)+p2mup3*p2mup4*F(892))
       Fa16(186) = mup3mup4*F(893)-4*(p1mup4*(p4mup3*F(860)+p5mup3*F(8
     &   63)+p6mup3*F(864)-p2mup3*F(870))+p1mup3*(p5mup4*F(499)-p2mup
     &   4*F(855)+p1mup4*F(856)+p4mup4*F(894)+p6mup4*F(895))-p4mup4*p
     &   5mup3*F(899)-p4mup4*p6mup3*F(900)-p5mup3*p5mup4*F(901)-p5mup
     &   3*p6mup4*F(902)-p5mup4*p6mup3*F(903)-p6mup3*p6mup4*F(904)+p2
     &   mup4*p6mup3*F(906)-p4mup3*(p4mup4*F(896)+p5mup4*F(897)+p6mup
     &   4*F(898)-p2mup4*F(907))+p2mup4*p5mup3*F(908)+p2mup3*(p4mup4*
     &   F(889)+p5mup4*F(905)+p6mup4*F(906)+p2mup4*F(909)))
       Fa16(187) = 2*(p1mup4*F(914)+p6mup4*F(915)-p4mup4*F(916)+p5mup4
     &   *F(917))+p2mup4*F(918)
       Fa16(188) = p4mup3*F(922)+2*(p5mup3*F(919)+p1mup3*F(920)+p6mup3
     &   *F(921)+p2mup3*F(923))
       Fa16(189) = -(mup3mup4*F(924))+2*(p5mup4*p6mup3*F(587)-p2mup4*p
     &   6mup3*F(862)+p1mup3*p2mup4*F(866)+p2mup3*p2mup4*F(868)-p2mup
     &   3*p4mup4*F(883)+p2mup3*p5mup4*F(925)+p2mup3*p6mup4*F(926)+p4
     &   mup3*(p5mup4*F(585)-p2mup4*F(859)+p4mup4*F(929)+p6mup4*F(930
     &   ))+p4mup4*p6mup3*F(932)+p5mup3*(p5mup4*F(586)-p2mup4*F(861)+
     &   p4mup4*F(931)+p6mup4*F(933))+p6mup3*p6mup4*F(934)-p1mup3*p5m
     &   up4*F(935)-p1mup3*p6mup4*F(936)-p1mup3*p4mup4*F(937)+p1mup4*
     &   (p5mup3*F(514)+p2mup3*F(866)+p4mup3*F(927)+p6mup3*F(928)-p1m
     &   up3*F(938)))
       Fa16(190) = mup3mup4*F(939)-4*(p2mup3*(p5mup4*F(583)+p4mup4*F(8
     &   72)+p6mup4*F(876)-p1mup4*F(878)-p2mup4*F(879))-p4mup3*(p5mup
     &   4*F(585)-p2mup4*F(883)+p4mup4*F(929)+p6mup4*F(930)-p1mup4*F(
     &   937))+p1mup3*(p5mup4*F(584)-p2mup4*F(878)+p4mup4*F(940)+p6mu
     &   p4*F(941)+p1mup4*F(942))-p5mup3*(p5mup4*F(586)-p2mup4*F(592)
     &   +p4mup4*F(931)+p6mup4*F(933)-p1mup4*F(943))-p6mup3*(p5mup4*F
     &   (587)-p2mup4*F(625)+p4mup4*F(932)+p6mup4*F(934)-p1mup4*F(944
     &   )))
       Fa16(191) = 2*mup3mup4*F(945)+p6mup3*p6mup4*F(959)+4*(p5mup4*p6
     &   mup3*F(657)-p1mup4*p6mup3*F(665)+p2mup4*(p5mup3*F(588)+p6mup
     &   3*F(617)-p1mup3*F(872)+p4mup3*F(881)-p2mup3*F(884))+p1mup4*p
     &   2mup3*F(946)-p1mup3*p5mup4*F(948)+p1mup3*p6mup4*F(949)-p1mup
     &   3*p1mup4*F(950)-p2mup3*p5mup4*F(951)-p2mup3*p6mup4*F(952)+p4
     &   mup3*p5mup4*F(954)+p4mup3*p6mup4*F(955)+p5mup3*p5mup4*F(957)
     &   +p5mup3*p6mup4*F(958)-p1mup4*p4mup3*F(960)-p1mup4*p5mup3*F(9
     &   61)+p4mup4*(p6mup3*F(656)-p1mup3*F(947)+p4mup3*F(953)+p5mup3
     &   *F(956)-p2mup3*F(962)))
       Fa16(192) = 2*mup3mup4*F(963)+p4mup3*p6mup4*F(967)+4*(p5mup4*p6
     &   mup3*F(542)+p1mup3*p5mup4*F(582)+p2mup3*p5mup4*F(591)+p2mup3
     &   *p4mup4*F(617)+p2mup4*p6mup3*F(618)+p2mup3*p6mup4*F(618)-p1m
     &   up4*p2mup3*F(625)-p2mup3*p2mup4*F(628)-p1mup3*p2mup4*F(876)+
     &   p1mup4*p6mup3*F(934)-p1mup3*p1mup4*F(941)+p1mup3*p4mup4*F(96
     &   4)+p1mup3*p6mup4*F(965)+p4mup3*(p5mup4*F(536)+p2mup4*F(882)+
     &   p1mup4*F(930)+p4mup4*F(966))+p4mup4*p6mup3*F(969)+p5mup3*(p5
     &   mup4*F(540)+p2mup4*F(590)+p1mup4*F(933)+p4mup4*F(968)+p6mup4
     &   *F(970))+p6mup3*p6mup4*F(971))
       Return
       End
