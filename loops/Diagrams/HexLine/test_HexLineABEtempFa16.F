c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe16(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa16)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup5, p2mup2, p2mup3, 
     -          p2mup4, p2mup5, p3mup2, p3mup3, p3mup4, p3mup5, 
     -          p4mup2, p4mup3, p4mup4, p4mup5, p5mup2, p5mup3, 
     -          p5mup4, p5mup5, p6mup2, p6mup3, p6mup4, p6mup5
       Complex*16   mup2mup3, mup2mup4, mup2mup5, mup3mup4, mup3mup5, 
     -          mup4mup5
        common/test_FhlFunctionsAbe/F
       COMMON/Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45
     -   ,s56,s16,s123,s234,s345
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
      Complex*16 Fa(445),F(5930)
      Real*8 P(11674) 
        COMMON/test_PFunctionsAbe/P
       Complex*16 Fa16(166:176)
       COMMON/test_FaFunctionsAbe/Fa
       Fa16(166) = -8*p3mup2*s12*(p6mup3*F(813)+3*((p1mup3+p2mup3)*F(8
     -   14)-p3mup3*F(815)+p4mup3*F(816)))+(p1mup2+p2mup2)*((p1mup3+p
     -   2mup3)*F(876)+p3mup3*F(877)+p6mup3*F(878)+p4mup3*F(879))
       Fa16(167) = p3mup4*F(880)+4*(p4mup4*F(881)-(p1mup4+p2mup4)*F(88
     -   2))+p6mup4*F(883)
       Fa16(168) = p3mup2*(p3mup4*F(884)+p4mup4*F(885)+(p1mup4+p2mup4)
     -   *F(889)-8*p6mup4*F(890))-p4mup2*(p6mup4*F(888)+8*((p1mup4+p2
     -   mup4)*F(886)-p3mup4*F(887)-p4mup4*F(891)))
       Fa16(169) = 4*mup3mup4*F(892)+p3mup3*p4mup4*F(894)+p4mup3*(p6mu
     -   p4*F(895)+(p1mup4+p2mup4)*F(896)+p3mup4*F(897))+p3mup4*p6mup
     -   3*F(898)+p4mup4*p6mup3*F(899)+p1mup4*p3mup3*F(901)+p2mup4*p3
     -   mup3*F(901)+p3mup3*p3mup4*F(902)+p6mup3*p6mup4*F(904)+p4mup3
     -   *p4mup4*F(905)+p1mup3*(p4mup4*F(893)+p3mup4*F(900)+8*p6mup4*
     -   F(903)+(p1mup4+p2mup4)*F(906))+p2mup3*(p4mup4*F(893)+p3mup4*
     -   F(900)+8*p6mup4*F(903)+(p1mup4+p2mup4)*F(906))+p1mup4*p6mup3
     -   *F(907)+p2mup4*p6mup3*F(907)-8*p3mup3*p6mup4*F(908)
       Fa16(170) = p1mup2*(8*p6mup4*F(890)+p3mup4*F(909)+p4mup4*F(910)
     -   +(p1mup4+p2mup4)*F(911))+p2mup2*(8*p6mup4*F(890)+p3mup4*F(90
     -   9)+p4mup4*F(910)+(p1mup4+p2mup4)*F(911))-8*p4mup2*s12*(p6mup
     -   4*F(912)+3*((p1mup4+p2mup4)*F(913)+p3mup4*F(914)+p4mup4*F(91
     -   5)))
       Fa16(171) = (p1mup3+p2mup3)*F(921)+p3mup3*F(922)+p6mup3*F(923)+
     -   p4mup3*F(924)
       Fa16(172) = 4*p6mup4*F(890)+p3mup4*F(925)+p4mup4*F(926)+(p1mup4
     -   +p2mup4)*F(927)
       Fa16(173) = p1mup2*p3mup3*p4mup5*F(931)+p2mup2*p3mup3*p4mup5*F(
     -   931)+p1mup2*p4mup3*p6mup5*F(932)+p2mup2*p4mup3*p6mup5*F(932)
     -   +8*(mup3mup5*(p1mup2+p2mup2)*F(929)+mup3mup5*p4mup2*s12*F(93
     -   3)+mup3mup5*p3mup2*s12*F(938))+32*(p1mup5*p3mup2*p4mup3*s12*
     -   F(950)+p2mup5*p3mup2*p4mup3*s12*F(950)+p3mup2*p3mup5*p4mup3*
     -   s12*F(954))+p1mup3*(p2mup2*p4mup5*F(930)+p1mup2*(p4mup5*F(93
     -   0)+16*(p3mup5*F(936)+p6mup5*F(940)-(p1mup5+p2mup5)*F(952)))-
     -   16*(p3mup2*p4mup5*s12*F(814)-p2mup2*p3mup5*F(936)-p2mup2*p6m
     -   up5*F(940)+p4mup2*p6mup5*s12*F(945)+p1mup5*p4mup2*s12*F(948)
     -   +p2mup5*p4mup2*s12*F(948)-p3mup5*p4mup2*s12*F(949)+p1mup5*p2
     -   mup2*F(952)+p2mup2*p2mup5*F(952)+p4mup2*p4mup5*s12*F(957)))+
     -   p2mup3*(p2mup2*p4mup5*F(930)+p1mup2*(p4mup5*F(930)+16*(p3mup
     -   5*F(936)+p6mup5*F(940)-(p1mup5+p2mup5)*F(952)))-16*(p3mup2*p
     -   4mup5*s12*F(814)-p2mup2*p3mup5*F(936)-p2mup2*p6mup5*F(940)+p
     -   4mup2*p6mup5*s12*F(945)+p1mup5*p4mup2*s12*F(948)+p2mup5*p4mu
     -   p2*s12*F(948)-p3mup5*p4mup2*s12*F(949)+p1mup5*p2mup2*F(952)+
     -   p2mup2*p2mup5*F(952)+p4mup2*p4mup5*s12*F(957)))+16*(p3mup2*p
     -   3mup3*p4mup5*s12*F(815)+p1mup2*p1mup5*p4mup3*F(934)+p1mup5*p
     -   2mup2*p4mup3*F(934)+p1mup2*p2mup5*p4mup3*F(934)+p2mup2*p2mup
     -   5*p4mup3*F(934)-p1mup2*p3mup5*p4mup3*F(935)-p2mup2*p3mup5*p4
     -   mup3*F(935)-p1mup2*p1mup5*p3mup3*F(937)-p1mup5*p2mup2*p3mup3
     -   *F(937)-p1mup2*p2mup5*p3mup3*F(937)-p2mup2*p2mup5*p3mup3*F(9
     -   37)+p1mup2*p4mup3*p4mup5*F(939)+p2mup2*p4mup3*p4mup5*F(939)-
     -   p1mup2*p1mup5*p6mup3*F(941)-p1mup5*p2mup2*p6mup3*F(941)-p1mu
     -   p2*p2mup5*p6mup3*F(941)-p2mup2*p2mup5*p6mup3*F(941)-p1mup2*p
     -   4mup5*p6mup3*F(942)-p2mup2*p4mup5*p6mup3*F(942)+p1mup2*p3mup
     -   3*p6mup5*F(943)+p2mup2*p3mup3*p6mup5*F(943)+p1mup2*p3mup5*p6
     -   mup3*F(944)+p2mup2*p3mup5*p6mup3*F(944)+p3mup2*p4mup3*p6mup5
     -   *s12*F(946)-p3mup3*p4mup2*p6mup5*s12*F(947)-p1mup5*p3mup3*p4
     -   mup2*s12*F(951)-p2mup5*p3mup3*p4mup2*s12*F(951)+p1mup2*p3mup
     -   3*p3mup5*F(953)+p2mup2*p3mup3*p3mup5*F(953)-p3mup3*p3mup5*p4
     -   mup2*s12*F(955)-p4mup2*p4mup3*p6mup5*s12*F(956)+p3mup2*p4mup
     -   3*p4mup5*s12*F(958)-p3mup3*p4mup2*p4mup5*s12*F(959)+p4mup2*p
     -   6mup3*p6mup5*s12*F(960)+p1mup2*p6mup3*p6mup5*F(961)+p2mup2*p
     -   6mup3*p6mup5*F(961)+p4mup2*p4mup5*p6mup3*s12*F(962)+p1mup5*p
     -   4mup2*p4mup3*s12*F(963)+p2mup5*p4mup2*p4mup3*s12*F(963)+p4mu
     -   p2*p4mup3*p4mup5*s12*F(964)+p3mup5*p4mup2*p4mup3*s12*F(965)-
     -   p1mup5*p4mup2*p6mup3*s12*F(966)-p2mup5*p4mup2*p6mup3*s12*F(9
     -   66)-p3mup5*p4mup2*p6mup3*s12*F(967)+p3mup2*p4mup5*p6mup3*s12
     -   *F(968))
       Fa16(174) = p3mup5*F(970)+4*(p4mup5*F(969)-(p1mup5+p2mup5)*F(97
     -   1)-p6mup5*F(972))
       Fa16(175) = p4mup2*((p1mup5+p2mup5)*F(974)-8*p3mup5*F(975)+p6mu
     -   p5*F(976)+p4mup5*F(979))+8*p3mup2*(p4mup5*F(973)-(p1mup5+p2m
     -   up5)*F(977)+p6mup5*F(978)+p3mup5*F(980))
       Fa16(176) = -4*mup3mup5*F(929)+p1mup3*p4mup5*F(981)+p2mup3*p4mu
     -   p5*F(981)+p3mup3*p4mup5*F(982)+p4mup3*(p6mup5*F(983)-8*(p4mu
     -   p5*F(939)-(p1mup5+p2mup5)*F(984)))-8*((p1mup3+p2mup3)*p6mup5
     -   *F(940)-(p1mup5+p2mup5)*p6mup3*F(941)-p4mup5*p6mup3*F(942)+p
     -   3mup3*p6mup5*F(943)-(p1mup3+p2mup3)*(p1mup5+p2mup5)*F(952)-p
     -   3mup5*(p4mup3*F(935)-(p1mup3+p2mup3)*F(936)-p6mup3*F(944)-p3
     -   mup3*F(953))+p6mup3*p6mup5*F(961)-(p1mup5+p2mup5)*p3mup3*F(9
     -   85))
       Return
       End
