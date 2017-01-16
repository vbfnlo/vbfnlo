c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa15(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   15)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBoxCross2FhlFunctions/F
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(392),F(1853)
      Real*8 P(1355) 
        COMMON/HexBoxCross2PFunctions/P
       Complex*16 Fa15(267:285)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa15(267) = p1mup6*F(845)+p4mup6*F(846)-p5mup6*F(847)-p6mup6*F(
     &   848)+p2mup6*F(849)
       Fa15(268) = mup3mup6*F(851)+2*(-(p6mup3*p6mup6*F(854))-p5mup6*p
     &   6mup3*F(855)+p2mup3*(p6mup6*F(850)+p5mup6*F(857)+p4mup6*F(85
     &   8)+p1mup6*F(859)+p2mup6*F(860))+p2mup6*p6mup3*F(861)+p1mup3*
     &   p6mup6*F(863)-p1mup6*p6mup3*F(865)+p4mup3*p6mup6*F(866)-p5mu
     &   p3*(-(p6mup6*F(853))+p2mup6*F(856)+p1mup6*F(864)+p5mup6*F(86
     &   7))-p2mup6*p4mup3*F(868)+p1mup6*p4mup3*F(869)+p4mup3*p5mup6*
     &   F(871)+p1mup3*p2mup6*F(872)-p4mup6*(p6mup3*F(852)+p5mup3*F(8
     &   62)-p4mup3*F(870)+p1mup3*F(873))-p1mup3*p1mup6*F(874)-p1mup3
     &   *p5mup6*F(875))
       Fa15(269) = mup3mup6*F(876)+2*(-(p5mup3*p6mup6*F(880))-p6mup3*p
     &   6mup6*F(881)-p4mup3*p4mup6*F(882)-p4mup3*p6mup6*F(883)-p5mup
     &   6*p6mup3*F(884)-p4mup3*p5mup6*F(885)-p5mup3*p5mup6*F(886)+p1
     &   mup3*p2mup6*F(887)+p2mup6*p5mup3*F(888)+p2mup6*p6mup3*F(889)
     &   +p2mup6*p4mup3*F(890)+p2mup3*(p5mup6*F(877)+p4mup6*F(878)+p6
     &   mup6*F(879)+p1mup6*F(891)+p2mup6*F(892))+p1mup3*p4mup6*F(893
     &   )+p4mup6*p5mup3*F(894)+p4mup6*p6mup3*F(895)-p1mup3*p5mup6*F(
     &   896)-p1mup6*p5mup3*F(897)-p1mup6*p6mup3*F(898)-p1mup6*p4mup3
     &   *F(899)+p1mup3*p6mup6*F(900))
       Fa15(270) = p4mup6*F(903)+p5mup6*F(904)+p1mup6*F(905)+p2mup6*F(
     &   906)+p6mup6*F(907)
       Fa15(271) = p6mup3*F(908)-p2mup3*F(909)+p1mup3*F(910)-p4mup3*F(
     &   911)-p5mup3*F(912)
       Fa15(272) = p1mup3*F(913)-p4mup3*F(914)-p5mup3*F(915)-p6mup3*F(
     &   916)-p2mup3*F(917)
       Fa15(273) = mup3mup6*F(918)+2*(p1mup3*p6mup6*F(887)-p4mup3*p6mu
     &   p6*F(921)+p5mup6*p6mup3*F(922)-p1mup6*p6mup3*F(925)-p2mup6*p
     &   6mup3*F(927)-p6mup3*p6mup6*F(928)-p1mup6*p4mup3*F(929)-p4mup
     &   6*(-(p6mup3*F(919))+p5mup3*F(923)+p4mup3*F(930))-p2mup6*p4mu
     &   p3*F(931)+p2mup3*p6mup6*F(932)-p4mup3*p5mup6*F(933)-p5mup3*(
     &   p6mup6*F(920)+p1mup6*F(924)+p2mup6*F(926)+p5mup6*F(934)))
       Fa15(274) = p6mup3*F(936)+p1mup3*F(937)-p4mup3*F(938)+p5mup3*F(
     &   939)+p2mup3*F(940)
       Fa15(275) = mup3mup6*F(956)-2*(p1mup6*p6mup3*F(942)+p4mup6*p6mu
     &   p3*F(947)-p5mup3*p6mup6*F(949)-p4mup3*p6mup6*F(950)-p1mup3*(
     &   p6mup6*F(941)-p2mup6*F(944)-p4mup6*F(946)-p1mup6*F(948)-p5mu
     &   p6*F(953))+p5mup6*p6mup3*F(954)-p2mup3*(-(p6mup6*F(943))+p2m
     &   up6*F(945)+p5mup6*F(951)+p4mup6*F(952)+p1mup6*F(955))-p6mup3
     &   *p6mup6*F(957)-p2mup6*p6mup3*F(958))
       Fa15(276) = mup3mup6*F(959)+2*(p1mup3*p2mup6*F(960)-p1mup6*p2mu
     &   p3*F(961)-p1mup3*p4mup6*F(962)+p1mup6*p4mup3*F(963)-p2mup3*p
     &   5mup6*F(964)-p2mup3*p6mup6*F(965)+p1mup3*p5mup6*F(966)+p1mup
     &   3*p6mup6*F(967)-p2mup6*p5mup3*F(968)-p2mup6*p6mup3*F(969)-p2
     &   mup6*p4mup3*F(970)+p2mup3*p2mup6*F(971)-p2mup3*p4mup6*F(972)
     &   -p4mup6*p5mup3*F(973)-p4mup6*p6mup3*F(974)-p4mup3*p4mup6*F(9
     &   75)+p4mup3*p5mup6*F(976)+p1mup3*p1mup6*F(977)-p1mup6*p5mup3*
     &   F(978)-p5mup3*p6mup6*F(979)-p6mup3*p6mup6*F(980)-p4mup3*p6mu
     &   p6*F(981)-p5mup6*p6mup3*F(982)+p1mup6*p6mup3*F(983)+p5mup3*p
     &   5mup6*F(984))
       Fa15(277) = mup3mup6*F(1000)-2*(p1mup6*p6mup3*F(986)+p4mup6*p6m
     &   up3*F(992)-p1mup3*(p6mup6*F(985)+p5mup6*F(988)-p2mup6*F(989)
     &   -p4mup6*F(991)-p1mup6*F(993))-p5mup3*p6mup6*F(994)-p4mup3*(p
     &   5mup6*F(950)+p6mup6*F(995))-p2mup3*(-(p6mup6*F(987))+p2mup6*
     &   F(990)+p5mup6*F(996)+p4mup6*F(997)+p1mup6*F(998))-p5mup3*p5m
     &   up6*F(999)+p5mup6*p6mup3*F(1001)+p6mup3*p6mup6*F(1002)+p1mup
     &   6*p5mup3*F(1003)-p2mup6*p6mup3*F(1004)+p2mup6*p5mup3*F(1005)
     &   -p4mup6*p5mup3*F(1006))
       Fa15(278) = mup3mup6*F(1007)+2*(-(p1mup6*p2mup3*F(1009))+p1mup6
     &   *p4mup3*F(1011)-p2mup3*p5mup6*F(1012)-p2mup3*p6mup6*F(1013)-
     &   p2mup6*(p5mup3*F(1015)+p6mup3*F(1016)+p4mup3*F(1017)-p2mup3*
     &   F(1018))-p2mup3*p4mup6*F(1019)-p4mup6*p5mup3*F(1020)-p4mup6*
     &   p6mup3*F(1021)-p4mup3*p4mup6*F(1022)+p4mup3*p5mup6*F(1023)-p
     &   1mup6*p5mup3*F(1024)+p5mup3*p5mup6*F(1025)+p1mup6*p6mup3*F(1
     &   026)+p5mup3*p6mup6*F(1027)+p6mup3*p6mup6*F(1028)+p4mup3*p6mu
     &   p6*F(1029)+p1mup3*(p2mup6*F(1008)-p4mup6*F(1010)+p5mup6*F(10
     &   14)+p1mup6*F(1030)-p6mup6*F(1031))+p5mup6*p6mup3*F(1032))
       Fa15(279) = mup3mup6*F(1037)+2*(-(p1mup6*p6mup3*F(1039))-p6mup3
     &   *p6mup6*F(1044)-p4mup3*p6mup6*F(1045)-p2mup3*p6mup6*F(1046)+
     &   p1mup3*(p6mup6*F(1038)+p4mup6*F(1040)+p1mup6*F(1041)+p2mup6*
     &   F(1042)+p5mup6*F(1047))-p4mup6*p6mup3*F(1048)+p1mup6*p2mup3*
     &   F(1049)-p1mup6*p4mup3*F(1050)+p2mup3*p2mup6*F(1051)-p2mup6*p
     &   4mup3*F(1052)+p2mup6*p6mup3*F(1053)-p2mup3*p5mup6*F(1054)-p2
     &   mup3*p4mup6*F(1055)+p4mup3*p4mup6*F(1056)+p4mup3*p5mup6*F(10
     &   57)-p5mup6*p6mup3*F(1058)-p5mup3*(p6mup6*F(1043)+p1mup6*F(10
     &   59)+p4mup6*F(1060)+p2mup6*F(1061)+p5mup6*F(1062)))
       Fa15(280) = mup3mup6*F(1063)-2*(p1mup6*p6mup3*F(1065)-p1mup3*(p
     &   6mup6*F(1064)+p5mup6*F(1066)+p4mup6*F(1067)+p1mup6*F(1068)+p
     &   2mup6*F(1069))+p6mup3*p6mup6*F(1071)+p4mup3*p6mup6*F(1072)+p
     &   5mup6*p6mup3*F(1074)+p4mup3*p5mup6*F(1075)+p4mup6*p6mup3*F(1
     &   077)+p1mup6*p4mup3*F(1079)+p2mup6*p4mup3*F(1082)+p2mup3*(p5m
     &   up6*F(1046)+p6mup6*F(1073)-p1mup6*F(1078)-p2mup6*F(1081)+p4m
     &   up6*F(1084))-p4mup3*p4mup6*F(1085)+p5mup3*(p6mup6*F(1070)+p5
     &   mup6*F(1076)+p1mup6*F(1080)+p2mup6*F(1083)-p4mup6*F(1086))+p
     &   2mup6*p6mup3*F(1087))
       Fa15(281) = 4*mup3mup6*F(1092)+p4mup6*(-8*(p5mup3*F(1099)+p6mup
     &   3*F(1100))+p4mup3*F(1103))+p1mup3*(p1mup6*F(1093)+p2mup6*F(1
     &   095)+p4mup6*F(1098)-8*(p6mup6*F(1094)+p5mup6*F(1104)))-8*(-(
     &   p5mup3*p6mup6*F(1102))+p5mup6*(p6mup3*F(1105)+p4mup3*F(1106)
     &   )+p5mup3*p5mup6*F(1107)-p2mup6*p5mup3*F(1109)-p2mup6*p6mup3*
     &   F(1110)-p2mup6*p4mup3*F(1111)-p6mup3*p6mup6*F(1112)-p4mup3*p
     &   6mup6*F(1113)+p2mup3*(p6mup6*F(1096)+p1mup6*F(1097)-p2mup6*F
     &   (1108)+p5mup6*F(1114)+p4mup6*F(1115))-p1mup6*(p5mup3*F(1101)
     &   -p4mup3*F(1116)-p6mup3*F(1117)))
       Fa15(282) = -(p6mup3*F(1119))+p2mup3*F(1120)+p4mup3*F(1121)+p5m
     &   up3*F(1122)+p1mup3*F(1123)
       Fa15(283) = p6mup3*F(1124)-p2mup3*F(1125)+p4mup3*F(1126)-p1mup3
     &   *F(1127)+p5mup3*F(1128)
       Fa15(284) = p1mup3*F(1129)-p6mup3*F(1130)+p2mup3*F(1131)+p5mup3
     &   *F(1132)-p4mup3*F(1133)
       Fa15(285) = p2mup6*F(1139)+8*(p1mup6*F(1137)+p6mup6*F(1138)-p4m
     &   up6*F(1140)+p5mup6*F(1141))
       Return
       End
