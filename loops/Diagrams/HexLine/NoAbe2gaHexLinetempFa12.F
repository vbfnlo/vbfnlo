c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa12(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa12)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
        common/NoAbe2gaHexLineFhlFunctions/F
       COMMON/NoAbe2gaHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,
     -   s12,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(872),F(9706)
      Real*8 P(9687) 
        COMMON/NoAbe2gaHexLinePFunctions/P
       Complex*16 Fa12(265:288)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa12(265) = p1mup2*p3mup5*p3mup6*F(825)+p2mup2*p3mup5*p3mup6*F(
     -   826)-p3mup6*p4mup5*p5mup2*F(884)-p3mup5*p4mup6*p5mup2*F(884)
     -   -p4mup5*p4mup6*p5mup2*F(884)+p3mup2*p3mup5*p3mup6*F(885)+p3m
     -   up5*p3mup6*p4mup2*F(885)+p3mup2*p3mup6*p4mup5*F(885)+p3mup2*
     -   p3mup5*p4mup6*F(885)+p3mup2*p4mup5*p4mup6*F(885)-(p3mup5+p4m
     -   up5)*(p3mup6+p4mup6)*p6mup2*F(886)
       Fa12(266) = p2mup2*F(888)+p5mup2*F(889)+2*(p1mup2*F(745)+(p3mup
     -   2+p4mup2)*F(747)+p6mup2*F(890))
       Fa12(267) = -4*mup2mup6*mup3mup5*F(891)+mup3mup6*(p6mup2*p6mup5
     -   *F(905)+p2mup2*p6mup5*F(912)+p2mup5*p6mup2*F(919)+p5mup2*(4*
     -   p5mup5*F(895)+p6mup5*F(896)+p2mup5*F(914)+(p3mup5+p4mup5)*F(
     -   924))+p3mup5*p6mup2*F(926)+p4mup5*p6mup2*F(926)+p1mup5*(p5mu
     -   p2*F(893)+p6mup2*F(903)-4*(p1mup2*F(899)-p2mup2*F(908)+(p3mu
     -   p2+p4mup2)*F(931)))+p2mup5*p3mup2*F(932)+p2mup5*p4mup2*F(932
     -   )+p3mup2*p6mup5*F(934)+p4mup2*p6mup5*F(934)+4*(-(p1mup2*p5mu
     -   p5*F(900))+p1mup2*p6mup5*F(901)+p5mup5*p6mup2*F(904)+p2mup2*
     -   p2mup5*F(909)+p2mup2*p5mup5*F(910)+p1mup2*p2mup5*F(916)+p1mu
     -   p2*p3mup5*F(925)+p1mup2*p4mup5*F(925)-p3mup2*p5mup5*F(933)-p
     -   4mup2*p5mup5*F(933)+p2mup2*p3mup5*F(935)+p2mup2*p4mup5*F(935
     -   )+mup2mup5*F(939)-(p3mup2+p4mup2)*(p3mup5+p4mup5)*F(940)))+m
     -   up3mup5*(p5mup6*p6mup2*F(906)-p1mup6*p6mup2*F(907)-p2mup2*p5
     -   mup6*F(911)+p2mup6*p6mup2*F(920)-p1mup6*p2mup2*F(922)+p5mup2
     -   *(p5mup6*F(892)+p1mup6*F(894)+4*p6mup6*F(897)+p2mup6*F(915)-
     -   (p3mup6+p4mup6)*F(927))+p1mup2*(p1mup6*F(898)+p5mup6*F(902)+
     -   4*p6mup6*F(917)+p2mup6*F(918)+(p3mup6+p4mup6)*F(928))+p3mup6
     -   *p6mup2*F(929)+p4mup6*p6mup2*F(929)+4*(p2mup2*p6mup6*F(913)+
     -   p6mup2*p6mup6*F(921)+p1mup6*p3mup2*F(923)+p1mup6*p4mup2*F(92
     -   3)-p2mup2*p2mup6*F(930)+p2mup2*p3mup6*F(936)+p2mup2*p4mup6*F
     -   (936)-p3mup2*p5mup6*F(937)-p4mup2*p5mup6*F(937)-p3mup2*p6mup
     -   6*F(938)-p4mup2*p6mup6*F(938)+p3mup2*p3mup6*F(941)+p3mup6*p4
     -   mup2*F(941)+p3mup2*p4mup6*F(941)+p4mup2*p4mup6*F(941)+p2mup6
     -   *(p3mup2+p4mup2)*F(942)))
       Fa12(268) = mup3mup5*(p2mup2*F(943)+4*(p5mup2*F(944)+p1mup2*F(9
     -   45)+p6mup2*F(946)-(p3mup2+p4mup2)*F(947)))
       Fa12(269) = 4*mup2mup6*mup3mup5*F(948)+mup3mup6*p1mup5*p2mup2*F
     -   (949)+mup3mup6*(p1mup5*p6mup2*F(961)+p5mup5*p6mup2*F(962)+p2
     -   mup5*p6mup2*F(969)+p1mup5*p5mup2*F(972)+p2mup2*(p5mup5*F(951
     -   )+p6mup5*F(953)-4*p2mup5*F(965)+(p3mup5+p4mup5)*F(977))+p3mu
     -   p5*p6mup2*F(978)+p4mup5*p6mup2*F(978)+p3mup2*p5mup5*F(988)+p
     -   4mup2*p5mup5*F(988)+p3mup2*p6mup5*F(989)+p4mup2*p6mup5*F(989
     -   )+p3mup5*p5mup2*F(990)+p4mup5*p5mup2*F(990)+p2mup5*p5mup2*F(
     -   992)+p1mup2*(p5mup5*F(958)-4*(p6mup5*F(959)+p1mup5*F(975))+(
     -   p3mup5+p4mup5)*F(991)+p2mup5*F(994))-4*(p5mup2*p5mup5*F(973)
     -   +p5mup2*p6mup5*F(974)+p6mup2*p6mup5*F(976)+p1mup5*(p3mup2+p4
     -   mup2)*F(986)+p2mup5*p3mup2*F(987)+p2mup5*p4mup2*F(987)-mup2m
     -   up5*F(996)+(p3mup2+p4mup2)*(p3mup5+p4mup5)*F(998)))+mup3mup5
     -   *(p3mup2*p6mup6*F(995)+p4mup2*p6mup6*F(995)+4*(p5mup2*p6mup6
     -   *F(956)+p1mup2*p6mup6*F(960)+p6mup2*p6mup6*F(964)+p2mup6*p5m
     -   up2*F(967)+p1mup2*p2mup6*F(968)-p2mup6*p6mup2*F(970)+p1mup6*
     -   (p2mup2*F(950)+p5mup2*F(955)+p1mup2*F(957)+p6mup2*F(963)+(p3
     -   mup2+p4mup2)*F(971))+p3mup6*p5mup2*F(979)+p4mup6*p5mup2*F(97
     -   9)-p2mup2*(p5mup6*F(952)-p6mup6*F(954)-p2mup6*F(966)+(p3mup6
     -   +p4mup6)*F(980))-p1mup2*p3mup6*F(981)-p1mup2*p4mup6*F(981)-p
     -   3mup6*p6mup2*F(982)-p4mup6*p6mup2*F(982)-p5mup2*p5mup6*F(983
     -   )-p1mup2*p5mup6*F(984)-p5mup6*p6mup2*F(985)+p2mup6*p3mup2*F(
     -   993)+p2mup6*p4mup2*F(993)+p3mup2*p5mup6*F(997)+p4mup2*p5mup6
     -   *F(997)+(p3mup2+p4mup2)*(p3mup6+p4mup6)*F(999)))
       Fa12(270) = 4*(mup2mup6*mup3mup5*F(1000)+mup3mup6*p1mup5*p2mup2
     -   *F(1001))+mup3mup6*(p1mup2*p5mup5*F(1014)+p5mup5*p6mup2*F(10
     -   19)+p1mup2*p2mup5*F(1027)+p2mup5*p6mup2*F(1028)+p5mup2*(-4*(
     -   p5mup5*F(1010)-p6mup5*F(1011))-p2mup5*F(1024)-(p3mup5+p4mup5
     -   )*F(1034))+p1mup2*p3mup5*F(1036)+p1mup2*p4mup5*F(1036)+p3mup
     -   5*p6mup2*F(1037)+p4mup5*p6mup2*F(1037)+p1mup5*(p5mup2*F(1008
     -   )+p6mup2*F(1018)+4*(p1mup2*F(1013)-(p3mup2+p4mup2)*F(1040)))
     -   +p3mup2*p6mup5*F(1043)+p4mup2*p6mup5*F(1043)+4*(p1mup2*p6mup
     -   5*F(1015)+p6mup2*p6mup5*F(1020)+p2mup2*(p5mup5*F(1003)+p6mup
     -   5*F(1005)+p2mup5*F(1023)-(p3mup5+p4mup5)*F(1035))-p2mup5*p3m
     -   up2*F(1041)-p2mup5*p4mup2*F(1041)-p3mup2*p5mup5*F(1042)-p4mu
     -   p2*p5mup5*F(1042)+mup2mup5*F(1048)-(p3mup2+p4mup2)*(p3mup5+p
     -   4mup5)*F(1050)))+mup3mup5*(-(p5mup2*p6mup6*F(1012))-p1mup2*p
     -   6mup6*F(1017)+p2mup2*(-(p6mup6*F(1006))-4*(p5mup6*F(1004)-p2
     -   mup6*F(1025)+(p3mup6+p4mup6)*F(1039)))-p3mup2*p6mup6*F(1049)
     -   -p4mup2*p6mup6*F(1049)-4*(p5mup2*p5mup6*F(1007)+p1mup2*p5mup
     -   6*F(1016)+p5mup6*p6mup2*F(1021)+p6mup2*p6mup6*F(1022)-p2mup6
     -   *p5mup2*F(1026)-p1mup2*p2mup6*F(1031)-p1mup6*(p2mup2*F(1002)
     -   -p5mup2*F(1009)+(p3mup2+p4mup2)*F(1029)-p1mup2*F(1030)+p6mup
     -   2*F(1032))-p2mup6*p6mup2*F(1033)-p3mup6*p5mup2*F(1038)-p4mup
     -   6*p5mup2*F(1038)-p1mup2*p3mup6*F(1044)-p1mup2*p4mup6*F(1044)
     -   -p3mup6*p6mup2*F(1045)-p4mup6*p6mup2*F(1045)-p2mup6*p3mup2*F
     -   (1046)-p2mup6*p4mup2*F(1046)-p3mup2*p5mup6*F(1047)-p4mup2*p5
     -   mup6*F(1047)-(p3mup2+p4mup2)*(p3mup6+p4mup6)*F(1051)))
       Fa12(271) = mup5mup6*F(1052)+p2mup5*p5mup6*F(1059)+p5mup5*p6mup
     -   6*F(1062)+p5mup6*p6mup5*F(1063)+p2mup5*p6mup6*F(1067)+p3mup5
     -   *p5mup6*F(1069)+p4mup5*p5mup6*F(1069)+p3mup5*p6mup6*F(1070)+
     -   p4mup5*p6mup6*F(1070)+p1mup6*(p2mup5*F(1056)+p5mup5*F(1057)+
     -   p6mup5*F(1058)+(p3mup5+p4mup5)*F(1071))-4*(p2mup6*p3mup5*F(1
     -   072)+p2mup6*p4mup5*F(1072))+p1mup5*(p1mup6*F(1053)+p2mup6*F(
     -   1054)+p6mup6*F(1055)+p5mup6*F(1065)-2*(p3mup6+p4mup6)*F(1074
     -   ))+p2mup5*p3mup6*F(1075)+p2mup5*p4mup6*F(1075)+2*(p2mup6*p5m
     -   up5*F(1060)+p2mup6*p6mup5*F(1061)+p6mup5*p6mup6*F(1064)-p2mu
     -   p5*p2mup6*F(1066)+p5mup5*p5mup6*F(1068)+p3mup6*p5mup5*F(1073
     -   )+p4mup6*p5mup5*F(1073)+p3mup6*p6mup5*F(1076)+p4mup6*p6mup5*
     -   F(1076))+8*(p3mup5+p4mup5)*(p3mup6+p4mup6)*F(1077)
       Fa12(272) = mup3mup6*(p1mup5*F(1080)+p5mup5*F(1082)+p6mup5*F(10
     -   84)+p2mup5*F(1086)+(p3mup5+p4mup5)*F(1088))+mup3mup5*(p6mup6
     -   *F(1085)+2*(p1mup6*F(1081)-p5mup6*F(1083)-p2mup6*F(1087)-(p3
     -   mup6+p4mup6)*F(1089)))
       Fa12(273) = mup3mup6*(p1mup5*F(1090)+p5mup5*F(1092)+p6mup5*F(10
     -   94)+p2mup5*F(1096)+(p3mup5+p4mup5)*F(1098))+mup3mup5*(-(p6mu
     -   p6*F(1095))+2*(p1mup6*F(1091)+p5mup6*F(1093)+p2mup6*F(1097)+
     -   (p3mup6+p4mup6)*F(1099)))
       Fa12(274) = mup5mup6*F(1101)+p2mup5*p6mup6*F(1107)+p2mup6*p6mup
     -   5*F(1109)-p6mup5*p6mup6*F(1110)+p2mup5*p2mup6*F(1112)+p3mup5
     -   *p6mup6*F(1113)+p4mup5*p6mup6*F(1113)+p2mup5*p5mup6*F(1114)+
     -   2*(p2mup6*p5mup5*F(1108)+p5mup5*p5mup6*F(1115))-p5mup5*p6mup
     -   6*F(1116)+p1mup6*(p2mup5*F(1104)+p5mup5*F(1105)+p6mup5*F(110
     -   6)-(p3mup5+p4mup5)*F(1117))+p3mup6*p5mup5*F(1118)+p4mup6*p5m
     -   up5*F(1118)+p2mup5*p3mup6*F(1120)+p2mup5*p4mup6*F(1120)+p3mu
     -   p6*p6mup5*F(1121)+p4mup6*p6mup5*F(1121)+p1mup5*(p1mup6*F(110
     -   2)+p6mup6*F(1103)+p2mup6*F(1111)+2*(p3mup6+p4mup6)*F(1119)+p
     -   5mup6*F(1122))+p5mup6*p6mup5*F(1123)+8*(p3mup5*p3mup6*F(1124
     -   )+p3mup6*p4mup5*F(1124)+p3mup5*p4mup6*F(1124)+p4mup5*p4mup6*
     -   F(1124))+4*(p3mup5*p5mup6*F(1125)+p4mup5*p5mup6*F(1125))+p2m
     -   up6*(p3mup5+p4mup5)*F(1126)
       Fa12(275) = 8*(p3mup6*p4mup5*F(886)+p4mup5*p4mup6*F(886))+mup5m
     -   up6*F(1127)+p1mup5*(p1mup6*F(1128)+p2mup6*F(1129)+p5mup6*F(1
     -   130)+p6mup6*F(1131))+p1mup6*p5mup5*F(1132)+p1mup6*p6mup5*F(1
     -   133)+p2mup5*p5mup6*F(1134)+p2mup6*p5mup5*F(1136)+p2mup6*p6mu
     -   p5*F(1137)+p5mup5*p5mup6*F(1138)+p1mup6*p2mup5*F(1141)+p2mup
     -   5*p2mup6*F(1142)+4*p4mup5*p6mup6*F(1144)+p2mup6*p4mup5*F(114
     -   6)+p3mup5*(8*(p3mup6+p4mup6)*F(886)+4*p6mup6*F(1144)+2*(p5mu
     -   p6*F(1143)+p1mup6*F(1145))+p2mup6*F(1146))+2*(p2mup5*p6mup6*
     -   F(1135)+p5mup5*p6mup6*F(1139)-p6mup5*p6mup6*F(1140)+p4mup5*p
     -   5mup6*F(1143)+p1mup6*p4mup5*F(1145)-p3mup6*p5mup5*F(1147)-p4
     -   mup6*p5mup5*F(1147)-p1mup5*p3mup6*F(1148)-p1mup5*p4mup6*F(11
     -   48)+p3mup6*p6mup5*F(1149)+p4mup6*p6mup5*F(1149))+p5mup6*p6mu
     -   p5*F(1150)+p2mup5*(p3mup6+p4mup6)*F(1151)
       Fa12(276) = p1mup6*F(1154)+p2mup6*F(1155)+p5mup6*F(1156)+p6mup6
     -   *F(1157)+(p3mup6+p4mup6)*F(1158)
       Fa12(277) = p2mup2*p5mup6*F(1167)+p6mup2*p6mup6*F(1172)+p1mup6*
     -   p2mup2*F(1174)+p5mup2*(p1mup6*F(1161)+p6mup6*F(1162)-2*(p5mu
     -   p6*F(1160)-p2mup6*F(1169)-(p3mup6+p4mup6)*F(1176)))+p1mup2*(
     -   -(p6mup6*F(1170))+2*(p1mup6*F(1163)+p5mup6*F(1164)+p2mup6*F(
     -   1171)+(p3mup6+p4mup6)*F(1177)))+p2mup2*p2mup6*F(1179)+p2mup2
     -   *p3mup6*F(1180)+p2mup2*p4mup6*F(1180)+p3mup2*p6mup6*F(1182)+
     -   p4mup2*p6mup6*F(1182)+2*(mup2mup6*F(1159)+p5mup6*p6mup2*F(11
     -   65)+p1mup6*p6mup2*F(1166)-p2mup2*p6mup6*F(1168)+p2mup6*p6mup
     -   2*F(1173)+p1mup6*p3mup2*F(1175)+p1mup6*p4mup2*F(1175)+p3mup6
     -   *p6mup2*F(1178)+p4mup6*p6mup2*F(1178)+p3mup2*p5mup6*F(1181)+
     -   p4mup2*p5mup6*F(1181)+p3mup2*p3mup6*F(1183)+p3mup6*p4mup2*F(
     -   1183)+p3mup2*p4mup6*F(1183)+p4mup2*p4mup6*F(1183)+p2mup6*(p3
     -   mup2+p4mup2)*F(1184))
       Fa12(278) = p1mup6*p5mup2*F(1188)+p5mup2*p6mup6*F(1189)+p1mup2*
     -   (2*(p3mup6+p4mup6)*F(981)-p1mup6*F(1190)+p6mup6*F(1191))+p6m
     -   up2*p6mup6*F(1192)+p2mup6*(2*(p6mup2*F(970)-(p3mup2+p4mup2)*
     -   F(993))+p5mup2*F(1194)+p1mup2*F(1195))+p2mup2*(-(p1mup6*F(11
     -   86))+p6mup6*F(1187)+2*(p5mup6*F(952)-p2mup6*F(966)+(p3mup6+p
     -   4mup6)*F(1198)))+p5mup2*p5mup6*F(1199)+p3mup2*p6mup6*F(1201)
     -   +p4mup2*p6mup6*F(1201)+2*(p3mup6*p6mup2*F(982)+p4mup6*p6mup2
     -   *F(982)+p5mup6*p6mup2*F(985)+mup2mup6*F(1185)-p1mup6*p6mup2*
     -   F(1193)+p1mup6*p3mup2*F(1196)+p1mup6*p4mup2*F(1196)+p3mup6*p
     -   5mup2*F(1197)+p4mup6*p5mup2*F(1197)+p1mup2*p5mup6*F(1200)+p3
     -   mup2*p5mup6*F(1202)+p4mup2*p5mup6*F(1202)+(p3mup2+p4mup2)*(p
     -   3mup6+p4mup6)*F(1203))
       Fa12(279) = p2mup2*p6mup6*F(1207)+p5mup2*(2*(p5mup6*F(1007)+p1m
     -   up6*F(1009)-p2mup6*F(1026)-(p3mup6+p4mup6)*F(1038))+p6mup6*F
     -   (1208))+p1mup2*p6mup6*F(1210)+p1mup6*p6mup2*F(1215)+p2mup6*p
     -   6mup2*F(1216)+p3mup6*p6mup2*F(1218)+p4mup6*p6mup2*F(1218)+p3
     -   mup2*p6mup6*F(1221)+p4mup2*p6mup6*F(1221)+2*(p6mup2*p6mup6*F
     -   (1022)+mup2mup6*F(1204)-p2mup2*(p2mup6*F(1025)+p1mup6*F(1205
     -   ))+p1mup6*(p3mup2+p4mup2)*F(1211)+p1mup2*p1mup6*F(1212)+p1mu
     -   p2*p2mup6*F(1213)+p2mup2*p3mup6*F(1214)+p2mup2*p4mup6*F(1214
     -   )-p1mup2*p3mup6*F(1217)-p1mup2*p4mup6*F(1217)+p2mup6*p3mup2*
     -   F(1219)+p2mup6*p4mup2*F(1219)+p5mup6*(p6mup2*F(1021)+p2mup2*
     -   F(1206)+p1mup2*F(1209)+(p3mup2+p4mup2)*F(1220))+(p3mup2+p4mu
     -   p2)*(p3mup6+p4mup6)*F(1222))
       Fa12(280) = p1mup5*F(1225)+p2mup5*F(1226)+p5mup5*F(1227)+p6mup5
     -   *F(1228)+(p3mup5+p4mup5)*F(1229)
       Fa12(281) = mup5mup6*(p2mup2*F(1232)+(p3mup2+p4mup2)*F(1233)+p5
     -   mup2*F(1234))
       Fa12(282) = p3mup2*(4*(p1mup5*p1mup6*F(1249)+p1mup6*p2mup5*F(12
     -   50)+p1mup5*p2mup6*F(1251)+p2mup5*p2mup6*F(1252)+p1mup5*p5mup
     -   6*F(1255)+p2mup5*p5mup6*F(1256))+p5mup5*(p6mup6*F(1236)-p1mu
     -   p6*F(1247)-2*(p2mup6*F(1253)+p5mup6*F(1257)))-2*(p1mup5*p6mu
     -   p6*F(1244)+p2mup5*p6mup6*F(1245)-p6mup5*(p6mup6*F(1237)-2*(p
     -   1mup6*F(1248)-p2mup6*F(1254)-p5mup6*F(1258)))))+p4mup2*(4*(p
     -   1mup5*p1mup6*F(1249)+p1mup6*p2mup5*F(1250)+p1mup5*p2mup6*F(1
     -   251)+p2mup5*p2mup6*F(1252)+p1mup5*p5mup6*F(1255)+p2mup5*p5mu
     -   p6*F(1256))+p5mup5*(p6mup6*F(1236)-p1mup6*F(1247)-2*(p2mup6*
     -   F(1253)+p5mup6*F(1257)))-2*(p1mup5*p6mup6*F(1244)+p2mup5*p6m
     -   up6*F(1245)-p6mup5*(p6mup6*F(1237)-2*(p1mup6*F(1248)-p2mup6*
     -   F(1254)-p5mup6*F(1258)))))-2*(p2mup2*(p2mup6*(p5mup5*F(1242)
     -   +2*(p2mup5*F(1241)+p6mup5*F(1243)))+2*p1mup5*(-(p5mup6*F(123
     -   9))+p2mup6*F(1240)+p1mup6*F(1246)))+p5mup2*p5mup6*(-2*(p6mup
     -   5*F(1238)+p1mup5*F(1259)-p2mup5*F(1260))+p5mup5*F(1261)))
       Fa12(283) = -2*mup5mup6*p6mup2*F(1263)+mup2mup6*(p6mup5*F(1264)
     -   -2*(p5mup5*F(1262)-p1mup5*F(1265)-p2mup5*F(1266)))+4*(-(p2mu
     -   p2*(-(p5mup5*p6mup6*F(1270))+2*(p2mup5*p6mup6*F(1268)+p6mup5
     -   *p6mup6*F(1269)+p5mup5*p5mup6*F(1274))+p1mup6*p5mup5*F(1278)
     -   +4*(p2mup5*p5mup6*F(1273)+p5mup6*p6mup5*F(1275)+p1mup6*p6mup
     -   5*F(1279)-p1mup6*p2mup5*F(1280))))-2*p1mup5*(p2mup2*p6mup6*F
     -   (1267)-p5mup2*(p6mup6*F(1276)+2*(p2mup6*F(1284)+p1mup6*F(128
     -   7))))+p5mup2*(-(p5mup5*p6mup6*F(1271))+p1mup6*p5mup5*F(1282)
     -   +2*(p6mup5*p6mup6*F(1272)+p2mup5*p6mup6*F(1277)+p2mup6*p5mup
     -   5*F(1286))-4*(p1mup6*p6mup5*F(1281)+p2mup6*p6mup5*F(1283)-p2
     -   mup5*p2mup6*F(1285)-p1mup6*p2mup5*F(1288))))
       Fa12(284) = mup2mup6*p4mup5*F(1289)+p1mup2*(p5mup5*(-2*(p6mup6*
     -   F(1290)-p1mup6*F(1292))+4*(p2mup6*F(1296)+p5mup6*F(1299)))+4
     -   *(p2mup5*p6mup6*F(1295)+p4mup5*p6mup6*F(1298)+p6mup5*(-(p6mu
     -   p6*F(1291))+2*(p1mup6*F(1293)+p2mup6*F(1297)+p5mup6*F(1300))
     -   )+p1mup5*(p6mup6*F(1294)+2*(p5mup6*F(1301)-p2mup6*F(1305)-p1
     -   mup6*F(1308)))+2*(p2mup5*p5mup6*F(1302)-p4mup5*p5mup6*F(1303
     -   )-p2mup5*p2mup6*F(1306)-p4mup5*p4mup6*F(1307)-p1mup6*p2mup5*
     -   F(1309))))-4*p2mup2*p4mup5*(p6mup6*F(1304)+2*(p4mup6*F(1310)
     -   -p5mup6*F(1311)))
       Fa12(285) = mup2mup6*p3mup5*F(1289)+mup2mup5*p6mup6*F(1313)-2*(
     -   mup2mup5*p1mup6*F(1312)+mup2mup5*p2mup6*F(1314)+mup2mup5*p5m
     -   up6*F(1315))+4*(-(p2mup2*p3mup5*p6mup6*F(1304))+p3mup5*p5mup
     -   2*p6mup6*F(1317)+p4mup5*p5mup2*p6mup6*F(1317)-p2mup5*p6mup2*
     -   p6mup6*F(1319)+p1mup2*(p3mup5*p6mup6*F(1298)-2*(p3mup5*p5mup
     -   6*F(1303)+p1mup6*p3mup5*F(1318)+p1mup6*p4mup5*F(1318)+p2mup6
     -   *(p3mup5+p4mup5)*F(1320))))+12*(p3mup2*p3mup5*p6mup6*F(1321)
     -   +p3mup5*p4mup2*p6mup6*F(1321)+p3mup2*p4mup5*p6mup6*F(1321)+p
     -   4mup2*p4mup5*p6mup6*F(1321))+8*(p2mup2*p3mup5*p5mup6*F(1311)
     -   +p1mup6*p2mup5*p6mup2*F(1316)-p1mup6*p2mup2*p3mup5*F(1322)-p
     -   1mup6*p2mup2*p4mup5*F(1322)+p2mup2*p2mup6*p3mup5*F(1323)+p2m
     -   up2*p2mup6*p4mup5*F(1323)+p3mup5*p5mup2*p5mup6*F(1324)+p4mup
     -   5*p5mup2*p5mup6*F(1324)-p2mup6*p3mup2*p3mup5*F(1325)-p2mup6*
     -   p3mup5*p4mup2*F(1325)-p2mup6*p3mup2*p4mup5*F(1325)-p2mup6*p4
     -   mup2*p4mup5*F(1325)+p3mup2*p3mup5*p5mup6*F(1326)+p3mup5*p4mu
     -   p2*p5mup6*F(1326)+p3mup2*p4mup5*p5mup6*F(1326)+p4mup2*p4mup5
     -   *p5mup6*F(1326)-p2mup6*p3mup5*p5mup2*F(1327)-p2mup6*p4mup5*p
     -   5mup2*F(1327)-p1mup6*p3mup2*p3mup5*F(1328)-p1mup6*p3mup5*p4m
     -   up2*F(1328)-p1mup6*p3mup2*p4mup5*F(1328)-p1mup6*p4mup2*p4mup
     -   5*F(1328)+p1mup6*p3mup5*p5mup2*F(1329)+p1mup6*p4mup5*p5mup2*
     -   F(1329)+p2mup5*p5mup6*p6mup2*F(9593))
       Fa12(286) = p4mup6*(p1mup2*(-(p5mup5*F(1330))+2*(p3mup5*F(1307)
     -   -p1mup5*F(1332)-p2mup5*F(1333)))+p2mup2*(p5mup5*F(1331)+2*(p
     -   3mup5*F(1310)-p1mup5*F(1334)-p2mup5*F(1335))))
       Fa12(287) = p6mup2*(-(p5mup5*(p5mup6*F(1338)+p2mup6*F(1339)))+2
     -   *p2mup5*p2mup6*F(9589)+(p3mup5+p4mup5)*(-(p6mup6*F(1337))+2*
     -   (p1mup6*F(1336)+p5mup6*F(9595)+p2mup6*F(9600))))
       Fa12(288) = mup2mup5*p3mup6*F(1340)+mup2mup5*p4mup6*F(1340)+p1m
     -   up6*p5mup5*p6mup2*F(1341)-p5mup5*p6mup2*p6mup6*F(1344)+2*(p3
     -   mup6*p5mup2*p5mup5*F(1257)+p4mup6*p5mup2*p5mup5*F(1257)-p1mu
     -   p2*p3mup6*p5mup5*F(1330)+p2mup2*p3mup6*p5mup5*F(1331)-p1mup5
     -   *(2*(p4mup6*p5mup2*F(1255)+p3mup6*(p5mup2*F(1255)+p1mup2*F(1
     -   332)+p2mup2*F(1334)))+p6mup2*(p6mup6*F(1342)-2*(p2mup6*F(134
     -   3)+p5mup6*F(1346)+p1mup6*F(1349))))+p3mup2*p3mup6*p5mup5*F(1
     -   351)+p3mup6*p4mup2*p5mup5*F(1351)+p3mup2*p4mup6*p5mup5*F(135
     -   1)+p4mup2*p4mup6*p5mup5*F(1351))+4*(p1mup2*p3mup6*p4mup5*F(1
     -   307)+p2mup2*p3mup6*p4mup5*F(1310)-p2mup5*(p3mup6*p5mup2*F(12
     -   56)+p4mup6*p5mup2*F(1256)+p1mup2*p3mup6*F(1333)+p2mup2*p3mup
     -   6*F(1335))-p1mup2*p3mup6*p6mup5*F(1345)-p1mup2*p4mup6*p6mup5
     -   *F(1345)+p3mup6*p5mup2*p6mup5*F(1347)+p4mup6*p5mup2*p6mup5*F
     -   (1347)+p2mup2*p3mup6*p6mup5*F(1348)+p2mup2*p4mup6*p6mup5*F(1
     -   348)-p1mup5*(p3mup2+p4mup2)*(p3mup6+p4mup6)*F(1350)+p2mup5*(
     -   p3mup2+p4mup2)*(p3mup6+p4mup6)*F(1352)+(p3mup2+p4mup2)*(p3mu
     -   p6+p4mup6)*p6mup5*F(1353))
       Return
       End
