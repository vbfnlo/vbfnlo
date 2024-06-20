c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim3(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,
     -   p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mu
     -   p3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup
     -   3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa3)
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
       Complex*16 Fa3(199:297)
       COMMON/FaPenOddGSim/Fa
       Fa3(199) = (mup3mup5*F(1138))/2.+p2mup5*p3mup3*F(1142)+p2mup5*p
     -   4mup3*F(1144)+p3mup3*p4mup5*F(1146)+p3mup5*p4mup3*F(1147)+p2
     -   mup3*(p3mup5*F(1141)+p4mup5*F(1143)+2*(p2mup5*F(1139)+p1mup5
     -   *F(1148)))+p1mup5*p4mup3*F(1149)+2*(p1mup3*p2mup5*F(1140)+p1
     -   mup3*p4mup5*F(1145)+p4mup3*p4mup5*F(1150)+(p1mup3*p3mup5*F(3
     -   44)+p1mup5*p3mup3*F(348)+p3mup3*p3mup5*F(351)+2*p1mup3*p1mup
     -   5*F(634))*P(7))
       Fa3(200) = (mup3mup5*F(1151))/2.+p1mup5*(2*p2mup3*F(800)+p3mup3
     -   *F(1153)+p4mup3*F(1155))+p2mup3*p4mup5*F(1156)+p3mup3*p4mup5
     -   *F(1157)+p3mup5*p4mup3*F(1158)+4*p2mup3*p2mup5*F(1159)+p1mup
     -   3*(p3mup5*F(1152)+p4mup5*F(1154)+2*(p1mup5*F(799)+p2mup5*F(1
     -   160)))-2*(p4mup3*p4mup5*F(813)-p2mup3*p3mup5*F(1161)-p2mup5*
     -   p3mup3*F(1162)-p2mup5*p4mup3*F(1163)-p3mup3*p3mup5*F(1164))
       Fa3(201) = (mup3mup5*F(1165))/2.+p1mup3*p2mup5*F(1166)+p1mup5*p
     -   2mup3*F(1167)+p1mup3*p4mup5*F(1168)+p1mup5*p4mup3*F(1169)+p2
     -   mup3*p4mup5*F(1170)+p2mup5*p4mup3*F(1171)+p2mup3*p2mup5*F(11
     -   72)+2*(p3mup3*p4mup5*F(822)+p3mup5*p4mup3*F(823)+p2mup3*p3mu
     -   p5*F(1173)+p2mup5*p3mup3*F(1174)+p4mup3*p4mup5*F(1175)+(p1mu
     -   p3*p3mup5*F(293)+p1mup5*p3mup3*F(310)+2*p3mup3*p3mup5*F(312)
     -   +p1mup3*p1mup5*F(828))*P(7))
       Fa3(202) = (mup3mup5*F(1176))/2.+p1mup5*p2mup3*F(1178)+p1mup5*p
     -   3mup3*F(1180)+p2mup3*p3mup5*F(1182)+p2mup5*(2*p4mup3*F(840)+
     -   p2mup3*F(1181)+p3mup3*F(1183))-2*(p1mup5*p4mup3*F(841)-p3mup
     -   3*p4mup5*F(842)+p3mup5*p4mup3*F(843)-p3mup3*p3mup5*F(1184)-p
     -   2mup3*p4mup5*F(1185))+p1mup3*(p2mup5*F(1177)+p3mup5*F(1179)-
     -   2*(p1mup5*F(832)-p4mup5*F(1186)))+4*p4mup3*p4mup5*F(1187)
       Fa3(203) = p1mup3*p3mup4*F(1189)+p1mup3*p4mup4*F(1192)+p1mup4*(
     -   -2*p2mup3*F(847)+p1mup3*F(1188)+p3mup3*F(1190)+p4mup3*F(1193
     -   ))+p3mup3*p3mup4*F(1195)-2*(p2mup4*(p4mup3*F(860)+p1mup3*F(1
     -   191)+2*p2mup3*F(1194))-p2mup3*p3mup4*F(1196))+p2mup4*p3mup3*
     -   F(1197)+p2mup3*p4mup4*F(1198)+p3mup3*p4mup4*F(1199)+p3mup4*p
     -   4mup3*F(1200)+p4mup3*p4mup4*F(1201)
       Fa3(204) = p1mup3*p4mup4*F(1206)+p1mup4*(2*p3mup3*F(873)+p2mup3
     -   *F(1203)+p1mup3*F(1204)+p4mup3*F(1207))+p2mup3*p3mup4*F(1209
     -   )+p2mup3*p4mup4*F(1210)+p2mup4*(2*p3mup3*F(874)+p1mup3*F(120
     -   2)+p2mup3*F(1208)+p4mup3*F(1211))+4*p3mup3*p3mup4*F(1212)+2*
     -   (p3mup3*p4mup4*F(875)+p3mup4*p4mup3*F(876)+p1mup3*p3mup4*F(1
     -   205)+p4mup3*p4mup4*F(1213))
       Fa3(205) = 8*(p1mup4*p1mup5*p2mup3*F(634)-p1mup3*p1mup4*p2mup5*
     -   F(882)+p2mup3*p2mup4*p2mup5*F(883))+mup4mup5*p4mup3*F(1215)+
     -   mup3mup5*p4mup4*F(1216)+mup3mup4*p4mup5*F(1217)+2*(mup3mup5*
     -   p3mup4*F(1214)+mup3mup4*p3mup5*F(273)*P(7))
       Fa3(206) = mup3mup4*p3mup5*F(885)-mup4mup5*p4mup3*F(886)+mup3mu
     -   p5*p4mup4*F(887)+mup3mup4*p4mup5*F(888)-2*(p1mup3*p1mup5*p2m
     -   up4*F(828)-p1mup4*(p1mup5*p2mup3*F(684)+p1mup3*p2mup5*F(889)
     -   ))-mup3mup5*p3mup4*F(1218)
       Fa3(207) = p1mup2*(p2mup4*p4mup5*F(892)+p4mup4*(2*(p4mup5*F(51)
     -   +p1mup5*F(755))+p2mup5*F(893))+2*(p2mup5*p3mup4*F(891)+p1mup
     -   4*(p4mup5*F(754)+p2mup5*F(890)+2*p3mup5*F(896)))+4*(p1mup5*p
     -   2mup4*F(750)+p1mup5*p3mup4*F(752)+p2mup4*p3mup5*F(897)+p3mup
     -   4*p3mup5*F(900)))+p3mup2*(p2mup4*p4mup5*F(903)+4*(p2mup4*p3m
     -   up5*F(102)+p1mup4*p3mup5*F(907)+p1mup4*p1mup5*F(910)+p1mup5*
     -   p2mup4*F(911))+p2mup5*p4mup4*F(912)+2*(p1mup4*p2mup5*F(670)+
     -   p4mup4*p4mup5*F(906)+p3mup5*p4mup4*F(908)+p3mup4*(p2mup5*F(6
     -   72)+p4mup5*F(904)+2*(p3mup5*F(100)+p1mup5*F(913)))))+p4mup2*
     -   (2*p1mup4*p4mup5*F(58)+p1mup5*p2mup4*F(895)+p2mup4*p3mup5*F(
     -   898)+p2mup4*p2mup5*F(901)+p2mup4*p4mup5*F(902)+4*(p1mup4*p1m
     -   up5*F(758)+p3mup4*p3mup5*F(905))+p4mup4*(2*(p4mup5*F(57)+p1m
     -   up5*F(60))+p2mup5*F(909))+p1mup4*p2mup5*F(1219)+p2mup5*p3mup
     -   4*F(1220))
       Fa3(208) = p1mup5*(-4*p2mup4*p3mup2*F(263)+8*p3mup4*(p1mup2*F(6
     -   34)+p3mup2*F(947)))+4*(-(p1mup2*p3mup5*(2*p1mup4*F(938)+p2mu
     -   p4*F(939)))+p3mup4*(p1mup2*p2mup5*F(940)+2*(p3mup5*p4mup2*F(
     -   316)+p1mup2*p3mup5*F(943))+p2mup5*p4mup2*F(948))+p3mup2*(p1m
     -   up4*p2mup5*F(346)-p2mup5*p4mup4*F(362)-p2mup4*(p4mup5*F(357)
     -   +p3mup5*F(941))+p3mup4*(p2mup5*F(942)+2*(p3mup5*F(306)+p4mup
     -   5*F(945)))+2*(p1mup4*p3mup5*F(944)+p3mup5*p4mup4*F(949))))+p
     -   2mup4*p3mup5*p4mup2*F(1221)
       Fa3(209) = p3mup3*F(1222)+p4mup3*F(1223)
       Fa3(210) = p2mup4*p3mup3*F(1226)+p2mup3*(2*(p2mup4*F(952)+p1mup
     -   4*F(956))+p3mup4*F(1225)+p4mup4*F(1227))+p2mup4*p4mup3*F(122
     -   8)+p3mup3*p3mup4*F(1230)+p3mup3*p4mup4*F(1233)+p3mup4*p4mup3
     -   *F(1234)+2*(p1mup3*p2mup4*F(1224)+p1mup4*(p3mup3*F(962)+p4mu
     -   p3*F(966)+2*p1mup3*F(1229))+p1mup3*p3mup4*F(1231)+p1mup3*p4m
     -   up4*F(1232)+p4mup3*p4mup4*F(1235))
       Fa3(211) = 4*p4mup3*p4mup4*F(983)+p1mup3*p2mup4*F(1236)+p1mup3*
     -   p1mup4*F(1238)+p1mup3*p3mup4*F(1239)+p2mup3*(2*p4mup4*F(977)
     -   +p1mup4*F(1237)+p2mup4*F(1241)+p3mup4*F(1242))+p3mup3*(2*p4m
     -   up4*F(981)+p1mup4*F(1240)+p2mup4*F(1243)+p3mup4*F(1244))-2*(
     -   p2mup4*p4mup3*F(978)-p1mup4*p4mup3*F(980)-p3mup4*p4mup3*F(98
     -   2)-p1mup3*p4mup4*F(1245))
       Fa3(212) = 2*mup4mup5*p4mup3*F(986)+8*(p1mup3*p1mup5*p2mup4*F(6
     -   34)-p1mup3*p1mup4*p2mup5*F(938)+p2mup3*p2mup4*p2mup5*F(989))
     -   +mup3mup5*p3mup4*F(1246)+mup3mup4*p3mup5*F(1247)+mup3mup4*p4
     -   mup5*F(1248)+mup3mup5*p4mup4*F(1249)
       Fa3(213) = mup3mup5*p3mup2*F(1252)+mup2mup5*p3mup3*F(1253)+mup2
     -   mup3*p3mup5*F(1254)+mup3mup5*p4mup2*F(1255)+mup2mup5*p4mup3*
     -   F(1256)+mup2mup3*p4mup5*F(1257)-2*(mup2mup5*p2mup3*F(1250)-m
     -   up3mup5*p1mup2*F(1258)-p1mup5*(-2*(p3mup2*p4mup3*F(263)-p4mu
     -   p2*(p3mup3*F(348)+2*(p2mup3*F(255)+p4mup3*F(364)+p1mup3*F(63
     -   4))))+mup2mup3*F(1259))-mup2mup3*p2mup5*F(1260))+p1mup3*(-2*
     -   mup2mup5*F(1251)+4*(p3mup5*p4mup2*F(344)-p3mup2*p4mup5*F(126
     -   3))+8*(p1mup2*p4mup5*F(633)+p2mup5*p4mup2*F(638)+p4mup2*p4mu
     -   p5*F(1270)))+4*(p3mup3*p3mup5*p4mup2*F(351)-p2mup5*p3mup2*p4
     -   mup3*F(379)+p1mup2*p3mup3*p4mup5*F(1261)-p1mup2*p3mup5*p4mup
     -   3*F(1262)-p2mup3*p3mup2*p4mup5*F(1265)+p2mup3*p3mup5*p4mup2*
     -   F(1266)+p2mup5*p3mup3*p4mup2*F(1267)-p3mup2*p3mup5*p4mup3*F(
     -   1268)-p3mup2*p4mup3*p4mup5*F(1271)+p3mup3*p4mup2*p4mup5*F(12
     -   72)+p3mup5*p4mup2*p4mup3*F(1273)+p3mup2*p3mup3*p4mup5*F(1274
     -   ))+8*(p1mup2*p2mup3*p4mup5*F(635)+p2mup3*p2mup5*p4mup2*F(126
     -   4)+p2mup3*p4mup2*p4mup5*F(1269)+p4mup2*p4mup3*p4mup5*F(1275)
     -   +p1mup2*p4mup3*p4mup5*F(1276)+p2mup5*p4mup2*p4mup3*F(1277))
       Fa3(214) = mup3mup5*p1mup2*F(1278)+mup3mup5*p4mup2*F(1281)+mup2
     -   mup5*p4mup3*F(1282)+mup2mup3*p4mup5*F(1283)+2*(mup2mup5*p2mu
     -   p3*F(1284)+mup2mup3*p2mup5*F(1285)+mup3mup5*p3mup2*F(1286)+m
     -   up2mup5*p3mup3*F(1287)+mup2mup3*p3mup5*F(1288))+4*(p1mup2*p2
     -   mup5*p4mup3*F(680)+p1mup2*p2mup3*p4mup5*F(1290)+p1mup2*p3mup
     -   5*p4mup3*F(1292)+p1mup2*p4mup3*p4mup5*F(1293)-p1mup2*p3mup3*
     -   p4mup5*F(1300))+p1mup3*(mup2mup5*F(1279)+4*(p1mup2*p4mup5*F(
     -   1289)-p4mup2*(p3mup5*F(293)+p1mup5*F(828)+p2mup5*F(1267)+p4m
     -   up5*F(1294))+p3mup2*p4mup5*F(1301)))+p1mup5*(mup2mup3*F(1280
     -   )+p4mup2*p4mup3*F(1296)-4*(p3mup3*p4mup2*F(310)-p1mup2*p4mup
     -   3*F(678)-p2mup3*p4mup2*F(1291)-p3mup2*p4mup3*F(1302)))-8*(p2
     -   mup3*p3mup5*p4mup2*F(311)-p3mup2*p3mup5*p4mup3*F(311)+p3mup3
     -   *p3mup5*p4mup2*F(312)-p2mup3*p4mup2*p4mup5*F(1295)-p3mup2*p4
     -   mup3*p4mup5*F(1297)+p3mup3*p4mup2*p4mup5*F(1298)+p2mup3*p3mu
     -   p2*p4mup5*F(1299)+p2mup3*p2mup5*p4mup2*F(1303)-p2mup5*p3mup2
     -   *p4mup3*F(1303)+p2mup5*p3mup3*p4mup2*F(1304)-p2mup5*p4mup2*p
     -   4mup3*F(1305)-p3mup2*p3mup3*p4mup5*F(1306)+p3mup5*p4mup2*p4m
     -   up3*F(1307))+16*p4mup2*p4mup3*p4mup5*F(1308)
       Fa3(215) = mup3mup5*p1mup2*F(1309)+mup2mup3*p1mup5*F(1311)+mup2
     -   mup5*p2mup3*F(1312)+mup3mup5*p3mup2*F(1313)+mup2mup5*p3mup3*
     -   F(1314)+mup2mup3*p3mup5*F(1315)-p1mup3*(-(mup2mup5*F(1310))+
     -   4*(4*p4mup2*p4mup5*F(153)+p1mup2*p3mup5*F(1320)+p2mup5*(4*p4
     -   mup2*F(63)+2*p1mup2*F(648)-p3mup2*F(1324))-p1mup5*p3mup2*F(1
     -   328)-p3mup2*p3mup5*F(1330)-p3mup2*p4mup5*F(1334)+p3mup5*p4mu
     -   p2*F(1338)+2*(p1mup2*p1mup5*F(1326)+p1mup2*p4mup5*F(1331)+p1
     -   mup5*p4mup2*F(1343))))-4*(p1mup2*p1mup5*p3mup3*F(1321)+p1mup
     -   2*p2mup3*p3mup5*F(1322)+p1mup5*p2mup3*p3mup2*F(1325)+p1mup2*
     -   p3mup3*p3mup5*F(1329)+p1mup2*p3mup3*p4mup5*F(1333)+p1mup2*p3
     -   mup5*p4mup3*F(1337)-p1mup5*p3mup2*p4mup3*F(1339)-p1mup5*p3mu
     -   p2*p3mup3*F(1340)+p1mup5*p3mup3*p4mup2*F(1345))-16*(p2mup3*p
     -   4mup2*p4mup5*F(157)+p4mup2*p4mup3*p4mup5*F(175)+p3mup2*p3mup
     -   3*p3mup5*F(1351))+2*(mup3mup5*p4mup2*F(1317)-mup2mup5*p4mup3
     -   *F(1318)-mup2mup3*p4mup5*F(1319)-p2mup5*(8*(p4mup2*p4mup3*F(
     -   171)+p2mup3*p4mup2*F(784))+mup2mup3*F(1316)+2*p1mup2*p3mup3*
     -   F(1323)+4*(p3mup3*p4mup2*F(115)+p1mup2*p2mup3*F(1327)+p2mup3
     -   *p3mup2*F(1347)+p1mup2*p4mup3*F(1349)-p3mup2*p4mup3*F(1350)+
     -   p3mup2*p3mup3*F(1352))))-8*(p3mup3*p4mup2*p4mup5*F(162)+p1mu
     -   p2*p1mup5*p2mup3*F(649)+p2mup3*p3mup2*p4mup5*F(665)+p3mup3*p
     -   3mup5*p4mup2*F(667)+p1mup2*p2mup3*p4mup5*F(1332)+p2mup3*p3mu
     -   p5*p4mup2*F(1335)-p3mup2*p3mup5*p4mup3*F(1336)+p1mup2*p4mup3
     -   *p4mup5*F(1341)-p3mup2*p4mup3*p4mup5*F(1342)+p1mup5*p2mup3*p
     -   4mup2*F(1344)+p2mup3*p3mup2*p3mup5*F(1346)+p1mup2*p1mup5*p4m
     -   up3*F(1348)+p3mup2*p3mup3*p4mup5*F(1353)+p1mup5*p4mup2*p4mup
     -   3*F(1354)+p3mup5*p4mup2*p4mup3*F(1355))
       Fa3(216) = 8*p2mup4*p2mup5*p3mup2*F(1264)+mup4mup5*p3mup2*F(135
     -   6)+mup2mup5*p3mup4*F(1357)+mup2mup4*p3mup5*F(1358)+mup4mup5*
     -   p4mup2*F(1361)+mup2mup5*p4mup4*F(1362)+mup2mup4*p4mup5*F(136
     -   3)+p1mup5*(4*(p3mup2*p4mup4*F(262)+p3mup4*p4mup2*F(263))+mup
     -   2mup4*F(1365))+mup2mup4*p2mup5*F(1366)+2*(mup2mup5*p2mup4*F(
     -   1359)+mup4mup5*p1mup2*F(1364)+p1mup4*(mup2mup5*F(1360)+2*(p3
     -   mup2*p4mup5*F(260)-p3mup5*p4mup2*F(1368))))+4*(p1mup2*p3mup4
     -   *p4mup5*F(1262)+p1mup2*p3mup5*p4mup4*F(1367)+p3mup4*p4mup2*p
     -   4mup5*F(1369)+p3mup5*p4mup2*p4mup4*F(1370))
       Fa3(217) = mup4mup5*p1mup2*F(1371)+mup2mup5*p1mup4*F(1372)+mup2
     -   mup5*p2mup4*F(1374)+mup2mup5*p3mup4*F(1375)+mup4mup5*p4mup2*
     -   F(1376)+mup2mup5*p4mup4*F(1377)+mup2mup4*p4mup5*F(1378)+2*(m
     -   up2mup4*p2mup5*F(1379)+mup4mup5*p3mup2*F(1380)+mup2mup4*p3mu
     -   p5*F(1381))+p1mup5*(-8*p1mup2*p1mup4*F(1326)+mup2mup4*F(1373
     -   )-4*(p3mup2*p4mup4*F(912)+p3mup4*p4mup2*F(1385)))-8*(p3mup5*
     -   p4mup2*p4mup4*F(936)+p2mup4*p2mup5*(2*p3mup2*F(771)+p1mup2*F
     -   (1327))+p3mup4*p4mup2*p4mup5*F(1386))-4*(p1mup2*p3mup4*p4mup
     -   5*F(1382)+p1mup2*p3mup5*p4mup4*F(1383)-p1mup4*p3mup5*p4mup2*
     -   F(1384)+p1mup4*p3mup2*p4mup5*F(1387))
       Fa3(218) = mup4mup5*p1mup2*F(1388)+p1mup5*(-4*(p3mup2*p4mup4*F(
     -   313)+p3mup4*p4mup2*F(948))+mup2mup4*F(1390))+mup4mup5*p3mup2
     -   *F(1391)+mup2mup5*p3mup4*F(1392)+mup2mup4*p3mup5*F(1393)-2*(
     -   mup2mup5*p2mup4*F(1394)-mup2mup4*p2mup5*F(1395)-mup4mup5*p4m
     -   up2*F(1396)-mup2mup5*p4mup4*F(1397)-mup2mup4*p4mup5*F(1398))
     -   -4*(p1mup2*p3mup4*p4mup5*F(1399)+p1mup2*p3mup5*p4mup4*F(1400
     -   ))+p1mup4*(mup2mup5*F(1389)-4*(p3mup2*p4mup5*F(298)-p3mup5*p
     -   4mup2*F(1401)))-8*(p3mup4*p4mup2*p4mup5*F(314)+p2mup4*p2mup5
     -   *p3mup2*F(692)-p3mup5*p4mup2*p4mup4*F(1402))
       Fa3(219) = p1mup2*(2*p1mup4*F(1403)+p2mup4*F(1404)+p3mup4*F(140
     -   5)+p4mup4*F(1406))+p1mup4*p4mup2*F(1407)+p2mup4*p4mup2*F(140
     -   8)+p3mup2*p3mup4*F(1409)+p3mup2*p4mup4*F(1412)+2*(p1mup4*p3m
     -   up2*F(1410)+p2mup4*p3mup2*F(1411)+p3mup4*p4mup2*F(1413))+p4m
     -   up2*p4mup4*F(1414)
       Fa3(220) = mup2mup5*p2mup4*F(1415)+mup2mup4*p2mup5*F(1416)+mup4
     -   mup5*p4mup2*F(1420)+mup2mup5*(2*(p1mup4*F(1418)+p3mup4*F(141
     -   9))+p4mup4*F(1421))+mup2mup4*p4mup5*F(1422)-p1mup2*(16*p1mup
     -   4*p1mup5*F(12)-mup4mup5*F(1417)+8*(p3mup5*p4mup4*F(916)+p2mu
     -   p4*p2mup5*F(1425)+p3mup4*p4mup5*F(1426)))+2*(mup4mup5*p3mup2
     -   *F(1423)+mup2mup4*p1mup5*F(1424)+mup2mup4*p3mup5*F(1427))-8*
     -   (p2mup4*p2mup5*p3mup2*F(671)+p3mup5*p4mup2*p4mup4*F(1429)+p1
     -   mup4*p3mup2*p4mup5*F(1430)+2*p1mup4*p3mup5*p4mup2*F(1431)+p1
     -   mup5*p3mup2*p4mup4*F(1432)+p3mup4*p4mup2*(p4mup5*F(1428)+2*p
     -   1mup5*F(1433)))
       Fa3(221) = mup2mup5*p1mup4*F(1435)+mup2mup4*p1mup5*F(1436)+mup2
     -   mup5*p2mup4*F(1437)+mup2mup4*p2mup5*F(1438)+mup4mup5*p3mup2*
     -   F(1439)+mup4mup5*p4mup2*F(1442)+2*(mup2mup5*p3mup4*F(1440)-m
     -   up2mup5*p4mup4*F(1443)-mup2mup4*p4mup5*F(1444)+p3mup5*(4*p1m
     -   up4*p4mup2*F(1429)-mup2mup4*F(1441)+8*(p4mup2*p4mup4*F(120)+
     -   p1mup2*p4mup4*F(1446))))+p1mup2*(8*p1mup4*p1mup5*F(136)+mup4
     -   mup5*F(1434)+4*p2mup4*p2mup5*F(1445)+16*p3mup4*p4mup5*F(1448
     -   ))+8*(p2mup4*p2mup5*p3mup2*F(665)+p1mup4*p3mup2*p4mup5*F(142
     -   9)-p1mup5*p3mup4*p4mup2*F(1447)+p1mup5*p3mup2*p4mup4*F(1449)
     -   )+16*p3mup4*p4mup2*p4mup5*F(1450)
       Fa3(222) = p1mup4*p3mup2*F(1452)+p1mup2*(p3mup4*F(1451)+2*(p1mu
     -   p4*F(1453)+p2mup4*F(1454))+p4mup4*F(1455))+p3mup2*(p3mup4*F(
     -   1457)+p2mup4*F(1458)+p4mup4*F(1459))+p4mup2*(p1mup4*F(1456)+
     -   p3mup4*F(1460)-2*(p2mup4*F(1461)-p4mup4*F(1462)))
       Fa3(223) = -8*(p3mup4*p4mup2*(p4mup5*F(361)+p1mup5*F(364))-p3mu
     -   p5*p4mup2*p4mup4*F(365))+p1mup4*(-8*p3mup5*p4mup2*F(359)+2*m
     -   up2mup5*F(630))+mup2mup4*p1mup5*F(631)+mup2mup4*p4mup5*F(124
     -   8)+mup2mup5*p4mup4*F(1249)+mup2mup5*p2mup4*F(1463)+mup2mup4*
     -   p2mup5*F(1464)+p1mup2*(-8*(p3mup4*p4mup5*F(354)+p3mup5*p4mup
     -   4*F(359))+2*mup4mup5*F(1465))+mup4mup5*p3mup2*F(1466)+mup2mu
     -   p5*p3mup4*F(1467)+mup2mup4*p3mup5*F(1468)+2*mup4mup5*p4mup2*
     -   F(1469)+4*p2mup4*p2mup5*p3mup2*F(1470)
       Fa3(224) = p2mup3*p3mup2*F(1472)-2*(p1mup2*(2*p1mup3*F(1471)+p2
     -   mup3*F(1473)+p3mup3*F(1474))-p1mup3*p3mup2*F(1475)-p3mup2*p3
     -   mup3*F(1476)-p4mup2*p4mup3*F(1477))
       Fa3(225) = p1mup2*(-2*(p1mup3*F(1478)+p2mup3*F(1479))+p3mup3*F(
     -   1480))+p1mup3*p3mup2*F(1481)+p4mup2*p4mup3*F(1483)+2*(p2mup3
     -   *p3mup2*F(1482)+p3mup2*p3mup3*F(1484))
       Fa3(226) = p1mup2*(p2mup3*F(1486)-2*(p1mup3*F(1485)+p3mup3*F(14
     -   88)))+p4mup2*p4mup3*F(1490)+2*p3mup2*(p2mup3*F(1487)+p1mup3*
     -   F(1489)+2*p3mup3*F(1491))
       Fa3(227) = mup2mup5*p2mup3*F(1492)+mup2mup3*p2mup5*F(1493)+mup3
     -   mup5*p4mup2*F(1495)+mup2mup5*p4mup3*F(1496)+mup2mup3*p4mup5*
     -   F(1497)-4*(p1mup3*p2mup5*p4mup2*F(258)-p1mup5*p2mup3*p4mup2*
     -   F(263)-p2mup3*p3mup5*p4mup2*F(273)-p2mup3*p4mup2*p4mup5*F(27
     -   5)-p2mup5*p3mup3*p4mup2*F(282)-p1mup2*p2mup5*p4mup3*F(940)-p
     -   2mup5*p3mup2*p4mup3*F(1302)-p1mup2*p2mup3*p4mup5*F(1499)-p2m
     -   up3*p2mup5*p4mup2*F(1500)-p2mup5*p4mup2*p4mup3*F(1502)-p2mup
     -   3*p3mup2*p4mup5*F(1505))+8*(p1mup2*p3mup5*p4mup3*F(266)-p1mu
     -   p2*p3mup3*p4mup5*F(279)+p3mup5*p4mup2*p4mup3*F(945)+p1mup5*p
     -   4mup3*(p4mup2*F(364)+p1mup2*F(634)+p3mup2*F(947))+p3mup2*p3m
     -   up5*p4mup3*F(1025)+p1mup2*p4mup3*p4mup5*F(1501)-p1mup3*p4mup
     -   5*(p3mup2*F(279)-p1mup2*F(882)+p4mup2*F(1503))+p3mup2*p4mup3
     -   *p4mup5*F(1504)+p4mup2*p4mup3*p4mup5*F(1506)+p3mup2*p3mup3*p
     -   4mup5*F(1507)-p3mup3*p4mup2*p4mup5*F(1508))+2*(mup3mup5*p1mu
     -   p2*F(1494)+mup3mup5*p3mup2*F(1498)+(-(mup2mup3*(p1mup5*F(263
     -   )+p3mup5*F(273)))+mup2mup5*(p1mup3*F(258)+p3mup3*F(600)))*P(
     -   7))
       Fa3(228) = p1mup2*(2*p1mup3*F(1509)+p2mup3*F(1510)+p3mup3*F(151
     -   1))+p1mup3*p3mup2*F(1512)+p2mup3*p3mup2*F(1513)+p3mup2*p3mup
     -   3*F(1514)+4*p4mup2*p4mup3*F(1515)
       Fa3(229) = mup2mup5*p2mup3*F(1516)+mup2mup3*p2mup5*F(1517)+mup3
     -   mup5*p3mup2*F(1520)+mup2mup5*p3mup3*F(1521)+mup2mup3*p3mup5*
     -   F(1522)-4*(p2mup5*p3mup3*p4mup2*F(1345)-p1mup3*p2mup5*p3mup2
     -   *F(1529)+p2mup3*p3mup2*p3mup5*F(1531)+p2mup3*p2mup5*p3mup2*F
     -   (1532)+p2mup3*p3mup2*p4mup5*F(1535)+p2mup3*p3mup5*p4mup2*F(1
     -   537)+p2mup5*p3mup2*p4mup3*F(1538)+p2mup5*p3mup2*p3mup3*F(154
     -   3))-8*(p1mup3*p3mup2*p3mup5*F(770)+p1mup3*p2mup5*p4mup2*F(13
     -   43)+p2mup3*p2mup5*p4mup2*F(1344)+p2mup5*p4mup2*p4mup3*F(1354
     -   )+p3mup2*p3mup3*p4mup5*F(1539)+p3mup3*p3mup5*p4mup2*F(1540)+
     -   p1mup3*p3mup5*p4mup2*F(1541)+p3mup2*p3mup3*p3mup5*F(1542)+p3
     -   mup3*p4mup2*p4mup5*F(1545)+p3mup2*p3mup5*p4mup3*F(1548)+p3mu
     -   p5*p4mup2*p4mup3*F(1550))+2*(mup2mup5*p1mup3*F(1519)+mup3mup
     -   5*p4mup2*F(1524)+mup2mup5*p4mup3*F(1525)+mup2mup3*p4mup5*F(1
     -   526)-p1mup2*(-(mup3mup5*F(1518))+2*(p2mup3*p3mup5*F(1527)+p3
     -   mup3*(p2mup5*F(1528)+2*(p3mup5*F(47)+p1mup5*F(753)+p4mup5*F(
     -   1536))))+4*(p2mup5*p4mup3*F(1343)+p2mup3*p2mup5*F(1425)+p1mu
     -   p3*(p3mup5*F(751)+p2mup5*F(890)+2*(p1mup5*F(12)+p4mup5*F(153
     -   3)))+p3mup5*p4mup3*F(1541))+8*(p1mup5*p2mup3*F(750)+p1mup5*p
     -   4mup3*F(758)+p2mup3*p4mup5*F(1534)+p4mup3*p4mup5*F(1544)))-p
     -   1mup5*(-(mup2mup3*F(1523))+2*p2mup3*p3mup2*F(1530)+4*(p3mup2
     -   *p3mup3*F(1546)+p3mup3*p4mup2*F(1547))+8*(p2mup3*p4mup2*F(61
     -   )+p1mup3*p4mup2*F(758)+p1mup3*p3mup2*F(910)+p4mup2*p4mup3*F(
     -   1549)+p3mup2*p4mup3*F(1552))))-16*(p2mup3*p4mup2*p4mup5*F(15
     -   3)+p1mup3*p3mup2*p4mup5*F(1431)+p1mup3*p4mup2*p4mup5*F(1544)
     -   +p3mup2*p4mup3*p4mup5*F(1551)+p4mup2*p4mup3*p4mup5*F(1553))
       Fa3(230) = mup2mup3*p1mup5*F(1556)+mup2mup5*p2mup3*F(1557)+mup2
     -   mup3*p2mup5*F(1558)+mup3mup5*p3mup2*F(1559)+mup3mup5*p4mup2*
     -   F(1562)+2*(mup2mup5*p3mup3*F(1560)+mup2mup3*p3mup5*F(1561)+m
     -   up2mup5*p4mup3*F(1563)+mup2mup3*p4mup5*F(1564))-p1mup2*(-(mu
     -   p3mup5*F(1554))+4*(4*(p3mup5*p4mup3*F(905)+p4mup3*p4mup5*F(1
     -   448))+2*p1mup5*p4mup3*F(1541)+p1mup3*(2*(p1mup5*F(79)+p3mup5
     -   *F(766)+p4mup5*F(1536))+p2mup5*F(1565))-p1mup5*p2mup3*F(1566
     -   )+p2mup3*p3mup5*F(1567)+p2mup3*p2mup5*F(1568)+p2mup3*p4mup5*
     -   F(1569)+p2mup5*p4mup3*F(1570)+p3mup3*(2*p1mup5*F(770)+4*(p3m
     -   up5*F(100)+p4mup5*F(905))+p2mup5*F(1575))))-4*(p1mup5*p2mup3
     -   *p4mup2*F(1572)+p1mup5*p2mup3*p3mup2*F(1577))-p1mup3*(-(mup2
     -   mup5*F(1555))+4*(p2mup5*p4mup2*F(1571)+p2mup5*p3mup2*F(1576)
     -   )+8*(p3mup2*p4mup5*F(1539)+p4mup2*p4mup5*F(1545)+p1mup5*(p4m
     -   up2*F(922)+p3mup2*F(1574))+p3mup5*(p4mup2*F(1539)+p3mup2*F(1
     -   578))))-8*(p2mup3*(p3mup5*p4mup2*F(98)+p3mup2*p4mup5*F(98)+p
     -   4mup2*p4mup5*F(162)+p2mup5*p4mup2*F(665))+p2mup5*p3mup3*p4mu
     -   p2*F(1353)+p2mup5*p3mup2*p4mup3*F(1353)+p2mup5*p4mup2*p4mup3
     -   *F(1355)+p1mup5*p3mup3*p4mup2*F(1548)+p1mup5*p3mup2*p4mup3*F
     -   (1548)+p1mup5*p4mup2*p4mup3*F(1550)+p2mup3*p3mup2*p3mup5*F(1
     -   573)+p2mup3*p2mup5*p3mup2*F(1579)+p1mup5*p3mup2*p3mup3*F(158
     -   0)+p2mup5*p3mup2*p3mup3*F(1581))-16*(p3mup3*p3mup5*p4mup2*F(
     -   933)+p3mup2*p3mup5*p4mup3*F(933)+p3mup2*p3mup3*p4mup5*F(933)
     -   +p3mup5*p4mup2*p4mup3*F(937)+p3mup3*p4mup2*p4mup5*F(937)+p3m
     -   up2*p4mup3*p4mup5*F(937)+p4mup2*p4mup3*p4mup5*F(1450)+p3mup2
     -   *p3mup3*p3mup5*F(1582))
       Fa3(231) = p2mup4*p3mup2*F(1584)+p1mup2*(p3mup4*F(1583)-4*p1mup
     -   4*F(1585)+p2mup4*F(1586)+p4mup4*F(1587))+p3mup2*(2*(p3mup4*F
     -   (1589)+p1mup4*F(1590))+p4mup4*F(1591))+p4mup2*(p2mup4*F(1588
     -   )+p3mup4*F(1592)+2*(p1mup4*F(1593)+p4mup4*F(1594)))
       Fa3(232) = p1mup4*p3mup2*F(1598)+p2mup4*p3mup2*F(1599)+p1mup2*(
     -   p2mup4*F(1596)+p3mup4*F(1597)+2*(p1mup4*F(1595)+p4mup4*F(160
     -   1)))+2*(p3mup2*p3mup4*F(1600)+p1mup4*p4mup2*F(1602)+p2mup4*p
     -   4mup2*F(1603)+p3mup2*p4mup4*F(1604))+p3mup4*p4mup2*F(1605)+4
     -   *p4mup2*p4mup4*F(1606)
       Fa3(233) = -(mup2mup3*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F
     -   (383)-p3mup5*F(385)))+mup2mup5*(2*(p1mup3*F(376)+p2mup3*F(37
     -   7))+p3mup3*F(378)+p4mup3*F(1012))-mup3mup5*(p3mup2*F(1607)-2
     -   *p1mup2*F(1608)+p4mup2*F(1609))
       Fa3(234) = mup3mup5*p1mup2*F(380)-mup2mup5*p1mup3*F(380)+mup2mu
     -   p3*p1mup5*F(380)+2*(mup3mup5*p4mup2*F(356)-mup2mup5*p4mup3*F
     -   (356)+mup2mup3*p4mup5*F(356)-mup2mup5*p2mup3*F(381)+mup2mup3
     -   *p2mup5*F(381)+(mup3mup5*p3mup2-mup2mup5*p3mup3+mup2mup3*p3m
     -   up5)*F(382))
       Fa3(235) = mup3mup5*(p1mup2*F(383)+2*(p3mup2*F(4465)+p4mup2*F(4
     -   467)))-mup2mup5*(p1mup3*F(383)+2*(p3mup3*F(4465)+p4mup3*F(44
     -   67)+p2mup3*F(4469)))+mup2mup3*(p1mup5*F(383)+2*(p3mup5*F(446
     -   5)+p4mup5*F(4467)+p2mup5*F(4469)))
       Fa3(236) = -(mup2mup5*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F
     -   (383)-p3mup4*F(385)))+mup2mup4*(2*(p1mup5*F(263)+p2mup5*F(37
     -   9))+p4mup5*F(383)+p3mup5*F(386))+mup4mup5*(2*p1mup2*F(1262)+
     -   p3mup2*F(1607)+p4mup2*F(1609))
       Fa3(237) = mup4mup5*p1mup2*F(380)-mup2mup5*p1mup4*F(380)+mup2mu
     -   p4*p1mup5*F(380)+2*(mup4mup5*p4mup2*F(356)-mup2mup5*p4mup4*F
     -   (356)+mup2mup4*p4mup5*F(356)-mup2mup5*p2mup4*F(381)+mup2mup4
     -   *p2mup5*F(381)+(mup4mup5*p3mup2-mup2mup5*p3mup4+mup2mup4*p3m
     -   up5)*F(382))
       Fa3(238) = mup4mup5*(p1mup2*F(383)+2*(p3mup2*F(4465)+p4mup2*F(4
     -   467)))-mup2mup5*(p1mup4*F(383)+2*(p3mup4*F(4465)+p4mup4*F(44
     -   67)+p2mup4*F(4469)))+mup2mup4*(p1mup5*F(383)+2*(p3mup5*F(446
     -   5)+p4mup5*F(4467)+p2mup5*F(4469)))
       Fa3(239) = p1mup2*(4*(p3mup4*F(270)+p4mup4*F(943))+2*p1mup4*F(1
     -   005)+p2mup4*F(1006))+2*(2*(p3mup2*p4mup4*F(306)+p3mup4*(p4mu
     -   p2*F(306)+p3mup2*F(312))+p4mup2*p4mup4*F(316))+p1mup4*p3mup2
     -   *F(351)+p1mup4*p4mup2*F(685)+p2mup4*p3mup2*F(1007)+p2mup4*p4
     -   mup2*F(1008))
       Fa3(240) = (mup2mup5*p2mup4-mup2mup4*p2mup5)*F(396)+2*mup4mup5*
     -   (p4mup2*F(356)+p3mup2*F(382)+p1mup2*F(1610))
       Fa3(241) = p1mup4*p3mup2*F(998)+p1mup4*p4mup2*F(1002)-p1mup2*(2
     -   *(p1mup4*F(678)+p2mup4*F(680))+p4mup4*F(1001)-p3mup4*F(1034)
     -   )+2*(p3mup4*p4mup2*F(1003)+p3mup2*p4mup4*F(1003)-p2mup4*p4mu
     -   p2*F(1004)+p2mup4*p3mup2*F(1037))+4*(p3mup2*p3mup4*F(1039)-p
     -   4mup2*p4mup4*F(4472))
       Fa3(242) = (mup2mup5*p2mup4-mup2mup4*p2mup5)*F(4464)+2*mup4mup5
     -   *(p1mup2*F(357)+p3mup2*F(4465)+p4mup2*F(4467))
       Fa3(243) = 4*p1mup3*(p3mup2*F(266)+p1mup2*F(1611))+p2mup3*p3mup
     -   2*F(1612)+2*(p1mup2*p2mup3*F(1613)+p1mup2*p3mup3*F(1614)+p3m
     -   up2*p3mup3*F(1615)+p4mup2*p4mup3*F(1616))
       Fa3(244) = 4*p4mup2*p4mup3*F(314)+p1mup2*(2*(p1mup3*F(1617)+p2m
     -   up3*F(1618))+p3mup3*F(1619))+p3mup2*(4*p3mup3*F(311)+p1mup3*
     -   F(1620)+2*p2mup3*F(1621))
       Fa3(245) = 4*p4mup2*p4mup3*F(316)+p1mup2*(4*p3mup3*F(270)+2*p1m
     -   up3*F(1005)+p2mup3*F(1006))+2*p3mup2*(2*p3mup3*F(312)+p1mup3
     -   *F(351)+p2mup3*F(1007))
       Fa3(246) = (mup2mup5*p2mup3-mup2mup3*p2mup5)*F(396)+2*mup3mup5*
     -   (p4mup2*F(356)+p3mup2*F(382)+p1mup2*F(1610))
       Fa3(247) = p1mup2*(2*p1mup3*F(684)+4*p3mup3*F(945)+p2mup3*F(100
     -   9))+2*p3mup2*(p1mup3*F(995)+p2mup3*F(1010)+2*p3mup3*F(4471))
     -   +4*p4mup2*p4mup3*F(4476)
       Fa3(248) = (mup2mup5*p2mup3-mup2mup3*p2mup5)*F(4464)+2*mup3mup5
     -   *(p1mup2*F(357)+p3mup2*F(4465)+p4mup2*F(4467))
       Fa3(249) = 4*p1mup4*(p4mup2*F(364)+p1mup2*F(634)+p3mup2*F(947))
     -   +p2mup4*p3mup2*F(990)+p2mup4*p4mup2*F(993)+2*(p1mup2*(p4mup4
     -   *F(684)+p2mup4*F(940)-p3mup4*F(991))+p3mup2*p4mup4*F(995)+p4
     -   mup2*p4mup4*F(996)-p3mup2*p3mup4*F(1029)-p3mup4*p4mup2*F(103
     -   1))
       Fa3(250) = p1mup2*(2*p1mup4*F(684)+4*(p4mup4*F(361)+p3mup4*F(94
     -   5))+p2mup4*F(1009))+2*(p1mup4*p3mup2*F(995)+p1mup4*p4mup2*F(
     -   996)+p2mup4*p3mup2*F(1010)+p2mup4*p4mup2*F(1011)+2*(p3mup2*p
     -   3mup4*F(4471)+p3mup4*p4mup2*F(4473)+p3mup2*p4mup4*F(4473)+p4
     -   mup2*p4mup4*F(4476)))
       Fa3(251) = mup2mup5*(2*(p1mup3*F(263)+p2mup3*F(379))+p4mup3*F(3
     -   83)+p3mup3*F(418))-mup2mup3*(2*(p1mup5*F(263)+p2mup5*F(379))
     -   +p4mup5*F(383)+p3mup5*F(418))-mup3mup5*(-2*p1mup2*F(1608)+p4
     -   mup2*F(1609)+p3mup2*F(1622))
       Fa3(252) = mup3mup5*(p1mup2*F(418)+2*(p4mup2*F(4465)+p3mup2*F(4
     -   466)))-mup2mup5*(p1mup3*F(418)+2*(p4mup3*F(4465)+p3mup3*F(44
     -   66)+p2mup3*F(4468)))+mup2mup3*(p1mup5*F(418)+2*(p4mup5*F(446
     -   5)+p3mup5*F(4466)+p2mup5*F(4468)))
       Fa3(253) = -(mup2mup5*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F
     -   (383)+p3mup4*F(418)))+mup2mup4*(2*(p1mup5*F(263)+p2mup5*F(37
     -   9))+p4mup5*F(383)+p3mup5*F(418))+mup4mup5*(2*p1mup2*F(1262)+
     -   p4mup2*F(1609)+p3mup2*F(1622))
       Fa3(254) = mup4mup5*(p1mup2*F(418)+2*(p4mup2*F(4465)+p3mup2*F(4
     -   466)))-mup2mup5*(p1mup4*F(418)+2*(p4mup4*F(4465)+p3mup4*F(44
     -   66)+p2mup4*F(4468)))+mup2mup4*(p1mup5*F(418)+2*(p4mup5*F(446
     -   5)+p3mup5*F(4466)+p2mup5*F(4468)))
       Fa3(255) = p1mup2*(2*p1mup4*F(828)+p2mup4*F(1023)+4*(p4mup4*F(9
     -   45)+p3mup4*F(1025)))+2*(p1mup4*p4mup2*F(995)+p2mup4*p4mup2*F
     -   (1010)+p1mup4*p3mup2*F(1019)+p2mup4*p3mup2*F(1024)+2*(p3mup4
     -   *p4mup2*F(4471)+p3mup2*p4mup4*F(4471)+p4mup2*p4mup4*F(4473)+
     -   p3mup2*p3mup4*F(4475)))
       Fa3(256) = (mup2mup5*p2mup4-mup2mup4*p2mup5)*F(4463)+2*mup4mup5
     -   *(p1mup2*F(273)+p4mup2*F(4465)+p3mup2*F(4466))
       Fa3(257) = p1mup2*(2*(p1mup4*F(678)+p2mup4*F(680))+p4mup4*F(100
     -   1)+p3mup4*F(1020))-p1mup4*(p4mup2*F(1002)+p3mup2*F(1021))+2*
     -   (p2mup4*p4mup2*F(1004)+p2mup4*p3mup2*F(1022)+2*(p3mup2*p3mup
     -   4*F(4470)+p4mup2*p4mup4*F(4472)+(p3mup4*p4mup2+p3mup2*p4mup4
     -   )*F(4474)))
       Fa3(258) = 4*p1mup3*p3mup2*F(947)+p2mup3*p3mup2*F(1018)+2*(p1mu
     -   p2*(2*p1mup3*F(634)+p3mup3*F(828)+p2mup3*F(940))+p4mup2*p4mu
     -   p3*F(996)+p3mup2*p3mup3*F(1019))
       Fa3(259) = p1mup2*(2*(p1mup3*F(678)+p2mup3*F(680))+p3mup3*F(102
     -   0))-p1mup3*p3mup2*F(1021)+2*p2mup3*p3mup2*F(1022)+4*(p3mup2*
     -   p3mup3*F(4470)+p4mup2*p4mup3*F(4472))
       Fa3(260) = p1mup2*(2*p1mup3*F(828)+p2mup3*F(1023)+4*p3mup3*F(10
     -   25))+4*p4mup2*p4mup3*F(4473)+2*p3mup2*(p1mup3*F(1019)+p2mup3
     -   *F(1024)+2*p3mup3*F(4475))
       Fa3(261) = (mup2mup5*p2mup3-mup2mup3*p2mup5)*F(4463)+2*mup3mup5
     -   *(p1mup2*F(273)+p4mup2*F(4465)+p3mup2*F(4466))
       Fa3(262) = 4*p1mup4*(p4mup2*F(364)+p1mup2*F(634)+p3mup2*F(947))
     -   +p2mup4*p4mup2*F(993)+p2mup4*p3mup2*F(1018)+2*(p1mup2*(p4mup
     -   4*F(684)+p3mup4*F(828)+p2mup4*F(940))+p3mup4*p4mup2*F(995)+p
     -   3mup2*p4mup4*F(995)+p4mup2*p4mup4*F(996)+p3mup2*p3mup4*F(101
     -   9))
       Fa3(263) = mup2mup3*(2*(p1mup5*F(263)+p2mup5*F(379))+p3mup5*F(3
     -   86)+p4mup5*F(447))-mup2mup5*(2*(p1mup3*F(376)+p2mup3*F(377))
     -   +p3mup3*F(378)+p4mup3*F(448))+mup3mup5*(2*p1mup2*F(1262)+p3m
     -   up2*F(1607)+p4mup2*F(1623))
       Fa3(264) = mup3mup5*p1mup2*F(443)-mup2mup5*p1mup3*F(443)+mup2mu
     -   p3*p1mup5*F(443)-2*(mup2mup5*p2mup3*F(444)-mup2mup3*p2mup5*F
     -   (444)-mup3mup5*p3mup2*F(445)+mup2mup5*p3mup3*F(445)-mup2mup3
     -   *p3mup5*F(445)-(mup3mup5*p4mup2-mup2mup5*p4mup3+mup2mup3*p4m
     -   up5)*F(459))
       Fa3(265) = mup2mup5*(2*(p1mup4*F(376)+p2mup4*F(377))+p3mup4*F(3
     -   78)+p4mup4*F(448))+mup2mup4*(-2*(p1mup5*F(263)+p2mup5*F(379)
     -   )+p3mup5*F(385)+p4mup5*F(485))-mup4mup5*(p3mup2*F(1607)-2*p1
     -   mup2*F(1608)+p4mup2*F(1623))
       Fa3(266) = mup4mup5*p1mup2*F(443)-mup2mup5*p1mup4*F(443)+mup2mu
     -   p4*p1mup5*F(443)-2*(mup2mup5*p2mup4*F(444)-mup2mup4*p2mup5*F
     -   (444)-mup4mup5*p3mup2*F(445)+mup2mup5*p3mup4*F(445)-mup2mup4
     -   *p3mup5*F(445)-(mup4mup5*p4mup2-mup2mup5*p4mup4+mup2mup4*p4m
     -   up5)*F(459))
       Fa3(267) = p1mup2*(2*(p1mup4*F(1617)+p2mup4*F(1618))+p3mup4*F(1
     -   619)+p4mup4*F(1624))+p3mup2*(4*p3mup4*F(311)+p1mup4*F(1620)+
     -   2*(p2mup4*F(1621)+p4mup4*F(1627)))+p4mup2*(4*p4mup4*F(314)+p
     -   1mup4*F(1625)+2*(p2mup4*F(1626)+p3mup4*F(1628)))
       Fa3(268) = (mup2mup5*p2mup4-mup2mup4*p2mup5)*F(477)+2*mup4mup5*
     -   (p3mup2*F(445)+p4mup2*F(459)+p1mup2*F(1629))
       Fa3(269) = 4*p1mup3*(p3mup2*F(266)+p1mup2*F(1611))+p2mup3*p3mup
     -   2*F(1612)+2*(p1mup2*p2mup3*F(1613)+p1mup2*p3mup3*F(1614)+p3m
     -   up2*p3mup3*F(1615)+p4mup2*p4mup3*F(1630))
       Fa3(270) = 4*p4mup2*p4mup3*F(314)+p1mup2*(2*(p1mup3*F(1617)+p2m
     -   up3*F(1618))+p3mup3*F(1619))+p3mup2*(4*p3mup3*F(311)-p1mup3*
     -   F(1631)-2*p2mup3*F(1632))
       Fa3(271) = 4*p4mup2*p4mup3*F(315)+p1mup2*(2*p1mup3*F(1042)+p2mu
     -   p3*F(1043)+4*p3mup3*F(1045))+2*p3mup2*(p2mup3*F(1044)+p1mup3
     -   *F(1046)+2*p3mup3*F(1049))
       Fa3(272) = (mup2mup5*p2mup3-mup2mup3*p2mup5)*F(477)+2*mup3mup5*
     -   (p3mup2*F(445)+p4mup2*F(459)+p1mup2*F(1629))
       Fa3(273) = 4*p1mup4*(p3mup2*F(266)+p4mup2*F(354)+p1mup2*F(1611)
     -   )+p2mup4*(p3mup2*F(1612)+2*p1mup2*F(1613)+p4mup2*F(1634))+2*
     -   (p1mup2*p3mup4*F(1614)+p3mup2*p3mup4*F(1615)+p4mup2*p4mup4*F
     -   (1630)+p1mup2*p4mup4*F(1633)+p3mup2*p4mup4*F(1635)+p3mup4*p4
     -   mup2*F(1636))
       Fa3(274) = p1mup2*p2mup4*F(1043)+2*(p2mup4*p3mup2*F(1044)+p1mup
     -   4*(p4mup2*F(276)+p1mup2*F(1042)+p3mup2*F(1046))+p2mup4*p4mup
     -   2*F(1047))+4*(p3mup2*p4mup4*F(308)+p4mup2*p4mup4*F(315)+p1mu
     -   p2*p4mup4*F(1048)+p3mup4*(p4mup2*F(308)+p1mup2*F(1045)+p3mup
     -   2*F(1049)))
       Fa3(275) = mup2mup5*(2*(p1mup3*F(376)+p2mup3*F(377))+p3mup3*F(3
     -   78)+p4mup3*F(448))+mup2mup3*(-2*(p1mup5*F(263)+p2mup5*F(379)
     -   )+p3mup5*F(385)+p4mup5*F(485))-mup3mup5*(p3mup2*F(1607)-2*p1
     -   mup2*F(1608)+p4mup2*F(1623))
       Fa3(276) = mup2mup4*(2*(p1mup5*F(263)+p2mup5*F(379))+p3mup5*F(3
     -   86)+p4mup5*F(447))-mup2mup5*(2*(p1mup4*F(376)+p2mup4*F(377))
     -   +p3mup4*F(378)+p4mup4*F(448))+mup4mup5*(2*p1mup2*F(1262)+p3m
     -   up2*F(1607)+p4mup2*F(1623))
       Fa3(277) = p1mup2*(2*(p1mup4*F(1617)+p2mup4*F(1618))+p3mup4*F(1
     -   619)+p4mup4*F(1624))+p3mup2*(4*p3mup4*F(311)-p1mup4*F(1631)+
     -   2*(p4mup4*F(1627)-p2mup4*F(1632)))+p4mup2*(4*p4mup4*F(314)-p
     -   1mup4*F(1637)-2*(p2mup4*F(1638)+p3mup4*F(1639)))
       Fa3(278) = -(mup2mup5*(2*(p1mup3*F(263)+p2mup3*F(379))+p4mup3*F
     -   (383)+p3mup3*F(418)))+mup2mup3*(2*(p1mup5*F(263)+p2mup5*F(37
     -   9))+p4mup5*F(383)+p3mup5*F(418))+mup3mup5*(2*p1mup2*F(1262)+
     -   p4mup2*F(1609)+p3mup2*F(1622))
       Fa3(279) = mup2mup5*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(3
     -   83)+p3mup4*F(418))-mup2mup4*(2*(p1mup5*F(263)+p2mup5*F(379))
     -   +p4mup5*F(383)+p3mup5*F(418))-mup4mup5*(-2*p1mup2*F(1608)+p4
     -   mup2*F(1609)+p3mup2*F(1622))
       Fa3(280) = mup2mup3*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F(3
     -   83)+p3mup5*F(386))-mup2mup5*(2*(p1mup3*F(376)+p2mup3*F(377))
     -   +p3mup3*F(378)+p4mup3*F(1012))+mup3mup5*(2*p1mup2*F(1262)+p3
     -   mup2*F(1607)+p4mup2*F(1609))
       Fa3(281) = -(mup2mup4*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F
     -   (383)-p3mup5*F(385)))+mup2mup5*(2*(p1mup4*F(263)+p2mup4*F(37
     -   9))+p4mup4*F(383)+p3mup4*F(386))-mup4mup5*(p3mup2*F(1607)-2*
     -   p1mup2*F(1608)+p4mup2*F(1609))
       Fa3(282) = p1mup2*(2*(p1mup4*F(678)+p2mup4*F(680))+p3mup4*F(997
     -   )+p4mup4*F(1001))-p1mup4*(p3mup2*F(998)+p4mup2*F(1002))+2*(p
     -   2mup4*p3mup2*F(999)-p3mup4*p4mup2*F(1003)-p3mup2*p4mup4*F(10
     -   03)+p2mup4*p4mup2*F(1004)+2*(p3mup2*p3mup4*F(1000)+p4mup2*p4
     -   mup4*F(4472)))
       Fa3(283) = 4*p1mup4*(p4mup2*F(364)+p1mup2*F(634)+p3mup2*F(947))
     -   +p2mup4*p3mup2*F(990)+p2mup4*p4mup2*F(993)+2*(p1mup2*(p4mup4
     -   *F(684)+p2mup4*F(940)-p3mup4*F(991))+p3mup2*p3mup4*F(992)+p3
     -   mup4*p4mup2*F(994)+p3mup2*p4mup4*F(995)+p4mup2*p4mup4*F(996)
     -   )
       Fa3(284) = -4*(p3mup3*p3mup5*p4mup2*F(351)-p2mup5*p3mup2*p4mup3
     -   *F(379)+p1mup2*p3mup3*p4mup5*F(1261)-p1mup2*p3mup5*p4mup3*F(
     -   1262)-p2mup3*p3mup2*p4mup5*F(1265)+p2mup3*p3mup5*p4mup2*F(12
     -   66)+p2mup5*p3mup3*p4mup2*F(1267)-p3mup2*p3mup5*p4mup3*F(1268
     -   )-p3mup2*p4mup3*p4mup5*F(1271)+p3mup3*p4mup2*p4mup5*F(1272)+
     -   p3mup5*p4mup2*p4mup3*F(1273))-8*(p1mup2*p2mup3*p4mup5*F(635)
     -   +p2mup3*p2mup5*p4mup2*F(1264)+p2mup3*p4mup2*p4mup5*F(1269)+p
     -   4mup2*p4mup3*p4mup5*F(1275)+p1mup2*p4mup3*p4mup5*F(1276)+p2m
     -   up5*p4mup2*p4mup3*F(1277))+2*(mup2mup5*p2mup3*F(1250)-p1mup3
     -   *(-(mup2mup5*F(1251))+2*(p3mup5*p4mup2*F(344)-p3mup2*p4mup5*
     -   F(1263))+4*(p1mup2*p4mup5*F(633)+p2mup5*p4mup2*F(638)+p4mup2
     -   *p4mup5*F(1270)))+mup3mup5*p1mup2*F(1364))+p1mup5*(4*(p3mup2
     -   *p4mup3*F(263)-p4mup2*(p3mup3*F(348)+2*(p2mup3*F(255)+p4mup3
     -   *F(364)+p1mup3*F(634))))+mup2mup3*F(1365))+mup2mup3*p2mup5*F
     -   (1366)+mup3mup5*p3mup2*F(1640)+mup2mup5*p3mup3*F(1641)+mup2m
     -   up3*p3mup5*F(1642)+mup3mup5*p4mup2*F(1643)+mup2mup5*p4mup3*F
     -   (1644)+mup2mup3*p4mup5*F(1645)+p3mup2*p3mup3*p4mup5*F(1646)
       Fa3(285) = -4*(p1mup2*p2mup5*p4mup3*F(680)+p1mup2*p2mup3*p4mup5
     -   *F(1290)+p1mup2*p3mup5*p4mup3*F(1292)+p1mup2*p4mup3*p4mup5*F
     -   (1293)-p1mup2*p3mup3*p4mup5*F(1300))+mup3mup5*p1mup2*F(1647)
     -   +mup3mup5*p4mup2*F(1650)+mup2mup5*p4mup3*F(1651)+mup2mup3*p4
     -   mup5*F(1652)-2*(mup2mup5*p3mup3*F(1287)-mup2mup5*p2mup3*F(16
     -   53)-mup2mup3*p2mup5*F(1654)-mup3mup5*p3mup2*F(1655)-mup2mup3
     -   *p3mup5*F(1656))+8*(p2mup3*p3mup5*p4mup2*F(311)-p3mup2*p3mup
     -   5*p4mup3*F(311)+p3mup3*p3mup5*p4mup2*F(312)-p3mup2*p4mup3*p4
     -   mup5*F(1297)+p3mup3*p4mup2*p4mup5*F(1298)+p2mup3*p3mup2*p4mu
     -   p5*F(1299)+p2mup3*p2mup5*p4mup2*F(1303)-p2mup5*p3mup2*p4mup3
     -   *F(1303)+p2mup5*p3mup3*p4mup2*F(1304)-p2mup5*p4mup2*p4mup3*F
     -   (1305)-p3mup2*p3mup3*p4mup5*F(1306)+p3mup5*p4mup2*p4mup3*F(1
     -   307)+p2mup3*p4mup2*p4mup5*F(1658))+p1mup5*(mup2mup3*F(1649)-
     -   4*(p4mup3*(p1mup2*F(678)+p3mup2*F(1302))-p4mup2*(p3mup3*F(31
     -   0)+p2mup3*F(1657)+p4mup3*F(1659))))+p1mup3*(-4*(p1mup2*p4mup
     -   5*F(1289)-p4mup2*(p3mup5*F(293)+p1mup5*F(828)+p2mup5*F(1267)
     -   +p4mup5*F(1294)))+mup2mup5*F(1648)+p3mup2*p4mup5*F(1660))+16
     -   *p4mup2*p4mup3*p4mup5*F(1661)
       Fa3(286) = 16*(p2mup3*p4mup2*p4mup5*F(157)+p4mup2*p4mup3*p4mup5
     -   *F(175)+p3mup2*p3mup3*p3mup5*F(1351))+mup3mup5*p1mup2*F(1662
     -   )+mup2mup3*p1mup5*F(1664)+mup2mup5*p2mup3*F(1665)+mup3mup5*p
     -   3mup2*F(1666)+mup2mup5*p3mup3*F(1667)+mup2mup3*p3mup5*F(1668
     -   )+p1mup3*(mup2mup5*F(1663)+4*(4*p4mup2*p4mup5*F(153)+p1mup2*
     -   p3mup5*F(1320)+p2mup5*(4*p4mup2*F(63)+2*p1mup2*F(648)-p3mup2
     -   *F(1324))-p1mup5*p3mup2*F(1328)-p3mup2*p4mup5*F(1334)+p3mup5
     -   *p4mup2*F(1338)+2*(p1mup2*p1mup5*F(1326)+p1mup2*p4mup5*F(133
     -   1)+p1mup5*p4mup2*F(1343))+p3mup2*p3mup5*F(1672)))+4*(p1mup2*
     -   p1mup5*p3mup3*F(1321)+p1mup2*p2mup3*p3mup5*F(1322)+p1mup5*p2
     -   mup3*p3mup2*F(1325)+p1mup2*p3mup3*p3mup5*F(1329)+p1mup2*p3mu
     -   p3*p4mup5*F(1333)+p1mup2*p3mup5*p4mup3*F(1337)-p1mup5*p3mup2
     -   *p4mup3*F(1339)+p1mup5*p3mup3*p4mup2*F(1345)+p1mup5*p3mup2*p
     -   3mup3*F(1674))+8*(p3mup3*p4mup2*p4mup5*F(162)+p1mup2*p1mup5*
     -   p2mup3*F(649)+p2mup3*p3mup2*p4mup5*F(665)+p3mup3*p3mup5*p4mu
     -   p2*F(667)+p1mup2*p2mup3*p4mup5*F(1332)+p2mup3*p3mup5*p4mup2*
     -   F(1335)+p1mup2*p4mup3*p4mup5*F(1341)+p1mup5*p2mup3*p4mup2*F(
     -   1344)+p2mup3*p3mup2*p3mup5*F(1346)+p1mup2*p1mup5*p4mup3*F(13
     -   48)+p3mup2*p3mup3*p4mup5*F(1353)+p1mup5*p4mup2*p4mup3*F(1354
     -   )+p3mup5*p4mup2*p4mup3*F(1355)+p3mup2*p3mup5*p4mup3*F(1673)+
     -   p3mup2*p4mup3*p4mup5*F(1675))+2*(mup2mup5*p4mup3*F(1318)+mup
     -   3mup5*p4mup2*F(1670)+mup2mup3*p4mup5*F(1671)+p2mup5*(2*p3mup
     -   3*(p1mup2*F(1323)+2*(p4mup2*F(115)+p3mup2*F(1352)))+mup2mup3
     -   *F(1669)+4*(2*(p4mup2*p4mup3*F(171)+p2mup3*p4mup2*F(784))+p1
     -   mup2*p2mup3*F(1327)+p2mup3*p3mup2*F(1347)+p1mup2*p4mup3*F(13
     -   49)+p3mup2*p4mup3*F(1676))))
       Fa3(287) = p1mup5*(-4*(p3mup2*p4mup4*F(262)+p3mup4*p4mup2*F(263
     -   ))+2*mup2mup4*F(1259))-8*p2mup4*p2mup5*p3mup2*F(1264)-4*(p1m
     -   up2*p3mup4*p4mup5*F(1262)+p1mup2*p3mup5*p4mup4*F(1367)+p3mup
     -   4*p4mup2*p4mup5*F(1369))+mup4mup5*p3mup2*F(1677)+mup2mup5*p3
     -   mup4*F(1678)+mup2mup4*p3mup5*F(1679)+2*(mup4mup5*p1mup2*F(12
     -   58)+mup2mup4*p2mup5*F(1260)+mup2mup5*p2mup4*F(1680))+p1mup4*
     -   (-4*(p3mup2*p4mup5*F(260)-p3mup5*p4mup2*F(1368))+2*mup2mup5*
     -   F(1681))+mup4mup5*p4mup2*F(1682)+mup2mup5*p4mup4*F(1683)+mup
     -   2mup4*p4mup5*F(1684)+p3mup5*p4mup2*p4mup4*F(1685)
       Fa3(288) = 8*(p3mup5*p4mup2*p4mup4*F(936)+p3mup4*p4mup2*p4mup5*
     -   F(1386))+4*(p1mup2*p3mup4*p4mup5*F(1382)+p1mup2*p3mup5*p4mup
     -   4*F(1383)-p1mup4*p3mup5*p4mup2*F(1384)+p1mup4*p3mup2*p4mup5*
     -   F(1387))+mup4mup5*p1mup2*F(1686)+mup2mup5*p1mup4*F(1687)+p1m
     -   up5*(8*p1mup2*p1mup4*F(1326)+4*(p3mup2*p4mup4*F(912)+p3mup4*
     -   p4mup2*F(1385))+mup2mup4*F(1688))+p2mup4*(8*p2mup5*(2*p3mup2
     -   *F(771)+p1mup2*F(1327))+mup2mup5*F(1689))+mup2mup5*p3mup4*F(
     -   1690)+mup4mup5*p4mup2*F(1691)+mup2mup5*p4mup4*F(1692)+mup2mu
     -   p4*p4mup5*F(1693)+2*(mup2mup4*p2mup5*F(1694)+mup4mup5*p3mup2
     -   *F(1695)+mup2mup4*p3mup5*F(1696))
       Fa3(289) = 4*(p1mup2*p3mup4*p4mup5*F(1399)+p1mup2*p3mup5*p4mup4
     -   *F(1400))+8*(p3mup4*p4mup2*p4mup5*F(314)+p2mup4*p2mup5*p3mup
     -   2*F(692)-p3mup5*p4mup2*p4mup4*F(1402))+mup4mup5*p1mup2*F(169
     -   7)+p1mup5*(4*(p3mup2*p4mup4*F(313)+p3mup4*p4mup2*F(948))+mup
     -   2mup4*F(1699))+mup4mup5*p3mup2*F(1700)+mup2mup5*p3mup4*F(170
     -   1)+mup2mup4*p3mup5*F(1702)+2*(mup2mup5*p2mup4*F(1394)-mup2mu
     -   p4*p2mup5*F(1703)+mup4mup5*p4mup2*F(1704)+mup2mup5*p4mup4*F(
     -   1705)-mup2mup4*p4mup5*F(1706))+p1mup4*(4*p3mup2*p4mup5*F(298
     -   )+mup2mup5*F(1698)+p3mup5*p4mup2*F(1707))
       Fa3(290) = -2*(p1mup4*p3mup2*F(1410)+p2mup4*p3mup2*F(1411)+p3mu
     -   p4*p4mup2*F(1413))+p1mup2*(p2mup4*F(1708)-2*(p1mup4*F(1403)+
     -   p3mup4*F(1709))+p4mup4*F(1710))+p1mup4*p4mup2*F(1711)+p2mup4
     -   *p4mup2*F(1712)-4*p3mup2*p3mup4*F(1713)+p3mup2*p4mup4*F(1714
     -   )+p4mup2*p4mup4*F(1715)
       Fa3(291) = 8*(p2mup4*p2mup5*p3mup2*F(671)+p3mup4*p4mup2*p4mup5*
     -   F(1428)+p3mup5*p4mup2*p4mup4*F(1429)+p1mup4*p3mup2*p4mup5*F(
     -   1430)+p1mup5*p3mup2*p4mup4*F(1432))+16*(p1mup4*p3mup5*p4mup2
     -   *F(1431)+p1mup5*p3mup4*p4mup2*F(1433))+mup2mup4*p2mup5*F(171
     -   7)+mup4mup5*p4mup2*F(1720)+mup2mup5*(p2mup4*F(1716)-2*(p1mup
     -   4*F(1418)+p3mup4*F(1719))+p4mup4*F(1721))+mup2mup4*p4mup5*F(
     -   1722)-2*(mup2mup4*p1mup5*F(1424)+mup2mup4*p3mup5*F(1427)-p1m
     -   up2*(8*p1mup4*p1mup5*F(12)+4*(p3mup5*p4mup4*F(916)+p2mup4*p2
     -   mup5*F(1425)+p3mup4*p4mup5*F(1426))+mup4mup5*F(1718))-mup4mu
     -   p5*p3mup2*F(1723))
       Fa3(292) = -4*p1mup2*(2*p1mup4*p1mup5*F(136)+p2mup4*p2mup5*F(14
     -   45)+4*p3mup4*p4mup5*F(1448))-8*(p2mup4*p2mup5*p3mup2*F(665)+
     -   p1mup4*p3mup2*p4mup5*F(1429)-p1mup5*p3mup4*p4mup2*F(1447)+p1
     -   mup5*p3mup2*p4mup4*F(1449))-16*p3mup4*p4mup2*p4mup5*F(1450)+
     -   mup4mup5*p1mup2*F(1724)+mup2mup5*p1mup4*F(1725)+mup2mup4*p1m
     -   up5*F(1726)+mup2mup5*p2mup4*F(1727)+mup2mup4*p2mup5*F(1728)+
     -   mup4mup5*p3mup2*F(1729)+2*(mup2mup5*p4mup4*F(1443)+mup2mup4*
     -   p4mup5*F(1444)-p3mup5*(4*p1mup4*p4mup2*F(1429)-mup2mup4*F(14
     -   41)+8*(p4mup2*p4mup4*F(120)+p1mup2*p4mup4*F(1446)))+mup2mup5
     -   *p3mup4*F(1730))+mup4mup5*p4mup2*F(1731)
       Fa3(293) = p1mup4*p3mup2*F(1733)+p1mup2*(p3mup4*F(1732)+2*p1mup
     -   4*F(1734)+p4mup4*F(1736))+p1mup4*p4mup2*F(1737)+p3mup2*p3mup
     -   4*F(1738)+p2mup4*(2*(p4mup2*F(1461)+p1mup2*F(1735))+p3mup2*F
     -   (1739))+p3mup2*p4mup4*F(1740)+p3mup4*p4mup2*F(1741)+p4mup2*p
     -   4mup4*F(1742)
       Fa3(294) = 8*(p1mup2*p3mup4*p4mup5*F(354)+p1mup4*p3mup5*p4mup2*
     -   F(359)+p1mup2*p3mup5*p4mup4*F(359)+p3mup4*p4mup2*p4mup5*F(36
     -   1)+p1mup5*p3mup4*p4mup2*F(364)-p3mup5*p4mup2*p4mup4*F(365))-
     -   4*p2mup4*p2mup5*p3mup2*F(1470)+mup4mup5*p3mup2*F(1746)+mup2m
     -   up5*(p1mup4*F(341)+p4mup4*F(988)+p2mup4*F(1743)+p3mup4*F(174
     -   7))+mup2mup4*(p1mup5*F(342)+p4mup5*F(987)+p2mup5*F(1744)+p3m
     -   up5*F(1748))+2*(mup4mup5*p1mup2*F(1745)+mup4mup5*p4mup2*F(17
     -   49))
       Fa3(295) = p2mup3*p3mup2*F(1751)+p1mup2*(4*p1mup3*F(1750)+2*p2m
     -   up3*F(1752)+p3mup3*F(1753))+2*(p1mup3*p3mup2*F(1754)+p3mup2*
     -   p3mup3*F(1755)+p4mup2*p4mup3*F(1756))
       Fa3(296) = p1mup2*(2*(p1mup3*F(1757)+p2mup3*F(1758))+p3mup3*F(1
     -   759))+p1mup3*p3mup2*F(1760)+2*p2mup3*p3mup2*F(1761)+p4mup2*p
     -   4mup3*F(1762)+p3mup2*p3mup3*F(1763)
       Fa3(297) = p1mup2*(p2mup3*F(1765)+2*(p1mup3*F(1764)+p3mup3*F(17
     -   67)))+p4mup2*p4mup3*F(1769)+2*p3mup2*(p2mup3*F(1766)+p1mup3*
     -   F(1768)+2*p3mup3*F(1770))
       Return
       End
