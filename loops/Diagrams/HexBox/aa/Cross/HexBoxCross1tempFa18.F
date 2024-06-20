c************** Calling the Fa functions*************************
       subroutine HexBoxCross1Fa18(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,
     &   p3mup4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa
     &   18)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup4, p2mup3, 
     -          p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, 
     -          p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4

        common/HexBoxCross1FhlFunctions/F
       COMMON/HexBoxCross1Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
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
      Complex*16 Fa(272),F(1461)
      Real*8 P(1825) 
        COMMON/HexBoxCross1PFunctions/P
       Complex*16 Fa18(222:234)
       COMMON/HexBoxCross1FaFunctions/Fa
       Fa18(222) = p1mup3*F(1152)+2*(p4mup3*F(1151)+p5mup3*F(1153)+p6m
     &   up3*F(1154)+p2mup3*F(1155))
       Fa18(223) = p4mup3*F(1157)+2*(p5mup3*F(1156)+p1mup3*F(1158)+p6m
     &   up3*F(1159)+p2mup3*F(1160))
       Fa18(224) = p5mup3*F(1163)+2*(p4mup3*F(1161)+p1mup3*F(1162)+p6m
     &   up3*F(1164)-p2mup3*F(1165))
       Fa18(225) = 2*(mup3mup4*F(1168)+2*(-(p5mup4*p6mup3*F(559))-p2mu
     &   p4*p6mup3*F(600)+p2mup3*p5mup4*F(883)+p1mup3*(p5mup4*F(556)-
     &   p2mup4*F(597)+p6mup4*F(1169)+p4mup4*F(1170)-p1mup4*F(1171))+
     &   p2mup3*p4mup4*F(1172)+p2mup3*p6mup4*F(1173)-p4mup4*p6mup3*F(
     &   1178)-p6mup3*p6mup4*F(1179)-p2mup3*p2mup4*F(1180)-p4mup3*(p5
     &   mup4*F(557)+p2mup4*F(599)+p4mup4*F(1174)+p6mup4*F(1177)+p1mu
     &   p4*F(1181))-p5mup3*(p5mup4*F(558)+p2mup4*F(560)+p6mup4*F(117
     &   5)+p4mup4*F(1176)+p1mup4*F(1182))-p1mup4*p6mup3*F(1183)-p1mu
     &   p4*p2mup3*F(1184)))
       Fa18(226) = mup3mup4*F(1185)-4*(p1mup3*(p5mup4*F(1186)+p6mup4*F
     &   (1187)+p4mup4*F(1188)+p2mup4*F(1192)-p1mup4*F(1194))+p2mup3*
     &   (p4mup4*F(1189)-p6mup4*F(1190)+p5mup4*F(1191)+p1mup4*F(1192)
     &   +p2mup4*F(1201))+p4mup3*(p4mup4*F(997)-p1mup4*F(1003)+p5mup4
     &   *F(1193)+p6mup4*F(1197)-p2mup4*F(1203))+p5mup3*(p5mup4*F(623
     &   )-p1mup4*F(626)+p6mup4*F(1195)+p4mup4*F(1196)-p2mup4*F(1204)
     &   )+p6mup3*(p5mup4*F(1198)+p4mup4*F(1199)+p6mup4*F(1200)-p1mup
     &   4*F(1202)-p2mup4*F(1205)))
       Fa18(227) = mup3mup4*F(1206)-4*(p2mup3*(p5mup4*F(635)+p6mup4*F(
     &   1072)+p4mup4*F(1073)-p1mup4*F(1076)-p2mup4*F(1078))+p1mup3*(
     &   p5mup4*F(638)-p2mup4*F(1076)+p6mup4*F(1207)+p4mup4*F(1208)+p
     &   1mup4*F(1209))-p4mup3*(p5mup4*F(646)+p4mup4*F(700)-p1mup4*F(
     &   707)-p2mup4*F(1004)+p6mup4*F(1212))-p5mup3*(p5mup4*F(647)-p2
     &   mup4*F(654)+p6mup4*F(1210)+p4mup4*F(1211)-p1mup4*F(1214))-p6
     &   mup3*(p5mup4*F(648)+p6mup4*F(850)-p2mup4*F(1081)+p4mup4*F(12
     &   13)-p1mup4*F(1215)))
       Fa18(228) = mup3mup4*F(1216)+4*(p4mup3*(p5mup4*F(647)-p2mup4*F(
     &   654)+p6mup4*F(1210)+p4mup4*F(1211)-p1mup4*F(1214))-p1mup3*(p
     &   5mup4*F(1220)+p6mup4*F(1221)+p4mup4*F(1222)-p2mup4*F(1223)+p
     &   1mup4*F(1224))+p5mup3*p5mup4*F(1225)+p5mup3*p6mup4*F(1226)+p
     &   4mup4*p5mup3*F(1227)+p5mup4*p6mup3*F(1228)+p4mup4*p6mup3*F(1
     &   229)+p6mup3*p6mup4*F(1230)-p2mup3*(p4mup4*F(1217)+p6mup4*F(1
     &   218)+p5mup4*F(1219)-p1mup4*F(1223)-p2mup4*F(1231))-p1mup4*p5
     &   mup3*F(1232)-p1mup4*p6mup3*F(1233)-p2mup4*p5mup3*F(1234)-p2m
     &   up4*p6mup3*F(1235))
       Fa18(229) = mup3mup4*F(1236)-4*(-(p5mup4*p6mup3*F(651))-p1mup4*
     &   p2mup3*F(991)-p2mup3*p5mup4*F(999)-p2mup3*p6mup4*F(1000)+p2m
     &   up3*p4mup4*F(1005)-p4mup4*p6mup3*F(1079)-p6mup3*p6mup4*F(108
     &   0)+p1mup4*p6mup3*F(1081)+p1mup3*(p5mup4*F(635)+p6mup4*F(1072
     &   )+p4mup4*F(1073)-p1mup4*F(1076)+p2mup4*F(1189))-p4mup3*(p5mu
     &   p4*F(649)+p4mup4*F(1001)+p6mup4*F(1002)-p1mup4*F(1004)-p2mup
     &   4*F(1237))-p5mup3*(p5mup4*F(650)+p4mup4*F(652)+p6mup4*F(653)
     &   -p1mup4*F(654)-p2mup4*F(1238))+p2mup4*p6mup3*F(1239)+p2mup3*
     &   p2mup4*F(1240))
       Fa18(230) = 2*mup3mup4*F(1241)+p6mup3*p6mup4*F(1251)+4*(p4mup3*
     &   (p4mup4*F(704)+p5mup4*F(705)+p6mup4*F(706))-p2mup3*(p5mup4*F
     &   (702)+p6mup4*F(703)+p4mup4*F(708))+p1mup4*(p2mup3*F(691)-p4m
     &   up3*F(709))+p2mup4*(p5mup3*F(652)+p4mup3*F(1001)-p2mup3*F(10
     &   05)-p1mup3*F(1073)+p6mup3*F(1079))-p1mup3*(p5mup4*F(1242)-p6
     &   mup4*F(1243)+p4mup4*F(1244)+p1mup4*F(1245))+p5mup3*p5mup4*F(
     &   1246)+p5mup3*p6mup4*F(1247)+p4mup4*p5mup3*F(1248)+p5mup4*p6m
     &   up3*F(1249)+p4mup4*p6mup3*F(1250)-p1mup4*p5mup3*F(1252)-p1mu
     &   p4*p6mup3*F(1253))
       Fa18(231) = -2*mup3mup4*F(1254)+4*(p4mup3*(p5mup4*F(1246)+p6mup
     &   4*F(1247)+p4mup4*F(1248)-p1mup4*F(1252))+p1mup4*p2mup3*F(125
     &   5)-p1mup3*(p5mup4*F(1256)+p6mup4*F(1257)+p4mup4*F(1258))-p1m
     &   up3*p1mup4*F(1259)+p2mup4*(p4mup3*F(652)-p2mup3*F(655)-p1mup
     &   3*F(1217)+p5mup3*F(1260)+p6mup3*F(1261))-p2mup3*p6mup4*F(126
     &   2)-p2mup3*p5mup4*F(1263)+p5mup3*p5mup4*F(1264)+p5mup3*p6mup4
     &   *F(1265)+p4mup4*p5mup3*F(1266)+p5mup4*p6mup3*F(1267)+p4mup4*
     &   p6mup3*F(1268)+p6mup3*p6mup4*F(1269)-p1mup4*p5mup3*F(1270)-p
     &   1mup4*p6mup3*F(1271)-p2mup3*p4mup4*F(1272))
       Fa18(232) = mup3mup4*F(1273)+4*(-(p1mup4*p5mup3*F(1195))-p1mup4
     &   *p6mup3*F(1200)+p1mup4*p2mup3*F(1205)-p2mup3*p4mup4*F(1239)-
     &   p1mup3*(p1mup4*F(1187)-p2mup4*F(1190)+p5mup4*F(1274)+p6mup4*
     &   F(1275)+p4mup4*F(1276))+p5mup3*p5mup4*F(1278)+p5mup3*p6mup4*
     &   F(1279)+p4mup4*p5mup3*F(1280)+p5mup4*p6mup3*F(1281)+p4mup4*p
     &   6mup3*F(1282)+p6mup3*p6mup4*F(1283)+p4mup3*(p4mup4*F(989)+p6
     &   mup4*F(996)-p1mup4*F(1197)+p5mup4*F(1277)-p2mup4*F(1284))-p2
     &   mup4*p5mup3*F(1285)-p2mup4*p6mup3*F(1286)-p2mup3*p6mup4*F(12
     &   86)-p2mup3*p5mup4*F(1287)-p2mup3*p2mup4*F(1288))
       Fa18(233) = 2*mup3mup4*F(1289)+p4mup3*p6mup4*F(1294)+4*(p4mup3*
     &   p5mup4*F(642)+p5mup4*p6mup3*F(645)+p2mup3*p5mup4*F(651)+p4mu
     &   p3*p4mup4*F(694)+p1mup4*p6mup3*F(850)+p2mup4*p4mup3*F(1002)+
     &   p2mup3*p4mup4*F(1079)+p2mup4*p6mup3*F(1080)+p2mup3*p6mup4*F(
     &   1080)-p1mup4*p2mup3*F(1081)-p2mup3*p2mup4*F(1082)+p1mup4*p4m
     &   up3*F(1212)+p1mup3*(p5mup4*F(633)-p2mup4*F(1072)-p1mup4*F(12
     &   07)+p6mup4*F(1290)+p4mup4*F(1291))+p5mup3*(p5mup4*F(640)+p2m
     &   up4*F(653)+p1mup4*F(1210)+p6mup4*F(1292)+p4mup4*F(1293))+p4m
     &   up4*p6mup3*F(1295)+p6mup3*p6mup4*F(1296))
       Fa18(234) = 2*mup3mup4*F(1297)+p5mup3*p6mup4*F(1302)+4*(p1mup4*
     &   p5mup3*F(1226)+p1mup4*p6mup3*F(1230)-p1mup4*p2mup3*F(1235)+p
     &   2mup3*p4mup4*F(1261)+p4mup3*(p5mup4*F(640)+p2mup4*F(653)+p1m
     &   up4*F(1210)+p6mup4*F(1292)+p4mup4*F(1293))+p1mup3*(-(p2mup4*
     &   F(1218))-p1mup4*F(1221)+p5mup4*F(1298)+p6mup4*F(1299)+p4mup4
     &   *F(1300))+p5mup3*p5mup4*F(1301)+p4mup4*p5mup3*F(1303)+p5mup4
     &   *p6mup3*F(1304)+p4mup4*p6mup3*F(1305)+p6mup3*p6mup4*F(1306)+
     &   p2mup4*p5mup3*F(1307)+p2mup4*p6mup3*F(1308)+p2mup3*p6mup4*F(
     &   1308)+p2mup3*p5mup4*F(1309)-p2mup3*p2mup4*F(1310))
       Return
       End
