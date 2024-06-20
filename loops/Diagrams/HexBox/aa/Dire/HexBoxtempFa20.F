c************** Calling the Fa functions*************************
       subroutine HexBoxFa20(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa20)
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
       Complex*16 Fa20(229:259)
       COMMON/HexBoxFaFunctions/Fa
       Fa20(229) = p6mup4*F(1192)+p4mup4*F(1193)+p5mup4*F(1194)+p2mup4
     &   *F(1195)-p1mup4*F(1196)
       Fa20(230) = p5mup3*F(1124)-p2mup3*F(1197)-p1mup3*F(1198)+p6mup3
     &   *F(1199)+p4mup3*F(1200)
       Fa20(231) = p5mup3*F(1130)+p2mup3*F(1184)+p1mup3*F(1201)+p6mup3
     &   *F(1202)-p4mup3*F(1203)
       Fa20(232) = p2mup4*F(1205)-p6mup4*F(1206)+p4mup4*F(1207)+p5mup4
     &   *F(1208)-p1mup4*F(1209)
       Fa20(233) = p2mup3*F(1210)-p4mup3*F(1211)-p5mup3*F(1212)+p1mup3
     &   *F(1213)+p6mup3*F(1214)
       Fa20(234) = p6mup4*F(1215)+p4mup4*F(1216)+p5mup4*F(1217)-p1mup4
     &   *F(1218)-p2mup4*F(1219)
       Fa20(235) = p6mup4*F(1220)+p4mup4*F(1221)+p5mup4*F(1222)+p1mup4
     &   *F(1223)-p2mup4*F(1224)
       Fa20(236) = p4mup4*F(1170)-p1mup4*F(1225)+p6mup4*F(1226)+p5mup4
     &   *F(1227)-p2mup4*F(1228)
       Fa20(237) = p6mup4*F(1229)+p4mup4*F(1230)+p5mup4*F(1231)-p2mup4
     &   *F(1232)-p1mup4*F(1233)
       Fa20(238) = -(p1mup4*F(1167))+p6mup4*F(1168)+p5mup4*F(1169)-p2m
     &   up4*F(1170)+p4mup4*F(1171)
       Fa20(239) = p6mup4*F(1234)+p4mup4*F(1235)+p5mup4*F(1236)-p2mup4
     &   *F(1237)-p1mup4*F(1238)
       Fa20(240) = -(p1mup4*F(1239))+p6mup4*F(1240)+p5mup4*F(1241)-p2m
     &   up4*F(1242)+p4mup4*F(1243)
       Fa20(241) = p2mup4*F(784)-p6mup4*F(1245)-p4mup4*F(1246)-p5mup4*
     &   F(1247)+p1mup4*F(1248)
       Fa20(242) = p2mup4*F(1251)+p6mup4*F(1252)+p4mup4*F(1253)+p5mup4
     &   *F(1254)-p1mup4*F(1255)
       Fa20(243) = p6mup4*F(1256)+p4mup4*F(1257)+p5mup4*F(1258)-p1mup4
     &   *F(1259)-p2mup4*F(1260)
       Fa20(244) = p1mup4*F(1261)+p2mup4*F(1262)+p6mup4*F(1263)-p4mup4
     &   *F(1264)+p5mup4*F(1265)
       Fa20(245) = p1mup4*F(1267)-p2mup4*F(1268)+p6mup4*F(1269)+p4mup4
     &   *F(1270)+p5mup4*F(1271)
       Fa20(246) = p2mup4*F(1270)+p1mup4*F(1272)-p6mup4*F(1273)-p4mup4
     &   *F(1274)+p5mup4*F(1275)
       Fa20(247) = p2mup4*F(1276)+p6mup4*F(1277)+p4mup4*F(1278)+p5mup4
     &   *F(1279)+p1mup4*F(1280)
       Fa20(248) = p1mup3*F(1187)+p6mup3*F(1189)+p4mup3*F(1190)+p5mup3
     &   *F(1191)+p2mup3*F(1283)
       Fa20(249) = p5mup3*F(1125)-p2mup3*F(1182)+p6mup3*F(1200)-p1mup3
     &   *F(1284)+p4mup3*F(1285)
       Fa20(250) = p1mup4*F(783)+p2mup4*F(784)-p4mup4*F(797)-p5mup4*F(
     &   798)-p6mup4*F(1313)
       Fa20(251) = p2mup3*F(793)-p1mup3*F(826)+p4mup3*F(831)+p5mup3*F(
     &   834)+p6mup3*F(835)
       Fa20(252) = -(p1mup3*F(1289))+p6mup3*F(1290)+p4mup3*F(1291)+p5m
     &   up3*F(1292)-p2mup3*F(1293)
       Fa20(253) = 4*(p5mup3*F(1130)+p2mup3*F(1184)+p1mup3*F(1201)+p6m
     &   up3*F(1202))+p4mup3*F(1294)
       Fa20(254) = p5mup3*F(1296)+4*(p4mup3*F(1130)+p6mup3*F(1131)+p1m
     &   up3*F(1295)+p2mup3*F(1297))
       Fa20(255) = 2*(p5mup3*F(1130)+p2mup3*F(1184)+p1mup3*F(1201)+p6m
     &   up3*F(1202))+p4mup3*F(1298)
       Fa20(256) = 2*(p4mup3*F(1130)+p6mup3*F(1131)+p1mup3*F(1295)+p2m
     &   up3*F(1297))+p5mup3*F(1299)
       Fa20(257) = -(p1mup3*F(1300))-p2mup3*F(1301)+p6mup3*F(1302)+p4m
     &   up3*F(1303)+p5mup3*F(1304)
       Fa20(258) = -(p1mup3*F(783))-p2mup3*F(784)+p4mup3*F(797)+p5mup3
     &   *F(798)+p6mup3*F(1313)
       Fa20(259) = p2mup3*F(792)+p1mup3*F(827)+p4mup3*F(838)+p5mup3*F(
     &   839)+p6mup3*F(840)
       Return
       End
