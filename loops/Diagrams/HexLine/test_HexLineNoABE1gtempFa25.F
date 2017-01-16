c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g25(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa25)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/test_FhlFunctionsNoAbe1g/F
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
      Complex*16 Fa(587),F(6395)
      Real*8 P(7129) 
        COMMON/test_PFunctionsNoAbe1g/P
       Complex*16 Fa25(337:350)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa25(337) = p5mup4*(-16*p1mup2*(-(p4mup6*(p4mup3*F(1331)+p5mup3
     -   *F(1336)))+p5mup6*(p3mup3*F(1322)+p4mup3*F(1332)+p5mup3*F(13
     -   37)))+p1mup2*p3mup6*p5mup3*F(1408)+p2mup2*(16*(p4mup3*p4mup6
     -   *F(1331)-p4mup3*p5mup6*F(1332)+p4mup6*p5mup3*F(1336)-p5mup3*
     -   p5mup6*F(1337)+p2mup3*p2mup6*F(1345)+p2mup3*p4mup6*F(1356)-p
     -   2mup3*p5mup6*F(1358)+p2mup6*p4mup3*F(1360)-p3mup3*(p4mup6*F(
     -   1319)+p5mup6*F(1322)+p2mup6*F(1351)+p3mup6*F(1364)))+p2mup6*
     -   p5mup3*F(1407)+p3mup6*(16*(p4mup3*F(1330)+p2mup3*F(1346))+p5
     -   mup3*F(1408))))
       Fa25(338) = 4*p3mup2*s12*(p4mup3*F(1389)+p3mup3*F(1390)+p5mup3*
     -   F(1391)+(p1mup3+p2mup3)*F(1392))+(p1mup2+p2mup2)*(4*((p1mup3
     -   +p2mup3)*F(1409)-p3mup3*F(1410))+p4mup3*F(1411)+p5mup3*F(141
     -   2))
       Fa25(339) = p5mup4*(p2mup2*(p4mup3*F(1413)+p5mup3*F(1414)+p3mup
     -   3*F(1415)+p2mup3*F(1420))-8*s12*(p2mup3*p4mup2*F(1416)-p3mup
     -   2*p4mup3*F(1417)+p3mup3*p4mup2*F(1418)+p4mup2*p4mup3*F(1419)
     -   +p3mup2*p3mup3*F(1421)+p3mup2*p5mup3*F(1422)+p4mup2*p5mup3*F
     -   (1423)+p2mup3*p3mup2*F(1424)))
       Fa25(340) = p1mup2*p4mup3*p5mup4*F(1413)+p1mup2*p5mup3*p5mup4*F
     -   (1414)+p1mup2*p3mup3*p5mup4*F(1415)+p1mup2*p1mup3*p5mup4*F(1
     -   420)+p1mup3*p2mup2*p5mup4*F(1420)+p1mup2*p2mup3*p5mup4*F(142
     -   0)+8*(-(p1mup3*p4mup2*p5mup4*s12*F(1416))-p1mup3*p3mup2*p5mu
     -   p4*s12*F(1424)+p1mup2*p1mup3*p4mup4*F(1425)+p1mup3*p2mup2*p4
     -   mup4*F(1425)+p1mup2*p2mup3*p4mup4*F(1425)+p2mup2*p2mup3*p4mu
     -   p4*F(1425)+(p1mup2+p2mup2)*p3mup3*p4mup4*F(1426))
       Fa25(341) = (p1mup2+p2mup2)*(p3mup4*(p5mup3*F(1428)+8*(p1mup3+p
     -   2mup3)*F(1430)+p4mup3*F(1431))+p4mup4*(p5mup3*F(1429)+8*p4mu
     -   p3*F(1432))+(p1mup4+p2mup4)*p5mup3*F(1433)+8*(p3mup3*p3mup4*
     -   F(1427)-mup3mup4*F(1434)))
       Fa25(342) = p1mup2*(8*p1mup3*p2mup4*F(1436)+p1mup4*(p4mup3*F(14
     -   35)+8*(p2mup3*F(1436)+p3mup3*F(1446))))+p2mup2*(8*p1mup3*p1m
     -   up4*F(1436)+p2mup4*(p4mup3*F(1435)+8*(p2mup3*F(1436)+p3mup3*
     -   F(1446))))-8*s12*(p3mup2*(p2mup3*p4mup4*F(1440)+p1mup3*p1mup
     -   4*F(1447)+p2mup4*(p5mup3*F(1443)+p2mup3*F(1447)+p3mup3*F(145
     -   0)-p4mup3*F(1453))+p3mup4*(p5mup3*F(1437)+p3mup3*F(1444)+p2m
     -   up3*F(1449)-p4mup3*F(1454))+p4mup4*(p5mup3*F(1445)+p3mup3*F(
     -   1451)+p4mup3*F(1456))+mup3mup4*F(1462))+p4mup2*(p1mup3*p1mup
     -   4*F(1448)+p2mup3*(p4mup4*F(1438)+p3mup4*F(1441)+p2mup4*F(144
     -   8))+p3mup3*(p4mup4*F(1439)+p2mup4*F(1442)+p3mup4*F(1452))+p2
     -   mup4*p4mup3*F(1455)+p3mup4*p4mup3*F(1457)+p4mup3*p4mup4*F(14
     -   58)+p4mup4*p5mup3*F(1459)+p2mup4*p5mup3*F(1460)+p3mup4*p5mup
     -   3*F(1461)+mup3mup4*F(1463)))
       Fa25(343) = p1mup2*p2mup4*(p4mup3*F(1435)+8*(p2mup3*F(1436)+p3m
     -   up3*F(1446)))+8*p1mup3*(p2mup2*p2mup4*F(1436)-s12*(p4mup2*p4
     -   mup4*F(1438)+p3mup2*p4mup4*F(1440)+p3mup4*p4mup2*F(1441)+p2m
     -   up4*p3mup2*F(1447)+p2mup4*p4mup2*F(1448)+p3mup2*p3mup4*F(144
     -   9)))+p1mup4*(p2mup2*(p4mup3*F(1435)+8*(p2mup3*F(1436)+p3mup3
     -   *F(1446)))+8*(p1mup2*p1mup3*F(1436)-s12*(p3mup3*p4mup2*F(144
     -   2)+p3mup2*p5mup3*F(1443)+p2mup3*p3mup2*F(1447)+p2mup3*p4mup2
     -   *F(1448)+p3mup2*p3mup3*F(1450)-p3mup2*p4mup3*F(1453)+p4mup2*
     -   p4mup3*F(1455)+p4mup2*p5mup3*F(1460))))
       Fa25(344) = mup4mup6*(-(p4mup3*F(1220))+p5mup3*F(1223)+(p1mup3+
     -   p2mup3)*F(1226)+p3mup3*F(1227))+mup3mup6*p5mup4*F(1377)
       Fa25(345) = mup3mup6*(p3mup4*F(1228)-p4mup4*F(1229))
       Fa25(346) = 4*p5mup6*(-(mup3mup4*F(1242))+2*(p2mup3*p2mup4*F(12
     -   54)+p2mup3*p3mup4*F(1255)+p2mup3*p4mup4*F(1256)-p2mup4*p4mup
     -   3*F(1258)-p2mup4*p3mup3*F(1265)-p3mup3*p3mup4*F(1266)-p3mup3
     -   *p4mup4*F(1267)+p3mup4*p4mup3*F(1271)+p4mup3*p4mup4*F(1278)+
     -   p2mup4*p5mup3*F(1283)+p1mup4*(p2mup3*F(1254)-p4mup3*F(1258)-
     -   p3mup3*F(1265)+p5mup3*F(1283))+p3mup4*p5mup3*F(1284)+p3mup3*
     -   p5mup4*F(1284)-p4mup4*p5mup3*F(1291)+p4mup3*p5mup4*F(1321)-p
     -   5mup3*p5mup4*F(1335)+p2mup3*p5mup4*F(1353)+p1mup3*((p1mup4+p
     -   2mup4)*F(1254)+p3mup4*F(1255)+p4mup4*F(1256)+p5mup4*F(1353))
     -   ))
       Fa25(347) = mup3mup6*(p1mup4+p2mup4)*F(1464)+8*p5mup4*(p4mup3*(
     -   p4mup6*F(1318)+(p1mup6+p2mup6)*F(1350)+p3mup6*F(1363))+p4mup
     -   6*(p5mup3*F(1333)+(p1mup3+p2mup3)*F(1343)+p3mup3*F(1366))+p5
     -   mup3*(p3mup6*F(1367)+p2mup6*F(1465)))
       Fa25(348) = -(mup3mup4*(p3mup6*F(1240)+p4mup6*F(1380)))+mup3mup
     -   4*(p1mup6+p2mup6)*F(1381)+2*((p1mup3*p1mup6+p2mup3*p2mup6)*(
     -   (p1mup4+p2mup4)*F(1244)+p3mup4*F(1245)+p4mup4*F(1246))-p2mup
     -   4*p2mup6*p4mup3*F(1251)+p2mup3*(p2mup4*p3mup6*F(1247)+p3mup4
     -   *p3mup6*F(1248)+p3mup6*p4mup4*F(1249)+p2mup4*p4mup6*F(1250)+
     -   p3mup4*p4mup6*F(1252)+p4mup4*p4mup6*F(1253))-p2mup4*p3mup6*p
     -   4mup3*F(1257)+p2mup4*p2mup6*p3mup3*F(1260)+p2mup6*p3mup3*p3m
     -   up4*F(1261)+p2mup6*p3mup4*p4mup3*F(1262)+p2mup4*p3mup3*p4mup
     -   6*F(1263)+p3mup3*p3mup4*p4mup6*F(1264)+p2mup4*p3mup3*p3mup6*
     -   F(1268)+p3mup3*p3mup4*p3mup6*F(1269)+p3mup4*p3mup6*p4mup3*F(
     -   1270)-p2mup6*p4mup4*p5mup3*F(1273)-p3mup6*p4mup4*p5mup3*F(12
     -   74)-p4mup4*p4mup6*p5mup3*F(1275)-p2mup4*p4mup6*p5mup3*F(1276
     -   )+p2mup6*p3mup3*p4mup4*F(1279)-p2mup6*p4mup3*p4mup4*F(1280)+
     -   p3mup6*p4mup3*p4mup4*F(1281)-p2mup4*p2mup6*p5mup3*F(1285)+p3
     -   mup3*p4mup4*p4mup6*F(1287)+p2mup4*p3mup6*p5mup3*F(1289)+p1mu
     -   p4*(p2mup3*p3mup6*F(1247)+p2mup3*p4mup6*F(1250)-p2mup6*p4mup
     -   3*F(1251)-p3mup6*p4mup3*F(1257)+p2mup6*p3mup3*F(1260)+p3mup3
     -   *p4mup6*F(1263)+p3mup3*p3mup6*F(1268)-p4mup6*p5mup3*F(1276)-
     -   p2mup6*p5mup3*F(1285)+p3mup6*p5mup3*F(1289))+p3mup4*p3mup6*p
     -   5mup3*F(1290)+p2mup6*p3mup4*p5mup3*F(1382)+p3mup4*p4mup6*p5m
     -   up3*F(1383)+p1mup6*p5mup3*p5mup4*F(1465))
       Fa25(349) = p1mup3*p2mup4*p2mup6*F(1244)+p1mup3*p2mup6*p3mup4*F
     -   (1245)+p1mup3*p2mup6*p4mup4*F(1246)+p1mup3*p2mup4*p3mup6*F(1
     -   247)+p1mup3*p3mup4*p3mup6*F(1248)+p1mup3*p3mup6*p4mup4*F(124
     -   9)+p1mup3*p2mup4*p4mup6*F(1250)+p1mup3*p3mup4*p4mup6*F(1252)
     -   +p1mup3*p4mup4*p4mup6*F(1253)-p2mup4*p4mup3*p4mup6*F(1259)+p
     -   3mup4*p4mup3*p4mup6*F(1272)+p4mup3*p4mup4*p4mup6*F(1282)+p1m
     -   up4*(p1mup3*p2mup6*F(1244)+p1mup3*p3mup6*F(1247)+p1mup3*p4mu
     -   p6*F(1250)-p4mup3*p4mup6*F(1259)+p1mup6*(p2mup3*F(1244)-p4mu
     -   p3*F(1251)+p3mup3*F(1260)-p5mup3*F(1285)))+p3mup3*p3mup6*p4m
     -   up4*F(1288)+p1mup6*(p2mup3*p2mup4*F(1244)+p2mup3*p3mup4*F(12
     -   45)+p2mup3*p4mup4*F(1246)-p2mup4*p4mup3*F(1251)+p2mup4*p3mup
     -   3*F(1260)+p3mup3*p3mup4*F(1261)+p3mup4*p4mup3*F(1262)-p4mup4
     -   *p5mup3*F(1273)+p3mup3*p4mup4*F(1279)-p4mup3*p4mup4*F(1280)-
     -   p2mup4*p5mup3*F(1285)+p3mup4*p5mup3*F(1382))
       Fa25(350) = p5mup4*(p3mup3*p3mup6*F(1290)+(p1mup3+p2mup3)*((p1m
     -   up6+p2mup6)*F(1340)+p3mup6*F(1342))+(p1mup6+p2mup6)*p3mup3*F
     -   (1347))
       Return
       End
