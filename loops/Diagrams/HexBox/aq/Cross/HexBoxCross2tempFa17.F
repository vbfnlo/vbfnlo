c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa17(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   17)
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
       Complex*16 Fa17(305:323)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa17(305) = -(p1mup3*F(1318))-p2mup3*F(1319)+p4mup3*F(1320)+2*p
     &   5mup3*F(1321)+p6mup3*F(1852)
       Fa17(306) = p2mup3*F(1262)-p4mup3*F(1263)+p5mup3*F(1265)+p1mup3
     &   *F(1266)-p6mup3*F(1322)
       Fa17(307) = -(p2mup3*F(1323))+p4mup3*F(1324)-p5mup3*F(1325)+p6m
     &   up3*F(1326)+p1mup3*F(1327)
       Fa17(308) = p1mup3*F(1292)-p2mup3*F(1306)+p4mup3*F(1328)
       Fa17(309) = -(p6mup3*p6mup6*F(1335))-p4mup3*p6mup6*F(1336)+p2mu
     &   p3*p6mup6*F(1337)+p1mup3*(p6mup6*F(1330)+p4mup6*F(1331)+p1mu
     &   p6*F(1332)+p2mup6*F(1333)+p5mup6*F(1338))+p1mup6*p6mup3*F(13
     &   40)+p1mup6*p4mup3*F(1341)-p1mup6*p2mup3*F(1342)+p2mup6*p6mup
     &   3*F(1344)+p2mup6*p4mup3*F(1345)-p2mup3*p2mup6*F(1346)-p4mup3
     &   *p5mup6*F(1347)-p4mup6*p6mup3*F(1348)-p4mup3*p4mup6*F(1349)+
     &   p2mup3*p5mup6*F(1350)+p2mup3*p4mup6*F(1351)-p5mup6*p6mup3*F(
     &   1352)-p5mup3*(p6mup6*F(1334)-p1mup6*F(1339)-p2mup6*F(1343)+p
     &   4mup6*F(1347)+p5mup6*F(1353))-mup3mup6*F(1354)
       Fa17(310) = mup3mup6*F(1355)-8*(p1mup3*p4mup6*F(1356)-p6mup3*(p
     &   6mup6*F(1071)+p2mup6*F(1087)+p4mup6*F(1360)+p5mup6*F(1365)-p
     &   1mup6*F(1369))+p1mup3*(p6mup6*F(1357)+p5mup6*F(1358)+p1mup6*
     &   F(1361)+p2mup6*F(1371))-p2mup3*(p5mup6*F(1046)+p6mup6*F(1073
     &   )-p2mup6*F(1081)-p1mup6*F(1366)-p4mup6*F(1372))-p5mup3*(p5mu
     &   p6*F(1043)+p4mup6*F(1359)+p6mup6*F(1362)-p1mup6*F(1368)+p2mu
     &   p6*F(1373))-p4mup3*(p4mup6*F(1363)+p6mup6*F(1364)+p5mup6*F(1
     &   367)-p1mup6*F(1370)+p2mup6*F(1374)))
       Fa17(311) = mup3mup6*F(1378)-8*(p2mup3*(p2mup6*F(1374)+p1mup6*F
     &   (1383))-p1mup3*(p6mup6*F(1380)+p5mup6*F(1381)-p2mup6*F(1382)
     &   -p1mup6*F(1385))-p5mup3*p6mup6*F(1386)-p6mup3*p6mup6*F(1387)
     &   -p4mup3*p6mup6*F(1389)-p5mup6*p6mup3*F(1390)-p4mup3*p5mup6*F
     &   (1391)-p5mup3*p5mup6*F(1392)-p1mup6*p5mup3*F(1393)-p1mup6*p6
     &   mup3*F(1394)-p1mup6*p4mup3*F(1395)-p2mup6*p5mup3*F(1396)-p2m
     &   up6*p6mup3*F(1397)-p2mup6*p4mup3*F(1398)-p2mup3*p5mup6*F(139
     &   9)-p4mup6*(p6mup3*F(1022)+p1mup3*F(1379)+p5mup3*F(1384)+p4mu
     &   p3*F(1388)+p2mup3*F(1400))-p2mup3*p6mup6*F(1401))
       Fa17(312) = p5mup3*F(1407)+2*(p1mup3*F(1405)+p4mup3*F(1406)+p6m
     &   up3*F(1408)+p2mup3*F(1409))
       Fa17(313) = p1mup3*F(1410)+p4mup3*F(1411)+p5mup3*F(1412)+2*p6mu
     &   p3*F(1413)+p2mup3*F(1414)
       Fa17(314) = -(p1mup3*F(1415))+p4mup3*F(1416)+p5mup3*F(1417)+p6m
     &   up3*F(1418)+p2mup3*F(1419)
       Fa17(315) = p6mup6*F(1424)-p1mup6*F(1425)-p2mup6*F(1426)+p4mup6
     &   *F(1427)+p5mup6*F(1428)
       Fa17(316) = 2*p2mup6*F(1319)+p4mup6*F(1429)-p5mup6*F(1430)-p6mu
     &   p6*F(1431)+p1mup6*F(1432)
       Fa17(317) = p1mup6*F(1435)+p6mup6*F(1436)+p2mup6*F(1437)+p4mup6
     &   *F(1438)+p5mup6*F(1439)
       Fa17(318) = p1mup6*F(1440)+p6mup6*F(1441)+p5mup6*F(1442)+p4mup6
     &   *F(1443)+p2mup6*F(1444)
       Fa17(319) = p1mup6*F(1445)+p5mup6*F(1446)+2*p6mup6*F(1447)+p4mu
     &   p6*F(1448)+p2mup6*F(1449)
       Fa17(320) = 2*(p1mup3*F(887)+p2mup3*F(932))-p4mup3*F(1452)-p5mu
     &   p3*F(1453)-p6mup3*F(1454)
       Fa17(321) = p1mup3*F(1455)-p4mup3*F(1456)-p5mup3*F(1457)-p6mup3
     &   *F(1458)-p2mup3*F(1459)
       Fa17(322) = p1mup3*F(1460)+p2mup3*F(1461)-p6mup3*F(1462)-p4mup3
     &   *F(1463)-p5mup3*F(1464)
       Fa17(323) = p6mup3*F(928)-p1mup3*F(1465)+p4mup3*F(1466)+p5mup3*
     &   F(1467)-p2mup3*F(1468)
       Return
       End
