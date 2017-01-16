c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa13(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa13)
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
       Complex*16 Fa13(289:312)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa13(289) = p3mup6*p5mup5*p6mup2*F(1356)+p4mup6*p5mup5*p6mup2*F
     -   (1356)-2*(-(p1mup5*(p3mup6+p4mup6)*p6mup2*F(1354))+p1mup6*p6
     -   mup2*p6mup5*F(1355)-p3mup6*p6mup2*p6mup5*F(1357)-p4mup6*p6mu
     -   p2*p6mup5*F(1357)+p5mup6*p6mup2*p6mup5*F(1358)+p2mup6*p6mup2
     -   *p6mup5*F(1359)-p3mup6*p4mup2*p4mup5*F(1360)-p4mup2*p4mup5*p
     -   4mup6*F(1360)+p3mup5*(p3mup6*p5mup2*F(1326)-p4mup2*p4mup6*F(
     -   1360))+p2mup5*(p3mup6+p4mup6)*p6mup2*F(9600))
       Fa13(290) = -(p1mup2*p3mup5*p3mup6*F(1307))-p2mup2*p3mup5*p3mup
     -   6*F(1310)+p3mup6*p4mup5*p5mup2*F(1326)+p3mup5*p4mup6*p5mup2*
     -   F(1326)+p4mup5*p4mup6*p5mup2*F(1326)-p3mup2*p3mup5*p3mup6*F(
     -   1360)-p3mup5*p3mup6*p4mup2*F(1360)-p3mup2*p3mup6*p4mup5*F(13
     -   60)-p3mup2*p3mup5*p4mup6*F(1360)-p3mup2*p4mup5*p4mup6*F(1360
     -   )+(p3mup5+p4mup5)*(p3mup6+p4mup6)*p6mup2*F(9592)
       Fa13(291) = mup5mup6*(p5mup2*F(1362)+p2mup2*F(1363)+(p3mup2+p4m
     -   up2)*F(1364))
       Fa13(292) = -8*p2mup2*(p2mup6*(p5mup5*F(1366)+2*(p6mup5*F(1367)
     -   -p2mup5*F(1374)-p1mup5*F(1385)))-2*p1mup5*(p5mup6*F(1368)+p1
     -   mup6*F(1386)))+(p3mup2+p4mup2)*(-4*(p1mup5*p6mup6*F(1376)+p2
     -   mup5*p6mup6*F(1380)-p5mup5*(p5mup6*F(1371)+2*(p1mup6*F(1377)
     -   +p2mup6*F(1381))+p6mup6*F(1383)))+p6mup5*p6mup6*F(1384)+16*(
     -   p5mup6*(p6mup5*F(1372)-p1mup5*F(1375)-p2mup5*F(1379))+p6mup5
     -   *(p1mup6*F(1378)+p2mup6*F(1382))-p1mup5*p1mup6*F(1387)-p1mup
     -   5*p2mup6*F(1388)-p1mup6*p2mup5*F(1389)-p2mup5*p2mup6*F(1390)
     -   ))+4*p5mup2*p5mup6*(p5mup5*F(1369)-4*(p1mup5*F(1370)+p2mup5*
     -   F(1373)+p6mup5*F(9593)))
       Fa13(293) = mup2mup6*(p5mup5*F(1392)+p6mup5*F(1393)+2*(p1mup5*F
     -   (1391)-p2mup5*F(1394)))+2*mup5mup6*p6mup2*F(1395)+p5mup2*p6m
     -   up5*p6mup6*F(1410)-8*(p1mup6*p2mup2*p5mup5*F(1399)-p2mup6*p5
     -   mup2*p5mup5*F(1411))-4*(p1mup6*p5mup2*p5mup5*F(1397)-p2mup2*
     -   p2mup5*p6mup6*F(1402)-p2mup2*p5mup5*p5mup6*F(1403)-p2mup2*p5
     -   mup5*p6mup6*F(1404)+p2mup2*p6mup5*p6mup6*F(1405)-p5mup2*p5mu
     -   p5*p6mup6*F(1408)+p2mup5*p5mup2*p6mup6*F(1413)-p1mup5*(p2mup
     -   2*p6mup6*F(1398)+p5mup2*(p6mup6*F(1396)-4*(p2mup6*F(1406)+p1
     -   mup6*F(1415)))))-16*(p1mup6*p2mup2*p6mup5*F(1400)-p2mup2*p2m
     -   up5*p5mup6*F(1401)-p1mup6*p2mup5*p5mup2*F(1407)-p1mup6*p5mup
     -   2*p6mup5*F(1409)-p2mup6*p5mup2*p6mup5*F(1412)-p1mup6*p2mup2*
     -   p2mup5*F(1414)-p2mup5*p2mup6*p5mup2*F(1416)+p2mup2*p5mup6*p6
     -   mup5*F(9589))
       Fa13(294) = mup2mup6*p4mup5*F(1417)+2*(p1mup2*(-(p1mup5*p6mup6*
     -   F(1420))-p2mup5*p6mup6*F(1425)+p5mup5*(p6mup6*F(1418)-2*(p1m
     -   up6*F(1421)+p2mup6*F(1426))+p5mup6*F(1428))-2*p4mup5*p6mup6*
     -   F(1429)+p6mup5*(p6mup6*F(1419)-4*(p1mup6*F(1422)+p2mup6*F(14
     -   27)+p5mup6*F(1430)))+4*(p1mup5*p5mup6*F(1423)+p2mup5*p5mup6*
     -   F(1424)+p4mup5*p5mup6*F(1431)+p1mup6*p2mup5*F(1432)+p2mup5*p
     -   2mup6*F(1433)+p4mup5*p4mup6*F(1434)+p1mup5*p1mup6*F(1437)+p1
     -   mup5*p2mup6*F(1438)))+2*p2mup2*p4mup5*(p6mup6*F(1436)-2*(p5m
     -   up6*F(1435)-p4mup6*F(1439))))
       Fa13(295) = mup2mup6*p3mup5*F(1417)+mup2mup5*p6mup6*F(1441)+2*(
     -   mup2mup5*p1mup6*F(1440)-mup2mup5*p2mup6*F(1442)-mup2mup5*p5m
     -   up6*F(1443))-4*(p1mup2*p3mup5*(p6mup6*F(1429)-2*p5mup6*F(143
     -   1))-p2mup2*p3mup5*p6mup6*F(1436)+p3mup5*p5mup2*p6mup6*F(1446
     -   )+p4mup5*p5mup2*p6mup6*F(1446)-p3mup2*p3mup5*p6mup6*F(1454)-
     -   p3mup5*p4mup2*p6mup6*F(1454)-p3mup2*p4mup5*p6mup6*F(1454)-p4
     -   mup2*p4mup5*p6mup6*F(1454))+8*(-(p2mup2*p3mup5*p5mup6*F(1435
     -   ))+p1mup6*p3mup5*p5mup2*F(1444)+p1mup6*p4mup5*p5mup2*F(1444)
     -   +p3mup5*p5mup2*p5mup6*F(1445)+p4mup5*p5mup2*p5mup6*F(1445)-p
     -   1mup6*p2mup2*p3mup5*F(1447)-p1mup6*p2mup2*p4mup5*F(1447)-p2m
     -   up2*p2mup6*p3mup5*F(1448)-p2mup2*p2mup6*p4mup5*F(1448)-p1mup
     -   2*p1mup6*(p3mup5+p4mup5)*F(1449)-p1mup2*p2mup6*(p3mup5+p4mup
     -   5)*F(1450)-p1mup6*p3mup2*p3mup5*F(1451)-p1mup6*p3mup5*p4mup2
     -   *F(1451)-p1mup6*p3mup2*p4mup5*F(1451)-p1mup6*p4mup2*p4mup5*F
     -   (1451)+p3mup2*p3mup5*p5mup6*F(1452)+p3mup5*p4mup2*p5mup6*F(1
     -   452)+p3mup2*p4mup5*p5mup6*F(1452)+p4mup2*p4mup5*p5mup6*F(145
     -   2)+p2mup6*p3mup2*p3mup5*F(1453)+p2mup6*p3mup5*p4mup2*F(1453)
     -   +p2mup6*p3mup2*p4mup5*F(1453)+p2mup6*p4mup2*p4mup5*F(1453)-p
     -   2mup6*p3mup5*p5mup2*F(1455)-p2mup6*p4mup5*p5mup2*F(1455)-p2m
     -   up5*p5mup6*p6mup2*F(9586))
       Fa13(296) = p4mup6*(p1mup2*(-(p5mup5*F(1456))+2*(p3mup5*F(1434)
     -   +p1mup5*F(1458)+p2mup5*F(1459)))+p2mup2*(p5mup5*F(1457)+2*(p
     -   3mup5*F(1439)+p1mup5*F(1460)+p2mup5*F(1461))))
       Fa13(297) = p6mup2*(-4*(p3mup5+p4mup5)*p5mup6*F(9582)+p5mup5*(p
     -   5mup6*F(9578)+4*p2mup6*F(9586)))
       Fa13(298) = mup2mup5*p3mup6*F(1462)+mup2mup5*p4mup6*F(1462)+2*(
     -   -(p1mup2*p3mup6*(-(p5mup5*F(1456))+2*(p4mup5*F(1434)+p1mup5*
     -   F(1458)+p2mup5*F(1459))))-p2mup2*p3mup6*(p5mup5*F(1457)+2*(p
     -   4mup5*F(1439)+p1mup5*F(1460)+p2mup5*F(1461)))+p3mup6*p5mup2*
     -   p5mup5*F(1465)+p4mup6*p5mup2*p5mup5*F(1465)+p3mup2*p3mup6*p5
     -   mup5*F(1468)+p3mup6*p4mup2*p5mup5*F(1468)+p3mup2*p4mup6*p5mu
     -   p5*F(1468)+p4mup2*p4mup6*p5mup5*F(1468)+p5mup5*p6mup2*p6mup6
     -   *F(9578))+4*(p1mup5*p3mup6*p5mup2*F(1463)+p1mup5*p4mup6*p5mu
     -   p2*F(1463)+p1mup2*(p3mup6+p4mup6)*p6mup5*F(1464)+p3mup6*p5mu
     -   p2*p6mup5*F(1466)+p4mup6*p5mup2*p6mup5*F(1466)+p2mup2*(p3mup
     -   6+p4mup6)*p6mup5*F(1467)-p2mup5*p3mup6*p5mup2*F(1469)-p2mup5
     -   *p4mup6*p5mup2*F(1469)-p1mup6*p5mup5*p6mup2*F(9585)+p1mup5*p
     -   5mup6*p6mup2*F(9585)+p3mup2*p3mup6*p6mup5*F(9591)+p3mup6*p4m
     -   up2*p6mup5*F(9591)+p3mup2*p4mup6*p6mup5*F(9591)+p4mup2*p4mup
     -   6*p6mup5*F(9591)-p2mup5*p3mup2*p3mup6*F(9610)-p2mup5*p3mup6*
     -   p4mup2*F(9610)-p2mup5*p3mup2*p4mup6*F(9610)-p2mup5*p4mup2*p4
     -   mup6*F(9610)-p1mup5*(p3mup2+p4mup2)*(p3mup6+p4mup6)*F(9612))
       Fa13(299) = (p3mup6+p4mup6)*(p5mup5*p6mup2*F(9582)+p4mup2*p4mup
     -   5*F(9607))+p3mup5*(p3mup6*p5mup2*F(1470)+p4mup2*p4mup6*F(960
     -   7))
       Fa13(300) = p1mup2*p3mup5*p3mup6*F(1434)+p2mup2*p3mup5*p3mup6*F
     -   (1439)+(p3mup6*p4mup5+(p3mup5+p4mup5)*p4mup6)*p5mup2*F(1470)
     -   +(p3mup5*p3mup6*p4mup2+p3mup2*(p3mup5+p4mup5)*(p3mup6+p4mup6
     -   ))*F(9607)
       Fa13(301) = mup3mup6*(p2mup2*F(1471)+p5mup2*F(1472)+4*(p1mup2*F
     -   (1473)+p6mup2*F(1474)-(p3mup2+p4mup2)*F(1475)))
       Fa13(302) = mup3mup6*(-(p2mup2*F(1476))+p5mup2*F(1477)+p1mup2*F
     -   (1478)+p6mup2*F(1479)+(p3mup2+p4mup2)*F(1480))
       Fa13(303) = p5mup2*p5mup5*F(1482)+p5mup2*p6mup5*F(1483)+p1mup2*
     -   p5mup5*F(1485)+2*(p4mup2*p5mup5*F(933)-mup2mup5*F(939)+p3mup
     -   5*p4mup2*F(940)+p4mup2*p4mup5*F(940)+p1mup2*p6mup5*F(1486))+
     -   p5mup5*p6mup2*F(1488)+p6mup2*p6mup5*F(1489)+p1mup5*(p5mup2*F
     -   (1481)+p1mup2*F(1484)+p6mup2*F(1487)+2*((p3mup2+p4mup2)*F(93
     -   1)+p2mup2*F(1490)))+p2mup2*p2mup5*F(1491)+p2mup2*p5mup5*F(14
     -   92)+p2mup2*p6mup5*F(1493)+p2mup5*p5mup2*F(1494)+p1mup2*p2mup
     -   5*F(1495)+p2mup5*p6mup2*F(1496)+p3mup5*p5mup2*F(1497)+p4mup5
     -   *p5mup2*F(1497)+p3mup5*p6mup2*F(1498)+p4mup5*p6mup2*F(1498)+
     -   p1mup2*p3mup5*F(1499)+p1mup2*p4mup5*F(1499)+p2mup5*p4mup2*F(
     -   1500)+p4mup2*p6mup5*F(1501)+p3mup2*(2*(p5mup5*F(933)+(p3mup5
     -   +p4mup5)*F(940))+p2mup5*F(1500)+p6mup5*F(1501))+p2mup2*(p3mu
     -   p5+p4mup5)*F(1502)
       Fa13(304) = mup3mup5*(-(p1mup2*F(1473))-p6mup2*F(1474)+(p3mup2+
     -   p4mup2)*F(1475)+p2mup2*F(1503)+p5mup2*F(1504))
       Fa13(305) = mup3mup5*(-(p1mup2*F(1478))+p2mup2*F(1505)+p5mup2*F
     -   (1506)+p6mup2*F(1507)+(p3mup2+p4mup2)*F(1508))
       Fa13(306) = 2*(p1mup2*F(1473)+p6mup2*F(1474)-(p3mup2+p4mup2)*F(
     -   1475))+p2mup2*F(1509)+p5mup2*F(1510)
       Fa13(307) = -(p2mup2*F(1476))+p5mup2*F(1477)+p1mup2*F(1478)+p6m
     -   up2*F(1479)+(p3mup2+p4mup2)*F(1480)
       Fa13(308) = p2mup2*p6mup5*F(1513)+p1mup2*p5mup5*F(1514)+p5mup5*
     -   p6mup2*F(1517)+p2mup5*p6mup2*F(1519)+p5mup2*p5mup5*F(1521)+p
     -   5mup2*p6mup5*F(1522)+p1mup5*(2*((p3mup2+p4mup2)*F(986)+p2mup
     -   2*F(1511))+p6mup2*F(1516)+p5mup2*F(1520)+p1mup2*F(1523))+p3m
     -   up5*p6mup2*F(1524)+p4mup5*p6mup2*F(1524)+p6mup2*p6mup5*F(152
     -   6)+p2mup5*p4mup2*F(1527)+2*(p3mup5*p4mup2*F(998)+p4mup2*p4mu
     -   p5*F(998)-p2mup2*p5mup5*F(1512)+p1mup2*p6mup5*F(1515)+p2mup2
     -   *p2mup5*F(1518)+p2mup2*p3mup5*F(1525)+p2mup2*p4mup5*F(1525)+
     -   p4mup2*p5mup5*F(1528))+p4mup2*p6mup5*F(1529)+p3mup2*(p2mup5*
     -   F(1527)+2*((p3mup5+p4mup5)*F(998)+p5mup5*F(1528))+p6mup5*F(1
     -   529))+p3mup5*p5mup2*F(1530)+p4mup5*p5mup2*F(1530)+p1mup2*p3m
     -   up5*F(1531)+p1mup2*p4mup5*F(1531)+p2mup5*p5mup2*F(1532)+p1mu
     -   p2*p2mup5*F(1533)+mup2mup5*F(1534)
       Fa13(309) = p2mup2*p6mup5*F(1537)+p5mup2*p5mup5*F(1539)-p5mup2*
     -   p6mup5*F(1540)+p1mup2*p5mup5*F(1542)-2*(p2mup2*p5mup5*F(1536
     -   )+p1mup2*p6mup5*F(1543))+p1mup5*(2*p2mup2*F(1535)+p5mup2*F(1
     -   538)+p1mup2*F(1541)+p6mup2*F(1544))+p5mup5*p6mup2*F(1545)+p6
     -   mup2*p6mup5*F(1546)+p2mup5*(2*((p3mup2+p4mup2)*F(1041)+p2mup
     -   2*F(1547))+p5mup2*F(1548)+p1mup2*F(1549)+p6mup2*F(1550))+(p3
     -   mup2+p4mup2)*(p1mup5*F(1555)+p5mup5*F(1556)+p6mup5*F(1557))+
     -   mup2mup5*F(1558)+(p3mup5+p4mup5)*(p5mup2*F(1551)+p6mup2*F(15
     -   52)+p1mup2*F(1553)+2*p2mup2*F(1554)+(p3mup2+p4mup2)*F(1559))
       Fa13(310) = mup5mup6*(-(p2mup2*F(1560))-p5mup2*F(1561)+(p3mup2+
     -   p4mup2)*F(1562))
       Fa13(311) = p2mup2*p2mup6*p5mup5*F(1565)+2*p5mup2*p5mup5*p5mup6
     -   *F(1579)+4*(p1mup5*p2mup2*p2mup6*F(1564)+p2mup2*p2mup6*p6mup
     -   5*F(1566)-p1mup5*p2mup2*p5mup6*F(1567)-p1mup5*p5mup2*p5mup6*
     -   F(1577)-p2mup5*p5mup2*p5mup6*F(1578)+p1mup5*p1mup6*p2mup2*F(
     -   1581)+p2mup2*p2mup5*p2mup6*F(1582)+p5mup2*p5mup6*p6mup5*F(15
     -   86))+p3mup2*(p2mup5*p6mup6*F(1573)-p5mup5*(p1mup6*F(1570)+p2
     -   mup6*F(1574))+2*p5mup5*p5mup6*F(1576)-p5mup5*p6mup6*F(1580)-
     -   p6mup5*(p6mup6*F(1357)+4*(p1mup6*F(1571)+p2mup6*F(1575)-p5mu
     -   p6*F(1585)))+p1mup5*(p6mup6*F(1569)+4*(p5mup6*F(1568)-p2mup6
     -   *F(1583)+p1mup6*F(1587)))+4*(p2mup5*p5mup6*F(1572)+p1mup6*p2
     -   mup5*F(1584)+p2mup5*p2mup6*F(1588)))+p4mup2*(p2mup5*p6mup6*F
     -   (1573)-p5mup5*(p1mup6*F(1570)+p2mup6*F(1574))+2*p5mup5*p5mup
     -   6*F(1576)-p5mup5*p6mup6*F(1580)-p6mup5*(p6mup6*F(1357)+4*(p1
     -   mup6*F(1571)+p2mup6*F(1575)-p5mup6*F(1585)))+p1mup5*(p6mup6*
     -   F(1569)+4*(p5mup6*F(1568)-p2mup6*F(1583)+p1mup6*F(1587)))+4*
     -   (p2mup5*p5mup6*F(1572)+p1mup6*p2mup5*F(1584)+p2mup5*p2mup6*F
     -   (1588)))
       Fa13(312) = mup5mup6*p6mup2*F(1591)+mup2mup6*(p1mup5*F(1589)+p2
     -   mup5*F(1592)-2*(p5mup5*F(1590)+p6mup5*F(1593)))-4*(p2mup2*(p
     -   1mup5*p6mup6*F(1594)-p1mup6*p5mup5*F(1595)+p2mup5*p6mup6*F(1
     -   597)+2*p5mup5*p5mup6*F(1600)+p5mup5*p6mup6*F(1601)-p6mup5*p6
     -   mup6*F(1602)+4*(p1mup6*p2mup5*F(1596)-p1mup6*p6mup5*F(1598)+
     -   p2mup5*p5mup6*F(1599)-p5mup6*p6mup5*F(1607)))-p5mup2*(p1mup5
     -   *p6mup6*F(1604)-p5mup5*(p2mup6*F(1603)+p1mup6*F(1605))+p2mup
     -   5*p6mup6*F(1606)+p6mup5*(p6mup6*F(1358)+4*(p2mup6*F(1608)+p1
     -   mup6*F(1609)))-p5mup5*p6mup6*F(1610)+4*(p1mup5*p2mup6*F(1611
     -   )-p1mup6*p2mup5*F(1612)+p1mup5*p1mup6*F(1613)+p2mup5*p2mup6*
     -   F(1614))))
       Return
       End
