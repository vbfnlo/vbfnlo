c************** Calling the Fa functions*************************
       subroutine HexBox2Fa19(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa19)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBox2FhlFunctions/F
       COMMON/HexBox2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,
     &   s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(385),F(1646)
      Real*8 P(1185) 
        COMMON/HexBox2PFunctions/P
       Complex*16 Fa19(343:361)
       COMMON/HexBox2FaFunctions/Fa
       Fa19(343) = p4mup6*F(1415)+p5mup6*F(1416)+p6mup6*F(1417)-p2mup6
     &   *F(1418)+p1mup6*F(1419)
       Fa19(344) = p1mup6*F(1414)+p4mup6*F(1420)+p5mup6*F(1421)+p6mup6
     &   *F(1422)-p2mup6*F(1423)
       Fa19(345) = p6mup3*F(1230)+p1mup3*F(1366)+p4mup3*F(1367)+p5mup3
     &   *F(1368)+p2mup3*F(1425)
       Fa19(346) = p4mup6*F(1430)-p5mup6*F(1431)+p6mup6*F(1432)-p2mup6
     &   *F(1433)+p1mup6*F(1434)
       Fa19(347) = -(p1mup6*F(1435))+p4mup6*F(1436)+p5mup6*F(1437)+p6m
     &   up6*F(1438)-p2mup6*F(1439)
       Fa19(348) = p2mup6*F(1440)-p4mup6*F(1441)+p5mup6*F(1442)-p6mup6
     &   *F(1443)+p1mup6*F(1444)
       Fa19(349) = p1mup3*F(1232)+p2mup3*F(1446)
       Fa19(350) = p4mup6*F(1265)+p2mup6*F(1391)-p5mup6*F(1448)+p6mup6
     &   *F(1449)+p1mup6*F(1450)
       Fa19(351) = -(p1mup3*F(1451))+2*(p4mup3*F(1270)+p5mup3*F(1271)+
     &   p2mup3*F(1452)+p6mup3*F(1644))
       Fa19(352) = p6mup3*F(1455)+p4mup3*F(1456)+p1mup3*F(1457)-p5mup3
     &   *F(1458)-p2mup3*F(1459)
       Fa19(353) = -(p1mup3*F(1460))+p6mup3*F(1461)+p4mup3*F(1462)+p5m
     &   up3*F(1463)-p2mup3*F(1464)
       Fa19(354) = -(p1mup3*F(1467))-p2mup3*F(1468)+p6mup3*F(1469)+p4m
     &   up3*F(1470)+p5mup3*F(1471)
       Fa19(355) = -(p1mup3*F(1212))+p6mup3*F(1474)+p4mup3*F(1475)+p5m
     &   up3*F(1476)-p2mup3*F(1477)
       Fa19(356) = p1mup3*F(1212)-p6mup3*F(1474)-p4mup3*F(1475)-p5mup3
     &   *F(1476)+p2mup3*F(1477)
       Fa19(357) = p4mup6*F(1481)+p5mup6*F(1482)+p6mup6*F(1483)+p2mup6
     &   *F(1484)+p1mup6*F(1485)
       Fa19(358) = -(p2mup6*F(1486))+p4mup6*F(1487)+p5mup6*F(1488)+p6m
     &   up6*F(1489)+p1mup6*F(1490)
       Fa19(359) = -(p1mup6*F(883))-p2mup6*F(888)+p5mup6*F(889)+p6mup6
     &   *F(890)+p4mup6*F(891)
       Fa19(360) = -(p1mup6*F(1491))-p2mup6*F(1492)+p4mup6*F(1493)+p5m
     &   up6*F(1494)+p6mup6*F(1495)
       Fa19(361) = p4mup6*F(1309)+p5mup6*F(1310)+p6mup6*F(1311)-p1mup6
     &   *F(1312)-p2mup6*F(1496)
       Return
       End
