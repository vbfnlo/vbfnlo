c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g27(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa27)
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
       Complex*16 Fa27(365:378)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa27(365) = 4*p4mup2*(mup4mup6*F(1498)+2*(p2mup4*p2mup6*F(1499)
     -   +p2mup4*p4mup6*F(1500)+p2mup4*p5mup6*F(1501)+p2mup4*p3mup6*F
     -   (1502)+p1mup4*((p1mup6+p2mup6)*F(1499)+p4mup6*F(1500)+p5mup6
     -   *F(1501)+p3mup6*F(1502))-p2mup6*p3mup4*F(1503)-p2mup6*p4mup4
     -   *F(1504)-p3mup4*p4mup6*F(1505)-p4mup4*p4mup6*F(1506)-p3mup4*
     -   p5mup6*F(1507)-p4mup4*p5mup6*F(1508)-p3mup4*p3mup6*F(1509)-p
     -   3mup6*p4mup4*F(1510)-p3mup6*p5mup4*F(1511)-p2mup6*p5mup4*F(1
     -   512)+p1mup6*(p2mup4*F(1499)-p3mup4*F(1503)-p4mup4*F(1504)-p5
     -   mup4*F(1512))-p4mup6*p5mup4*F(1513)-p5mup4*p5mup6*F(1514)))
       Fa27(366) = -2*mup4mup6*F(1498)-4*(p2mup4*p2mup6*F(1499)+p2mup4
     -   *p4mup6*F(1500)+p2mup4*p5mup6*F(1501)+p2mup4*p3mup6*F(1502)+
     -   p1mup4*((p1mup6+p2mup6)*F(1499)+p4mup6*F(1500)+p5mup6*F(1501
     -   )+p3mup6*F(1502))-p2mup6*p4mup4*F(1504)-p3mup4*p4mup6*F(1505
     -   )-p4mup4*p4mup6*F(1506)-p3mup4*p5mup6*F(1507)-p4mup4*p5mup6*
     -   F(1508)-p3mup4*p3mup6*F(1509)-p3mup6*p4mup4*F(1510)-p3mup6*p
     -   5mup4*F(1511))+p2mup6*p3mup4*F(1533)+p5mup4*p5mup6*F(1534)+p
     -   2mup6*p5mup4*F(1535)+p1mup6*(-4*(p2mup4*F(1499)-p4mup4*F(150
     -   4))+p3mup4*F(1533)+p5mup4*F(1535))+p4mup6*p5mup4*F(1536)
       Fa27(367) = p4mup6*F(1537)+2*(p1mup6+p2mup6)*F(1538)+p5mup6*F(1
     -   539)+p3mup6*F(1540)
       Fa27(368) = 2*(p4mup4*F(1525)+(p1mup4+p2mup4)*F(1528))+p3mup4*F
     -   (1542)+p5mup4*F(1543)
       Fa27(369) = 4*p3mup2*s12*(p4mup4*F(1529)+p3mup4*F(1530)+p5mup4*
     -   F(1531)+(p1mup4+p2mup4)*F(1532))+(p1mup2+p2mup2)*(p4mup4*F(1
     -   544)+p5mup4*F(1545)+p3mup4*F(1546)-(p1mup4+p2mup4)*F(1547))
       Fa27(370) = p5mup4*F(1548)+p4mup4*F(1549)-(p1mup4+p2mup4)*F(155
     -   0)+p3mup4*F(1551)
       Fa27(371) = p5mup3*p5mup4*F(1422)+p1mup3*p5mup4*F(1424)+p2mup3*
     -   p5mup4*F(1424)+p3mup4*p5mup3*F(1437)+p1mup3*p4mup4*F(1440)+p
     -   2mup3*p4mup4*F(1440)+p1mup4*p5mup3*F(1443)+p2mup4*p5mup3*F(1
     -   443)+p4mup4*p5mup3*F(1445)+p1mup3*p1mup4*F(1447)+p1mup4*p2mu
     -   p3*F(1447)+p1mup3*p2mup4*F(1447)+p2mup3*p2mup4*F(1447)+(p1mu
     -   p3+p2mup3)*p3mup4*F(1449)+p3mup3*(p5mup4*F(1421)+p3mup4*F(14
     -   44)+(p1mup4+p2mup4)*F(1450)+p4mup4*F(1451))+mup3mup4*F(1462)
     -   -p4mup3*(p5mup4*F(1417)+(p1mup4+p2mup4)*F(1453)-p4mup4*F(145
     -   6)-p3mup4*F(1552))
       Fa27(372) = p3mup3*p5mup4*F(1418)+p4mup3*p5mup4*F(1419)+p5mup3*
     -   p5mup4*F(1423)+p3mup3*p4mup4*F(1439)+p1mup4*p3mup3*F(1442)+p
     -   2mup4*p3mup3*F(1442)+p1mup3*(p5mup4*F(1416)+p4mup4*F(1438)+p
     -   3mup4*F(1441)+(p1mup4+p2mup4)*F(1448))+p2mup3*(p5mup4*F(1416
     -   )+p4mup4*F(1438)+p3mup4*F(1441)+(p1mup4+p2mup4)*F(1448))+p3m
     -   up3*p3mup4*F(1452)+p1mup4*p4mup3*F(1455)+p2mup4*p4mup3*F(145
     -   5)+p3mup4*p4mup3*F(1457)+p4mup3*p4mup4*F(1458)+p4mup4*p5mup3
     -   *F(1459)+p1mup4*p5mup3*F(1460)+p2mup4*p5mup3*F(1460)+p3mup4*
     -   p5mup3*F(1461)+mup3mup4*F(1463)
       Fa27(373) = 4*p3mup2*(-(mup3mup6*F(1553))+2*(-(p2mup3*p4mup6*F(
     -   1554))-p1mup6*p2mup3*F(1555)-p2mup3*p2mup6*F(1555)-p2mup3*p3
     -   mup6*F(1556)-p2mup3*p5mup6*F(1557)-p1mup3*(p4mup6*F(1554)+(p
     -   1mup6+p2mup6)*F(1555)+p3mup6*F(1556)+p5mup6*F(1557))+p3mup3*
     -   p3mup6*F(1558)+p1mup6*p3mup3*F(1559)+p2mup6*p3mup3*F(1559)+p
     -   3mup3*p4mup6*F(1560)+p3mup3*p5mup6*F(1561)+p4mup3*(p4mup6*F(
     -   1505)+p3mup6*F(1509)+(p1mup6+p2mup6)*F(1562)+p5mup6*F(1563))
     -   +p3mup6*p5mup3*F(1564)+p1mup6*p5mup3*F(1565)+p2mup6*p5mup3*F
     -   (1565)+p4mup6*p5mup3*F(1566)+p5mup3*p5mup6*F(1567)))
       Fa27(374) = 4*(p1mup2+p2mup2)*(-(mup3mup6*F(1553))+2*(-(p2mup3*
     -   p4mup6*F(1554))-p1mup6*p2mup3*F(1555)-p2mup3*p2mup6*F(1555)-
     -   p2mup3*p3mup6*F(1556)-p2mup3*p5mup6*F(1557)-p1mup3*(p4mup6*F
     -   (1554)+(p1mup6+p2mup6)*F(1555)+p3mup6*F(1556)+p5mup6*F(1557)
     -   )+p3mup3*p3mup6*F(1558)+p1mup6*p3mup3*F(1559)+p2mup6*p3mup3*
     -   F(1559)+p3mup3*p4mup6*F(1560)+p3mup3*p5mup6*F(1561)+p4mup3*(
     -   p4mup6*F(1505)+p3mup6*F(1509)+(p1mup6+p2mup6)*F(1562)+p5mup6
     -   *F(1563))+p3mup6*p5mup3*F(1564)+p1mup6*p5mup3*F(1565)+p2mup6
     -   *p5mup3*F(1565)+p4mup6*p5mup3*F(1566)+p5mup3*p5mup6*F(1567))
     -   )
       Fa27(375) = mup3mup4*(p3mup6*F(1073)+(p1mup6+p2mup6)*F(1074)+p4
     -   mup6*F(1083)+p5mup6*F(1087))+2*(mup4mup6*((p1mup3+p2mup3)*F(
     -   1076)-p3mup3*F(1077)+p5mup3*F(1078)+p4mup3*F(1376))-mup3mup6
     -   *((p1mup4+p2mup4)*F(1089)-p5mup4*F(1090)-p4mup4*F(1091)+p3mu
     -   p4*F(1378)))
       Fa27(376) = 2*(p4mup3*F(1389)+p3mup3*F(1390)+p5mup3*F(1391))+(p
     -   1mup3+p2mup3)*F(1568)
       Fa27(377) = 4*p4mup2*(-(mup3mup6*F(1553))+2*(-(p2mup3*p4mup6*F(
     -   1554))-p1mup6*p2mup3*F(1555)-p2mup3*p2mup6*F(1555)-p2mup3*p3
     -   mup6*F(1556)-p2mup3*p5mup6*F(1557)-p1mup3*(p4mup6*F(1554)+(p
     -   1mup6+p2mup6)*F(1555)+p3mup6*F(1556)+p5mup6*F(1557))+p3mup3*
     -   p3mup6*F(1558)+p1mup6*p3mup3*F(1559)+p2mup6*p3mup3*F(1559)+p
     -   3mup3*p4mup6*F(1560)+p3mup3*p5mup6*F(1561)+p4mup3*(p4mup6*F(
     -   1505)+p3mup6*F(1509)+(p1mup6+p2mup6)*F(1562)+p5mup6*F(1563))
     -   +p3mup6*p5mup3*F(1564)+p1mup6*p5mup3*F(1565)+p2mup6*p5mup3*F
     -   (1565)+p4mup6*p5mup3*F(1566)+p5mup3*p5mup6*F(1567)))
       Fa27(378) = -2*mup3mup6*F(1553)+4*(-(p2mup3*p4mup6*F(1554))-p1m
     -   up6*p2mup3*F(1555)-p2mup3*p2mup6*F(1555)-p2mup3*p3mup6*F(155
     -   6)-p2mup3*p5mup6*F(1557)-p1mup3*(p4mup6*F(1554)+(p1mup6+p2mu
     -   p6)*F(1555)+p3mup6*F(1556)+p5mup6*F(1557))+p3mup3*p3mup6*F(1
     -   558)+p1mup6*p3mup3*F(1559)+p2mup6*p3mup3*F(1559)+p3mup3*p4mu
     -   p6*F(1560)+p3mup3*p5mup6*F(1561)+p4mup3*(p4mup6*F(1505)+p3mu
     -   p6*F(1509)+(p1mup6+p2mup6)*F(1562)+p5mup6*F(1563))+p3mup6*p5
     -   mup3*F(1564)+p1mup6*p5mup3*F(1565)+p2mup6*p5mup3*F(1565)+p4m
     -   up6*p5mup3*F(1566)+p5mup3*p5mup6*F(1567))
       Return
       End
