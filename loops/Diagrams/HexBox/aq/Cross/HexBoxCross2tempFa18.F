c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa18(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   18)
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
       Complex*16 Fa18(324:342)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa18(324) = p6mup6*F(1469)+p1mup6*F(1470)+p4mup6*F(1471)+p5mup6
     &   *F(1472)+p2mup6*F(1473)
       Fa18(325) = p5mup6*F(1474)+p6mup6*F(1475)+p1mup6*F(1476)+p4mup6
     &   *F(1477)+p2mup6*F(1478)
       Fa18(326) = p5mup6*F(1479)+p6mup6*F(1480)+p2mup6*F(1481)+p1mup6
     &   *F(1482)+p4mup6*F(1483)
       Fa18(327) = -(p6mup6*F(1484))+p2mup6*F(1485)+p1mup6*F(1486)+p4m
     &   up6*F(1487)+p5mup6*F(1488)
       Fa18(328) = -(p5mup6*F(1489))-p6mup6*F(1490)+p2mup6*F(1491)+p1m
     &   up6*F(1492)+p4mup6*F(1493)
       Fa18(329) = p6mup6*F(1494)-2*(p1mup6*F(1003)+p2mup6*F(1005)-p4m
     &   up6*F(1006)-p5mup6*F(1495))
       Fa18(330) = p2mup6*F(1496)-p4mup6*F(1497)+p5mup6*F(1498)+p6mup6
     &   *F(1499)+p1mup6*F(1500)
       Fa18(331) = p1mup6*F(898)+p2mup6*F(1501)-p4mup6*F(1502)+p5mup6*
     &   F(1503)+p6mup6*F(1504)
       Fa18(332) = p2mup6*F(1505)-p1mup6*F(1506)-p4mup6*F(1507)+p5mup6
     &   *F(1508)+p6mup6*F(1509)
       Fa18(333) = p4mup6*F(1514)+p5mup6*F(1515)+p6mup6*F(1516)+p1mup6
     &   *F(1517)-p2mup6*F(1518)
       Fa18(334) = -(p1mup3*F(1460))-p2mup3*F(1461)+p4mup3*F(1463)+p5m
     &   up3*F(1464)+p6mup3*F(1519)
       Fa18(335) = p1mup6*F(929)+p4mup6*F(930)+p2mup6*F(931)+p5mup6*F(
     &   933)+p6mup6*F(1521)
       Fa18(336) = p1mup6*F(899)+p4mup6*F(1522)+p5mup6*F(1523)+p6mup6*
     &   F(1524)-p2mup6*F(1525)
       Fa18(337) = p6mup6*F(1526)-p1mup6*F(1527)-p2mup6*F(1528)+p4mup6
     &   *F(1529)+p5mup6*F(1530)
       Fa18(338) = p1mup6*F(1531)+p4mup6*F(1532)+p5mup6*F(1533)+p6mup6
     &   *F(1534)-p2mup6*F(1535)
       Fa18(339) = p5mup6*F(958)+p6mup6*F(1536)-p1mup6*F(1537)-p2mup6*
     &   F(1538)+p4mup6*F(1539)
       Fa18(340) = p1mup6*F(986)+p4mup6*F(992)-p2mup6*F(1004)+p5mup6*F
     &   (1540)+p6mup6*F(1541)
       Fa18(341) = p2mup3*F(1545)+p5mup3*F(1546)+p6mup3*F(1547)+p1mup3
     &   *F(1548)+p4mup3*F(1549)
       Fa18(342) = 2*(p5mup3*F(949)+p4mup3*F(950))+p1mup3*F(1552)-p2mu
     &   p3*F(1553)+p6mup3*F(1554)
       Return
       End
