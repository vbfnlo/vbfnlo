c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe25(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa25)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup5, p2mup2, p2mup3, 
     -          p2mup4, p2mup5, p3mup2, p3mup3, p3mup4, p3mup5, 
     -          p4mup2, p4mup3, p4mup4, p4mup5, p5mup2, p5mup3, 
     -          p5mup4, p5mup5, p6mup2, p6mup3, p6mup4, p6mup5
       Complex*16   mup2mup3, mup2mup4, mup2mup5, mup3mup4, mup3mup5, 
     -          mup4mup5
        common/test_FhlFunctionsAbe/F
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
      Complex*16 Fa(445),F(5930)
      Real*8 P(11674) 
        COMMON/test_PFunctionsAbe/P
       Complex*16 Fa25(265:275)
       COMMON/test_FaFunctionsAbe/Fa
       Fa25(265) = -16*(p1mup2*p3mup4*p4mup3*F(1464)+p2mup2*p3mup4*p4m
     -   up3*F(1465)-p1mup3*p3mup2*p3mup4*F(1466)-p2mup3*p3mup2*p3mup
     -   4*F(1467)-p1mup2*p2mup4*p4mup3*F(1468)+p1mup4*p2mup2*p4mup3*
     -   F(1469))+p1mup2*p4mup3*p4mup4*F(1470)+p2mup2*p4mup3*p4mup4*F
     -   (1471)+mup2mup4*(p2mup3*F(1457)+8*(p1mup3*F(1456)+p4mup3*F(1
     -   458)-p3mup3*F(1461)-p6mup3*F(1474)))-8*(-(mup2mup3*p6mup4*F(
     -   1462))+mup3mup4*(p1mup2*F(1459)-p2mup2*F(1460)-p6mup2*F(1463
     -   ))-2*(-(p1mup3*p4mup2*p6mup4*F(1472))+p4mup2*p6mup3*p6mup4*F
     -   (1473)+p1mup2*p2mup4*p6mup3*F(1475)+p1mup2*p3mup4*p6mup3*F(1
     -   476)+p1mup2*p4mup4*p6mup3*F(1477)+p1mup4*p2mup2*p6mup3*F(147
     -   8)+p2mup2*p4mup4*p6mup3*F(1479)-p4mup3*p6mup2*p6mup4*F(1480)
     -   -p2mup2*p2mup4*p4mup3*F(1481)+p1mup4*p3mup2*p3mup3*F(1482)+p
     -   2mup4*p3mup2*p3mup3*F(1483)+p3mup2*p3mup3*p3mup4*F(1484)+p1m
     -   up2*p6mup3*p6mup4*F(1485)+p1mup3*p6mup2*p6mup4*F(1486)-p6mup
     -   2*p6mup3*p6mup4*F(1487)+p2mup2*p2mup4*p6mup3*F(1488)+p2mup2*
     -   p3mup4*p6mup3*F(1489)-p1mup2*p1mup4*p4mup3*F(1490)+p1mup2*p1
     -   mup4*p6mup3*F(1491)))
       Fa25(266) = 8*(mup3mup4*p4mup2*F(1494)-mup2mup3*(p2mup4*F(1492)
     -   -p3mup4*F(1493)-p4mup4*F(1495)+p1mup4*F(1496)))+p2mup4*p4mup
     -   2*p4mup3*F(1523)+p3mup3*p4mup2*p4mup4*F(1524)+p3mup4*p4mup2*
     -   p4mup3*F(1525)+p1mup4*p4mup2*p4mup3*F(1526)+p4mup3*p4mup4*p6
     -   mup2*F(1527)+p3mup4*p4mup3*p6mup2*F(1529)+p2mup4*p4mup2*p6mu
     -   p3*F(1533)+p2mup3*p4mup2*p4mup4*F(1544)+p3mup2*(p2mup4*p4mup
     -   3*F(1507)+p2mup4*p6mup3*F(1508)+8*mup3mup4*F(1517)+p3mup4*(-
     -   16*p6mup3*F(1497)+p4mup3*F(1521))-16*(p1mup4*p4mup3*F(1513)+
     -   p4mup4*p6mup3*F(1515)+p3mup3*p4mup4*F(1520)-p2mup3*(p2mup4*F
     -   (1501)+p4mup4*F(1504)-p1mup4*F(1531))-p1mup3*(p2mup4*F(1499)
     -   +p4mup4*F(1511)-p1mup4*F(1535))-p1mup4*p6mup3*F(1536))+p4mup
     -   3*p4mup4*F(1545))+p1mup3*p3mup4*p6mup2*F(1547)+p1mup3*p4mup2
     -   *p4mup4*F(1548)+p4mup2*p4mup4*p6mup3*F(1550)+16*(p3mup3*p3mu
     -   p4*p6mup2*F(1498)+p2mup3*p4mup2*(p1mup4*F(1500)+p2mup4*F(150
     -   2))+p2mup3*p2mup4*p6mup2*F(1503)+p2mup3*p3mup4*p6mup2*F(1505
     -   )+p2mup3*p4mup4*p6mup2*F(1506)+p2mup4*p3mup3*p4mup2*F(1509)+
     -   p2mup4*p3mup3*p6mup2*F(1510)-p1mup3*p3mup4*p4mup2*F(1512)+p1
     -   mup4*p3mup3*p4mup2*F(1514)+p3mup3*p4mup4*p6mup2*F(1516)-p1mu
     -   p3*p2mup4*p4mup2*F(1518)+p2mup3*p3mup4*p4mup2*F(1519)+p3mup3
     -   *p3mup4*p4mup2*F(1522)+p2mup4*p4mup3*p6mup2*F(1528)+p1mup3*p
     -   2mup4*p6mup2*F(1530)+p1mup4*p2mup3*p6mup2*F(1532)-p2mup4*p6m
     -   up2*p6mup3*F(1534)+p1mup4*p3mup3*p6mup2*F(1537)-p3mup4*p6mup
     -   2*p6mup3*F(1538)-p1mup3*p1mup4*p4mup2*F(1539)+p1mup3*p4mup4*
     -   p6mup2*F(1540)+p1mup4*p4mup2*p6mup3*F(1541)+p1mup4*p4mup3*p6
     -   mup2*F(1542)-p4mup4*p6mup2*p6mup3*F(1543)+p3mup4*p4mup2*p6mu
     -   p3*F(1546)-p4mup2*p4mup3*p4mup4*F(1549)+p1mup3*p1mup4*p6mup2
     -   *F(1551)-p1mup4*p6mup2*p6mup3*F(1552))
       Fa25(267) = 4*mup2mup3*F(1553)+p2mup2*(p4mup3*F(1565)-8*(p2mup3
     -   *F(1555)-p3mup3*F(1557)-p1mup3*F(1563)+p6mup3*F(1569)))+p1mu
     -   p2*(8*(p2mup3*F(1554)+p3mup3*F(1556)-p1mup3*F(1571)+p6mup3*F
     -   (1572))+p4mup3*F(1575))+8*(p3mup2*p6mup3*F(1560)-p3mup3*p6mu
     -   p2*F(1561)-p3mup2*p3mup3*F(1564)+p2mup3*(p3mup2*F(1559)-p6mu
     -   p2*F(1562)+p4mup2*F(1567))+p3mup3*p4mup2*F(1568)-p4mup2*p6mu
     -   p3*F(1570)+p1mup3*(p3mup2*F(1558)+p4mup2*F(1566)-p6mup2*F(15
     -   73))+p6mup2*p6mup3*F(1574)-p3mup2*p4mup3*F(1576)+p4mup2*p4mu
     -   p3*F(1577)+p4mup3*p6mup2*F(1578))
       Fa25(268) = p1mup4*F(1579)+p2mup4*F(1580)+p3mup4*F(1581)+4*p4mu
     -   p4*F(1582)+p6mup4*F(1583)
       Fa25(269) = 4*mup2mup4*F(1584)+p1mup2*p6mup4*F(1601)-8*(p3mup4*
     -   p4mup2*F(1588)-p1mup2*p3mup4*F(1590)-p1mup2*p4mup4*F(1591)+p
     -   1mup4*p2mup2*F(1592)-p2mup2*p6mup4*F(1593)+p1mup4*p4mup2*F(1
     -   594)+p3mup4*p6mup2*F(1596)-p4mup2*p6mup4*F(1597)+p4mup4*p6mu
     -   p2*F(1598)+p3mup2*(p3mup4*F(1586)+p4mup4*F(1587)-p1mup4*F(15
     -   99)-p6mup4*F(1600))-p1mup4*p6mup2*F(1602)+p2mup2*p3mup4*F(16
     -   03)+p2mup2*p4mup4*F(1604)-p1mup2*p1mup4*F(1606)+p2mup4*(p3mu
     -   p2*F(1585)-p1mup2*F(1589)+p6mup2*F(1595)+p4mup2*F(1605)+p2mu
     -   p2*F(1607))+p4mup2*p4mup4*F(1608)+p6mup2*p6mup4*F(1609))
       Fa25(270) = 4*mup3mup4*F(1610)+p3mup3*p3mup4*F(1615)+p4mup3*p4m
     -   up4*F(1627)+p1mup3*p1mup4*F(1634)+8*(p3mup3*p4mup4*F(1616)+p
     -   3mup4*p4mup3*F(1617)+p1mup3*p3mup4*F(1619)+p1mup3*p4mup4*F(1
     -   620)-p4mup3*p6mup4*F(1621)+p1mup4*p4mup3*F(1622)+p3mup4*p6mu
     -   p3*F(1624)+p4mup4*p6mup3*F(1625)-p1mup4*p3mup3*F(1626)+p2mup
     -   4*(p3mup3*F(1614)+p1mup3*F(1618)+p6mup3*F(1623)+p4mup3*F(162
     -   8))-p1mup3*p6mup4*F(1630)+p2mup3*(p2mup4*F(1611)+p3mup4*F(16
     -   12)+p4mup4*F(1613)-p1mup4*F(1629)-p6mup4*F(1631))-p3mup3*p6m
     -   up4*F(1632)-p1mup4*p6mup3*F(1633)+p6mup3*p6mup4*F(1635))
       Fa25(271) = 4*mup2mup4*F(1636)+p1mup2*(p3mup4*F(1638)+8*(p2mup4
     -   *F(1637)+p4mup4*F(1639)-p1mup4*F(1656))+p6mup4*F(1657))-8*(p
     -   1mup4*p4mup2*F(1643)+p2mup4*p4mup2*F(1644)+p2mup4*p6mup2*F(1
     -   645)-p3mup4*p6mup2*F(1646)-p4mup2*p6mup4*F(1647)+p4mup4*p6mu
     -   p2*F(1648)-p2mup4*p3mup2*F(1649)-p3mup2*p3mup4*F(1650)-p3mup
     -   2*p4mup4*F(1651)+p2mup2*(p2mup4*F(1640)-p3mup4*F(1641)+p4mup
     -   4*F(1642)+p1mup4*F(1652)-p6mup4*F(1653))-p3mup4*p4mup2*F(165
     -   4)+p4mup2*p4mup4*F(1655)-p1mup4*p6mup2*F(1658)+p6mup2*p6mup4
     -   *F(1659)+p1mup4*p3mup2*F(1660)-p3mup2*p6mup4*F(1661))
       Fa25(272) = p1mup2*F(1664)+4*(p4mup2*F(1663)+p6mup2*F(1665)+p3m
     -   up2*F(1666)+p2mup2*F(1667))
       Fa25(273) = p3mup2*F(1670)-p1mup2*F(1671)-p6mup2*F(1672)+p2mup2
     -   *F(1673)+p4mup2*F(1674)
       Fa25(274) = p2mup3*F(1675)+p3mup3*F(1676)+p1mup3*F(1677)+p6mup3
     -   *F(1678)+p4mup3*F(1679)
       Fa25(275) = -(p1mup2*F(1680))+p2mup2*F(1681)+p6mup2*F(1682)+p3m
     -   up2*F(1683)+p4mup2*F(1684)
       Return
       End
