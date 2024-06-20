c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g28(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa28)
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
       Complex*16 Fa28(379:392)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa28(379) = 2*(p1mup6+p2mup6)*F(1107)+p4mup6*F(1569)+p5mup6*F(1
     -   570)+p3mup6*F(1571)
       Fa28(380) = -2*(p4mup4*F(1529)+p5mup4*F(1531))+(p1mup4+p2mup4)*
     -   F(1572)+p3mup4*F(1573)
       Fa28(381) = 6*(p3mup6*F(1073)+(p1mup6+p2mup6)*F(1074)+p4mup6*F(
     -   1083))+2*p5mup6*F(1084)
       Fa28(382) = p4mup2*(p3mup3*p5mup4*F(1575)+p3mup3*p4mup4*F(1577)
     -   +p1mup4*p3mup3*F(1579)+p2mup4*p3mup3*F(1579)+p1mup3*(p5mup4*
     -   F(1574)+p4mup4*F(1576)+p3mup4*F(1578)+(p1mup4+p2mup4)*F(1580
     -   ))+p2mup3*(p5mup4*F(1574)+p4mup4*F(1576)+p3mup4*F(1578)+(p1m
     -   up4+p2mup4)*F(1580))+p3mup3*p3mup4*F(1581)+p4mup3*p5mup4*F(1
     -   582)-p4mup4*p5mup3*F(1583)-p1mup4*p5mup3*F(1584)-p2mup4*p5mu
     -   p3*F(1584)-p3mup4*p5mup3*F(1585)-p5mup3*p5mup4*F(1586)+mup3m
     -   up4*F(1587)-p1mup4*p4mup3*F(1588)-p2mup4*p4mup3*F(1588)-2*p3
     -   mup4*p4mup3*F(1589)-p4mup3*p4mup4*F(1590))
       Fa28(383) = p3mup2*(p3mup3*p5mup4*F(1575)+p3mup3*p4mup4*F(1577)
     -   +p1mup4*p3mup3*F(1579)+p2mup4*p3mup3*F(1579)+p1mup3*(p5mup4*
     -   F(1574)+p4mup4*F(1576)+p3mup4*F(1578)+(p1mup4+p2mup4)*F(1580
     -   ))+p2mup3*(p5mup4*F(1574)+p4mup4*F(1576)+p3mup4*F(1578)+(p1m
     -   up4+p2mup4)*F(1580))+p3mup3*p3mup4*F(1581)+p4mup3*p5mup4*F(1
     -   582)-p4mup4*p5mup3*F(1583)-p1mup4*p5mup3*F(1584)-p2mup4*p5mu
     -   p3*F(1584)-p3mup4*p5mup3*F(1585)-p5mup3*p5mup4*F(1586)+mup3m
     -   up4*F(1587)-p1mup4*p4mup3*F(1588)-p2mup4*p4mup3*F(1588)-2*p3
     -   mup4*p4mup3*F(1589)-p4mup3*p4mup4*F(1590))
       Fa28(384) = (p1mup2+p2mup2)*(p3mup3*p5mup4*F(1575)+p3mup3*p4mup
     -   4*F(1577)+p1mup4*p3mup3*F(1579)+p2mup4*p3mup3*F(1579)+p1mup3
     -   *(p5mup4*F(1574)+p4mup4*F(1576)+p3mup4*F(1578)+(p1mup4+p2mup
     -   4)*F(1580))+p2mup3*(p5mup4*F(1574)+p4mup4*F(1576)+p3mup4*F(1
     -   578)+(p1mup4+p2mup4)*F(1580))+p3mup3*p3mup4*F(1581)+p4mup3*p
     -   5mup4*F(1582)-p4mup4*p5mup3*F(1583)-p1mup4*p5mup3*F(1584)-p2
     -   mup4*p5mup3*F(1584)-p3mup4*p5mup3*F(1585)-p5mup3*p5mup4*F(15
     -   86)+mup3mup4*F(1587)-p1mup4*p4mup3*F(1588)-p2mup4*p4mup3*F(1
     -   588)-2*p3mup4*p4mup3*F(1589)-p4mup3*p4mup4*F(1590))
       Fa28(385) = p3mup3*p5mup4*F(1575)+p3mup3*p4mup4*F(1577)+p1mup4*
     -   p3mup3*F(1579)+p2mup4*p3mup3*F(1579)+p1mup3*(p5mup4*F(1574)+
     -   p4mup4*F(1576)+p3mup4*F(1578)+(p1mup4+p2mup4)*F(1580))+p2mup
     -   3*(p5mup4*F(1574)+p4mup4*F(1576)+p3mup4*F(1578)+(p1mup4+p2mu
     -   p4)*F(1580))+p3mup3*p3mup4*F(1581)+p4mup3*p5mup4*F(1582)-p4m
     -   up4*p5mup3*F(1583)-p1mup4*p5mup3*F(1584)-p2mup4*p5mup3*F(158
     -   4)-p3mup4*p5mup3*F(1585)-p5mup3*p5mup4*F(1586)+mup3mup4*F(15
     -   87)-2*p3mup4*p4mup3*F(1589)-p4mup3*p4mup4*F(1590)+(p1mup4+p2
     -   mup4)*p4mup3*F(1591)
       Fa28(386) = mup4mup6*F(1592)+p2mup4*p5mup6*F(1598)+p3mup4*(p3mu
     -   p6*F(1593)+4*p5mup6*F(1596)+p1mup6*F(1599)+p2mup6*F(1600)+p4
     -   mup6*F(1601))+p1mup4*p4mup6*F(1605)+p2mup4*p2mup6*F(1607)+p5
     -   mup4*p5mup6*F(1610)+p3mup6*(p2mup4*F(1595)+p1mup4*F(1609)+4*
     -   (p4mup4*F(1594)+p5mup4*F(1613)))+4*(p4mup4*p5mup6*F(1597)-p1
     -   mup6*p4mup4*F(1602)+p2mup6*p4mup4*F(1603)+p4mup4*p4mup6*F(16
     -   04)+p1mup6*p2mup4*F(1606)+p2mup4*p4mup6*F(1608)-p1mup6*p5mup
     -   4*F(1611)+p2mup6*p5mup4*F(1612)+p4mup6*p5mup4*F(1614))+p1mup
     -   4*p5mup6*F(1615)+p1mup4*p1mup6*F(1616)+p1mup4*p2mup6*F(1617)
       Fa28(387) = mup3mup6*F(1618)+p1mup6*p4mup3*F(1624)+p3mup3*(p4mu
     -   p6*F(1619)+p3mup6*F(1627)-4*(p1mup6*F(1629)+p2mup6*F(1631)+p
     -   5mup6*F(1632)))+p1mup3*p1mup6*F(1638)-4*(p3mup6*p4mup3*F(162
     -   0)+p4mup3*p4mup6*F(1622)+p1mup3*p4mup6*F(1623)+p2mup6*p4mup3
     -   *F(1625)+p4mup3*p5mup6*F(1626)-p1mup3*p3mup6*F(1628)+p3mup6*
     -   p5mup3*F(1633)-p1mup6*p5mup3*F(1634)+p2mup6*p5mup3*F(1635)+p
     -   5mup3*p5mup6*F(1636)+p4mup6*p5mup3*F(1637)-p1mup3*p2mup6*F(1
     -   639)-p1mup3*p5mup6*F(1642))+p2mup3*(p4mup6*F(1621)+p3mup6*F(
     -   1630)+p1mup6*F(1640)-4*(p2mup6*F(1641)+p5mup6*F(1643)))
       Fa28(388) = p2mup6*F(1644)+2*p4mup6*F(1645)+p5mup6*F(1646)+p1mu
     -   p6*F(1647)+p3mup6*F(1648)
       Fa28(389) = mup2mup6*F(1649)+p1mup2*p4mup6*F(1652)+p2mup2*(p4mu
     -   p6*F(1653)+4*p5mup6*F(1654)+p1mup6*F(1655)+p2mup6*F(1656))+p
     -   3mup2*(p4mup6*F(1650)+p5mup6*F(1661)+4*p1mup6*F(1662)+p2mup6
     -   *F(1663)+p3mup6*F(1664))+4*(p4mup2*p4mup6*F(1657)+p2mup6*p4m
     -   up2*F(1658)-p1mup6*p4mup2*F(1659)+p4mup2*p5mup6*F(1660)+p4mu
     -   p6*p5mup2*F(1666)+p2mup6*p5mup2*F(1667)+p1mup6*p5mup2*F(1669
     -   )+p5mup2*p5mup6*F(1670))+p3mup6*(p1mup2*F(1665)+4*(p4mup2*F(
     -   1651)+p5mup2*F(1668)+p2mup2*F(1671)))+p1mup2*p5mup6*F(1672)+
     -   p1mup2*p1mup6*F(1673)+p1mup2*p2mup6*F(1674)
       Fa28(390) = 2*mup4mup6*F(1676)+p2mup4*p5mup6*F(1679)+p3mup4*(p2
     -   mup6*F(1681)+4*(p5mup6*F(1677)+p1mup6*F(1680)+p3mup6*F(1682)
     -   +p4mup6*F(1683)))+p2mup4*p2mup6*F(1689)+p2mup4*p3mup6*F(1690
     -   )+p2mup4*p4mup6*F(1691)+p4mup4*(p4mup6*F(1686)+4*(p5mup6*F(1
     -   678)-p1mup6*F(1684)+p3mup6*F(1685)+p2mup6*F(1692)))+p1mup4*p
     -   2mup6*F(1693)+4*(p1mup4*p4mup6*F(1687)+p1mup6*p2mup4*F(1688)
     -   +p5mup4*p5mup6*F(1694)+p1mup6*p5mup4*F(1695)+p3mup6*p5mup4*F
     -   (1696)+p4mup6*p5mup4*F(1697)+p2mup6*p5mup4*F(1698)+p1mup4*p5
     -   mup6*F(1699)+p1mup4*p1mup6*F(1700)+p1mup4*p3mup6*F(1701))
       Fa28(391) = mup3mup6*F(1702)+p1mup3*p2mup6*F(1710)+p2mup3*p2mup
     -   6*F(1712)+p2mup6*p3mup3*F(1714)+p2mup6*p5mup3*F(1716)-4*(p1m
     -   up3*p4mup6*F(1705)-p2mup3*p4mup6*F(1708)+p4mup3*(p4mup6*F(17
     -   03)+p3mup6*F(1704)-p1mup6*F(1706)+p5mup6*F(1707)+p2mup6*F(17
     -   09))-p1mup6*p2mup3*F(1711)+p2mup3*p3mup6*F(1713)+p2mup3*p5mu
     -   p6*F(1715)+p3mup3*p4mup6*F(1718)+p4mup6*p5mup3*F(1725)+(p3mu
     -   p3*p3mup6*F(1717)+p1mup3*p1mup6*F(1719)+p1mup3*p3mup6*F(1720
     -   )+p1mup6*p3mup3*F(1721)+p1mup3*p5mup6*F(1722)+p3mup3*p5mup6*
     -   F(1723)+p3mup6*p5mup3*F(1724)+p1mup6*p5mup3*F(1726)+p5mup3*p
     -   5mup6*F(1727))*P(1979))
       Fa28(392) = 2*mup2mup6*F(1728)+p4mup2*p4mup6*F(1731)+p3mup6*p4m
     -   up2*F(1732)+p2mup6*p4mup2*F(1736)+p3mup2*p4mup6*F(1737)+p2mu
     -   p2*(p4mup6*F(1730)+p2mup6*F(1735)+4*(p5mup6*F(1738)-p1mup6*F
     -   (1741))+p3mup6*F(1742))+p4mup6*p5mup2*F(1743)+p1mup6*p3mup2*
     -   F(1746)+p3mup2*p3mup6*F(1747)-4*(p1mup6*p4mup2*F(1733)-p4mup
     -   2*p5mup6*F(1734)-p2mup6*p3mup2*F(1739)-p2mup6*p5mup2*F(1744)
     -   -p3mup2*p5mup6*F(1745)-p5mup2*p5mup6*F(1748))+p3mup6*p5mup2*
     -   F(1749)+p1mup6*p5mup2*F(1750)+p1mup2*(p4mup6*F(1729)+4*p2mup
     -   6*F(1740)+p5mup6*F(1751)+p1mup6*F(1752)+p3mup6*F(1753))
       Return
       End
