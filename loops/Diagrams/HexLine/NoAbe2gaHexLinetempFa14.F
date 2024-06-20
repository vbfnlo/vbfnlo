c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa14(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa14)
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
       Complex*16 Fa14(313:336)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa14(313) = mup2mup6*p4mup5*F(1615)+4*(p1mup2*(-(p1mup5*p6mup6*
     -   F(1617))-p2mup5*p6mup6*F(1620)+p2mup6*p5mup5*F(1621)-p5mup5*
     -   p6mup6*F(1627)+p6mup5*p6mup6*F(1628)+2*(p5mup5*p5mup6*F(1626
     -   )+p4mup5*p6mup6*F(1632))+p1mup6*(p5mup5*F(1616)-4*(p2mup5*F(
     -   1618)-p6mup5*F(1622)-p1mup5*F(1634)))+4*(p1mup5*p2mup6*F(161
     -   9)-p1mup5*p5mup6*F(1623)-p2mup5*p5mup6*F(1624)+p2mup6*p6mup5
     -   *F(1625)+p5mup6*p6mup5*F(1629)+p4mup5*p5mup6*F(1630)+p2mup5*
     -   p2mup6*F(1635)+p4mup5*p4mup6*F(1636)))+2*p2mup2*p4mup5*(p6mu
     -   p6*F(1633)+2*(p5mup6*F(1631)+p4mup6*F(9613))))
       Fa14(314) = mup2mup6*p3mup5*F(1615)+2*(mup2mup5*p6mup6*F(1638)+
     -   2*(mup2mup5*p1mup6*F(1637)-mup2mup5*p2mup6*F(1639)-mup2mup5*
     -   p5mup6*F(1640))+4*(p2mup2*(p3mup5*p6mup6*F(1633)-2*(-(p3mup5
     -   *p5mup6*F(1631))+p1mup6*p3mup5*F(1641)+p1mup6*p4mup5*F(1641)
     -   -p2mup6*(p3mup5+p4mup5)*F(1645)))+p1mup2*(p3mup5*p6mup6*F(16
     -   32)-2*(-(p3mup5*p5mup6*F(1630))+p1mup6*p3mup5*F(1642)+p1mup6
     -   *p4mup5*F(1642)-p2mup6*(p3mup5+p4mup5)*F(1646)))+p3mup5*p5mu
     -   p2*p6mup6*F(1652)+p4mup5*p5mup2*p6mup6*F(1652)+p3mup2*p3mup5
     -   *p6mup6*F(1653)+p3mup5*p4mup2*p6mup6*F(1653)+p3mup2*p4mup5*p
     -   6mup6*F(1653)+p4mup2*p4mup5*p6mup6*F(1653)+p2mup5*p6mup2*p6m
     -   up6*F(9578))-8*(p1mup6*(p3mup5+p4mup5)*p5mup2*F(1643)+p2mup6
     -   *p3mup5*p5mup2*F(1644)+p2mup6*p4mup5*p5mup2*F(1644)+p1mup6*p
     -   2mup5*p6mup2*F(1647)+p1mup6*p3mup2*p3mup5*F(1648)+p1mup6*p3m
     -   up5*p4mup2*F(1648)+p1mup6*p3mup2*p4mup5*F(1648)+p1mup6*p4mup
     -   2*p4mup5*F(1648)+p2mup6*p3mup2*p3mup5*F(1649)+p2mup6*p3mup5*
     -   p4mup2*F(1649)+p2mup6*p3mup2*p4mup5*F(1649)+p2mup6*p4mup2*p4
     -   mup5*F(1649)-p3mup2*p3mup5*p5mup6*F(1650)-p3mup5*p4mup2*p5mu
     -   p6*F(1650)-p3mup2*p4mup5*p5mup6*F(1650)-p4mup2*p4mup5*p5mup6
     -   *F(1650)+p3mup5*p5mup2*p5mup6*F(1651)+p4mup5*p5mup2*p5mup6*F
     -   (1651)+p2mup5*p5mup6*p6mup2*F(9580)))
       Fa14(315) = p4mup6*(p1mup2*(-(p5mup5*F(1656))+2*(p3mup5*F(1636)
     -   +p1mup5*F(1654)+p2mup5*F(1658)))+p2mup2*(p5mup5*F(1657)+2*(p
     -   1mup5*F(1655)+p2mup5*F(1659)+p3mup5*F(9613))))
       Fa14(316) = -(p2mup6*p5mup5*p6mup2*F(1660))+4*(p1mup6+p2mup6)*(
     -   p3mup5+p4mup5)*p6mup2*F(9582)
       Fa14(317) = p4mup6*(-(mup2mup5*F(1661))+2*(-2*(p1mup2*p6mup5*F(
     -   1664)-p2mup2*p6mup5*F(1665))+p5mup2*(2*(p6mup5*F(1585)+p1mup
     -   5*F(1662)+p2mup5*F(1663))+p5mup5*F(1669))+(p3mup2+p4mup2)*(-
     -   (p5mup5*F(1668))+2*(p1mup5*F(1666)+p2mup5*F(1667)-p6mup5*F(1
     -   670)))))+p6mup2*(-(p1mup6*p5mup5*F(1660))+2*p1mup5*(p6mup6*F
     -   (9578)+2*(p2mup6*F(1647)-p5mup6*F(9580))))+p3mup6*(p1mup2*(-
     -   2*p5mup5*F(1656)+4*(p4mup5*F(1636)+p1mup5*F(1654)+p2mup5*F(1
     -   658)))-mup2mup5*F(1661)+2*(p2mup2*p5mup5*F(1657)+(p3mup2+p4m
     -   up2)*(2*(p1mup5*F(1666)+p2mup5*F(1667))-p5mup5*F(1668))+p5mu
     -   p2*(2*(p6mup5*F(1585)+p1mup5*F(1662)+p2mup5*F(1663))+p5mup5*
     -   F(1669)))+4*(p1mup5*p2mup2*F(1655)+p2mup2*p2mup5*F(1659)+p6m
     -   up5*(-(p1mup2*F(1664))+p2mup2*F(1665)-(p3mup2+p4mup2)*F(1670
     -   ))+p2mup2*p4mup5*F(9613)))
       Fa14(318) = p1mup5*(p3mup6+p4mup6)*p6mup2*F(9582)+p2mup5*(p3mup
     -   6+p4mup6)*p6mup2*F(9582)-p4mup2*(p3mup6*p4mup5+(p3mup5+p4mup
     -   5)*p4mup6)*F(9608)-p3mup5*p3mup6*p5mup2*F(9611)
       Fa14(319) = p1mup2*p3mup5*p3mup6*F(1636)+p3mup5*p3mup6*p4mup2*F
     -   (9608)+p3mup2*(p3mup5+p4mup5)*(p3mup6+p4mup6)*F(9608)+p3mup6
     -   *p4mup5*p5mup2*F(9611)+p3mup5*p4mup6*p5mup2*F(9611)+p4mup5*p
     -   4mup6*p5mup2*F(9611)+p2mup2*p3mup5*p3mup6*F(9613)
       Fa14(320) = mup3mup6*(p2mup2*F(1671)+p5mup2*F(1672)+(p3mup2+p4m
     -   up2)*F(1673)+p1mup2*F(1674))
       Fa14(321) = mup3mup5*(-(p1mup2*F(1674))+p2mup2*F(1675)+p5mup2*F
     -   (1676)+(p3mup2+p4mup2)*F(1677))
       Fa14(322) = p2mup2*F(1671)+p5mup2*F(1672)+(p3mup2+p4mup2)*F(167
     -   3)+p1mup2*F(1674)
       Fa14(323) = 2*(mup5mup6*F(1678)+p5mup5*p6mup6*F(1690)-p6mup5*p6
     -   mup6*F(1691)+p1mup6*(-(p5mup5*F(1680))+4*(p2mup5*F(1682)-p6m
     -   up5*F(1684)+(p3mup5+p4mup5)*F(1692))))+p3mup5*p6mup6*F(1695)
     -   +p4mup5*p6mup6*F(1695)-4*(p2mup6*p5mup5*F(1687)-p5mup5*p5mup
     -   6*F(1689)+p3mup6*p5mup5*F(1698)+p4mup6*p5mup5*F(1698))+p1mup
     -   5*(p6mup6*F(1679)-8*(p2mup6*F(1681)-p5mup6*F(1683)+(p3mup6+p
     -   4mup6)*F(1693)-p1mup6*F(1701)))+p2mup5*(p6mup6*F(1686)+8*(p5
     -   mup6*F(1685)-(p3mup6+p4mup6)*F(1700)-p2mup6*F(1702)))+8*(-(p
     -   2mup6*p6mup5*F(1688))+p3mup5*p5mup6*F(1694)+p4mup5*p5mup6*F(
     -   1694)-p5mup6*p6mup5*F(1696)+p2mup6*p3mup5*F(1697)+p2mup6*p4m
     -   up5*F(1697)-p3mup6*p6mup5*F(1699)-p4mup6*p6mup5*F(1699)+(p3m
     -   up5+p4mup5)*(p3mup6+p4mup6)*F(9592))
       Fa14(324) = p1mup5*F(1703)+p5mup5*F(1704)+p6mup5*F(1705)+p2mup5
     -   *F(1706)+(p3mup5+p4mup5)*F(1707)
       Fa14(325) = p1mup5*F(1708)+p5mup5*F(1709)+p6mup5*F(1710)+p2mup5
     -   *F(1711)+(p3mup5+p4mup5)*F(1712)
       Fa14(326) = p1mup6*F(1713)+p5mup6*F(1714)+p6mup6*F(1715)+p2mup6
     -   *F(1716)+(p3mup6+p4mup6)*F(1717)
       Fa14(327) = p1mup6*F(1718)+p5mup6*F(1719)+p6mup6*F(1720)+p2mup6
     -   *F(1721)+(p3mup6+p4mup6)*F(1722)
       Fa14(328) = mup3mup6*(p1mup5*F(1726)+p5mup5*F(1728)+2*p6mup5*F(
     -   1730)+p2mup5*F(1732)+(p3mup5+p4mup5)*F(1734))+2*mup3mup5*(p1
     -   mup6*F(1727)+p5mup6*F(1729)-p6mup6*F(1731)+p2mup6*F(1733)+(p
     -   3mup6+p4mup6)*F(1735))
       Fa14(329) = mup3mup6*(p3mup2*p6mup5*F(1753)+p4mup2*p6mup5*F(175
     -   3)+4*(-(p2mup2*p5mup5*F(1737))-p1mup2*p5mup5*F(1741)+p5mup2*
     -   p5mup5*F(1746)+p5mup5*p6mup2*F(1749)+p3mup2*p5mup5*F(1752)+p
     -   4mup2*p5mup5*F(1752)+4*(p6mup2*p6mup5*F(1361)-p2mup2*p6mup5*
     -   F(1739)-p1mup2*p6mup5*F(1742)+p5mup2*p6mup5*F(1747)+mup2mup5
     -   *F(1759))+8*(p2mup2*p2mup5*F(1760)+p1mup2*p2mup5*F(1762)))+1
     -   6*(p4mup5*p6mup2*F(1757)+p3mup5*(p6mup2*F(1757)+2*(p4mup2*F(
     -   1321)+p5mup2*F(1756)-p2mup2*F(1764)-p1mup2*F(1765)))+p1mup5*
     -   p6mup2*F(1773)+p2mup5*p6mup2*F(1776))+32*(p4mup2*p4mup5*F(13
     -   21)+p3mup2*(p3mup5+p4mup5)*F(1321)+p4mup5*p5mup2*F(1756)-p2m
     -   up2*p4mup5*F(1764)-p1mup2*p4mup5*F(1765)+p1mup5*p2mup2*F(176
     -   7)+p1mup2*p1mup5*F(1770)-p1mup5*p5mup2*F(1771)-p2mup5*p5mup2
     -   *F(1774)+p1mup5*p3mup2*F(1777)+p1mup5*p4mup2*F(1777)-p2mup5*
     -   (p3mup2+p4mup2)*F(1778)))+4*mup3mup5*(mup2mup6*F(1736)-p5mup
     -   2*p6mup6*F(1748)+p6mup2*p6mup6*F(1751)-p3mup2*p6mup6*F(1755)
     -   -p4mup2*p6mup6*F(1755)+p2mup2*(4*p5mup6*F(1738)+p6mup6*F(174
     -   0)+8*((p3mup6+p4mup6)*F(1382)-p2mup6*F(1761)-p1mup6*F(1768))
     -   )+p1mup2*(4*p5mup6*F(1743)+p6mup6*F(1744)-8*(p2mup6*F(1763)-
     -   (p3mup6+p4mup6)*F(1766)+p1mup6*F(1769)))-4*(p5mup2*p5mup6*F(
     -   1745)+p5mup6*p6mup2*F(1750)+p3mup2*p5mup6*F(1754)+p4mup2*p5m
     -   up6*F(1754)+p3mup6*p5mup2*F(1758)+p4mup6*p5mup2*F(1758)-p1mu
     -   p6*p5mup2*F(1772)-p2mup6*p5mup2*F(1775))+8*(p3mup6*p6mup2*F(
     -   9582)+p4mup6*p6mup2*F(9582)+p1mup6*p6mup2*F(9585)+p2mup6*p6m
     -   up2*F(9586)-p3mup2*p3mup6*F(9591)-p3mup6*p4mup2*F(9591)-p3mu
     -   p2*p4mup6*F(9591)-p4mup2*p4mup6*F(9591)-p1mup6*p3mup2*F(9597
     -   )-p1mup6*p4mup2*F(9597)-p2mup6*(p3mup2+p4mup2)*F(9605)))
       Fa14(330) = p1mup6*F(1779)+p5mup6*F(1780)+p6mup6*F(1781)-p2mup6
     -   *F(1782)+(p3mup6+p4mup6)*F(1783)
       Fa14(331) = 2*(p5mup2*p6mup6*F(1748)+p3mup2*p6mup6*F(1755)+p4mu
     -   p2*p6mup6*F(1755)-p2mup2*(4*p5mup6*F(1738)+p6mup6*F(1740)+8*
     -   ((p3mup6+p4mup6)*F(1382)-p2mup6*F(1761)-p1mup6*F(1768)))-p1m
     -   up2*(4*p5mup6*F(1743)+p6mup6*F(1744)-8*(p2mup6*F(1763)-(p3mu
     -   p6+p4mup6)*F(1766)+p1mup6*F(1769)))+mup2mup6*F(1784))+p6mup2
     -   *p6mup6*F(1785)+8*(p5mup2*p5mup6*F(1745)+p5mup6*p6mup2*F(175
     -   0)+p3mup2*p5mup6*F(1754)+p4mup2*p5mup6*F(1754)+p3mup6*p5mup2
     -   *F(1758)+p4mup6*p5mup2*F(1758)+p1mup6*p5mup2*F(1786)+p2mup6*
     -   p5mup2*F(1787))-16*(p3mup6*p6mup2*F(9582)+p4mup6*p6mup2*F(95
     -   82)+p1mup6*p6mup2*F(9585)+p2mup6*p6mup2*F(9586)-p3mup2*p3mup
     -   6*F(9591)-p3mup6*p4mup2*F(9591)-p3mup2*p4mup6*F(9591)-p4mup2
     -   *p4mup6*F(9591)-p1mup6*p3mup2*F(9597)-p1mup6*p4mup2*F(9597)-
     -   p2mup6*(p3mup2+p4mup2)*F(9605))
       Fa14(332) = p1mup5*F(1788)+p5mup5*F(1789)+p6mup5*F(1790)+p2mup5
     -   *F(1791)+(p3mup5+p4mup5)*F(1792)
       Fa14(333) = -16*(p4mup2*p4mup5*F(1321)+p4mup5*p5mup2*F(1756)-p2
     -   mup2*p4mup5*F(1764)-p1mup2*p4mup5*F(1765)+p1mup5*p2mup2*F(17
     -   67)+p1mup2*p1mup5*F(1770)-p1mup5*p5mup2*F(1771)-p2mup5*p5mup
     -   2*F(1774)+p1mup5*p4mup2*F(1777)-p2mup5*p4mup2*F(1778))+p5mup
     -   5*p6mup2*F(1794)-8*(p4mup5*p6mup2*F(1757)+p3mup5*(p6mup2*F(1
     -   757)+2*(p4mup2*F(1321)+p5mup2*F(1756)-p2mup2*F(1764)-p1mup2*
     -   F(1765)))-p4mup2*p6mup5*F(1795))-2*(-(p2mup2*p5mup5*F(1737))
     -   -p1mup2*p5mup5*F(1741)+p4mup2*p5mup5*F(1752)+4*(p6mup2*p6mup
     -   5*F(1361)-p2mup2*p6mup5*F(1739)-p1mup2*p6mup5*F(1742)+p5mup2
     -   *p6mup5*F(1747)+mup2mup5*F(1759))+8*(p2mup2*p2mup5*F(1760)+p
     -   1mup2*p2mup5*F(1762))-p5mup2*p5mup5*F(1793)+p3mup2*(p5mup5*F
     -   (1752)+8*((p3mup5+p4mup5)*F(1321)+p1mup5*F(1777)-p2mup5*F(17
     -   78))-4*p6mup5*F(1795)))+p1mup5*p6mup2*F(1796)+p2mup5*p6mup2*
     -   F(1797)
       Fa14(334) = mup5mup6*(p5mup3*F(1798)+p2mup3*F(1799)+p4mup3*F(18
     -   00)-2*p3mup3*F(1801))
       Fa14(335) = mup3mup6*((p1mup5+p2mup5)*F(1802)+p5mup5*F(1804)+p6
     -   mup5*F(1805))+mup3mup5*(p6mup6*F(1806)+4*((p1mup6+p2mup6)*F(
     -   1803)+p5mup6*F(1807)))
       Fa14(336) = -2*mup5mup6*p6mup3*F(1808)+mup3mup6*p4mup5*F(1809)
       Return
       End
