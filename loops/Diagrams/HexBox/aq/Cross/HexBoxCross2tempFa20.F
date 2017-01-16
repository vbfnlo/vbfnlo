c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa20(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   20)
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
       Complex*16 Fa20(362:392)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa20(362) = p2mup3*F(1460)+p1mup3*F(1593)-p4mup3*F(1594)-p5mup3
     &   *F(1595)+p6mup3*F(1644)
       Fa20(363) = p5mup3*F(949)+p4mup3*F(950)+p2mup3*F(1464)+p1mup3*F
     &   (1602)+p6mup3*F(1645)
       Fa20(364) = p1mup3*p6mup6*F(1646)+p1mup3*p5mup6*F(1647)-p1mup3*
     &   p2mup6*F(1648)-p1mup3*p4mup6*F(1649)+p4mup3*p6mup6*F(1650)+p
     &   4mup3*p5mup6*F(1651)+p5mup3*p6mup6*F(1652)+p5mup3*p5mup6*F(1
     &   653)-mup3mup6*F(1654)+p1mup3*p1mup6*F(1655)+p1mup6*p5mup3*F(
     &   1656)+p1mup6*p4mup3*F(1658)+p2mup6*p5mup3*F(1659)+p2mup6*p4m
     &   up3*F(1661)+p4mup6*p5mup3*F(1663)+p6mup3*(-(p6mup6*F(1027))+
     &   p5mup6*F(1652)+p1mup6*F(1657)+p2mup6*F(1660)+p4mup6*F(1664))
     &   +p4mup3*p4mup6*F(1665)+p2mup3*(p5mup6*F(1197)-p2mup6*F(1373)
     &   -p1mup6*F(1622)+p6mup6*F(1662)-p4mup6*F(1666))
       Fa20(365) = p5mup3*p5mup6*F(1672)+p4mup3*p4mup6*F(1678)+4*(-(p2
     &   mup3*p2mup6*F(1087))-p1mup3*p2mup6*F(1611)+p1mup3*p1mup6*F(1
     &   614)+p1mup6*p5mup3*F(1618)-p1mup6*p2mup3*F(1623)+p2mup6*p5mu
     &   p3*F(1662)+p2mup3*p5mup6*F(1662)+p1mup3*p6mup6*F(1667)+p1mup
     &   3*p5mup6*F(1668)-p1mup3*p4mup6*F(1669)+p5mup3*p6mup6*F(1671)
     &   -mup3mup6*F(1673)+p2mup3*p6mup6*F(1674)+p4mup3*(p6mup6*F(102
     &   1)+p1mup6*F(1619)+p5mup6*F(1670)+p2mup6*F(1675))+p4mup6*p5mu
     &   p3*F(1676)+p6mup3*(-(p6mup6*F(1028))+p1mup6*F(1031)+p5mup6*F
     &   (1671)+p2mup6*F(1674)+p4mup6*F(1677))-p2mup3*p4mup6*F(1679))
       Fa20(366) = p1mup3*F(1680)+p2mup3*F(1681)+p4mup3*F(1682)+p5mup3
     &   *F(1683)+p6mup3*F(1684)
       Fa20(367) = p5mup6*p6mup3*F(1622)+p4mup6*p6mup3*F(1624)+mup3mup
     &   6*F(1680)+p1mup3*p6mup6*F(1685)-p2mup3*p6mup6*F(1686)+p1mup3
     &   *p2mup6*F(1687)-p1mup3*p4mup6*F(1688)+p6mup3*p6mup6*F(1690)+
     &   p4mup3*p6mup6*F(1691)-p1mup3*p5mup6*F(1692)+p2mup3*p2mup6*F(
     &   1693)-p2mup3*p5mup6*F(1694)-p2mup3*p4mup6*F(1695)-p1mup6*(p4
     &   mup3*F(872)+p6mup3*F(1288)-p2mup3*F(1687)+p1mup3*F(1689)+p5m
     &   up3*F(1696))-p2mup6*p6mup3*F(1698)-p2mup6*p4mup3*F(1699)+p4m
     &   up3*p5mup6*F(1700)+p5mup3*(p5mup6*F(1059)+p6mup6*F(1284)-p2m
     &   up6*F(1697)+p4mup6*F(1700))+p4mup3*p4mup6*F(1701)
       Fa20(368) = -(p4mup6*p6mup3*F(1661))+mup3mup6*F(1683)+p1mup3*p6
     &   mup6*F(1702)+p1mup3*p4mup6*F(1703)+p1mup3*p1mup6*F(1704)+p1m
     &   up3*p2mup6*F(1705)-p5mup3*p6mup6*F(1706)-p6mup3*p6mup6*F(170
     &   7)-p4mup3*p6mup6*F(1708)+p1mup6*p5mup3*F(1710)+p1mup6*p6mup3
     &   *F(1711)+p1mup6*p4mup3*F(1712)+p2mup6*p5mup3*F(1713)+p2mup6*
     &   p6mup3*F(1714)+p2mup6*p4mup3*F(1715)-p4mup6*p5mup3*F(1717)-p
     &   4mup3*p4mup6*F(1718)-p5mup6*(p6mup3*F(1659)-p1mup3*F(1709)+p
     &   4mup3*F(1717)+p5mup3*F(1719))+p2mup3*(p5mup6*F(1061)-p1mup6*
     &   F(1697)+p6mup6*F(1716)+p4mup6*F(1720)+p2mup6*F(1721))
       Fa20(369) = p2mup3*p5mup6*F(1373)+p2mup3*p4mup6*F(1374)-p4mup3*
     &   p4mup6*F(1400)+p1mup3*p5mup6*F(1610)+mup3mup6*F(1684)+p1mup3
     &   *p1mup6*F(1685)-p1mup3*p2mup6*F(1686)+p1mup6*p4mup3*F(1691)-
     &   p1mup6*p2mup3*F(1698)+p1mup3*p6mup6*F(1722)+p1mup3*p4mup6*F(
     &   1723)-p4mup3*p6mup6*F(1725)+p2mup3*p6mup6*F(1726)-p6mup3*(p6
     &   mup6*F(1016)+p5mup6*F(1662)+p4mup6*F(1675)-p1mup6*F(1690)-p2
     &   mup6*F(1726))+p2mup6*p4mup3*F(1727)-p4mup3*p5mup6*F(1728)-p5
     &   mup3*(p5mup6*F(1197)-p1mup6*F(1284)-p2mup6*F(1716)+p6mup6*F(
     &   1724)+p4mup6*F(1728))+p2mup3*p2mup6*F(1729)
       Fa20(370) = p6mup6*F(1730)+p1mup6*F(1731)+p4mup6*F(1732)+p5mup6
     &   *F(1733)+p2mup6*F(1734)
       Fa20(371) = p2mup6*F(945)+p5mup6*F(951)+p4mup6*F(952)+p1mup6*F(
     &   955)-p6mup6*F(1735)
       Fa20(372) = p2mup6*F(990)+p4mup6*F(997)+p1mup6*F(998)-p6mup6*F(
     &   1736)+p5mup6*F(1737)
       Fa20(373) = p2mup3*F(1682)-p1mup3*F(1738)+p4mup3*F(1739)+p5mup3
     &   *F(1740)+p6mup3*F(1741)
       Fa20(374) = -(p6mup3*p6mup6*F(1619))+p2mup6*p6mup3*F(1624)+p5mu
     &   p6*p6mup3*F(1625)+p4mup6*p6mup3*F(1627)-p1mup3*p2mup6*F(1688
     &   )-p2mup3*p2mup6*F(1695)+p2mup6*p4mup3*F(1701)+p2mup3*p6mup6*
     &   F(1723)-mup3mup6*F(1738)+p1mup3*p6mup6*F(1742)+p1mup3*p4mup6
     &   *F(1743)+p1mup6*(p4mup3*F(873)+p5mup3*F(875)+p6mup3*F(1289)-
     &   p2mup3*F(1688)+p1mup3*F(1744))-p4mup3*p6mup6*F(1745)+p1mup3*
     &   p5mup6*F(1746)+p2mup3*p5mup6*F(1747)+p2mup3*p4mup6*F(1748)+p
     &   4mup3*p5mup6*F(1749)+p4mup3*p4mup6*F(1750)+p5mup3*(-(p6mup6*
     &   F(1285))+p2mup6*F(1700)+p4mup6*F(1749)+p5mup6*F(1751))
       Fa20(375) = p5mup3*p6mup6*F(1651)+p6mup3*(p6mup6*F(1650)+p1mup6
     &   *F(1658)+p2mup6*F(1661)+p5mup6*F(1663)+p4mup6*F(1665))-p1mup
     &   3*p2mup6*F(1703)+p2mup6*p5mup3*F(1717)+p2mup6*p4mup3*F(1718)
     &   -mup3mup6*F(1740)+p1mup3*p6mup6*F(1752)-p1mup3*p4mup6*F(1753
     &   )-p1mup3*p1mup6*F(1754)+p4mup3*p6mup6*F(1755)-p1mup3*p5mup6*
     &   F(1756)+p1mup6*p5mup3*F(1757)+p1mup6*p4mup3*F(1758)+p4mup6*p
     &   5mup3*F(1759)+p4mup3*p5mup6*F(1759)+p4mup3*p4mup6*F(1760)+p5
     &   mup3*p5mup6*F(1761)-p2mup3*(p5mup6*F(1060)+p1mup6*F(1700)+p2
     &   mup6*F(1720)-p6mup6*F(1728)+p4mup6*F(1762))
       Fa20(376) = p4mup6*p6mup3*F(1678)+p4mup3*p6mup6*F(1765)-4*(p2mu
     &   p3*p2mup6*F(1374)-p2mup6*p4mup3*F(1400)+p1mup6*p2mup3*F(1624
     &   )+p2mup3*p5mup6*F(1666)-p2mup3*p6mup6*F(1675)-p6mup3*(p6mup6
     &   *F(1021)+p1mup6*F(1619)+p2mup6*F(1675)+p5mup6*F(1676))+p1mup
     &   3*p2mup6*F(1723)+mup3mup6*F(1741)+p1mup3*p1mup6*F(1742)-p1mu
     &   p6*p4mup3*F(1745)-p1mup3*p6mup6*F(1763)-p1mup3*p4mup6*F(1764
     &   )+p1mup3*p5mup6*F(1766)-p4mup3*p5mup6*F(1767)-p5mup3*(p5mup6
     &   *F(1201)+p1mup6*F(1285)+p6mup6*F(1670)+p2mup6*F(1728)+p4mup6
     &   *F(1767))-p4mup3*p4mup6*F(1768)+p2mup3*p4mup6*F(1769))
       Fa20(377) = p2mup6*F(1308)-p4mup6*F(1309)-p6mup6*F(1770)+p1mup6
     &   *F(1771)-p5mup6*F(1772)
       Fa20(378) = p6mup6*F(1773)-p1mup6*F(1774)-p2mup6*F(1775)+p4mup6
     &   *F(1776)+p5mup6*F(1777)
       Fa20(379) = p4mup6*F(923)+p1mup6*F(924)+p5mup6*F(934)+p6mup6*F(
     &   1778)-p2mup6*F(1779)
       Fa20(380) = p1mup6*F(897)+p5mup6*F(1780)+p6mup6*F(1781)-p2mup6*
     &   F(1782)-p4mup6*F(1783)
       Fa20(381) = -4*(p6mup6*F(1784)-p1mup6*F(1785)-p2mup6*F(1786))+p
     &   4mup6*F(1787)+p5mup6*F(1788)
       Fa20(382) = p1mup6*F(925)+p6mup6*F(1789)+p4mup6*F(1790)+p5mup6*
     &   F(1791)-p2mup6*F(1792)
       Fa20(383) = p1mup6*F(898)+p6mup6*F(1504)+p5mup6*F(1793)-p2mup6*
     &   F(1794)-p4mup6*F(1795)
       Fa20(384) = p4mup6*F(1005)+p6mup6*F(1796)-p1mup6*F(1797)-p2mup6
     &   *F(1798)+p5mup6*F(1799)
       Fa20(385) = p5mup6*F(1508)+p6mup6*F(1509)-p1mup6*F(1800)-p2mup6
     &   *F(1801)+p4mup6*F(1802)
       Fa20(386) = -(p1mup3*F(1803))+p4mup3*F(1804)+p5mup3*F(1805)+p6m
     &   up3*F(1806)-p2mup3*F(1807)
       Fa20(387) = p1mup3*p2mup6*F(1680)+p2mup3*p2mup6*F(1681)+p2mup6*
     &   p4mup3*F(1682)+p2mup3*p4mup6*F(1682)+p2mup3*p5mup6*F(1683)+p
     &   2mup6*p6mup3*F(1684)-p1mup3*p4mup6*F(1738)+p4mup3*p4mup6*F(1
     &   739)+p4mup3*p5mup6*F(1740)+p4mup6*p6mup3*F(1741)-p1mup3*p6mu
     &   p6*F(1811)-p2mup3*p6mup6*F(1813)+p5mup3*(p5mup6*F(1317)+p2mu
     &   p6*F(1683)+p4mup6*F(1740)+p6mup6*F(1814))+p6mup3*p6mup6*F(18
     &   15)+p4mup3*p6mup6*F(1816)-p1mup3*p5mup6*F(1817)+p5mup6*p6mup
     &   3*F(1820)+p1mup6*(p5mup3*F(1312)+p1mup3*F(1812)+p6mup3*F(181
     &   8)+p4mup3*F(1819)+p2mup3*F(1821))+mup3mup6*F(1853)
       Fa20(388) = p2mup6*F(1319)-p4mup6*F(1320)-p5mup6*F(1321)-p6mup6
     &   *F(1431)+p1mup6*F(1432)
       Fa20(389) = p2mup6*F(923)-p6mup6*F(1822)+p1mup6*F(1823)-p4mup6*
     &   F(1824)+p5mup6*F(1825)
       Fa20(390) = p2mup6*F(1790)-p6mup6*F(1826)+p1mup6*F(1827)-p4mup6
     &   *F(1828)+p5mup6*F(1829)
       Fa20(391) = -(p1mup6*F(1003))-p2mup6*F(1005)+p4mup6*F(1006)+p5m
     &   up6*F(1495)+p6mup6*F(1830)
       Fa20(392) = -(p2mup3*F(1319))+p4mup3*F(1320)+p5mup3*F(1321)-p1m
     &   up3*F(1834)+p6mup3*F(1852)
       Return
       End
