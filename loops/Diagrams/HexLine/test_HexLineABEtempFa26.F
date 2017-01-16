c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineAbe26(p1mup2,p1mup3,p1mup4,p1m
     -   up5,p2mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup2,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mu
     -   p2,p6mup3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,
     -   mup3mup5,mup4mup5,Fa26)
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
       Complex*16 Fa26(276:286)
       COMMON/test_FaFunctionsAbe/Fa
       Fa26(276) = p2mup4*F(1685)+p3mup4*F(1686)+p4mup4*F(1687)+p1mup4
     -   *F(1688)+p6mup4*F(1689)
       Fa26(277) = -16*p1mup2*(p1mup3*(p4mup5*F(1698)-p2mup5*F(1701)-p
     -   3mup5*F(1702)+p6mup5*F(1703)+p1mup5*F(1707))+p2mup3*(p4mup5*
     -   F(1691)+p6mup5*F(1699)+p1mup5*F(1704)+p2mup5*F(1708)-p3mup5*
     -   F(1712))+p3mup3*(p4mup5*F(1692)+p1mup5*F(1705)-p6mup5*F(1706
     -   )+p2mup5*F(1709)+p3mup5*F(1713)))-p2mup2*(-(p1mup3*(p6mup5*F
     -   (1700)+16*(p2mup5*F(1693)+p3mup5*F(1695)-p4mup5*F(1697)-p1mu
     -   p5*F(1710))))+16*(p1mup5*p2mup3*F(1694)+p1mup5*p3mup3*F(1696
     -   )+p2mup3*p2mup5*F(1711)-p2mup3*p3mup5*F(1714)))
       Fa26(278) = p2mup2*(p2mup3*(p4mup5*F(1715)-p6mup5*F(1717))+p3mu
     -   p3*(p4mup5*F(1716)-p6mup5*F(1718)+p2mup5*F(1719)+p3mup5*F(17
     -   20)))
       Fa26(279) = 8*mup3mup5*(p1mup2*F(1724)+p2mup2*F(1725)-p6mup2*F(
     -   1727))+p6mup5*(p2mup3*p4mup2*F(1728)+p3mup3*p4mup2*F(1730)-1
     -   6*p4mup3*(p3mup2*F(1729)+p4mup2*F(1731)))+mup2mup5*(p4mup3*F
     -   (1723)+8*(p1mup3*F(1721)+p2mup3*F(1722)+p3mup3*F(1726)+p6mup
     -   3*F(1736)))+p6mup5*(p1mup2*(p4mup3*F(1732)+16*p6mup3*F(1746)
     -   )-8*(mup2mup3*F(1735)-2*(p2mup2*p4mup3*F(1737)-p2mup3*p3mup2
     -   *F(1738)-p1mup3*(-(p4mup2*F(1733))+p3mup2*F(1739)+p6mup2*F(1
     -   740))-p3mup2*p3mup3*F(1741)+p2mup3*p6mup2*F(1742)+p3mup3*p6m
     -   up2*F(1743)-p4mup3*p6mup2*F(1744)-p6mup3*(p4mup2*F(1734)-p2m
     -   up2*F(1745)-p3mup2*F(1747)+p6mup2*F(1748)))))
       Fa26(280) = mup3mup5*(p4mup2*F(1749)-p3mup2*F(1750))+2*p6mup3*(
     -   -(p2mup2*(p1mup5*F(1752)+p4mup5*F(1753)))-p1mup5*p4mup2*F(17
     -   55)-p2mup2*p3mup5*F(1757)+p1mup5*p3mup2*F(1759)+p3mup2*p4mup
     -   5*F(1760)-p3mup5*p4mup2*F(1761)+p4mup5*p6mup2*F(1763)+p1mup5
     -   *p6mup2*F(1765)+p3mup2*p3mup5*F(1766)+p4mup2*p4mup5*F(1767)+
     -   p2mup5*(p4mup2*F(1754)-p2mup2*F(1756)-p3mup2*F(1758)-p6mup2*
     -   F(1768))+p1mup2*(p2mup5*F(1751)-p4mup5*F(1762)+p1mup5*F(1764
     -   )+p3mup5*F(1769))-p3mup5*p6mup2*F(1770))
       Fa26(281) = 8*mup2mup3*(p1mup5*F(1771)+p2mup5*F(1772)-p3mup5*F(
     -   1773)-p4mup5*F(1774))+16*p4mup3*(p1mup2*(p3mup5*F(1775)-p2mu
     -   p5*F(1777)+p1mup5*F(1780))+p2mup2*(p3mup5*F(1776)+p1mup5*F(1
     -   778)+p4mup5*F(1779)+p2mup5*F(1781)))+p1mup2*p4mup3*p4mup5*F(
     -   1782)
       Fa26(282) = p2mup3*p4mup2*p4mup5*F(1799)+p1mup5*(p4mup2*p4mup3*
     -   F(1798)+16*(p3mup3*p4mup2*F(1787)+p1mup3*p4mup2*F(1804)+p4mu
     -   p3*p6mup2*F(1808)+p2mup3*(p4mup2*F(1783)+p6mup2*F(1806)-p3mu
     -   p2*F(1818))-p3mup3*p6mup2*F(1820)+p3mup2*(p4mup3*F(1786)-p3m
     -   up3*F(1816)-p1mup3*F(1823))+p1mup3*p6mup2*F(1826)))+16*(p2mu
     -   p3*p2mup5*p4mup2*F(1788)-p2mup5*p3mup2*p4mup3*F(1789)+p2mup5
     -   *p3mup3*p4mup2*F(1790)+p2mup3*p3mup2*p4mup5*F(1792)+p2mup3*p
     -   3mup5*p4mup2*F(1793)-p3mup2*p3mup3*p4mup5*F(1794)-p3mup2*p3m
     -   up5*p4mup3*F(1795)+p3mup3*p3mup5*p4mup2*F(1796)-p2mup5*p4mup
     -   2*p4mup3*F(1800)-p3mup2*p4mup3*p4mup5*F(1801)+p3mup3*p4mup2*
     -   p4mup5*F(1802)-p3mup5*p4mup2*p4mup3*F(1803)+p2mup3*p4mup5*p6
     -   mup2*F(1809)-p3mup3*p4mup5*p6mup2*F(1810)+p3mup5*p4mup3*p6mu
     -   p2*F(1811)+p4mup2*p4mup3*p4mup5*F(1812)+p2mup5*p4mup3*p6mup2
     -   *F(1813)+p4mup3*p4mup5*p6mup2*F(1814)+p1mup3*(p3mup5*p4mup2*
     -   F(1785)+p2mup5*p4mup2*F(1791)+p4mup2*p4mup5*F(1797)-p2mup5*p
     -   6mup2*F(1805)+p3mup2*(-(p4mup5*F(1784))+p3mup5*F(1815)+p2mup
     -   5*F(1817))-p6mup2*(p4mup5*F(1807)+p3mup5*F(1819)))-p2mup3*p3
     -   mup5*p6mup2*F(1821)-p2mup5*p3mup3*p6mup2*F(1822)+p2mup3*p2mu
     -   p5*p3mup2*F(1824)-p2mup5*p3mup2*p3mup3*F(1825)+p2mup3*p2mup5
     -   *p6mup2*F(1827)+p2mup3*p3mup2*p3mup5*F(1828)-p3mup2*p3mup3*p
     -   3mup5*F(1829)+p3mup3*p3mup5*p6mup2*F(1830))
       Fa26(283) = p2mup5*F(1832)+p4mup5*F(1833)+4*(p1mup5*F(1831)+p6m
     -   up5*F(1834))+p3mup5*F(1835)
       Fa26(284) = 4*mup2mup5*F(1836)+p1mup2*(p3mup5*F(1845)+8*(p4mup5
     -   *F(1839)+p6mup5*F(1846)+p1mup5*F(1855))+p2mup5*F(1856))+8*(p
     -   2mup2*(p1mup5*F(1840)-p6mup5*F(1841)-p2mup5*F(1848))+p4mup2*
     -   (p3mup5*F(1838)+p1mup5*F(1842)+p6mup5*F(1843)+p2mup5*F(1850)
     -   )+p4mup5*(p2mup2*F(1849)+p4mup2*F(1851))-p2mup2*p3mup5*F(185
     -   7)-p6mup2*(p1mup5*F(1847)+p4mup5*F(1852)+p2mup5*F(1858)-p3mu
     -   p5*F(1859)-p6mup5*F(1860))+p3mup2*(p4mup5*F(1837)-p1mup5*F(1
     -   844)+p2mup5*F(1853)+p6mup5*F(1854)+p3mup5*F(1861)))
       Fa26(285) = 4*mup3mup5*F(1862)+8*(p1mup3*p4mup5*F(1866)-p4mup3*
     -   (p3mup5*F(1865)+p6mup5*F(1867)+p1mup5*F(1868))-p1mup3*p3mup5
     -   *F(1869)-p4mup3*p4mup5*F(1871)+p1mup3*p6mup5*F(1872)-p1mup5*
     -   p6mup3*F(1873)+p4mup5*p6mup3*F(1874)+p2mup5*p4mup3*F(1875)-p
     -   1mup3*p2mup5*F(1876)+p2mup5*p6mup3*F(1882)-p3mup5*p6mup3*F(1
     -   883)+p1mup3*p1mup5*F(1884)+p2mup3*(p4mup5*F(1863)+p1mup5*F(1
     -   877)-p3mup5*F(1878)+p6mup5*F(1879)+p2mup5*F(1885))+p3mup3*(p
     -   4mup5*F(1864)+p1mup5*F(1870)+p2mup5*F(1880)+p6mup5*F(1881)+p
     -   3mup5*F(1886))-p6mup3*p6mup5*F(1887))
       Fa26(286) = -4*mup2mup5*F(1888)+8*(-(p3mup2*p4mup5*F(1894))+p4m
     -   up2*(p1mup5*F(1891)+p2mup5*F(1892)-p6mup5*F(1893)+p3mup5*F(1
     -   898))+p4mup2*p4mup5*F(1899)+p4mup5*p6mup2*F(1900)+p1mup5*p3m
     -   up2*F(1901)+p2mup5*p3mup2*F(1902)+p3mup2*p6mup5*F(1903)+p2mu
     -   p2*(p4mup5*F(1890)+p1mup5*F(1896)+p6mup5*F(1897)+p3mup5*F(19
     -   04)+p2mup5*F(1907))+p1mup5*p6mup2*F(1908)+p2mup5*p6mup2*F(19
     -   09)-p6mup2*p6mup5*F(1910)+p3mup2*p3mup5*F(1911)+p1mup2*(p4mu
     -   p5*F(1889)+p2mup5*F(1895)+p1mup5*F(1905)+p6mup5*F(1906)+p3mu
     -   p5*F(1912))+p3mup5*p6mup2*F(1913))
       Return
       End
