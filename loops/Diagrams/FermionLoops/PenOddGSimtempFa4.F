c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim4(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,
     -   p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mu
     -   p3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup
     -   3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa4)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenOddGSim/F
      Complex*16 EpsF(126)
       common/EPSFPenOddGSim/EPSF
       COMMON/InvariantsPenOddGSim/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s3
     -   4,s45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       Real*8   p1sq 
        Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s234,m1,m1sq,m3,m3sq,m4,m4sq
      Complex*16 Fa(3196),F(4476)
      Real*8 P(1383) 
        COMMON/PPenOddGSim/P
       Complex*16 Fa4(298:396)
       COMMON/FaPenOddGSim/Fa
       Fa4(298) = 4*(p1mup3*p2mup5*p4mup2*F(258)-p1mup5*p2mup3*p4mup2*
     -   F(263)-p2mup3*p3mup5*p4mup2*F(273)-p2mup3*p4mup2*p4mup5*F(27
     -   5)-p2mup5*p3mup3*p4mup2*F(282)-p1mup2*p2mup5*p4mup3*F(940)-p
     -   2mup5*p3mup2*p4mup3*F(1302)-p1mup2*p2mup3*p4mup5*F(1499)-p2m
     -   up3*p2mup5*p4mup2*F(1500)-p2mup5*p4mup2*p4mup3*F(1502))-8*(p
     -   1mup2*p3mup5*p4mup3*F(266)-p1mup2*p3mup3*p4mup5*F(279)+p3mup
     -   5*p4mup2*p4mup3*F(945)+p1mup5*p4mup3*(p4mup2*F(364)+p1mup2*F
     -   (634)+p3mup2*F(947))+p3mup2*p3mup5*p4mup3*F(1025)+p1mup2*p4m
     -   up3*p4mup5*F(1501)-p1mup3*p4mup5*(p3mup2*F(279)-p1mup2*F(882
     -   )+p4mup2*F(1503))+p3mup2*p4mup3*p4mup5*F(1504)+p4mup2*p4mup3
     -   *p4mup5*F(1506)+p3mup2*p3mup3*p4mup5*F(1507)-p3mup3*p4mup2*p
     -   4mup5*F(1508))+mup2mup5*p2mup3*F(1771)+mup2mup3*p2mup5*F(177
     -   2)+mup3mup5*p1mup2*F(1773)+mup3mup5*p4mup2*F(1774)+mup2mup5*
     -   p4mup3*F(1775)+mup2mup3*p4mup5*F(1776)+p2mup3*p3mup2*p4mup5*
     -   F(1778)+2*(mup3mup5*p3mup2*F(1777)+(mup2mup3*p1mup5*F(263)+m
     -   up2mup5*p1mup3*F(264)+mup2mup3*p3mup5*F(273)+mup2mup5*p3mup3
     -   *F(282))*P(7))
       Fa4(299) = p1mup2*(-2*p1mup3*F(1779)+p2mup3*F(1780)+p3mup3*F(17
     -   81))+p1mup3*p3mup2*F(1782)+p2mup3*p3mup2*F(1783)+p3mup2*p3mu
     -   p3*F(1784)-4*p4mup2*p4mup3*F(1785)
       Fa4(300) = 8*(p1mup3*p3mup2*p3mup5*F(770)+p1mup3*p2mup5*p4mup2*
     -   F(1343)+p2mup3*p2mup5*p4mup2*F(1344)+p2mup5*p4mup2*p4mup3*F(
     -   1354)+p3mup2*p3mup3*p4mup5*F(1539)+p3mup3*p3mup5*p4mup2*F(15
     -   40)+p1mup3*p3mup5*p4mup2*F(1541)+p3mup2*p3mup3*p3mup5*F(1542
     -   )+p3mup3*p4mup2*p4mup5*F(1545)+p3mup2*p3mup5*p4mup3*F(1548)+
     -   p3mup5*p4mup2*p4mup3*F(1550))+16*(p2mup3*p4mup2*p4mup5*F(153
     -   )+p1mup3*p3mup2*p4mup5*F(1431)+p1mup3*p4mup2*p4mup5*F(1544)+
     -   p3mup2*p4mup3*p4mup5*F(1551)+p4mup2*p4mup3*p4mup5*F(1553))+m
     -   up2mup5*p2mup3*F(1786)+mup2mup3*p2mup5*F(1787)+mup2mup5*p1mu
     -   p3*F(1789)+mup3mup5*p3mup2*F(1790)+mup2mup5*p3mup3*F(1791)+m
     -   up2mup3*p3mup5*F(1792)-2*(mup2mup5*p4mup3*F(1525)-p1mup2*(2*
     -   (p2mup3*p3mup5*F(1527)+p3mup3*(p2mup5*F(1528)+2*(p3mup5*F(47
     -   )+p1mup5*F(753)+p4mup5*F(1536))))+4*(p2mup5*p4mup3*F(1343)+p
     -   2mup3*p2mup5*F(1425)+p1mup3*(p3mup5*F(751)+p2mup5*F(890)+2*(
     -   p1mup5*F(12)+p4mup5*F(1533)))+p3mup5*p4mup3*F(1541))+8*(p1mu
     -   p5*p2mup3*F(750)+p1mup5*p4mup3*F(758)+p2mup3*p4mup5*F(1534)+
     -   p4mup3*p4mup5*F(1544))+mup3mup5*F(1788))-p1mup5*(2*p2mup3*p3
     -   mup2*F(1530)+4*(p3mup2*p3mup3*F(1546)+p3mup3*p4mup2*F(1547))
     -   +8*(p2mup3*p4mup2*F(61)+p1mup3*p4mup2*F(758)+p1mup3*p3mup2*F
     -   (910)+p4mup2*p4mup3*F(1549)+p3mup2*p4mup3*F(1552))+mup2mup3*
     -   F(1793))-mup3mup5*p4mup2*F(1794)-mup2mup3*p4mup5*F(1795))+4*
     -   (p2mup5*p3mup3*p4mup2*F(1345)+p2mup3*p3mup2*p3mup5*F(1531)+p
     -   2mup3*p2mup5*p3mup2*F(1532)+p2mup3*p3mup2*p4mup5*F(1535)+p2m
     -   up3*p3mup5*p4mup2*F(1537)+p2mup5*p3mup2*p4mup3*F(1538)+p2mup
     -   5*p3mup2*p3mup3*F(1543)+p1mup3*p2mup5*p3mup2*F(1796))
       Fa4(301) = -2*(mup2mup5*p3mup3*F(1560)+mup2mup3*p3mup5*F(1561)+
     -   mup2mup5*p4mup3*F(1563)+mup2mup3*p4mup5*F(1564))+4*(p1mup5*p
     -   2mup3*p4mup2*F(1572)+p1mup2*(4*(p3mup5*p4mup3*F(905)+p4mup3*
     -   p4mup5*F(1448))+2*p1mup5*p4mup3*F(1541)+p1mup3*(2*(p1mup5*F(
     -   79)+p3mup5*F(766)+p4mup5*F(1536))+p2mup5*F(1565))-p1mup5*p2m
     -   up3*F(1566)+p2mup3*p3mup5*F(1567)+p2mup3*p2mup5*F(1568)+p2mu
     -   p3*p4mup5*F(1569)+p2mup5*p4mup3*F(1570)+p3mup3*(2*p1mup5*F(7
     -   70)+4*(p3mup5*F(100)+p4mup5*F(905))+p2mup5*F(1575)))+p1mup5*
     -   p2mup3*p3mup2*F(1577))+8*(p2mup3*(p3mup5*p4mup2*F(98)+p3mup2
     -   *p4mup5*F(98)+p4mup2*p4mup5*F(162)+p2mup5*p4mup2*F(665))+p2m
     -   up5*p3mup3*p4mup2*F(1353)+p2mup5*p3mup2*p4mup3*F(1353)+p2mup
     -   5*p4mup2*p4mup3*F(1355)+p1mup5*p3mup3*p4mup2*F(1548)+p1mup5*
     -   p3mup2*p4mup3*F(1548)+p1mup5*p4mup2*p4mup3*F(1550)+p2mup3*p3
     -   mup2*p3mup5*F(1573)+p2mup3*p2mup5*p3mup2*F(1579)+p1mup5*p3mu
     -   p2*p3mup3*F(1580)+p2mup5*p3mup2*p3mup3*F(1581))+16*(p3mup3*p
     -   3mup5*p4mup2*F(933)+p3mup2*p3mup5*p4mup3*F(933)+p3mup2*p3mup
     -   3*p4mup5*F(933)+p3mup5*p4mup2*p4mup3*F(937)+p3mup3*p4mup2*p4
     -   mup5*F(937)+p3mup2*p4mup3*p4mup5*F(937)+p4mup2*p4mup3*p4mup5
     -   *F(1450)+p3mup2*p3mup3*p3mup5*F(1582))+mup3mup5*p1mup2*F(179
     -   7)+p1mup3*(4*(p2mup5*p4mup2*F(1571)+p2mup5*p3mup2*F(1576))+8
     -   *(p3mup2*p4mup5*F(1539)+p4mup2*p4mup5*F(1545)+p1mup5*(p4mup2
     -   *F(922)+p3mup2*F(1574))+p3mup5*(p4mup2*F(1539)+p3mup2*F(1578
     -   )))+mup2mup5*F(1798))+mup2mup3*p1mup5*F(1799)+mup2mup5*p2mup
     -   3*F(1800)+mup2mup3*p2mup5*F(1801)+mup3mup5*p3mup2*F(1802)+mu
     -   p3mup5*p4mup2*F(1803)
       Fa4(302) = p2mup4*p3mup2*F(1805)+p1mup2*(4*p1mup4*F(1585)+p3mup
     -   4*F(1804)+p2mup4*F(1806)+p4mup4*F(1807))+p2mup4*p4mup2*F(180
     -   8)+p3mup2*p4mup4*F(1810)+p3mup4*p4mup2*F(1811)-2*(p1mup4*p3m
     -   up2*F(1590)+p3mup2*p3mup4*F(1809)-p1mup4*p4mup2*F(1812)-p4mu
     -   p2*p4mup4*F(1813))
       Fa4(303) = p3mup2*(p1mup4*F(1817)+p2mup4*F(1818)-2*(p4mup4*F(16
     -   04)-p3mup4*F(1819)))+p1mup2*(p2mup4*F(1815)+p3mup4*F(1816)+2
     -   *(p1mup4*F(1814)+p4mup4*F(1820)))+p4mup2*(2*(p1mup4*F(1821)+
     -   p2mup4*F(1822))+p3mup4*F(1823)+4*p4mup4*F(1824))
       Fa4(304) = mup1mup2*(p3mup5*F(1254)+2*(p1mup5*F(1259)+p2mup5*F(
     -   1260))+p4mup5*F(1684))
       Fa4(305) = mup1mup2*(p1mup5*F(1280)+2*(p2mup5*F(1285)+p3mup5*F(
     -   1288))+p4mup5*F(1825))
       Fa4(306) = mup1mup2*(p1mup5*F(1699)-2*(p2mup5*F(1703)+p4mup5*F(
     -   1706))+p3mup5*F(1826))
       Fa4(307) = mup1mup3*(p4mup5*F(1363)+p1mup5*F(1365)+p2mup5*F(136
     -   6)+p3mup5*F(1642))
       Fa4(308) = mup1mup3*(p1mup5*F(1649)+2*(p2mup5*F(1654)+p3mup5*F(
     -   1656))+p4mup5*F(1827))
       Fa4(309) = mup1mup3*(p1mup5*F(1390)+2*(p2mup5*F(1395)+p4mup5*F(
     -   1398))+p3mup5*F(1828))
       Fa4(310) = -4*p3mup1*(p3mup4*p4mup5*F(271)+p3mup5*p4mup4*F(272)
     -   +p4mup4*p4mup5*F(276)+p1mup5*(p4mup4*F(262)+2*(p2mup4*F(255)
     -   +p3mup4*F(257)+p1mup4*F(634)))+p1mup4*(p4mup5*F(260)+2*(p3mu
     -   p5*F(266)+p2mup5*F(638)))+p2mup4*p4mup5*F(639)+p2mup5*p4mup4
     -   *F(641)+2*(p2mup4*p3mup5*F(265)+p3mup4*p3mup5*F(270)+p2mup5*
     -   p3mup4*F(637)+p2mup4*p2mup5*F(1264)))+mup4mup5*p3mup1*F(1829
     -   )
       Fa4(311) = p3mup1*(4*(p1mup4*(p3mup5*F(272)+p4mup5*F(298)+p2mup
     -   5*F(641)+p1mup5*F(684))+p1mup5*(p3mup4*F(295)+p4mup4*F(313)+
     -   p2mup4*F(1069)))+mup4mup5*F(1830)+8*(p2mup4*p3mup5*F(304)+p3
     -   mup4*p3mup5*F(306)+p3mup5*p4mup4*F(308)+p2mup4*p4mup5*F(314)
     -   +p4mup4*p4mup5*F(315)+p3mup4*p4mup5*F(316)+p2mup5*p4mup4*F(6
     -   91)+p2mup4*p2mup5*F(692)+p2mup5*p3mup4*F(1831)))
       Fa4(312) = p2mup3*p4mup2*F(1833)+2*(p1mup2*p4mup3*F(1832)+p1mup
     -   3*p4mup2*F(1834))+p3mup2*p4mup3*F(1835)+p3mup3*p4mup2*F(1836
     -   )
       Fa4(313) = p1mup2*p4mup3*F(1837)+p1mup3*p4mup2*F(1838)+p2mup3*p
     -   4mup2*F(1839)+p3mup2*p4mup3*F(1840)+p3mup3*p4mup2*F(1841)
       Fa4(314) = p1mup2*p4mup3*F(1842)+p1mup3*p4mup2*F(1843)+p2mup3*p
     -   4mup2*F(1844)+p3mup2*p4mup3*F(1845)+2*p3mup3*p4mup2*F(1846)
       Fa4(315) = 2*(p1mup2*p4mup3*F(1847)+p1mup3*p4mup2*F(1848)+p2mup
     -   3*p4mup2*F(1849)+p3mup2*p4mup3*F(1850))+p3mup3*p4mup2*F(1851
     -   )
       Fa4(316) = mup2mup4*(2*(p1mup3*F(1852)+p2mup3*F(1853))+p3mup3*F
     -   (1855)+p4mup3*F(1860))+mup2mup3*(p3mup4*F(1856)+2*(p1mup4*F(
     -   1857)+p2mup4*F(1858))+p4mup4*F(1861))+mup3mup4*(p3mup2*F(185
     -   4)+p4mup2*F(1859)+2*p1mup2*F(1862))+4*(p1mup4*p3mup2*p4mup3*
     -   F(1882)-p1mup4*p3mup3*p4mup2*F(1883)-p2mup3*p3mup2*p4mup4*F(
     -   1889)-p2mup3*p3mup4*p4mup2*F(1890)+p2mup4*p3mup2*p4mup3*F(18
     -   91)-p2mup4*p3mup3*p4mup2*F(1892)+p3mup2*p3mup4*p4mup3*F(1893
     -   )-p3mup3*p3mup4*p4mup2*F(1894)+p3mup2*p4mup3*p4mup4*F(1895)-
     -   p3mup4*p4mup2*p4mup3*F(1896)+p3mup2*p3mup3*p4mup4*F(1900)-p1
     -   mup3*(p3mup2*(2*(p3mup4*F(1867)+p1mup4*F(1872)+p2mup4*F(1873
     -   ))+p4mup4*F(1880))+p4mup2*(p3mup4*F(1881)+2*(p1mup4*F(1876)+
     -   p2mup4*F(1887)-p4mup4*F(1901))))+p3mup3*p4mup2*p4mup4*F(1905
     -   )+p1mup2*(p3mup3*p4mup4*F(1878)+p3mup4*p4mup3*F(1879)+2*(p3m
     -   up3*p3mup4*F(1866)+p1mup4*p3mup3*F(1869)+p2mup4*p3mup3*F(187
     -   1)+p1mup3*(2*(p1mup4*F(1863)+p2mup4*F(1864)+p3mup4*F(1868))+
     -   p4mup4*F(1874))+p2mup3*p4mup4*F(1875)+p4mup3*p4mup4*F(1906))
     -   +4*(p1mup4*p2mup3*F(1864)+p2mup3*p2mup4*F(1865)+p2mup3*p3mup
     -   4*F(1870)+p1mup4*p4mup3*F(1907)+p2mup4*p4mup3*F(1908)))-2*(p
     -   1mup4*p2mup3*p4mup2*F(1877)+p2mup3*p2mup4*p4mup2*F(1884)+p2m
     -   up3*p3mup2*p3mup4*F(1885)+p1mup4*p2mup3*p3mup2*F(1886)+p2mup
     -   3*p2mup4*p3mup2*F(1888)-p3mup2*p3mup3*p3mup4*F(1897)-p1mup4*
     -   p3mup2*p3mup3*F(1898)-p2mup4*p3mup2*p3mup3*F(1899)-p1mup4*p4
     -   mup2*p4mup3*F(1902)-p2mup3*p4mup2*p4mup4*F(1903)-p2mup4*p4mu
     -   p2*p4mup3*F(1904)+p4mup2*p4mup3*p4mup4*F(1909)))
       Fa4(317) = mup2mup3*(p2mup4*F(1912)-2*(p1mup4*F(1913)+p3mup4*F(
     -   1914))+p4mup4*F(1917))+mup3mup4*(p1mup2*F(1910)+p4mup2*F(191
     -   5)+2*p3mup2*F(1918))+16*p1mup2*p1mup3*(p1mup4*F(1919)+p2mup4
     -   *F(1920))+mup2mup4*(p2mup3*F(1911)+p4mup3*F(1916)+2*(p1mup3*
     -   F(1931)+p3mup3*F(1954)))+4*(-(p1mup4*p2mup3*p4mup2*F(1930))-
     -   p2mup3*p2mup4*p4mup2*F(1932)-p2mup3*p3mup4*p4mup2*F(1934)+p2
     -   mup3*p4mup2*p4mup4*F(1937)-p2mup4*p4mup2*p4mup3*F(1938)+p2mu
     -   p4*p3mup2*p4mup3*F(1942)-p1mup3*(-(p2mup4*p4mup2*F(1929))+4*
     -   (p1mup4*p4mup2*F(1939)+p3mup4*p4mup2*F(1940)+p3mup2*p3mup4*F
     -   (1945)+p1mup4*p3mup2*F(1947)+p2mup4*p3mup2*F(1948))+2*(p4mup
     -   2*p4mup4*F(1936)+p3mup2*p4mup4*F(1950)))+p1mup2*(p2mup3*p4mu
     -   p4*F(1927)-p2mup4*p4mup3*F(1928)+4*(p1mup3*p3mup4*F(1923)-p3
     -   mup3*p3mup4*F(1945)-p2mup4*p3mup3*F(1948))+2*(p2mup3*p2mup4*
     -   F(1922)+p2mup3*p3mup4*F(1924)+p1mup3*p4mup4*F(1925)-p3mup4*p
     -   4mup3*F(1933)-p4mup3*p4mup4*F(1935)+p1mup4*(p2mup3*F(1921)-p
     -   4mup3*F(1926)-2*p3mup3*F(1947))-p3mup3*p4mup4*F(1950)))+p2mu
     -   p3*p3mup2*p4mup4*F(1952)-p2mup4*p3mup3*p4mup2*F(1953)+4*(p3m
     -   up2*p3mup3*p3mup4*F(1956)+p1mup4*p3mup2*p3mup3*F(1957)+p2mup
     -   4*p3mup2*p3mup3*F(1958)-p1mup4*p3mup3*p4mup2*F(1960)-p3mup3*
     -   p3mup4*p4mup2*F(1961))+2*(p1mup4*p3mup2*p4mup3*F(1941)+p3mup
     -   2*p3mup4*p4mup3*F(1943)-p3mup4*p4mup2*p4mup3*F(1944)-p2mup3*
     -   p3mup2*p3mup4*F(1946)-p1mup4*p2mup3*p3mup2*F(1949)+p2mup3*p2
     -   mup4*p3mup2*F(1951)+p1mup4*p4mup2*p4mup3*F(1955)+p3mup2*p3mu
     -   p3*p4mup4*F(1959)-p3mup2*p4mup3*p4mup4*F(1962)-p3mup3*p4mup2
     -   *p4mup4*F(1963)-p4mup2*p4mup3*p4mup4*F(1964)))
       Fa4(318) = mup3mup4*(p3mup2*F(1969)+2*(p1mup2*F(1965)+p4mup2*F(
     -   1973)))+mup2mup4*(p2mup3*F(1966)+p1mup3*F(1968)+p3mup3*F(197
     -   0)-2*p4mup3*F(1974))+mup2mup3*(p2mup4*F(1967)+p3mup4*F(1971)
     -   -2*(p1mup4*F(1972)-p4mup4*F(1975)))+4*(p1mup4*p2mup3*p3mup2*
     -   F(1985)+p2mup3*p3mup2*p3mup4*F(1987)+p2mup3*p2mup4*p3mup2*F(
     -   1988)+p2mup4*p3mup2*p3mup3*F(1992)+p2mup3*p3mup4*p4mup2*F(19
     -   95)+p2mup4*p3mup2*p4mup3*F(1996)-p2mup3*p3mup2*p4mup4*F(2003
     -   )+p2mup4*p3mup3*p4mup2*F(2004)+p1mup2*(p2mup3*p3mup4*F(1982)
     -   +p2mup4*p3mup3*F(1983)+2*(p2mup3*p2mup4*F(1979)+p1mup4*p3mup
     -   3*F(1981)+p3mup3*p3mup4*F(1986)+p3mup4*p4mup3*F(1994)+p1mup3
     -   *(p2mup4*F(1977)+p3mup4*F(1980)+2*(p1mup4*F(1976)-p4mup4*F(1
     -   999)))-p2mup4*p4mup3*F(2001)-p3mup3*p4mup4*F(2002))+4*(p1mup
     -   4*p2mup3*F(1978)-p1mup4*p4mup3*F(1999)-p2mup3*p4mup4*F(2000)
     -   +p4mup3*p4mup4*F(2006)))-p1mup3*(p2mup4*p3mup2*F(1984)+2*(p3
     -   mup2*p3mup4*F(1990)-p3mup4*p4mup2*F(1994)+p2mup4*p4mup2*F(20
     -   01))+4*(p1mup4*p3mup2*F(1993)+p1mup4*p4mup2*F(1999)-p4mup2*p
     -   4mup4*F(2006)+p3mup2*p4mup4*F(2011)))-2*(p3mup2*p3mup3*p3mup
     -   4*F(1989)+p1mup4*p3mup2*p3mup3*F(1991)-p3mup3*p3mup4*p4mup2*
     -   F(1997)-p2mup3*p2mup4*p4mup2*F(1998)+p1mup4*p3mup3*p4mup2*F(
     -   2002)-p3mup2*p3mup4*p4mup3*F(2005)-p2mup4*p4mup2*p4mup3*F(20
     -   08)-p3mup3*p4mup2*p4mup4*F(2009)-p3mup4*p4mup2*p4mup3*F(2010
     -   )+p3mup2*p3mup3*p4mup4*F(2012))-4*(p1mup4*p2mup3*p4mup2*F(20
     -   00)-p1mup4*p4mup2*p4mup3*F(2006)-p2mup3*p4mup2*p4mup4*F(2007
     -   )+p1mup4*p3mup2*p4mup3*F(2011)-p4mup2*p4mup3*p4mup4*F(2013)-
     -   p3mup2*p4mup3*p4mup4*F(2014)))
       Fa4(319) = mup2mup3*(p1mup4*F(2017)+p2mup4*F(2018)+p3mup4*F(201
     -   9)+p4mup4*F(2022))+mup3mup4*(p1mup2*F(2015)+p4mup2*F(2020)+2
     -   *p3mup2*F(2024))+mup2mup4*(p1mup3*F(2016)+p4mup3*F(2021)+2*(
     -   p2mup3*F(2023)-p3mup3*F(2025)))+4*(-(p1mup4*p2mup3*p4mup2*F(
     -   2037))+p1mup4*p4mup2*p4mup3*F(2043)+p1mup4*p3mup2*p4mup3*F(2
     -   044)+p1mup4*p3mup3*p4mup2*F(2045)+p1mup2*(p1mup3*(2*(p1mup4*
     -   F(2026)+p2mup4*F(2027)+p3mup4*F(2030))+p4mup4*F(2032))-p1mup
     -   4*p4mup3*F(2033)+p2mup3*p4mup4*F(2034)-p2mup4*p4mup3*F(2035)
     -   -p3mup4*p4mup3*F(2039)-p4mup3*p4mup4*F(2042)+2*(p1mup4*p2mup
     -   3*F(2028)+p2mup3*p2mup4*F(2029)+p2mup3*p3mup4*F(2031)-p3mup3
     -   *p3mup4*F(2051)-p1mup4*p3mup3*F(2054)-p2mup4*p3mup3*F(2055))
     -   -p3mup3*p4mup4*F(2059))+p1mup3*(p2mup4*p4mup2*F(2036)+p1mup4
     -   *p4mup2*F(2038)+p3mup4*p4mup2*F(2040)-2*(p3mup2*p3mup4*F(205
     -   2)-p1mup4*p3mup2*F(2056)+p2mup4*p3mup2*F(2057))+p3mup2*p4mup
     -   4*F(2060)-p4mup2*p4mup4*F(2063))+2*(p2mup3*p2mup4*p4mup2*F(2
     -   041)-p2mup3*p3mup4*p4mup2*F(2046)+p2mup4*p3mup2*p4mup3*F(204
     -   7)-p2mup4*p3mup3*p4mup2*F(2048)+p3mup2*p3mup4*p4mup3*F(2049)
     -   -p3mup3*p3mup4*p4mup2*F(2050)-p2mup3*p3mup2*p4mup4*F(2062)+p
     -   2mup3*p4mup2*p4mup4*F(2064)-p2mup4*p4mup2*p4mup3*F(2065)+p3m
     -   up2*p3mup3*p4mup4*F(2069)-p3mup2*p4mup3*p4mup4*F(2070)+p3mup
     -   4*p4mup2*p4mup3*F(2071)-p3mup3*p4mup2*p4mup4*F(2072))-4*(p2m
     -   up3*p3mup2*p3mup4*F(2053)+p1mup4*p2mup3*p3mup2*F(2058)+p2mup
     -   3*p2mup4*p3mup2*F(2061)-p3mup2*p3mup3*p3mup4*F(2066)-p1mup4*
     -   p3mup2*p3mup3*F(2067)-p2mup4*p3mup2*p3mup3*F(2068)-p4mup2*p4
     -   mup3*p4mup4*F(2073)))
       Fa4(320) = mup3mup4*(p1mup2*F(2074)+p3mup2*F(2078)+2*p4mup2*F(2
     -   082))+mup2mup4*(p1mup3*F(2075)+p2mup3*F(2077)+p3mup3*F(2079)
     -   +2*p4mup3*F(2083))+mup2mup3*(p1mup4*F(2076)+p3mup4*F(2080)+2
     -   *(p2mup4*F(2081)+p4mup4*F(2084)))+4*(p1mup4*p2mup3*p3mup2*F(
     -   2094)-p1mup4*p3mup2*p4mup3*F(2102)-p1mup4*p3mup2*p3mup3*F(21
     -   08)+p1mup4*p3mup3*p4mup2*F(2118)+p1mup2*(p1mup4*p3mup3*F(209
     -   0)+p2mup3*p3mup4*F(2091)+p2mup4*p3mup3*F(2092)+p3mup3*p3mup4
     -   *F(2095)+p3mup4*p4mup3*F(2099)+p1mup3*(p3mup4*F(2089)+2*(p1m
     -   up4*F(2085)+p2mup4*F(2086)-p4mup4*F(2110)))-p3mup3*p4mup4*F(
     -   2117)+2*(p1mup4*p2mup3*F(2087)+p2mup3*p2mup4*F(2088)-p1mup4*
     -   p4mup3*F(2111)-p2mup3*p4mup4*F(2112)-p2mup4*p4mup3*F(2113)-p
     -   4mup3*p4mup4*F(2124)))+p1mup3*(-(p2mup4*p3mup2*F(2093))-p3mu
     -   p2*p3mup4*F(2096)-p1mup4*p3mup2*F(2097)+p3mup2*p4mup4*F(2100
     -   )+p3mup4*p4mup2*F(2101)+2*p1mup4*p4mup2*F(2114)-4*(p2mup4*p4
     -   mup2*F(2115)-p4mup2*p4mup4*F(2126)))-2*(p2mup3*p3mup2*p3mup4
     -   *F(2098)-p2mup3*p2mup4*p3mup2*F(2103)+p2mup3*p3mup2*p4mup4*F
     -   (2104)-p2mup3*p3mup4*p4mup2*F(2105)+p3mup2*p3mup4*p4mup3*F(2
     -   106)-p3mup3*p3mup4*p4mup2*F(2107)+p1mup4*p2mup3*p4mup2*F(211
     -   6)-p2mup4*p3mup2*p4mup3*F(2119)+p2mup4*p3mup3*p4mup2*F(2120)
     -   -p2mup4*p3mup2*p3mup3*F(2122)+p3mup2*p3mup3*p4mup4*F(2123)-p
     -   1mup4*p4mup2*p4mup3*F(2125)-p3mup3*p4mup2*p4mup4*F(2128)-p3m
     -   up4*p4mup2*p4mup3*F(2129)+p3mup2*p4mup3*p4mup4*F(2130))+4*(-
     -   (p2mup3*p2mup4*p4mup2*F(2109))+p3mup2*p3mup3*p3mup4*F(2121)+
     -   p2mup4*p4mup2*p4mup3*F(2127)+p2mup3*p4mup2*p4mup4*F(2127)+p4
     -   mup2*p4mup3*p4mup4*F(2131)))
       Fa4(321) = mup3mup4*(p1mup2*F(2132)+p3mup2*F(2137)+p4mup2*F(214
     -   0))+mup2mup4*(p1mup3*F(2133)+p2mup3*F(2135)-2*(p3mup3*F(2138
     -   )+p4mup3*F(2141)))+mup2mup3*(p1mup4*F(2134)+p2mup4*F(2136)+2
     -   *(p3mup4*F(2139)+p4mup4*F(2142)))+4*(-(p1mup3*p2mup4*p3mup2*
     -   F(2156))+p1mup3*p2mup4*p4mup2*F(2163)+p2mup3*(-(p1mup4*p3mup
     -   2*F(2157))+p1mup4*p4mup2*F(2164)+2*(p2mup4*p4mup2*F(2149)-p3
     -   mup2*p3mup4*F(2152)-p2mup4*p3mup2*F(2167)+p3mup4*p4mup2*F(21
     -   68)+p3mup2*p4mup4*F(2175)-p4mup2*p4mup4*F(2180)))-2*(p1mup3*
     -   p3mup2*p3mup4*F(2151)+p1mup3*p1mup4*p3mup2*F(2155)-p1mup3*p1
     -   mup4*p4mup2*F(2162)-p1mup3*p3mup4*p4mup2*F(2166)-p1mup4*p3mu
     -   p2*p3mup3*F(2169)-p2mup4*p3mup2*p3mup3*F(2170)-p1mup3*p3mup2
     -   *p4mup4*F(2173)-p1mup4*p3mup3*p4mup2*F(2174)-p1mup4*p3mup2*p
     -   4mup3*F(2174)-p2mup4*p3mup3*p4mup2*F(2176)-p2mup4*p3mup2*p4m
     -   up3*F(2176)+p1mup3*p4mup2*p4mup4*F(2178)-p1mup4*p4mup2*p4mup
     -   3*F(2179)-p2mup4*p4mup2*p4mup3*F(2181))+p1mup2*(-(p1mup4*p2m
     -   up3*F(2144))+p2mup3*p2mup4*F(2146)+p2mup3*p3mup4*F(2148)+p2m
     -   up4*p3mup3*F(2154)+p1mup3*(p2mup4*F(2143)+2*(p1mup4*F(2145)+
     -   p3mup4*F(2147)+p4mup4*F(2158)))+2*(p1mup4*p3mup3*F(2153)+p1m
     -   up4*p4mup3*F(2159))+p2mup3*p4mup4*F(2160)+p2mup4*p4mup3*F(21
     -   61)-4*(p3mup3*p3mup4*F(2150)-p3mup4*p4mup3*F(2165)-p3mup3*p4
     -   mup4*F(2172)-p4mup3*p4mup4*F(2182)))+4*(-(p3mup2*p3mup3*p3mu
     -   p4*F(2171))+p3mup3*p3mup4*p4mup2*F(2177)+p3mup2*p3mup4*p4mup
     -   3*F(2177)+p3mup4*p4mup2*p4mup3*F(2183)+p3mup2*p3mup3*p4mup4*
     -   F(2184)+p4mup2*p4mup3*p4mup4*F(2185)+(p3mup3*p4mup2+p3mup2*p
     -   4mup3)*p4mup4*F(2186)))
       Fa4(322) = p2mup1*(p1mup4*p4mup5*F(260)+p1mup5*p4mup4*F(262)+p3
     -   mup4*p4mup5*F(271)+p3mup5*p4mup4*F(272)+p4mup4*p4mup5*F(276)
     -   +p2mup4*p4mup5*F(639)+p2mup5*p4mup4*F(641)+2*(p1mup5*p3mup4*
     -   F(257)+p1mup4*p3mup5*F(266)+p3mup4*p3mup5*F(270)+p2mup4*p2mu
     -   p5*F(1264)))
       Fa4(323) = p2mup1*(p1mup4*(p3mup5*F(272)+p4mup5*F(298))+p3mup4*
     -   (2*(p3mup5*F(306)+p4mup5*F(316))-p1mup5*F(607))+p4mup4*(p1mu
     -   p5*F(313)+2*(p3mup5*F(308)+p4mup5*F(315)+p2mup5*F(691)))+2*(
     -   p2mup4*p4mup5*F(314)+p2mup4*p2mup5*F(692)))
       Fa4(324) = mup4mup5*p3mup1*F(2187)+4*(p1mup5*(-(p2mup4*p3mup1)+
     -   p2mup1*p3mup4)*F(263)+p1mup4*p2mup5*p3mup1*F(346)-p1mup4*p2m
     -   up1*p3mup5*F(346)+p2mup1*p3mup4*p3mup5*F(356)+p2mup1*p3mup4*
     -   p4mup5*F(357)-p2mup5*p3mup1*p4mup4*F(362)+p2mup1*p3mup5*p4mu
     -   p4*F(362)-p2mup4*p3mup1*(p3mup5*F(273)+p4mup5*F(357)+p2mup5*
     -   F(1470))+p2mup5*p3mup1*p3mup4*F(2188))
       Fa4(325) = mup1mup3*(p4mup5*F(1217)+p2mup5*F(1772)+2*(p1mup5*F(
     -   263)+p3mup5*F(273))*P(7))
       Fa4(326) = mup1mup2*(p4mup5*F(881)+p2mup5*F(1493)-2*(p1mup5*F(2
     -   63)+p3mup5*F(273))*P(7))
       Fa4(327) = mup1mup3*(p1mup5*F(631)+p3mup5*F(1247)+p4mup5*F(1248
     -   )+p2mup5*F(1464))
       Fa4(328) = mup1mup3*(p1mup5*F(603)-2*(p3mup5*F(885)+p4mup5*F(88
     -   8))+p2mup5*F(2189))
       Fa4(329) = mup1mup2*(p1mup5*F(342)+p3mup5*F(985)+p4mup5*F(987)+
     -   p2mup5*F(1744))
       Fa4(330) = mup1mup2*(p1mup5*F(286)+2*(p3mup5*F(885)+p4mup5*F(88
     -   8))+p2mup5*F(2190))
       Fa4(331) = mup1mup2*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F(3
     -   83)-p3mup5*F(385))
       Fa4(332) = mup1mup2*(p1mup5*F(380)+2*(p4mup5*F(356)+p2mup5*F(38
     -   1)+p3mup5*F(382)))
       Fa4(333) = mup1mup2*(p1mup5*F(383)+2*(p3mup5*F(4465)+p4mup5*F(4
     -   467)+p2mup5*F(4469)))
       Fa4(334) = mup1mup3*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F(3
     -   83)+p3mup5*F(386))
       Fa4(335) = mup1mup3*(p1mup5*F(380)+2*(p4mup5*F(356)+p2mup5*F(38
     -   1)+p3mup5*F(382)))
       Fa4(336) = mup1mup3*(p1mup5*F(383)+2*(p3mup5*F(4465)+p4mup5*F(4
     -   467)+p2mup5*F(4469)))
       Fa4(337) = 4*p1mup3*p4mup2*F(354)+p2mup3*p4mup2*F(2192)+2*(p3mu
     -   p3*p4mup2*F(1636)+p1mup2*p4mup3*F(2191)+p3mup2*p4mup3*F(2193
     -   ))
       Fa4(338) = p1mup2*p4mup3*F(2194)+p1mup3*p4mup2*F(2195)+2*(p2mup
     -   3*p4mup2*F(2196)-p3mup2*p4mup3*F(2197)+p3mup3*p4mup2*F(2198)
     -   )
       Fa4(339) = p1mup3*p4mup2*F(685)+2*(p3mup3*p4mup2*F(306)+p3mup2*
     -   p4mup3*F(306)+p1mup2*p4mup3*F(943))+p2mup3*p4mup2*F(1008)
       Fa4(340) = p1mup3*p4mup2*F(996)+p2mup3*p4mup2*F(1011)+2*(p1mup2
     -   *p4mup3*F(361)+(p3mup3*p4mup2+p3mup2*p4mup3)*F(4473))
       Fa4(341) = mup2mup3*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F(3
     -   83)+p3mup4*F(386))-mup2mup4*(2*(p1mup3*F(376)+p2mup3*F(377))
     -   +p3mup3*F(378)+p4mup3*F(1012))+mup3mup4*(2*p1mup2*F(1262)+p3
     -   mup2*F(1607)+p4mup2*F(1609))
       Fa4(342) = (mup2mup4*p2mup3-mup2mup3*p2mup4)*F(396)+2*mup3mup4*
     -   (p4mup2*F(356)+p3mup2*F(382)+p1mup2*F(1610))
       Fa4(343) = (mup2mup4*p2mup3-mup2mup3*p2mup4)*F(4464)+2*mup3mup4
     -   *(p1mup2*F(357)+p3mup2*F(4465)+p4mup2*F(4467))
       Fa4(344) = mup3mup4*p1mup2*F(380)-mup2mup4*p1mup3*F(380)+mup2mu
     -   p3*p1mup4*F(380)+2*(mup3mup4*p4mup2*F(356)-mup2mup4*p4mup3*F
     -   (356)+mup2mup3*p4mup4*F(356)-mup2mup4*p2mup3*F(381)+mup2mup3
     -   *p2mup4*F(381)+(mup3mup4*p3mup2-mup2mup4*p3mup3+mup2mup3*p3m
     -   up4)*F(382))
       Fa4(345) = mup3mup4*(p1mup2*F(383)+2*(p3mup2*F(4465)+p4mup2*F(4
     -   467)))-mup2mup4*(p1mup3*F(383)+2*(p3mup3*F(4465)+p4mup3*F(44
     -   67)+p2mup3*F(4469)))+mup2mup3*(p1mup4*F(383)+2*(p3mup4*F(446
     -   5)+p4mup4*F(4467)+p2mup4*F(4469)))
       Fa4(346) = mup1mup2*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F(3
     -   83)+p3mup5*F(418))
       Fa4(347) = mup1mup2*(p1mup5*F(418)+2*(p4mup5*F(4465)+p3mup5*F(4
     -   466)+p2mup5*F(4468)))
       Fa4(348) = mup1mup3*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F(3
     -   83)+p3mup5*F(418))
       Fa4(349) = mup1mup3*(p1mup5*F(418)+2*(p4mup5*F(4465)+p3mup5*F(4
     -   466)+p2mup5*F(4468)))
       Fa4(350) = 4*p1mup3*p4mup2*F(364)+p2mup3*p4mup2*F(993)+2*(p1mup
     -   2*p4mup3*F(684)+(p3mup3*p4mup2+p3mup2*p4mup3)*F(995))
       Fa4(351) = p1mup2*p4mup3*F(1001)-p1mup3*p4mup2*F(1002)+2*p2mup3
     -   *p4mup2*F(1004)+4*(p3mup3*p4mup2+p3mup2*p4mup3)*F(4474)
       Fa4(352) = p1mup3*p4mup2*F(995)+p2mup3*p4mup2*F(1010)+2*(p1mup2
     -   *p4mup3*F(945)+(p3mup3*p4mup2+p3mup2*p4mup3)*F(4471))
       Fa4(353) = -(mup2mup4*(2*(p1mup3*F(263)+p2mup3*F(379))+p4mup3*F
     -   (383)+p3mup3*F(418)))+mup2mup3*(2*(p1mup4*F(263)+p2mup4*F(37
     -   9))+p4mup4*F(383)+p3mup4*F(418))+mup3mup4*(2*p1mup2*F(1262)+
     -   p4mup2*F(1609)+p3mup2*F(1622))
       Fa4(354) = (mup2mup4*p2mup3-mup2mup3*p2mup4)*F(4463)+2*mup3mup4
     -   *(p1mup2*F(273)+p4mup2*F(4465)+p3mup2*F(4466))
       Fa4(355) = mup3mup4*(p1mup2*F(418)+2*(p4mup2*F(4465)+p3mup2*F(4
     -   466)))-mup2mup4*(p1mup3*F(418)+2*(p4mup3*F(4465)+p3mup3*F(44
     -   66)+p2mup3*F(4468)))+mup2mup3*(p1mup4*F(418)+2*(p4mup4*F(446
     -   5)+p3mup4*F(4466)+p2mup4*F(4468)))
       Fa4(356) = mup1mup2*(2*(p1mup5*F(263)+p2mup5*F(379))+p3mup5*F(3
     -   86)+p4mup5*F(447))
       Fa4(357) = mup1mup2*(p1mup5*F(443)+2*(p2mup5*F(444)+p3mup5*F(44
     -   5)+p4mup5*F(459)))
       Fa4(358) = mup1mup3*(-2*(p1mup5*F(263)+p2mup5*F(379))+p3mup5*F(
     -   385)+p4mup5*F(485))
       Fa4(359) = mup1mup3*(p1mup5*F(443)+2*(p2mup5*F(444)+p3mup5*F(44
     -   5)+p4mup5*F(459)))
       Fa4(360) = 4*p1mup3*p4mup2*F(354)+p2mup3*p4mup2*F(1634)+2*(p1mu
     -   p2*p4mup3*F(1633)+p3mup2*p4mup3*F(1635)+p3mup3*p4mup2*F(1636
     -   ))
       Fa4(361) = -(p1mup2*p4mup3*F(1624))-2*p3mup2*p4mup3*F(1627)+p4m
     -   up2*(p1mup3*F(1637)+2*(p2mup3*F(1638)+p3mup3*F(1639)))
       Fa4(362) = p1mup3*p4mup2*F(276)+p2mup3*p4mup2*F(1047)+2*(p3mup3
     -   *p4mup2*F(308)+p3mup2*p4mup3*F(308)+p1mup2*p4mup3*F(1048))
       Fa4(363) = mup2mup4*(2*(p1mup3*F(376)+p2mup3*F(377))+p3mup3*F(3
     -   78)+p4mup3*F(448))-mup2mup3*(2*(p1mup4*F(376)+p2mup4*F(377))
     -   +p3mup4*F(378)+p4mup4*F(448))-mup3mup4*(p3mup2*F(1607)-2*p1m
     -   up2*F(1608)+p4mup2*F(1623))
       Fa4(364) = (mup2mup4*p2mup3-mup2mup3*p2mup4)*F(477)+2*mup3mup4*
     -   (p3mup2*F(445)+p4mup2*F(459)+p1mup2*F(1629))
       Fa4(365) = mup3mup4*p1mup2*F(443)-mup2mup4*p1mup3*F(443)+mup2mu
     -   p3*p1mup4*F(443)-2*(mup2mup4*p2mup3*F(444)-mup2mup3*p2mup4*F
     -   (444)-mup3mup4*p3mup2*F(445)+mup2mup4*p3mup3*F(445)-mup2mup3
     -   *p3mup4*F(445)-(mup3mup4*p4mup2-mup2mup4*p4mup3+mup2mup3*p4m
     -   up4)*F(459))
       Fa4(366) = mup1mup2*(-2*(p1mup5*F(263)+p2mup5*F(379))+p3mup5*F(
     -   385)+p4mup5*F(485))
       Fa4(367) = mup1mup3*(2*(p1mup5*F(263)+p2mup5*F(379))+p3mup5*F(3
     -   86)+p4mup5*F(447))
       Fa4(368) = p1mup2*p4mup3*F(1624)+p1mup3*p4mup2*F(1625)+2*(p2mup
     -   3*p4mup2*F(1626)+p3mup2*p4mup3*F(1627)+p3mup3*p4mup2*F(1628)
     -   )
       Fa4(369) = -(mup2mup4*(2*(p1mup3*F(376)+p2mup3*F(377))+p3mup3*F
     -   (378)+p4mup3*F(448)))+mup2mup3*(2*(p1mup4*F(376)+p2mup4*F(37
     -   7))+p3mup4*F(378)+p4mup4*F(448))+mup3mup4*(2*p1mup2*F(1262)+
     -   p3mup2*F(1607)+p4mup2*F(1623))
       Fa4(370) = mup2mup4*(2*(p1mup3*F(263)+p2mup3*F(379))+p4mup3*F(3
     -   83)+p3mup3*F(418))-mup2mup3*(2*(p1mup4*F(263)+p2mup4*F(379))
     -   +p4mup4*F(383)+p3mup4*F(418))-mup3mup4*(-2*p1mup2*F(1608)+p4
     -   mup2*F(1609)+p3mup2*F(1622))
       Fa4(371) = mup1mup2*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F(3
     -   83)+p3mup5*F(386))
       Fa4(372) = mup1mup3*(2*(p1mup5*F(263)+p2mup5*F(379))+p4mup5*F(3
     -   83)-p3mup5*F(385))
       Fa4(373) = -(mup2mup3*(2*(p1mup4*F(263)+p2mup4*F(379))+p4mup4*F
     -   (383)-p3mup4*F(385)))+mup2mup4*(2*(p1mup3*F(376)+p2mup3*F(37
     -   7))+p3mup3*F(378)+p4mup3*F(1012))-mup3mup4*(p3mup2*F(1607)-2
     -   *p1mup2*F(1608)+p4mup2*F(1609))
       Fa4(374) = mup1mup2*(p4mup5*F(1363)+p1mup5*F(1365)+p2mup5*F(136
     -   6)+p3mup5*F(1642))
       Fa4(375) = mup1mup2*(p1mup5*F(1649)+2*(p2mup5*F(1654)+p3mup5*F(
     -   1656))+p4mup5*F(1827))
       Fa4(376) = mup1mup2*(p1mup5*F(1390)+2*(p2mup5*F(1395)+p4mup5*F(
     -   1398))+p3mup5*F(1828))
       Fa4(377) = mup1mup3*(p3mup5*F(1254)+2*(p1mup5*F(1259)+p2mup5*F(
     -   1260))+p4mup5*F(1684))
       Fa4(378) = mup1mup3*(p1mup5*F(1280)+2*(p2mup5*F(1285)+p3mup5*F(
     -   1288))+p4mup5*F(1825))
       Fa4(379) = mup1mup3*(p1mup5*F(1699)-2*(p2mup5*F(1703)+p4mup5*F(
     -   1706))+p3mup5*F(1826))
       Fa4(380) = 4*p3mup1*(p3mup4*p4mup5*F(271)+p3mup5*p4mup4*F(272)+
     -   p4mup4*p4mup5*F(276)+p1mup5*(p4mup4*F(262)+2*(p2mup4*F(255)+
     -   p3mup4*F(257)+p1mup4*F(634)))+p1mup4*(p4mup5*F(260)+2*(p3mup
     -   5*F(266)+p2mup5*F(638)))+p2mup4*p4mup5*F(639)+p2mup5*p4mup4*
     -   F(641)+2*(p2mup4*p3mup5*F(265)+p3mup4*p3mup5*F(270)+p2mup5*p
     -   3mup4*F(637)+p2mup4*p2mup5*F(1264)))+mup4mup5*p3mup1*F(2199)
       Fa4(381) = -(p3mup1*(4*(p1mup5*p4mup4*F(313)-p1mup5*p3mup4*F(60
     -   7)+p1mup4*(p3mup5*F(272)+p4mup5*F(298)+p2mup5*F(641)+p1mup5*
     -   F(684)))+p2mup4*(-4*p1mup5*F(681)+8*(p3mup5*F(304)+p4mup5*F(
     -   314)+p2mup5*F(692)))+8*(p3mup4*p3mup5*F(306)+p3mup5*p4mup4*F
     -   (308)+p4mup4*p4mup5*F(315)+p3mup4*p4mup5*F(316)+p2mup5*p4mup
     -   4*F(691)+p2mup5*p3mup4*F(1831))-mup4mup5*F(2200)))
       Fa4(382) = 2*p1mup2*p4mup3*F(2201)+p2mup3*p4mup2*F(2202)+p1mup3
     -   *p4mup2*F(2203)+p3mup2*p4mup3*F(2204)+p3mup3*p4mup2*F(2205)
       Fa4(383) = p1mup2*p4mup3*F(2206)+p1mup3*p4mup2*F(2207)+p2mup3*p
     -   4mup2*F(2208)+p3mup2*p4mup3*F(2209)+p3mup3*p4mup2*F(2210)
       Fa4(384) = p1mup2*p4mup3*F(2211)+p1mup3*p4mup2*F(2212)+p2mup3*p
     -   4mup2*F(2213)+p3mup2*p4mup3*F(2214)-2*p3mup3*p4mup2*F(2215)
       Fa4(385) = -2*(p1mup2*p4mup3*F(2216)+p1mup3*p4mup2*F(2217)+p2mu
     -   p3*p4mup2*F(2218)-p3mup2*p4mup3*F(2219))+p3mup3*p4mup2*F(222
     -   0)
       Fa4(386) = 4*(p1mup3*p3mup2*p4mup4*F(1880)+p1mup3*p3mup4*p4mup2
     -   *F(1881)-p1mup4*p3mup2*p4mup3*F(1882)+p1mup4*p3mup3*p4mup2*F
     -   (1883)+p2mup3*p3mup2*p4mup4*F(1889)+p2mup3*p3mup4*p4mup2*F(1
     -   890)-p2mup4*p3mup2*p4mup3*F(1891)+p2mup4*p3mup3*p4mup2*F(189
     -   2)-p3mup2*p3mup4*p4mup3*F(1893)+p3mup3*p3mup4*p4mup2*F(1894)
     -   -p3mup2*p4mup3*p4mup4*F(1895)+p3mup4*p4mup2*p4mup3*F(1896)-p
     -   3mup2*p3mup3*p4mup4*F(1900)-p3mup3*p4mup2*p4mup4*F(1905)-p1m
     -   up2*(p3mup3*p4mup4*F(1878)+p3mup4*p4mup3*F(1879)+2*(p3mup3*p
     -   3mup4*F(1866)+p1mup4*p3mup3*F(1869)+p2mup4*p3mup3*F(1871)+p1
     -   mup3*(2*(p1mup4*F(1863)+p2mup4*F(1864)+p3mup4*F(1868))+p4mup
     -   4*F(1874))+p2mup3*p4mup4*F(1875)+p4mup3*p4mup4*F(1906))+4*(p
     -   1mup4*p2mup3*F(1864)+p2mup3*p2mup4*F(1865)+p2mup3*p3mup4*F(1
     -   870)+p1mup4*p4mup3*F(1907)+p2mup4*p4mup3*F(1908))))+8*(p1mup
     -   3*p3mup2*p3mup4*F(1867)+p1mup3*p1mup4*p3mup2*F(1872)+p1mup3*
     -   p2mup4*p3mup2*F(1873)+p1mup3*p1mup4*p4mup2*F(1876)+p1mup4*p2
     -   mup3*p4mup2*F(1877)+p2mup3*p2mup4*p4mup2*F(1884)+p2mup3*p3mu
     -   p2*p3mup4*F(1885)+p1mup4*p2mup3*p3mup2*F(1886)+p1mup3*p2mup4
     -   *p4mup2*F(1887)+p2mup3*p2mup4*p3mup2*F(1888)-p3mup2*p3mup3*p
     -   3mup4*F(1897)-p1mup4*p3mup2*p3mup3*F(1898)-p2mup4*p3mup2*p3m
     -   up3*F(1899)-p1mup3*p4mup2*p4mup4*F(1901)-p1mup4*p4mup2*p4mup
     -   3*F(1902)-p2mup3*p4mup2*p4mup4*F(1903)-p2mup4*p4mup2*p4mup3*
     -   F(1904)+p4mup2*p4mup3*p4mup4*F(1909))+mup3mup4*p3mup2*F(2221
     -   )+mup2mup3*p3mup4*F(2223)+mup3mup4*p4mup2*F(2226)+mup2mup4*(
     -   -2*(p1mup3*F(1852)+p2mup3*F(1853))+p3mup3*F(2222)+p4mup3*F(2
     -   227))+mup2mup3*p4mup4*F(2228)-2*(mup2mup3*p1mup4*F(2224)-mup
     -   2mup3*p2mup4*F(2225)-mup3mup4*p1mup2*F(2229))
       Fa4(387) = 16*(p1mup3*p1mup4*p4mup2*F(1939)+p1mup3*p3mup4*p4mup
     -   2*F(1940)+p1mup3*p3mup2*p3mup4*F(1945)+p1mup3*p1mup4*p3mup2*
     -   F(1947)+p1mup3*p2mup4*p3mup2*F(1948)-p3mup2*p3mup3*p3mup4*F(
     -   1956)-p1mup4*p3mup2*p3mup3*F(1957)-p2mup4*p3mup2*p3mup3*F(19
     -   58)+p1mup4*p3mup3*p4mup2*F(1960)+p3mup3*p3mup4*p4mup2*F(1961
     -   ))+mup2mup4*p2mup3*F(2231)+mup3mup4*p4mup2*F(2234)+mup2mup4*
     -   p4mup3*F(2235)+mup2mup3*(p2mup4*F(2232)+2*(p1mup4*F(1913)+p3
     -   mup4*F(2233))+p4mup4*F(2236))+8*(p1mup3*p4mup2*p4mup4*F(1936
     -   )-p1mup4*p3mup2*p4mup3*F(1941)-p3mup2*p3mup4*p4mup3*F(1943)+
     -   p3mup4*p4mup2*p4mup3*F(1944)+p2mup3*p3mup2*p3mup4*F(1946)+p1
     -   mup4*p2mup3*p3mup2*F(1949)+p1mup3*p3mup2*p4mup4*F(1950)-p1mu
     -   p4*p4mup2*p4mup3*F(1955)-p3mup2*p3mup3*p4mup4*F(1959)+p3mup2
     -   *p4mup3*p4mup4*F(1962)+p3mup3*p4mup2*p4mup4*F(1963)+p4mup2*p
     -   4mup3*p4mup4*F(1964)+p2mup3*p2mup4*p3mup2*F(2241))-4*(p1mup3
     -   *p2mup4*p4mup2*F(1929)-p1mup4*p2mup3*p4mup2*F(1930)-p2mup3*p
     -   2mup4*p4mup2*F(1932)-p2mup3*p3mup4*p4mup2*F(1934)-p2mup4*p4m
     -   up2*p4mup3*F(1938)-p2mup4*p3mup3*p4mup2*F(1953)-p2mup3*p4mup
     -   2*p4mup4*F(2239)-p2mup4*p3mup2*p4mup3*F(2240)-p2mup3*p3mup2*
     -   p4mup4*F(2242))-2*(p1mup2*(2*(p2mup3*p4mup4*F(1927)-p2mup4*p
     -   4mup3*F(1928))-8*(p3mup3*p3mup4*F(1945)+p1mup4*p3mup3*F(1947
     -   )+p2mup4*p3mup3*F(1948))+4*(p1mup4*p2mup3*F(1921)+p2mup3*p2m
     -   up4*F(1922)+p2mup3*p3mup4*F(1924)+p1mup3*(2*(p1mup4*F(1919)+
     -   p2mup4*F(1920)+p3mup4*F(1923))+p4mup4*F(1925))-p1mup4*p4mup3
     -   *F(1926)-p3mup4*p4mup3*F(1933)-p4mup3*p4mup4*F(1935)-p3mup3*
     -   p4mup4*F(1950))-mup3mup4*F(2230))-mup3mup4*p3mup2*F(2237)-mu
     -   p2mup4*p1mup3*F(2238)-mup2mup4*p3mup3*F(2243))
       Fa4(388) = mup3mup4*p3mup2*F(2248)+mup2mup4*(p2mup3*F(2245)+2*(
     -   p4mup3*F(1974)+p1mup3*F(2247))+p3mup3*F(2249))+2*(mup3mup4*p
     -   1mup2*F(2244)+mup3mup4*p4mup2*F(2251))+mup2mup3*(p2mup4*F(22
     -   46)+p3mup4*F(2250)+2*(p1mup4*F(1972)+p4mup4*F(2252)))-8*(p1m
     -   up2*p3mup4*p4mup3*F(1994)+p3mup3*p3mup4*p4mup2*F(1997)-p1mup
     -   3*p2mup4*p4mup2*F(2001)-p1mup2*p2mup4*p4mup3*F(2001)-p1mup4*
     -   p3mup3*p4mup2*F(2002)-p1mup2*p3mup3*p4mup4*F(2002)+p3mup2*p3
     -   mup4*p4mup3*F(2005)+p2mup4*p4mup2*p4mup3*F(2008)+p3mup3*p4mu
     -   p2*p4mup4*F(2009)+p3mup4*p4mup2*p4mup3*F(2010)-p3mup2*p3mup3
     -   *p4mup4*F(2012)-p2mup3*p2mup4*p4mup2*F(2254))-4*(p1mup2*(4*p
     -   1mup4*p2mup3*F(1978)+p2mup3*p3mup4*F(1982)+p2mup4*p3mup3*F(1
     -   983)+2*(p2mup3*p2mup4*F(1979)+p1mup3*(2*p1mup4*F(1976)+p2mup
     -   4*F(1977)+p3mup4*F(1980))+p1mup4*p3mup3*F(1981)+p3mup3*p3mup
     -   4*F(1986)))+p3mup2*(p1mup4*p2mup3*F(1985)+p2mup3*p3mup4*F(19
     -   87)+p2mup3*p2mup4*F(1988)-2*(p3mup3*p3mup4*F(1989)+p1mup4*p3
     -   mup3*F(1991))-p1mup3*(p2mup4*F(1984)+2*p3mup4*F(1990)+4*p1mu
     -   p4*F(1993)))+2*p1mup3*p3mup4*p4mup2*F(1994)+p2mup3*p3mup4*p4
     -   mup2*F(1995)+p2mup4*p3mup2*p4mup3*F(1996)-p2mup3*p3mup2*p4mu
     -   p4*F(2003)-p2mup4*p3mup2*p3mup3*F(2253)-p2mup4*p3mup3*p4mup2
     -   *F(2255))+16*(p1mup3*p1mup4*p4mup2*F(1999)+p1mup2*p1mup4*p4m
     -   up3*F(1999)+p1mup2*p1mup3*p4mup4*F(1999)+p1mup4*p2mup3*p4mup
     -   2*F(2000)+p1mup2*p2mup3*p4mup4*F(2000)-p1mup4*p4mup2*p4mup3*
     -   F(2006)-p1mup3*p4mup2*p4mup4*F(2006)-p1mup2*p4mup3*p4mup4*F(
     -   2006)-p2mup3*p4mup2*p4mup4*F(2007)+p1mup4*p3mup2*p4mup3*F(20
     -   11)+p1mup3*p3mup2*p4mup4*F(2011)+p4mup2*p4mup3*p4mup4*F(2256
     -   )+p3mup2*p4mup3*p4mup4*F(2257))
       Fa4(389) = mup3mup4*p1mup2*F(2258)+mup2mup3*p1mup4*F(2260)+mup2
     -   mup3*p3mup4*F(2262)+mup3mup4*p4mup2*F(2263)+mup2mup3*p4mup4*
     -   F(2265)+mup2mup4*(p1mup3*F(2259)+p4mup3*F(2264)+2*(p3mup3*F(
     -   2025)+p2mup3*F(2266)))+2*mup3mup4*p3mup2*F(2267)+p2mup4*(16*
     -   p2mup3*p3mup2*F(2061)+mup2mup3*F(2261)+8*(p1mup3*p3mup2*F(20
     -   57)+p4mup2*p4mup3*F(2065)+p3mup3*(p4mup2*F(2048)+p1mup2*F(20
     -   55)-2*p3mup2*F(2068))+p2mup3*p4mup2*F(2268)))+8*(p3mup3*p3mu
     -   p4*p4mup2*F(2050)+p1mup2*p3mup3*p3mup4*F(2051)+p1mup3*p3mup2
     -   *p3mup4*F(2052)+p1mup2*p1mup4*p3mup3*F(2054)+p2mup3*p3mup2*p
     -   4mup4*F(2062)-p2mup3*p4mup2*p4mup4*F(2064)-p3mup2*p3mup3*p4m
     -   up4*F(2069)+p3mup2*p4mup3*p4mup4*F(2070)-p3mup4*p4mup2*p4mup
     -   3*F(2071)+p3mup3*p4mup2*p4mup4*F(2072)+p1mup3*p1mup4*p3mup2*
     -   F(2271))-4*(-(p1mup4*p2mup3*p4mup2*F(2037))+p1mup3*p4mup2*(p
     -   2mup4*F(2036)+p1mup4*F(2038)+p3mup4*F(2040))+p1mup2*(2*(p1mu
     -   p4*p2mup3*F(2028)+p2mup3*p2mup4*F(2029)+p2mup3*p3mup4*F(2031
     -   ))+p1mup3*(2*(p1mup4*F(2026)+p2mup4*F(2027)+p3mup4*F(2030))+
     -   p4mup4*F(2032))-p1mup4*p4mup3*F(2033)+p2mup3*p4mup4*F(2034)-
     -   p2mup4*p4mup3*F(2035)-p3mup4*p4mup3*F(2039)-p4mup3*p4mup4*F(
     -   2042))+p1mup4*p4mup2*p4mup3*F(2043)-2*(p2mup3*p3mup4*p4mup2*
     -   F(2046)-p2mup4*p3mup2*p4mup3*F(2047)-p3mup2*p3mup4*p4mup3*F(
     -   2049))-p1mup2*p3mup3*p4mup4*F(2059)-p1mup3*p4mup2*p4mup4*F(2
     -   063)-p1mup4*p3mup2*p4mup3*F(2269)-p1mup4*p3mup3*p4mup2*F(227
     -   0)-p1mup3*p3mup2*p4mup4*F(2272))+16*(p2mup3*p3mup2*p3mup4*F(
     -   2053)+p1mup4*p2mup3*p3mup2*F(2058)-p3mup2*p3mup3*p3mup4*F(20
     -   66)-p1mup4*p3mup2*p3mup3*F(2067)+p4mup2*p4mup3*p4mup4*F(2273
     -   ))
       Fa4(390) = mup3mup4*p1mup2*F(2274)+mup3mup4*p3mup2*F(2278)+mup2
     -   mup4*(-2*p4mup3*F(2083)+p1mup3*F(2275)+p2mup3*F(2277)+p3mup3
     -   *F(2279))+2*mup3mup4*p4mup2*F(2281)+mup2mup3*(p1mup4*F(2276)
     -   +p3mup4*F(2280)-2*(p2mup4*F(2081)-p4mup4*F(2282)))+8*(p2mup3
     -   *p3mup2*p4mup4*F(2104)-p2mup3*p3mup4*p4mup2*F(2105)+p3mup2*p
     -   3mup4*p4mup3*F(2106)-p3mup3*p3mup4*p4mup2*F(2107)+p1mup2*p1m
     -   up3*p4mup4*F(2110)+p1mup2*p1mup4*p4mup3*F(2111)+p1mup2*p2mup
     -   3*p4mup4*F(2112)+p1mup2*p2mup4*p4mup3*F(2113)+p1mup4*p2mup3*
     -   p4mup2*F(2116)-p2mup4*p3mup2*p4mup3*F(2119)+p2mup4*p3mup3*p4
     -   mup2*F(2120)-p2mup4*p3mup2*p3mup3*F(2122)+p3mup2*p3mup3*p4mu
     -   p4*F(2123)+p1mup2*p4mup3*p4mup4*F(2124)-p1mup4*p4mup2*p4mup3
     -   *F(2125)-p3mup3*p4mup2*p4mup4*F(2128)-p3mup4*p4mup2*p4mup3*F
     -   (2129)+p3mup2*p4mup3*p4mup4*F(2130)+p2mup3*p2mup4*p3mup2*F(2
     -   284)+p1mup3*p1mup4*p4mup2*F(2285))-4*(p1mup2*(2*(p1mup3*p1mu
     -   p4*F(2085)+p1mup3*p2mup4*F(2086)+p1mup4*p2mup3*F(2087)+p2mup
     -   3*p2mup4*F(2088))+p1mup3*p3mup4*F(2089)+p1mup4*p3mup3*F(2090
     -   )+p2mup3*p3mup4*F(2091)+p2mup4*p3mup3*F(2092)+p3mup3*p3mup4*
     -   F(2095))-p3mup2*(p1mup3*p2mup4*F(2093)-p1mup4*p2mup3*F(2094)
     -   +p1mup3*p3mup4*F(2096)+p1mup3*p1mup4*F(2097)+2*p2mup3*p3mup4
     -   *F(2098))+p1mup2*p3mup4*p4mup3*F(2099)+p1mup3*p3mup4*p4mup2*
     -   F(2101)-p1mup4*p3mup2*p4mup3*F(2102)-p1mup4*p3mup2*p3mup3*F(
     -   2108)-p1mup2*p3mup3*p4mup4*F(2117)-p1mup3*p3mup2*p4mup4*F(22
     -   83)-p1mup4*p3mup3*p4mup2*F(2286))+16*(p2mup3*p2mup4*p4mup2*F
     -   (2109)+p1mup3*p2mup4*p4mup2*F(2115)-p3mup2*p3mup3*p3mup4*F(2
     -   121)-p1mup3*p4mup2*p4mup4*F(2126)-p2mup4*p4mup2*p4mup3*F(212
     -   7)-p2mup3*p4mup2*p4mup4*F(2127)+p4mup2*p4mup3*p4mup4*F(2287)
     -   )
       Fa4(391) = mup3mup4*p1mup2*F(2288)+mup2mup4*(2*(p3mup3*F(2138)+
     -   p4mup3*F(2141))+p1mup3*F(2289)+p2mup3*F(2291))+mup3mup4*p3mu
     -   p2*F(2293)+mup2mup3*(p1mup4*F(2290)+p2mup4*F(2292)-2*(p4mup4
     -   *F(2142)-p3mup4*F(2294)))+mup3mup4*p4mup2*F(2295)-4*(-(p1mup
     -   3*p2mup4*p3mup2*F(2156))+p1mup3*p2mup4*p4mup2*F(2163)+p1mup2
     -   *(-(p1mup4*p2mup3*F(2144))+p2mup3*p2mup4*F(2146)+p1mup3*(p2m
     -   up4*F(2143)+2*(p1mup4*F(2145)+p3mup4*F(2147)))+p2mup3*p3mup4
     -   *F(2148)+2*p1mup4*p4mup3*F(2159)+p2mup4*p4mup3*F(2161)-4*(p3
     -   mup3*p3mup4*F(2150)-p3mup4*p4mup3*F(2165)))+p2mup3*(-(p1mup4
     -   *p3mup2*F(2157))+p1mup4*p4mup2*F(2164)+2*(p2mup4*p4mup2*F(21
     -   49)-p3mup2*p3mup4*F(2152)-p2mup4*p3mup2*F(2167)+p3mup4*p4mup
     -   2*F(2168)))-2*(p1mup3*p3mup2*p3mup4*F(2151)+p1mup3*p1mup4*p3
     -   mup2*F(2155)-p1mup3*p1mup4*p4mup2*F(2162)-p1mup3*p3mup4*p4mu
     -   p2*F(2166)-p1mup4*p3mup2*p3mup3*F(2169)-p2mup4*p3mup2*p3mup3
     -   *F(2170))-4*p3mup2*p3mup3*p3mup4*F(2171)-p1mup2*p2mup4*p3mup
     -   3*F(2297)-p1mup2*p2mup3*p4mup4*F(2299))-8*(p1mup4*p3mup3*p4m
     -   up2*F(2174)+p1mup4*p3mup2*p4mup3*F(2174)+p2mup4*p3mup3*p4mup
     -   2*F(2176)+p2mup4*p3mup2*p4mup3*F(2176)-p1mup3*p4mup2*p4mup4*
     -   F(2178)+p1mup4*p4mup2*p4mup3*F(2179)-p2mup3*p4mup2*p4mup4*F(
     -   2180)+p2mup4*p4mup2*p4mup3*F(2181)-p1mup2*p1mup4*p3mup3*F(22
     -   96)-p1mup2*p1mup3*p4mup4*F(2298)-p1mup3*p3mup2*p4mup4*F(2300
     -   )-p2mup3*p3mup2*p4mup4*F(2301))-16*(p1mup2*p3mup3*p4mup4*F(2
     -   172)+p3mup3*p3mup4*p4mup2*F(2177)+p3mup2*p3mup4*p4mup3*F(217
     -   7)+p1mup2*p4mup3*p4mup4*F(2182)+p3mup4*p4mup2*p4mup3*F(2183)
     -   -p3mup2*p3mup3*p4mup4*F(2302)-p4mup2*p4mup3*p4mup4*F(2303)-(
     -   p3mup3*p4mup2+p3mup2*p4mup3)*p4mup4*F(2304))
       Fa4(392) = p2mup1*(p1mup4*(p3mup5*F(272)+p4mup5*F(298))+p1mup5*
     -   (p3mup4*F(295)+p4mup4*F(313))+2*(p3mup4*p3mup5*F(306)+p3mup5
     -   *p4mup4*F(308)+p2mup4*p4mup5*F(314)+p4mup4*p4mup5*F(315)+p3m
     -   up4*p4mup5*F(316)+p2mup5*p4mup4*F(691)+p2mup4*p2mup5*F(692))
     -   )
       Fa4(393) = 4*(p1mup5*(p2mup4*p3mup1-p2mup1*p3mup4)*F(263)-p1mup
     -   4*(p2mup5*p3mup1-p2mup1*p3mup5)*F(346)-p2mup1*p3mup4*p3mup5*
     -   F(356)-p2mup1*p3mup4*p4mup5*F(357)+p2mup5*p3mup1*p4mup4*F(36
     -   2)-p2mup1*p3mup5*p4mup4*F(362)+p2mup4*p3mup1*(p3mup5*F(273)+
     -   p4mup5*F(357)+p2mup5*F(1470))-p2mup5*p3mup1*p3mup4*F(2188))+
     -   mup4mup5*p3mup1*F(2305)
       Fa4(394) = mup1mup3*(p4mup5*F(881)+p2mup5*F(1493)-2*(p1mup5*F(2
     -   63)+p3mup5*F(273))*P(7))
       Fa4(395) = mup1mup2*(p4mup5*F(1217)+p2mup5*F(1772)+2*(p1mup5*F(
     -   263)+p3mup5*F(273))*P(7))
       Fa4(396) = mup1mup3*(p1mup5*F(342)+p3mup5*F(985)+p4mup5*F(987)+
     -   p2mup5*F(1744))
       Return
       End
