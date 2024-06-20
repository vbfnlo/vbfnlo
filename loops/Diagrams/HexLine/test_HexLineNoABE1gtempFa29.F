c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g29(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa29)
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
       Complex*16 Fa29(393:406)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa29(393) = mup4mup6*F(1755)+p3mup4*(p1mup6*F(1757)+p2mup6*F(17
     -   58)-4*(p5mup6*F(1756)+p3mup6*F(1759)+p4mup6*F(1760)))+p2mup4
     -   *p2mup6*F(1765)+4*(p2mup4*p5mup6*F(1762)+p1mup6*p2mup4*F(176
     -   4)+p2mup4*p3mup6*F(1766)+p2mup4*p4mup6*F(1767)+p4mup4*p5mup6
     -   *F(1768)-p1mup6*p4mup4*F(1769)+p2mup6*p4mup4*F(1770)+p3mup6*
     -   p4mup4*F(1771)+p4mup4*p4mup6*F(1772)+p5mup4*p5mup6*F(1773)-p
     -   1mup6*p5mup4*F(1774)+p2mup6*p5mup4*F(1775)+p3mup6*p5mup4*F(1
     -   776)+p4mup6*p5mup4*F(1777)-p1mup4*(p3mup6*F(1761)+p2mup6*F(1
     -   763)+p5mup6*F(1778)+p1mup6*F(1779)+p4mup6*F(1780)))
       Fa29(394) = mup3mup6*F(1781)+p2mup3*p2mup6*F(1785)+p1mup3*p2mup
     -   6*F(1792)+p2mup6*p3mup3*F(1793)+p2mup6*p4mup3*F(1794)+p5mup3
     -   *p5mup6*F(1800)+4*(p2mup3*p4mup6*F(1786)+p1mup3*p3mup6*F(178
     -   7)+p2mup3*p3mup6*F(1789)+p2mup3*p5mup6*F(1790)+p3mup3*(p3mup
     -   6*F(1782)+p4mup6*F(1783)+p1mup6*F(1788)+p5mup6*F(1791))-p3mu
     -   p6*p5mup3*F(1795)-p4mup6*p5mup3*F(1796)-p2mup6*p5mup3*F(1798
     -   )-p4mup3*p5mup6*F(1799)+p1mup3*p4mup6*F(1802)+p1mup3*p5mup6*
     -   F(1803)-p4mup3*p4mup6*F(1804)+p3mup6*p4mup3*F(1805))+p1mup6*
     -   (p2mup3*F(1784)+p1mup3*F(1801)-4*(p5mup3*F(1797)-p4mup3*F(18
     -   06)))
       Fa29(395) = mup2mup6*F(1807)+p2mup2*p3mup6*F(1813)+p2mup6*p3mup
     -   2*F(1816)+p1mup6*p2mup2*F(1818)+p4mup2*p4mup6*F(1820)+p3mup6
     -   *p4mup2*F(1822)+p1mup6*p4mup2*F(1823)-4*(p3mup2*(p4mup6*F(18
     -   08)+p5mup6*F(1809)-p1mup6*F(1810)+p3mup6*F(1811))-p2mup2*p4m
     -   up6*F(1814)-p2mup2*p5mup6*F(1815)-p2mup2*p2mup6*F(1819)-p2mu
     -   p6*p4mup2*F(1821)-p4mup2*p5mup6*F(1824)-p4mup6*p5mup2*F(1825
     -   )-p2mup6*p5mup2*F(1826)+p1mup6*p5mup2*F(1827)-p5mup2*p5mup6*
     -   F(1828)-p3mup6*p5mup2*F(1829)-p1mup2*(p3mup6*F(1812)+p2mup6*
     -   F(1817)-p4mup6*F(1830)-p5mup6*F(1831)+p1mup6*F(1832)))
       Fa29(396) = mup3mup6*p5mup4*F(1838)+mup4mup6*(p4mup3*F(1835)+p1
     -   mup3*F(1836)+4*(p2mup3*F(1834)+p3mup3*F(1837))+p5mup3*F(1839
     -   ))
       Fa29(397) = mup3mup6*(p3mup4*F(1840)+p4mup4*F(1841))
       Fa29(398) = 4*p5mup6*(mup3mup4*F(1842)+2*(-(p3mup4*p4mup3*F(184
     -   8))+p2mup4*p4mup3*F(1850)+p4mup3*p4mup4*F(1851)+p3mup4*p5mup
     -   3*F(1854)+p4mup3*p5mup4*F(1855)+p4mup4*p5mup3*F(1856)-p2mup4
     -   *p5mup3*F(1857)+p5mup3*p5mup4*F(1858)+p1mup3*p1mup4*F(1859)+
     -   p1mup3*p2mup4*F(1860)+p2mup3*(p3mup4*F(1844)+p4mup4*F(1846)+
     -   p2mup4*F(1849)+p5mup4*F(1852)+p1mup4*F(1861))+p1mup3*p3mup4*
     -   F(1862)+p3mup3*(p3mup4*F(1843)+p2mup4*F(1845)+p4mup4*F(1847)
     -   +p5mup4*F(1853)+p1mup4*F(1863))+p1mup3*p4mup4*F(1864)-p1mup4
     -   *p4mup3*F(1865)+p1mup3*p5mup4*F(1866)+p1mup4*p5mup3*F(1867))
     -   )
       Fa29(399) = 8*mup3mup6*(p1mup4*F(1868)+p2mup4*F(1869))+p5mup4*(
     -   p4mup6*p5mup3*F(1877)+8*(p1mup3*p4mup6*F(1871)+p4mup3*(p4mup
     -   6*F(1870)+p2mup6*F(1872)+p3mup6*F(1873)+p1mup6*F(1874))+p2mu
     -   p3*p4mup6*F(1875)+p3mup3*p4mup6*F(1876)-p2mup6*p5mup3*F(1878
     -   ))+p3mup6*p5mup3*F(1879))
       Fa29(400) = 4*mup3mup4*(p3mup6*F(1880)+p4mup6*F(1881)+p1mup6*F(
     -   1882)+p2mup6*F(1883))+p2mup6*p3mup3*p4mup4*F(1888)+p2mup4*p3
     -   mup6*p5mup3*F(1912)+p3mup4*p3mup6*p5mup3*F(1913)+p3mup6*p4mu
     -   p4*p5mup3*F(1914)+p2mup4*p4mup6*p5mup3*F(1915)+p3mup4*p4mup6
     -   *p5mup3*F(1916)+p1mup6*p5mup3*p5mup4*F(1919)+8*(p2mup4*p2mup
     -   6*p3mup3*F(1886)+p2mup6*p3mup3*p3mup4*F(1887)+p2mup6*p3mup4*
     -   p4mup3*F(1892)+p1mup4*p3mup3*p4mup6*F(1893)+p1mup4*p3mup6*p4
     -   mup3*F(1894)+p3mup3*p4mup4*p4mup6*F(1895)-p3mup4*p3mup6*p4mu
     -   p3*F(1896)+p1mup4*p2mup6*p4mup3*F(1901)+p2mup4*p2mup6*p4mup3
     -   *F(1903)+p2mup4*p3mup3*p4mup6*F(1904)+p2mup4*p3mup6*p4mup3*F
     -   (1905)+p3mup3*p3mup4*p4mup6*F(1906)+p2mup6*p4mup3*p4mup4*F(1
     -   907)-p1mup4*p2mup6*p5mup3*F(1908)-p2mup4*p2mup6*p5mup3*F(190
     -   9)-p2mup6*p3mup4*p5mup3*F(1910)-p2mup6*p4mup4*p5mup3*F(1911)
     -   -p4mup4*p4mup6*p5mup3*F(1917)+p3mup6*p4mup3*p4mup4*F(1918)+p
     -   1mup3*p1mup4*p1mup6*F(1920)+p1mup3*p1mup6*p2mup4*F(1921)+p1m
     -   up3*p1mup6*p3mup4*F(1922)+p1mup3*p1mup6*p4mup4*F(1923)+p1mup
     -   4*p2mup6*p3mup3*F(1925)+p2mup3*(p2mup4*p3mup6*F(1884)+p3mup4
     -   *p3mup6*F(1885)+p1mup4*p4mup6*F(1889)+p2mup4*p4mup6*F(1890)+
     -   p3mup4*p4mup6*F(1891)+p2mup4*p2mup6*F(1897)+p2mup6*p3mup4*F(
     -   1898)+p2mup6*p4mup4*F(1899)+p3mup6*p4mup4*F(1900)+p4mup4*p4m
     -   up6*F(1902)+p1mup4*p3mup6*F(1924)+p1mup4*p2mup6*F(1926))+p1m
     -   up4*p3mup3*p3mup6*F(1927)+p2mup4*p3mup3*p3mup6*F(1928)+p3mup
     -   3*p3mup4*p3mup6*F(1929))+p1mup4*p3mup6*p5mup3*F(1930)+p1mup4
     -   *p4mup6*p5mup3*F(1931)
       Fa29(401) = p1mup6*p2mup3*p2mup4*F(1934)+p1mup6*p2mup3*p3mup4*F
     -   (1935)+p3mup3*(8*p3mup6*p4mup4*F(1932)+p1mup6*(p2mup4*F(1940
     -   )+8*(p3mup4*F(1941)+p4mup4*F(1942)+p1mup4*F(1960))))+8*(p3mu
     -   p4*p4mup3*p4mup6*F(1933)+p1mup3*(p3mup4*p3mup6*F(1938)+p3mup
     -   6*p4mup4*F(1939)+p2mup6*p3mup4*F(1946)+p2mup6*p4mup4*F(1947)
     -   +p2mup4*(p3mup6*F(1937)+p2mup6*F(1945)+p4mup6*F(1949))+p3mup
     -   4*p4mup6*F(1950)+p4mup4*p4mup6*F(1951)+p1mup4*p2mup6*F(1957)
     -   +p1mup4*p3mup6*F(1959)+p1mup4*p4mup6*F(1961))+p4mup3*p4mup6*
     -   (p2mup4*F(1944)+p1mup4*F(1956)+p4mup4*F(1963)))+p1mup6*(p2mu
     -   p4*p4mup3*F(1948)+p2mup4*p5mup3*F(1953)+p3mup4*p5mup3*F(1954
     -   )+p4mup4*p5mup3*F(1955)+p1mup4*p2mup3*F(1958)+8*(p2mup3*p4mu
     -   p4*F(1936)+p3mup4*p4mup3*F(1943)+p4mup3*p4mup4*F(1952)+p1mup
     -   4*p4mup3*F(1962))+p1mup4*p5mup3*F(1964))
       Fa29(402) = p5mup4*(p2mup6*(p3mup3*F(1965)+p2mup3*F(1966)+p1mup
     -   3*F(1967))+p3mup6*(p1mup3*F(1969)+p2mup3*F(1971)+p3mup3*F(19
     -   72))+p1mup6*(p2mup3*F(1968)+p3mup3*F(1970)+p1mup3*F(1973)))
       Fa29(403) = mup4mup6*(p1mup2*F(1974)+p2mup2*F(1975)+p3mup2*F(19
     -   76)+p4mup2*F(1977)+p5mup2*F(1978))+mup2mup6*p5mup4*F(1979)
       Fa29(404) = p5mup6*(p2mup2*(p3mup4*F(1980)-p4mup4*F(1981)-p5mup
     -   4*F(1982))+p1mup2*(p3mup4*F(1983)+p4mup4*F(1984)-p5mup4*F(19
     -   85)))
       Fa29(405) = p2mup2*(p3mup4*(-(p3mup6*F(1987))+p1mup6*F(1991)+p2
     -   mup6*F(1992)-p4mup6*F(1993))+p1mup6*p4mup4*F(1995)+p2mup6*p4
     -   mup4*F(1996)+p4mup4*p4mup6*F(1997)+p1mup6*p5mup4*F(2003)-p2m
     -   up6*p5mup4*F(2004)-p3mup6*(p4mup4*F(1989)+p5mup4*F(2005))-p4
     -   mup6*p5mup4*F(2006)+p1mup4*p5mup6*F(2009)+p2mup4*p5mup6*F(20
     -   10))+p1mup2*(p3mup6*p4mup4*F(1988)-p2mup6*p4mup4*F(1994)+p4m
     -   up4*p4mup6*F(1999)-p2mup6*p5mup4*F(2000)-p3mup6*p5mup4*F(200
     -   1)-p4mup6*p5mup4*F(2002)+p1mup4*p5mup6*F(2007)+p2mup4*p5mup6
     -   *F(2008)+p3mup4*(p3mup6*F(1986)-p2mup6*F(1990)-p4mup6*F(1998
     -   )+p1mup6*F(2011))+p1mup6*p4mup4*F(2012)-p1mup6*p5mup4*F(2013
     -   ))
       Fa29(406) = p1mup2*(p1mup4*(p1mup6*F(2030)+p2mup6*F(2031)+p3mup
     -   6*F(2032)+p4mup6*F(2033))-p2mup4*(p3mup6*F(2014)+p2mup6*F(20
     -   16)+p4mup6*F(2019)-p1mup6*F(2034)))+p2mup2*(p1mup4*(p2mup6*F
     -   (2017)+p3mup6*F(2018)+p4mup6*F(2020)+p1mup6*F(2035))+p2mup4*
     -   (-(p3mup6*F(2015))+p1mup6*F(2036)+p2mup6*F(2037)-p4mup6*F(20
     -   38)))+p5mup4*(-(p3mup2*(p1mup6*F(2024)+p5mup6*F(2026)+p4mup6
     -   *F(2027)))-p3mup6*(p4mup2*F(2023)+p5mup2*F(2028)+p3mup2*F(20
     -   29))+p2mup6*(p3mup2*F(2021)-p4mup2*F(2022)+p5mup2*F(2039))+p
     -   1mup6*(p4mup2*F(2025)+p5mup2*F(2040)))
       Return
       End
