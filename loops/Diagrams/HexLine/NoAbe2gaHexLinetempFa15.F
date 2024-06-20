c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa15(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa15)
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
       Complex*16 Fa15(337:360)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa15(337) = p4mup6*(mup3mup5*F(1811)-p3mup3*p5mup5*F(1814)+p4mu
     -   p3*(p5mup5*F(1812)-p4mup5*F(1817)-2*(p3mup5*F(1816)+p2mup5*F
     -   (1819)))-p2mup3*p5mup5*F(1820)+p5mup3*(p5mup5*F(1813)+p4mup5
     -   *F(1818)-2*(p3mup5*F(1815)+p2mup5*F(1821)))+p2mup3*p4mup5*F(
     -   1822)-p3mup3*p4mup5*F(1823)+2*(p2mup3*p3mup5*F(1824)-p2mup5*
     -   p3mup3*F(1825)+p2mup3*p2mup5*F(1826)+p3mup3*p3mup5*F(1827)))
       Fa15(338) = 4*p4mup6*(p1mup3*(-(p5mup5*F(1820))+p4mup5*F(1822)+
     -   2*p3mup5*F(1824))-2*(p1mup5*(p4mup3*F(1819)+p5mup3*F(1821)+p
     -   3mup3*F(1825))-(p1mup5*p2mup3+p1mup3*(p1mup5+p2mup5))*F(1826
     -   )))
       Fa15(339) = -(p4mup3*(-2*p3mup6*p5mup5*F(1832)-p5mup5*p5mup6*F(
     -   1833)+p4mup5*(p5mup6*F(1830)+p3mup6*F(1837)+p1mup6*F(1839))+
     -   4*(p3mup5*p3mup6*F(1838)+p2mup5*p3mup6*F(1842)+p2mup5*p5mup6
     -   *F(1843)+p3mup5*p5mup6*F(1844))))+p2mup6*(p4mup3*(2*p5mup5*F
     -   (1829)-p4mup5*F(1839)-4*(p2mup5*F(1834)+p3mup5*F(1842)))-p4m
     -   up5*(p5mup3*F(1841)+p3mup3*F(1846))+(p1mup3+p2mup3)*p4mup5*F
     -   (1847))+p4mup5*(-(p5mup3*(p5mup6*F(1831)+p3mup6*F(1836)+p1mu
     -   p6*F(1841)))+(p1mup3+p2mup3)*(p5mup6*F(1840)+p3mup6*F(1845)+
     -   p1mup6*F(1847))+p3mup3*(p5mup6*F(1835)-p1mup6*F(1846)+p3mup6
     -   *F(1848)))
       Fa15(340) = p1mup3*p5mup5*p5mup6*F(1849)+p2mup3*p5mup5*p5mup6*F
     -   (1849)+p5mup3*p5mup5*p5mup6*F(1852)+4*(p1mup3*p2mup5*p2mup6*
     -   F(1853)+p2mup3*p2mup5*p2mup6*F(1853)-p1mup3*p3mup5*p5mup6*F(
     -   1854)-p2mup3*p3mup5*p5mup6*F(1854)-p1mup3*p2mup5*p5mup6*F(18
     -   55)-p2mup3*p2mup5*p5mup6*F(1855)-p1mup3*p3mup5*p3mup6*F(1857
     -   )-p2mup3*p3mup5*p3mup6*F(1857)+p1mup3*p2mup6*p3mup5*F(1858)+
     -   p2mup3*p2mup6*p3mup5*F(1858)+p1mup3*p2mup5*p3mup6*F(1858)+p2
     -   mup3*p2mup5*p3mup6*F(1858)+p2mup5*p2mup6*p3mup3*F(1860)+p2mu
     -   p6*p3mup3*p3mup5*F(1861)+p2mup5*p3mup3*p3mup6*F(1861)+p3mup3
     -   *p3mup5*p3mup6*F(1863)+p2mup5*p2mup6*p5mup3*F(1864)+p2mup5*p
     -   3mup3*p5mup6*F(1865)+p1mup5*(p2mup6*(-(p4mup3*F(1834))+(p1mu
     -   p3+p2mup3)*F(1853)+p3mup3*F(1860)+p5mup3*F(1864))-p5mup6*(p4
     -   mup3*F(1843)+(p1mup3+p2mup3)*F(1855)-p5mup3*F(1862)-p3mup3*F
     -   (1865))+p3mup6*(-(p4mup3*F(1842))+(p1mup3+p2mup3)*F(1858)+p3
     -   mup3*F(1861)+p5mup3*F(1865))))-2*(p1mup3*p2mup6*p5mup5*F(185
     -   0)+p2mup3*p2mup6*p5mup5*F(1850)-p1mup3*p3mup6*p5mup5*F(1851)
     -   -p2mup3*p3mup6*p5mup5*F(1851)-p1mup6*(p4mup3*(p5mup5*F(1829)
     -   -2*((p1mup5+p2mup5)*F(1834)+p3mup5*F(1842)))-p2mup3*p5mup5*F
     -   (1850)+p5mup3*p5mup5*F(1856)+p1mup3*(-(p5mup5*F(1850))+2*((p
     -   1mup5+p2mup5)*F(1853)+p3mup5*F(1858)))-p3mup3*p5mup5*F(1859)
     -   +2*(p1mup5*p2mup3*F(1853)+p2mup3*p2mup5*F(1853)+p2mup3*p3mup
     -   5*F(1858)+p1mup5*p3mup3*F(1860)+p2mup5*p3mup3*F(1860)+p3mup3
     -   *p3mup5*F(1861)+p1mup5*p5mup3*F(1864)+p2mup5*p5mup3*F(1864)+
     -   p3mup5*p5mup3*F(1865))))
       Fa15(341) = -(p3mup3*p5mup5*p5mup6*F(1866))+2*(p2mup6*(p5mup3*p
     -   5mup5*F(1856)-p3mup3*p5mup5*F(1859)+2*p3mup5*p5mup3*F(1865))
     -   +p3mup6*p5mup3*p5mup5*F(1867)+p3mup3*p3mup6*p5mup5*F(1868))+
     -   4*(p2mup5*p5mup3*(p5mup6*F(1862)+p3mup6*F(1865))+p3mup5*p3mu
     -   p6*p5mup3*F(1869)+p3mup3*p3mup5*p5mup6*F(1869)+p3mup5*p5mup3
     -   *p5mup6*F(1870))
       Fa15(342) = p6mup3*(-(p5mup5*p5mup6*F(1871))+(p3mup6*p4mup5+2*p
     -   3mup5*p4mup6)*F(1872)+((p1mup6+p2mup6)*p4mup5+2*(p1mup5+p2mu
     -   p5)*p4mup6)*F(1873))
       Fa15(343) = -(p2mup3*p5mup5*p6mup6*F(1876))+p4mup5*p5mup6*p6mup
     -   3*F(1880)-p4mup6*p5mup5*p6mup3*F(1881)-3*p4mup5*p4mup6*p6mup
     -   3*F(1886)-2*(p4mup3*p4mup6*p6mup5*F(1882)-p2mup3*p4mup6*p6mu
     -   p5*F(1887))+6*(p1mup5*p2mup3*p6mup6*F(1884)+p2mup3*p2mup5*p6
     -   mup6*F(1884)-p2mup3*p3mup5*p6mup6*F(1885)+p2mup3*p4mup5*p6mu
     -   p6*F(1888))+4*(p1mup5*p5mup6*p6mup3*F(1874)+p2mup5*p5mup6*p6
     -   mup3*F(1874)-p1mup6*p5mup5*p6mup3*F(1875)-p2mup6*p5mup5*p6mu
     -   p3*F(1875)-p2mup3*p5mup6*p6mup5*F(1877)+p3mup5*p5mup6*p6mup3
     -   *F(1878)-p3mup6*p5mup5*p6mup3*F(1879)+p5mup3*p5mup6*p6mup5*F
     -   (1883)+p1mup3*p1mup6*p6mup5*F(1884)+p1mup6*p2mup3*p6mup5*F(1
     -   884)+p1mup3*p2mup6*p6mup5*F(1884)+p2mup3*p2mup6*p6mup5*F(188
     -   4)-p2mup3*p3mup6*p6mup5*F(1885)-p3mup3*p3mup6*p6mup5*F(1889)
     -   )
       Fa15(344) = p6mup5*(p1mup3*(2*(p5mup6*F(1877)+p3mup6*F(1885))-p
     -   4mup6*F(1887))+p4mup6*p5mup3*F(1894)-p3mup3*p4mup6*F(1896)+2
     -   *((p1mup6+p2mup6)*(p4mup3*F(1890)-p5mup3*F(1891))+p3mup6*p4m
     -   up3*F(1892)-p3mup6*p5mup3*F(1893)+p4mup3*p5mup6*F(1895)-p3mu
     -   p3*p5mup6*F(1897)+(p1mup6+p2mup6)*p3mup3*F(1898)))
       Fa15(345) = p6mup3*p6mup5*(p4mup6*F(1899)-2*p5mup6*F(9579))
       Fa15(346) = p6mup6*(3*p2mup3*p6mup5*F(1900)+p2mup5*p5mup3*F(190
     -   1))
       Fa15(347) = p6mup6*(p1mup3*(p5mup5*F(1876)-6*((p1mup5+p2mup5)*F
     -   (1884)-p3mup5*F(1885)+p4mup5*F(1888)))-p4mup3*p5mup5*F(1902)
     -   -p5mup3*p5mup5*F(1903)+2*(p1mup5*p5mup3*F(1901)+p1mup5*p4mup
     -   3*F(1904)+p2mup5*p4mup3*F(1904)+p3mup5*p4mup3*F(1905)+p3mup5
     -   *p5mup3*F(1906)+p4mup5*p5mup3*F(1907)+p4mup3*p4mup5*F(1908))
     -   +p3mup3*(-(p5mup5*F(1909))+6*(p3mup5*F(1889)+(p1mup5+p2mup5)
     -   *F(1898)+p4mup5*F(1910))))
       Fa15(348) = p6mup5*p6mup6*(-(p4mup3*F(1912))+p5mup3*F(1913)+3*(
     -   p1mup3*F(1900)+p3mup3*F(1914)))
       Fa15(349) = p2mup2*(p5mup3*((p1mup6+p2mup6)*F(1918)-4*(p5mup6*F
     -   (1917)-p6mup6*F(1919)))+4*(p1mup6*p4mup3*F(1921)+p2mup6*p4mu
     -   p3*F(1921)+p4mup3*p5mup6*F(1922)+p4mup3*p6mup6*F(1923)+p1mup
     -   6*p3mup3*F(1925)+p2mup6*p3mup3*F(1925)+p3mup3*p5mup6*F(1927)
     -   +p3mup3*p6mup6*F(1928)))+4*(mup3mup6*p2mup2*F(1915)+p2mup3*(
     -   p2mup2*((p1mup6+p2mup6)*F(1916)+p5mup6*F(1920)-p6mup6*F(1924
     -   ))+p5mup2*s12*(-4*p6mup6*F(1900)+(p1mup6+p2mup6)*F(1933)+p5m
     -   up6*F(1936))))+p5mup2*s12*(mup3mup6*F(1926)-4*(p5mup3*(-2*p6
     -   mup6*F(1929)+(p1mup6+p2mup6)*F(1930)+p5mup6*F(1931))+p4mup3*
     -   (p5mup6*F(1932)+(p1mup6+p2mup6)*F(1937)-4*p6mup6*F(1938))+p3
     -   mup3*(-(p5mup6*F(1934))-4*p6mup6*F(1935)+(p1mup6+p2mup6)*F(1
     -   939))))
       Fa15(350) = p1mup2*(p2mup6*p5mup3*F(1918)+p1mup6*(p5mup3*F(1918
     -   )+4*(p2mup3*F(1916)+p4mup3*F(1921)+p3mup3*F(1925)))+4*(mup3m
     -   up6*F(1915)-p5mup3*p5mup6*F(1917)+p5mup3*p6mup6*F(1919)+p2mu
     -   p6*p4mup3*F(1921)+p4mup3*p5mup6*F(1922)+p4mup3*p6mup6*F(1923
     -   )+p2mup3*(p2mup6*F(1916)+p5mup6*F(1920)-p6mup6*F(1924))+p2mu
     -   p6*p3mup3*F(1925)+p3mup3*p5mup6*F(1927)+p3mup3*p6mup6*F(1928
     -   )))+4*p1mup3*((p1mup2+p2mup2)*((p1mup6+p2mup6)*F(1916)+p5mup
     -   6*F(1920)-p6mup6*F(1924))+p5mup2*s12*(-4*p6mup6*F(1900)+(p1m
     -   up6+p2mup6)*F(1933)+p5mup6*F(1936)))
       Fa15(351) = (p1mup2+p2mup2)*(p6mup3*p6mup6*F(1942)+4*(p5mup6*p6
     -   mup3*F(1941)-p4mup6*p5mup3*F(1943)+p4mup6*p6mup3*F(1944)))
       Fa15(352) = p4mup3*p6mup2*(p5mup6*F(1945)+p1mup6*F(1946)+p6mup6
     -   *F(1947))
       Fa15(353) = p2mup2*p2mup3*p3mup6*F(1948)-p2mup2*p3mup6*p5mup3*F
     -   (1949)+p2mup2*p3mup6*p6mup3*F(1950)+p2mup2*p2mup3*p4mup6*F(1
     -   951)+p1mup3*p2mup2*(p3mup6*F(1948)+p4mup6*F(1951))-p2mup2*p3
     -   mup6*p4mup3*F(1952)-p2mup2*p4mup3*p4mup6*F(1953)+p1mup2*(p1m
     -   up3*p3mup6*F(1948)+p2mup3*p3mup6*F(1948)-p3mup6*p5mup3*F(194
     -   9)+p3mup6*p6mup3*F(1950)+p1mup3*p4mup6*F(1951)+p2mup3*p4mup6
     -   *F(1951)-p3mup6*p4mup3*F(1952)-p4mup3*p4mup6*F(1953))+p3mup6
     -   *p4mup3*p5mup2*s12*F(1954)-p4mup3*p4mup6*p5mup2*s12*F(1955)
       Fa15(354) = p1mup6*p2mup2*p6mup3*F(1956)+s12*(p2mup3*p5mup6*p6m
     -   up2*F(1958)-p5mup3*(p5mup6*p6mup2*F(1957)+p3mup6*p5mup2*F(19
     -   60)-p3mup6*p6mup2*F(1961)+p4mup6*p5mup2*F(1962)-p4mup6*p6mup
     -   2*F(1963))+p3mup3*p5mup6*p6mup2*F(1964)-p1mup3*p3mup6*p5mup2
     -   *F(1965)-p2mup3*p3mup6*p5mup2*F(1965)+p3mup6*p4mup3*p6mup2*F
     -   (1966)+p4mup3*p4mup6*p6mup2*F(1967)+p1mup3*p4mup6*p5mup2*F(1
     -   968)+p2mup3*p4mup6*p5mup2*F(1968)-p3mup3*p3mup6*p5mup2*F(196
     -   9))+p2mup2*p3mup3*p3mup6*F(1970)+p2mup2*p3mup3*p4mup6*F(1971
     -   )+p1mup2*(p1mup6*p6mup3*F(1956)+p3mup3*p3mup6*F(1970)+p3mup3
     -   *p4mup6*F(1971))+p3mup3*p4mup6*p5mup2*s12*F(1972)+p2mup3*p4m
     -   up6*p6mup2*s12*F(9598)+p3mup3*p3mup6*p6mup2*s12*F(9599)+p2mu
     -   p3*p3mup6*p6mup2*s12*F(9604)+p2mup6*((p1mup2+p2mup2)*p6mup3*
     -   F(1956)+p6mup2*s12*(p4mup3*F(1946)+p5mup3*F(1959)+p2mup3*F(9
     -   603)+p3mup3*F(9604)))+p3mup3*p4mup6*p6mup2*s12*F(9606)
       Fa15(355) = p6mup3*(p5mup6*p6mup2*F(1973)+(p1mup6+p2mup6)*p5mup
     -   2*F(1974)+p3mup6*p5mup2*F(1975)+p5mup2*p5mup6*F(1976)+p4mup6
     -   *p5mup2*F(1977)+p3mup6*p6mup2*F(9581)+p2mup6*p6mup2*F(9584)+
     -   p4mup6*p6mup2*F(9587))-p1mup3*p6mup2*(p5mup6*F(1958)+p4mup6*
     -   F(9598)+p2mup6*F(9603)+p3mup6*F(9604))
       Fa15(356) = p1mup6*p6mup2*(p5mup3*F(1959)-p6mup3*F(9584)+p1mup3
     -   *F(9603)+p3mup3*F(9604))
       Fa15(357) = p6mup2*p6mup6*((p1mup3+p2mup3)*F(1900)-p6mup3*F(191
     -   1)-p3mup3*F(1935)-p5mup3*F(1978))
       Fa15(358) = p5mup3*F(1981)+2*(p3mup3*F(1979)+p4mup3*F(1980)+p6m
     -   up3*F(1982)+(p1mup3+p2mup3)*F(1983))
       Fa15(359) = mup5mup6*((p2mup2*F(1984))/s12+p5mup2*F(1985))
       Fa15(360) = p1mup5*p1mup6*p2mup2*F(1986)+p1mup5*p2mup2*p2mup6*F
     -   (1986)+p2mup2*p2mup5*p2mup6*F(1986)+p2mup2*p2mup6*p5mup5*F(1
     -   987)+p2mup2*p2mup6*p6mup5*F(1988)+p1mup5*p2mup2*p5mup6*F(198
     -   9)+8*p5mup2*p5mup5*p5mup6*s12*F(1990)+16*(p5mup2*p5mup6*p6mu
     -   p5*s12*F(1991)-(p1mup5+p2mup5)*p5mup2*p5mup6*s12*F(1992))
       Return
       End
