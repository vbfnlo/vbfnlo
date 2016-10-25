c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenGSim6(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa6)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenEvenGSim/F
      Complex*16 Fa(123),F(3806)
      Real*8 P(1320) 
        COMMON/PPenEvenGSim/P
       Complex*16 Fa6(41:48)
       COMMON/FaPenEvenGSim/Fa
       Fa6(41) = mup1mup3*mup2mup5*p3mup4*F(1836)+mup1mup3*mup2mup4*p3
     -   mup5*F(1836)+mup1mup4*mup2mup5*p3mup3*F(1837)+mup1mup4*mup2m
     -   up3*p3mup5*F(1837)+mup2mup3*mup4mup5*p1mup1*F(1877)+mup2mup4
     -   *mup3mup5*p1mup1*F(1878)+mup2mup5*mup3mup4*p1mup1*F(1879)+mu
     -   p1mup3*mup4mup5*p1mup2*F(1880)+mup1mup4*mup3mup5*p1mup2*F(18
     -   81)+mup1mup4*mup2mup5*p1mup3*F(1884)-mup1mup3*mup2mup5*p1mup
     -   4*F(1884)-mup1mup4*mup2mup3*p1mup5*F(1887)+mup1mup3*mup2mup4
     -   *p1mup5*F(1887)-mup2mup4*mup3mup5*p2mup1*F(1888)+mup2mup3*mu
     -   p4mup5*p2mup1*F(1888)+mup2mup5*mup3mup4*p2mup1*F(1889)+mup1m
     -   up4*mup3mup5*p2mup2*F(1890)-mup1mup3*mup4mup5*p2mup2*F(1890)
     -   +mup1mup4*mup2mup5*p2mup3*F(1893)+mup1mup3*mup2mup5*p2mup4*F
     -   (1896)-mup1mup4*mup2mup3*p2mup5*F(1899)+mup1mup3*mup2mup4*p2
     -   mup5*F(1899)+mup2mup3*mup4mup5*p3mup1*F(1900)+mup2mup4*mup3m
     -   up5*p3mup1*F(1901)+mup2mup5*mup3mup4*p3mup1*F(1902)+mup1mup3
     -   *mup4mup5*p3mup2*F(1903)+mup1mup4*mup3mup5*p3mup2*F(1904)+mu
     -   p2mup3*mup4mup5*p4mup1*F(1906)+mup2mup4*mup3mup5*p4mup1*F(19
     -   07)+mup2mup5*mup3mup4*p4mup1*F(1908)+mup1mup3*mup4mup5*p4mup
     -   2*F(1909)+mup1mup4*mup3mup5*p4mup2*F(1910)+mup1mup4*mup2mup5
     -   *p4mup3*F(1913)+mup1mup3*mup2mup5*p4mup4*F(1914)+mup1mup5*(m
     -   up2mup3*(p3mup4*F(1837)+p1mup4*F(1885)+p2mup4*F(1897))+mup3m
     -   up4*(p1mup2*F(1882)+p2mup2*F(1891)+p3mup2*F(1905)+p4mup2*F(1
     -   911))+mup2mup4*(p3mup3*F(1836)-p1mup3*F(1885)+p2mup3*F(1894)
     -   +p4mup3*F(1914)))+mup1mup2*(mup3mup4*(p3mup5*F(1835)+p1mup5*
     -   F(1886)+p2mup5*F(1898))+mup4mup5*(p3mup3*F(1838)-p1mup3*F(18
     -   83)+p2mup3*F(1892)+p4mup3*F(1912))+mup3mup5*(p3mup4*F(1835)+
     -   p1mup4*F(1883)+p2mup4*F(1895)+p4mup4*F(1915)))+8*(-(mup4mup5
     -   *p3mup3*p4mup1*p4mup2*F(1848))+mup3mup5*p3mup4*p4mup1*p4mup2
     -   *F(1848)+mup3mup4*p3mup5*p4mup1*p4mup2*F(1848)-(p3mup2*p4mup
     -   1+p3mup1*p4mup2)*(mup4mup5*p4mup3-mup3mup5*p4mup4-mup3mup4*p
     -   4mup5)*F(1916))
       Fa6(42) = (mup1mup5*mup2mup3*p4mup4*F(1917)+p4mup5*(mup1mup2*mu
     -   p3mup4*F(1918)+mup1mup3*mup2mup4*F(1919)+mup1mup4*mup2mup3*F
     -   (1920))+16*(p1mup2*(p1mup4*p1mup5*(p2mup1*p3mup3*F(1935)-p2m
     -   up3*p3mup1*F(1936))+p1mup3*(p1mup4*p2mup1*p3mup5*F(1932)+p1m
     -   up4*p2mup5*p3mup1*F(1933)+2*(p1mup5*p2mup1*p3mup4*F(1934)+p1
     -   mup5*p2mup4*p3mup1*F(1937))))+p1mup3*p1mup4*p1mup5*(-(p2mup1
     -   *p3mup2*F(1938))+p2mup2*p3mup1*F(1939))+p1mup1*(p1mup2*(p1mu
     -   p4*p2mup3*p3mup5*F(1923)+p1mup4*p2mup5*p3mup3*F(1924)+4*p1mu
     -   p5*p2mup3*p3mup4*F(1925)+2*(p1mup3*(p2mup4*p3mup5*F(1921)+p2
     -   mup5*p3mup4*F(1922))+p1mup5*p2mup4*p3mup3*F(1926)))+p1mup4*p
     -   1mup5*(p2mup2*p3mup3*F(1930)+p2mup3*p3mup2*F(1931))+p1mup3*(
     -   p1mup4*p2mup2*p3mup5*F(1927)-p1mup4*p2mup5*p3mup2*F(1928)+2*
     -   p1mup5*p2mup2*p3mup4*F(1929)+4*p1mup5*p2mup4*p3mup2*F(1940))
     -   )))/2.
       Fa6(43) = mup1mup4*mup2mup3*p4mup5*F(1843)+(mup1mup5*mup2mup3*p
     -   4mup4+mup1mup3*mup2mup4*p4mup5)*F(1844)+mup1mup2*mup3mup4*p4
     -   mup5*F(1874)
       Fa6(44) = mup1mup2*mup3mup4*p4mup5*F(1842)+(mup1mup5*mup2mup3*p
     -   4mup4+mup1mup3*mup2mup4*p4mup5)*F(1843)+mup1mup4*mup2mup3*p4
     -   mup5*F(1844)
       Fa6(45) = mup1mup5*mup2mup3*p4mup4*F(1913)+p4mup5*(mup1mup4*mup
     -   2mup3*F(1913)+mup1mup3*mup2mup4*F(1914)+mup1mup2*mup3mup4*F(
     -   1915))
       Fa6(46) = mup2mup3*p2mup4*p2mup5*F(1943)+mup2mup3*p3mup4*p3mup5
     -   *F(1946)+mup2mup3*p2mup4*p3mup5*F(1951)+mup2mup3*p2mup5*p3mu
     -   p4*F(1952)+mup4mup5*p2mup2*(p4mup3*F(1953)-2*p2mup3*F(1955)+
     -   p3mup3*F(1957))+mup3mup4*p2mup2*p3mup5*F(1958)+mup4mup5*p2mu
     -   p3*p3mup2*F(1959)+mup2mup4*p2mup3*p3mup5*F(1960)+mup3mup4*p2
     -   mup5*p3mup2*F(1961)+mup2mup4*p2mup5*p3mup3*F(1962)+mup3mup5*
     -   (p2mup4*p3mup2*F(1949)+p2mup2*(2*p2mup4*F(1941)+p3mup4*F(194
     -   7)+p4mup4*F(1954))+2*(p3mup2*p3mup4*F(1944)+p4mup2*p4mup4*F(
     -   1964)))+mup2mup5*(p2mup3*p3mup4*F(1948)+p2mup4*p3mup3*F(1950
     -   )-2*(p2mup3*p2mup4*F(1942)+p3mup3*p3mup4*F(1945)+p4mup3*p4mu
     -   p4*F(1966)))+2*(mup2mup4*p2mup3*p2mup5*F(1956)+mup3mup4*p2mu
     -   p2*p2mup5*F(1963)+mup3mup4*p4mup2*p4mup5*F(1965)+mup2mup4*p4
     -   mup3*p4mup5*F(1967)+mup2mup3*p4mup4*p4mup5*F(1968)-mup4mup5*
     -   p3mup2*p3mup3*F(1969)+mup3mup4*p3mup2*p3mup5*F(1970)+mup4mup
     -   5*p4mup2*p4mup3*F(1971)+mup2mup4*p3mup3*p3mup5*F(1981))+4*(p
     -   2mup2*p3mup3*p3mup4*p4mup5*F(1973)+p2mup2*p3mup3*p3mup5*p4mu
     -   p4*F(1974)+p2mup3*p3mup2*p3mup4*p4mup5*F(1975)+p2mup3*p3mup2
     -   *p3mup5*p4mup4*F(1976)-p2mup4*p3mup2*p3mup5*p4mup3*F(1977)+p
     -   2mup4*p3mup3*p3mup5*p4mup2*F(1978)+p2mup5*p3mup2*p3mup4*p4mu
     -   p3*F(1979)-p2mup5*p3mup3*p3mup4*p4mup2*F(1980)+p2mup2*p3mup4
     -   *p3mup5*p4mup3*F(1982)+p2mup3*p3mup4*p3mup5*p4mup2*F(1983)+p
     -   2mup2*p3mup3*p4mup4*p4mup5*F(1984)+p2mup2*p3mup4*p4mup3*p4mu
     -   p5*F(1985)+p2mup3*p3mup2*p4mup4*p4mup5*F(1986)+p2mup3*p3mup4
     -   *p4mup2*p4mup5*F(1987)-p2mup4*p3mup2*p4mup3*p4mup5*F(1988)+p
     -   2mup4*p3mup3*p4mup2*p4mup5*F(1989)+p2mup4*p3mup2*p3mup3*p4mu
     -   p5*F(1990)+p2mup5*p3mup2*p3mup3*p4mup4*F(1991)+p2mup2*p3mup5
     -   *p4mup3*p4mup4*F(1995)+p2mup3*p3mup5*p4mup2*p4mup4*F(1996)+p
     -   2mup4*p3mup5*p4mup2*p4mup3*F(1997)-p2mup5*p3mup2*p4mup3*p4mu
     -   p4*F(1998)+p2mup5*p3mup3*p4mup2*p4mup4*F(1999)+p2mup5*p3mup4
     -   *p4mup2*p4mup3*F(2000))+8*(p3mup2*p3mup3*p3mup4*p4mup5*F(199
     -   2)+p3mup2*p3mup3*p3mup5*p4mup4*F(1993)+p2mup2*p4mup3*p4mup4*
     -   p4mup5*F(2001)+p2mup4*p4mup2*p4mup3*p4mup5*F(2003)+p3mup3*p3
     -   mup4*p4mup2*p4mup5*F(2004)+p3mup4*p3mup5*p4mup2*p4mup3*F(200
     -   5)+p3mup3*p4mup2*p4mup4*p4mup5*F(2006)+p3mup4*p4mup2*p4mup3*
     -   p4mup5*F(2007)+p3mup2*p3mup4*p3mup5*p4mup3*F(2008)+p3mup2*p3
     -   mup3*p4mup4*p4mup5*F(2010)+p3mup3*p3mup5*p4mup2*p4mup4*F(201
     -   1)+p3mup2*p3mup4*p4mup3*p4mup5*F(2011)+p2mup5*p4mup2*p4mup3*
     -   p4mup4*F(2012)+p3mup5*p4mup2*p4mup3*p4mup4*F(2013)+p3mup2*p3
     -   mup5*p4mup3*p4mup4*F(2014))+16*(p2mup2*p2mup3*p2mup4*p2mup5*
     -   F(1972)+p3mup2*p3mup3*p3mup4*p3mup5*F(1994)+p2mup3*p4mup2*p4
     -   mup4*p4mup5*F(2002)+p3mup3*p3mup4*p3mup5*p4mup2*F(2009)+p3mu
     -   p2*p4mup3*p4mup4*p4mup5*F(2015))
       Fa6(47) = mup3mup5*p2mup2*p3mup4*F(2019)+mup3mup4*p2mup2*p3mup5
     -   *F(2020)+mup2mup5*(p2mup4*p3mup3+p2mup3*p3mup4)*F(2022)+mup2
     -   mup3*p2mup4*p3mup5*F(2022)+mup2mup4*p2mup5*p3mup3*F(2023)+mu
     -   p2mup3*p2mup5*p3mup4*F(2023)+mup2mup4*p2mup3*p3mup5*F(2023)+
     -   mup3mup5*p2mup4*p3mup2*F(2024)+mup3mup4*p2mup5*p3mup2*F(2025
     -   )+mup3mup5*p2mup2*p4mup4*F(2030)+mup4mup5*(p2mup3*p3mup2*F(2
     -   021)+p2mup2*(2*p2mup3*F(2016)+p3mup3*F(2019)+p4mup3*F(2030))
     -   +2*(p3mup2*p3mup3*F(2028)+p4mup2*p4mup3*F(2031)))+2*(mup2mup
     -   5*p2mup3*p2mup4*F(2017)+mup2mup4*p2mup3*p2mup5*F(2018)+mup3m
     -   up5*p2mup2*p2mup4*F(2026)+mup3mup4*p2mup2*p2mup5*F(2027)+mup
     -   3mup5*p3mup2*p3mup4*F(2028)+mup3mup4*p3mup2*p3mup5*F(2029)+m
     -   up3mup5*p4mup2*p4mup4*F(2031)+mup3mup4*p4mup2*p4mup5*F(2032)
     -   )
       Fa6(48) = mup3mup5*p2mup4*p3mup2*F(2034)+mup3mup5*p2mup2*p3mup4
     -   *F(2034)+mup2mup4*p2mup5*p3mup3*F(2035)+mup2mup5*p2mup3*p3mu
     -   p4*F(2035)+mup2mup3*p2mup4*p3mup5*F(2035)+mup2mup5*p2mup4*p3
     -   mup3*F(2036)+mup2mup3*p2mup5*p3mup4*F(2036)+mup2mup4*p2mup3*
     -   p3mup5*F(2036)+mup3mup5*p2mup2*p4mup4*F(2037)+mup4mup5*(p2mu
     -   p3*p3mup2*F(2033)+p2mup2*(2*p2mup3*F(2026)+p4mup3*F(2030)+p3
     -   mup3*F(2033))+2*(p4mup2*p4mup3*F(2031)+p3mup2*p3mup3*F(2038)
     -   ))+mup3mup4*(p2mup5*p3mup2*F(2033)+p2mup2*p3mup5*F(2033)+2*(
     -   p2mup2*p2mup5*F(2026)+p4mup2*p4mup5*F(2031)+p3mup2*p3mup5*F(
     -   2038)))+2*(mup3mup5*p2mup2*p2mup4*F(2027)+mup3mup5*p4mup2*p4
     -   mup4*F(2032)+mup3mup5*p3mup2*p3mup4*F(2039))
       Return
       End
