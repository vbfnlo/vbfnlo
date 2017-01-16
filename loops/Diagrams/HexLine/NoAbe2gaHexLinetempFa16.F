c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa16(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa16)
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
       Complex*16 Fa16(361:384)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa16(361) = p1mup6*p2mup2*p2mup5*F(1986)+p1mup6*p2mup2*p5mup5*F
     -   (1987)+p1mup6*p2mup2*p6mup5*F(1988)+p2mup2*p2mup5*p5mup6*F(1
     -   989)+p1mup5*p2mup2*p6mup6*F(1993)+p2mup2*p2mup5*p6mup6*F(199
     -   3)-4*p2mup2*p5mup5*p5mup6*F(1994)+p2mup2*p5mup5*p6mup6*F(199
     -   5)+p2mup2*p5mup6*p6mup5*F(1996)+p2mup2*p6mup5*p6mup6*F(1997)
     -   +mup5mup6*p6mup2*s12*F(1998)-8*(p1mup5*p5mup2*p6mup6*s12*F(1
     -   874)+p2mup5*p5mup2*p6mup6*s12*F(1874)+p5mup2*p6mup5*p6mup6*s
     -   12*F(1999)-p1mup6*p5mup2*p5mup5*s12*F(2000)-p2mup6*p5mup2*p5
     -   mup5*s12*F(2000))+p5mup2*p5mup5*p6mup6*s12*F(2002)+16*(p1mup
     -   6*p5mup2*p6mup5*s12*F(2001)+p2mup6*p5mup2*p6mup5*s12*F(2001)
     -   -(p1mup5+p2mup5)*(p1mup6+p2mup6)*p5mup2*s12*F(2003))
       Fa16(362) = -4*p2mup2*p4mup5*(p5mup6*F(2004)+p4mup6*F(2005)-p6m
     -   up6*F(2006))+p1mup2*(p1mup6*p2mup5*F(1986)+p2mup5*p2mup6*F(1
     -   986)+p1mup6*p5mup5*F(1987)+p2mup6*p5mup5*F(1987)+p1mup6*p6mu
     -   p5*F(1988)+p2mup6*p6mup5*F(1988)+p2mup5*p5mup6*F(1989)+p2mup
     -   5*p6mup6*F(1993)+p1mup5*((p1mup6+p2mup6)*F(1986)+p5mup6*F(19
     -   89)+p6mup6*F(1993))+p5mup5*p6mup6*F(1995)+p5mup6*p6mup5*F(19
     -   96)+p6mup5*p6mup6*F(1997)-4*(p5mup5*p5mup6*F(1994)+p4mup5*p5
     -   mup6*F(2004)+p4mup5*p4mup6*F(2005)-p4mup5*p6mup6*F(2006)))
       Fa16(363) = 4*((p1mup2+p2mup2)*(p1mup6+p2mup6)*p4mup5*F(2007)+p
     -   4mup5*p5mup2*s12*(2*p6mup6*F(2011)-(p1mup6+p2mup6)*F(2013)-p
     -   5mup6*F(2014)))+8*p2mup5*p6mup2*s12*(p6mup6*F(1875)+2*(p5mup
     -   6*F(2012)-p1mup6*F(2015)))+p3mup5*((p1mup2+p2mup2)*((p1mup6+
     -   p2mup6)*F(2008)+p6mup6*F(2009)-4*p5mup6*F(2010))-8*p5mup2*s1
     -   2*(p6mup6*F(1878)+2*((p1mup6+p2mup6)*F(2016)+p5mup6*F(2017))
     -   ))
       Fa16(364) = (p1mup2+p2mup2)*p4mup6*(p5mup5*F(2018)+(p1mup5+p2mu
     -   p5)*F(2019)-8*p3mup5*F(2020))
       Fa16(365) = p6mup2*(p4mup5*p5mup6*F(2021)+p4mup5*p6mup6*F(2022)
     -   +p5mup5*p5mup6*F(2023)+2*(p2mup6*p5mup5*F(2024)-p3mup5*(p6mu
     -   p6*F(1879)+2*p5mup6*F(2025)))+4*(p2mup5*p2mup6*F(2015)+(p1mu
     -   p6+p2mup6)*p3mup5*F(2026))+(p1mup6+p2mup6)*p4mup5*F(2027))
       Fa16(366) = p1mup2*p3mup6*p5mup5*F(2029)+p2mup2*p3mup6*p5mup5*F
     -   (2029)+p1mup2*p1mup5*p3mup6*F(2030)+p1mup5*p2mup2*p3mup6*F(2
     -   030)+(p1mup2+p2mup2)*p2mup5*p3mup6*F(2030)+p1mup2*p3mup6*p6m
     -   up5*F(2031)+p2mup2*p3mup6*p6mup5*F(2031)-16*(p1mup5*p3mup6*p
     -   5mup2*s12*F(2016)+p2mup5*p3mup6*p5mup2*s12*F(2016)-p3mup6*p5
     -   mup2*p6mup5*s12*F(2035))-4*(p1mup2*p3mup6*p4mup5*F(2032)+p2m
     -   up2*p3mup6*p4mup5*F(2032)-p5mup5*p6mup2*p6mup6*s12*F(2033)-p
     -   4mup6*p5mup2*p5mup5*s12*F(2036))+8*(p1mup5*p6mup2*s12*(p6mup
     -   6*F(1875)+2*(p5mup6*F(2012)-(p1mup6+p2mup6)*F(2015)))-p1mup6
     -   *p5mup5*p6mup2*s12*F(2024)+(p1mup2+p2mup2)*p4mup6*p6mup5*F(2
     -   028)+p3mup6*p5mup2*p5mup5*s12*F(2034)+p4mup6*p5mup2*p6mup5*s
     -   12*F(2037)+(p1mup5+p2mup5)*p4mup6*p5mup2*s12*F(2038))
       Fa16(367) = p6mup2*(4*(p1mup6*p6mup5*F(1875)+p2mup6*p6mup5*F(18
     -   75))+p4mup6*p5mup5*F(2039)-2*(-(p4mup6*p6mup5*F(2040))+(p1mu
     -   p5+p2mup5)*p4mup6*F(2041)+p3mup6*(-2*(p6mup5*F(1879)-(p1mup5
     -   +p2mup5)*F(2026))+p5mup5*F(2042))-2*p5mup6*p6mup5*F(2043)))-
     -   4*p3mup5*p3mup6*p5mup2*F(2044)
       Fa16(368) = p1mup2*p3mup5*p3mup6*F(2045)+p2mup2*p3mup5*p3mup6*F
     -   (2045)-4*s12*(p4mup6*(p4mup5*p5mup2*F(2048)+p4mup5*p6mup2*F(
     -   2049)+2*(p3mup5*p6mup2*F(2047)-p3mup5*p5mup2*F(2051)))+p3mup
     -   6*(p4mup5*p5mup2*F(2046)+4*p3mup5*p6mup2*F(2050)+p4mup5*p6mu
     -   p2*F(2052)))
       Fa16(369) = 2*(p1mup2+p2mup2)*F(2053)+s12*(p5mup2*F(2054)+2*p6m
     -   up2*F(2055))
       Fa16(370) = mup5mup6*p5mup2*(p3mup3*F(2056)+p5mup3*F(2057))
       Fa16(371) = mup5mup6*p5mup2*((p1mup3+p2mup3)*F(2058)+p4mup3*F(2
     -   059))
       Fa16(372) = mup5mup6*p6mup2*(p5mup3*F(2060)+p3mup3*F(2061)+p4mu
     -   p3*F(2062)-(p1mup3+p2mup3)*F(2063))
       Fa16(373) = p5mup2*(mup3mup6*(p6mup5*F(2066)-4*(p5mup5*F(2065)-
     -   (p1mup5+p2mup5)*F(2068)))+mup3mup5*((p1mup6+p2mup6)*F(2064)+
     -   p6mup6*F(2067)+p5mup6*F(2069)))
       Fa16(374) = -4*mup5mup6*p5mup2*p6mup3*F(2070)+mup3mup6*p6mup2*(
     -   (p1mup5+p2mup5)*F(2071)+p5mup5*F(2072)+p6mup5*F(2073))-8*mup
     -   3mup5*p5mup6*p6mup2*F(2074)
       Fa16(375) = 4*(mup3mup5*p6mup2*(2*(p1mup6+p2mup6)*F(2075)+p6mup
     -   6*F(2076))-mup3mup6*p5mup2*(p4mup5*F(2077)+p3mup5*F(2078)))
       Fa16(376) = mup3mup6*p6mup2*(p3mup5*F(2079)+4*p4mup5*F(2080))
       Fa16(377) = mup3mup5*(p3mup6*p5mup2*F(2081)+p4mup6*p5mup2*F(208
     -   2)-8*(p3mup6*p6mup2*F(2083)+p4mup6*p6mup2*F(2084)))
       Fa16(378) = p5mup2*(p2mup6*p5mup3*p5mup5*F(2104)+p4mup3*(-2*(p4
     -   mup6*p5mup5*F(2089)-p3mup6*p5mup5*F(2093)+p5mup5*p5mup6*F(20
     -   99))+p2mup6*p5mup5*F(2102)-4*(p4mup5*(p5mup6*F(2088)-p3mup6*
     -   F(2092)+p4mup6*F(2097)-p2mup6*F(2101))+p3mup5*p5mup6*F(2106)
     -   +p2mup5*p4mup6*F(2109)-p2mup5*p5mup6*F(2111)-p3mup5*p4mup6*F
     -   (2114)))-2*(p4mup6*p5mup3*p5mup5*F(2091)-p3mup6*p5mup3*p5mup
     -   5*F(2095)+p2mup3*(p4mup6*p5mup5*F(2085)+p5mup5*p5mup6*F(2086
     -   )-2*(p4mup5*p5mup6*F(2087)+p4mup5*p4mup6*F(2096)))+p5mup3*p5
     -   mup5*p5mup6*F(2100)+p3mup3*p5mup5*p5mup6*F(2107)-p3mup3*p4mu
     -   p6*p5mup5*F(2115))-4*(p4mup5*p5mup3*p5mup6*F(2090)-p3mup6*p4
     -   mup5*p5mup3*F(2094)+p4mup5*p4mup6*p5mup3*F(2098)-p2mup6*p4mu
     -   p5*p5mup3*F(2103)+p3mup3*p4mup5*p5mup6*F(2105)+p3mup5*p5mup3
     -   *p5mup6*F(2108)-p3mup3*p4mup5*p4mup6*F(2110)-p2mup5*p5mup3*p
     -   5mup6*F(2112)+p2mup5*p4mup6*p5mup3*F(2113)-p3mup5*p4mup6*p5m
     -   up3*F(2116)))
       Fa16(379) = p5mup2*(p3mup3*p3mup6*(p5mup5*F(2119)-2*p4mup5*F(21
     -   21))+2*(p1mup3*p1mup5*p1mup6*F(2117)+p3mup5*(p3mup6*p4mup3*F
     -   (2118)+p3mup6*p5mup3*F(2120)-p3mup3*p3mup6*F(2122)+p3mup3*p5
     -   mup6*F(2123)-p3mup3*p4mup6*F(2124))))
       Fa16(380) = p5mup2*(p1mup3*p1mup6*(-(p5mup5*F(2125))+4*(p2mup5*
     -   F(2117)+p3mup5*F(2127)+p4mup5*F(2128)))+4*p1mup5*(p1mup3*(p2
     -   mup6*F(2117)+p3mup6*F(2126)+p4mup6*F(2129)+p5mup6*F(2130))+p
     -   1mup6*(p2mup3*F(2117)+p4mup3*F(2131)+p5mup3*F(2132)-p3mup3*F
     -   (2133))))
       Fa16(381) = p5mup2*(p1mup6*p2mup3*p2mup5*F(2117)+p1mup5*p2mup3*
     -   p2mup6*F(2117)+p1mup3*p2mup5*p2mup6*F(2117)+p1mup5*p2mup3*p3
     -   mup6*F(2126)+p1mup3*p2mup5*p3mup6*F(2126)+p1mup6*p2mup3*p3mu
     -   p5*F(2127)+p1mup3*p2mup6*p3mup5*F(2127)-(p1mup6*p2mup5+p1mup
     -   5*p2mup6)*p3mup3*F(2133))
       Fa16(382) = p5mup2*(p1mup3*(-(p2mup6*p5mup5*F(2125))+4*(p2mup6*
     -   p4mup5*F(2128)+p2mup5*p4mup6*F(2129)+p2mup5*p5mup6*F(2130)+p
     -   3mup5*p3mup6*F(2134)+p3mup5*p4mup6*F(2135)+p3mup6*p4mup5*F(2
     -   136)))+4*p1mup5*(p2mup3*p4mup6*F(2129)+p2mup3*p5mup6*F(2130)
     -   +p2mup6*p4mup3*F(2131)+p2mup6*p5mup3*F(2132)-p3mup3*p4mup6*F
     -   (2137)-p3mup6*p4mup3*F(2138)+p3mup3*p3mup6*F(2141))+p1mup6*(
     -   p2mup3*(-(p5mup5*F(2125))+4*p4mup5*F(2128))+4*(p2mup5*p4mup3
     -   *F(2131)+p2mup5*p5mup3*F(2132)+p3mup3*p4mup5*F(2139)+p3mup5*
     -   p4mup3*F(2140)-p3mup3*p3mup5*F(2142))))
       Fa16(383) = p5mup2*(-2*p1mup3*(p4mup6*p5mup5*F(2085)+p5mup5*p5m
     -   up6*F(2086)+p3mup6*p5mup5*F(2143)-2*(p4mup5*p5mup6*F(2087)+p
     -   4mup5*p4mup6*F(2096)+p3mup5*p5mup6*F(2145)))+4*(p2mup3*(p2mu
     -   p5*p2mup6*F(2117)+p2mup5*p3mup6*F(2126)+p2mup6*p3mup5*F(2127
     -   )+p2mup5*p4mup6*F(2129))-p2mup5*p2mup6*p3mup3*F(2133)+p1mup5
     -   *(-(p4mup3*p4mup6*F(2109))+p4mup3*p5mup6*F(2111)+p5mup3*p5mu
     -   p6*F(2112)-p4mup6*p5mup3*F(2113)+p3mup3*p5mup6*F(2146)+p3mup
     -   6*p5mup3*F(2147)))+p1mup6*(p4mup3*p5mup5*F(2102)+p5mup3*p5mu
     -   p5*F(2104)-p3mup3*p5mup5*F(2144)+4*(p4mup3*p4mup5*F(2101)+p4
     -   mup5*p5mup3*F(2103)+p3mup5*p5mup3*F(2148))))
       Fa16(384) = -8*(p2mup6*p3mup3*p5mup2*p5mup5*F(2144)-p5mup3*p5mu
     -   p5*p5mup6*p6mup2*F(2149)+p2mup3*(p2mup6*p5mup2*(p5mup5*F(212
     -   5)-4*p4mup5*F(2128))-2*(2*(p2mup5*p5mup2*p5mup6*F(2130)+p3mu
     -   p5*(p3mup6*p5mup2*F(2134)+p4mup6*p5mup2*F(2135)+p5mup2*p5mup
     -   6*F(2145)+p1mup6*p6mup2*F(2150)))+p3mup6*(-(p5mup2*p5mup5*F(
     -   2143))+2*(p4mup5*p5mup2*F(2136)+p1mup5*p6mup2*F(2152))))))+p
     -   3mup3*p3mup6*p5mup5*p6mup2*F(2155)+32*(p2mup6*p3mup3*p4mup5*
     -   p5mup2*F(2139)+p2mup6*p3mup5*p4mup3*p5mup2*F(2140)-p2mup6*p3
     -   mup3*p3mup5*p5mup2*F(2142)+p2mup6*p3mup5*p5mup2*p5mup3*F(214
     -   8)+p1mup3*p1mup6*p3mup5*p6mup2*F(2150)+p1mup3*p2mup6*p3mup5*
     -   p6mup2*F(2150)-p4mup3*p4mup5*p4mup6*p6mup2*F(2151)+p1mup3*p1
     -   mup5*p3mup6*p6mup2*F(2152)-p3mup3*p3mup5*p5mup6*p6mup2*F(215
     -   3)+p3mup5*p3mup6*p5mup3*p6mup2*F(2154)-p1mup5*(p1mup6+p2mup6
     -   )*p3mup3*p6mup2*F(2156)+p2mup5*(p2mup6*p5mup2*(p4mup3*F(2131
     -   )+p5mup3*F(2132))+p3mup6*(-(p4mup3*p5mup2*F(2138))+p5mup2*p5
     -   mup3*F(2147)+p1mup3*p6mup2*F(2152))+p3mup3*(-(p4mup6*p5mup2*
     -   F(2137))+p3mup6*p5mup2*F(2141)+p5mup2*p5mup6*F(2146)-p1mup6*
     -   p6mup2*F(2156))))
       Return
       End
