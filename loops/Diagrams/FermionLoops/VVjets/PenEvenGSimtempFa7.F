c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenGSim7(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa7)
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
       Complex*16 Fa7(49:56)
       COMMON/FaPenEvenGSim/Fa
       Fa7(49) = mup3mup4*p2mup2*p3mup5*F(2019)+mup4mup5*p2mup2*p3mup3
     -   *F(2020)+mup2mup5*p2mup3*p3mup4*F(2022)+mup2mup4*(p2mup5*p3m
     -   up3+p2mup3*p3mup5)*F(2022)+mup2mup5*p2mup4*p3mup3*F(2023)+mu
     -   p2mup3*p2mup5*p3mup4*F(2023)+mup2mup3*p2mup4*p3mup5*F(2023)+
     -   mup3mup4*p2mup5*p3mup2*F(2024)+mup4mup5*p2mup3*p3mup2*F(2041
     -   )+mup4mup5*p2mup2*p4mup3*F(2042)+2*(mup2mup4*p2mup3*p2mup5*F
     -   (2017)+mup2mup3*p2mup4*p2mup5*F(2018)+mup3mup4*p2mup2*p2mup5
     -   *F(2026)+mup3mup4*p3mup2*p3mup5*F(2028)+mup4mup5*p3mup2*p3mu
     -   p3*F(2029)+mup4mup5*p2mup2*p2mup3*F(2040)+mup4mup5*p4mup2*p4
     -   mup3*F(2044)+mup3mup4*p4mup2*p4mup5*F(2045))+mup3mup5*(p2mup
     -   4*p3mup2*F(2021)+p2mup2*(2*p2mup4*F(2016)+p3mup4*F(2019)+p4m
     -   up4*F(2043))+2*(p3mup2*p3mup4*F(2028)+p4mup2*p4mup4*F(2045))
     -   )
       Fa7(50) = (mup2mup5*mup3mup4*F(2048))/2.+mup3mup4*p2mup2*p4mup5
     -   *F(2049)+mup4mup5*p2mup3*p4mup2*F(2050)+mup2mup5*p2mup3*p4mu
     -   p4*F(2051)+mup3mup5*p2mup4*p4mup2*F(2053)+mup2mup5*p2mup4*p4
     -   mup3*F(2054)+mup3mup4*p2mup5*p4mup2*F(2056)+mup4mup5*p3mup2*
     -   p4mup3*F(2059)+mup3mup5*p3mup2*p4mup4*F(2060)+mup3mup4*p3mup
     -   2*p4mup5*F(2061)+mup4mup5*p3mup3*p4mup2*F(2062)+mup2mup5*p3m
     -   up3*p4mup4*F(2063)+mup3mup5*p3mup4*p4mup2*F(2065)+mup2mup5*p
     -   3mup4*p4mup3*F(2066)+mup3mup4*p3mup5*p4mup2*F(2068)+mup2mup4
     -   *((mup3mup5*F(2047))/2.+p2mup3*p4mup5*F(2052)+p2mup5*p4mup3*
     -   F(2057)+p3mup3*p4mup5*F(2064)+p3mup5*p4mup3*F(2069))+mup2mup
     -   3*((mup4mup5*F(2046))/2.+p2mup4*p4mup5*F(2055)+p2mup5*p4mup4
     -   *F(2058)+p3mup4*p4mup5*F(2067)+p3mup5*p4mup4*F(2070))+4*(p2m
     -   up2*p2mup3*p3mup4*p4mup5*F(2078)+p2mup2*p2mup3*p3mup5*p4mup4
     -   *F(2079)+p2mup2*p2mup4*p3mup3*p4mup5*F(2080)-p2mup2*p2mup4*p
     -   3mup5*p4mup3*F(2081)+p2mup2*p2mup5*p3mup3*p4mup4*F(2082)+p2m
     -   up2*p2mup5*p3mup4*p4mup3*F(2083)+p2mup3*p2mup4*p3mup2*p4mup5
     -   *F(2088)+p2mup3*p2mup4*p3mup5*p4mup2*F(2089)+p2mup3*p2mup5*p
     -   3mup2*p4mup4*F(2090)-p2mup3*p2mup5*p3mup4*p4mup2*F(2091)+p1m
     -   up2*p3mup3*p3mup4*p4mup5*F(2094)+p1mup2*p3mup3*p3mup5*p4mup4
     -   *F(2095)+p1mup3*p3mup2*p3mup4*p4mup5*F(2096)+p1mup3*p3mup2*p
     -   3mup5*p4mup4*F(2097)-p1mup4*p3mup2*p3mup5*p4mup3*F(2098)+p1m
     -   up4*p3mup3*p3mup5*p4mup2*F(2099)+p1mup5*p3mup2*p3mup4*p4mup3
     -   *F(2100)+p1mup5*p3mup3*p3mup4*p4mup2*F(2101)-p2mup4*p2mup5*p
     -   3mup2*p4mup3*F(2103)+p2mup4*p2mup5*p3mup3*p4mup2*F(2104)+p1m
     -   up2*p3mup4*p3mup5*p4mup3*F(2107)+p1mup3*p3mup4*p3mup5*p4mup2
     -   *F(2108)+p1mup2*p3mup3*p4mup4*p4mup5*F(2109)+p1mup2*p3mup4*p
     -   4mup3*p4mup5*F(2110)+p1mup3*p3mup2*p4mup4*p4mup5*F(2111)+p1m
     -   up3*p3mup4*p4mup2*p4mup5*F(2112)+p1mup4*p3mup2*p4mup3*p4mup5
     -   *F(2113)+p1mup4*p3mup3*p4mup2*p4mup5*F(2114)+p1mup4*p3mup2*p
     -   3mup3*p4mup5*F(2118)+p1mup5*p3mup2*p3mup3*p4mup4*F(2119)+p1m
     -   up2*p3mup5*p4mup3*p4mup4*F(2122)+p1mup3*p3mup5*p4mup2*p4mup4
     -   *F(2123)+p1mup4*p3mup5*p4mup2*p4mup3*F(2124)+p1mup5*p3mup2*p
     -   4mup3*p4mup4*F(2125)+p1mup5*p3mup3*p4mup2*p4mup4*F(2126))+8*
     -   (p2mup2*p2mup3*p2mup4*p3mup5*F(2071)+p2mup2*p2mup4*p2mup5*p3
     -   mup3*F(2073)+p2mup2*p2mup3*p2mup4*p4mup5*F(2074)+p2mup2*p2mu
     -   p3*p2mup5*p4mup4*F(2075)+p2mup2*p2mup3*p3mup4*p3mup5*F(2076)
     -   +p2mup2*p2mup5*p3mup3*p3mup4*F(2077)+p2mup3*p2mup4*p2mup5*p3
     -   mup2*F(2085)+p2mup3*p2mup5*p3mup2*p3mup4*F(2086)+p2mup2*p2mu
     -   p4*p3mup3*p3mup5*F(2086)+p2mup3*p2mup4*p3mup2*p3mup5*F(2087)
     -   +p2mup2*p3mup3*p3mup4*p3mup5*F(2092)+p2mup3*p3mup2*p3mup4*p3
     -   mup5*F(2093)+p2mup3*p2mup4*p2mup5*p4mup2*F(2102)+p2mup2*p2mu
     -   p3*p4mup4*p4mup5*F(2105)+p2mup3*p2mup4*p4mup2*p4mup5*F(2106)
     -   +p2mup4*p2mup5*p3mup2*p3mup3*F(2115)+p2mup5*p3mup2*p3mup3*p3
     -   mup4*F(2117)+p2mup3*p2mup5*p4mup2*p4mup4*F(2120)+p2mup2*p2mu
     -   p4*p4mup3*p4mup5*F(2120)+p2mup4*p2mup5*p4mup2*p4mup3*F(2121)
     -   +p2mup2*p2mup5*p4mup3*p4mup4*F(2128))+16*(p2mup2*p2mup3*p2mu
     -   p5*p3mup4*F(2072)+p1mup2*p1mup3*p1mup4*p1mup5*F(2084)+p1mup3
     -   *p1mup4*p1mup5*p2mup2*F(2084)+p1mup2*p1mup4*p1mup5*p2mup3*F(
     -   2084)+p1mup4*p1mup5*p2mup2*p2mup3*F(2084)+p1mup2*p1mup3*p1mu
     -   p5*p2mup4*F(2084)+p1mup3*p1mup5*p2mup2*p2mup4*F(2084)+p1mup2
     -   *p1mup5*p2mup3*p2mup4*F(2084)+p1mup5*p2mup2*p2mup3*p2mup4*F(
     -   2084)+p1mup2*p1mup3*p1mup4*p2mup5*F(2084)+p1mup3*p1mup4*p2mu
     -   p2*p2mup5*F(2084)+p1mup2*p1mup4*p2mup3*p2mup5*F(2084)+p1mup4
     -   *p2mup2*p2mup3*p2mup5*F(2084)+p1mup2*p1mup3*p2mup4*p2mup5*F(
     -   2084)+p1mup3*p2mup2*p2mup4*p2mup5*F(2084)+p1mup2*p2mup3*p2mu
     -   p4*p2mup5*F(2084)+p2mup4*p3mup2*p3mup3*p3mup5*F(2116)+p2mup2
     -   *p2mup4*p2mup5*p4mup3*F(2127)+p4mup2*p4mup3*p4mup4*p4mup5*F(
     -   2129))
       Fa7(51) = mup2mup3*mup4mup5*F(2130)+mup3mup5*(2*p2mup4*p4mup2*F
     -   (2030)+mup2mup4*F(2131)+4*(p3mup4*p4mup2*F(2136)+p3mup2*p4mu
     -   p4*F(2138)))+2*((mup2mup5*p2mup4*p4mup3+mup2mup5*p2mup3*p4mu
     -   p4+mup2mup3*p2mup4*p4mup5)*F(2134)+(mup2mup4*p2mup5*p4mup3+m
     -   up2mup3*p2mup5*p4mup4+mup2mup4*p2mup3*p4mup5)*F(2135)+mup4mu
     -   p5*(p2mup3*p4mup2*F(2133)+2*(p3mup3*p4mup2*F(2136)+p3mup2*p4
     -   mup3*F(2138))))+mup3mup4*(2*(p2mup5*p4mup2*F(2037)+p2mup2*p4
     -   mup5*F(2037))+mup2mup5*F(2132)+4*(p3mup5*p4mup2*F(2137)+p3mu
     -   p2*p4mup5*F(2139)))
       Fa7(52) = 2*(mup3mup5*p2mup4*p4mup2*F(2037)+(mup2mup5*p2mup4*p4
     -   mup3+mup2mup3*p2mup5*p4mup4+mup2mup4*p2mup3*p4mup5)*F(2134)+
     -   (mup2mup4*p2mup5*p4mup3+mup2mup5*p2mup3*p4mup4+mup2mup3*p2mu
     -   p4*p4mup5)*F(2135)+2*mup3mup5*(p3mup4*p4mup2+p3mup2*p4mup4)*
     -   F(2139))+mup4mup5*(2*p2mup3*p4mup2*F(2030)+4*(p3mup3*p4mup2+
     -   p3mup2*p4mup3)*F(2138)+mup2mup3*F(2140))+mup3mup4*(2*(p2mup5
     -   *p4mup2+p2mup2*p4mup5)*F(2030)+4*(p3mup5*p4mup2+p3mup2*p4mup
     -   5)*F(2138)+mup2mup5*F(2140))+mup2mup4*mup3mup5*F(2141)
       Fa7(53) = 4*mup3mup5*p3mup4*p4mup2*F(2136)+mup2mup4*mup3mup5*F(
     -   2143)+mup4mup5*(mup2mup3*F(2142)+2*p2mup3*p4mup2*F(2145)+4*(
     -   p3mup3*p4mup2*F(2137)+p3mup2*p4mup3*F(2150)))+2*(mup2mup4*p2
     -   mup5*p4mup3*F(2146)+mup2mup4*p2mup3*p4mup5*F(2146)+mup3mup5*
     -   p2mup4*p4mup2*F(2147)+mup2mup3*(p2mup5*p4mup4+p2mup4*p4mup5)
     -   *F(2148)+mup2mup5*(p2mup3*p4mup4*F(2146)+p2mup4*p4mup3*F(214
     -   8))+2*mup3mup5*p3mup2*p4mup4*F(2151))+mup3mup4*(mup2mup5*F(2
     -   144)+2*(p2mup2*p4mup5*F(2043)+p2mup5*p4mup2*F(2149))+4*(p3mu
     -   p5*p4mup2*F(2136)+p3mup2*p4mup5*F(2151)))
       Fa7(54) = mup2mup3*p1mup4*p1mup5*F(2156)+mup2mup3*p1mup5*p2mup4
     -   *F(2158)+mup2mup3*p1mup4*p3mup5*F(2163)+mup2mup3*p1mup5*p3mu
     -   p4*F(2164)+mup4mup5*p1mup2*p4mup3*F(2165)+mup3mup4*p1mup2*p4
     -   mup5*F(2167)+mup4mup5*p1mup3*p4mup2*F(2168)+mup2mup4*p1mup3*
     -   p4mup5*F(2170)+mup3mup5*(p1mup2*(2*(p1mup4*F(2152)-p2mup4*F(
     -   2153))+p3mup4*F(2159)+p4mup4*F(2166))+p1mup4*(2*p2mup2*F(215
     -   2)+p3mup2*F(2161)+p4mup2*F(2171)))+mup2mup5*(p1mup3*(-2*(p1m
     -   up4*F(2154)+p2mup4*F(2155))+p3mup4*F(2160)+p4mup4*F(2169))+p
     -   1mup4*(-2*p2mup3*F(2154)+p3mup3*F(2162)+p4mup3*F(2172)))+mup
     -   2mup3*p1mup4*p4mup5*F(2173)+mup3mup4*p1mup5*p4mup2*F(2174)+m
     -   up2mup4*p1mup5*p4mup3*F(2175)+mup2mup3*p1mup5*p4mup4*F(2176)
     -   +mup4mup5*p1mup2*p3mup3*F(2179)+mup3mup4*p1mup2*p3mup5*F(218
     -   0)+mup4mup5*p1mup3*p3mup2*F(2181)+mup2mup4*p1mup3*p3mup5*F(2
     -   182)+mup3mup4*p1mup5*p3mup2*F(2183)+mup2mup4*p1mup5*p3mup3*F
     -   (2184)+2*(mup2mup3*p1mup4*p2mup5*F(2157)+mup4mup5*p1mup2*p1m
     -   up3*F(2177)+mup4mup5*p1mup3*p2mup2*F(2177)+mup4mup5*p1mup2*p
     -   2mup3*F(2177)+mup2mup4*p1mup3*p1mup5*F(2178)+mup2mup4*p1mup5
     -   *p2mup3*F(2178)+mup2mup4*p1mup3*p2mup5*F(2178)+mup3mup4*p1mu
     -   p2*p1mup5*F(2185)+mup3mup4*p1mup5*p2mup2*F(2185)+mup3mup4*p1
     -   mup2*p2mup5*F(2185))+4*p1mup5*p3mup4*p4mup2*p4mup3*F(2193)+1
     -   6*(p1mup2*p1mup3*p1mup5*p3mup4*F(2187)+p1mup3*p1mup5*p2mup2*
     -   p3mup4*F(2187)+p1mup2*p1mup5*p2mup3*p3mup4*F(2187)+p1mup5*p2
     -   mup2*p2mup3*p3mup4*F(2187)+p1mup2*p1mup3*p2mup5*p3mup4*F(218
     -   7)+p1mup3*p2mup2*p2mup5*p3mup4*F(2187)+p1mup2*p2mup3*p2mup5*
     -   p3mup4*F(2187)+p1mup3*p4mup2*p4mup4*p4mup5*F(2195)+p1mup2*p1
     -   mup4*p1mup5*p4mup3*F(2197)+p1mup4*p1mup5*p2mup2*p4mup3*F(219
     -   7)+p1mup2*p1mup5*p2mup4*p4mup3*F(2197)+p1mup2*p1mup4*p2mup5*
     -   p4mup3*F(2197)+p1mup4*p2mup2*p2mup5*p4mup3*F(2197)+p1mup2*p2
     -   mup4*p2mup5*p4mup3*F(2197))+8*(p1mup2*p1mup3*p1mup4*p3mup5*F
     -   (2186)+p1mup3*p1mup4*p2mup2*p3mup5*F(2186)+p1mup2*p1mup4*p2m
     -   up3*p3mup5*F(2186)+p1mup4*p2mup2*p2mup3*p3mup5*F(2186)+p1mup
     -   2*p1mup3*p2mup4*p3mup5*F(2186)+p1mup3*p2mup2*p2mup4*p3mup5*F
     -   (2186)+p1mup2*p2mup3*p2mup4*p3mup5*F(2186)+p1mup2*p1mup4*p1m
     -   up5*p3mup3*F(2188)+p1mup4*p1mup5*p2mup2*p3mup3*F(2188)+p1mup
     -   2*p1mup5*p2mup4*p3mup3*F(2188)+p1mup5*p2mup2*p2mup4*p3mup3*F
     -   (2188)+p1mup2*p1mup4*p2mup5*p3mup3*F(2188)+p1mup4*p2mup2*p2m
     -   up5*p3mup3*F(2188)+p1mup2*p2mup4*p2mup5*p3mup3*F(2188)+p1mup
     -   2*p1mup3*p1mup4*p4mup5*F(2189)+p1mup3*p1mup4*p2mup2*p4mup5*F
     -   (2189)+p1mup2*p1mup4*p2mup3*p4mup5*F(2189)+p1mup4*p2mup2*p2m
     -   up3*p4mup5*F(2189)+p1mup2*p1mup3*p2mup4*p4mup5*F(2189)+p1mup
     -   3*p2mup2*p2mup4*p4mup5*F(2189)+p1mup2*p2mup3*p2mup4*p4mup5*F
     -   (2189)+p1mup2*p1mup3*p1mup5*p4mup4*F(2190)+p1mup3*p1mup5*p2m
     -   up2*p4mup4*F(2190)+p1mup2*p1mup5*p2mup3*p4mup4*F(2190)+p1mup
     -   2*p1mup3*p2mup5*p4mup4*F(2190)+p1mup3*p2mup2*p2mup5*p4mup4*F
     -   (2190)+p1mup2*p2mup3*p2mup5*p4mup4*F(2190)+p1mup3*p1mup4*p1m
     -   up5*p3mup2*F(2191)+p1mup4*p1mup5*p2mup3*p3mup2*F(2191)+p1mup
     -   3*p1mup5*p2mup4*p3mup2*F(2191)+p1mup5*p2mup3*p2mup4*p3mup2*F
     -   (2191)+p1mup3*p1mup4*p2mup5*p3mup2*F(2191)+p1mup4*p2mup3*p2m
     -   up5*p3mup2*F(2191)+p1mup3*p2mup4*p2mup5*p3mup2*F(2191)+p1mup
     -   3*p1mup4*p1mup5*p4mup2*F(2192)+p1mup4*p1mup5*p2mup3*p4mup2*F
     -   (2192)+p1mup3*p1mup5*p2mup4*p4mup2*F(2192)+p1mup3*p1mup4*p2m
     -   up5*p4mup2*F(2192)+p1mup4*p2mup3*p2mup5*p4mup2*F(2192)+p1mup
     -   3*p2mup4*p2mup5*p4mup2*F(2192)+p1mup2*p4mup3*p4mup4*p4mup5*F
     -   (2194)+p1mup4*p4mup2*p4mup3*p4mup5*F(2196)+p1mup5*p4mup2*p4m
     -   up3*p4mup4*F(2198))
       Fa7(55) = p1mup2*(mup3mup5*p3mup4-mup3mup4*p3mup5)*F(2201)+mup2
     -   mup5*p1mup4*p3mup3*F(2203)-mup2mup4*p1mup5*p3mup3*F(2203)+mu
     -   p2mup5*p1mup3*p3mup4*F(2203)-mup2mup3*p1mup5*p3mup4*F(2203)-
     -   mup2mup4*p1mup3*p3mup5*F(2203)+mup2mup3*p1mup4*p3mup5*F(2203
     -   )+mup3mup5*p1mup4*p3mup2*F(2204)+mup3mup4*p1mup5*p3mup2*F(22
     -   05)-2*(mup2mup4*p1mup3*p1mup5*F(2200)-mup2mup5*p1mup4*p2mup3
     -   *F(2200)+mup2mup4*p1mup5*p2mup3*F(2200)-mup2mup5*p1mup3*(p1m
     -   up4+p2mup4)*F(2200)+mup2mup4*p1mup3*p2mup5*F(2200)-mup3mup5*
     -   p1mup2*p1mup4*F(2206)+mup3mup4*p1mup2*p1mup5*F(2206)-mup3mup
     -   5*p1mup4*p2mup2*F(2206)+mup3mup4*p1mup5*p2mup2*F(2206)-mup3m
     -   up5*p1mup2*p2mup4*F(2206)+mup3mup4*p1mup2*p2mup5*F(2206))+mu
     -   p3mup5*p1mup4*p4mup2*F(2207)-mup3mup4*p1mup5*p4mup2*F(2207)+
     -   mup3mup5*p1mup2*p4mup4*F(2207)-mup3mup4*p1mup2*p4mup5*F(2207
     -   )+mup4mup5*(p1mup2*(2*(p1mup3+p2mup3)*F(2199)+p3mup3*F(2201)
     -   +p4mup3*F(2207))+p1mup3*(2*p2mup2*F(2199)+p3mup2*F(2202)+p4m
     -   up2*F(2208)))+(mup2mup5*p1mup4*p4mup3-mup2mup4*p1mup5*p4mup3
     -   +mup2mup5*p1mup3*p4mup4-mup2mup3*p1mup5*p4mup4-mup2mup4*p1mu
     -   p3*p4mup5+mup2mup3*p1mup4*p4mup5)*F(3789)
       Fa7(56) = 2*(mup3mup4*p1mup2*p1mup5*F(2206)+mup3mup4*p1mup5*p2m
     -   up2*F(2206)+mup3mup4*p1mup2*p2mup5*F(2206))+mup3mup4*p1mup5*
     -   p4mup2*F(2207)+mup3mup4*p1mup2*p4mup5*F(2207)+mup3mup4*p1mup
     -   5*p3mup2*F(2209)+mup3mup4*p1mup2*p3mup5*F(2209)+mup4mup5*(p1
     -   mup3*(2*p2mup2*F(2206)+p4mup2*F(2207)+p3mup2*F(2209))+p1mup2
     -   *(2*(p1mup3+p2mup3)*F(2206)+p4mup3*F(2207)+p3mup3*F(2209)))-
     -   mup3mup5*(p1mup4*(2*p2mup2*F(2206)+p4mup2*F(2207)+p3mup2*F(2
     -   209))+p1mup2*(2*(p1mup4+p2mup4)*F(2206)+p4mup4*F(2207)+p3mup
     -   4*F(2209)))+mup2mup5*p1mup4*p3mup3*F(3788)-mup2mup4*p1mup5*p
     -   3mup3*F(3788)-mup2mup5*p1mup3*p3mup4*F(3788)+mup2mup3*p1mup5
     -   *p3mup4*F(3788)+mup2mup4*p1mup3*p3mup5*F(3788)-mup2mup3*p1mu
     -   p4*p3mup5*F(3788)+(mup2mup5*p1mup4*p4mup3-mup2mup4*p1mup5*p4
     -   mup3-mup2mup5*p1mup3*p4mup4+mup2mup3*p1mup5*p4mup4+mup2mup4*
     -   p1mup3*p4mup5-mup2mup3*p1mup4*p4mup5)*F(3789)
       Return
       End
