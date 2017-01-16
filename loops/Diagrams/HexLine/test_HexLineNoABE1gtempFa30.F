c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g30(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa30)
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
       Complex*16 Fa30(407:420)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa30(407) = mup2mup4*(p3mup6*F(2043)+4*(p1mup6*F(2041)+p2mup6*F
     -   (2042)-p4mup6*F(2044)-p5mup6*F(2045)))+8*p5mup4*(-(p4mup2*(p
     -   4mup6*F(2046)+p5mup6*F(2047)))+p5mup2*(p4mup6*F(2048)+p5mup6
     -   *F(2049)))
       Fa30(408) = p1mup4*p2mup6*p3mup2*F(2100)-8*(mup2mup6*(p1mup4*F(
     -   2050)+p2mup4*F(2051)-p3mup4*F(2052)-p4mup4*F(2053))+p2mup6*p
     -   3mup2*p3mup4*F(2055)-p2mup6*p3mup2*p4mup4*F(2056)-p2mup6*p3m
     -   up4*p4mup2*F(2057)+p1mup4*p3mup2*p4mup6*F(2058)-p1mup4*p3mup
     -   6*p4mup2*F(2059)+p3mup2*p4mup4*p4mup6*F(2060)-p3mup4*p3mup6*
     -   p4mup2*F(2061)+p3mup2*p3mup4*p5mup6*F(2063)+p3mup2*p4mup4*p5
     -   mup6*F(2064)-p3mup4*p4mup2*p5mup6*F(2065)+p3mup2*p3mup6*p4mu
     -   p4*F(2066)-p3mup4*p4mup2*p4mup6*F(2067)+p1mup6*p3mup2*p3mup4
     -   *F(2069)+p1mup6*p3mup2*p4mup4*F(2070)+p1mup6*p3mup4*p4mup2*F
     -   (2071)-p1mup4*p2mup6*p4mup2*F(2072)+p3mup2*p3mup4*p4mup6*F(2
     -   079)+p3mup2*p3mup4*p3mup6*F(2081)-p3mup6*p4mup2*p4mup4*F(208
     -   2)+p1mup4*p4mup2*p4mup6*F(2083)-p1mup4*p2mup6*p5mup2*F(2084)
     -   -p2mup6*p4mup4*p5mup2*F(2086)+p3mup4*p3mup6*p5mup2*F(2088)+p
     -   3mup6*p4mup4*p5mup2*F(2089)-p1mup4*p4mup6*p5mup2*F(2090)-p4m
     -   up4*p4mup6*p5mup2*F(2092)-p4mup4*p5mup2*p5mup6*F(2094)+p2mup
     -   4*(p3mup6*p4mup2*F(2075)+p4mup2*p5mup6*F(2076)+p4mup2*p4mup6
     -   *F(2077)+p1mup6*p4mup2*F(2078)+p3mup2*(p5mup6*F(2062)+p1mup6
     -   *F(2068)+p4mup6*F(2074)+p3mup6*F(2080))+p2mup6*(p3mup2*F(205
     -   4)+p4mup2*F(2073)-p5mup2*F(2085))+p3mup6*p5mup2*F(2087)-p4mu
     -   p6*p5mup2*F(2091)-p5mup2*p5mup6*F(2093)+p1mup6*p5mup2*F(2095
     -   ))+p1mup6*p3mup4*p5mup2*F(2096)+p1mup6*p4mup4*p5mup2*F(2097)
     -   -p3mup4*p4mup6*p5mup2*F(2098)+p1mup4*p3mup2*p3mup6*F(2099)-p
     -   2mup6*p3mup4*p5mup2*F(2101)-p1mup4*p3mup6*p5mup2*F(2102)-p3m
     -   up4*p5mup2*p5mup6*F(2103)+p1mup4*p3mup2*p5mup6*F(2104)-p1mup
     -   4*p4mup2*p5mup6*F(2105)+p1mup4*p1mup6*p3mup2*F(2106)-p1mup4*
     -   p1mup6*p4mup2*F(2107)-p2mup6*p4mup2*p4mup4*F(2108)-p4mup2*p4
     -   mup4*p5mup6*F(2109)-p4mup2*p4mup4*p4mup6*F(2110)-p1mup6*p4mu
     -   p2*p4mup4*F(2111)-p1mup4*p5mup2*p5mup6*F(2112)-p1mup4*p1mup6
     -   *p5mup2*F(2113))
       Fa30(409) = 4*mup2mup3*(p2mup6*F(2120)+p3mup6*F(2122)+p4mup6*F(
     -   2124))+mup3mup6*(p1mup2*F(2114)+p2mup2*F(2115)+p3mup2*F(2121
     -   )+p4mup2*F(2123)+p5mup2*F(2125))+p2mup2*p2mup3*p4mup6*F(2132
     -   )+p2mup3*p2mup6*p4mup2*F(2133)+p2mup2*p3mup3*p3mup6*F(2144)+
     -   mup2mup6*(p2mup3*F(2116)+p4mup3*F(2117)+4*p1mup3*F(2118)+p3m
     -   up3*F(2119)+p5mup3*F(2148))+p2mup2*p3mup6*p5mup3*F(2171)-8*(
     -   p3mup3*p4mup2*p4mup6*F(2127)+p2mup2*p3mup3*p4mup6*F(2128)-p2
     -   mup2*p3mup6*p4mup3*F(2129)+p2mup6*p3mup3*p4mup2*F(2131)+p3mu
     -   p3*p3mup6*p4mup2*F(2135)+p2mup2*p2mup3*p2mup6*F(2136)-p2mup2
     -   *p4mup3*p4mup6*F(2137)+p2mup3*p4mup2*p4mup6*F(2138)+p2mup3*p
     -   3mup6*p4mup2*F(2140)+p2mup2*p2mup3*p3mup6*F(2141)-p2mup2*p2m
     -   up6*p3mup3*F(2142)-p2mup2*p2mup6*p4mup3*F(2143)-p2mup6*p4mup
     -   2*p4mup3*F(2146)-p3mup6*p4mup2*p4mup3*F(2147)-p4mup2*p4mup6*
     -   p5mup3*F(2149)+p4mup3*p4mup6*p5mup2*F(2150)+p2mup3*p3mup6*p5
     -   mup2*F(2151)+p2mup6*p3mup3*p5mup2*F(2153)+p2mup2*p4mup6*p5mu
     -   p3*F(2154)+p2mup3*p4mup6*p5mup2*F(2155)-p2mup6*p4mup2*p5mup3
     -   *F(2156)+p2mup6*p4mup3*p5mup2*F(2157)+p3mup3*p4mup6*p5mup2*F
     -   (2158)-p3mup6*p4mup2*p5mup3*F(2159)-p3mup6*p4mup3*p5mup2*F(2
     -   160)-p2mup2*p2mup6*p5mup3*F(2161)+p2mup3*p2mup6*p5mup2*F(216
     -   2)+p4mup2*p4mup3*p4mup6*F(2166)+p2mup6*p5mup2*p5mup3*F(2167)
     -   +p3mup6*p5mup2*p5mup3*F(2168)-p4mup6*p5mup2*p5mup3*F(2169)+p
     -   3mup3*p3mup6*p5mup2*F(2173))+p3mup2*(p4mup3*(p2mup6*F(2130)+
     -   8*(p4mup6*F(2126)+p3mup6*F(2134)))+p3mup3*p4mup6*F(2145)+p3m
     -   up3*p3mup6*F(2163)+8*(p2mup6*p5mup3*F(2152)-p2mup6*p3mup3*F(
     -   2165)-p2mup3*(p4mup6*F(2139)+p3mup6*F(2164)+p2mup6*F(2170))+
     -   p3mup6*p5mup3*F(2172)+p4mup6*p5mup3*F(2174)))
       Fa30(410) = 4*mup2mup3*p1mup6*F(2175)+p1mup2*p3mup3*p4mup6*F(21
     -   80)+p1mup3*p2mup2*p2mup6*F(2186)+p1mup2*p3mup3*p3mup6*F(2221
     -   )+p1mup6*(p3mup2*p4mup3*F(2177)+p2mup2*p4mup3*F(2185)+p4mup2
     -   *p4mup3*F(2191)+p2mup2*p5mup3*F(2196)+p2mup3*(p4mup2*F(2176)
     -   +8*(p2mup2*F(2184)+p5mup2*F(2197)+p3mup2*F(2204)-p1mup2*F(22
     -   11)))-8*(p4mup2*p5mup3*F(2194)-p4mup3*p5mup2*F(2195)-p5mup2*
     -   p5mup3*F(2203)+p3mup2*p5mup3*F(2206)+p1mup2*p1mup3*F(2209)+p
     -   1mup3*p2mup2*F(2213)-p3mup3*(p4mup2*F(2178)+p2mup2*F(2189)+p
     -   3mup2*F(2201)+p5mup2*F(2207)-p1mup2*F(2215))+p1mup2*p4mup3*F
     -   (2219)+p1mup3*p4mup2*F(2220)+p1mup3*p3mup2*F(2223)-p1mup2*p5
     -   mup3*F(2224)+p1mup3*p5mup2*F(2225)))-8*(p1mup2*p2mup3*p4mup6
     -   *F(2179)-p1mup2*p3mup6*p4mup3*F(2181)-p1mup3*p3mup2*p4mup6*F
     -   (2182)+p1mup3*p3mup6*p4mup2*F(2183)-p1mup2*p2mup6*p4mup3*F(2
     -   187)+p1mup3*p2mup6*p4mup2*F(2188)+p1mup3*p2mup2*p3mup6*F(219
     -   0)+p1mup3*p2mup2*p4mup6*F(2192)+p1mup3*p4mup2*p4mup6*F(2193)
     -   -p1mup2*p2mup6*p5mup3*F(2198)+p1mup3*p2mup6*p5mup2*F(2199)+p
     -   1mup3*p4mup6*p5mup2*F(2200)+p1mup3*p3mup2*p3mup6*F(2202)+p1m
     -   up3*p2mup6*p3mup2*F(2205)+p1mup3*p3mup6*p5mup2*F(2208)+p1mup
     -   2*p1mup3*p2mup6*F(2210)+p1mup2*p2mup3*p2mup6*F(2212)+p1mup2*
     -   p1mup3*p3mup6*F(2214)+p1mup2*p2mup3*p3mup6*F(2216)+p1mup2*p2
     -   mup6*p3mup3*F(2217)-p1mup2*p1mup3*p4mup6*F(2218)-p1mup2*p4mu
     -   p3*p4mup6*F(2222)-p1mup2*p3mup6*p5mup3*F(2226)-p1mup2*p4mup6
     -   *p5mup3*F(2227))
       Fa30(411) = p5mup6*(p1mup2*(-(p1mup3*F(2230))+p2mup3*F(2231)+p3
     -   mup3*F(2232))+p2mup2*(p2mup3*F(2228)+p3mup3*F(2229)-p1mup3*F
     -   (2233)))
       Fa30(412) = p5mup6*(-(p4mup3*(p3mup2*F(2236)+p2mup2*F(2237)))+p
     -   4mup2*(p2mup3*F(2234)+p3mup3*F(2235)+p4mup3*F(2238)-p5mup3*F
     -   (2239)+p1mup3*F(2240))-p1mup2*p4mup3*F(2241))
       Fa30(413) = 4*p5mup6*(mup2mup3*F(2242)-2*(p2mup3*p5mup2*F(2244)
     -   +p4mup3*p5mup2*F(2245)-p3mup2*p3mup3*F(2246)+p2mup3*p3mup2*F
     -   (2247)+p3mup3*p5mup2*F(2248)+p1mup3*p3mup2*F(2249)+p1mup3*p5
     -   mup2*F(2252)-p5mup3*(p2mup2*F(2243)+p1mup2*F(2250)+p3mup2*F(
     -   2251)+p5mup2*F(2253))))
       Fa30(414) = p5mup4*F(2254)+p4mup4*F(2255)+p2mup4*F(2256)+p1mup4
     -   *F(2257)+p3mup4*F(2258)
       Fa30(415) = p5mup3*F(2259)+p3mup3*F(2262)+2*(p2mup3*F(2260)+p4m
     -   up3*F(2261)+p1mup3*F(2263))
       Fa30(416) = p2mup2*F(2264)+p4mup2*F(2265)+p5mup2*F(2266)+p3mup2
     -   *F(2267)+p1mup2*F(2268)
       Fa30(417) = mup4mup6*(p2mup2*p2mup3*F(2271)-p3mup2*(p4mup3*F(22
     -   69)+p5mup3*F(2273))+p3mup3*(p4mup2*F(2270)+p5mup2*F(2274))-p
     -   4mup2*(p4mup3*F(2272)+p5mup3*F(2275))+p5mup2*p5mup3*F(2276))
       Fa30(418) = mup4mup6*(p2mup3*(p4mup2*F(2277)+8*(p5mup2*F(2280)+
     -   p3mup2*F(2284)+p1mup2*F(2288)))+8*(p4mup3*p5mup2*F(2282)-mup
     -   2mup3*F(2283)+p3mup2*p3mup3*F(2285)+p2mup2*(p3mup3*F(2278)+p
     -   1mup3*F(2279)+p4mup3*F(2281)+p5mup3*F(2286))+p1mup2*p1mup3*F
     -   (2287)+p1mup2*p3mup3*F(2289)))
       Fa30(419) = mup4mup6*(p1mup2*p4mup3*F(2290)+p1mup3*p4mup2*F(229
     -   1)+p1mup2*p5mup3*F(2292)+p1mup3*p5mup2*F(2293)+p1mup3*p3mup2
     -   *F(2294))
       Fa30(420) = mup3mup6*(p3mup2*(2*(p3mup4*F(2295)+p4mup4*F(2296))
     -   +p5mup4*F(2299))+2*(p1mup2*(p3mup4*F(2297)+p4mup4*F(2298)+p5
     -   mup4*F(2300))+p4mup2*(p5mup4*F(2303)+p3mup4*F(2304)-p4mup4*F
     -   (2305))-p2mup2*(p3mup4*F(2301)+p4mup4*F(2302)+p5mup4*F(2306)
     -   )))
       Return
       End
