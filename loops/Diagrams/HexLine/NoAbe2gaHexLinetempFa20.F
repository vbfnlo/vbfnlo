c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa20(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa20)
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
       Complex*16 Fa20(457:480)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa20(457) = p4mup6*(p1mup3*(p5mup5*F(2085)-2*(p4mup5*F(2096)+(p
     -   1mup5+p2mup5)*F(2129)+p3mup5*F(2135)))+2*p1mup5*(p4mup3*F(21
     -   09)+p5mup3*F(2113)-p2mup3*F(2129)+p3mup3*F(2137)))
       Fa20(458) = 2*(p3mup6*p4mup3*p5mup5*F(2093)-p4mup3*p5mup5*p5mup
     -   6*F(2099))+p2mup6*p4mup3*p5mup5*F(2102)-4*(p4mup3*p4mup5*p5m
     -   up6*F(2088)+p4mup5*p5mup3*p5mup6*F(2090)-p3mup6*p4mup3*p4mup
     -   5*F(2092)-p3mup6*p4mup5*p5mup3*F(2094)-p1mup6*p4mup3*p4mup5*
     -   F(2101)-p2mup6*p4mup3*p4mup5*F(2101)-p1mup6*p4mup5*p5mup3*F(
     -   2103)-p2mup6*p4mup5*p5mup3*F(2103)+p3mup3*p4mup5*p5mup6*F(21
     -   05)+p3mup5*p4mup3*p5mup6*F(2106)-p2mup5*p4mup3*p5mup6*F(2111
     -   )-p3mup5*p3mup6*p4mup3*F(2118)+p3mup3*p3mup6*p4mup5*F(2121)-
     -   p2mup5*p2mup6*p4mup3*F(2131)-p1mup3*p4mup5*(p5mup6*F(2087)+(
     -   p1mup6+p2mup6)*F(2128)+p3mup6*F(2136))-p2mup3*p4mup5*(p5mup6
     -   *F(2087)+(p1mup6+p2mup6)*F(2128)+p3mup6*F(2136))+p2mup5*p3mu
     -   p6*p4mup3*F(2138)-p1mup6*p3mup3*p4mup5*F(2139)-p2mup6*p3mup3
     -   *p4mup5*F(2139)-p2mup6*p3mup5*p4mup3*F(2140))
       Fa20(459) = -2*p5mup3*p5mup5*p5mup6*F(2100)+p1mup6*p4mup3*p5mup
     -   5*F(2102)+p1mup6*p5mup3*p5mup5*F(2104)-p1mup6*p3mup3*p5mup5*
     -   F(2144)+p1mup3*(-(p5mup5*((p1mup6+p2mup6)*F(2125)+2*(p5mup6*
     -   F(2086)+p3mup6*F(2143))))+4*((p1mup5+p2mup5)*((p1mup6+p2mup6
     -   )*F(2117)+p3mup6*F(2126))+(p1mup5+p2mup5)*p5mup6*F(2130)+p3m
     -   up5*((p1mup6+p2mup6)*F(2127)+p3mup6*F(2134)+p5mup6*F(2145)))
     -   )+p2mup3*(-(p5mup5*((p1mup6+p2mup6)*F(2125)+2*(p5mup6*F(2086
     -   )+p3mup6*F(2143))))+4*((p1mup5+p2mup5)*((p1mup6+p2mup6)*F(21
     -   17)+p3mup6*F(2126))+(p1mup5+p2mup5)*p5mup6*F(2130)+p3mup5*((
     -   p1mup6+p2mup6)*F(2127)+p3mup6*F(2134)+p5mup6*F(2145))))+4*(p
     -   1mup5*p4mup3*p5mup6*F(2111)+p1mup5*p5mup3*p5mup6*F(2112)-p3m
     -   up3*p3mup5*p3mup6*F(2122)+p1mup5*p1mup6*p4mup3*F(2131)+p1mup
     -   6*p2mup5*p4mup3*F(2131)+p1mup5*p2mup6*p4mup3*F(2131)+p1mup5*
     -   p1mup6*p5mup3*F(2132)+p1mup6*p2mup5*p5mup3*F(2132)+p1mup5*p2
     -   mup6*p5mup3*F(2132)+p2mup5*p2mup6*p5mup3*F(2132)-p1mup5*p1mu
     -   p6*p3mup3*F(2133)-p1mup6*p2mup5*p3mup3*F(2133)-p1mup5*p2mup6
     -   *p3mup3*F(2133)-p2mup5*p2mup6*p3mup3*F(2133)-p1mup5*p3mup6*p
     -   4mup3*F(2138)+p1mup6*p3mup5*p4mup3*F(2140)+p1mup5*p3mup3*p3m
     -   up6*F(2141)+p2mup5*p3mup3*p3mup6*F(2141)-p1mup6*p3mup3*p3mup
     -   5*F(2142)-p2mup6*p3mup3*p3mup5*F(2142)+p1mup5*p3mup3*p5mup6*
     -   F(2146)+p2mup5*p3mup3*p5mup6*F(2146)+p1mup5*p3mup6*p5mup3*F(
     -   2147)+p1mup6*p3mup5*p5mup3*F(2148))
       Fa20(460) = -2*(p5mup6*(p3mup3*p5mup5*F(2107)+2*(p3mup5*p5mup3*
     -   F(2108)-p2mup5*p5mup3*F(2112)-p3mup3*p3mup5*F(2123)))-p3mup6
     -   *(p3mup3*p5mup5*F(2119)+p5mup3*(p5mup5*F(2095)+2*(p3mup5*F(2
     -   120)+p2mup5*F(2147)))))+p2mup6*(p5mup3*p5mup5*F(2104)-p3mup3
     -   *p5mup5*F(2144)+4*p3mup5*p5mup3*F(2148))
       Fa20(461) = p6mup3*(p5mup5*p5mup6*F(2218)-2*(p1mup6*p3mup5*F(22
     -   08)+(p1mup6+p2mup6)*p4mup5*F(2210)-p3mup6*p4mup5*F(2214)-p1m
     -   up5*(p3mup6*F(2207)+p4mup6*F(2209)+(p1mup6+p2mup6)*F(2216))-
     -   p2mup5*(p3mup6*F(2207)+p4mup6*F(2209)+(p1mup6+p2mup6)*F(2216
     -   ))+p3mup5*p4mup6*F(9606)))
       Fa20(462) = -(p1mup6*p5mup5*p6mup3*F(2206))-p2mup6*p5mup5*p6mup
     -   3*F(2206)-p2mup3*p5mup5*p6mup6*F(2256)-2*(p3mup6*p5mup5*p6mu
     -   p3*F(2215)+p4mup6*p5mup5*p6mup3*F(2217)-p1mup5*p2mup3*p6mup6
     -   *F(2257)-p2mup3*p2mup5*p6mup6*F(2257)-p2mup3*p3mup5*p6mup6*F
     -   (2258)-p4mup5*(2*(p4mup6*p6mup3*F(1967)+p5mup6*p6mup3*F(2212
     -   ))+p2mup3*p6mup6*F(2259)))+4*(p2mup6*p3mup5*p6mup3*F(2208)-p
     -   1mup5*p5mup6*p6mup3*F(2211)-p2mup5*p5mup6*p6mup3*F(2211)+p1m
     -   up6*(p1mup3+p2mup3)*p6mup5*F(2231)+p1mup3*p2mup6*p6mup5*F(22
     -   31)+p2mup3*p2mup6*p6mup5*F(2231)+p2mup3*p3mup6*p6mup5*F(2232
     -   )+p2mup3*p4mup6*p6mup5*F(2233)+p2mup3*p5mup6*p6mup5*F(2234)+
     -   p4mup3*p4mup6*p6mup5*F(2238)+p3mup5*p5mup6*p6mup3*F(9590)-p3
     -   mup3*p3mup6*p6mup5*F(9590)+p3mup5*p3mup6*p6mup3*F(9599)-p5mu
     -   p3*p5mup6*p6mup5*F(9602))
       Fa20(463) = p6mup5*(-(p1mup3*(p3mup6*F(2232)+p4mup6*F(2233)+p5m
     -   up6*F(2234)))-(p1mup6+p2mup6)*(p3mup3*F(2235)-p4mup3*F(2236)
     -   -p5mup3*F(2237))+p4mup3*p5mup6*F(2239)+p3mup6*p4mup3*F(2240)
     -   +p3mup6*p5mup3*F(9594)+p3mup3*p5mup6*F(9594)+p4mup6*p5mup3*F
     -   (9596)+p3mup3*p4mup6*F(9601))
       Fa20(464) = (p1mup6+p2mup6)*p6mup3*p6mup5
       Fa20(465) = p6mup6*(p2mup5*p5mup3*F(2262)+p2mup3*p6mup5*F(2281)
     -   )
       Fa20(466) = p6mup6*(p1mup3*(-(p5mup5*F(2256))+2*((p1mup5+p2mup5
     -   )*F(2257)+p3mup5*F(2258)+p4mup5*F(2259)))+p4mup3*p5mup5*F(22
     -   67)+p5mup3*p5mup5*F(2269)-p3mup3*p5mup5*F(2275)-6*(p3mup5*p5
     -   mup3*F(2276)+p3mup3*p3mup5*F(2279))+2*(p1mup5*p3mup3*F(2260)
     -   +p2mup5*p3mup3*F(2260)-p1mup5*p5mup3*F(2262)+p4mup3*p4mup5*F
     -   (2266)+p4mup5*p5mup3*F(2268)-p3mup3*p4mup5*F(2273)-p3mup5*p4
     -   mup3*F(2274)+(p1mup5+p2mup5)*p4mup3*F(2471)))
       Fa20(467) = p6mup3*p6mup6*((p1mup5+p2mup5)*F(2280)+p4mup5*F(228
     -   5)+p5mup5*F(2286)-3*p3mup5*F(9581))
       Fa20(468) = p6mup5*p6mup6*(-(p1mup3*F(2281))+p4mup3*F(2291)+3*(
     -   p5mup3*F(9583)+p3mup3*F(9588)))
       Fa20(469) = mup5mup6*(-(p3mup3*F(2061))+p2mup3*F(2063)-p5mup3*F
     -   (2472)+p4mup3*F(2473))
       Fa20(470) = -2*mup3mup5*(p6mup6*F(2404)+2*((p1mup6+p2mup6)*F(24
     -   03)-p5mup6*F(2405)))+mup3mup6*((p1mup5+p2mup5)*F(2474)+p5mup
     -   5*F(2475)+p6mup5*F(2476))
       Fa20(471) = 2*mup5mup6*p6mup3*F(2222)+mup3mup6*p4mup5*F(2477)
       Fa20(472) = p4mup6*(2*(p5mup3*(p5mup5*F(2166)+2*(p2mup5*F(2168)
     -   -p3mup5*F(2181)+p4mup5*F(2184)))+p4mup3*(p5mup5*F(2163)+2*(p
     -   4mup5*F(2151)-p3mup5*F(2173)-p2mup5*F(2185)))-p2mup3*(-(p5mu
     -   p5*F(2167))+2*(p4mup5*F(2183)+p2mup5*F(2189)+p3mup5*F(2193))
     -   )-p3mup3*(p5mup5*F(2179)+2*(p4mup5*F(2197)+p2mup5*F(2199)-p3
     -   mup5*F(2220))))+mup3mup5*F(2409))
       Fa20(473) = p4mup6*(p1mup3*(-(p5mup5*F(2167))+2*(p4mup5*F(2183)
     -   +(p1mup5+p2mup5)*F(2189)+p3mup5*F(2193)))+2*p1mup5*(-(p5mup3
     -   *F(2168))+p4mup3*F(2185)+p2mup3*F(2189)+p3mup3*F(2199)))
       Fa20(474) = p4mup3*(p5mup5*(p5mup6*F(2164)+2*(p3mup6*F(2182)+p2
     -   mup6*F(2187)))+4*(-(p3mup5*p5mup6*F(2180))-p2mup5*p5mup6*F(2
     -   186)+p2mup5*p2mup6*F(2198)-p2mup6*p3mup5*F(2201)-p2mup5*p3mu
     -   p6*F(2204)+p4mup5*(p5mup6*F(2162)-(p1mup6+p2mup6)*F(2196)-p3
     -   mup6*F(2205))-p3mup5*p3mup6*F(2221)))+4*p4mup5*(p2mup3*p5mup
     -   6*F(2172)-p3mup3*p5mup6*F(2178)+p1mup6*p2mup3*F(2192)+p2mup3
     -   *p2mup6*F(2192)+p2mup3*p3mup6*F(2194)+p1mup3*(p5mup6*F(2172)
     -   +(p1mup6+p2mup6)*F(2192)+p3mup6*F(2194))-p1mup6*p3mup3*F(220
     -   0)-p2mup6*p3mup3*F(2200)+p5mup3*(p5mup6*F(2165)+(p1mup6+p2mu
     -   p6)*F(2169)+p3mup6*F(2203))+p3mup3*p3mup6*F(9615))
       Fa20(475) = -(p2mup3*p5mup5*p5mup6*F(2170))-2*(p1mup6*p2mup3*p5
     -   mup5*F(2158)+p2mup3*p2mup6*p5mup5*F(2158)+p2mup3*p3mup6*p5mu
     -   p5*F(2159)+p1mup6*p3mup3*p5mup5*F(2177)-p1mup6*p4mup3*p5mup5
     -   *F(2187))+p5mup3*(p5mup5*(p5mup6*F(2149)-p1mup6*F(2157))+4*(
     -   p1mup6*p3mup5*F(2161)-p1mup5*(p3mup6*F(2160)+p5mup6*F(2171))
     -   -(p1mup5+p2mup5)*(p1mup6+p2mup6)*F(2191)))+p1mup3*(-2*(p1mup
     -   6*p5mup5*F(2158)+p2mup6*p5mup5*F(2158)+p3mup6*p5mup5*F(2159)
     -   )-p5mup5*p5mup6*F(2170)+4*(p1mup6*(p1mup5+p2mup5)*F(2188)+(p
     -   1mup5+p2mup5)*p2mup6*F(2188)+p1mup5*p5mup6*F(2190)+p2mup5*p5
     -   mup6*F(2190)+p3mup5*p5mup6*F(2195)+p3mup5*p3mup6*F(2213)))+4
     -   *((p1mup3+p2mup3)*((p1mup6+p2mup6)*p3mup5*F(2150)+(p1mup5+p2
     -   mup5)*p3mup6*F(2152))-(p1mup5+p2mup5)*(p1mup6+p2mup6)*p3mup3
     -   *F(2156)-p1mup5*p4mup3*p5mup6*F(2186)+p1mup5*p1mup6*p2mup3*F
     -   (2188)+p1mup6*p2mup3*p2mup5*F(2188)+p1mup5*p2mup3*p2mup6*F(2
     -   188)+p2mup3*p2mup5*p2mup6*F(2188)+p1mup5*p2mup3*p5mup6*F(219
     -   0)+p2mup3*p2mup5*p5mup6*F(2190)+p2mup3*p3mup5*p5mup6*F(2195)
     -   +p1mup5*p1mup6*p4mup3*F(2198)+p1mup6*p2mup5*p4mup3*F(2198)+p
     -   1mup5*p2mup6*p4mup3*F(2198)-p1mup6*p3mup5*p4mup3*F(2201)-p1m
     -   up5*p3mup3*p5mup6*F(2202)-p2mup5*p3mup3*p5mup6*F(2202)-p1mup
     -   5*p3mup6*p4mup3*F(2204)+p2mup3*p3mup5*p3mup6*F(2213)+p1mup6*
     -   p3mup3*p3mup5*F(2219)+p2mup6*p3mup3*p3mup5*F(2219)+p3mup3*p3
     -   mup5*p3mup6*F(9609)+(p1mup5+p2mup5)*p3mup3*p3mup6*F(9614))
       Fa20(476) = p5mup3*(p2mup6*p5mup5*F(2157)+4*(p2mup5*p3mup6*F(21
     -   60)+p2mup5*p5mup6*F(2171)-p3mup5*(p3mup6*F(2154)+p2mup6*F(21
     -   61)-p5mup6*F(2175)))-2*p3mup6*p5mup5*F(2176))+p3mup3*(4*p3mu
     -   p5*p5mup6*F(2153)+p5mup5*(p5mup6*F(2174)+2*(p2mup6*F(2177)+p
     -   3mup6*F(2479))))
       Fa20(477) = p6mup3*(-(p5mup5*p5mup6*F(2223))-4*(-((p1mup6+p2mup
     -   6)*p4mup5*F(9598))+p3mup5*(p4mup6*F(2214)-p1mup6*F(9604))+(p
     -   1mup5+p2mup5)*(p4mup6*F(2228)-(p1mup6+p2mup6)*F(9603)-p3mup6
     -   *F(9604))-p3mup6*p4mup5*F(9606)))
       Fa20(478) = -(p2mup3*p5mup5*p6mup6*F(2263))+p4mup5*(4*(p4mup6*p
     -   6mup3*F(1967)-p5mup6*p6mup3*F(2224))+p2mup3*p6mup6*F(2272))+
     -   2*(p2mup3*p3mup5*p6mup6*F(2242)+p4mup6*(p5mup5*p6mup3*F(2225
     -   )+2*(p4mup3*p6mup5*F(2243)+p2mup3*p6mup5*F(2246)))+p1mup5*p2
     -   mup3*p6mup6*F(2265)+p2mup3*p2mup5*p6mup6*F(2265)-p3mup6*p5mu
     -   p5*p6mup3*F(2277)-p1mup6*p5mup5*p6mup3*F(2480)-p2mup6*p5mup5
     -   *p6mup3*F(2480))+4*(p3mup5*p5mup6*p6mup3*F(2226)+p1mup5*p5mu
     -   p6*p6mup3*F(2229)+p2mup5*p5mup6*p6mup3*F(2229)+p6mup5*((p1mu
     -   p3+p2mup3)*(p1mup6+p2mup6)*F(2241)+p2mup3*p3mup6*F(2242))+p2
     -   mup3*p5mup6*p6mup5*F(2247)-p5mup3*p5mup6*p6mup5*F(2255)+p3mu
     -   p5*p3mup6*p6mup3*F(9599)+p3mup3*p3mup6*p6mup5*F(9599)+p2mup6
     -   *p3mup5*p6mup3*F(9604))
       Fa20(479) = p6mup5*(-(p4mup6*p5mup3*F(2244))-p1mup3*(p3mup6*F(2
     -   242)+p4mup6*F(2246)+p5mup6*F(2247))+(p1mup6+p2mup6)*p5mup3*F
     -   (2248)+p4mup3*p5mup6*F(2249)+p3mup6*p5mup3*F(2250)+p3mup6*p4
     -   mup3*F(2252)+(p1mup6+p2mup6)*p4mup3*F(2253)+p3mup3*(p4mup6*F
     -   (2214)-(p1mup6+p2mup6)*F(2251)+p5mup6*F(2254)))
       Fa20(480) = p6mup3*p6mup5*(-(p4mup6*F(2245))+p3mup6*F(9581)+(p1
     -   mup6+p2mup6)*F(9584)+p5mup6*F(9588))
       Return
       End
