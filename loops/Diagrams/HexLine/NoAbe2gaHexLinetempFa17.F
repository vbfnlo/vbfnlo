c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa17(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa17)
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
       Complex*16 Fa17(385:408)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa17(385) = p6mup2*(-(p1mup6*p5mup3*p5mup5*F(2157))-p2mup6*p5mu
     -   p3*p5mup5*F(2157)+p4mup3*p5mup5*p5mup6*F(2164)-p1mup3*p5mup5
     -   *p5mup6*F(2170)-p3mup3*p5mup5*p5mup6*F(2174)-2*(p1mup3*p1mup
     -   6*p5mup5*F(2158)+p1mup3*p2mup6*p5mup5*F(2158)+p1mup3*p3mup6*
     -   p5mup5*F(2159)+p4mup3*p4mup6*p5mup5*F(2163)+p4mup6*p5mup3*p5
     -   mup5*F(2166)+p1mup3*p4mup6*p5mup5*F(2167)-p3mup6*p5mup3*p5mu
     -   p5*F(2176)+p1mup6*p3mup3*p5mup5*F(2177)+p2mup6*p3mup3*p5mup5
     -   *F(2177)-p3mup3*p4mup6*p5mup5*F(2179)-p3mup6*p4mup3*p5mup5*F
     -   (2182)-p1mup6*p4mup3*p5mup5*F(2187)-p2mup6*p4mup3*p5mup5*F(2
     -   187))+p2mup3*(-2*(p1mup6*p5mup5*F(2158)+p3mup6*p5mup5*F(2159
     -   )+p4mup6*p5mup5*F(2167))-p5mup5*p5mup6*F(2170)+p2mup6*(-2*p5
     -   mup5*F(2158)+4*(p3mup5*F(2150)+p1mup5*F(2188)+p4mup5*F(2192)
     -   ))+4*(p4mup5*p5mup6*F(2172)+p4mup5*p4mup6*F(2183)+p1mup5*p1m
     -   up6*F(2188)+p1mup5*p4mup6*F(2189)+p1mup5*p5mup6*F(2190)+p2mu
     -   p5*(p3mup6*F(2152)+p4mup6*F(2189)+p5mup6*F(2190))+p1mup6*p4m
     -   up5*F(2192)+p3mup5*p4mup6*F(2193)+p3mup6*p4mup5*F(2194)+p3mu
     -   p5*p5mup6*F(2195)))+4*(-(p1mup5*p3mup6*p5mup3*F(2160))+p1mup
     -   6*p3mup5*p5mup3*F(2161)+p2mup6*p3mup5*p5mup3*F(2161)+p4mup3*
     -   p4mup5*p5mup6*F(2162)+p4mup5*p5mup3*p5mup6*F(2165)-p1mup5*p4
     -   mup6*p5mup3*F(2168)+p1mup6*p4mup5*p5mup3*F(2169)+p2mup6*p4mu
     -   p5*p5mup3*F(2169)-p1mup5*p5mup3*p5mup6*F(2171)+p1mup3*p4mup5
     -   *p5mup6*F(2172)+p3mup5*p4mup3*p4mup6*F(2173)-p3mup5*p5mup3*p
     -   5mup6*F(2175)-p3mup3*p4mup5*p5mup6*F(2178)-p3mup5*p4mup3*p5m
     -   up6*F(2180)+p3mup5*p4mup6*p5mup3*F(2181)+p1mup3*p4mup5*p4mup
     -   6*F(2183)-p4mup5*p4mup6*p5mup3*F(2184)+p1mup5*p4mup3*p4mup6*
     -   F(2185)-p1mup5*p4mup3*p5mup6*F(2186)+p1mup3*p1mup5*p1mup6*F(
     -   2188)+p1mup3*p1mup5*p2mup6*F(2188)+p1mup3*p1mup5*p4mup6*F(21
     -   89)+p1mup3*p1mup5*p5mup6*F(2190)-p1mup5*p1mup6*p5mup3*F(2191
     -   )-p1mup5*p2mup6*p5mup3*F(2191)+p1mup3*p1mup6*p4mup5*F(2192)+
     -   p1mup3*p2mup6*p4mup5*F(2192)+p1mup3*p3mup5*p4mup6*F(2193)+p1
     -   mup3*p3mup6*p4mup5*F(2194)+p1mup3*p3mup5*p5mup6*F(2195)-p1mu
     -   p6*p4mup3*p4mup5*F(2196)-p2mup6*p4mup3*p4mup5*F(2196)+p3mup3
     -   *p4mup5*p4mup6*F(2197)+p1mup5*p1mup6*p4mup3*F(2198)+p1mup5*p
     -   2mup6*p4mup3*F(2198)+p1mup5*p3mup3*p4mup6*F(2199)-p1mup6*p3m
     -   up3*p4mup5*F(2200)-p2mup6*p3mup3*p4mup5*F(2200)-p1mup6*p3mup
     -   5*p4mup3*F(2201)-p2mup6*p3mup5*p4mup3*F(2201)-p1mup5*p3mup3*
     -   p5mup6*F(2202)+p3mup6*p4mup5*p5mup3*F(2203)-p1mup5*p3mup6*p4
     -   mup3*F(2204)-p2mup5*(p3mup6*p5mup3*F(2160)+p4mup6*p5mup3*F(2
     -   168)+p5mup3*p5mup6*F(2171)-p4mup3*p4mup6*F(2185)+p4mup3*p5mu
     -   p6*F(2186)-p1mup3*p1mup6*F(2188)-p1mup3*(p4mup6*F(2189)+p5mu
     -   p6*F(2190))+p1mup6*p5mup3*F(2191)-p1mup6*p4mup3*F(2198)+p2mu
     -   p6*(p3mup3*F(2156)-p1mup3*F(2188)+p5mup3*F(2191)-p4mup3*F(21
     -   98))-p3mup3*p4mup6*F(2199)+p3mup3*p5mup6*F(2202)+p3mup6*p4mu
     -   p3*F(2204))-p3mup6*p4mup3*p4mup5*F(2205)+p3mup3*p3mup5*p3mup
     -   6*F(9609)))
       Fa17(386) = p2mup6*p5mup2*p5mup5*p6mup3*F(2206)+p1mup6*p5mup2*p
     -   6mup3*(p5mup5*F(2206)-4*(p3mup5*F(2208)-(p1mup5+p2mup5)*F(22
     -   16)))+2*(p3mup6*p5mup2*p5mup5*p6mup3*F(2215)+p4mup6*p5mup2*p
     -   5mup5*p6mup3*F(2217)+p5mup2*p5mup5*p5mup6*p6mup3*F(2218))+4*
     -   (p2mup3*p2mup5*(p1mup6+p2mup6)*p6mup2*F(2188)+p1mup5*p3mup6*
     -   p5mup2*p6mup3*F(2207)+p2mup5*p3mup6*p5mup2*p6mup3*F(2207)-p2
     -   mup6*p3mup5*p5mup2*p6mup3*F(2208)+p1mup5*p4mup6*p5mup2*p6mup
     -   3*F(2209)+p2mup5*p4mup6*p5mup2*p6mup3*F(2209)+p1mup5*p5mup2*
     -   p5mup6*p6mup3*F(2211)+p2mup5*p5mup2*p5mup6*p6mup3*F(2211)+p1
     -   mup3*p3mup5*p3mup6*p6mup2*F(2213)+p2mup3*p3mup5*p3mup6*p6mup
     -   2*F(2213)-p4mup5*p5mup2*p6mup3*(p4mup6*F(1967)+(p1mup6+p2mup
     -   6)*F(2210)+p5mup6*F(2212)-p3mup6*F(2214))+p1mup5*p2mup6*p5mu
     -   p2*p6mup3*F(2216)+p2mup5*p2mup6*p5mup2*p6mup3*F(2216)+p1mup6
     -   *p3mup3*p3mup5*p6mup2*F(2219)+p2mup6*p3mup3*p3mup5*p6mup2*F(
     -   2219)-p3mup3*p3mup5*p4mup6*p6mup2*F(2220)-p3mup5*p3mup6*p4mu
     -   p3*p6mup2*F(2221)-p3mup5*p5mup2*p5mup6*p6mup3*F(9590)-p3mup5
     -   *p3mup6*p5mup2*p6mup3*F(9599)-p3mup5*p4mup6*p5mup2*p6mup3*F(
     -   9606)+(p1mup5+p2mup5)*p3mup3*p3mup6*p6mup2*F(9614)+p3mup3*p3
     -   mup6*p4mup5*p6mup2*F(9615))
       Fa17(387) = p6mup2*p6mup3*(4*mup5mup6*F(2222)-8*p5mup5*p5mup6*F
     -   (2223)+16*p4mup6*p5mup5*F(2225)-p3mup6*p5mup5*F(2227)-p1mup6
     -   *p5mup5*F(2230)-p2mup6*p5mup5*F(2230)+32*(-(p1mup5*p4mup6*F(
     -   2228))-p2mup5*p4mup6*F(2228)+p1mup5*p5mup6*F(2229)+p2mup5*p5
     -   mup6*F(2229)+p1mup5*p1mup6*F(9603)+p1mup6*p2mup5*F(9603)+p1m
     -   up5*p2mup6*F(9603)+p2mup5*p2mup6*F(9603)+(p1mup5+p2mup5)*p3m
     -   up6*F(9604)+p3mup5*(-(p4mup6*F(2214))+p5mup6*F(2226)+p3mup6*
     -   F(9599)+(p1mup6+p2mup6)*F(9604))+p4mup5*(p4mup6*F(1967)-p5mu
     -   p6*F(2224)+(p1mup6+p2mup6)*F(9598)+p3mup6*F(9606))))
       Fa17(388) = p5mup2*p6mup5*(p2mup3*(p2mup6*F(2231)+p3mup6*F(2232
     -   ))+p1mup3*((p1mup6+p2mup6)*F(2231)+p3mup6*F(2232)+p4mup6*F(2
     -   233)+p5mup6*F(2234))+p1mup6*(p2mup3*F(2231)+p3mup3*F(2235)-p
     -   4mup3*F(2236)-p5mup3*F(2237)))
       Fa17(389) = p5mup2*p6mup5*(-(p2mup3*(p4mup6*F(2233)+p5mup6*F(22
     -   34)))+p2mup6*(-(p3mup3*F(2235))+p4mup3*F(2236)+p5mup3*F(2237
     -   ))-p4mup3*p4mup6*F(2238)+p4mup3*p5mup6*F(2239)+p3mup6*p4mup3
     -   *F(2240)+p3mup3*p3mup6*F(9590)+p3mup6*p5mup3*F(9594)+p3mup3*
     -   p5mup6*F(9594)+p4mup6*p5mup3*F(9596)+p3mup3*p4mup6*F(9601)+p
     -   5mup3*p5mup6*F(9602))
       Fa17(390) = p6mup5*(p1mup6*p2mup3*p6mup2*F(2241)+p2mup3*p2mup6*
     -   p6mup2*F(2241)+p2mup3*p3mup6*p6mup2*F(2242)+p4mup3*p4mup6*p6
     -   mup2*F(2243)+p4mup6*p5mup3*p6mup2*F(2244)+p4mup6*p6mup2*p6mu
     -   p3*F(2245)+p2mup3*p4mup6*p6mup2*F(2246)+p2mup3*p5mup6*p6mup2
     -   *F(2247)+p1mup3*p6mup2*((p1mup6+p2mup6)*F(2241)+p3mup6*F(224
     -   2)+p4mup6*F(2246)+p5mup6*F(2247))-p1mup6*p5mup3*p6mup2*F(224
     -   8)-p2mup6*p5mup3*p6mup2*F(2248)-p4mup3*p5mup6*p6mup2*F(2249)
     -   -p3mup6*p5mup3*p6mup2*F(2250)-p3mup6*p4mup3*p6mup2*F(2252)-p
     -   1mup6*p4mup3*p6mup2*F(2253)-p2mup6*p4mup3*p6mup2*F(2253)-p5m
     -   up3*p5mup6*p6mup2*F(2255)+p1mup6*p5mup2*p6mup3*F(9579)+p2mup
     -   6*p5mup2*p6mup3*F(9579)-p3mup6*p6mup2*p6mup3*F(9581)-p1mup6*
     -   p6mup2*p6mup3*F(9584)-p2mup6*p6mup2*p6mup3*F(9584)-p5mup6*p6
     -   mup2*p6mup3*F(9588)+p3mup3*p6mup2*(-(p4mup6*F(2214))+(p1mup6
     -   +p2mup6)*F(2251)-p5mup6*F(2254)+p3mup6*F(9599)))
       Fa17(391) = p5mup2*p6mup6*(p1mup3*(-(p5mup5*F(2256))+2*((p1mup5
     -   +p2mup5)*F(2257)+p3mup5*F(2258)+p4mup5*F(2259)))+2*p1mup5*(p
     -   2mup3*F(2257)+p3mup3*F(2260)-p4mup3*F(2261)-p5mup3*F(2262)))
       Fa17(392) = p6mup6*(p4mup3*p5mup2*p5mup5*F(2267)+p5mup2*p5mup3*
     -   p5mup5*F(2269)-p4mup3*p4mup5*p6mup2*F(2270)+p5mup3*p5mup5*p6
     -   mup2*F(2271)+p1mup3*p6mup2*(-(p5mup5*F(2263))+2*(p3mup5*F(22
     -   42)+(p1mup5+p2mup5)*F(2265))+p4mup5*F(2272))+p2mup3*(2*p3mup
     -   5*(p6mup2*F(2242)+p5mup2*F(2258))+p5mup2*(-(p5mup5*F(2256))+
     -   2*(p2mup5*F(2257)+p4mup5*F(2259)))+p6mup2*(-(p5mup5*F(2263))
     -   +2*(p1mup5+p2mup5)*F(2265)+p4mup5*F(2272)))-p3mup3*p5mup2*p5
     -   mup5*F(2275)-p3mup3*p5mup5*p6mup2*F(2277)-2*(p3mup5*p5mup3*p
     -   6mup2*F(2250)-p2mup5*p3mup3*p5mup2*F(2260)+p2mup5*p4mup3*p5m
     -   up2*F(2261)+p2mup5*p5mup2*p5mup3*F(2262)+p1mup5*p5mup3*p6mup
     -   2*F(2264)+p2mup5*p5mup3*p6mup2*F(2264)-p4mup3*p4mup5*p5mup2*
     -   F(2266)-p4mup5*p5mup2*p5mup3*F(2268)+p3mup3*p4mup5*p5mup2*F(
     -   2273)+p3mup5*p4mup3*p5mup2*F(2274)+p1mup5*p4mup3*p6mup2*F(22
     -   78)+p2mup5*p4mup3*p6mup2*F(2278))-6*(p3mup5*p5mup2*p5mup3*F(
     -   2276)+p3mup3*p3mup5*p5mup2*F(2279)))
       Fa17(393) = p6mup6*(p4mup3*p5mup5*p6mup2*F(2283)-p4mup5*p5mup3*
     -   p6mup2*F(2284)-p4mup5*p6mup2*p6mup3*F(2287)+p5mup5*p6mup2*p6
     -   mup3*F(2288)+p3mup3*p4mup5*p6mup2*F(2292)+6*(p5mup2*p5mup3*p
     -   6mup5*F(9583)+p3mup3*p5mup2*p6mup5*F(9588))-2*(-((p1mup5+p2m
     -   up5)*p5mup2*p6mup3*F(2280))+(p1mup3+p2mup3)*p5mup2*p6mup5*F(
     -   2281)+p1mup3*p6mup2*p6mup5*F(2282)+p2mup3*p6mup2*p6mup5*F(22
     -   82)-p4mup5*p5mup2*p6mup3*F(2285)-p5mup2*p5mup5*p6mup3*F(2286
     -   )-p4mup3*p6mup2*p6mup5*F(2289)-p5mup3*p6mup2*p6mup5*F(2290)-
     -   p4mup3*p5mup2*p6mup5*F(2291)-p6mup2*p6mup3*p6mup5*F(9577)+p3
     -   mup3*p6mup2*p6mup5*F(9581)+p1mup5*p6mup2*p6mup3*F(9584)+p2mu
     -   p5*p6mup2*p6mup3*F(9584)+p3mup5*(p4mup3*p6mup2*F(2252)+(3*p5
     -   mup2+p6mup2)*p6mup3*F(9581)-p3mup3*p6mup2*F(9599))-(p1mup5+p
     -   2mup5)*p3mup3*p6mup2*F(9604)))
       Fa17(394) = p5mup2*(p5mup3*F(2293)+p4mup3*F(2295)-4*((p1mup3+p2
     -   mup3)*F(2294)-p6mup3*F(2299)+p3mup3*F(2301)))+4*p6mup2*(p4mu
     -   p3*F(2296)+p5mup3*F(2297)+(p1mup3+p2mup3)*F(2298)+p6mup3*F(2
     -   300)-p3mup3*F(2302))
       Fa17(395) = p4mup5*(p2mup2*(-4*p4mup3*F(2303)+p5mup3*F(2304)-8*
     -   (p2mup3*F(2305)-p3mup3*F(2306)))+8*p5mup2*s12*(p5mup3*F(2307
     -   )+p4mup3*F(2308)+p2mup3*F(2309)-p3mup3*F(2310)))
       Fa17(396) = p4mup5*(p1mup2*(-4*p4mup3*F(2303)+p5mup3*F(2304)-8*
     -   ((p1mup3+p2mup3)*F(2305)-p3mup3*F(2306)))-8*p1mup3*(p2mup2*F
     -   (2305)-p5mup2*s12*F(2309)))
       Fa17(397) = -4*mup3mup5*((p1mup2+p2mup2)*F(2311)-2*s12*(p5mup2*
     -   F(2312)+p6mup2*F(2313)))
       Fa17(398) = (p1mup2+p2mup2)*p4mup5*p6mup3
       Fa17(399) = p4mup3*(p1mup2*(p1mup5+p2mup5)*F(2316)+p2mup2*(p5mu
     -   p5*F(2315)+(p1mup5+p2mup5)*F(2316)+p3mup5*F(2317))+p5mup2*s1
     -   2*(p5mup5*F(2318)-2*(p3mup5*F(2319)+p2mup5*F(2320))))
       Fa17(400) = (p1mup2+p2mup2)*(4*(p5mup3*p5mup5*F(2321)+p1mup3*p5
     -   mup5*F(2322)+p2mup3*p5mup5*F(2322))+p5mup5*p6mup3*F(2323)+8*
     -   (p2mup3*p2mup5*F(2324)+p2mup5*p5mup3*F(2325)))+4*(p1mup2*p4m
     -   up3*p5mup5*F(2315)+p1mup2*p3mup5*p4mup3*F(2317)+p5mup3*p5mup
     -   5*p6mup2*s12*F(2329)-p2mup3*p5mup2*p5mup5*s12*F(2330)-p2mup3
     -   *p5mup5*p6mup2*s12*F(2331)+p5mup2*p5mup3*p5mup5*s12*F(2334)+
     -   p4mup3*p5mup5*p6mup2*s12*F(2335)-p1mup2*p3mup3*p5mup5*F(2337
     -   )-p2mup2*p3mup3*p5mup5*F(2337)+p3mup3*p5mup5*p6mup2*s12*F(23
     -   38))+p3mup3*p5mup2*p5mup5*s12*F(2340)-8*(p1mup5*p4mup3*p5mup
     -   2*s12*F(2320)+(p1mup2+p2mup2)*(p1mup3+p2mup3)*p3mup5*F(2326)
     -   +p2mup5*p5mup3*p6mup2*s12*F(2332)+p4mup5*p5mup3*p6mup2*s12*F
     -   (2336)-p4mup3*p4mup5*p6mup2*s12*F(2339)+p2mup5*p5mup2*p5mup3
     -   *s12*F(2341)+p2mup5*p4mup3*p6mup2*s12*F(2342)-p2mup3*p4mup5*
     -   p6mup2*s12*F(2343)-p2mup2*p2mup5*p3mup3*F(2345)+p2mup5*p3mup
     -   3*p5mup2*s12*F(2346)+p3mup5*(p1mup2*p5mup3*F(2327)+p2mup2*p5
     -   mup3*F(2327)+p1mup2*p6mup3*F(2328)+p2mup2*p6mup3*F(2328)+p5m
     -   up3*p6mup2*s12*F(2333)-(p1mup2+p2mup2)*p3mup3*F(2344)+p5mup2
     -   *s12*(p5mup3*F(2276)+p3mup3*F(2279)+p2mup3*F(2346))))+p3mup3
     -   *p3mup5*p6mup2*s12*F(2347)
       Fa17(401) = 4*p1mup3*(-(p5mup2*s12*(p5mup5*F(2330)+2*p3mup5*F(2
     -   346)))+2*(p1mup5+p2mup5)*(p2mup2*F(2324)+p5mup2*s12*F(2353))
     -   )+8*(p2mup2*p2mup5*p6mup3*F(2348)+p1mup2*(p1mup5*p2mup3*F(23
     -   24)+p1mup3*(p1mup5+p2mup5)*F(2324)+p1mup5*p5mup3*F(2325)+p1m
     -   up5*p3mup3*F(2345)+p2mup5*p3mup3*F(2345)+(p1mup5+p2mup5)*p6m
     -   up3*F(2348))+p1mup5*(p2mup2*(p2mup3*F(2324)+p5mup3*F(2325)+p
     -   3mup3*F(2345)+p6mup3*F(2348))-p5mup2*s12*(p5mup3*F(2341)+p3m
     -   up3*F(2346)-p2mup3*F(2353)))+s12*(p2mup3*p2mup5*p6mup2*F(234
     -   9)-p2mup3*p3mup5*p6mup2*F(2350)-p3mup5*p4mup3*p6mup2*F(2351)
     -   +p3mup3*p4mup5*p6mup2*F(2352)+p2mup3*p2mup5*p5mup2*F(2353)+p
     -   2mup5*p3mup3*p6mup2*F(2354)))
       Fa17(402) = p6mup5*(p2mup2*p4mup3*F(2357)+p1mup2*(p4mup3*F(2357
     -   )+2*((p1mup3+p2mup3)*F(2355)+p5mup3*F(2356)-p3mup3*F(2361)))
     -   +2*(p2mup2*p2mup3*F(2355)+p2mup2*p5mup3*F(2356)+p2mup3*p5mup
     -   2*s12*F(2358)+p1mup3*(p2mup2*F(2355)+p5mup2*s12*F(2358))+p5m
     -   up2*p5mup3*s12*F(2359)+p4mup3*p5mup2*s12*F(2360)-p2mup2*p3mu
     -   p3*F(2361)+p3mup3*p5mup2*s12*F(2362)))
       Fa17(403) = 8*(p1mup2*p6mup3*p6mup5*F(2363)+p2mup2*p6mup3*p6mup
     -   5*F(2363))-4*s12*(p5mup2*p5mup5*p6mup3*F(2364)+2*(p2mup3*p6m
     -   up2*p6mup5*F(2365)-p2mup5*p5mup2*p6mup3*F(2366)))
       Fa17(404) = 4*p6mup2*(p1mup3*(-(p5mup5*F(2331))+2*(p4mup5*F(234
     -   3)+(p1mup5+p2mup5)*F(2349)-p3mup5*F(2350)))+2*(-(p1mup5*p5mu
     -   p3*F(2332))+p1mup5*p3mup3*F(2354)+p5mup3*p6mup5*F(2367)+p4mu
     -   p3*p6mup5*F(2368)+p3mup3*p6mup5*F(2369)))
       Fa17(405) = 4*p6mup3*(2*(p1mup5*p5mup2*F(2366)+p5mup2*(p3mup5*F
     -   (2372)+p4mup5*F(2374)))+p6mup2*(p5mup5*F(2370)-2*(p2mup5*F(2
     -   371)+p3mup5*F(2373)+p4mup5*F(2375))))
       Fa17(406) = p1mup3*p6mup2*p6mup5*F(2365)+p6mup3*(p1mup5*p6mup2*
     -   F(2371)-p6mup5*(p6mup2*F(2376)+p5mup2*F(2377)))
       Fa17(407) = mup5mup6*p2mup2*(4*p2mup3*F(2058)+p4mup3*F(2378))
       Fa17(408) = mup5mup6*(p2mup2*(4*(p3mup3*F(2056)+p1mup3*F(2058))
     -   -p5mup3*F(2379))+p1mup2*(4*(p3mup3*F(2056)+(p1mup3+p2mup3)*F
     -   (2058))+p4mup3*F(2378)-p5mup3*F(2379)))
       Return
       End
