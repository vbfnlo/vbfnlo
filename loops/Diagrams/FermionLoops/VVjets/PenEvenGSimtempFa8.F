c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenGSim8(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa8)
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
       Complex*16 Fa8(57:64)
       COMMON/FaPenEvenGSim/Fa
       Fa8(57) = mup3mup4*p1mup2*p3mup5*F(2201)-mup2mup5*p1mup4*p3mup3
     -   *F(2203)+mup2mup4*p1mup5*p3mup3*F(2203)+mup2mup5*p1mup3*p3mu
     -   p4*F(2203)-mup2mup3*p1mup5*p3mup4*F(2203)+mup2mup4*p1mup3*p3
     -   mup5*F(2203)-mup2mup3*p1mup4*p3mup5*F(2203)+mup3mup4*p1mup5*
     -   p3mup2*F(2204)+2*(mup2mup4*p1mup3*p1mup5*F(2200)-mup2mup3*p1
     -   mup4*p1mup5*F(2200)+mup2mup4*p1mup5*p2mup3*F(2200)-mup2mup3*
     -   p1mup5*p2mup4*F(2200)+mup2mup4*p1mup3*p2mup5*F(2200)-mup2mup
     -   3*p1mup4*p2mup5*F(2200)+mup3mup4*p1mup2*p1mup5*F(2206)+mup3m
     -   up4*p1mup5*p2mup2*F(2206)+mup3mup4*p1mup2*p2mup5*F(2206))+mu
     -   p3mup4*p1mup2*p4mup5*F(2210)-mup4mup5*(p1mup2*(2*(p1mup3+p2m
     -   up3)*F(2199)+p3mup3*F(2201)+p4mup3*F(2210))+p1mup3*(2*p2mup2
     -   *F(2199)+p3mup2*F(2202)+p4mup2*F(2211)))+mup3mup5*(p1mup2*(2
     -   *(p1mup4+p2mup4)*F(2199)+p3mup4*F(2201)+p4mup4*F(2210))+p1mu
     -   p4*(2*p2mup2*F(2199)+p3mup2*F(2202)+p4mup2*F(2211)))-mup2mup
     -   5*p1mup4*p4mup3*F(2212)+mup2mup4*p1mup5*p4mup3*F(2212)+mup2m
     -   up5*p1mup3*p4mup4*F(2212)-mup2mup3*p1mup5*p4mup4*F(2212)+mup
     -   2mup4*p1mup3*p4mup5*F(2212)-mup2mup3*p1mup4*p4mup5*F(2212)+m
     -   up3mup4*p1mup5*p4mup2*F(2213)
       Fa8(58) = p1mup3*p2mup2*p3mup4*p4mup5*F(2216)+p1mup3*p2mup2*p3m
     -   up5*p4mup4*F(2217)+p1mup4*p2mup2*p3mup3*p4mup5*F(2218)+p1mup
     -   4*p2mup2*p3mup5*p4mup3*F(2219)+p1mup3*p1mup4*p3mup2*p4mup5*F
     -   (2224)+p1mup4*p2mup3*p3mup2*p4mup5*F(2224)+p1mup3*p2mup4*p3m
     -   up2*p4mup5*F(2224)+p1mup3*p1mup4*p3mup5*p4mup2*F(2225)+p1mup
     -   4*p2mup3*p3mup5*p4mup2*F(2225)+p1mup3*p2mup4*p3mup5*p4mup2*F
     -   (2225)+p1mup3*p2mup5*p3mup2*p4mup4*F(2226)-p1mup3*p2mup5*p3m
     -   up4*p4mup2*F(2227)-p1mup4*p2mup5*p3mup2*p4mup3*F(2230)+p1mup
     -   4*p2mup5*p3mup3*p4mup2*F(2231)+4*p1mup4*p3mup2*p3mup3*p3mup5
     -   *F(2235)+2*(p1mup3*p2mup2*p3mup4*p3mup5*F(2214)+p1mup3*p2mup
     -   5*p3mup2*p3mup4*F(2222)+p1mup4*p2mup2*p3mup3*p3mup5*F(2222)+
     -   p1mup3*p1mup4*p3mup2*p3mup5*F(2223)+p1mup4*p2mup3*p3mup2*p3m
     -   up5*F(2223)+p1mup3*p2mup4*p3mup2*p3mup5*F(2223)+p1mup3*p3mup
     -   2*p3mup4*p3mup5*F(2229)+p1mup3*p2mup2*p4mup4*p4mup5*F(2232)+
     -   p1mup3*p1mup4*p4mup2*p4mup5*F(2233)+p1mup4*p2mup3*p4mup2*p4m
     -   up5*F(2233)+p1mup3*p2mup4*p4mup2*p4mup5*F(2233)+p1mup4*p2mup
     -   5*p3mup2*p3mup3*F(2234)+p1mup3*p2mup5*p4mup2*p4mup4*F(2237)+
     -   p1mup4*p2mup2*p4mup3*p4mup5*F(2237)+p1mup4*p2mup5*p4mup2*p4m
     -   up3*F(2238))+p1mup5*(4*p2mup2*p2mup4*p4mup3*F(2197)+p1mup2*p
     -   3mup3*p4mup4*F(2220)+p2mup2*p3mup3*p4mup4*F(2220)+p1mup2*p3m
     -   up4*p4mup3*F(2221)+p2mup2*p3mup4*p4mup3*F(2221)+p1mup3*p3mup
     -   2*p4mup4*F(2226)-p1mup3*p3mup4*p4mup2*F(2227)-p1mup4*p3mup2*
     -   p4mup3*F(2230)-p2mup4*p3mup2*p4mup3*F(2230)+p1mup4*p3mup3*p4
     -   mup2*F(2231)+p2mup4*p3mup3*p4mup2*F(2231)+p2mup3*(p3mup2*p4m
     -   up4*F(2226)-p3mup4*p4mup2*F(2227)+2*(p2mup4*p4mup2*F(2192)+p
     -   3mup2*p3mup4*F(2222)+p4mup2*p4mup4*F(2237)))+2*(p2mup2*p2mup
     -   3*p4mup4*F(2190)+p1mup2*p3mup3*p3mup4*F(2215)+p2mup2*p3mup3*
     -   p3mup4*F(2215)+p1mup3*p3mup2*p3mup4*F(2222)+p1mup4*p3mup2*p3
     -   mup3*F(2234)+p2mup4*p3mup2*p3mup3*F(2234)+p3mup2*p3mup3*p3mu
     -   p4*F(2236)+p1mup3*p4mup2*p4mup4*F(2237)+p1mup4*p4mup2*p4mup3
     -   *F(2238)+p2mup4*p4mup2*p4mup3*F(2238)+(p1mup2+p2mup2)*p4mup3
     -   *p4mup4*F(2239)))+p1mup2*(p1mup4*p3mup3*p4mup5*F(2218)+p2mup
     -   4*p3mup3*p4mup5*F(2218)+p1mup4*p3mup5*p4mup3*F(2219)+p2mup4*
     -   p3mup5*p4mup3*F(2219)+p2mup5*p3mup3*p4mup4*F(2220)+p2mup5*p3
     -   mup4*p4mup3*F(2221)+p1mup3*(p3mup4*p4mup5*F(2216)+p3mup5*p4m
     -   up4*F(2217)+2*(p3mup4*p3mup5*F(2214)+p4mup4*p4mup5*F(2232)))
     -   +p2mup3*(p3mup4*p4mup5*F(2216)+p3mup5*p4mup4*F(2217)+2*(p3mu
     -   p4*p3mup5*F(2214)+p4mup4*p4mup5*F(2232)))+2*(p2mup5*p3mup3*p
     -   3mup4*F(2215)+p1mup4*p3mup3*p3mup5*F(2222)+p2mup4*p3mup3*p3m
     -   up5*F(2222)+p3mup3*p3mup4*p3mup5*F(2228)+p1mup4*p4mup3*p4mup
     -   5*F(2237)+p2mup4*p4mup3*p4mup5*F(2237)+p2mup5*p4mup3*p4mup4*
     -   F(2239)))
       Fa8(59) = -2*p2mup1*p2mup3*(p1mup4*p1mup5*F(2251)+p1mup5*p2mup4
     -   *F(2253)+2*(p1mup4*p2mup5*F(2252)+p2mup4*p2mup5*F(2254)))+p1
     -   mup3*(-2*p2mup1*(p1mup5*p2mup4*F(2246)+p1mup4*p1mup5*F(2248)
     -   +p1mup4*p2mup5*F(2249)+p2mup4*p2mup5*F(2250))+p3mup1*(p3mup4
     -   *p4mup5*F(2257)+p3mup5*p4mup4*F(2258)))+p1mup1*(2*(p1mup3*(p
     -   1mup4*p2mup5*F(2241)+2*(p1mup4*p1mup5*F(2240)+p1mup5*p2mup4*
     -   F(2242))+p2mup4*p2mup5*F(2243))+p2mup3*(p1mup4*p1mup5*F(2244
     -   )+p1mup5*p2mup4*F(2245)-p1mup4*p2mup5*F(2246)-p2mup4*p2mup5*
     -   F(2247)))+p3mup3*p3mup4*p4mup5*F(2255)+p3mup3*p3mup5*p4mup4*
     -   F(2256)+p3mup4*p3mup5*p4mup3*F(2259))
       Fa8(60) = mup3mup4*p1mup1*p2mup5*F(2263)+mup1mup5*p1mup3*p2mup4
     -   *F(2267)+mup1mup4*p1mup3*p2mup5*F(2268)+mup1mup5*p1mup4*p2mu
     -   p3*F(2270)+mup1mup3*p1mup4*p2mup5*F(2271)+mup3mup4*p1mup5*p2
     -   mup1*F(2272)+mup1mup4*p1mup5*p2mup3*F(2273)+mup1mup3*p1mup5*
     -   p2mup4*F(2274)+mup4mup5*(p1mup1*p2mup3*F(2261)+p1mup3*p2mup1
     -   *F(2266)-2*(p1mup1*p1mup3*F(2260)-p2mup1*p2mup3*F(2277)+p3mu
     -   p1*p3mup3*F(2287)-p4mup1*p4mup3*F(2289)))+mup3mup5*(p1mup1*p
     -   2mup4*F(2262)+p1mup4*p2mup1*F(2269)-2*(p1mup1*p1mup4*F(2275)
     -   -p2mup1*p2mup4*F(2280)-p3mup1*p3mup4*F(2282)+p4mup1*p4mup4*F
     -   (2290)))-2*(mup1mup5*p1mup3*p1mup4*F(2264)-mup1mup4*p1mup3*p
     -   1mup5*F(2265)-mup1mup3*p1mup4*p1mup5*F(2276)-mup1mup5*p2mup3
     -   *p2mup4*F(2278)+mup1mup4*p2mup3*p2mup5*F(2279)+mup1mup3*p2mu
     -   p4*p2mup5*F(2281)-mup1mup5*p3mup3*p3mup4*F(2283)-mup1mup3*p3
     -   mup4*p3mup5*F(2284)-mup3mup4*p1mup1*p1mup5*F(2285)-mup3mup4*
     -   p2mup1*p2mup5*F(2286)-mup3mup4*p3mup1*p3mup5*F(2288)-mup3mup
     -   4*p4mup1*p4mup5*F(2291)-mup1mup5*p4mup3*p4mup4*F(2292)+mup1m
     -   up4*p4mup3*p4mup5*F(2293)-mup1mup4*p3mup3*p3mup5*F(2298))-4*
     -   (p1mup4*p3mup1*p3mup5*p4mup3*F(2294)-p1mup4*p3mup3*p3mup5*p4
     -   mup1*F(2295)-p1mup5*p3mup1*p3mup4*p4mup3*F(2296)+p1mup5*p3mu
     -   p3*p3mup4*p4mup1*F(2297)+p2mup1*p3mup3*p3mup4*p4mup5*F(2299)
     -   +p2mup1*p3mup3*p3mup5*p4mup4*F(2300)+p2mup3*p3mup1*p3mup4*p4
     -   mup5*F(2301)+p2mup3*p3mup1*p3mup5*p4mup4*F(2302)-p2mup4*p3mu
     -   p1*p3mup5*p4mup3*F(2303)+p2mup4*p3mup3*p3mup5*p4mup1*F(2304)
     -   +p2mup5*p3mup1*p3mup4*p4mup3*F(2305)-p2mup5*p3mup3*p3mup4*p4
     -   mup1*F(2306)-p1mup3*p3mup4*p3mup5*p4mup1*F(2307)-p1mup1*p3mu
     -   p3*p4mup4*p4mup5*F(2308)-p1mup1*p3mup4*p4mup3*p4mup5*F(2309)
     -   -p1mup3*p3mup1*p4mup4*p4mup5*F(2310)-p1mup3*p3mup4*p4mup1*p4
     -   mup5*F(2311)+p1mup4*p3mup1*p4mup3*p4mup5*F(2312)-p1mup4*p3mu
     -   p3*p4mup1*p4mup5*F(2313)+p2mup1*p3mup4*p3mup5*p4mup3*F(2314)
     -   +p2mup3*p3mup4*p3mup5*p4mup1*F(2315)+p2mup1*p3mup3*p4mup4*p4
     -   mup5*F(2316)+p2mup1*p3mup4*p4mup3*p4mup5*F(2317)+p2mup3*p3mu
     -   p1*p4mup4*p4mup5*F(2318)+p2mup3*p3mup4*p4mup1*p4mup5*F(2319)
     -   -p2mup4*p3mup1*p4mup3*p4mup5*F(2320)+p2mup4*p3mup3*p4mup1*p4
     -   mup5*F(2321)-p1mup4*p3mup1*p3mup3*p4mup5*F(2322)-p1mup5*p3mu
     -   p1*p3mup3*p4mup4*F(2323)-p2mup4*p3mup1*p3mup3*p4mup5*F(2324)
     -   -p2mup5*p3mup1*p3mup3*p4mup4*F(2325)-p1mup1*p3mup5*p4mup3*p4
     -   mup4*F(2328)-p1mup3*p3mup5*p4mup1*p4mup4*F(2329)-p1mup4*p3mu
     -   p5*p4mup1*p4mup3*F(2330)+p1mup5*p3mup1*p4mup3*p4mup4*F(2331)
     -   -p1mup5*p3mup3*p4mup1*p4mup4*F(2332)-p1mup5*p3mup4*p4mup1*p4
     -   mup3*F(2333)-p2mup4*p3mup5*p4mup1*p4mup3*F(2334)-p2mup5*p3mu
     -   p4*p4mup1*p4mup3*F(2335)+p2mup1*p3mup5*p4mup3*p4mup4*F(2345)
     -   +p2mup3*p3mup5*p4mup1*p4mup4*F(2346)-p2mup5*p3mup1*p4mup3*p4
     -   mup4*F(2347)+p2mup5*p3mup3*p4mup1*p4mup4*F(2348))+8*(p3mup1*
     -   p3mup3*p3mup4*p4mup5*F(2326)+p3mup1*p3mup3*p3mup5*p4mup4*F(2
     -   327)+p1mup1*p4mup3*p4mup4*p4mup5*F(2336)+p1mup4*p4mup1*p4mup
     -   3*p4mup5*F(2338)+p2mup3*p4mup1*p4mup4*p4mup5*F(2339)+p2mup4*
     -   p4mup1*p4mup3*p4mup5*F(2340)+p3mup3*p3mup4*p4mup1*p4mup5*F(2
     -   341)+p3mup4*p3mup5*p4mup1*p4mup3*F(2342)+p3mup3*p4mup1*p4mup
     -   4*p4mup5*F(2343)+p3mup4*p4mup1*p4mup3*p4mup5*F(2344)-p2mup1*
     -   p4mup3*p4mup4*p4mup5*F(2349)+p3mup1*p3mup4*p3mup5*p4mup3*F(2
     -   350)+p3mup1*p3mup3*p4mup4*p4mup5*F(2352)+p3mup3*p3mup5*p4mup
     -   1*p4mup4*F(2353)+p3mup1*p3mup4*p4mup3*p4mup5*F(2353)+p1mup5*
     -   p4mup1*p4mup3*p4mup4*F(2354)+p2mup5*p4mup1*p4mup3*p4mup4*F(2
     -   355)+p3mup5*p4mup1*p4mup3*p4mup4*F(2356)+p3mup1*p3mup5*p4mup
     -   3*p4mup4*F(2357))+16*(p1mup3*p4mup1*p4mup4*p4mup5*F(2337)+p3
     -   mup3*p3mup4*p3mup5*p4mup1*F(2351)+p3mup1*p4mup3*p4mup4*p4mup
     -   5*F(2358))
       Fa8(61) = mup4mup5*p1mup1*p2mup3*F(2360)+mup3mup5*p1mup1*p2mup4
     -   *F(2361)-mup3mup4*p1mup1*p2mup5*F(2361)+mup4mup5*p1mup3*p2mu
     -   p1*F(2364)-mup1mup4*(p1mup5*p2mup3+p1mup3*p2mup5)*F(2365)+mu
     -   p1mup5*(-2*(p2mup3*p2mup4*F(2200)-p1mup3*p1mup4*F(2362))+(p1
     -   mup4*p2mup3+p1mup3*p2mup4)*F(2365))-mup3mup5*p1mup4*p2mup1*F
     -   (2366)+mup3mup4*p1mup5*p2mup1*F(2366)-mup1mup3*p1mup5*p2mup4
     -   *F(2367)+mup1mup3*p1mup4*p2mup5*F(2367)+2*(mup4mup5*p1mup1*p
     -   1mup3*F(2359)+mup1mup4*(p2mup3*p2mup5*F(2200)+p1mup3*p1mup5*
     -   F(2363))-mup4mup5*p2mup1*p2mup3*F(2368)+mup3mup5*p1mup1*p1mu
     -   p4*F(2369)+mup3mup4*p1mup1*p1mup5*F(2370)+mup3mup5*p2mup1*p2
     -   mup4*F(2371)+mup3mup4*p2mup1*p2mup5*F(2372)+mup4mup5*p3mup1*
     -   p3mup3*F(2373)+mup3mup5*p3mup1*p3mup4*F(2373)+mup3mup4*p3mup
     -   1*p3mup5*F(2374)+mup4mup5*p4mup1*p4mup3*F(2375)+mup3mup5*p4m
     -   up1*p4mup4*F(2375)+mup3mup4*p4mup1*p4mup5*F(2376))
       Fa8(62) = mup4mup5*(p1mup3*p2mup1*F(2377)+p1mup1*p2mup3*F(2377)
     -   +2*(p1mup1*p1mup3*F(2369)+p4mup1*p4mup3*F(2375)-p2mup1*p2mup
     -   3*F(2379)+p3mup1*p3mup3*F(2380)))+mup3mup4*(p1mup5*p2mup1*F(
     -   2377)+p1mup1*p2mup5*F(2377)+2*(p1mup1*p1mup5*F(2369)+p4mup1*
     -   p4mup5*F(2375)-p2mup1*p2mup5*F(2379)+p3mup1*p3mup5*F(2380)))
     -   +mup3mup5*(p1mup4*p2mup1*F(2378)+p1mup1*p2mup4*F(2378)+2*(p1
     -   mup1*p1mup4*F(2370)+p4mup1*p4mup4*F(2376)+p2mup1*p2mup4*F(23
     -   79)+p3mup1*p3mup4*F(2381)))+(mup1mup5*p1mup4*p2mup3-mup1mup4
     -   *p1mup5*p2mup3-mup1mup5*p1mup3*p2mup4+mup1mup3*p1mup5*p2mup4
     -   +mup1mup4*p1mup3*p2mup5-mup1mup3*p1mup4*p2mup5)*F(3787)
       Fa8(63) = mup3mup5*p1mup1*p2mup4*F(2360)+mup3mup4*p1mup1*p2mup5
     -   *F(2361)+mup3mup5*p1mup4*p2mup1*F(2364)-mup4mup5*(p1mup1*p2m
     -   up3*F(2360)+p1mup3*p2mup1*F(2364))-mup1mup3*(p1mup5*p2mup4+p
     -   1mup4*p2mup5)*F(2365)+mup1mup4*(-2*(p2mup3*p2mup5*F(2200)-p1
     -   mup3*p1mup5*F(2362))+(p1mup5*p2mup3+p1mup3*p2mup5)*F(2365))-
     -   mup3mup4*p1mup5*p2mup1*F(2366)-mup1mup5*p1mup4*p2mup3*F(2367
     -   )+mup1mup5*p1mup3*p2mup4*F(2367)+2*(mup3mup5*p1mup1*p1mup4*F
     -   (2359)+mup1mup3*(p2mup4*p2mup5*F(2200)+p1mup4*p1mup5*F(2363)
     -   )+mup4mup5*p2mup1*p2mup3*F(2368)-mup3mup5*p2mup1*p2mup4*F(23
     -   68)+mup3mup4*p1mup1*p1mup5*F(2369)+mup3mup4*p2mup1*p2mup5*F(
     -   2371)+mup3mup5*p3mup1*p3mup4*F(2373)+mup3mup4*p3mup1*p3mup5*
     -   F(2373)+mup4mup5*p3mup1*p3mup3*F(2374)+mup4mup5*p1mup1*p1mup
     -   3*F(2382)+mup4mup5*p4mup1*p4mup3*F(2383)+p4mup1*(mup3mup5*p4
     -   mup4+mup3mup4*p4mup5)*F(2384))
       Fa8(64) = (mup1mup5*mup3mup4*F(2387))/2.+mup4mup5*p1mup1*p3mup3
     -   *F(2388)+mup3mup5*p1mup1*p3mup4*F(2389)+mup3mup4*p1mup1*p3mu
     -   p5*F(2390)+mup4mup5*p1mup3*p3mup1*F(2391)+mup1mup5*p1mup3*p3
     -   mup4*F(2392)+mup3mup5*p1mup4*p3mup1*F(2394)+mup1mup5*p1mup4*
     -   p3mup3*F(2395)+mup3mup4*p1mup5*p3mup1*F(2397)+mup4mup5*p1mup
     -   1*p4mup3*F(2400)+mup3mup5*p1mup1*p4mup4*F(2401)+mup3mup4*p1m
     -   up1*p4mup5*F(2402)+mup4mup5*p1mup3*p4mup1*F(2403)+mup1mup5*p
     -   1mup3*p4mup4*F(2404)+mup3mup5*p1mup4*p4mup1*F(2406)+mup1mup5
     -   *p1mup4*p4mup3*F(2407)+mup3mup4*p1mup5*p4mup1*F(2409)+mup4mu
     -   p5*p2mup1*p3mup3*F(2412)+mup3mup5*p2mup1*p3mup4*F(2413)+mup3
     -   mup4*p2mup1*p3mup5*F(2414)+mup4mup5*p2mup3*p3mup1*F(2415)+mu
     -   p1mup5*p2mup3*p3mup4*F(2416)+mup3mup5*p2mup4*p3mup1*F(2418)+
     -   mup1mup5*p2mup4*p3mup3*F(2419)+mup3mup4*p2mup5*p3mup1*F(2421
     -   )+mup4mup5*p2mup1*p4mup3*F(2424)+mup3mup5*p2mup1*p4mup4*F(24
     -   25)+mup3mup4*p2mup1*p4mup5*F(2426)+mup4mup5*p2mup3*p4mup1*F(
     -   2427)+mup1mup5*p2mup3*p4mup4*F(2428)+mup3mup5*p2mup4*p4mup1*
     -   F(2430)+mup1mup5*p2mup4*p4mup3*F(2431)+mup3mup4*p2mup5*p4mup
     -   1*F(2433)+mup4mup5*p3mup1*p4mup3*F(2436)+mup3mup5*p3mup1*p4m
     -   up4*F(2437)+mup3mup4*p3mup1*p4mup5*F(2438)+mup4mup5*p3mup3*p
     -   4mup1*F(2439)+mup1mup5*p3mup3*p4mup4*F(2440)+mup3mup5*p3mup4
     -   *p4mup1*F(2442)+mup1mup5*p3mup4*p4mup3*F(2443)+mup3mup4*p3mu
     -   p5*p4mup1*F(2445)+mup1mup4*((mup3mup5*F(2386))/2.+p1mup3*p3m
     -   up5*F(2393)+p1mup5*p3mup3*F(2398)+p1mup3*p4mup5*F(2405)+p1mu
     -   p5*p4mup3*F(2410)+p2mup3*p3mup5*F(2417)+p2mup5*p3mup3*F(2422
     -   )+p2mup3*p4mup5*F(2429)+p2mup5*p4mup3*F(2434)+p3mup3*p4mup5*
     -   F(2441)+p3mup5*p4mup3*F(2446))+8*(p1mup1*p1mup3*p1mup4*p3mup
     -   5*F(2448)+p1mup1*p1mup3*p2mup4*p3mup5*F(2450)+p1mup1*p1mup3*
     -   p2mup5*p3mup4*F(2451)+p1mup1*p1mup4*p1mup5*p3mup3*F(2452)+p1
     -   mup1*p1mup5*p2mup3*p3mup4*F(2453)+p1mup1*p1mup5*p2mup4*p3mup
     -   3*F(2454)+p1mup1*p1mup3*p1mup4*p4mup5*F(2455)+p1mup1*p1mup3*
     -   p1mup5*p4mup4*F(2456)-p1mup1*p1mup4*p2mup3*p3mup5*F(2457)-p1
     -   mup1*p1mup4*p2mup5*p3mup3*F(2458)-p1mup1*p2mup3*p2mup4*p3mup
     -   5*F(2459)-p1mup1*p2mup3*p2mup5*p3mup4*F(2460)-p1mup1*p2mup4*
     -   p2mup5*p3mup3*F(2461)-p1mup3*p1mup4*p2mup1*p3mup5*F(2462)-p1
     -   mup3*p1mup5*p2mup1*p3mup4*F(2463)-p1mup3*p2mup1*p2mup4*p3mup
     -   5*F(2464)-p1mup3*p2mup1*p2mup5*p3mup4*F(2465)-p1mup4*p1mup5*
     -   p2mup1*p3mup3*F(2466)-p1mup4*p2mup1*p2mup3*p3mup5*F(2467)-p1
     -   mup4*p2mup1*p2mup5*p3mup3*F(2468)-p1mup5*p2mup1*p2mup3*p3mup
     -   4*F(2469)-p1mup5*p2mup1*p2mup4*p3mup3*F(2470)-p2mup1*p2mup3*
     -   p2mup4*p3mup5*F(2471)-p2mup1*p2mup4*p2mup5*p3mup3*F(2473)+p1
     -   mup3*p1mup4*p1mup5*p3mup1*F(2474)+p1mup3*p1mup4*p2mup5*p3mup
     -   1*F(2475)+p1mup3*p1mup5*p2mup4*p3mup1*F(2476)+p1mup3*p2mup4*
     -   p2mup5*p3mup1*F(2477)-p1mup4*p1mup5*p2mup3*p3mup1*F(2478)-p1
     -   mup4*p2mup3*p2mup5*p3mup1*F(2479)-p1mup5*p2mup3*p2mup4*p3mup
     -   1*F(2480)-p2mup3*p2mup4*p2mup5*p3mup1*F(2481))+mup1mup3*((mu
     -   p4mup5*F(2385))/2.+p1mup4*p3mup5*F(2396)+p1mup5*p3mup4*F(239
     -   9)+p1mup4*p4mup5*F(2408)+p1mup5*p4mup4*F(2411)+p2mup4*p3mup5
     -   *F(2420)+p2mup5*p3mup4*F(2423)+p2mup4*p4mup5*F(2432)+p2mup5*
     -   p4mup4*F(2435)+p3mup4*p4mup5*F(2444)+p3mup5*p4mup4*F(2447)+2
     -   *p4mup4*p4mup5*F(2482))+16*(p1mup1*p1mup3*p1mup5*p3mup4*F(24
     -   49)-p2mup1*p2mup3*p2mup5*p3mup4*F(2472)+p3mup1*p3mup3*p3mup4
     -   *p3mup5*F(2483)+p4mup1*p4mup3*p4mup4*p4mup5*F(2484))
       Return
       End
