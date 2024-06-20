c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa21(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa21)
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
       Complex*16 Fa21(481:504)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa21(481) = p6mup6*(p2mup5*p5mup3*F(2264)+p2mup3*p6mup5*F(2282)
     -   )
       Fa21(482) = p6mup6*(-4*(p4mup3*p4mup5*F(2270)-p5mup3*p5mup5*F(2
     -   271)-p1mup3*(-(p5mup5*F(2263))+2*(p3mup5*F(2242)+(p1mup5+p2m
     -   up5)*F(2265))+p4mup5*F(2272))-p4mup3*p5mup5*F(2283)-p3mup3*p
     -   4mup5*F(2292)-p4mup5*p5mup3*F(2413))+p3mup3*p5mup5*F(2481)-8
     -   *(p1mup5*p5mup3*F(2264)+p1mup5*p4mup3*F(2278)+p2mup5*p4mup3*
     -   F(2278)-p3mup5*(-(p4mup3*F(2252))+p5mup3*F(2482)+p3mup3*F(95
     -   99))-(p1mup5+p2mup5)*p3mup3*F(9604)))
       Fa21(483) = p6mup3*p6mup6*(-4*p4mup5*F(2287)+p5mup5*F(2483)-8*(
     -   p3mup5*F(9581)+(p1mup5+p2mup5)*F(9584)))
       Fa21(484) = p6mup5*p6mup6*(-(p1mup3*F(2282))+p4mup3*F(2289)+p5m
     -   up3*F(2290)+p6mup3*F(9577)-p3mup3*F(9581))
       Fa21(485) = 2*(mup3mup5*F(2311)-p5mup3*p5mup5*F(2321))+(p1mup3+
     -   p2mup3)*p5mup5*F(2484)-p5mup5*p6mup3*F(2485)+p2mup5*p5mup3*F
     -   (2486)+p2mup5*p6mup3*F(2487)+(p1mup3+p2mup3)*p6mup5*F(2488)+
     -   p5mup3*p6mup5*F(2489)+p4mup3*p4mup5*F(2490)+4*((p1mup3+p2mup
     -   3)*(-((p1mup5+p2mup5)*F(2324))+p3mup5*F(2326))-p6mup3*p6mup5
     -   *F(2363)+p4mup5*p5mup3*F(2491)+p1mup3*p4mup5*F(2492)+p2mup3*
     -   p4mup5*F(2492)+p4mup5*p6mup3*F(2493))+p3mup5*p5mup3*F(2494)+
     -   p3mup5*p6mup3*F(2495)+p4mup3*p5mup5*F(2496)+p2mup5*p4mup3*F(
     -   2497)+p1mup5*(p5mup3*F(2486)+p6mup3*F(2487)+p4mup3*F(2497))+
     -   p4mup3*p6mup5*F(2498)+p3mup5*p4mup3*F(2499)+p3mup3*(-4*(p4mu
     -   p5*F(2306)+p3mup5*F(2344)+(p1mup5+p2mup5)*F(2345))+p5mup5*F(
     -   2500)+p6mup5*F(2501))
       Fa21(486) = 2*((p1mup3+p2mup3)*F(2294)-p6mup3*F(2299)+p3mup3*F(
     -   2301))+p4mup3*F(2502)+p5mup3*F(2503)
       Fa21(487) = -2*((p1mup3+p2mup3)*F(2298)+p6mup3*F(2300)-p3mup3*F
     -   (2302)+p4mup3*F(2416))+p5mup3*F(2504)
       Fa21(488) = mup5mup6*F(2505)+p4mup5*p4mup6*F(2508)+2*(p1mup6*p4
     -   mup5*F(2013)+p2mup6*p4mup5*F(2013)+p4mup5*p5mup6*F(2014)+p3m
     -   up6*p4mup5*F(2046)+p4mup6*p5mup5*F(2509))+4*(p3mup5*p6mup6*F
     -   (1878)+p6mup5*p6mup6*F(1999)-p4mup5*p6mup6*F(2011)-p5mup5*(p
     -   5mup6*F(1990)+(p1mup6+p2mup6)*F(2000)+p3mup6*F(2034))+p4mup6
     -   *p6mup5*F(2510))+p5mup5*p6mup6*F(2511)+8*(p1mup6*p3mup5*F(20
     -   16)+p2mup6*p3mup5*F(2016)+p3mup5*p5mup6*F(2017)+p3mup5*p3mup
     -   6*F(2044)+p1mup6*p6mup5*F(2506)+p2mup6*p6mup5*F(2506)+p3mup6
     -   *p6mup5*F(2507)+p5mup6*p6mup5*F(2512))+p1mup5*(4*p6mup6*F(18
     -   74)+8*(p5mup6*F(1992)+(p1mup6+p2mup6)*F(2003)+p3mup6*F(2016)
     -   )+p4mup6*F(2513))+p2mup5*(4*p6mup6*F(1874)+8*(p5mup6*F(1992)
     -   +(p1mup6+p2mup6)*F(2003)+p3mup6*F(2016))+p4mup6*F(2513))+p3m
     -   up5*p4mup6*F(2514)
       Fa21(489) = 2*(p4mup5*p5mup6*F(2021)+p4mup5*p6mup6*F(2022)+p5mu
     -   p5*p5mup6*F(2023)-p5mup5*p6mup6*F(2033))-4*(p1mup5*p6mup6*F(
     -   1875)+p2mup5*p6mup6*F(1875)+p3mup5*p6mup6*F(1879)-p6mup5*p6m
     -   up6*F(1911)+p4mup6*p6mup5*F(2040)-p1mup5*p4mup6*F(2041)-p2mu
     -   p5*p4mup6*F(2041)-p3mup6*p5mup5*F(2042)-p3mup5*p4mup6*F(2047
     -   ))-8*(p3mup6*p6mup5*F(1879)+p1mup5*p5mup6*F(2012)+p2mup5*p5m
     -   up6*F(2012)+p3mup5*p5mup6*F(2025)-p1mup5*p3mup6*F(2026)-p2mu
     -   p5*p3mup6*F(2026)+p5mup6*p6mup5*F(2043)-p3mup5*p3mup6*F(2050
     -   ))+mup5mup6*F(2515)+p4mup6*p5mup5*F(2516)+p4mup5*p4mup6*F(25
     -   17)+p1mup6*(4*p5mup5*F(2024)-8*(p6mup5*F(1875)-(p1mup5+p2mup
     -   5)*F(2015)-p3mup5*F(2026))+p4mup5*F(2518))+p2mup6*(4*p5mup5*
     -   F(2024)-8*(p6mup5*F(1875)-(p1mup5+p2mup5)*F(2015)-p3mup5*F(2
     -   026))+p4mup5*F(2518))+p3mup6*p4mup5*F(2519)
       Fa21(490) = mup3mup6*F(2520)+2*(p4mup3*p5mup6*F(1932)+p1mup6*p4
     -   mup3*F(1937)+p2mup6*p4mup3*F(1937)+p1mup6*p3mup3*F(1939)+p2m
     -   up6*p3mup3*F(1939)-4*(p3mup3*p6mup6*F(1935)-p4mup3*p6mup6*F(
     -   1947))+p3mup6*p4mup3*F(1954)+p5mup3*(-2*p6mup6*F(1929)+(p1mu
     -   p6+p2mup6)*F(1930)+p5mup6*F(1931)+p3mup6*F(1960)+p4mup6*F(19
     -   62))+p3mup3*p3mup6*F(1969)+p6mup3*(-4*p6mup6*F(1911)+(p1mup6
     -   +p2mup6)*F(1974)+p3mup6*F(1975)+p5mup6*F(1976)+p4mup6*F(1977
     -   ))+p3mup3*p5mup6*F(2522)+p4mup3*p4mup6*F(2523)+p1mup3*(4*p6m
     -   up6*F(1900)-p5mup6*F(1936)+p3mup6*F(1965)+(p1mup6+p2mup6)*F(
     -   2521)+p4mup6*F(2524))+p2mup3*(4*p6mup6*F(1900)-p5mup6*F(1936
     -   )+p3mup6*F(1965)+(p1mup6+p2mup6)*F(2521)+p4mup6*F(2524))+p3m
     -   up3*p4mup6*F(2525))
       Fa21(491) = p3mup3*p5mup6*F(2528)+2*(p6mup3*p6mup6*F(1911)+p3mu
     -   p3*p6mup6*F(1935)-p4mup3*p6mup6*F(1947)-p5mup3*(-(p5mup6*F(1
     -   957))+(p1mup6+p2mup6)*F(1959)+p3mup6*F(1961)+p4mup6*F(1963))
     -   -p4mup3*p4mup6*F(1967)+p5mup6*p6mup3*F(1973)+p5mup3*p6mup6*F
     -   (1978)+p3mup6*p4mup3*F(2214)+p1mup6*p4mup3*F(2228)+p2mup6*p4
     -   mup3*F(2228)+mup3mup6*F(2526)+p4mup3*p5mup6*F(2527)+p3mup6*p
     -   6mup3*F(9581)+p1mup6*p6mup3*F(9584)+p2mup6*p6mup3*F(9584)+p4
     -   mup6*p6mup3*F(9587)-p3mup3*p3mup6*F(9599)-p1mup6*p3mup3*F(96
     -   04)-p2mup6*p3mup3*F(9604)-p1mup3*(p6mup6*F(1900)+p5mup6*F(19
     -   58)+p4mup6*F(9598)+(p1mup6+p2mup6)*F(9603)+p3mup6*F(9604))-p
     -   2mup3*(p6mup6*F(1900)+p5mup6*F(1958)+p4mup6*F(9598)+(p1mup6+
     -   p2mup6)*F(9603)+p3mup6*F(9604))-p3mup3*p4mup6*F(9606))
       Fa21(492) = (p1mup6+p2mup6)*F(2529)+p5mup6*F(2530)+p6mup6*F(253
     -   1)+p3mup6*F(2532)+p4mup6*F(2533)
       Fa21(493) = p5mup2*(-(p5mup6*F(2534))+(p1mup6+p2mup6)*F(2535)+p
     -   6mup6*F(2536)+p3mup6*F(2540)+p4mup6*F(2542))-p6mup2*(p6mup6*
     -   F(2538)+2*(p5mup6*F(2537)-(p1mup6+p2mup6)*F(2539)+p3mup6*F(2
     -   541)+p4mup6*F(2543)))
       Fa21(494) = p1mup6*p2mup2*F(2544)+p2mup2*p2mup6*F(2544)+p2mup2*
     -   p5mup6*F(2545)+p2mup2*p6mup6*F(2546)+p2mup2*p3mup6*F(2547)+p
     -   2mup2*p4mup6*F(2548)+p1mup2*((p1mup6+p2mup6)*F(2544)+p5mup6*
     -   F(2545)+p6mup6*F(2546)+p3mup6*F(2547)+p4mup6*F(2548))+2*(p5m
     -   up6*p6mup2*s12*F(2549)+p6mup2*p6mup6*s12*F(2550)+p1mup6*p6mu
     -   p2*s12*F(2551)+p2mup6*p6mup2*s12*F(2551)+p3mup6*p6mup2*s12*F
     -   (2552)+p4mup6*p6mup2*s12*F(2553))
       Fa21(495) = p1mup2*(2*(p3mup6*F(2541)+p4mup6*F(2543)+(p1mup6+p2
     -   mup6)*F(2554))+p5mup6*F(2555)+p6mup6*F(2556))+p2mup2*(2*(p3m
     -   up6*F(2541)+p4mup6*F(2543)+(p1mup6+p2mup6)*F(2554))+p5mup6*F
     -   (2555)+p6mup6*F(2556))+2*p5mup2*s12*(-(p5mup6*F(2549))+p6mup
     -   6*F(2557)+(p1mup6+p2mup6)*F(2558)+p3mup6*F(2559)+p4mup6*F(25
     -   60))
       Fa21(496) = p1mup3*p6mup5*F(2562)+p2mup3*p6mup5*F(2562)+4*(-(mu
     -   p3mup5*F(2312))+p1mup5*p4mup3*F(2320)+p2mup5*p4mup3*F(2320)+
     -   p1mup5*p5mup3*F(2341)+p2mup5*p5mup3*F(2341)+p1mup5*p3mup3*F(
     -   2346)+p2mup5*p3mup3*F(2346)-p1mup3*p1mup5*F(2353)-p1mup5*p2m
     -   up3*F(2353)-p1mup3*p2mup5*F(2353)-p2mup3*p2mup5*F(2353)-p5mu
     -   p3*p6mup5*F(2359)-p4mup3*p6mup5*F(2360)-p3mup3*p6mup5*F(2362
     -   )-p6mup3*p6mup5*F(2377)+p1mup5*p6mup3*F(2561)+p2mup5*p6mup3*
     -   F(2561)+p3mup5*(p5mup3*F(2276)+p3mup3*F(2279)+p4mup3*F(2319)
     -   +(p1mup3+p2mup3)*F(2346)+p6mup3*F(2563))-p4mup5*(p5mup3*F(23
     -   07)+p4mup3*F(2308)+(p1mup3+p2mup3)*F(2309)-p3mup3*F(2310)-p6
     -   mup3*F(2566)))+2*(p5mup5*((p1mup3+p2mup3)*F(2330)+p6mup3*F(2
     -   364))+p5mup3*p5mup5*F(2564)+p4mup3*p5mup5*F(2565)+p3mup3*p5m
     -   up5*F(2567))
       Fa21(497) = 2*(p2mup3*p5mup5*F(2331)-p4mup3*p5mup5*F(2335)+p1mu
     -   p3*(p5mup5*F(2331)-2*(p4mup5*F(2343)+(p1mup5+p2mup5)*F(2349)
     -   -p3mup5*F(2350)-p6mup5*F(2365))))+p5mup3*(-2*p5mup5*F(2329)+
     -   4*((p1mup5+p2mup5)*F(2332)+p3mup5*F(2333)+p4mup5*F(2336)-p6m
     -   up5*F(2367)))+p5mup5*p6mup3*F(2568)+p3mup3*p5mup5*F(2569)+p3
     -   mup3*p6mup5*F(2570)+p3mup3*p4mup5*F(2571)-4*(mup3mup5*F(2313
     -   )+p4mup3*p4mup5*F(2339)-p1mup5*p4mup3*F(2342)-p2mup5*p4mup3*
     -   F(2342)+p2mup3*p4mup5*F(2343)+p1mup5*p2mup3*F(2349)+p2mup3*p
     -   2mup5*F(2349)-p2mup3*p3mup5*F(2350)-p3mup5*p4mup3*F(2351)-p2
     -   mup3*p6mup5*F(2365)+p4mup3*p6mup5*F(2368)-p1mup5*p6mup3*F(23
     -   71)-p2mup5*p6mup3*F(2371)-p3mup5*p6mup3*F(2373)-p4mup5*p6mup
     -   3*F(2375)+p6mup3*p6mup5*F(2376)-p3mup3*p3mup5*F(2572))+(p1mu
     -   p5+p2mup5)*p3mup3*F(2573)
       Fa21(498) = (p1mup5+p2mup5)*F(2574)+p5mup5*F(2575)+p6mup5*F(257
     -   6)+p3mup5*F(2577)+p4mup5*F(2578)
       Fa21(499) = p6mup2*(p4mup5*p6mup6*F(2588)-4*(mup5mup6*F(2581)-(
     -   p1mup6+p2mup6)*(p5mup5*F(2583)-4*(p3mup5*F(2584)+p4mup5*F(25
     -   86)))+p5mup5*p6mup6*F(2593))-8*(p2mup5*p6mup6*F(2582)+p4mup6
     -   *p5mup5*F(2590)-p5mup5*p5mup6*F(2596)-p1mup5*(-(p6mup6*F(258
     -   2))+2*(p4mup6*F(2589)+p5mup6*F(2592)+(p1mup6+p2mup6)*F(2597)
     -   ))+p6mup5*p6mup6*F(9579)-p3mup5*p6mup6*F(9588)-p3mup6*(p5mup
     -   5*F(2218)+2*((p1mup5+p2mup5)*F(2585)-p4mup5*F(2599)+p6mup5*F
     -   (9588)-p3mup5*F(9590))))+16*(p4mup5*p5mup6*F(2587)+p2mup5*p4
     -   mup6*F(2589)-p4mup6*p6mup5*F(2591)+p2mup5*p5mup6*F(2592)-p1m
     -   up6*p6mup5*F(2594)-p2mup6*p6mup5*F(2594)-p4mup5*p4mup6*F(259
     -   5)+p1mup6*p2mup5*F(2597)+p2mup5*p2mup6*F(2597)+p3mup5*p4mup6
     -   *F(2598)+p5mup6*p6mup5*F(9583)-p3mup5*p5mup6*F(9594)))
       Fa21(500) = p5mup2*(-16*p4mup5*p5mup6*F(2587)+mup5mup6*F(2600)+
     -   8*(p1mup5*(p6mup6*F(2582)-2*(p4mup6*F(2589)+p5mup6*F(2592)+(
     -   p1mup6+p2mup6)*F(2597)))-p3mup6*(p5mup5*F(2218)+2*((p1mup5+p
     -   2mup5)*F(2585)-p4mup5*F(2599)+p6mup5*F(9588)-p3mup5*F(9590))
     -   ))-4*((p1mup6+p2mup6)*(p5mup5*F(2583)-4*(p3mup5*F(2584)+p4mu
     -   p5*F(2586)))-p5mup5*p6mup6*F(2593)-p4mup5*p6mup6*F(2601)+2*(
     -   p5mup5*p5mup6*F(2596)-p2mup5*(p6mup6*F(2582)-2*(p4mup6*F(258
     -   9)+p5mup6*F(2592)+(p1mup6+p2mup6)*F(2597)))-p4mup6*(p5mup5*F
     -   (2590)+2*(p6mup5*F(2591)+p4mup5*F(2595)-p3mup5*F(2598)))-p6m
     -   up5*p6mup6*F(9579)+p3mup5*p6mup6*F(9588))-4*((p1mup6+p2mup6)
     -   *p6mup5*F(2594)-p5mup6*p6mup5*F(9583)+p3mup5*p5mup6*F(9594))
     -   ))
       Fa21(501) = p6mup2*(p6mup3*p6mup6*F(2608)-16*(-(p5mup6*p6mup3*F
     -   (2606))+p3mup3*p5mup6*F(2610)+p1mup3*p5mup6*F(2613)+p2mup3*p
     -   5mup6*F(2613))+4*(-(mup3mup6*F(2602))+p4mup3*(8*(p4mup6*F(19
     -   67)-p3mup6*F(2214)-(p1mup6+p2mup6)*F(2228))+4*p5mup6*F(2603)
     -   +p6mup6*F(2604))+p3mup3*p6mup6*F(2612)+p5mup3*(p6mup6*F(2607
     -   )+4*(p5mup6*F(2605)+p4mup6*F(2609)+p3mup6*F(2611)+(p1mup6+p2
     -   mup6)*F(2614)))+p1mup3*p6mup6*F(2615)+p2mup3*p6mup6*F(2615))
     -   -32*(p3mup6*p6mup3*F(9581)+p1mup6*p6mup3*F(9584)+p2mup6*p6mu
     -   p3*F(9584)+p4mup6*p6mup3*F(9587)-p1mup3*p4mup6*F(9598)-p2mup
     -   3*p4mup6*F(9598)-p3mup3*p3mup6*F(9599)-p1mup3*p1mup6*F(9603)
     -   -p1mup6*p2mup3*F(9603)-p1mup3*p2mup6*F(9603)-p2mup3*p2mup6*F
     -   (9603)-p1mup6*p3mup3*F(9604)-p2mup6*p3mup3*F(9604)-p1mup3*p3
     -   mup6*F(9604)-p2mup3*p3mup6*F(9604)-p3mup3*p4mup6*F(9606)))
       Fa21(502) = p5mup2*(16*(-(p5mup6*p6mup3*F(2606))+p3mup3*p5mup6*
     -   F(2610)+p1mup3*p5mup6*F(2613)+p2mup3*p5mup6*F(2613))+p3mup3*
     -   p6mup6*F(2617)+4*(mup3mup6*F(2602)-p4mup3*(8*(p4mup6*F(1967)
     -   -p3mup6*F(2214)-(p1mup6+p2mup6)*F(2228))+4*p5mup6*F(2603)+p6
     -   mup6*F(2604))+p6mup3*p6mup6*F(2616)+p5mup3*(-(p6mup6*F(2607)
     -   )-4*(p5mup6*F(2605)+p4mup6*F(2609)-p3mup6*F(2618)-(p1mup6+p2
     -   mup6)*F(2619))))+p1mup3*p6mup6*F(2620)+p2mup3*p6mup6*F(2620)
     -   +32*(p3mup6*p6mup3*F(9581)+p1mup6*p6mup3*F(9584)+p2mup6*p6mu
     -   p3*F(9584)+p4mup6*p6mup3*F(9587)-p1mup3*p4mup6*F(9598)-p2mup
     -   3*p4mup6*F(9598)-p3mup3*p3mup6*F(9599)-p1mup3*p1mup6*F(9603)
     -   -p1mup6*p2mup3*F(9603)-p1mup3*p2mup6*F(9603)-p2mup3*p2mup6*F
     -   (9603)-p1mup6*p3mup3*F(9604)-p2mup6*p3mup3*F(9604)-p1mup3*p3
     -   mup6*F(9604)-p2mup3*p3mup6*F(9604)-p3mup3*p4mup6*F(9606)))
       Fa21(503) = p1mup5*p5mup2*F(2621)+p2mup5*p5mup2*F(2621)+p5mup2*
     -   p5mup5*F(2622)+p5mup2*p6mup5*F(2623)+p5mup5*p6mup2*F(2625)+p
     -   6mup2*p6mup5*F(2626)+p3mup5*p5mup2*F(2627)+2*(p1mup5*p6mup2*
     -   F(2624)+p2mup5*p6mup2*F(2624)-p3mup5*p6mup2*F(2628)-p4mup5*p
     -   6mup2*F(2629))+p4mup5*p5mup2*F(2630)
       Fa21(504) = p6mup2*(32*p3mup5*(p5mup3*F(2276)+p3mup3*F(2279)+p4
     -   mup3*F(2319)+(p1mup3+p2mup3)*F(2346))-16*p4mup5*(2*(p5mup3*F
     -   (2307)+p4mup3*F(2308)+(p1mup3+p2mup3)*F(2309)-p3mup3*F(2310)
     -   )+p6mup3*F(2637))+p3mup5*p6mup3*F(2641)+4*(8*(p1mup5*(p4mup3
     -   *F(2320)+p5mup3*F(2341)+p3mup3*F(2346)-(p1mup3+p2mup3)*F(235
     -   3))+p2mup5*(p4mup3*F(2320)+p5mup3*F(2341)+p3mup3*F(2346)-(p1
     -   mup3+p2mup3)*F(2353)))-p5mup3*p5mup5*F(2631)-p4mup3*p5mup5*F
     -   (2632)-p5mup5*p6mup3*F(2635)+p3mup3*p5mup5*F(2639)+p1mup3*p5
     -   mup5*F(2642)+p2mup3*p5mup5*F(2642)-4*(p4mup3*p6mup5*F(2633)+
     -   p5mup3*p6mup5*F(2634)+p6mup3*p6mup5*F(2636)+mup3mup5*F(2638)
     -   -p3mup3*p6mup5*F(2640)-(p1mup3+p2mup3)*p6mup5*F(2643)))+(p1m
     -   up5+p2mup5)*p6mup3*F(2644))
       Return
       End
