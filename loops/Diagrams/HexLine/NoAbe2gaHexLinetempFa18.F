c************** Calling the Fa functions*************************
       subroutine NoAbe2gaHexLineFa18(p1mup2,p1mup3,p1mup5,p1mup6,p2mu
     -   p2,p2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p
     -   4mup3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup
     -   3,p6mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6
     -   ,mup5mup6,Fa18)
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
       Complex*16 Fa18(409:432)
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       Fa18(409) = mup5mup6*p6mup2*((p1mup3+p2mup3)*F(2380)+p4mup3*F(2
     -   381)+p5mup3*F(2382)+p3mup3*F(2383))
       Fa18(410) = p2mup2*(mup3mup6*(p6mup5*F(2386)+4*(p5mup5*F(2385)+
     -   (p1mup5+p2mup5)*F(2388)))+mup3mup5*((p1mup6+p2mup6)*F(2384)+
     -   p6mup6*F(2387)+p5mup6*F(2389)))
       Fa18(411) = 4*mup5mup6*(p1mup2+p2mup2)*p6mup3*F(2070)+mup3mup5*
     -   p1mup2*((p1mup6+p2mup6)*F(2384)+p6mup6*F(2387)+p5mup6*F(2389
     -   ))+16*mup3mup5*p5mup6*p6mup2*s12*F(2391)+mup3mup6*(p1mup2*(4
     -   *p5mup5*F(2385)+p6mup5*F(2386))+4*(p1mup2*p1mup5*F(2388)+p1m
     -   up2*p2mup5*F(2388)+p5mup5*p6mup2*s12*F(2390)-4*(p1mup5*p6mup
     -   2*s12*F(1874)+p2mup5*p6mup2*s12*F(1874)+p6mup2*p6mup5*s12*F(
     -   9579))))
       Fa18(412) = 2*mup3mup5*p6mup2*s12*(2*(p1mup6+p2mup6)*F(1875)-p6
     -   mup6*F(1911))+mup3mup6*(p1mup2+p2mup2)*(p4mup5*F(2392)+p3mup
     -   5*F(2393))
       Fa18(413) = mup3mup6*p6mup2*(p3mup5*F(1878)+p4mup5*F(2394))
       Fa18(414) = mup3mup5*((p1mup2+p2mup2)*p3mup6*F(2395)+(p1mup2+p2
     -   mup2)*p4mup6*F(2396)+16*(p3mup6*p6mup2*s12*F(1879)+p4mup6*p6
     -   mup2*s12*F(2397)))
       Fa18(415) = p2mup2*(p4mup3*(p4mup6*p5mup5*F(2089)+2*(p4mup5*p5m
     -   up6*F(2088)+p4mup5*p4mup6*F(2097))+p5mup5*p5mup6*F(2099))+p5
     -   mup3*(p4mup6*p5mup5*F(2091)+2*(p4mup5*p5mup6*F(2090)+p4mup5*
     -   p4mup6*F(2098))+p5mup5*p5mup6*F(2100))-2*p1mup3*p1mup5*p1mup
     -   6*F(2117))-2*(p1mup2*(p1mup5*p1mup6*p2mup3+p1mup3*p1mup6*p2m
     -   up5+p1mup3*p1mup5*p2mup6)*F(2117)+p1mup2*p1mup3*p1mup5*p3mup
     -   6*F(2126))
       Fa18(416) = p1mup3*(4*p1mup2*(p2mup5*p2mup6*F(2117)+p2mup5*p3mu
     -   p6*F(2126)+p2mup6*p3mup5*F(2127))+p1mup6*(-(p1mup2*p5mup5*F(
     -   2125))+4*(p2mup2*p2mup5*F(2117)+p1mup2*p3mup5*F(2127)+p1mup2
     -   *p4mup5*F(2128))))+4*(p2mup3*p2mup5*(p1mup2*p1mup6+p2mup2*p2
     -   mup6)*F(2117)+p1mup5*(p1mup6*p2mup2*p2mup3*F(2117)+p1mup3*p2
     -   mup2*p2mup6*F(2117)+p1mup2*p2mup3*p2mup6*F(2117)+p1mup2*p2mu
     -   p3*p3mup6*F(2126)+p1mup2*p1mup3*p4mup6*F(2129)+p1mup2*p1mup3
     -   *p5mup6*F(2130)+p1mup2*p1mup6*p4mup3*F(2131)+p1mup2*p1mup6*p
     -   5mup3*F(2132)-p1mup2*(p1mup6+p2mup6)*p3mup3*F(2133)))
       Fa18(417) = p2mup2*(4*(p1mup3*(p2mup5*p2mup6*F(2117)+p2mup5*p3m
     -   up6*F(2126)+p2mup6*p3mup5*F(2127))+p1mup5*(p2mup3*p2mup6*F(2
     -   117)+p1mup3*p3mup6*F(2126)+p2mup3*p3mup6*F(2126)+p1mup3*p4mu
     -   p6*F(2129)+p1mup3*p5mup6*F(2130)-p2mup6*p3mup3*F(2133)))+p1m
     -   up6*(p1mup3*(-(p5mup5*F(2125))+4*(p3mup5*F(2127)+p4mup5*F(21
     -   28)))+4*(p2mup3*p2mup5*F(2117)+p1mup5*(p4mup3*F(2131)+p5mup3
     -   *F(2132)-p3mup3*F(2133)))))+p1mup2*(p1mup3*(p4mup6*(-2*p5mup
     -   5*F(2085)+4*(p4mup5*F(2096)+p2mup5*F(2129)+p3mup5*F(2135)))-
     -   p5mup5*(p2mup6*F(2125)+2*(p5mup6*F(2086)+p3mup6*F(2143)))+4*
     -   (p4mup5*p5mup6*F(2087)+p2mup6*p4mup5*F(2128)+p2mup5*p5mup6*F
     -   (2130)+p3mup5*p3mup6*F(2134)+p3mup6*p4mup5*F(2136)+p3mup5*p5
     -   mup6*F(2145)))+4*(p2mup3*p2mup5*p2mup6*F(2117)+p2mup3*p2mup5
     -   *p3mup6*F(2126)+p2mup3*p2mup6*p3mup5*F(2127)-p2mup5*p2mup6*p
     -   3mup3*F(2133)+p1mup5*(p2mup3*p4mup6*F(2129)+p2mup3*p5mup6*F(
     -   2130)-p3mup3*p4mup6*F(2137)+p4mup3*(-(p4mup6*F(2109))+p5mup6
     -   *F(2111)+p2mup6*F(2131)-p3mup6*F(2138))+p3mup3*p3mup6*F(2141
     -   )+p3mup3*p5mup6*F(2146)+p5mup3*(p5mup6*F(2112)-p4mup6*F(2113
     -   )+p2mup6*F(2132)+p3mup6*F(2147))))+p1mup6*(p5mup3*p5mup5*F(2
     -   104)-p2mup3*p5mup5*F(2125)+p4mup3*(p5mup5*F(2102)+4*(p4mup5*
     -   F(2101)+p2mup5*F(2131)+p3mup5*F(2140)))-p3mup3*p5mup5*F(2144
     -   )+4*(p2mup3*p3mup5*F(2127)+p2mup5*p5mup3*F(2132)-p2mup5*p3mu
     -   p3*F(2133)+p4mup5*(p5mup3*F(2103)+p2mup3*F(2128)+p3mup3*F(21
     -   39))-p3mup3*p3mup5*F(2142)+p3mup5*p5mup3*F(2148))))
       Fa18(418) = p1mup2*(p2mup3*(-(p5mup5*(p2mup6*F(2125)+2*p3mup6*F
     -   (2143)))+4*(p2mup5*p4mup6*F(2129)+p2mup5*p5mup6*F(2130)+p3mu
     -   p5*p3mup6*F(2134)+p3mup5*p4mup6*F(2135)+p4mup5*(p5mup6*F(208
     -   7)+p4mup6*F(2096)+p2mup6*F(2128)+p3mup6*F(2136))+p3mup5*p5mu
     -   p6*F(2145)))+4*p2mup5*(-(p3mup3*p4mup6*F(2137))-p4mup3*(p4mu
     -   p6*F(2109)+p3mup6*F(2138))+p3mup3*p3mup6*F(2141)+p3mup3*p5mu
     -   p6*F(2146)+p3mup6*p5mup3*F(2147))+p2mup6*(-(p3mup3*p5mup5*F(
     -   2144))+4*(p2mup5*p5mup3*F(2132)+p3mup3*p4mup5*F(2139)+p4mup3
     -   *(p4mup5*F(2101)+p2mup5*F(2131)+p3mup5*F(2140))-p3mup3*p3mup
     -   5*F(2142)+p3mup5*p5mup3*F(2148))))+p2mup2*(p1mup3*(-(p5mup5*
     -   (p2mup6*F(2125)+2*p3mup6*F(2143)))+4*(p2mup5*p4mup6*F(2129)+
     -   p2mup5*p5mup6*F(2130)+p3mup5*p3mup6*F(2134)+p3mup5*p4mup6*F(
     -   2135)+p4mup5*(p4mup6*F(2096)+p2mup6*F(2128)+p3mup6*F(2136))+
     -   p3mup5*p5mup6*F(2145)))+4*p1mup5*(p2mup3*p4mup6*F(2129)+p2mu
     -   p3*p5mup6*F(2130)+p2mup6*p5mup3*F(2132)-p3mup3*p4mup6*F(2137
     -   )-p4mup3*(p4mup6*F(2109)-p2mup6*F(2131)+p3mup6*F(2138))+p3mu
     -   p3*p3mup6*F(2141)+p3mup3*p5mup6*F(2146)+p3mup6*p5mup3*F(2147
     -   ))+p1mup6*(-(p2mup3*p5mup5*F(2125))-p3mup3*p5mup5*F(2144)+4*
     -   (p2mup3*p3mup5*F(2127)+p2mup3*p4mup5*F(2128)+p2mup5*p5mup3*F
     -   (2132)-p2mup5*p3mup3*F(2133)+p3mup3*p4mup5*F(2139)+p4mup3*(p
     -   4mup5*F(2101)+p2mup5*F(2131)+p3mup5*F(2140))-p3mup3*p3mup5*F
     -   (2142)+p3mup5*p5mup3*F(2148))))
       Fa18(419) = 2*p1mup3*p2mup2*(p4mup6*p5mup5*F(2085)+p5mup5*p5mup
     -   6*F(2086)-2*p4mup5*p5mup6*F(2087))-p2mup2*(p1mup6*(p4mup3*p5
     -   mup5*F(2102)+4*p4mup5*p5mup3*F(2103)+p5mup3*p5mup5*F(2104))+
     -   4*p1mup5*(p4mup3*p5mup6*F(2111)+p5mup3*p5mup6*F(2112)-p4mup6
     -   *p5mup3*F(2113)))+p1mup2*(-(p2mup6*p4mup3*p5mup5*F(2102))-p2
     -   mup6*p5mup3*p5mup5*F(2104)+2*(p2mup3*p5mup5*(p4mup6*F(2085)+
     -   p5mup6*F(2086))+p3mup3*p5mup5*p5mup6*F(2107)-p3mup3*p4mup6*p
     -   5mup5*F(2115)-p3mup6*(p4mup3*p5mup5*F(2093)+p3mup3*p5mup5*F(
     -   2119)+2*(p4mup3*p4mup5*F(2092)+p4mup5*p5mup3*F(2094)+p3mup5*
     -   p4mup3*F(2118)+p3mup5*p5mup3*F(2120)-p3mup3*p4mup5*F(2121)-p
     -   3mup3*p3mup5*F(2122))))-4*(p2mup6*p4mup5*p5mup3*F(2103)-p3mu
     -   p3*p4mup5*p5mup6*F(2105)-p3mup5*p4mup3*p5mup6*F(2106)+p3mup3
     -   *p4mup5*p4mup6*F(2110)+p2mup5*p4mup3*p5mup6*F(2111)+p2mup5*p
     -   5mup3*p5mup6*F(2112)-p2mup5*p4mup6*p5mup3*F(2113)+p3mup5*p4m
     -   up3*p4mup6*F(2114)+p3mup5*p4mup6*p5mup3*F(2116)+p3mup3*p3mup
     -   5*p5mup6*F(2123)-p3mup3*p3mup5*p4mup6*F(2124)))
       Fa18(420) = 2*p1mup2*(p4mup3*(p4mup6*p5mup5*F(2089)+2*(p4mup5*p
     -   5mup6*F(2088)+p4mup5*p4mup6*F(2097))+p5mup5*p5mup6*F(2099))+
     -   p5mup3*(p4mup6*p5mup5*F(2091)-p3mup6*p5mup5*F(2095)+p5mup5*p
     -   5mup6*F(2100)+2*(p4mup5*(p5mup6*F(2090)+p4mup6*F(2098))+p3mu
     -   p5*p5mup6*F(2108)))-2*p1mup3*p1mup5*p1mup6*F(2117))+p2mup2*(
     -   p2mup3*(p2mup6*p5mup5*F(2125)+2*(p4mup6*(p5mup5*F(2085)-2*(p
     -   4mup5*F(2096)+p2mup5*F(2129)+p3mup5*F(2135)))+p3mup6*p5mup5*
     -   F(2143))-4*(p2mup6*p3mup5*F(2127)+p2mup5*(p3mup6*F(2126)+p5m
     -   up6*F(2130))+p3mup5*p3mup6*F(2134)+p4mup5*(p5mup6*F(2087)+p2
     -   mup6*F(2128)+p3mup6*F(2136))+p3mup5*p5mup6*F(2145)))+4*p2mup
     -   5*(p4mup6*p5mup3*F(2113)+p3mup3*p4mup6*F(2137)+p4mup3*(p4mup
     -   6*F(2109)-p5mup6*F(2111)+p3mup6*F(2138))-p3mup3*p3mup6*F(214
     -   1)-p3mup3*p5mup6*F(2146)-p3mup6*p5mup3*F(2147))-p2mup6*(p4mu
     -   p3*(p5mup5*F(2102)+4*(p4mup5*F(2101)+p2mup5*F(2131)+p3mup5*F
     -   (2140)))-p3mup3*p5mup5*F(2144)+4*(p2mup5*p5mup3*F(2132)-p2mu
     -   p5*p3mup3*F(2133)+p4mup5*(p5mup3*F(2103)+p3mup3*F(2139))-p3m
     -   up3*p3mup5*F(2142)+p3mup5*p5mup3*F(2148))))
       Fa18(421) = p2mup2*(-(p2mup6*p5mup3*p5mup5*F(2104))+2*(p2mup3*p
     -   5mup5*p5mup6*F(2086)+p3mup3*p5mup5*p5mup6*F(2107)-p3mup3*p4m
     -   up6*p5mup5*F(2115)-p3mup6*(p4mup3*p5mup5*F(2093)+p5mup3*p5mu
     -   p5*F(2095)+p3mup3*p5mup5*F(2119)+2*(p4mup3*p4mup5*F(2092)+p4
     -   mup5*p5mup3*F(2094)+p3mup5*p4mup3*F(2118)+p3mup5*p5mup3*F(21
     -   20)-p3mup3*p4mup5*F(2121)-p3mup3*p3mup5*F(2122))))+4*(p3mup3
     -   *p4mup5*p5mup6*F(2105)+p3mup5*p4mup3*p5mup6*F(2106)+p3mup5*p
     -   5mup3*p5mup6*F(2108)-p3mup3*p4mup5*p4mup6*F(2110)-p2mup5*p5m
     -   up3*p5mup6*F(2112)-p3mup5*p4mup3*p4mup6*F(2114)-p3mup5*p4mup
     -   6*p5mup3*F(2116)-p3mup3*p3mup5*p5mup6*F(2123)+p3mup3*p3mup5*
     -   p4mup6*F(2124)))
       Fa18(422) = (p1mup2+p2mup2)*p6mup3*(p2mup6*p5mup5*F(2206)+p1mup
     -   6*(p5mup5*F(2206)-4*(p3mup5*F(2208)-(p1mup5+p2mup5)*F(2216))
     -   )+2*(p3mup6*p5mup5*F(2215)+p4mup6*p5mup5*F(2217)+p5mup5*p5mu
     -   p6*F(2218))-4*(p2mup6*p3mup5*F(2208)-(p1mup5+p2mup5)*(p3mup6
     -   *F(2207)+p4mup6*F(2209)+p5mup6*F(2211))+p4mup5*(p4mup6*F(196
     -   7)+(p1mup6+p2mup6)*F(2210)+p5mup6*F(2212)-p3mup6*F(2214))-p1
     -   mup5*p2mup6*F(2216)-p2mup5*p2mup6*F(2216)+p3mup5*p5mup6*F(95
     -   90)+p3mup5*p3mup6*F(9599)+p3mup5*p4mup6*F(9606)))
       Fa18(423) = p6mup5*(p1mup3*(p1mup6*p2mup2*F(2231)+(p1mup2+p2mup
     -   2)*(p2mup6*F(2231)+p3mup6*F(2232)+p4mup6*F(2233)+p5mup6*F(22
     -   34)))+p2mup2*(p2mup3*(p2mup6*F(2231)+p3mup6*F(2232)+p4mup6*F
     -   (2233)+p5mup6*F(2234))+p2mup6*p3mup3*F(2235)-p2mup6*p4mup3*F
     -   (2236)-p2mup6*p5mup3*F(2237)+p1mup6*(p2mup3*F(2231)+p3mup3*F
     -   (2235)-p4mup3*F(2236)-p5mup3*F(2237))-p3mup6*p4mup3*F(2240)-
     -   p3mup3*p3mup6*F(9590)-p3mup6*p5mup3*F(9594)-p3mup3*p5mup6*F(
     -   9594)-p3mup3*p4mup6*F(9601))+p1mup2*(p2mup3*(p2mup6*F(2231)+
     -   p3mup6*F(2232)+p4mup6*F(2233)+p5mup6*F(2234))+p2mup6*p3mup3*
     -   F(2235)-p2mup6*p4mup3*F(2236)-p2mup6*p5mup3*F(2237)+p1mup6*(
     -   p2mup3*F(2231)+p3mup3*F(2235)-p4mup3*F(2236)-p5mup3*F(2237))
     -   +p4mup3*p4mup6*F(2238)-p4mup3*p5mup6*F(2239)-p3mup6*p4mup3*F
     -   (2240)-p3mup3*p3mup6*F(9590)-p3mup6*p5mup3*F(9594)-p3mup3*p5
     -   mup6*F(9594)-p4mup6*p5mup3*F(9596)-p3mup3*p4mup6*F(9601)-p5m
     -   up3*p5mup6*F(9602)))
       Fa18(424) = p6mup5*(-(p1mup2*p1mup3*p1mup6*F(2231))+p2mup2*(-(p
     -   4mup3*p4mup6*F(2238))+p4mup3*p5mup6*F(2239)+p4mup6*p5mup3*F(
     -   9596)+p5mup3*p5mup6*F(9602)))
       Fa18(425) = (p1mup2+p2mup2)*(p1mup6+p2mup6)*p6mup3*p6mup5
       Fa18(426) = p6mup6*(p2mup2*(p1mup3*(-(p5mup5*F(2256))+2*((p1mup
     -   5+p2mup5)*F(2257)+p3mup5*F(2258)+p4mup5*F(2259)))+p2mup3*(-(
     -   p5mup5*F(2256))+2*((p1mup5+p2mup5)*F(2257)+p3mup5*F(2258)+p4
     -   mup5*F(2259)))+2*(p1mup5+p2mup5)*(p3mup3*F(2260)-p4mup3*F(22
     -   61)-p5mup3*F(2262)))+p1mup2*(p1mup3*(-(p5mup5*F(2256))+2*(p2
     -   mup5*F(2257)+p3mup5*F(2258)+p4mup5*F(2259)))+p2mup3*(-(p5mup
     -   5*F(2256))+2*((p1mup5+p2mup5)*F(2257)+p3mup5*F(2258)+p4mup5*
     -   F(2259)))+p4mup3*p5mup5*F(2267)+p5mup3*p5mup5*F(2269)+2*(p1m
     -   up5*p3mup3*F(2260)+p2mup5*p3mup3*F(2260)-p1mup5*p4mup3*F(226
     -   1)-p2mup5*p4mup3*F(2261)-p1mup5*p5mup3*F(2262)-p2mup5*p5mup3
     -   *F(2262)+p4mup3*p4mup5*F(2266)+p4mup5*p5mup3*F(2268)-p3mup3*
     -   p4mup5*F(2273)-p3mup5*p4mup3*F(2274))-p3mup3*p5mup5*F(2275)-
     -   6*(p3mup5*p5mup3*F(2276)+p3mup3*p3mup5*F(2279))))
       Fa18(427) = p6mup6*(2*p1mup2*p1mup3*p1mup5*F(2257)+p2mup2*(p5mu
     -   p3*p5mup5*F(2269)+2*(p4mup5*p5mup3*F(2268)-p3mup3*p4mup5*F(2
     -   273))+p4mup3*(p5mup5*F(2267)+2*(p4mup5*F(2266)-p3mup5*F(2274
     -   )))-p3mup3*p5mup5*F(2275)-6*(p3mup5*p5mup3*F(2276)+p3mup3*p3
     -   mup5*F(2279))))
       Fa18(428) = (p1mup2+p2mup2)*p6mup6*(p6mup3*((p1mup5+p2mup5)*F(2
     -   280)+p4mup5*F(2285)+p5mup5*F(2286)-3*p3mup5*F(9581))+p6mup5*
     -   (-((p1mup3+p2mup3)*F(2281))+p4mup3*F(2291)+3*(p5mup3*F(9583)
     -   +p3mup3*F(9588))))
       Fa18(429) = mup5mup6*(-2*p5mup2*s12*(p5mup3*F(2382)+p3mup3*F(23
     -   83))+p2mup2*(p2mup3*F(2063)+p4mup3*F(2398)))
       Fa18(430) = mup5mup6*(-(p2mup2*(p5mup3*F(2060)+p3mup3*F(2061)-p
     -   1mup3*F(2063)))+p1mup2*(-(p5mup3*F(2060))-p3mup3*F(2061)+(p1
     -   mup3+p2mup3)*F(2063)+p4mup3*F(2398))+2*p5mup2*s12*(-(p4mup3*
     -   F(2381))+(p1mup3+p2mup3)*F(2399)))
       Fa18(431) = -4*mup3mup5*(4*(p1mup6*p5mup2*s12*F(1875)+p2mup6*p5
     -   mup2*s12*F(1875)+p5mup2*p5mup6*s12*F(2391))+p2mup2*p6mup6*F(
     -   2404)+2*(-(p5mup2*p6mup6*s12*F(1911))+p1mup6*p2mup2*F(2403)+
     -   p2mup2*p2mup6*F(2403)-p2mup2*p5mup6*F(2405)))+mup3mup6*(-4*p
     -   5mup2*p5mup5*s12*F(2390)+p1mup5*p2mup2*F(2400)+p2mup2*p2mup5
     -   *F(2400)+p2mup2*p5mup5*F(2401)+p2mup2*p6mup5*F(2402)+16*(p1m
     -   up5*p5mup2*s12*F(1874)+p2mup5*p5mup2*s12*F(1874)+p5mup2*p6mu
     -   p5*s12*F(9579)))
       Fa18(432) = 4*mup5mup6*p6mup3*((p1mup2+p2mup2)*F(2222)-2*p5mup2
     -   *s12*F(2376))+p1mup2*(mup3mup6*((p1mup5+p2mup5)*F(2400)+p5mu
     -   p5*F(2401)+p6mup5*F(2402))-4*mup3mup5*(p6mup6*F(2404)+2*((p1
     -   mup6+p2mup6)*F(2403)-p5mup6*F(2405))))
       Return
       End
