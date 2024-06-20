c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenG6(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup
     -   1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3
     -   mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3
     -   ,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,
     -   mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa6)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenEvenG/F
      Complex*16 Fa(150),F(4300)
      Real*8 P(2563) 
        COMMON/PPenEvenG/P
       Complex*16 Fa6(51:60)
       COMMON/FaPenEvenG/Fa
       Fa6(51) = mup2mup3*mup4mup5*F(2130)+mup3mup5*(2*p2mup4*p4mup2*F
     -   (2030)+mup2mup4*F(2131)+4*(p3mup4*p4mup2*F(2136)+p3mup2*p4mu
     -   p4*F(2138)))+2*((mup2mup5*p2mup4*p4mup3+mup2mup5*p2mup3*p4mu
     -   p4+mup2mup3*p2mup4*p4mup5)*F(2134)+(mup2mup4*p2mup5*p4mup3+m
     -   up2mup3*p2mup5*p4mup4+mup2mup4*p2mup3*p4mup5)*F(2135)+mup4mu
     -   p5*(p2mup3*p4mup2*F(2133)+2*(p3mup3*p4mup2*F(2136)+p3mup2*p4
     -   mup3*F(2138))))+mup3mup4*(2*(p2mup5*p4mup2*F(2037)+p2mup2*p4
     -   mup5*F(2037))+mup2mup5*F(2132)+4*(p3mup5*p4mup2*F(2137)+p3mu
     -   p2*p4mup5*F(2139)))
       Fa6(52) = 2*(mup3mup5*p2mup4*p4mup2*F(2037)+(mup2mup5*p2mup4*p4
     -   mup3+mup2mup3*p2mup5*p4mup4+mup2mup4*p2mup3*p4mup5)*F(2134)+
     -   (mup2mup4*p2mup5*p4mup3+mup2mup5*p2mup3*p4mup4+mup2mup3*p2mu
     -   p4*p4mup5)*F(2135)+2*mup3mup5*(p3mup4*p4mup2+p3mup2*p4mup4)*
     -   F(2139))+mup4mup5*(2*p2mup3*p4mup2*F(2030)+4*(p3mup3*p4mup2+
     -   p3mup2*p4mup3)*F(2138)+mup2mup3*F(2140))+mup3mup4*(2*(p2mup5
     -   *p4mup2+p2mup2*p4mup5)*F(2030)+4*(p3mup5*p4mup2+p3mup2*p4mup
     -   5)*F(2138)+mup2mup5*F(2140))+mup2mup4*mup3mup5*F(2141)
       Fa6(53) = 4*mup3mup5*p3mup4*p4mup2*F(2136)+mup2mup4*mup3mup5*F(
     -   2143)+mup4mup5*(mup2mup3*F(2142)+2*p2mup3*p4mup2*F(2145)+4*(
     -   p3mup3*p4mup2*F(2137)+p3mup2*p4mup3*F(2150)))+2*(mup2mup4*p2
     -   mup5*p4mup3*F(2146)+mup2mup4*p2mup3*p4mup5*F(2146)+mup3mup5*
     -   p2mup4*p4mup2*F(2147)+mup2mup3*(p2mup5*p4mup4+p2mup4*p4mup5)
     -   *F(2148)+mup2mup5*(p2mup3*p4mup4*F(2146)+p2mup4*p4mup3*F(214
     -   8))+2*mup3mup5*p3mup2*p4mup4*F(2151))+mup3mup4*(mup2mup5*F(2
     -   144)+2*(p2mup2*p4mup5*F(2043)+p2mup5*p4mup2*F(2149))+4*(p3mu
     -   p5*p4mup2*F(2136)+p3mup2*p4mup5*F(2151)))
       Fa6(54) = mup2mup3*p1mup4*p1mup5*F(2156)+mup2mup3*p1mup5*p2mup4
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
       Fa6(55) = p1mup2*(mup3mup5*p3mup4-mup3mup4*p3mup5)*F(2201)+mup2
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
     -   p3*p4mup5+mup2mup3*p1mup4*p4mup5)*F(4283)
       Fa6(56) = 2*(mup3mup4*p1mup2*p1mup5*F(2206)+mup3mup4*p1mup5*p2m
     -   up2*F(2206)+mup3mup4*p1mup2*p2mup5*F(2206))+mup3mup4*p1mup5*
     -   p4mup2*F(2207)+mup3mup4*p1mup2*p4mup5*F(2207)+mup3mup4*p1mup
     -   5*p3mup2*F(2209)+mup3mup4*p1mup2*p3mup5*F(2209)+mup4mup5*(p1
     -   mup3*(2*p2mup2*F(2206)+p4mup2*F(2207)+p3mup2*F(2209))+p1mup2
     -   *(2*(p1mup3+p2mup3)*F(2206)+p4mup3*F(2207)+p3mup3*F(2209)))-
     -   mup3mup5*(p1mup4*(2*p2mup2*F(2206)+p4mup2*F(2207)+p3mup2*F(2
     -   209))+p1mup2*(2*(p1mup4+p2mup4)*F(2206)+p4mup4*F(2207)+p3mup
     -   4*F(2209)))+mup2mup5*p1mup4*p3mup3*F(4282)-mup2mup4*p1mup5*p
     -   3mup3*F(4282)-mup2mup5*p1mup3*p3mup4*F(4282)+mup2mup3*p1mup5
     -   *p3mup4*F(4282)+mup2mup4*p1mup3*p3mup5*F(4282)-mup2mup3*p1mu
     -   p4*p3mup5*F(4282)+(mup2mup5*p1mup4*p4mup3-mup2mup4*p1mup5*p4
     -   mup3-mup2mup5*p1mup3*p4mup4+mup2mup3*p1mup5*p4mup4+mup2mup4*
     -   p1mup3*p4mup5-mup2mup3*p1mup4*p4mup5)*F(4283)
       Fa6(57) = mup3mup4*p1mup2*p3mup5*F(2201)-mup2mup5*p1mup4*p3mup3
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
       Fa6(58) = p1mup3*p2mup2*p3mup4*p4mup5*F(2216)+p1mup3*p2mup2*p3m
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
       Fa6(59) = -2*p2mup1*p2mup3*(p1mup4*p1mup5*F(2251)+p1mup5*p2mup4
     -   *F(2253)+2*(p1mup4*p2mup5*F(2252)+p2mup4*p2mup5*F(2254)))+p1
     -   mup3*(-2*p2mup1*(p1mup5*p2mup4*F(2246)+p1mup4*p1mup5*F(2248)
     -   +p1mup4*p2mup5*F(2249)+p2mup4*p2mup5*F(2250))+p3mup1*(p3mup4
     -   *p4mup5*F(2257)+p3mup5*p4mup4*F(2258)))+p1mup1*(2*(p1mup3*(p
     -   1mup4*p2mup5*F(2241)+2*(p1mup4*p1mup5*F(2240)+p1mup5*p2mup4*
     -   F(2242))+p2mup4*p2mup5*F(2243))+p2mup3*(p1mup4*p1mup5*F(2244
     -   )+p1mup5*p2mup4*F(2245)-p1mup4*p2mup5*F(2246)-p2mup4*p2mup5*
     -   F(2247)))+p3mup3*p3mup4*p4mup5*F(2255)+p3mup3*p3mup5*p4mup4*
     -   F(2256)+p3mup4*p3mup5*p4mup3*F(2259))
       Fa6(60) = mup3mup4*p1mup1*p2mup5*F(2263)+mup1mup5*p1mup3*p2mup4
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
       Return
       End
