c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim31(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa31)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenOddGSim/F
      Complex*16 EpsF(126)
       common/EPSFPenOddGSim/EPSF
       COMMON/InvariantsPenOddGSim/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s3
     -   4,s45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       Real*8   p1sq 
        Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s234,m1,m1sq,m3,m3sq,m4,m4sq
      Complex*16 Fa(3196),F(4476)
      Real*8 P(1383) 
        COMMON/PPenOddGSim/P
       Complex*16 Fa31(2971:3069)
       COMMON/FaPenOddGSim/Fa
       Fa31(2971) = mup1mup3*F(3658)+4*Fa(1945)+2*(p3mup1*p4mup3*F(403
     -   9)+Fa(2050))
       Fa31(2972) = mup1mup3*F(3660)+4*p3mup1*p4mup3*F(4041)+2*Fa(1962
     -   )
       Fa31(2973) = mup1mup4*F(3655)+2*Fa(1963)
       Fa31(2974) = mup1mup4*F(3656)+2*Fa(1964)
       Fa31(2975) = mup1mup4*F(3657)+2*Fa(1965)
       Fa31(2976) = mup1mup4*F(3661)+2*Fa(1966)
       Fa31(2977) = mup1mup4*F(3663)+2*Fa(1967)
       Fa31(2978) = mup3mup4*p4mup1*F(1056)+mup1mup4*p4mup3*F(4038)+Fa
     -   (1968)
       Fa31(2979) = 2*mup1mup4*p4mup3*F(4041)+Fa(1969)
       Fa31(2980) = 2*Fa(1563)+Fa(1970)+Fa(2114)
       Fa31(2981) = 2*mup1mup4*p4mup2*F(4040)+Fa(1971)+Fa(2052)
       Fa31(2982) = mup1mup4*F(3662)+2*Fa(1972)
       Fa31(2983) = Fa(1973)+Fa(2053)
       Fa31(2984) = 2*mup1mup4*p4mup2*F(2369)+Fa(1974)+Fa(2054)
       Fa31(2985) = mup1mup3*F(3717)+4*Fa(1959)+2*Fa(2055)
       Fa31(2986) = 4*p2mup1*p4mup3*F(362)+mup1mup3*F(3718)+2*Fa(1975)
       Fa31(2987) = mup1mup3*F(3664)+4*p2mup1*p4mup3*F(4040)+2*Fa(1976
     -   )
       Fa31(2988) = mup1mup2*F(3655)+4*Fa(1567)+2*Fa(2112)
       Fa31(2989) = mup1mup2*F(3657)+4*p3mup1*p4mup2*F(4040)+2*Fa(2049
     -   )
       Fa31(2990) = mup1mup3*F(3659)+2*Fa(2051)
       Fa31(2991) = mup1mup2*F(3661)+4*Fa(1568)+2*Fa(2115)
       Fa31(2992) = 2*mup1mup2*p1mup3*F(258)+Fa(2056)+Fa(2116)
       Fa31(2993) = mup1mup2*F(3662)+2*Fa(2057)+4*Fa(2117)
       Fa31(2994) = mup1mup2*F(3663)+4*p2mup1*p4mup2*F(4041)+2*Fa(2058
     -   )
       Fa31(2995) = 2*(mup1mup2*p1mup3*F(346)+mup1mup3*p4mup2*F(4442))
     -   +Fa(2059)
       Fa31(2996) = mup1mup2*p1mup3*F(4029)+2*mup1mup3*p4mup2*F(4041)+
     -   Fa(2060)
       Fa31(2997) = 2*mup1mup3*p4mup2*F(4041)+Fa(2061)
       Fa31(2998) = mup1mup2*p2mup4*F(1056)+2*mup1mup4*p4mup2*F(4040)+
     -   Fa(2062)
       Fa31(2999) = mup1mup2*F(3656)+2*Fa(2113)
       Fa31(3000) = mup3mup4*p2mup1*F(1056)+Fa(1977)
       Fa31(3001) = Fa(1678)+Fa(1783)
       Fa31(3002) = Fa(1679)+Fa(1784)
       Fa31(3003) = Fa(1680)+Fa(1785)
       Fa31(3004) = Fa(1681)+Fa(1786)
       Fa31(3005) = mup3mup4*F(3744)+2*Fa(1780)
       Fa31(3006) = mup3mup4*F(3745)+2*Fa(1781)
       Fa31(3007) = mup3mup4*F(3746)+2*Fa(1782)
       Fa31(3008) = mup3mup4*F(3717)+2*Fa(1787)
       Fa31(3009) = mup3mup4*F(3664)+2*Fa(1788)
       Fa31(3010) = mup3mup4*F(3718)+2*Fa(1789)
       Fa31(3011) = mup2mup3*F(3725)+4*Fa(1845)+2*Fa(1905)
       Fa31(3012) = mup2mup3*F(3726)+4*Fa(1846)+2*Fa(1906)
       Fa31(3013) = mup2mup3*F(3727)+2*Fa(1847)+4*Fa(1907)
       Fa31(3014) = -(mup2mup4*F(3728))+2*Fa(1848)
       Fa31(3015) = mup2mup4*F(3729)+2*Fa(1849)
       Fa31(3016) = mup2mup4*F(3730)+2*Fa(1850)
       Fa31(3017) = mup2mup4*F(3714)+2*Fa(1851)
       Fa31(3018) = mup2mup4*F(3716)+2*Fa(1852)
       Fa31(3019) = mup2mup4*F(3715)+2*Fa(1853)
       Fa31(3020) = mup2mup3*F(3711)+4*Fa(1854)+2*Fa(1908)
       Fa31(3021) = mup2mup3*F(3713)+2*Fa(1855)+4*Fa(1909)
       Fa31(3022) = mup2mup3*F(3712)+2*Fa(1856)+4*Fa(1910)
       Fa31(3023) = -(mup1mup3*F(3722))+4*Fa(1978)+2*(p3mup1*p4mup3*F(
     -   4038)+Fa(2064))
       Fa31(3024) = mup1mup3*F(3724)+4*p3mup1*p4mup3*F(4040)+2*Fa(1979
     -   )
       Fa31(3025) = mup1mup4*F(3719)+2*Fa(1980)
       Fa31(3026) = mup1mup4*F(3720)+2*Fa(1981)
       Fa31(3027) = mup1mup4*F(3721)+2*Fa(1982)
       Fa31(3028) = mup1mup4*F(3717)+2*Fa(1983)
       Fa31(3029) = mup1mup4*F(3664)+2*Fa(1984)
       Fa31(3030) = mup3mup4*p4mup1*F(1057)+mup1mup4*p4mup3*F(4039)+Fa
     -   (1985)
       Fa31(3031) = 2*mup1mup4*p4mup3*F(4040)+Fa(1986)
       Fa31(3032) = 2*Fa(1569)+Fa(1987)+Fa(2121)
       Fa31(3033) = 2*mup1mup4*p4mup2*F(4041)+Fa(1988)+Fa(2066)
       Fa31(3034) = mup1mup4*F(3718)+2*Fa(1989)
       Fa31(3035) = Fa(1990)+Fa(2067)
       Fa31(3036) = 2*mup1mup4*p4mup2*F(4442)+Fa(1991)+Fa(2068)
       Fa31(3037) = mup1mup3*F(3661)+4*Fa(1992)+2*Fa(2069)
       Fa31(3038) = mup1mup3*F(3662)+4*p2mup1*p4mup3*F(3969)+2*Fa(1993
     -   )
       Fa31(3039) = mup1mup3*F(3663)+4*p2mup1*p4mup3*F(4041)+2*Fa(1994
     -   )
       Fa31(3040) = mup1mup2*F(3719)+4*Fa(1561)+2*Fa(2119)
       Fa31(3041) = mup1mup2*F(3721)+4*p3mup1*p4mup2*F(4041)+2*Fa(2063
     -   )
       Fa31(3042) = mup1mup3*F(3723)+2*Fa(2065)
       Fa31(3043) = mup1mup2*F(3717)+4*Fa(1562)+2*Fa(2122)
       Fa31(3044) = 2*mup1mup2*p1mup3*F(264)+Fa(2070)+Fa(2123)
       Fa31(3045) = mup1mup2*F(3718)+2*Fa(2071)+4*Fa(2124)
       Fa31(3046) = mup1mup2*F(3664)+4*p2mup1*p4mup2*F(4040)+2*Fa(2072
     -   )
       Fa31(3047) = 2*(mup1mup3*p4mup2*F(2369)+mup1mup2*p1mup3*F(4024)
     -   )+Fa(2073)
       Fa31(3048) = mup1mup2*p1mup3*F(4028)+2*mup1mup3*p4mup2*F(4040)+
     -   Fa(2074)
       Fa31(3049) = 2*mup1mup3*p4mup2*F(4040)+Fa(2075)
       Fa31(3050) = mup1mup2*p2mup4*F(1057)+2*mup1mup4*p4mup2*F(4041)+
     -   Fa(2076)
       Fa31(3051) = mup1mup2*F(3720)+2*Fa(2120)
       Fa31(3052) = mup3mup4*p2mup1*F(1057)+Fa(1995)
       Fa31(3053) = Fa(1694)+Fa(1773)
       Fa31(3054) = Fa(1695)+M1*M3*Fa(1790)
       Fa31(3055) = Fa(1696)+Fa(1775)
       Fa31(3056) = Fa(1697)+Fa(1776)
       Fa31(3057) = -(mup3mup4*F(3776))+2*Fa(1770)
       Fa31(3058) = mup3mup4*F(3777)+2*Fa(1771)
       Fa31(3059) = mup3mup4*F(3778)+2*Fa(1772)
       Fa31(3060) = -(mup2mup3*F(3728))+4*Fa(1834)+2*Fa(1895)
       Fa31(3061) = mup2mup3*F(3729)+4*Fa(1823)+2*Fa(1896)
       Fa31(3062) = mup2mup3*F(3730)+2*Fa(1835)+4*Fa(1897)
       Fa31(3063) = mup2mup4*F(3725)+2*Fa(1836)
       Fa31(3064) = mup2mup4*F(3726)+2*Fa(1837)
       Fa31(3065) = mup2mup4*F(3727)+2*Fa(1838)
       Fa31(3066) = mup1mup3*F(3719)+4*Fa(1945)+2*(p3mup1*p4mup3*F(403
     -   9)+Fa(2050))
       Fa31(3067) = mup1mup3*F(3721)+4*p3mup1*p4mup3*F(4041)+2*Fa(1962
     -   )
       Fa31(3068) = -(mup1mup4*F(3722))+2*Fa(1963)
       Fa31(3069) = mup1mup4*F(3723)+2*Fa(1964)
       Return
       End
