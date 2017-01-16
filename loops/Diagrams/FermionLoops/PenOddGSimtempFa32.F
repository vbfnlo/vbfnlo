c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim32(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa32)
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
       Complex*16 Fa32(3070:3196)
       COMMON/FaPenOddGSim/Fa
       Fa32(3070) = mup1mup4*F(3724)+2*Fa(1965)
       Fa32(3071) = -(mup1mup2*F(3722))+4*Fa(1567)+2*Fa(2112)
       Fa32(3072) = mup1mup2*F(3724)+4*p3mup1*p4mup2*F(4040)+2*Fa(2049
     -   )
       Fa32(3073) = mup1mup3*F(3720)+2*Fa(2051)
       Fa32(3074) = mup1mup2*F(3723)+2*Fa(2113)
       Fa32(3075) = Fa(1710)+Fa(1783)
       Fa32(3076) = Fa(1711)+Fa(1791)
       Fa32(3077) = Fa(1712)+Fa(1785)
       Fa32(3078) = Fa(1713)+Fa(1786)
       Fa32(3079) = mup3mup4*F(3793)+2*Fa(1780)
       Fa32(3080) = mup3mup4*F(3794)+2*Fa(1781)
       Fa32(3081) = mup3mup4*F(3795)+2*Fa(1782)
       Fa32(3082) = mup2mup3*F(3668)+4*Fa(1845)+2*Fa(1905)
       Fa32(3083) = mup2mup3*F(3669)+4*Fa(1846)+2*Fa(1906)
       Fa32(3084) = mup2mup3*F(3670)+2*Fa(1847)+4*Fa(1907)
       Fa32(3085) = mup2mup4*F(3665)+2*Fa(1848)
       Fa32(3086) = mup2mup4*F(3666)+2*Fa(1849)
       Fa32(3087) = mup2mup4*F(3667)+2*Fa(1850)
       Fa32(3088) = mup1mup3*F(3655)+4*Fa(1978)+2*(p3mup1*p4mup3*F(403
     -   8)+Fa(2064))
       Fa32(3089) = mup1mup3*F(3657)+4*p3mup1*p4mup3*F(4040)+2*Fa(1979
     -   )
       Fa32(3090) = mup1mup4*F(3658)+2*Fa(1980)
       Fa32(3091) = mup1mup4*F(3659)+2*Fa(1981)
       Fa32(3092) = mup1mup4*F(3660)+2*Fa(1982)
       Fa32(3093) = mup1mup2*F(3658)+4*Fa(1561)+2*Fa(2119)
       Fa32(3094) = mup1mup2*F(3660)+4*p3mup1*p4mup2*F(4041)+2*Fa(2063
     -   )
       Fa32(3095) = mup1mup3*F(3656)+2*Fa(2065)
       Fa32(3096) = mup1mup2*F(3659)+2*Fa(2120)
       Fa32(3097) = Fa(1726)+Fa(1793)
       Fa32(3098) = Fa(1727)+Fa(1794)
       Fa32(3099) = Fa(1728)+Fa(1795)
       Fa32(3100) = Fa(1729)+Fa(1796)
       Fa32(3101) = mup3mup4*F(3597)-2*Fa(1760)
       Fa32(3102) = mup3mup4*F(3598)+2*Fa(1761)
       Fa32(3103) = mup3mup4*F(3599)+2*Fa(1792)
       Fa32(3104) = mup3mup4*F(3653)+2*Fa(1797)
       Fa32(3105) = mup3mup4*F(3603)+2*Fa(1798)
       Fa32(3106) = mup3mup4*F(3654)-2*Fa(1799)
       Fa32(3107) = mup2mup3*F(3607)+4*Fa(1857)+2*Fa(1915)
       Fa32(3108) = mup2mup3*F(3798)+4*Fa(1846)+2*Fa(1916)
       Fa32(3109) = mup2mup3*F(3799)+2*Fa(1858)+4*Fa(1917)
       Fa32(3110) = mup2mup4*F(3604)+2*Fa(1859)
       Fa32(3111) = mup2mup4*F(3800)-2*Fa(1826)
       Fa32(3112) = mup2mup4*F(3801)+2*Fa(1860)
       Fa32(3113) = mup2mup4*F(3839)-2*Fa(1828)
       Fa32(3114) = mup2mup4*F(3652)+2*Fa(1829)
       Fa32(3115) = mup2mup4*F(3840)+2*Fa(1861)
       Fa32(3116) = mup2mup3*F(3841)+4*Fa(1862)+2*Fa(1918)
       Fa32(3117) = mup2mup3*F(3842)-2*Fa(1832)-4*Fa(1889)
       Fa32(3118) = mup2mup3*F(3648)-2*Fa(1833)-4*Fa(1890)
       Fa32(3119) = mup1mup3*F(3594)+4*Fa(1978)+2*(p3mup1*p4mup3*F(401
     -   7)-Fa(2036))
       Fa32(3120) = mup1mup3*F(3596)-4*p3mup1*p4mup3*F(4020)-2*Fa(1946
     -   )
       Fa32(3121) = mup1mup4*F(3597)-2*Fa(1996)
       Fa32(3122) = mup1mup4*F(3598)-2*Fa(1948)
       Fa32(3123) = mup1mup4*F(3599)+2*Fa(1949)
       Fa32(3124) = mup1mup4*F(3653)+2*Fa(1950)
       Fa32(3125) = mup1mup4*F(3603)-2*Fa(1951)
       Fa32(3126) = -(mup1mup4*p4mup3*F(4017))+mup3mup4*p4mup1*F(4444)
     -   +Fa(1997)
       Fa32(3127) = 2*Fa(1569)-Fa(1954)+Fa(2127)
       Fa32(3128) = 2*mup1mup4*p4mup2*F(4020)-Fa(1955)+Fa(2038)
       Fa32(3129) = mup1mup4*F(3654)+2*Fa(1998)
       Fa32(3130) = Fa(1999)+Fa(2077)
       Fa32(3131) = 2*mup1mup4*p4mup2*F(4442)+Fa(2000)+Fa(2078)
       Fa32(3132) = mup1mup3*F(3600)+4*Fa(1992)-2*Fa(2041)
       Fa32(3133) = mup1mup3*F(3601)+4*p2mup1*p4mup3*F(4025)+2*Fa(1960
     -   )
       Fa32(3134) = mup1mup3*F(3602)+4*p2mup1*p4mup3*F(4020)+2*Fa(1961
     -   )
       Fa32(3135) = mup1mup2*F(3597)+4*Fa(1561)+2*Fa(2126)
       Fa32(3136) = mup1mup2*F(3599)+4*p3mup1*p4mup2*F(4020)+2*Fa(2035
     -   )
       Fa32(3137) = mup1mup3*F(3595)+2*Fa(2037)
       Fa32(3138) = mup1mup2*F(3653)+4*Fa(1562)-2*Fa(2109)
       Fa32(3139) = 2*mup1mup2*p1mup3*F(264)+Fa(2079)+Fa(2128)
       Fa32(3140) = mup1mup2*F(3654)+2*Fa(2080)+4*Fa(2129)
       Fa32(3141) = mup1mup2*F(3603)-4*p2mup1*p4mup2*F(4020)-2*Fa(2044
     -   )
       Fa32(3142) = -2*(mup1mup3*p4mup2*F(2369)-mup1mup2*p1mup3*F(4021
     -   ))+Fa(2081)
       Fa32(3143) = mup1mup2*p1mup3*F(4018)-2*mup1mup3*p4mup2*F(4020)+
     -   Fa(2082)
       Fa32(3144) = 2*mup1mup4*p4mup2*F(4020)+mup1mup2*p2mup4*F(4444)+
     -   Fa(2083)
       Fa32(3145) = mup1mup2*F(3598)-2*Fa(2107)
       Fa32(3146) = Fa(1742)+Fa(1803)
       Fa32(3147) = Fa(1743)+Fa(1804)
       Fa32(3148) = Fa(1744)+Fa(1805)
       Fa32(3149) = Fa(1745)+Fa(1806)
       Fa32(3150) = mup3mup4*F(3862)+Fa(1800)
       Fa32(3151) = mup3mup4*F(3863)+Fa(1801)
       Fa32(3152) = mup3mup4*F(3864)+Fa(1802)
       Fa32(3153) = mup3mup4*F(3592)+Fa(1807)
       Fa32(3154) = mup3mup4*F(3539)+Fa(1808)
       Fa32(3155) = mup3mup4*F(3593)+Fa(1809)
       Fa32(3156) = mup2mup3*F(3843)+4*Fa(1863)+Fa(1921)
       Fa32(3157) = mup2mup3*F(3844)+4*Fa(1864)+Fa(1922)
       Fa32(3158) = mup2mup3*F(3845)+Fa(1865)+4*Fa(1923)
       Fa32(3159) = mup2mup4*F(3846)+Fa(1866)
       Fa32(3160) = mup2mup4*F(3847)+Fa(1867)
       Fa32(3161) = mup2mup4*F(3848)+Fa(1868)
       Fa32(3162) = mup2mup4*F(3889)+Fa(1869)
       Fa32(3163) = mup2mup4*F(3890)-2*Fa(1870)
       Fa32(3164) = mup2mup4*F(3891)+Fa(1871)
       Fa32(3165) = mup2mup3*F(3892)+4*Fa(1872)+Fa(1924)
       Fa32(3166) = mup2mup3*F(3893)+Fa(1873)+4*Fa(1879)
       Fa32(3167) = mup2mup3*F(3894)+Fa(1874)+4*Fa(1880)
       Fa32(3168) = mup1mup3*F(3530)-2*p3mup1*p4mup3*F(3957)+4*Fa(2001
     -   )+Fa(2086)
       Fa32(3169) = mup1mup3*F(3532)+4*p3mup1*p4mup3*F(3962)+Fa(2002)
       Fa32(3170) = mup1mup4*F(3533)+Fa(2003)
       Fa32(3171) = mup1mup4*F(3534)+Fa(2004)
       Fa32(3172) = mup1mup4*F(3535)-2*Fa(2005)
       Fa32(3173) = mup1mup4*F(3592)+Fa(2006)
       Fa32(3174) = mup1mup4*F(3539)+Fa(2007)
       Fa32(3175) = 2*(mup1mup4*p4mup3*F(3957)-mup3mup4*p4mup1*F(4412)
     -   )+Fa(2008)
       Fa32(3176) = 4*mup1mup4*p4mup3*F(3962)+Fa(2009)
       Fa32(3177) = Fa(2010)+4*Fa(2088)-2*Fa(2102)
       Fa32(3178) = Fa(2011)-2*(2*mup1mup4*p4mup2*F(4408)+Fa(2089))
       Fa32(3179) = mup1mup4*F(3593)+Fa(2012)
       Fa32(3180) = Fa(2013)+Fa(2090)
       Fa32(3181) = 4*m1sq*mup1mup4*p4mup2*F(2369)+Fa(2014)+Fa(2091)
       Fa32(3182) = mup1mup3*F(3536)-4*Fa(2015)+2*Fa(2026)
       Fa32(3183) = mup1mup3*F(3537)-4*p2mup1*p4mup3*F(3968)+Fa(2016)
       Fa32(3184) = mup1mup3*F(3538)-4*p2mup1*p4mup3*F(3962)-2*Fa(2017
     -   )
       Fa32(3185) = mup1mup2*F(3533)-4*Fa(2084)+Fa(2130)
       Fa32(3186) = mup1mup2*F(3535)-4*p3mup1*p4mup2*F(4408)-2*Fa(2085
     -   )
       Fa32(3187) = mup1mup3*F(3531)-2*Fa(2087)
       Fa32(3188) = mup1mup2*F(3592)+4*Fa(2092)+Fa(2132)
       Fa32(3189) = 4*mup1mup2*p1mup3*F(3951)+Fa(2093)+2*Fa(2104)
       Fa32(3190) = mup1mup2*F(3593)+Fa(2094)+4*Fa(2133)
       Fa32(3191) = mup1mup2*F(3539)+4*p2mup1*p4mup2*F(4408)+Fa(2095)
       Fa32(3192) = 4*(mup1mup2*p1mup3*F(3965)+m1sq*mup1mup3*p4mup2*F(
     -   4442))+Fa(2096)
       Fa32(3193) = -2*mup1mup2*p1mup3*F(3963)+4*mup1mup3*p4mup2*F(440
     -   8)+Fa(2097)
       Fa32(3194) = 4*mup1mup3*p4mup2*F(4414)+Fa(2098)
       Fa32(3195) = -2*mup1mup2*p2mup4*F(4412)-4*mup1mup4*p4mup2*F(441
     -   4)+Fa(2099)
       Fa32(3196) = mup1mup2*F(3534)+Fa(2131)
       Return
       End
