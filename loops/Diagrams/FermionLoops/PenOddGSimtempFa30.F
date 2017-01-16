c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim30(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa30)
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
       Complex*16 Fa30(2872:2970)
       COMMON/FaPenOddGSim/Fa
       Fa30(2872) = mup1mup4*F(3531)-2*Fa(1931)
       Fa30(2873) = mup1mup4*F(3532)+Fa(1932)
       Fa30(2874) = mup1mup4*F(3536)+Fa(1933)
       Fa30(2875) = mup1mup4*F(3538)-2*Fa(1934)
       Fa30(2876) = -2*(mup1mup4*p4mup3*F(3957)-mup3mup4*p4mup1*F(4412
     -   ))+Fa(1935)
       Fa30(2877) = -4*mup1mup4*p4mup3*F(3962)+Fa(1936)
       Fa30(2878) = Fa(1937)-4*Fa(2022)+2*Fa(2102)
       Fa30(2879) = 4*mup1mup4*p4mup2*F(4408)-2*Fa(1938)+Fa(2023)
       Fa30(2880) = mup1mup4*F(3537)+Fa(1939)
       Fa30(2881) = Fa(1940)+Fa(2024)
       Fa30(2882) = 4*m1sq*mup1mup4*p4mup2*F(4442)+Fa(1941)+Fa(2025)
       Fa30(2883) = mup1mup3*F(3592)+4*Fa(1942)-2*Fa(2026)
       Fa30(2884) = mup1mup3*F(3593)+4*p2mup1*p4mup3*F(3968)+Fa(1943)
       Fa30(2885) = mup1mup3*F(3539)+4*p2mup1*p4mup3*F(3962)+Fa(1944)
       Fa30(2886) = mup1mup2*F(3530)+4*Fa(2018)+Fa(2100)
       Fa30(2887) = mup1mup2*F(3532)+4*p3mup1*p4mup2*F(4408)+Fa(2019)
       Fa30(2888) = mup1mup3*F(3534)+Fa(2021)
       Fa30(2889) = mup1mup2*F(3536)-4*Fa(2027)+Fa(2103)
       Fa30(2890) = Fa(2028)-2*(2*mup1mup2*p1mup3*F(3951)+Fa(2104))
       Fa30(2891) = mup1mup2*F(3537)+Fa(2029)+4*Fa(2105)
       Fa30(2892) = mup1mup2*F(3538)-4*p2mup1*p4mup2*F(4408)-2*Fa(2030
     -   )
       Fa30(2893) = 4*(m1sq*mup1mup3*p4mup2*F(2369)-mup1mup2*p1mup3*F(
     -   3965))+Fa(2031)
       Fa30(2894) = 2*mup1mup2*p1mup3*F(3963)-4*mup1mup3*p4mup2*F(4408
     -   )+Fa(2032)
       Fa30(2895) = 2*mup1mup3*p4mup2*F(4414)+Fa(2033)
       Fa30(2896) = 2*mup1mup2*p2mup4*F(4412)+4*mup1mup4*p4mup2*F(4414
     -   )+Fa(2034)
       Fa30(2897) = mup1mup2*F(3531)-2*Fa(2101)
       Fa30(2898) = -Fa(1646)+Fa(1763)
       Fa30(2899) = Fa(1647)+Fa(1764)
       Fa30(2900) = Fa(1648)+Fa(1765)
       Fa30(2901) = Fa(1649)+Fa(1766)
       Fa30(2902) = mup3mup4*F(3594)+2*Fa(1760)
       Fa30(2903) = mup3mup4*F(3595)-2*Fa(1761)
       Fa30(2904) = mup3mup4*F(3596)-2*Fa(1762)
       Fa30(2905) = mup3mup4*F(3600)+2*Fa(1767)
       Fa30(2906) = mup3mup4*F(3602)+2*Fa(1768)
       Fa30(2907) = mup3mup4*F(3601)+2*Fa(1769)
       Fa30(2908) = mup2mup3*F(3604)+4*Fa(1822)+2*Fa(1885)
       Fa30(2909) = mup2mup3*F(3605)+4*Fa(1823)+2*Fa(1886)
       Fa30(2910) = mup2mup3*F(3606)+2*Fa(1824)+4*Fa(1887)
       Fa30(2911) = mup2mup4*F(3607)-2*Fa(1825)
       Fa30(2912) = mup2mup4*F(3608)+2*Fa(1826)
       Fa30(2913) = mup2mup4*F(3609)+2*Fa(1827)
       Fa30(2914) = mup2mup4*F(3647)+2*Fa(1828)
       Fa30(2915) = mup2mup4*F(3648)-2*Fa(1829)
       Fa30(2916) = mup2mup4*F(3649)+2*Fa(1830)
       Fa30(2917) = mup2mup3*F(3650)+4*Fa(1831)+2*Fa(1888)
       Fa30(2918) = mup2mup3*F(3651)+2*Fa(1832)+4*Fa(1889)
       Fa30(2919) = mup2mup3*F(3652)+2*Fa(1833)+4*Fa(1890)
       Fa30(2920) = mup1mup3*F(3597)+4*Fa(1945)-2*(p3mup1*p4mup3*F(401
     -   7)-Fa(2036))
       Fa30(2921) = mup1mup3*F(3599)+4*p3mup1*p4mup3*F(4020)+2*Fa(1946
     -   )
       Fa30(2922) = mup1mup4*F(3594)+2*Fa(1947)
       Fa30(2923) = mup1mup4*F(3595)+2*Fa(1948)
       Fa30(2924) = mup1mup4*F(3596)-2*Fa(1949)
       Fa30(2925) = mup1mup4*F(3600)-2*Fa(1950)
       Fa30(2926) = mup1mup4*F(3602)+2*Fa(1951)
       Fa30(2927) = mup1mup4*p4mup3*F(4017)-mup3mup4*p4mup1*F(4444)+Fa
     -   (1952)
       Fa30(2928) = -2*mup1mup4*p4mup3*F(4020)+Fa(1953)
       Fa30(2929) = 2*Fa(1563)+Fa(1954)+Fa(2108)
       Fa30(2930) = -2*mup1mup4*p4mup2*F(4020)+Fa(1955)-Fa(2038)
       Fa30(2931) = mup1mup4*F(3601)+2*Fa(1956)
       Fa30(2932) = Fa(1957)+Fa(2039)
       Fa30(2933) = 2*mup1mup4*p4mup2*F(2369)+Fa(1958)+Fa(2040)
       Fa30(2934) = mup1mup3*F(3653)+4*Fa(1959)+2*Fa(2041)
       Fa30(2935) = mup1mup3*F(3654)-4*p2mup1*p4mup3*F(4025)-2*Fa(1960
     -   )
       Fa30(2936) = mup1mup3*F(3603)-4*p2mup1*p4mup3*F(4020)-2*Fa(1961
     -   )
       Fa30(2937) = mup1mup2*F(3594)+4*Fa(1567)+2*Fa(2106)
       Fa30(2938) = mup1mup2*F(3596)-4*p3mup1*p4mup2*F(4020)-2*Fa(2035
     -   )
       Fa30(2939) = mup1mup3*F(3598)-2*Fa(2037)
       Fa30(2940) = mup1mup2*F(3600)+4*Fa(1568)+2*Fa(2109)
       Fa30(2941) = 2*mup1mup2*p1mup3*F(258)+Fa(2042)+Fa(2110)
       Fa30(2942) = mup1mup2*F(3601)+2*Fa(2043)+4*Fa(2111)
       Fa30(2943) = mup1mup2*F(3602)+4*p2mup1*p4mup2*F(4020)+2*Fa(2044
     -   )
       Fa30(2944) = 2*(mup1mup2*p1mup3*F(4021)+mup1mup3*p4mup2*F(4442)
     -   )+Fa(2045)
       Fa30(2945) = -(mup1mup2*p1mup3*F(4018))+2*mup1mup3*p4mup2*F(402
     -   0)+Fa(2046)
       Fa30(2946) = 2*mup1mup3*p4mup2*F(4020)+Fa(2047)
       Fa30(2947) = -2*mup1mup4*p4mup2*F(4020)-mup1mup2*p2mup4*F(4444)
     -   +Fa(2048)
       Fa30(2948) = mup1mup2*F(3595)+2*Fa(2107)
       Fa30(2949) = Fa(1662)+Fa(1773)
       Fa30(2950) = Fa(1663)+Fa(1774)
       Fa30(2951) = Fa(1664)+Fa(1775)
       Fa30(2952) = Fa(1665)+Fa(1776)
       Fa30(2953) = mup3mup4*F(3684)+2*Fa(1770)
       Fa30(2954) = mup3mup4*F(3685)+2*Fa(1771)
       Fa30(2955) = mup3mup4*F(3686)+2*Fa(1772)
       Fa30(2956) = mup3mup4*F(3661)+2*Fa(1777)
       Fa30(2957) = mup3mup4*F(3663)+2*Fa(1778)
       Fa30(2958) = mup3mup4*F(3662)+2*Fa(1779)
       Fa30(2959) = mup2mup3*F(3665)+4*Fa(1834)+2*Fa(1895)
       Fa30(2960) = mup2mup3*F(3666)+4*Fa(1823)+2*Fa(1896)
       Fa30(2961) = mup2mup3*F(3667)+2*Fa(1835)+4*Fa(1897)
       Fa30(2962) = mup2mup4*F(3668)+2*Fa(1836)
       Fa30(2963) = mup2mup4*F(3669)+2*Fa(1837)
       Fa30(2964) = mup2mup4*F(3670)+2*Fa(1838)
       Fa30(2965) = mup2mup4*F(3711)+2*Fa(1839)
       Fa30(2966) = mup2mup4*F(3712)+2*Fa(1840)
       Fa30(2967) = mup2mup4*F(3713)+2*Fa(1841)
       Fa30(2968) = mup2mup3*F(3714)+4*Fa(1842)+2*Fa(1898)
       Fa30(2969) = mup2mup3*F(3715)+2*Fa(1843)+4*Fa(1899)
       Fa30(2970) = mup2mup3*F(3716)+2*Fa(1844)+4*Fa(1900)
       Return
       End
