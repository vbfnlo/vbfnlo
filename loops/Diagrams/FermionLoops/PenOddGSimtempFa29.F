c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim29(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa29)
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
       Complex*16 Fa29(2773:2871)
       COMMON/FaPenOddGSim/Fa
       Fa29(2773) = Fa(740)-Fa(1518)+2*Fa(1563)
       Fa29(2774) = mup1mup5*F(3155)+2*Fa(624)
       Fa29(2775) = -4*p3mup1*p4mup2*F(459)+mup1mup2*F(3156)-2*Fa(1566
     -   )
       Fa29(2776) = mup1mup2*F(3157)+4*Fa(1561)+2*Fa(1597)
       Fa29(2777) = mup1mup2*F(3158)+4*Fa(1562)+2*Fa(1596)
       Fa29(2778) = -4*p3mup1*p4mup2*F(459)+2*Fa(1572)+mup1mup2*F(477)
     -   *P(1)
       Fa29(2779) = mup1mup2*F(3155)+2*Fa(1594)
       Fa29(2780) = mup2mup4*F(3107)-2*Fa(1454)
       Fa29(2781) = mup2mup4*F(3207)+2*Fa(1439)
       Fa29(2782) = 2*Fa(1455)-mup2mup4*F(4464)*P(1)
       Fa29(2783) = mup2mup4*F(3208)+2*Fa(1441)
       Fa29(2784) = mup2mup4*F(3106)+2*Fa(1442)
       Fa29(2785) = -(mup1mup5*p4mup2*F(1609))+mup2mup5*p4mup1*F(4464)
     -   +Fa(1470)-2*Fa(1613)
       Fa29(2786) = 2*p3mup5*p4mup1*F(4464)+Fa(1524)
       Fa29(2787) = mup1mup4*F(3105)+2*Fa(1498)
       Fa29(2788) = 2*p2mup5*p4mup1*F(4464)+Fa(1526)
       Fa29(2789) = mup1mup4*F(3066)-2*Fa(1502)
       Fa29(2790) = mup1mup4*F(3064)+2*Fa(1503)
       Fa29(2791) = mup1mup4*F(3085)-2*Fa(1504)
       Fa29(2792) = Fa(730)+Fa(1505)+2*Fa(1569)
       Fa29(2793) = mup1mup5*F(3104)-2*Fa(595)
       Fa29(2794) = mup1mup2*F(3106)+4*p3mup1*p4mup2*F(4467)+2*Fa(1560
     -   )
       Fa29(2795) = mup1mup2*F(3107)+4*Fa(1567)+2*Fa(1598)
       Fa29(2796) = mup1mup2*F(3108)+4*Fa(1568)-2*Fa(1593)
       Fa29(2797) = 4*p3mup1*p4mup2*F(4467)+2*Fa(1573)-mup1mup2*F(4464
     -   )*P(1)
       Fa29(2798) = mup1mup2*F(3104)-2*Fa(1591)
       Fa29(2799) = mup2mup4*F(3062)-2*Fa(1454)
       Fa29(2800) = mup2mup4*F(3223)+2*Fa(1439)
       Fa29(2801) = mup2mup4*F(3224)+2*Fa(1441)
       Fa29(2802) = mup2mup4*F(3061)+2*Fa(1442)
       Fa29(2803) = 2*p3mup5*p4mup1*F(4464)+Fa(1527)
       Fa29(2804) = mup1mup4*F(3060)+2*Fa(1498)
       Fa29(2805) = 2*p2mup5*p4mup1*F(4464)+Fa(1529)
       Fa29(2806) = mup1mup4*F(3020)-2*Fa(1502)
       Fa29(2807) = mup1mup4*F(3018)+2*Fa(1503)
       Fa29(2808) = mup1mup4*F(3039)-2*Fa(1504)
       Fa29(2809) = mup1mup5*F(3059)-2*Fa(595)
       Fa29(2810) = mup1mup2*F(3061)+4*p3mup1*p4mup2*F(4467)+2*Fa(1560
     -   )
       Fa29(2811) = mup1mup2*F(3062)+4*Fa(1567)+2*Fa(1598)
       Fa29(2812) = mup1mup2*F(3063)+4*Fa(1568)-2*Fa(1593)
       Fa29(2813) = mup1mup2*F(3059)-2*Fa(1591)
       Fa29(2814) = mup2mup4*F(3299)+Fa(1457)
       Fa29(2815) = mup2mup4*F(3300)-2*Fa(1458)
       Fa29(2816) = mup2mup4*F(3301)+Fa(1459)
       Fa29(2817) = mup2mup4*F(3302)+Fa(1460)
       Fa29(2818) = mup2mup4*F(3303)+Fa(1461)
       Fa29(2819) = mup2mup4*F(3304)+Fa(1462)
       Fa29(2820) = -2*(mup1mup5*p4mup2*F(3481)-mup2mup5*p4mup1*F(3529
     -   ))+Fa(1473)-4*Fa(1617)
       Fa29(2821) = -4*mup1mup5*p4mup2*F(3490)-2*mup2mup5*p3mup1*F(352
     -   9)+Fa(1474)+Fa(1618)
       Fa29(2822) = Fa(1475)-2*(2*mup1mup5*p4mup2*F(3485)+Fa(1619))
       Fa29(2823) = -4*mup1mup5*p4mup2*F(3495)-2*Fa(1476)+Fa(1621)
       Fa29(2824) = Fa(1530)+Fa(1574)
       Fa29(2825) = mup1mup4*F(3010)-2*Fa(1532)
       Fa29(2826) = -2*mup4mup5*p4mup1*F(3529)+Fa(1533)+Fa(1614)
       Fa29(2827) = -2*Fa(1534)+Fa(1615)
       Fa29(2828) = 2*mup4mup5*p3mup1*F(3529)+Fa(1535)+Fa(1616)
       Fa29(2829) = Fa(1536)-2*Fa(1576)
       Fa29(2830) = mup1mup4*F(3011)+Fa(1539)
       Fa29(2831) = mup1mup4*F(3012)-2*Fa(1540)
       Fa29(2832) = Fa(1541)-2*Fa(1620)
       Fa29(2833) = 4*mup1mup4*p4mup2*F(3485)-2*Fa(1542)+Fa(1579)
       Fa29(2834) = mup1mup4*F(3008)+Fa(1543)
       Fa29(2835) = mup1mup4*F(3009)+2*Fa(1491)
       Fa29(2836) = 4*Fa(1582)-2*(Fa(1544)-Fa(1589))
       Fa29(2837) = mup1mup4*F(3000)+Fa(1545)
       Fa29(2838) = 4*mup1mup4*p4mup2*F(3490)+Fa(1546)+Fa(1583)
       Fa29(2839) = 4*mup1mup4*p4mup2*F(3495)+Fa(1547)+Fa(1584)
       Fa29(2840) = mup1mup5*F(2999)+Fa(1575)
       Fa29(2841) = mup1mup2*F(3013)-4*p3mup1*p4mup2*F(3485)-2*Fa(1577
     -   )
       Fa29(2842) = mup1mup2*F(3014)+4*p2mup1*p4mup2*F(3485)+Fa(1578)
       Fa29(2843) = mup1mup2*F(3015)+4*Fa(1580)+Fa(1600)
       Fa29(2844) = mup1mup2*F(3016)-4*Fa(1581)+Fa(1601)
       Fa29(2845) = mup1mup2*F(3017)+Fa(1585)-4*Fa(1590)
       Fa29(2846) = mup1mup2*F(2999)+Fa(1599)
       Fa29(2847) = Fa(1630)+Fa(1753)
       Fa29(2848) = Fa(1631)+Fa(1754)
       Fa29(2849) = Fa(1632)+Fa(1755)
       Fa29(2850) = Fa(1633)+Fa(1756)
       Fa29(2851) = mup3mup4*F(3559)+Fa(1750)
       Fa29(2852) = mup3mup4*F(3560)+Fa(1751)
       Fa29(2853) = mup3mup4*F(3561)+Fa(1752)
       Fa29(2854) = mup3mup4*F(3536)+Fa(1757)
       Fa29(2855) = mup3mup4*F(3538)+Fa(1758)
       Fa29(2856) = mup3mup4*F(3537)+Fa(1759)
       Fa29(2857) = mup2mup3*F(3540)+4*Fa(1810)+Fa(1875)
       Fa29(2858) = mup2mup3*F(3541)+Fa(1811)+Fa(1876)
       Fa29(2859) = mup2mup3*F(3542)+Fa(1812)-4*Fa(1877)
       Fa29(2860) = mup2mup4*F(3543)+Fa(1813)
       Fa29(2861) = mup2mup4*F(3544)+Fa(1814)
       Fa29(2862) = mup2mup4*F(3545)+Fa(1815)
       Fa29(2863) = mup2mup4*F(3586)+Fa(1816)
       Fa29(2864) = mup2mup4*F(3587)+Fa(1817)
       Fa29(2865) = mup2mup4*F(3588)+Fa(1818)
       Fa29(2866) = mup2mup3*F(3589)+Fa(1819)+Fa(1878)
       Fa29(2867) = mup2mup3*F(3590)-2*Fa(1820)-4*Fa(1879)
       Fa29(2868) = mup2mup3*F(3591)-2*Fa(1821)-4*Fa(1880)
       Fa29(2869) = mup1mup3*F(3533)+2*p3mup1*p4mup3*F(3957)-4*Fa(1928
     -   )+Fa(2020)
       Fa29(2870) = mup1mup3*F(3535)-4*p3mup1*p4mup3*F(3962)-2*Fa(1929
     -   )
       Fa29(2871) = mup1mup4*F(3530)+Fa(1930)
       Return
       End
