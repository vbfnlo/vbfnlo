c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim28(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa28)
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
       Complex*16 Fa28(2674:2772)
       COMMON/FaPenOddGSim/Fa
       Fa28(2674) = 4*mup1mup5*p4mup2*F(3485)-2*Fa(1468)+Fa(1607)
       Fa28(2675) = 4*mup1mup5*p4mup2*F(3495)+Fa(1469)+Fa(1609)
       Fa28(2676) = Fa(1477)+Fa(1548)
       Fa28(2677) = mup1mup4*F(2999)+Fa(1479)
       Fa28(2678) = Fa(1480)+2*(mup4mup5*p4mup1*F(3529)-Fa(1602))
       Fa28(2679) = Fa(1481)-2*Fa(1603)
       Fa28(2680) = -2*mup4mup5*p3mup1*F(3529)+Fa(1482)+Fa(1604)
       Fa28(2681) = Fa(1483)+Fa(1550)
       Fa28(2682) = mup1mup4*F(3013)-2*Fa(1486)
       Fa28(2683) = mup1mup4*F(3014)+Fa(1487)
       Fa28(2684) = -2*Fa(1488)+Fa(1608)
       Fa28(2685) = Fa(1489)-2*(2*mup1mup4*p4mup2*F(3485)+Fa(1553))
       Fa28(2686) = mup1mup4*F(3015)+Fa(1490)
       Fa28(2687) = mup1mup4*F(3016)-2*Fa(1491)
       Fa28(2688) = Fa(1492)-2*(2*Fa(1556)+Fa(1589))
       Fa28(2689) = mup1mup4*F(3017)+Fa(1493)
       Fa28(2690) = -4*mup1mup4*p4mup2*F(3490)+Fa(1494)+Fa(1557)
       Fa28(2691) = -4*mup1mup4*p4mup2*F(3495)-2*Fa(1495)+Fa(1558)
       Fa28(2692) = mup1mup5*F(3010)-2*Fa(1549)
       Fa28(2693) = mup1mup2*F(3011)+4*p3mup1*p4mup2*F(3485)+Fa(1551)
       Fa28(2694) = mup1mup2*F(3012)-4*p2mup1*p4mup2*F(3485)-2*Fa(1552
     -   )
       Fa28(2695) = mup1mup2*F(3008)-4*Fa(1554)-2*Fa(1587)
       Fa28(2696) = mup1mup2*F(3009)+4*Fa(1555)+Fa(1588)
       Fa28(2697) = mup1mup2*F(3000)+Fa(1559)+4*Fa(1590)
       Fa28(2698) = mup1mup2*F(3010)-2*Fa(1586)
       Fa28(2699) = mup2mup4*F(3018)+2*Fa(1438)
       Fa28(2700) = mup2mup4*F(3057)-2*Fa(1439)
       Fa28(2701) = 2*Fa(1440)+mup2mup4*F(4464)*P(1)
       Fa28(2702) = mup2mup4*F(3058)-2*Fa(1441)
       Fa28(2703) = mup2mup4*F(3020)-2*Fa(1442)
       Fa28(2704) = -(mup1mup5*p4mup2*F(1609))+mup2mup5*p4mup1*F(4464)
     -   +Fa(1470)+2*Fa(1611)
       Fa28(2705) = -(mup2mup5*p3mup1)+mup1mup5*p3mup2+mup1mup2*p3mup5
       Fa28(2706) = 2*p3mup5*p4mup1*F(4464)+Fa(1496)
       Fa28(2707) = mup1mup4*F(3059)-2*Fa(1498)
       Fa28(2708) = mup4mup5*p4mup1*F(4464)+Fa(1499)-Fa(1610)
       Fa28(2709) = mup4mup5*p3mup1+mup1mup5*p3mup4-mup1mup4*p3mup5
       Fa28(2710) = 2*p2mup5*p4mup1*F(4464)+Fa(1500)
       Fa28(2711) = mup1mup4*F(3061)+2*Fa(1502)
       Fa28(2712) = mup1mup4*F(3062)-2*Fa(1503)
       Fa28(2713) = mup1mup4*F(3063)+2*Fa(1504)
       Fa28(2714) = Fa(744)-Fa(1505)+2*Fa(1563)
       Fa28(2715) = 2*Fa(1506)-mup1mup4*P(1)
       Fa28(2716) = -(mup1mup2*p3mup4)+Fa(1564)
       Fa28(2717) = mup1mup5*F(3060)+2*Fa(595)
       Fa28(2718) = mup1mup2*F(3020)-4*p3mup1*p4mup2*F(4467)-2*Fa(1560
     -   )
       Fa28(2719) = mup1mup2*F(3018)+4*Fa(1561)+2*Fa(1592)
       Fa28(2720) = mup1mup2*F(3039)+4*Fa(1562)+2*Fa(1593)
       Fa28(2721) = -4*p3mup1*p4mup2*F(4467)+2*Fa(1565)+mup1mup2*F(446
     -   4)*P(1)
       Fa28(2722) = mup1mup2*F(3060)+2*Fa(1591)
       Fa28(2723) = mup2mup4*F(3064)+2*Fa(1438)
       Fa28(2724) = mup2mup4*F(3102)-2*Fa(1439)
       Fa28(2725) = mup2mup4*F(3103)-2*Fa(1441)
       Fa28(2726) = mup2mup4*F(3066)-2*Fa(1442)
       Fa28(2727) = 2*p3mup5*p4mup1*F(4464)+Fa(1507)
       Fa28(2728) = mup1mup4*F(3104)-2*Fa(1498)
       Fa28(2729) = 2*p2mup5*p4mup1*F(4464)+Fa(1509)
       Fa28(2730) = mup1mup4*F(3106)+2*Fa(1502)
       Fa28(2731) = mup1mup4*F(3107)-2*Fa(1503)
       Fa28(2732) = mup1mup4*F(3108)+2*Fa(1504)
       Fa28(2733) = mup1mup5*F(3105)+2*Fa(595)
       Fa28(2734) = mup1mup2*F(3066)-4*p3mup1*p4mup2*F(4467)-2*Fa(1560
     -   )
       Fa28(2735) = mup1mup2*F(3064)+4*Fa(1561)+2*Fa(1592)
       Fa28(2736) = mup1mup2*F(3085)+4*Fa(1562)+2*Fa(1593)
       Fa28(2737) = mup1mup2*F(3105)+2*Fa(1591)
       Fa28(2738) = mup2mup4*F(3109)-2*Fa(1445)
       Fa28(2739) = mup2mup4*F(3153)+2*Fa(1446)
       Fa28(2740) = 2*Fa(1447)-mup2mup4*F(477)*P(1)
       Fa28(2741) = mup2mup4*F(3154)+2*Fa(1448)
       Fa28(2742) = mup2mup4*F(3113)+2*Fa(1449)
       Fa28(2743) = mup2mup5*p4mup1*F(477)-mup1mup5*p4mup2*F(1623)+Fa(
     -   1471)-2*Fa(1613)
       Fa28(2744) = -2*p3mup5*p4mup1*F(477)+Fa(1510)
       Fa28(2745) = mup1mup4*F(3155)+2*Fa(1512)
       Fa28(2746) = mup4mup5*p4mup1*F(477)+Fa(1513)-Fa(1612)
       Fa28(2747) = 2*p2mup5*p4mup1*F(477)+Fa(1514)
       Fa28(2748) = mup1mup4*F(3156)-2*Fa(1515)
       Fa28(2749) = mup1mup4*F(3157)+2*Fa(1516)
       Fa28(2750) = mup1mup4*F(3158)-2*Fa(1517)
       Fa28(2751) = Fa(742)+Fa(1518)+2*Fa(1569)
       Fa28(2752) = 2*Fa(1519)+mup1mup4*P(1)
       Fa28(2753) = mup1mup2*p3mup4+Fa(1570)
       Fa28(2754) = mup1mup5*F(3110)-2*Fa(624)
       Fa28(2755) = 4*p3mup1*p4mup2*F(459)+mup1mup2*F(3113)+2*Fa(1566)
       Fa28(2756) = mup1mup2*F(3109)+4*Fa(1567)+2*Fa(1595)
       Fa28(2757) = mup1mup2*F(3134)+4*Fa(1568)-2*Fa(1596)
       Fa28(2758) = 4*p3mup1*p4mup2*F(459)+2*Fa(1571)-mup1mup2*F(477)*
     -   P(1)
       Fa28(2759) = mup1mup2*F(3110)-2*Fa(1594)
       Fa28(2760) = mup2mup4*F(3157)+2*Fa(1451)
       Fa28(2761) = mup2mup4*F(3191)-2*Fa(1446)
       Fa28(2762) = -2*Fa(1447)+mup2mup4*F(477)*P(1)
       Fa28(2763) = mup2mup4*F(3192)-2*Fa(1452)
       Fa28(2764) = mup2mup4*F(3156)-2*Fa(1449)
       Fa28(2765) = mup2mup5*p4mup1*F(477)-mup1mup5*p4mup2*F(1623)+Fa(
     -   1472)+2*Fa(1611)
       Fa28(2766) = 2*p3mup5*p4mup1*F(477)+Fa(1520)
       Fa28(2767) = mup1mup4*F(3110)-2*Fa(1512)
       Fa28(2768) = mup4mup5*p4mup1*F(477)+Fa(1522)-Fa(1612)
       Fa28(2769) = -2*p2mup5*p4mup1*F(477)+Fa(1523)
       Fa28(2770) = mup1mup4*F(3113)+2*Fa(1515)
       Fa28(2771) = mup1mup4*F(3109)-2*Fa(1516)
       Fa28(2772) = mup1mup4*F(3134)+2*Fa(1517)
       Return
       End
