c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenOddGSim16(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa16)
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
       Complex*16 Fa16(1486:1584)
       COMMON/FaPenOddGSim/Fa
       Fa16(1486) = p3mup1*(p1mup4*F(3482)+2*(p2mup4*F(3483)+p3mup4*F(
     -   3484)+p4mup4*F(3485)))
       Fa16(1487) = 2*p2mup1*(p1mup4*F(3482)+2*(p2mup4*F(3483)+p3mup4*
     -   F(3484)+p4mup4*F(3485)))
       Fa16(1488) = mup1mup4*(2*(p1mup5*F(2527)+p4mup5*F(2531)+p3mup5*
     -   F(2566))+p2mup5*F(3493))-mup4mup5*p2mup1*F(3529)
       Fa16(1489) = 2*mup1mup2*(p1mup4*F(3482)+2*(p2mup4*F(3483)+p3mup
     -   4*F(3484)+p4mup4*F(3485)))
       Fa16(1490) = 4*p1mup4*(p4mup1*F(3001)+p3mup1*F(3474))+2*p3mup1*
     -   (2*(p2mup4*F(3475)+p3mup4*F(3478))+p4mup4*F(3479))+p4mup1*(-
     -   4*(p4mup4*F(3007)-p2mup4*F(3486))+p3mup4*F(3487))
       Fa16(1491) = 2*p1mup4*(p4mup1*F(2525)+p2mup1*F(3474))+p2mup1*(2
     -   *(p2mup4*F(3475)+p3mup4*F(3478))+p4mup4*F(3479))+p4mup1*(-2*
     -   (p4mup4*F(2535)-p3mup4*F(2565))+p2mup4*F(3492))
       Fa16(1492) = 2*mup1mup2*(2*(p1mup4*F(3474)+p2mup4*F(3475)+p3mup
     -   4*F(3478))+p4mup4*F(3479))
       Fa16(1493) = 4*(p1mup4*p3mup1*F(2525)-p3mup1*p4mup4*F(2535)+p3m
     -   up1*p3mup4*F(2565)-p1mup4*p2mup1*F(3001)+p2mup1*p4mup4*F(300
     -   7)-p2mup1*p2mup4*F(3486))+2*p2mup4*p3mup1*F(3492)+p2mup1*p3m
     -   up4*F(3527)
       Fa16(1494) = mup1mup2*(4*(p1mup4*F(3001)-p4mup4*F(3007)+p2mup4*
     -   F(3486))+p3mup4*F(3487))
       Fa16(1495) = mup1mup2*(2*(p1mup4*F(2525)-p4mup4*F(2535)+p3mup4*
     -   F(2565))+p2mup4*F(3492))
       Fa16(1496) = 2*p3mup1*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(383))+mup1mup5*F(3018)
       Fa16(1497) = 2*p1mup5*p3mup1*F(383)-mup1mup5*F(3020)+4*p3mup1*(
     -   p3mup5*F(4465)+p4mup5*F(4467)+p2mup5*F(4469))
       Fa16(1498) = p4mup1*(p1mup4*F(383)+2*(p3mup4*F(4465)+p4mup4*F(4
     -   467)+p2mup4*F(4469)))
       Fa16(1499) = mup1mup4*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(383))
       Fa16(1500) = 2*p2mup1*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(383))-mup1mup5*F(3039)
       Fa16(1501) = 2*(p2mup5*p3mup1-p2mup1*p3mup5)+mup1mup5*P(1)
       Fa16(1502) = p3mup1*(p1mup4*F(383)+2*(p3mup4*F(4465)+p4mup4*F(4
     -   467)+p2mup4*F(4469)))
       Fa16(1503) = p3mup1*(2*(p1mup4*F(263)+p3mup4*F(273)+p2mup4*F(37
     -   9))+p4mup4*F(383))+p3mup4*p4mup1*F(4464)
       Fa16(1504) = p2mup1*(2*(p1mup4*F(263)+p3mup4*F(273)+p2mup4*F(37
     -   9))+p4mup4*F(383))+p2mup4*p4mup1*F(4464)
       Fa16(1505) = mup1mup2*(2*(p1mup4*F(263)+p3mup4*F(273)+p2mup4*F(
     -   379))+p4mup4*F(383))
       Fa16(1506) = -(p2mup4*p3mup1)+p2mup1*p3mup4
       Fa16(1507) = 2*p3mup1*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(383))+mup1mup5*F(3064)
       Fa16(1508) = 2*p1mup5*p3mup1*F(383)-mup1mup5*F(3066)+4*p3mup1*(
     -   p3mup5*F(4465)+p4mup5*F(4467)+p2mup5*F(4469))
       Fa16(1509) = 2*p2mup1*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(383))-mup1mup5*F(3085)
       Fa16(1510) = -4*p3mup1*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(37
     -   9))+2*p3mup1*p4mup5*F(485)+mup1mup5*F(3109)
       Fa16(1511) = 2*p1mup5*p3mup1*F(443)+4*p3mup1*(p2mup5*F(444)+p3m
     -   up5*F(445)+p4mup5*F(459))+mup1mup5*F(3113)
       Fa16(1512) = p4mup1*(p1mup4*F(443)+2*(p2mup4*F(444)+p3mup4*F(44
     -   5)+p4mup4*F(459)))
       Fa16(1513) = mup1mup4*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(447))
       Fa16(1514) = 2*p2mup1*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(447))+mup1mup5*F(3134)
       Fa16(1515) = p3mup1*(p1mup4*F(443)+2*(p2mup4*F(444)+p3mup4*F(44
     -   5)+p4mup4*F(459)))
       Fa16(1516) = p3mup4*p4mup1*F(477)+p3mup1*(p4mup4*F(448)+2*(p1mu
     -   p4*F(376)+p2mup4*F(377)+p3mup4*F(3135)))
       Fa16(1517) = 2*p1mup4*p2mup1*F(376)+p2mup4*p4mup1*F(477)+p2mup1
     -   *(p4mup4*F(448)+2*(p2mup4*F(377)+p3mup4*F(3135)))
       Fa16(1518) = mup1mup2*(p4mup4*F(448)+2*(p1mup4*F(376)+p2mup4*F(
     -   377)+p3mup4*F(3135)))
       Fa16(1519) = p2mup4*p3mup1-p2mup1*p3mup4
       Fa16(1520) = 2*p3mup1*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(447))+mup1mup5*F(3157)
       Fa16(1521) = 2*p1mup5*p3mup1*F(443)+4*p3mup1*(p2mup5*F(444)+p3m
     -   up5*F(445)+p4mup5*F(459))-mup1mup5*F(3156)
       Fa16(1522) = mup1mup4*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))-p4mup5*F(485))
       Fa16(1523) = -4*p2mup1*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(37
     -   9))+2*p2mup1*p4mup5*F(485)+mup1mup5*F(3158)
       Fa16(1524) = 2*p3mup1*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(383))-mup1mup5*F(3107)
       Fa16(1525) = 2*p1mup5*p3mup1*F(383)+mup1mup5*F(3106)+4*p3mup1*(
     -   p3mup5*F(4465)+p4mup5*F(4467)+p2mup5*F(4469))
       Fa16(1526) = 2*p2mup1*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(383))+mup1mup5*F(3108)
       Fa16(1527) = 2*p3mup1*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(383))-mup1mup5*F(3062)
       Fa16(1528) = 2*p1mup5*p3mup1*F(383)+mup1mup5*F(3061)+4*p3mup1*(
     -   p3mup5*F(4465)+p4mup5*F(4467)+p2mup5*F(4469))
       Fa16(1529) = 2*p2mup1*(2*(p1mup5*F(263)+p3mup5*F(273)+p2mup5*F(
     -   379))+p4mup5*F(383))+mup1mup5*F(3063)
       Fa16(1530) = mup1mup5*F(3015)+p3mup1*(4*(p1mup5*F(3312)+p2mup5*
     -   F(3313)+p3mup5*F(3476))+p4mup5*F(3525))
       Fa16(1531) = mup1mup5*F(3013)-2*p3mup1*(p1mup5*F(3482)+2*(p2mup
     -   5*F(3483)+p3mup5*F(3484)+p4mup5*F(3485)))
       Fa16(1532) = p4mup1*(p1mup4*F(3482)+2*(p2mup4*F(3483)+p3mup4*F(
     -   3484)+p4mup4*F(3485)))
       Fa16(1533) = -4*mup1mup4*(p1mup5*F(3312)+p2mup5*F(3313)+p3mup5*
     -   F(3476))+2*mup1mup4*p4mup5*F(3477)
       Fa16(1534) = mup1mup4*(p1mup5*F(3482)+2*(p2mup5*F(3483)+p3mup5*
     -   F(3484)+p4mup5*F(3485)))
       Fa16(1535) = mup1mup4*(4*(p1mup5*F(3002)+p4mup5*F(3003)+p2mup5*
     -   F(3382))+p3mup5*F(3488))
       Fa16(1536) = mup1mup5*F(3016)-4*p2mup1*(p1mup5*F(3312)+p2mup5*F
     -   (3313)+p3mup5*F(3476))+2*p2mup1*p4mup5*F(3477)
       Fa16(1537) = mup1mup5*F(3017)+4*(p1mup5*(p3mup1*F(2527)+p2mup1*
     -   F(3002))+p2mup1*p4mup5*F(3003)+p2mup1*p2mup5*F(3382))+p2mup1
     -   *p3mup5*F(3488)+2*p3mup1*(2*(p4mup5*F(2531)+p3mup5*F(2566))+
     -   p2mup5*F(3493))
       Fa16(1538) = mup1mup5*F(3014)+2*p2mup1*(p1mup5*F(3482)+2*(p2mup
     -   5*F(3483)+p3mup5*F(3484)+p4mup5*F(3485)))
       Fa16(1539) = 2*p3mup1*(p1mup4*F(3482)+2*(p2mup4*F(3483)+p3mup4*
     -   F(3484)+p4mup4*F(3485)))
       Fa16(1540) = p2mup1*(p1mup4*F(3482)+2*(p2mup4*F(3483)+p3mup4*F(
     -   3484)+p4mup4*F(3485)))
       Fa16(1541) = 2*(mup1mup4*(2*(p1mup5*F(2527)+p4mup5*F(2531)+p3mu
     -   p5*F(2566))+p2mup5*F(3493))-mup4mup5*p2mup1*F(3529))
       Fa16(1542) = mup1mup2*(p1mup4*F(3482)+2*(p2mup4*F(3483)+p3mup4*
     -   F(3484)+p4mup4*F(3485)))
       Fa16(1543) = -4*p1mup4*(p4mup1*F(3001)+p3mup1*F(3474))-2*p3mup1
     -   *(2*(p2mup4*F(3475)+p3mup4*F(3478))+p4mup4*F(3479))+p4mup1*(
     -   4*(p4mup4*F(3007)-p2mup4*F(3486))+p3mup4*F(3527))
       Fa16(1544) = mup1mup2*(2*(p1mup4*F(3474)+p2mup4*F(3475)+p3mup4*
     -   F(3478))+p4mup4*F(3479))
       Fa16(1545) = -4*(p1mup4*p3mup1*F(2525)-p3mup1*p4mup4*F(2535)+p3
     -   mup1*p3mup4*F(2565)-p1mup4*p2mup1*F(3001)+p2mup1*p4mup4*F(30
     -   07)-p2mup1*p2mup4*F(3486))+p2mup1*p3mup4*F(3487)-2*p2mup4*p3
     -   mup1*F(3492)
       Fa16(1546) = mup1mup2*(-4*(p1mup4*F(3001)-p4mup4*F(3007)+p2mup4
     -   *F(3486))+p3mup4*F(3527))
       Fa16(1547) = 2*mup1mup2*(2*(p1mup4*F(2525)-p4mup4*F(2535)+p3mup
     -   4*F(2565))+p2mup4*F(3492))
       Fa16(1548) = p4mup1*(4*(p1mup5*F(3002)+p4mup5*F(3003)+p2mup5*F(
     -   3382))+p3mup5*F(3488))
       Fa16(1549) = p4mup1*(p1mup5*F(3482)+2*(p2mup5*F(3483)+p3mup5*F(
     -   3484)+p4mup5*F(3485)))
       Fa16(1550) = 2*p4mup1*(2*(p1mup5*F(2527)+p4mup5*F(2531)+p3mup5*
     -   F(2566))+p2mup5*F(3493))
       Fa16(1551) = 2*p3mup1*(p1mup2*F(3482)+2*p3mup2*F(3484))
       Fa16(1552) = p2mup1*(p1mup2*F(3482)+2*p3mup2*F(3484))
       Fa16(1553) = mup1mup4*(p1mup2*F(3482)+2*p3mup2*F(3484))
       Fa16(1554) = p3mup1*(p1mup2*F(3310)-p3mup2*F(3480))
       Fa16(1555) = p2mup1*(p1mup2*F(3310)+p3mup2*F(3526))
       Fa16(1556) = mup1mup4*(p1mup2*F(3310)-p3mup2*F(3480))
       Fa16(1557) = -4*mup1mup4*p1mup2*F(3311)+2*(mup1mup4*p3mup2*F(34
     -   89)-mup2mup4*p3mup1*F(3529))
       Fa16(1558) = -4*mup1mup4*(p1mup2*F(3491)+p3mup2*F(3494))+2*mup2
     -   mup4*p2mup1*F(3529)
       Fa16(1559) = -2*p2mup1*p3mup2*F(3489)+4*(p1mup2*(p2mup1*F(3311)
     -   +p3mup1*F(3491))+p3mup1*p3mup2*F(3494))
       Fa16(1560) = p3mup1*(p1mup2*F(383)+2*p3mup2*F(4465))
       Fa16(1561) = p3mup1*(p1mup2*F(1262)+p3mup2*F(2368))
       Fa16(1562) = p2mup1*(p1mup2*F(1608)+p3mup2*F(3506))
       Fa16(1563) = mup1mup4*(p1mup2*F(1262)+p3mup2*F(2368))
       Fa16(1564) = mup2mup4*p3mup1-mup1mup4*p3mup2
       Fa16(1565) = p2mup1*p3mup2*F(4464)-2*p3mup1*(p1mup2*F(357)+p3mu
     -   p2*F(4465))
       Fa16(1566) = p3mup1*(p1mup2*F(443)+2*p3mup2*F(445))
       Fa16(1567) = p3mup1*(p1mup2*F(1608)+p3mup2*F(3506))
       Fa16(1568) = p2mup1*(p1mup2*F(1262)+p3mup2*F(2368))
       Fa16(1569) = mup1mup4*(p1mup2*F(1608)+p3mup2*F(3506))
       Fa16(1570) = -(mup2mup4*p3mup1)+mup1mup4*p3mup2
       Fa16(1571) = -(p2mup1*p3mup2*F(477))+2*p3mup1*(p3mup2*F(445)+p1
     -   mup2*F(1629))
       Fa16(1572) = p2mup1*p3mup2*F(477)-2*p3mup1*(p3mup2*F(445)+p1mup
     -   2*F(1629))
       Fa16(1573) = -(p2mup1*p3mup2*F(4464))+2*(p1mup2*p3mup1*F(357)+p
     -   3mup1*p3mup2*F(4465))
       Fa16(1574) = 2*p4mup1*(2*(p1mup5*F(3305)+p4mup5*F(3306)+p2mup5*
     -   F(3314))+p3mup5*F(3528))
       Fa16(1575) = 2*p4mup1*(p1mup5*F(3482)+2*(p2mup5*F(3483)+p3mup5*
     -   F(3484)+p4mup5*F(3485)))
       Fa16(1576) = p4mup1*(2*(p1mup5*F(2527)+p4mup5*F(2531)+p3mup5*F(
     -   2566))+p2mup5*F(3493))
       Fa16(1577) = p3mup1*(p1mup2*F(3482)+2*p3mup2*F(3484))
       Fa16(1578) = 2*p2mup1*(p1mup2*F(3482)+2*p3mup2*F(3484))
       Fa16(1579) = 2*mup1mup4*(p1mup2*F(3482)+2*p3mup2*F(3484))
       Fa16(1580) = p3mup1*(p1mup2*F(3310)+p3mup2*F(3526))
       Fa16(1581) = p2mup1*(p1mup2*F(3310)-p3mup2*F(3480))
       Fa16(1582) = mup1mup4*(p1mup2*F(3310)+p3mup2*F(3526))
       Fa16(1583) = 4*mup1mup4*p1mup2*F(3311)-2*(mup1mup4*p3mup2*F(348
     -   9)-mup2mup4*p3mup1*F(3529))
       Fa16(1584) = 4*mup1mup4*(p1mup2*F(3491)+p3mup2*F(3494))-2*mup2m
     -   up4*p2mup1*F(3529)
       Return
       End
