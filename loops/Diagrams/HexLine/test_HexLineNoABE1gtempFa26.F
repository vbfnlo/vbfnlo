c************** Calling the Fa functions*************************
       subroutine test_FaFunctionhexlineNoAbe1g26(p1mup2,p1mup3,p1mup4
     -   ,p1mup6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3m
     -   up6,p4mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,
     -   p6mup2,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3m
     -   up4,mup3mup6,mup4mup6,Fa26)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6
        common/test_FhlFunctionsNoAbe1g/F
       COMMON/Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45
     -   ,s56,s16,s123,s234,s345
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
      Complex*16 Fa(587),F(6395)
      Real*8 P(7129) 
        COMMON/test_PFunctionsNoAbe1g/P
       Complex*16 Fa26(351:364)
       COMMON/test_FaFunctionsNoAbe1g/Fa
       Fa26(351) = p4mup3*F(1466)+p5mup3*F(1467)+(p1mup3+p2mup3)*F(146
     -   8)+p3mup3*F(1469)
       Fa26(352) = mup4mup6*((p1mup3+p2mup3)*F(1470)+p3mup3*F(1471)+p4
     -   mup3*F(1472)+p5mup3*F(1473))+mup3mup6*p5mup4*F(1474)
       Fa26(353) = mup3mup6*(p3mup4*F(1475)+p4mup4*F(1476))
       Fa26(354) = p5mup6*(8*(p2mup3*p2mup4*F(1297)+p2mup4*p3mup3*F(12
     -   98)+p2mup3*p3mup4*F(1304)+p3mup3*p3mup4*F(1305)-p2mup3*p4mup
     -   4*F(1310)-p3mup3*p4mup4*F(1311)+p3mup4*p4mup3*F(1317)-p3mup3
     -   *p5mup4*F(1322)+p3mup4*p5mup3*F(1323)+p4mup3*p4mup4*F(1325)-
     -   p4mup4*p5mup3*F(1327)-p4mup3*p5mup4*F(1332)-p5mup3*p5mup4*F(
     -   1337)+p2mup4*p4mup3*F(1348)+p2mup4*p5mup3*F(1355)+p1mup4*(p2
     -   mup3*F(1297)+p3mup3*F(1298)+p4mup3*F(1348)+p5mup3*F(1355))-p
     -   2mup3*p5mup4*F(1358)+p1mup3*((p1mup4+p2mup4)*F(1297)+p3mup4*
     -   F(1304)-p4mup4*F(1310)-p5mup4*F(1358)))-mup3mup4*F(1477))
       Fa26(355) = mup3mup6*(p1mup4+p2mup4)*F(1478)-8*(p3mup3*p4mup6*p
     -   5mup4*F(1319)-p5mup4*(p4mup3*p4mup6*F(1331)+p4mup6*p5mup3*F(
     -   1336)+(p1mup3+p2mup3)*p4mup6*F(1356)+(p1mup6+p2mup6)*p4mup3*
     -   F(1360)+p2mup6*p5mup3*F(1479))-p3mup6*p5mup4*(p4mup3*F(1330)
     -   +p5mup3*F(1480)))
       Fa26(356) = 8*(p2mup3*p2mup4*p2mup6*F(1292)-p2mup3*p2mup6*p3mup
     -   4*F(1293)+p2mup4*p2mup6*p3mup3*F(1294)+p2mup3*p2mup4*p3mup6*
     -   F(1295)+p2mup4*p3mup3*p3mup6*F(1296)+p2mup3*p2mup4*p4mup6*F(
     -   1299)+p2mup4*p3mup3*p4mup6*F(1300)-p2mup6*p3mup3*p3mup4*F(13
     -   01)+p2mup3*p3mup4*p3mup6*F(1302)+p3mup3*p3mup4*p3mup6*F(1303
     -   )+p2mup3*p3mup4*p4mup6*F(1306)+p3mup3*p3mup4*p4mup6*F(1307)+
     -   p2mup4*p2mup6*p4mup3*F(1308)+p2mup3*p2mup6*p4mup4*F(1308)+p1
     -   mup3*p1mup6*((p1mup4+p2mup4)*F(1292)-p3mup4*F(1293)+p4mup4*F
     -   (1308))+p2mup3*p3mup6*p4mup4*F(1309)+p2mup3*p4mup4*p4mup6*F(
     -   1312)+p2mup6*p3mup3*p4mup4*F(1313)+p3mup3*p4mup4*p4mup6*F(13
     -   15)+p3mup4*p4mup6*p5mup3*F(1320)+p3mup6*p4mup4*p5mup3*F(1324
     -   )+p4mup4*p4mup6*p5mup3*F(1326)+p3mup6*p4mup3*p4mup4*F(1328)+
     -   p2mup4*p3mup6*p4mup3*F(1338)+p2mup6*p3mup4*p4mup3*F(1339)+p2
     -   mup4*p2mup6*p5mup3*F(1341)+p2mup4*p3mup6*p5mup3*F(1344)+p3mu
     -   p4*p3mup6*p4mup3*F(1349)+p2mup6*p3mup4*p5mup3*F(1352)+p2mup4
     -   *p4mup6*p5mup3*F(1354)+p1mup4*(p2mup3*p2mup6*F(1292)+p2mup6*
     -   p3mup3*F(1294)+p2mup3*p3mup6*F(1295)+p3mup3*p3mup6*F(1296)+p
     -   2mup3*p4mup6*F(1299)+p3mup3*p4mup6*F(1300)+p2mup6*p4mup3*F(1
     -   308)+p3mup6*p4mup3*F(1338)+p2mup6*p5mup3*F(1341)+p3mup6*p5mu
     -   p3*F(1344)+p4mup6*p5mup3*F(1354))+p2mup6*p4mup4*p5mup3*F(135
     -   7)+p2mup6*p4mup3*p4mup4*F(1359)+p3mup4*p3mup6*p5mup3*F(1365)
     -   +p1mup6*p5mup3*p5mup4*F(1479))-mup3mup4*p3mup6*F(1481)-mup3m
     -   up4*p4mup6*F(1482)+mup3mup4*(p1mup6+p2mup6)*F(1483)
       Fa26(357) = p1mup3*p2mup4*p2mup6*F(1292)-p1mup3*p2mup6*p3mup4*F
     -   (1293)+p1mup3*p2mup4*p3mup6*F(1295)+p1mup3*p2mup4*p4mup6*F(1
     -   299)+p1mup3*p3mup4*p3mup6*F(1302)+p1mup3*p3mup4*p4mup6*F(130
     -   6)+p1mup3*p2mup6*p4mup4*F(1308)+p1mup3*p3mup6*p4mup4*F(1309)
     -   +p2mup4*p4mup3*p4mup6*F(1312)+p1mup3*p4mup4*p4mup6*F(1312)+p
     -   3mup3*p3mup6*p4mup4*F(1314)+p3mup4*p4mup3*p4mup6*F(1316)+p4m
     -   up3*p4mup4*p4mup6*F(1329)+p1mup4*(p1mup3*p2mup6*F(1292)+p1mu
     -   p3*p3mup6*F(1295)+p1mup3*p4mup6*F(1299)+p4mup3*p4mup6*F(1312
     -   )+p1mup6*(p2mup3*F(1292)+p3mup3*F(1294)+p4mup3*F(1308)+p5mup
     -   3*F(1341)))+p1mup6*(p2mup3*p2mup4*F(1292)-p2mup3*p3mup4*F(12
     -   93)+p2mup4*p3mup3*F(1294)-p3mup3*p3mup4*F(1301)+p2mup4*p4mup
     -   3*F(1308)+p2mup3*p4mup4*F(1308)+p3mup3*p4mup4*F(1313)+p3mup4
     -   *p4mup3*F(1339)+p2mup4*p5mup3*F(1341)+p3mup4*p5mup3*F(1352)+
     -   p4mup4*p5mup3*F(1357)+p4mup3*p4mup4*F(1359))
       Fa26(358) = p5mup4*(p1mup6*p2mup3*F(1345)+p2mup3*p2mup6*F(1345)
     -   +p1mup3*(p1mup6+p2mup6)*F(1345)+p1mup3*p3mup6*F(1346)+p2mup3
     -   *p3mup6*F(1346)-p1mup6*p3mup3*F(1351)-p2mup6*p3mup3*F(1351)-
     -   p3mup3*p3mup6*F(1364))
       Fa26(359) = 2*(p1mup3+p2mup3)*F(1409)+p5mup3*F(1484)+p4mup3*F(1
     -   485)+p3mup3*F(1486)
       Fa26(360) = p4mup3*p5mup4*F(1487)-p3mup4*p5mup3*F(1488)-p4mup4*
     -   p5mup3*F(1489)-p1mup4*p4mup3*F(1490)-p2mup4*p4mup3*F(1490)-p
     -   3mup4*p4mup3*F(1491)+4*(p3mup3*p4mup4*F(1426)+p3mup3*p3mup4*
     -   F(1427)+p1mup3*(p4mup4*F(1425)+p3mup4*F(1430)+(p1mup4+p2mup4
     -   )*F(1436))+p2mup3*(p4mup4*F(1425)+p3mup4*F(1430)+(p1mup4+p2m
     -   up4)*F(1436))+p1mup4*p3mup3*F(1446)+p2mup4*p3mup3*F(1446)+p4
     -   mup3*p4mup4*F(1492))-p5mup3*p5mup4*F(1493)+p3mup3*p5mup4*F(1
     -   494)-(p1mup4+p2mup4)*p5mup3*F(1495)-mup3mup4*F(1496)+p1mup3*
     -   p5mup4*F(1497)+p2mup3*p5mup4*F(1497)
       Fa26(361) = 4*p3mup2*(-(mup4mup6*F(1498))-2*(p2mup4*p2mup6*F(14
     -   99)+p2mup4*p4mup6*F(1500)+p2mup4*p5mup6*F(1501)+p2mup4*p3mup
     -   6*F(1502)+p1mup4*((p1mup6+p2mup6)*F(1499)+p4mup6*F(1500)+p5m
     -   up6*F(1501)+p3mup6*F(1502))-p2mup6*p3mup4*F(1503)-p2mup6*p4m
     -   up4*F(1504)-p3mup4*p4mup6*F(1505)-p4mup4*p4mup6*F(1506)-p3mu
     -   p4*p5mup6*F(1507)-p4mup4*p5mup6*F(1508)-p3mup4*p3mup6*F(1509
     -   )-p3mup6*p4mup4*F(1510)-p3mup6*p5mup4*F(1511)-p2mup6*p5mup4*
     -   F(1512)+p1mup6*(p2mup4*F(1499)-p3mup4*F(1503)-p4mup4*F(1504)
     -   -p5mup4*F(1512))-p4mup6*p5mup4*F(1513)-p5mup4*p5mup6*F(1514)
     -   ))
       Fa26(362) = (p1mup2+p2mup2)*(mup4mup6*F(1498)+2*(p2mup4*p2mup6*
     -   F(1499)+p2mup4*p4mup6*F(1500)+p2mup4*p5mup6*F(1501)+p2mup4*p
     -   3mup6*F(1502)+p1mup4*((p1mup6+p2mup6)*F(1499)+p4mup6*F(1500)
     -   +p5mup6*F(1501)+p3mup6*F(1502))-p2mup6*p3mup4*F(1503)-p2mup6
     -   *p4mup4*F(1504)-p3mup4*p4mup6*F(1505)-p4mup4*p4mup6*F(1506)-
     -   p3mup4*p5mup6*F(1507)-p4mup4*p5mup6*F(1508)-p3mup4*p3mup6*F(
     -   1509)-p3mup6*p4mup4*F(1510)-p3mup6*p5mup4*F(1511)-p2mup6*p5m
     -   up4*F(1512)+p1mup6*(p2mup4*F(1499)-p3mup4*F(1503)-p4mup4*F(1
     -   504)-p5mup4*F(1512))-p4mup6*p5mup4*F(1513)-p5mup4*p5mup6*F(1
     -   514)))
       Fa26(363) = p4mup2*(p4mup4*F(1519)+p3mup4*F(1520)+p5mup4*F(1521
     -   )+(p1mup4+p2mup4)*F(1522))+p3mup2*(p4mup4*F(1517)+p3mup4*F(1
     -   518)+p5mup4*F(1523)+(p1mup4+p2mup4)*F(1524))
       Fa26(364) = p1mup2*(p3mup4*F(1526)+p5mup4*F(1527)+4*(p4mup4*F(1
     -   525)+(p1mup4+p2mup4)*F(1528)))+p2mup2*(p3mup4*F(1526)+p5mup4
     -   *F(1527)+4*(p4mup4*F(1525)+(p1mup4+p2mup4)*F(1528)))-4*p4mup
     -   2*s12*(p4mup4*F(1529)+p3mup4*F(1530)+p5mup4*F(1531)+(p1mup4+
     -   p2mup4)*F(1532))
       Return
       End
