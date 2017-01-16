c************** Calling the Fa functions*************************
       subroutine HexBox2Fa20(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa20)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBox2FhlFunctions/F
       COMMON/HexBox2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,
     &   s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(385),F(1646)
      Real*8 P(1185) 
        COMMON/HexBox2PFunctions/P
       Complex*16 Fa20(362:385)
       COMMON/HexBox2FaFunctions/Fa
       Fa20(362) = p1mup6*F(1497)-p2mup6*F(1498)+p4mup6*F(1499)+p5mup6
     &   *F(1500)+p6mup6*F(1501)
       Fa20(363) = 4*(p4mup6*F(1263)+p5mup6*F(1391)+p6mup6*F(1392)+p1m
     &   up6*F(1394))+p2mup6*F(1502)
       Fa20(364) = 4*(p2mup6*F(1263)+p5mup6*F(1265)+p6mup6*F(1266)+p1m
     &   up6*F(1267))+p4mup6*F(1503)
       Fa20(365) = p4mup6*F(1270)+p5mup6*F(1271)-p1mup6*F(1505)-p2mup6
     &   *F(1506)+p6mup6*F(1644)
       Fa20(366) = -(p1mup6*F(1509))-p2mup6*F(1510)+p4mup6*F(1511)+p5m
     &   up6*F(1512)+p6mup6*F(1513)
       Fa20(367) = p1mup6*F(1514)+p2mup6*F(1515)+4*(p4mup6*F(1516)+p5m
     &   up6*F(1517)+p6mup6*F(1518))
       Fa20(368) = p1mup6*p4mup3*F(1108)+p1mup6*p2mup3*F(1174)+p1mup6*
     &   p6mup3*F(1175)+p5mup3*(p6mup6*F(915)+p2mup6*F(1037)+p4mup6*F
     &   (1038)+p1mup6*F(1176))-p1mup3*p4mup6*F(1251)+p1mup3*p2mup6*F
     &   (1383)+p1mup3*p1mup6*F(1521)+p5mup6*(p6mup3*F(915)+p2mup3*F(
     &   1037)+p4mup3*F(1038)+p5mup3*F(1039)-p1mup3*F(1522))-p1mup3*p
     &   6mup6*F(1523)+p2mup3*p2mup6*F(1524)+p2mup6*p6mup3*F(1525)+p2
     &   mup3*p6mup6*F(1525)+p2mup6*p4mup3*F(1526)+p2mup3*p4mup6*F(15
     &   26)+p4mup6*p6mup3*F(1527)+p4mup3*p6mup6*F(1527)+p4mup3*p4mup
     &   6*F(1528)+p6mup3*p6mup6*F(1645)+mup3mup6*F(1646)
       Fa20(369) = mup3mup6*F(1037)+p1mup3*(p2mup6*F(1529)+p5mup6*F(15
     &   30)+p6mup6*F(1531)+p1mup6*F(1532)+p4mup6*F(1533))+p1mup6*p5m
     &   up3*F(1534)+p1mup6*p6mup3*F(1535)+p1mup6*p4mup3*F(1536)+p2mu
     &   p6*p4mup3*F(1537)+p2mup6*p5mup3*F(1538)+p2mup6*p6mup3*F(1539
     &   )-p4mup6*p5mup3*F(1540)-p4mup3*p5mup6*F(1540)-p4mup6*p6mup3*
     &   F(1541)-p4mup3*p6mup6*F(1541)-p5mup3*p5mup6*F(1542)-p5mup6*p
     &   6mup3*F(1543)-p5mup3*p6mup6*F(1543)-p6mup3*p6mup6*F(1544)-p4
     &   mup3*p4mup6*F(1545)+p2mup3*(p4mup6*F(943)+p5mup6*F(944)+p6mu
     &   p6*F(945)-p1mup6*F(1386)+p2mup6*F(1546))
       Fa20(370) = -(p4mup3*p5mup6*F(1071))-p5mup3*p5mup6*F(1072)-p5mu
     &   p6*p6mup3*F(1073)-p5mup3*p6mup6*F(1073)+p1mup6*p4mup3*F(1178
     &   )+p1mup6*p5mup3*F(1179)+p1mup6*p6mup3*F(1180)-p4mup3*p6mup6*
     &   F(1210)+p1mup3*p2mup6*F(1277)+p1mup3*p5mup6*F(1278)+p1mup3*p
     &   6mup6*F(1279)-p4mup6*(p5mup3*F(1071)+p6mup3*F(1210)+p4mup3*F
     &   (1211)-p1mup3*F(1280))+p1mup3*p1mup6*F(1281)-p6mup3*p6mup6*F
     &   (1359)+mup3mup6*F(1525)+p2mup6*p4mup3*F(1547)+p2mup6*p5mup3*
     &   F(1548)+p2mup6*p6mup3*F(1549)+p2mup3*(p5mup6*F(945)-p1mup6*F
     &   (1191)+p4mup6*F(1360)+p6mup6*F(1361)+p2mup6*F(1550))
       Fa20(371) = p2mup3*(p6mup6*F(911)+p1mup6*F(942)+p2mup6*F(943)+p
     &   5mup6*F(946)+p4mup6*F(948))+mup3mup6*F(1038)+p1mup3*p2mup6*F
     &   (1533)-p2mup6*p5mup3*F(1540)-p2mup6*p6mup3*F(1541)-p2mup6*p4
     &   mup3*F(1545)+p1mup3*p5mup6*F(1551)+p1mup3*p6mup6*F(1552)+p1m
     &   up3*p1mup6*F(1553)+p1mup3*p4mup6*F(1554)-p1mup6*p5mup3*F(155
     &   5)-p1mup6*p6mup3*F(1556)-p1mup6*p4mup3*F(1557)-p4mup6*p5mup3
     &   *F(1558)-p4mup3*p5mup6*F(1558)+p4mup6*p6mup3*F(1559)+p4mup3*
     &   p6mup6*F(1559)-p5mup3*p5mup6*F(1560)-p5mup6*p6mup3*F(1561)-p
     &   5mup3*p6mup6*F(1561)-p6mup3*p6mup6*F(1562)-p4mup3*p4mup6*F(1
     &   563)
       Fa20(372) = -(p5mup6*p6mup3*F(904))-p4mup3*p5mup6*F(905)+p2mup3
     &   *p5mup6*F(911)-p5mup3*(p5mup6*F(903)+p6mup6*F(904)+p4mup6*F(
     &   905)+p1mup6*F(1070)+p2mup6*F(1071))-p1mup6*p6mup3*F(1181)-p1
     &   mup6*p4mup3*F(1186)+p1mup6*p2mup3*F(1193)-p2mup6*p6mup3*F(12
     &   10)-p2mup6*p4mup3*F(1211)+p2mup3*p2mup6*F(1360)+mup3mup6*F(1
     &   527)+p1mup3*(p5mup6*F(896)+p2mup6*F(1280)+p1mup6*F(1282)+p6m
     &   up6*F(1564)+p4mup6*F(1565))+p4mup6*p6mup3*F(1566)+p4mup3*p6m
     &   up6*F(1566)-p6mup3*p6mup6*F(1567)-p4mup3*p4mup6*F(1568)+p2mu
     &   p3*p6mup6*F(1569)+p2mup3*p4mup6*F(1570)
       Fa20(373) = p2mup3*F(1515)+4*(p4mup3*F(1516)+p5mup3*F(1517)+p6m
     &   up3*F(1518)+p1mup3*F(1571))
       Fa20(374) = p1mup6*p2mup3*F(941)-p2mup3*p4mup6*F(1071)-p2mup3*p
     &   5mup6*F(1072)-p2mup3*p6mup6*F(1073)+p2mup3*p2mup6*F(1548)-p1
     &   mup6*p4mup3*F(1556)+p1mup3*p2mup6*F(1572)-p1mup3*p5mup6*F(15
     &   73)-p1mup3*p6mup6*F(1574)-p1mup3*p4mup6*F(1575)-p2mup6*p4mup
     &   3*F(1576)-p2mup6*p5mup3*F(1577)-p2mup6*p6mup3*F(1578)-p4mup6
     &   *p5mup3*F(1579)-p4mup3*p5mup6*F(1579)-p4mup6*p6mup3*F(1580)-
     &   p4mup3*p6mup6*F(1580)-p4mup3*p4mup6*F(1581)-p5mup3*p5mup6*F(
     &   1582)+p5mup6*p6mup3*F(1583)+p5mup3*p6mup6*F(1583)+p6mup3*p6m
     &   up6*F(1584)-p1mup3*p1mup6*F(1585)-p1mup6*p5mup3*F(1586)-p1mu
     &   p6*p6mup3*F(1587)+mup3mup6*F(1588)
       Fa20(375) = p1mup6*(p5mup3*F(1069)+p4mup3*F(1181)+p6mup3*F(1182
     &   )-p2mup3*F(1192)-p1mup3*F(1285))+p2mup3*(p5mup6*F(1073)+p4mu
     &   p6*F(1210)+p6mup6*F(1359)-p2mup6*F(1549))-p1mup3*p2mup6*F(15
     &   89)+p1mup3*p5mup6*F(1590)+p1mup3*p6mup6*F(1591)+p1mup3*p4mup
     &   6*F(1592)+p2mup6*p4mup3*F(1593)+p2mup6*p5mup3*F(1594)+p2mup6
     &   *p6mup3*F(1595)+p4mup6*p5mup3*F(1596)+p4mup3*p5mup6*F(1596)+
     &   p4mup6*p6mup3*F(1597)+p4mup3*p6mup6*F(1597)-p4mup3*p4mup6*F(
     &   1598)-p5mup3*p5mup6*F(1599)+p5mup6*p6mup3*F(1600)+p5mup3*p6m
     &   up6*F(1600)-p6mup3*p6mup6*F(1601)-mup3mup6*F(1602)
       Fa20(376) = p1mup3*F(1460)-p6mup3*F(1461)-p4mup3*F(1462)-p5mup3
     &   *F(1463)+p2mup3*F(1464)
       Fa20(377) = p1mup3*F(1603)-p6mup3*F(1604)-p4mup3*F(1605)-p5mup3
     &   *F(1606)+p2mup3*F(1607)
       Fa20(378) = p2mup3*F(991)+p1mup3*F(1159)+p6mup3*F(1160)+p5mup3*
     &   F(1162)-p4mup3*F(1608)
       Fa20(379) = p5mup3*F(1137)+p1mup3*F(1609)+p6mup3*F(1610)-p4mup3
     &   *F(1611)+p2mup3*F(1612)
       Fa20(380) = -(p1mup3*F(1269))+p4mup3*F(1270)+p5mup3*F(1271)-p2m
     &   up3*F(1506)+p6mup3*F(1644)
       Fa20(381) = -(p2mup6*F(1614))+p4mup6*F(1615)+p5mup6*F(1616)+p6m
     &   up6*F(1617)-p1mup6*F(1618)
       Fa20(382) = -(p1mup6*F(1622))-p2mup6*F(1623)+p4mup6*F(1624)+p5m
     &   up6*F(1625)+p6mup6*F(1626)
       Fa20(383) = p5mup3*F(1037)+p1mup3*F(1383)+p2mup3*F(1524)+p6mup3
     &   *F(1525)+p4mup3*F(1526)
       Fa20(384) = -(p2mup6*F(1624))-p1mup6*F(1627)+p4mup6*F(1628)+p5m
     &   up6*F(1629)+p6mup6*F(1630)
       Fa20(385) = p5mup3*F(1038)-p1mup3*F(1251)+p2mup3*F(1526)+p6mup3
     &   *F(1527)+p4mup3*F(1528)
       Return
       End
