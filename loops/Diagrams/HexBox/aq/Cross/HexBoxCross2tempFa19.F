c************** Calling the Fa functions*************************
       subroutine HexBoxCross2Fa19(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,
     &   p3mup6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa
     &   19)
       IMPLICIT NONE
      Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
        common/HexBoxCross2FhlFunctions/F
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
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
      Complex*16 Fa(392),F(1853)
      Real*8 P(1355) 
        COMMON/HexBoxCross2PFunctions/P
       Complex*16 Fa19(343:361)
       COMMON/HexBoxCross2FaFunctions/Fa
       Fa19(343) = p1mup3*F(1555)-p2mup3*F(1556)-p4mup3*F(1557)+p5mup3
     &   *F(1558)+p6mup3*F(1559)
       Fa19(344) = p1mup6*F(891)+p4mup6*F(1561)+p5mup6*F(1562)+p6mup6*
     &   F(1563)+p2mup6*F(1564)
       Fa19(345) = -(p4mup6*F(1565))-p5mup6*F(1566)+p6mup6*F(1567)+p2m
     &   up6*F(1568)+p1mup6*F(1569)
       Fa19(346) = p1mup6*F(1570)+p5mup6*F(1571)+p6mup6*F(1572)-p2mup6
     &   *F(1573)+p4mup6*F(1574)
       Fa19(347) = p1mup6*F(1575)+p4mup6*F(1576)+p5mup6*F(1577)+p6mup6
     &   *F(1578)-p2mup6*F(1579)
       Fa19(348) = p1mup6*F(1580)+p5mup6*F(1581)+p6mup6*F(1582)-p2mup6
     &   *F(1583)+p4mup6*F(1584)
       Fa19(349) = -(p1mup3*F(1585))-p2mup3*F(1586)+p4mup3*F(1587)+p5m
     &   up3*F(1588)+p6mup3*F(1589)
       Fa19(350) = p1mup3*F(887)+p2mup3*F(932)+p4mup3*F(1590)+p5mup3*F
     &   (1591)+p6mup3*F(1592)
       Fa19(351) = p2mup3*F(1460)+p1mup3*F(1593)-p4mup3*F(1594)-p5mup3
     &   *F(1595)-p6mup3*F(1596)
       Fa19(352) = p1mup3*F(1597)-p4mup3*F(1598)-p5mup3*F(1599)-p6mup3
     &   *F(1600)+p2mup3*F(1601)
       Fa19(353) = p5mup3*F(949)+p2mup3*F(1464)+p1mup3*F(1602)-p4mup3*
     &   F(1603)+p6mup3*F(1604)
       Fa19(354) = p6mup3*F(928)-p2mup3*F(1468)-p1mup3*F(1605)+p4mup3*
     &   F(1606)+p5mup3*F(1607)
       Fa19(355) = p2mup3*F(1519)-p4mup3*F(1557)+p6mup3*F(1559)+p1mup3
     &   *F(1608)+p5mup3*F(1609)
       Fa19(356) = -(p4mup3*p5mup6*F(1285))+p1mup6*p4mup3*F(1289)-p2mu
     &   p3*p2mup6*F(1371)+p2mup3*p5mup6*F(1610)+p2mup3*p6mup6*F(1611
     &   )+p1mup3*p5mup6*F(1612)+mup3mup6*F(1613)-p1mup3*p6mup6*F(161
     &   4)-p1mup6*p2mup3*F(1615)-p1mup3*p2mup6*F(1615)+p1mup3*p4mup6
     &   *F(1616)-p1mup3*p1mup6*F(1617)-p4mup3*p6mup6*F(1619)+p2mup3*
     &   p4mup6*F(1620)+p2mup6*p4mup3*F(1624)+p5mup3*(-(p5mup6*F(1194
     &   ))+p1mup6*F(1291)-p6mup6*F(1618)+p2mup6*F(1622)+p4mup6*F(162
     &   5))+p6mup3*(-(p6mup6*F(1031))-p5mup6*F(1618)+p1mup6*F(1621)+
     &   p2mup6*F(1623)+p4mup6*F(1626))+p4mup3*p4mup6*F(1627)
       Fa19(357) = -(p6mup6*F(1629))+p2mup6*F(1630)+p4mup6*F(1631)+p1m
     &   up6*F(1632)-p5mup6*F(1633)
       Fa19(358) = p1mup3*F(896)-p2mup3*F(1634)+p4mup3*F(1635)-p5mup3*
     &   F(1636)+p6mup3*F(1637)
       Fa19(359) = p1mup3*F(887)-p6mup3*F(1323)+p4mup3*F(1590)+p5mup3*
     &   F(1591)-p2mup3*F(1638)
       Fa19(360) = -(p1mup3*F(1597))+p4mup3*F(1598)+p5mup3*F(1599)-p2m
     &   up3*F(1601)+p6mup3*F(1639)
       Fa19(361) = 4*(p6mup3*F(957)+p4mup3*F(1558)+p1mup3*F(1640)-p2mu
     &   p3*F(1641))+p5mup3*F(1642)
       Return
       End
