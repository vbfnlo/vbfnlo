c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenGSim5(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4
     -   ,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5m
     -   up3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mu
     -   p3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa5)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenEvenGSim/F
      Complex*16 Fa(123),F(3806)
      Real*8 P(1320) 
        COMMON/PPenEvenGSim/P
       Complex*16 Fa5(33:40)
       COMMON/FaPenEvenGSim/Fa
       Fa5(33) = -(mup2mup3*p2mup5*p3mup1*p3mup4*F(1309))-mup2mup3*p2m
     -   up4*p3mup1*p3mup5*F(1309)+mup3mup5*p2mup1*p2mup4*p3mup2*F(14
     -   48)-mup2mup3*p2mup1*p2mup5*p3mup4*F(1449)-mup2mup3*p2mup1*p2
     -   mup4*p3mup5*F(1449)+mup3mup4*p2mup1*p2mup5*p3mup2*F(1571)-mu
     -   p4mup5*p2mup2*p3mup1*p3mup3*F(1575)+mup3mup5*p2mup2*p3mup1*p
     -   3mup4*F(1575)+mup3mup4*p2mup2*p3mup1*p3mup5*F(1575)-mup4mup5
     -   *p2mup3*p3mup1*p3mup2*F(1576)+mup3mup5*p2mup4*p3mup1*p3mup2*
     -   F(1576)+mup3mup4*p2mup5*p3mup1*p3mup2*F(1577)-mup4mup5*p2mup
     -   2*p3mup3*p4mup1*F(1581)-mup4mup5*p2mup1*p2mup2*p4mup3*F(1590
     -   )+mup3mup5*p2mup1*p2mup2*p4mup4*F(1590)+mup3mup4*p2mup1*p2mu
     -   p2*p4mup5*F(1590)-mup4mup5*p2mup1*p2mup3*p4mup2*F(1591)+mup3
     -   mup5*p2mup1*p2mup4*p4mup2*F(1591)-mup2mup3*p2mup1*p2mup5*p4m
     -   up4*F(1592)-mup2mup3*p2mup1*p2mup4*p4mup5*F(1592)-mup2mup5*(
     -   p2mup4*p3mup1*p3mup3*F(1309)-p2mup3*p3mup1*p3mup4*F(1309)+p2
     -   mup1*p2mup4*p3mup3*F(1449)+p2mup1*(p2mup4*p4mup3-p2mup3*p4mu
     -   p4)*F(1592))+mup2mup4*(p2mup5*p3mup1*p3mup3*F(1309)+p2mup3*p
     -   3mup1*p3mup5*F(1309)+p2mup1*p2mup5*p3mup3*F(1449)+2*(p2mup3*
     -   p2mup5*p3mup1*F(1574)+p2mup3*p2mup5*p4mup1*F(1580))+p2mup1*(
     -   p2mup5*p4mup3+p2mup3*p4mup5)*F(1592))+mup3mup4*p2mup1*p2mup5
     -   *p4mup2*F(1593)-mup4mup5*p2mup2*p3mup1*p4mup3*F(1594)+mup3mu
     -   p5*p2mup2*p3mup1*p4mup4*F(1594)+mup3mup4*p2mup2*p3mup1*p4mup
     -   5*F(1594)-2*(mup4mup5*p2mup2*p2mup3*p3mup1*F(1573)-mup3mup5*
     -   p2mup2*p2mup4*p3mup1*F(1573)+mup2mup3*p2mup4*p2mup5*p3mup1*F
     -   (1574)+mup4mup5*p2mup2*p2mup3*p4mup1*F(1579)-mup3mup5*p2mup2
     -   *p2mup4*p4mup1*F(1579)+mup2mup3*p2mup4*p2mup5*p4mup1*F(1580)
     -   -mup3mup4*p2mup2*p2mup5*p3mup1*F(1583)+mup4mup5*p2mup1*p3mup
     -   2*p3mup3*F(1584)-mup3mup5*p2mup1*p3mup2*p3mup4*F(1584)-mup3m
     -   up4*p2mup1*p3mup2*p3mup5*F(1584)-mup3mup4*p2mup2*p2mup5*p4mu
     -   p1*F(1587)+mup4mup5*p2mup1*p3mup3*p4mup2*F(1588)-mup3mup5*p2
     -   mup1*p3mup4*p4mup2*F(1588)-mup3mup4*p2mup1*p3mup5*p4mup2*F(1
     -   588)+p2mup1*p3mup2*(mup4mup5*p4mup3-mup3mup5*p4mup4-mup3mup4
     -   *p4mup5)*F(1595))
       Fa5(34) = mup3mup5*p2mup2*p3mup4*p4mup1*F(1598)+mup1mup3*p2mup2
     -   *p3mup4*p4mup5*F(1600)+mup3mup4*p2mup2*p3mup5*p4mup1*F(1601)
     -   +mup1mup3*p2mup2*p3mup5*p4mup4*F(1603)+mup4mup5*p2mup3*p3mup
     -   1*p4mup2*F(1604)+mup2mup5*p2mup3*p3mup1*p4mup4*F(1605)+mup2m
     -   up4*p2mup3*p3mup1*p4mup5*F(1606)+mup4mup5*p2mup3*p3mup2*p4mu
     -   p1*F(1607)+mup2mup5*p2mup3*p3mup4*p4mup1*F(1610)+mup1mup2*p2
     -   mup3*p3mup4*p4mup5*F(1612)+mup2mup4*p2mup3*p3mup5*p4mup1*F(1
     -   613)+mup1mup2*p2mup3*p3mup5*p4mup4*F(1615)+mup3mup5*p2mup4*p
     -   3mup1*p4mup2*F(1616)+mup2mup5*p2mup4*p3mup1*p4mup3*F(1617)+m
     -   up2mup3*p2mup4*p3mup1*p4mup5*F(1618)+mup3mup5*p2mup4*p3mup2*
     -   p4mup1*F(1619)+mup1mup3*p2mup4*p3mup2*p4mup5*F(1621)+mup2mup
     -   5*p2mup4*p3mup3*p4mup1*F(1622)+mup1mup2*p2mup4*p3mup3*p4mup5
     -   *F(1624)+mup2mup3*p2mup4*p3mup5*p4mup1*F(1625)+mup1mup3*p2mu
     -   p4*p3mup5*p4mup2*F(1626)+mup1mup2*p2mup4*p3mup5*p4mup3*F(162
     -   7)+mup3mup4*p2mup5*p3mup1*p4mup2*F(1628)+mup2mup4*p2mup5*p3m
     -   up1*p4mup3*F(1629)+mup2mup3*p2mup5*p3mup1*p4mup4*F(1630)+mup
     -   3mup4*p2mup5*p3mup2*p4mup1*F(1631)+mup1mup3*p2mup5*p3mup2*p4
     -   mup4*F(1633)+mup2mup4*p2mup5*p3mup3*p4mup1*F(1634)+mup1mup2*
     -   p2mup5*p3mup3*p4mup4*F(1636)+mup2mup3*p2mup5*p3mup4*p4mup1*F
     -   (1637)+mup1mup3*p2mup5*p3mup4*p4mup2*F(1638)+mup1mup2*p2mup5
     -   *p3mup4*p4mup3*F(1639)+mup2mup5*p3mup3*p3mup4*p4mup1*F(1661)
     -   +mup2mup4*p3mup3*p3mup5*p4mup1*F(1664)+mup1mup2*p2mup3*p4mup
     -   4*p4mup5*F(1685)+mup1mup5*(p2mup3*p3mup2*p4mup4*F(1608)+p2mu
     -   p3*p3mup4*p4mup2*F(1611)+p2mup4*p3mup2*p4mup3*F(1620)+p2mup4
     -   *p3mup3*p4mup2*F(1623)+p2mup2*(p3mup3*p4mup4*F(1596)+p3mup4*
     -   p4mup3*F(1599)+2*p4mup3*p4mup4*F(1681))+2*(p3mup2*p3mup3*p4m
     -   up4*F(1653)+p3mup2*p3mup4*p4mup3*F(1656)+p3mup3*p3mup4*p4mup
     -   2*F(1662)+p2mup3*p4mup2*p4mup4*F(1683)+p2mup4*p4mup2*p4mup3*
     -   F(1686)))+mup1mup4*(p2mup2*p3mup3*p4mup5*F(1597)+p2mup2*p3mu
     -   p5*p4mup3*F(1602)+p2mup3*p3mup2*p4mup5*F(1609)+p2mup3*p3mup5
     -   *p4mup2*F(1614)+p2mup5*p3mup2*p4mup3*F(1632)+p2mup5*p3mup3*p
     -   4mup2*F(1635)+p3mup3*p3mup5*p4mup2*F(1665)-2*(p3mup2*p3mup3*
     -   p4mup5*F(1654)-p3mup2*p3mup5*p4mup3*F(1659)+p2mup2*p4mup3*p4
     -   mup5*F(1682)-p2mup3*p4mup2*p4mup5*F(1684)-p2mup5*p4mup2*p4mu
     -   p3*F(1688)))+mup3mup5*p2mup1*p4mup2*p4mup4*F(1690)+mup3mup4*
     -   p2mup1*p4mup2*p4mup5*F(1691)+mup2mup3*p2mup1*p4mup4*p4mup5*F
     -   (1692)+2*(mup4mup5*p3mup1*p3mup2*p4mup3*F(1640)+mup3mup5*p3m
     -   up1*p3mup2*p4mup4*F(1641)+mup3mup4*p3mup1*p3mup2*p4mup5*F(16
     -   42)+mup4mup5*p3mup1*p3mup3*p4mup2*F(1643)-mup2mup5*p3mup1*p3
     -   mup3*p4mup4*F(1644)-mup2mup4*p3mup1*p3mup3*p4mup5*F(1645)+mu
     -   p3mup5*p3mup1*p3mup4*p4mup2*F(1646)-mup2mup5*p3mup1*p3mup4*p
     -   4mup3*F(1647)+mup2mup3*p3mup1*p3mup4*p4mup5*F(1648)-mup3mup4
     -   *p3mup1*p3mup5*p4mup2*F(1649)+mup2mup4*p3mup1*p3mup5*p4mup3*
     -   F(1650)+mup2mup3*p3mup1*p3mup5*p4mup4*F(1651)+mup4mup5*p3mup
     -   2*p3mup3*p4mup1*F(1652)+mup3mup5*p3mup2*p3mup4*p4mup1*F(1655
     -   )+mup1mup3*p3mup2*p3mup4*p4mup5*F(1657)+mup3mup4*p3mup2*p3mu
     -   p5*p4mup1*F(1658)+mup1mup3*p3mup2*p3mup5*p4mup4*F(1660)+mup1
     -   mup2*p3mup3*p3mup4*p4mup5*F(1663)+mup4mup5*p2mup2*p4mup1*p4m
     -   up3*F(1666)-mup3mup5*p2mup2*p4mup1*p4mup4*F(1667)-mup3mup4*p
     -   2mup2*p4mup1*p4mup5*F(1668)+mup4mup5*p2mup3*p4mup1*p4mup2*F(
     -   1669)-mup2mup5*p2mup3*p4mup1*p4mup4*F(1670)-mup2mup4*p2mup3*
     -   p4mup1*p4mup5*F(1671)-mup3mup5*p2mup4*p4mup1*p4mup2*F(1672)+
     -   mup2mup5*p2mup4*p4mup1*p4mup3*F(1673)+mup2mup3*p2mup4*p4mup1
     -   *p4mup5*F(1674)+mup3mup4*p2mup5*p4mup1*p4mup2*F(1675)-mup2mu
     -   p4*p2mup5*p4mup1*p4mup3*F(1676)+mup2mup3*p2mup5*p4mup1*p4mup
     -   4*F(1677)+mup4mup5*p2mup1*p4mup2*p4mup3*F(1678)+mup2mup5*p2m
     -   up1*p4mup3*p4mup4*F(1679)+mup2mup4*p2mup1*p4mup3*p4mup5*F(16
     -   80)+mup1mup2*p2mup4*p4mup3*p4mup5*F(1687)+mup1mup2*p2mup5*p4
     -   mup3*p4mup4*F(1689)-mup1mup3*p2mup2*p4mup4*p4mup5*F(1693)-mu
     -   p1mup3*p2mup4*p4mup2*p4mup5*F(1694)-mup1mup3*p2mup5*p4mup2*p
     -   4mup4*F(1695))
       Fa5(35) = -(mup2mup3*p2mup5*p3mup4*p4mup1*F(1453))+mup2mup3*p2m
     -   up4*p3mup5*p4mup1*F(1453)+mup3mup5*p2mup2*p3mup4*p4mup1*F(15
     -   81)-mup3mup4*p2mup2*p3mup5*p4mup1*F(1581)+mup3mup5*p2mup4*p3
     -   mup1*p4mup2*F(1582)-mup3mup4*p2mup5*p3mup1*p4mup2*F(1582)+mu
     -   p4mup5*p2mup3*p3mup2*p4mup1*F(1696)+mup3mup5*p2mup4*p3mup2*p
     -   4mup1*F(1697)+mup3mup4*p2mup5*p3mup2*p4mup1*F(1698)+mup4mup5
     -   *p2mup3*p3mup1*p4mup2*F(1699)+mup3mup5*p2mup4*p4mup1*p4mup2*
     -   F(1701)-mup3mup4*p2mup5*p4mup1*p4mup2*F(1701)+mup4mup5*p2mup
     -   2*p4mup1*p4mup3*F(1701)+mup3mup5*p2mup2*p4mup1*p4mup4*F(1701
     -   )-mup3mup4*p2mup2*p4mup1*p4mup5*F(1701)+mup4mup5*p2mup3*p4mu
     -   p1*p4mup2*F(1702)-mup2mup3*p2mup5*p3mup1*p4mup4*F(3792)+mup2
     -   mup3*p2mup4*p3mup1*p4mup5*F(3792)-mup2mup3*p2mup5*p4mup1*p4m
     -   up4*F(3794)+mup2mup3*p2mup4*p4mup1*p4mup5*F(3794)+mup2mup5*(
     -   p2mup4*p3mup3*p4mup1*F(1453)+p2mup3*p3mup4*p4mup1*F(1453)+p2
     -   mup4*p3mup1*p4mup3*F(3792)+p2mup3*p3mup1*p4mup4*F(3792)+p4mu
     -   p1*(p2mup4*p4mup3+p2mup3*p4mup4)*F(3794))-mup2mup4*(p2mup5*p
     -   3mup3*p4mup1*F(1453)+p2mup3*p3mup5*p4mup1*F(1453)+p2mup5*p3m
     -   up1*p4mup3*F(3792)+p2mup3*p3mup1*p4mup5*F(3792)+p4mup1*(p2mu
     -   p5*p4mup3+p2mup3*p4mup5)*F(3794))+2*(mup4mup5*p3mup2*p3mup3*
     -   p4mup1*F(1700)+mup3mup5*p3mup2*p3mup4*p4mup1*F(1700)-mup3mup
     -   4*p3mup2*p3mup5*p4mup1*F(1700)+mup4mup5*p3mup1*p3mup3*p4mup2
     -   *F(1700)+mup3mup5*p3mup1*p3mup4*p4mup2*F(1700)-mup3mup4*p3mu
     -   p1*p3mup5*p4mup2*F(1700)+mup4mup5*p3mup1*p3mup2*p4mup3*F(379
     -   9)+mup3mup5*p3mup1*p3mup2*p4mup4*F(3799)-mup3mup4*p3mup1*p3m
     -   up2*p4mup5*F(3799)+p2mup1*p4mup2*(mup4mup5*p4mup3+mup3mup5*p
     -   4mup4-mup3mup4*p4mup5)*F(3801))
       Fa5(36) = -(mup3mup4*p2mup5*p3mup2*p4mup1*F(1582))-mup3mup4*p2m
     -   up2*p3mup5*p4mup1*F(1582)-mup3mup4*p2mup5*p3mup1*p4mup2*F(15
     -   82)-mup4mup5*p2mup3*(p3mup2*p4mup1+p3mup1*p4mup2)*F(1582)+mu
     -   p3mup5*p2mup4*(p3mup2*p4mup1+p3mup1*p4mup2)*F(1582)-mup4mup5
     -   *p2mup3*p4mup1*p4mup2*F(1701)+mup3mup5*p2mup4*p4mup1*p4mup2*
     -   F(1701)-mup3mup4*p2mup5*p4mup1*p4mup2*F(1701)-mup4mup5*p2mup
     -   2*p4mup1*p4mup3*F(1701)-mup3mup4*p2mup2*p4mup1*p4mup5*F(1701
     -   )+mup3mup5*p2mup2*p4mup1*(p3mup4*F(1582)+p4mup4*F(1701))-mup
     -   2mup5*p2mup4*p3mup3*p4mup1*F(3792)+mup2mup4*p2mup5*p3mup3*p4
     -   mup1*F(3792)+mup2mup5*p2mup3*p3mup4*p4mup1*F(3792)-mup2mup3*
     -   p2mup5*p3mup4*p4mup1*F(3792)-mup2mup4*p2mup3*p3mup5*p4mup1*F
     -   (3792)+mup2mup3*p2mup4*p3mup5*p4mup1*F(3792)-mup2mup5*p2mup4
     -   *p3mup1*p4mup3*F(3792)+mup2mup4*p2mup5*p3mup1*p4mup3*F(3792)
     -   +mup2mup5*p2mup3*p3mup1*p4mup4*F(3792)-mup2mup3*p2mup5*p3mup
     -   1*p4mup4*F(3792)-mup2mup4*p2mup3*p3mup1*p4mup5*F(3792)+mup2m
     -   up3*p2mup4*p3mup1*p4mup5*F(3792)-mup2mup5*p2mup4*p4mup1*p4mu
     -   p3*F(3794)+mup2mup4*p2mup5*p4mup1*p4mup3*F(3794)+mup2mup5*p2
     -   mup3*p4mup1*p4mup4*F(3794)-mup2mup3*p2mup5*p4mup1*p4mup4*F(3
     -   794)-mup2mup4*p2mup3*p4mup1*p4mup5*F(3794)+mup2mup3*p2mup4*p
     -   4mup1*p4mup5*F(3794)-2*(mup3mup4*p3mup2*p3mup5*p4mup1*F(3799
     -   )+mup3mup4*p3mup1*p3mup5*p4mup2*F(3799)+mup4mup5*(p3mup2*p3m
     -   up3*p4mup1+p3mup1*p3mup3*p4mup2+p3mup1*p3mup2*p4mup3)*F(3799
     -   )-mup3mup5*(p3mup2*p3mup4*p4mup1+p3mup1*p3mup4*p4mup2+p3mup1
     -   *p3mup2*p4mup4)*F(3799)+mup3mup4*p3mup1*p3mup2*p4mup5*F(3799
     -   )+p2mup1*p4mup2*(mup4mup5*p4mup3-mup3mup5*p4mup4+mup3mup4*p4
     -   mup5)*F(3801))
       Fa5(37) = -(mup2mup3*p2mup5*p3mup4*p4mup1*F(1453))-mup2mup3*p2m
     -   up4*p3mup5*p4mup1*F(1453)-mup2mup3*p2mup5*p3mup1*p4mup4*F(14
     -   64)-mup2mup3*p2mup4*p3mup1*p4mup5*F(1464)-mup2mup3*p2mup5*p4
     -   mup1*p4mup4*F(1469)-mup2mup3*p2mup4*p4mup1*p4mup5*F(1469)+mu
     -   p2mup4*(p2mup5*p3mup3*p4mup1*F(1453)+p2mup3*p3mup5*p4mup1*F(
     -   1453)+p2mup5*p3mup1*p4mup3*F(1464)+p2mup3*p3mup1*p4mup5*F(14
     -   64)+p4mup1*(p2mup5*p4mup3+p2mup3*p4mup5)*F(1469))+mup2mup5*(
     -   -(p2mup4*(p3mup3*p4mup1*F(1453)+p3mup1*p4mup3*F(1464)+p4mup1
     -   *p4mup3*F(1469)))+p2mup3*(p3mup4*p4mup1*F(1453)+p3mup1*p4mup
     -   4*F(1464)+p4mup1*p4mup4*F(1469)))+mup3mup5*p2mup2*p3mup4*p4m
     -   up1*F(1581)+mup3mup4*p2mup2*p3mup5*p4mup1*F(1581)-mup4mup5*p
     -   2mup3*p3mup2*p4mup1*F(1696)+mup3mup5*p2mup4*p3mup2*p4mup1*F(
     -   1696)+mup3mup4*p2mup5*p3mup2*p4mup1*F(1697)-mup4mup5*p2mup3*
     -   p3mup1*p4mup2*F(1703)+mup3mup5*p2mup4*p3mup1*p4mup2*F(1703)+
     -   mup3mup4*p2mup5*p3mup1*p4mup2*F(1704)-mup4mup5*p2mup2*p4mup1
     -   *p4mup3*F(1705)+mup3mup5*p2mup2*p4mup1*p4mup4*F(1705)+mup3mu
     -   p4*p2mup2*p4mup1*p4mup5*F(1705)-mup4mup5*p2mup3*p4mup1*p4mup
     -   2*F(1706)+mup3mup5*p2mup4*p4mup1*p4mup2*F(1706)+mup3mup4*p2m
     -   up5*p4mup1*p4mup2*F(1707)+2*(-(mup4mup5*p3mup2*p3mup3*p4mup1
     -   *F(1700))+mup3mup5*p3mup2*p3mup4*p4mup1*F(1700)+mup3mup4*p3m
     -   up2*p3mup5*p4mup1*F(1700)-mup4mup5*p3mup1*p3mup3*p4mup2*F(17
     -   00)+mup3mup5*p3mup1*p3mup4*p4mup2*F(1700)+mup3mup4*p3mup1*p3
     -   mup5*p4mup2*F(1700)-mup4mup5*p3mup1*p3mup2*p4mup3*F(1708)+mu
     -   p3mup5*p3mup1*p3mup2*p4mup4*F(1708)+mup3mup4*p3mup1*p3mup2*p
     -   4mup5*F(1708)+p2mup1*p4mup2*(-(mup4mup5*p4mup3)+mup3mup5*p4m
     -   up4+mup3mup4*p4mup5)*F(1709))
       Fa5(38) = mup2mup5*mup3mup4*p1mup1*F(1712)+mup1mup3*mup4mup5*p1
     -   mup2*F(1713)+mup1mup4*mup3mup5*p1mup2*F(1714)+mup1mup5*mup3m
     -   up4*p1mup2*F(1715)+mup1mup2*mup4mup5*p1mup3*F(1716)+mup1mup4
     -   *mup2mup5*p1mup3*F(1717)+mup1mup2*mup3mup5*p1mup4*F(1719)+mu
     -   p1mup3*mup2mup5*p1mup4*F(1720)+mup1mup2*mup3mup4*p1mup5*F(17
     -   22)+mup1mup2*mup4mup5*p2mup3*F(1725)+mup1mup4*mup2mup5*p2mup
     -   3*F(1726)+mup2mup5*mup3mup4*p3mup1*F(1730)+mup1mup3*mup4mup5
     -   *p3mup2*F(1731)+mup1mup4*mup3mup5*p3mup2*F(1732)+mup1mup5*mu
     -   p3mup4*p3mup2*F(1733)+mup1mup2*mup4mup5*p3mup3*F(1734)+mup1m
     -   up4*mup2mup5*p3mup3*F(1735)+mup1mup2*mup3mup5*p3mup4*F(1737)
     -   +mup1mup3*mup2mup5*p3mup4*F(1738)+mup1mup2*mup3mup4*p3mup5*F
     -   (1740)+mup2mup5*mup3mup4*p2mup1*F(1745)+mup1mup3*mup4mup5*p2
     -   mup2*F(1746)+mup1mup4*mup3mup5*p2mup2*F(1747)+mup1mup5*mup3m
     -   up4*p2mup2*F(1748)+mup1mup2*mup3mup5*p2mup4*F(1749)+mup1mup3
     -   *mup2mup5*p2mup4*F(1750)+mup1mup2*mup3mup4*p2mup5*F(1752)+mu
     -   p2mup5*mup3mup4*p4mup1*F(1757)+mup1mup3*mup4mup5*p4mup2*F(17
     -   58)+mup1mup4*mup3mup5*p4mup2*F(1759)+mup1mup5*mup3mup4*p4mup
     -   2*F(1760)+mup1mup2*mup4mup5*p4mup3*F(1761)+mup1mup4*mup2mup5
     -   *p4mup3*F(1762)+mup1mup2*mup3mup5*p4mup4*F(1764)+mup1mup3*mu
     -   p2mup5*p4mup4*F(1765)+mup2mup4*(mup1mup3*(p1mup5*F(1723)+p3m
     -   up5*F(1741)+p2mup5*F(1753))+mup3mup5*(p1mup1*F(1711)+p3mup1*
     -   F(1729)+p2mup1*F(1744)+p4mup1*F(1756))+mup1mup5*(p1mup3*F(17
     -   18)+p2mup3*F(1727)+p3mup3*F(1736)+p4mup3*F(1763))-4*(p3mup3*
     -   p4mup1*p4mup5*F(1773)-p3mup5*p4mup1*p4mup3*F(1778))+2*p3mup1
     -   *p4mup3*p4mup5*F(1784))+mup2mup3*(mup1mup5*p1mup4*F(1721)+mu
     -   p1mup4*p1mup5*F(1724)+mup1mup5*p3mup4*F(1739)+mup1mup4*p3mup
     -   5*F(1742)+mup1mup5*p2mup4*F(1751)+mup1mup4*p2mup5*F(1754)+mu
     -   p4mup5*(p1mup1*F(1710)+p3mup1*F(1728)+p2mup1*F(1743)+p4mup1*
     -   F(1755))+4*(p3mup4*p4mup1*p4mup5*F(1776)-p3mup5*p4mup1*p4mup
     -   4*F(1779))+2*(p3mup4*p3mup5*p4mup1*F(1767)+p3mup1*p4mup4*p4m
     -   up5*F(1785)))+2*mup1mup3*p3mup2*p4mup4*p4mup5*F(1787)+4*(mup
     -   1mup2*p3mup3*p3mup5*p4mup4*F(1766)+mup4mup5*p3mup2*p4mup1*p4
     -   mup3*F(1768)+mup3mup5*p3mup2*p4mup1*p4mup4*F(1769)+mup3mup4*
     -   p3mup2*p4mup1*p4mup5*F(1770)+mup4mup5*p3mup3*p4mup1*p4mup2*F
     -   (1771)-mup2mup5*p3mup3*p4mup1*p4mup4*F(1772)+mup3mup5*p3mup4
     -   *p4mup1*p4mup2*F(1774)-mup2mup5*p3mup4*p4mup1*p4mup3*F(1775)
     -   -mup3mup4*p3mup5*p4mup1*p4mup2*F(1777)+mup1mup2*p3mup4*p3mup
     -   5*p4mup3*F(1780)+mup4mup5*p3mup1*p4mup2*p4mup3*F(1781)+mup3m
     -   up5*p3mup1*p4mup2*p4mup4*F(1782)+mup3mup4*p3mup1*p4mup2*p4mu
     -   p5*F(1783)+mup1mup4*p3mup2*p4mup3*p4mup5*F(1786)+mup1mup5*p3
     -   mup3*p4mup2*p4mup4*F(1788)+mup1mup4*p3mup3*p4mup2*p4mup5*F(1
     -   789)+mup1mup2*p3mup3*p4mup4*p4mup5*F(1790)+mup1mup5*p3mup4*p
     -   4mup2*p4mup3*F(1791)+mup1mup3*p3mup4*p4mup2*p4mup5*F(1792)+m
     -   up1mup2*p3mup4*p4mup3*p4mup5*F(1793)+mup1mup4*p3mup5*p4mup2*
     -   p4mup3*F(1794)-mup1mup3*p3mup5*p4mup2*p4mup4*F(1795)+mup1mup
     -   3*p3mup4*p3mup5*p4mup2*F(1796)+mup2mup5*p3mup1*p4mup3*p4mup4
     -   *F(1797)+mup1mup5*p3mup2*p4mup3*p4mup4*F(1798)+mup1mup2*p3mu
     -   p5*p4mup3*p4mup4*F(1799))+64*p1mup1*p1mup2*p1mup5*p2mup3*p2m
     -   up4*F(1801)+32*(p1mup1*p1mup2*p1mup4*p2mup3*p2mup5*F(1800)+p
     -   1mup1*p1mup2*p1mup3*p2mup4*p2mup5*F(1800)+p1mup1*p1mup4*p1mu
     -   p5*p2mup2*p2mup3*F(1802)+p1mup1*p1mup3*p1mup5*p2mup2*p2mup4*
     -   F(1802)+p1mup2*p1mup4*p1mup5*p2mup1*p2mup3*F(1803)+p1mup2*p1
     -   mup3*p1mup5*p2mup1*p2mup4*F(1803)+p1mup1*p1mup3*p1mup4*p2mup
     -   2*p2mup5*F(1803)+p1mup2*p1mup3*p1mup4*p2mup1*p2mup5*F(1804)+
     -   p1mup3*p1mup4*p1mup5*p2mup1*p2mup2*F(1805))
       Fa5(39) = mup2mup4*mup3mup5*p1mup1*F(1807)+mup1mup3*mup4mup5*p1
     -   mup2*F(1808)+mup1mup4*mup3mup5*p1mup2*F(1809)-mup1mup5*mup3m
     -   up4*p1mup2*F(1810)+mup1mup2*mup4mup5*p1mup3*F(1811)+mup1mup5
     -   *mup2mup4*p1mup3*F(1812)+mup1mup2*mup3mup5*p1mup4*F(1813)-mu
     -   p1mup3*mup2mup4*p1mup5*F(1814)+mup1mup2*mup3mup4*p1mup5*F(18
     -   16)+mup2mup4*mup3mup5*p2mup1*F(1818)+mup1mup3*mup4mup5*p2mup
     -   2*F(1819)+mup1mup4*mup3mup5*p2mup2*F(1820)+mup1mup5*mup3mup4
     -   *p2mup2*F(1821)+mup1mup2*mup4mup5*p2mup3*F(1822)+mup1mup5*mu
     -   p2mup4*p2mup3*F(1824)+mup1mup2*mup3mup5*p2mup4*F(1825)-mup1m
     -   up3*mup2mup4*p2mup5*F(1826)+mup1mup2*mup3mup4*p2mup5*F(1828)
     -   +mup2mup4*mup3mup5*p3mup1*F(1830)+mup1mup3*mup4mup5*p3mup2*F
     -   (1832)+mup1mup4*mup3mup5*p3mup2*F(1833)+mup1mup5*mup3mup4*p3
     -   mup2*F(1834)+mup1mup2*mup4mup5*p3mup3*F(1835)+mup1mup2*mup3m
     -   up5*p3mup4*F(1835)+mup1mup5*mup2mup4*p3mup3*F(1837)+mup1mup3
     -   *mup2mup4*p3mup5*F(1837)+mup1mup2*mup3mup4*p3mup5*F(1838)+mu
     -   p2mup3*(-(mup1mup5*p1mup4*F(1815))+mup1mup4*p1mup5*F(1815)-m
     -   up1mup5*p2mup4*F(1827)+mup1mup4*p2mup5*F(1827)+mup1mup4*p3mu
     -   p5*F(1836)+mup1mup5*p3mup4*F(1837)+mup4mup5*(p1mup1*F(1806)-
     -   p2mup1*F(1817)+p3mup1*F(1829)+p4mup1*F(1839)))+mup2mup4*mup3
     -   mup5*p4mup1*F(1840)+mup1mup2*mup4mup5*p4mup3*F(1842)+mup1mup
     -   2*mup3mup5*p4mup4*F(1842)+mup2mup5*(-(mup1mup4*p1mup3*F(1812
     -   ))+mup1mup3*p1mup4*F(1814)+mup1mup4*p2mup3*F(1823)+mup1mup3*
     -   p2mup4*F(1826)+mup1mup4*p3mup3*F(1836)+mup1mup3*p3mup4*F(183
     -   6)+mup3mup4*(-(p1mup1*F(1807))-p2mup1*F(1818)+p3mup1*F(1831)
     -   +p4mup1*F(1841))+(mup1mup4*p4mup3+mup1mup3*p4mup4)*F(1843))+
     -   mup1mup5*mup2mup4*p4mup3*F(1844)+mup1mup3*mup4mup5*p4mup2*F(
     -   1845)+mup1mup4*mup3mup5*p4mup2*F(1846)+mup1mup5*mup3mup4*p4m
     -   up2*F(1847)+8*(mup4mup5*p3mup3*p4mup1*p4mup2*F(1848)+mup3mup
     -   5*p3mup4*p4mup1*p4mup2*F(1848)-mup3mup4*p3mup5*p4mup1*p4mup2
     -   *F(1848)+(p3mup2*p4mup1+p3mup1*p4mup2)*(mup4mup5*p4mup3+mup3
     -   mup5*p4mup4-mup3mup4*p4mup5)*F(3802))
       Fa5(40) = mup1mup5*mup2mup4*p4mup3*F(1844)+mup1mup3*mup2mup5*p4
     -   mup4*F(1844)+mup2mup5*mup3mup4*p1mup1*F(1849)+mup2mup3*mup4m
     -   up5*p1mup1*F(1849)+mup2mup4*mup3mup5*p1mup1*F(1850)+mup1mup5
     -   *mup3mup4*p1mup2*F(1851)+mup1mup3*mup4mup5*p1mup2*F(1851)+mu
     -   p1mup5*mup2mup3*p1mup4*F(1854)+mup1mup3*mup2mup4*p1mup5*F(18
     -   54)+mup1mup5*mup2mup4*p1mup3*F(1855)+mup1mup3*mup2mup5*p1mup
     -   4*F(1855)+mup2mup5*mup3mup4*p2mup1*F(1857)+mup2mup3*mup4mup5
     -   *p2mup1*F(1857)-mup2mup4*mup3mup5*p2mup1*F(1858)+mup1mup5*mu
     -   p3mup4*p2mup2*F(1859)+mup1mup3*mup4mup5*p2mup2*F(1859)-mup1m
     -   up5*mup2mup4*p2mup3*F(1862)+mup1mup5*mup2mup3*p2mup4*F(1862)
     -   -mup1mup3*mup2mup5*p2mup4*F(1862)+mup1mup3*mup2mup4*p2mup5*F
     -   (1862)+mup2mup5*mup3mup4*p3mup1*F(1864)+mup2mup3*mup4mup5*p3
     -   mup1*F(1864)+mup2mup4*mup3mup5*p3mup1*F(1865)+mup1mup5*mup2m
     -   up3*p3mup4*F(1867)+mup1mup3*mup2mup4*p3mup5*F(1867)+mup1mup5
     -   *mup2mup4*p3mup3*F(1868)+mup1mup3*mup2mup5*p3mup4*F(1868)+mu
     -   p1mup5*mup3mup4*p3mup2*F(1870)+mup1mup3*mup4mup5*p3mup2*F(18
     -   70)+mup2mup5*mup3mup4*p4mup1*F(1872)+mup2mup3*mup4mup5*p4mup
     -   1*F(1872)+mup2mup4*mup3mup5*p4mup1*F(1873)+mup1mup2*(mup4mup
     -   5*p4mup3*F(1842)+mup4mup5*p1mup3*F(1853)+mup3mup4*p1mup5*F(1
     -   853)+mup4mup5*p2mup3*F(1861)+mup3mup4*p2mup5*F(1861)+(mup4mu
     -   p5*p3mup3+mup3mup4*p3mup5)*F(1866)+mup3mup5*(p1mup4*F(1856)+
     -   p2mup4*F(1863)+p3mup4*F(1869)+p4mup4*F(1874)))+mup1mup5*mup3
     -   mup4*p4mup2*F(1875)+mup1mup3*mup4mup5*p4mup2*F(1875)+mup1mup
     -   4*(mup2mup5*(p4mup3*F(1843)+p1mup3*F(1854)+p2mup3*F(1862)+p3
     -   mup3*F(1867))+mup2mup3*(p1mup5*F(1855)-p2mup5*F(1862)+p3mup5
     -   *F(1868))+mup3mup5*(p1mup2*F(1852)+p2mup2*F(1860)+p3mup2*F(1
     -   871)+p4mup2*F(1876)))+8*(mup4mup5*(p3mup3*p4mup1*p4mup2+p3mu
     -   p2*p4mup1*p4mup3+p3mup1*p4mup2*p4mup3)-mup3mup5*(p3mup4*p4mu
     -   p1*p4mup2+p3mup2*p4mup1*p4mup4+p3mup1*p4mup2*p4mup4)+mup3mup
     -   4*(p3mup5*p4mup1*p4mup2+p3mup2*p4mup1*p4mup5+p3mup1*p4mup2*p
     -   4mup5))*F(3802)
       Return
       End
