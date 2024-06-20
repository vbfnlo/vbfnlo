c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenG13(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mu
     -   p1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p
     -   3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup
     -   3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3
     -   ,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa13)
       IMPLICIT NONE
      Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
        common/FPenEvenG/F
      Complex*16 Fa(150),F(4300)
      Real*8 P(2563) 
        COMMON/PPenEvenG/P
       Complex*16 Fa13(121:130)
       COMMON/FaPenEvenG/Fa
       Fa13(121) = mup1mup2*mup3mup4*F(4275)+mup1mup3*mup2mup4*F(4276)
     -   +mup1mup4*mup2mup3*F(4277)+8*(p2mup3*p3mup2*p4mup1*p5mup4*F(
     -   1582)+p2mup3*p3mup1*p4mup2*p5mup4*F(1582)+p2mup3*p4mup1*p4mu
     -   p2*p5mup4*F(1701)-p2mup3*p3mup4*p4mup1*p5mup2*F(4286)-p2mup3
     -   *p3mup1*p4mup4*p5mup2*F(4286)-p2mup3*p4mup1*p4mup4*p5mup2*F(
     -   4288))+2*(mup1mup2*p4mup3*p5mup4*F(1842)+mup1mup4*p4mup3*p5m
     -   up2*F(1843)+mup2mup4*p4mup3*p5mup1*F(1844)+mup3mup4*p1mup1*p
     -   5mup2*F(1849)+mup2mup3*p1mup1*p5mup4*F(1849)+mup2mup4*p1mup1
     -   *p5mup3*F(1850)+mup3mup4*p1mup2*p5mup1*F(1851)+mup1mup3*p1mu
     -   p2*p5mup4*F(1851)+mup1mup4*p1mup2*p5mup3*F(1852)+mup1mup2*p1
     -   mup3*p5mup4*F(1853)+mup2mup3*p1mup4*p5mup1*F(1854)+mup1mup4*
     -   p1mup3*p5mup2*F(1854)+mup2mup4*p1mup3*p5mup1*F(1855)+mup1mup
     -   3*p1mup4*p5mup2*F(1855)+mup1mup2*p1mup4*p5mup3*F(1856)+mup3m
     -   up4*p2mup1*p5mup2*F(1857)+mup2mup3*p2mup1*p5mup4*F(1857)-mup
     -   2mup4*p2mup1*p5mup3*F(1858)+p2mup2*(-4*(p3mup4*p4mup1*p5mup3
     -   *F(1582)-p3mup3*p4mup1*p5mup4*F(1582)+p4mup1*p4mup4*p5mup3*F
     -   (1701)-p4mup1*p4mup3*p5mup4*F(1701))+mup3mup4*p5mup1*F(1859)
     -   +mup1mup3*p5mup4*F(1859)+mup1mup4*p5mup3*F(1860))+mup1mup2*p
     -   2mup3*p5mup4*F(1861)-mup2mup4*p2mup3*p5mup1*F(1862)+mup1mup4
     -   *p2mup3*p5mup2*F(1862)+mup3mup4*p3mup1*p5mup2*F(1864)+mup2mu
     -   p3*p3mup1*p5mup4*F(1864)+mup2mup4*p3mup1*p5mup3*F(1865)+mup1
     -   mup2*p3mup3*p5mup4*F(1866)+mup2mup3*p3mup4*p5mup1*F(1867)+mu
     -   p1mup4*p3mup3*p5mup2*F(1867)+mup2mup4*p3mup3*p5mup1*F(1868)+
     -   mup1mup3*p3mup4*p5mup2*F(1868)+mup1mup2*p3mup4*p5mup3*F(1869
     -   )+mup3mup4*p3mup2*p5mup1*F(1870)+mup1mup3*p3mup2*p5mup4*F(18
     -   70)+mup1mup4*p3mup2*p5mup3*F(1871)+mup3mup4*p4mup1*p5mup2*F(
     -   1872)+mup2mup3*p4mup1*p5mup4*F(1872)+mup2mup4*p4mup1*p5mup3*
     -   F(1873)+mup3mup4*p4mup2*p5mup1*F(1875)+mup1mup3*p4mup2*p5mup
     -   4*F(1875)+mup1mup4*p4mup2*p5mup3*F(1876)+p2mup4*((mup2mup3*p
     -   5mup1-mup1mup3*p5mup2)*F(1862)+mup1mup2*p5mup3*F(1863)-4*(p3
     -   mup2*p4mup1*p5mup3*F(1582)+p3mup1*p4mup2*p5mup3*F(1582)+p4mu
     -   p1*p4mup2*p5mup3*F(1701)-p3mup3*p4mup1*p5mup2*F(4286)-p3mup1
     -   *p4mup3*p5mup2*F(4286)-p4mup1*p4mup3*p5mup2*F(4288))))-16*(p
     -   3mup2*p3mup4*p4mup1*p5mup3*F(4293)+p3mup1*p3mup4*p4mup2*p5mu
     -   p3*F(4293)+p3mup1*p3mup2*p4mup4*p5mup3*F(4293)-p3mup2*p3mup3
     -   *p4mup1*p5mup4*F(4293)-p3mup1*p3mup3*p4mup2*p5mup4*F(4293)-p
     -   3mup1*p3mup2*p4mup3*p5mup4*F(4293)+p2mup1*p4mup2*p4mup4*p5mu
     -   p3*F(4295)-p2mup1*p4mup2*p4mup3*p5mup4*F(4295)-(-((p3mup4*p4
     -   mup1*p4mup2+p3mup2*p4mup1*p4mup4+p3mup1*p4mup2*p4mup4)*p5mup
     -   3)+(p3mup3*p4mup1*p4mup2+p3mup2*p4mup1*p4mup3+p3mup1*p4mup2*
     -   p4mup3)*p5mup4)*F(4296))
       Fa13(122) = 8*(p2mup2*p3mup4*p4mup1*p5mup3*F(1581)-p2mup2*p3mup
     -   3*p4mup1*p5mup4*F(1581)+p2mup2*p4mup1*p4mup4*p5mup3*F(1705)-
     -   p2mup2*p4mup1*p4mup3*p5mup4*F(1705))+2*(mup1mup2*p3mup4*p5mu
     -   p3*F(1835)+mup2mup4*p3mup3*p5mup1*F(1836)+mup1mup3*p3mup4*p5
     -   mup2*F(1836)+mup2mup3*p3mup4*p5mup1*F(1837)+mup1mup4*p3mup3*
     -   p5mup2*F(1837)+mup1mup2*p3mup3*p5mup4*F(1838)+mup2mup3*p1mup
     -   1*p5mup4*F(1877)+mup2mup4*p1mup1*p5mup3*F(1878)+mup3mup4*p1m
     -   up1*p5mup2*F(1879)+mup1mup3*p1mup2*p5mup4*F(1880)+mup1mup4*p
     -   1mup2*p5mup3*F(1881)+mup3mup4*p1mup2*p5mup1*F(1882)+mup1mup2
     -   *p1mup4*p5mup3*F(1883)-mup1mup2*p1mup3*p5mup4*F(1883)+mup1mu
     -   p4*p1mup3*p5mup2*F(1884)-mup1mup3*p1mup4*p5mup2*F(1884)-mup2
     -   mup4*p1mup3*p5mup1*F(1885)+mup2mup3*p1mup4*p5mup1*F(1885)-mu
     -   p2mup4*p2mup1*p5mup3*F(1888)+mup2mup3*p2mup1*p5mup4*F(1888)+
     -   mup3mup4*p2mup1*p5mup2*F(1889)+mup1mup4*p2mup2*p5mup3*F(1890
     -   )-mup1mup3*p2mup2*p5mup4*F(1890)+mup3mup4*p2mup2*p5mup1*F(18
     -   91)+p2mup3*(4*(p3mup4*p4mup1*p5mup2*F(1453)+p3mup1*p4mup4*p5
     -   mup2*F(1464)+p4mup1*p4mup4*p5mup2*F(1469)-p3mup2*p4mup1*p5mu
     -   p4*F(1696)-p3mup1*p4mup2*p5mup4*F(1703)-p4mup1*p4mup2*p5mup4
     -   *F(1706))+mup1mup2*p5mup4*F(1892)+mup1mup4*p5mup2*F(1893)+mu
     -   p2mup4*p5mup1*F(1894))+p2mup4*(-4*(p3mup3*p4mup1*p5mup2*F(14
     -   53)+p3mup1*p4mup3*p5mup2*F(1464)+p4mup1*p4mup3*p5mup2*F(1469
     -   )-p3mup2*p4mup1*p5mup3*F(1696)-p3mup1*p4mup2*p5mup3*F(1703)-
     -   p4mup1*p4mup2*p5mup3*F(1706))+mup1mup2*p5mup3*F(1895)+mup1mu
     -   p3*p5mup2*F(1896)+mup2mup3*p5mup1*F(1897))+mup2mup3*p3mup1*p
     -   5mup4*F(1900)+mup2mup4*p3mup1*p5mup3*F(1901)+mup3mup4*p3mup1
     -   *p5mup2*F(1902)+mup1mup3*p3mup2*p5mup4*F(1903)+mup1mup4*p3mu
     -   p2*p5mup3*F(1904)+mup3mup4*p3mup2*p5mup1*F(1905)+mup2mup3*p4
     -   mup1*p5mup4*F(1906)+mup2mup4*p4mup1*p5mup3*F(1907)+mup3mup4*
     -   p4mup1*p5mup2*F(1908)+mup1mup3*p4mup2*p5mup4*F(1909)+mup1mup
     -   4*p4mup2*p5mup3*F(1910)+mup3mup4*p4mup2*p5mup1*F(1911)+mup1m
     -   up2*p4mup3*p5mup4*F(1912)+mup1mup4*p4mup3*p5mup2*F(1913)+mup
     -   2mup4*p4mup3*p5mup1*F(1914))+16*(p3mup2*p3mup4*p4mup1*p5mup3
     -   *F(1700)+p3mup1*p3mup4*p4mup2*p5mup3*F(1700)-p3mup2*p3mup3*p
     -   4mup1*p5mup4*F(1700)-p3mup1*p3mup3*p4mup2*p5mup4*F(1700)+p3m
     -   up1*p3mup2*p4mup4*p5mup3*F(1708)-p3mup1*p3mup2*p4mup3*p5mup4
     -   *F(1708)+p2mup1*p4mup2*p4mup4*p5mup3*F(1709)-p2mup1*p4mup2*p
     -   4mup3*p5mup4*F(1709)+p3mup4*p4mup1*p4mup2*p5mup3*F(1848)-p3m
     -   up3*p4mup1*p4mup2*p5mup4*F(1848)+p3mup2*p4mup1*p4mup4*p5mup3
     -   *F(1916)+p3mup1*p4mup2*p4mup4*p5mup3*F(1916)-p3mup2*p4mup1*p
     -   4mup3*p5mup4*F(1916)-p3mup1*p4mup2*p4mup3*p5mup4*F(1916))+mu
     -   p1mup2*mup3mup4*F(4278)+mup1mup3*mup2mup4*F(4279)+mup1mup4*m
     -   up2mup3*F(4280)
       Fa13(123) = p4mup4*(mup1mup2*p5mup3*F(1764)+mup1mup3*p5mup2*F(1
     -   765)+mup2mup3*p5mup1*F(1917))
       Fa13(124) = p4mup4*(mup1mup2*p5mup3*F(1842)+mup1mup3*p5mup2*F(1
     -   843)+mup2mup3*p5mup1*F(1844))
       Fa13(125) = p4mup4*(mup2mup3*p5mup1*F(1843)+mup1mup3*p5mup2*F(1
     -   844)+mup1mup2*p5mup3*F(1874))
       Fa13(126) = p4mup4*(mup2mup3*p5mup1*F(1913)+mup1mup3*p5mup2*F(1
     -   914)+mup1mup2*p5mup3*F(1915))
       Fa13(127) = 4*(Fa(2)+Fa(3)+Fa(5)+Fa(6)+Fa(7))+8*(Fa(1)+Fa(4)+Fa
     -   (8)+Fa(9))+Fa(10)+Fa(14)+Fa(18)+Fa(22)+Fa(26)+Fa(30)+Fa(34)+
     -   Fa(38)/2.+Fa(42)
       Fa13(128) = 8*Fa(11)+2*(Fa(19)+Fa(23))+4*(Fa(15)+Fa(27)+Fa(31)+
     -   Fa(35))+Fa(39)+Fa(43)
       Fa13(129) = 8*Fa(12)+2*(Fa(20)+Fa(24))+4*(Fa(16)-Fa(28)+Fa(32)-
     -   Fa(36))+Fa(40)+Fa(44)
       Fa13(130) = 8*Fa(13)+2*(Fa(21)+Fa(25))-4*(Fa(17)-Fa(29)-Fa(33)-
     -   Fa(37))+Fa(41)+Fa(45)
       Return
       End
