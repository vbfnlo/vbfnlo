c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenGSim14(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p
     -   2mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup
     -   4,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5
     -   mup3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2m
     -   up3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa14)
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
       Complex*16 Fa14(104:112)
       COMMON/FaPenEvenGSim/Fa
       Fa14(104) = 4*(Fa(2)+Fa(3)+Fa(5)+Fa(6)+Fa(7))+8*(Fa(1)+Fa(4)+Fa
     -   (8)+Fa(9))+Fa(10)+Fa(14)+Fa(18)+Fa(22)+Fa(26)+Fa(30)+Fa(34)+
     -   Fa(38)/2.+Fa(42)
       Fa14(105) = 8*Fa(11)+2*(Fa(19)+Fa(23))+4*(Fa(15)+Fa(27)+Fa(31)+
     -   Fa(35))+Fa(39)+Fa(43)
       Fa14(106) = 8*Fa(12)+2*(Fa(20)+Fa(24))+4*(Fa(16)-Fa(28)+Fa(32)-
     -   Fa(36))+Fa(40)+Fa(44)
       Fa14(107) = 8*Fa(13)+2*(Fa(21)+Fa(25))-4*(Fa(17)-Fa(29)-Fa(33)-
     -   Fa(37))+Fa(41)+Fa(45)
       Fa14(108) = Fa(46)+Fa(50)+Fa(54)+4*Fa(58)
       Fa14(109) = Fa(51)+2*(Fa(47)+Fa(55))
       Fa14(110) = Fa(52)+2*(Fa(48)+Fa(56))
       Fa14(111) = Fa(53)+2*(Fa(49)+Fa(57))
       Fa14(112) = Fa(60)+Fa(64)+4*(Fa(59)+Fa(68)+Fa(69))
       Return
       End
