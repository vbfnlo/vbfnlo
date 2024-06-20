c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenGSim15(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p
     -   2mup1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup
     -   4,p3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5
     -   mup3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2m
     -   up3,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa15)
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
       Complex*16 Fa15(113:123)
       COMMON/FaPenEvenGSim/Fa
       Fa15(113) = 2*Fa(61)+Fa(65)
       Fa15(114) = 2*Fa(62)+Fa(66)
       Fa15(115) = 2*Fa(63)+Fa(67)
       Fa15(116) = Fa(70)+Fa(73)+Fa(77)+Fa(81)
       Fa15(117) = -4*Fa(71)+2*Fa(74)+Fa(78)+8*Fa(82)
       Fa15(118) = 4*Fa(71)+2*Fa(75)+Fa(79)-8*Fa(83)
       Fa15(119) = 4*Fa(72)+2*Fa(76)+Fa(80)+8*Fa(84)
       Fa15(120) = -2*Fa(89)+Fa(93)+Fa(97)+4*(Fa(85)+Fa(101))
       Fa15(121) = 4*Fa(90)+2*Fa(94)+Fa(98)+8*(Fa(86)-Fa(102))
       Fa15(122) = 4*Fa(91)+2*Fa(95)+Fa(99)+8*(Fa(87)+Fa(103))
       Fa15(123) = 4*Fa(92)+2*Fa(96)+Fa(100)+8*(Fa(88)+Fa(102))
       Return
       End
