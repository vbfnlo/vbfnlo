c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenG14(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mu
     -   p1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p
     -   3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup
     -   3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3
     -   ,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa14)
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
       Complex*16 Fa14(131:140)
       COMMON/FaPenEvenG/Fa
       Fa14(131) = Fa(46)+Fa(50)+Fa(54)+4*Fa(58)
       Fa14(132) = Fa(51)+2*(Fa(47)+Fa(55))
       Fa14(133) = Fa(52)+2*(Fa(48)+Fa(56))
       Fa14(134) = Fa(53)+2*(Fa(49)+Fa(57))
       Fa14(135) = Fa(60)+Fa(64)+4*(Fa(59)+Fa(68)+Fa(69))
       Fa14(136) = 2*Fa(61)+Fa(65)
       Fa14(137) = 2*Fa(62)+Fa(66)
       Fa14(138) = 2*Fa(63)+Fa(67)
       Fa14(139) = Fa(70)+Fa(73)+Fa(77)+Fa(81)
       Fa14(140) = -4*Fa(71)+2*Fa(74)+Fa(78)+8*Fa(82)
       Return
       End
