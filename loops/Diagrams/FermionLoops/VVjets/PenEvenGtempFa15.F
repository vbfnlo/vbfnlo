c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c***********************************************************************************
c************** Calling the Fa functions*************************
       subroutine FaPenEvenG15(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mu
     -   p1,p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p
     -   3mup5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup
     -   3,p5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3
     -   ,mup2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa15)
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
       Complex*16 Fa15(141:150)
       COMMON/FaPenEvenG/Fa
       Fa15(141) = 4*Fa(71)+2*Fa(75)+Fa(79)-8*Fa(83)
       Fa15(142) = 4*Fa(72)+2*Fa(76)+Fa(80)+8*Fa(84)
       Fa15(143) = -2*Fa(89)+Fa(93)+Fa(97)+4*(Fa(85)+Fa(101))
       Fa15(144) = 4*Fa(90)+2*Fa(94)+Fa(98)+8*(Fa(86)-Fa(102))
       Fa15(145) = 4*Fa(91)+2*Fa(95)+Fa(99)+8*(Fa(87)+Fa(103))
       Fa15(146) = 4*Fa(92)+2*Fa(96)+Fa(100)+8*(Fa(88)+Fa(102))
       Fa15(147) = Fa(104)+Fa(105)+Fa(106)+Fa(107)/2.+Fa(111)/2.+Fa(11
     -   5)+Fa(119)/4.+Fa(123)/2.
       Fa15(148) = 4*Fa(116)+Fa(120)+2*(Fa(108)+Fa(112)+Fa(124))
       Fa15(149) = 4*Fa(117)+Fa(121)+2*(Fa(109)+Fa(113)+Fa(125))
       Fa15(150) = 4*Fa(118)+Fa(122)+2*(Fa(110)+Fa(114)+Fa(126))
       Return
       End
