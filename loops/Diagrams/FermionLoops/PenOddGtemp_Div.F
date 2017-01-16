       subroutine PenOddG_Div(cpl,MM1,MM3,MM4,p1,p2,p3,p4,p5, mup1, mup2,m
     -   up3,mup4,mup5,musq,comp,result,Div)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/02/2012
c Modified:27/2/2012
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16   mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
       Real*8 dotrr
       Complex*16 B0finG,C0finG,D0finG,E0finG
       EXTERNAL dotrr,B0finG,C0finG,D0finG,E0finG
        Real*8   p1sq, p1p2, p1p3, p1p4, p1p5 
       Real*8   p2sq, p2p3, p2p4, p2p5 
       Real*8   p3sq, p3p4, p3p5 
       Real*8   p4sq, p4p5 
       Real*8   p5sq 
       Real*8   s12, s13, s14, s15 
       Real*8   s23, s24, s25 
       Real*8   s34, s35 
       Real*8   s45 
       Real*8   s234


       Real*8  B012R,B013R,B014R,B015R 
       Real*8  B023R,B024R,B025R 
       Real*8  B034R,B035R 
       Real*8  B045R 
       Real*8  B012I,B013I,B014I,B015I 
       Real*8  B023I,B024I,B025I 
       Real*8  B034I,B035I 
       Real*8  B045I 
       Real*8  B12R(3,3),B13R(3,3),B14R(3,3),B15R(3,3) 
       Real*8  B23R(3,3),B24R(3,3),B25R(3,3) 
       Real*8  B34R(3,3),B35R(3,3) 
       Real*8  B45R(3,3) 
       Real*8  B12I(3,3),B13I(3,3),B14I(3,3),B15I(3,3) 
       Real*8  B23I(3,3),B24I(3,3),B25I(3,3) 
       Real*8  B34I(3,3),B35I(3,3) 
       Real*8  B45I(3,3)


       Real*8 C0123R,C0124R,C0125R 
       Real*8 C0134R,C0135R 
       Real*8 C0145R 
       Real*8 C0234R,C0235R 
       Real*8 C0245R 
       Real*8 C0345R 
       Real*8 C0123I,C0124I,C0125I 
       Real*8 C0134I,C0135I 
       Real*8 C0145I 
       Real*8 C0234I,C0235I 
       Real*8 C0245I 
       Real*8 C0345I  
       Real*8 C123R(9,4),C124R(9,4),C125R(9,4) 
       Real*8 C134R(9,4),C135R(9,4) 
       Real*8 C145R(9,4) 
       Real*8 C234R(9,4),C235R(9,4) 
       Real*8 C245R(9,4) 
       Real*8 C345R(9,4) 
       Real*8 C123I(9,4),C124I(9,4),C125I(9,4) 
       Real*8 C134I(9,4),C135I(9,4) 
       Real*8 C145I(9,4) 
       Real*8 C234I(9,4),C235I(9,4) 
       Real*8 C245I(9,4) 
       Real*8 C345I(9,4)


       Real*8 D01234R,D01235R 
       Real*8 D01245R 
       Real*8 D01345R 
       Real*8 D02345R  
       Real*8 D01234I,D01235I 
       Real*8 D01245I 
       Real*8 D01345I 
       Real*8 D02345I 
       Real*8 D1234R(34,5),D1235R(34,5) 
       Real*8 D1245R(34,5) 
       Real*8 D1345R(34,5) 
       Real*8 D2345R(34,5) 
       Real*8 D1234I(34,5),D1235I(34,5) 
       Real*8 D1245I(34,5) 
       Real*8 D1345I(34,5) 
       Real*8 D2345I(34,5) 

       Real*8 E012345R 
       Real*8 E012345I 
       Real*8 E12345R(80,5) 
       Real*8 E12345I(80,5) 
       Logical PrintB,PrintC,PrintD,PrintE 
       Integer comp,Div 
       Real*8  MM1,MM3,MM4,M1,M3,M4,m1sq,m3sq,m4sq 
       Real*8 cpl(10) 
      Complex*16 F(7459),Fa(4390)
      Real*8 P(1989)
      Complex*16 dotrc,dotcc,result(1,6)
       Real*8 musq
      EXTERNAL   dotrc,dotcc

      Complex*16 EpsF(126)
       common/EPSFPenOddG/EPSF
      Complex*16 EPSCCCC,EPSCCCR,EPSCCRR,EPSCRRR
       Real*8 EPSRRRR
      EXTERNAL   EPSCCCC,EPSCCCR,EPSCCRR,EPSCRRR,EPSRRRR
       COMMON/FaPenOddG/Fa
       COMMON/FPenOddG/F
       COMMON/PPenOddG/P
      Save/FPenOddG/
      Save/PPenOddG/
       COMMON/InvariantsPenOddG/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s
     -   45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       COMMON/EVALFunctionsPenOddG/ E012345R,E12345R,E012345I,E12345I
       COMMON/DVALFunctionsPenOddG/ D01234R,D1234R,D01234I,D1234I, D01
     -   235R,D1235R,D01235I,D1235I, D01245R,D1245R,D01245I,D1245I, D
     -   01345R,D1345R,D01345I,D1345I, D02345R,D2345R,D02345I,D2345I
       COMMON/CVALFunctionsPenOddG/ C0123R,C123R,C0123I,C123I, C0124R,
     -   C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0134R,C134R,
     -   C0134I,C134I, C0135R,C135R,C0135I,C135I, C0145R,C145R,C0145I
     -   ,C145I, C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I
     -   , C0245R,C245R,C0245I,C245I, C0345R,C345R,C0345I,C345I 
       COMMON/BVALFunctionsPenOddG/ B012R,B012I, B013R,B013I, B014R,B0
     -   14I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,B025I, B03
     -   4R,B034I, B035R,B035I, B045R,B045I, B12R,B12I, B13R,B13I, B1
     -   4R,B14I, B15R,B15I, B23R,B23I, B24R,B24I, B25R,B25I, B34R,B3
     -   4I, B35R,B35I, B45R,B45I 

       parameter(printB=.false., printC=.true.,printD=.false.,printE=.false.)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the scalar products. Not inlcueded the contraction of the
c       moments with the external currents  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p1p5 = dotrr(p1,p5)
       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)
       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)
       p5sq = dotrr(p5,p5)
       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s15 = (p1sq +p5sq+ 2*p1p5) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s25 = (p2sq +p5sq+ 2*p2p5) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s35 = (p3sq +p5sq+ 2*p3p5) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
c       Write(*,'(a5,F20.10)')," p1sq ", p1sq 
c       Write(*,'(a5,F20.10)')," p1p2 ", p1p2
c       Write(*,'(a5,F20.10)')," p1p3 ", p1p3
c       Write(*,'(a5,F20.10)')," p1p4 ", p1p4
c       Write(*,'(a5,F20.10)')," p1p5 ", p1p5
c       Write(*,'(a5,F20.10)')," p2sq ", p2sq 
c       Write(*,'(a5,F20.10)')," p2p3 ", p2p3
c       Write(*,'(a5,F20.10)')," p2p4 ", p2p4
c       Write(*,'(a5,F20.10)')," p2p5 ", p2p5
c       Write(*,'(a5,F20.10)')," p3sq ", p3sq 
c       Write(*,'(a5,F20.10)')," p3p4 ", p3p4
c       Write(*,'(a5,F20.10)')," p3p5 ", p3p5
c       Write(*,'(a5,F20.10)')," p4sq ", p4sq 
c       Write(*,'(a5,F20.10)')," p4p5 ", p4p5
c       Write(*,'(a5,F20.10)')," p5sq ", p5sq 
      s234=s15
       m1=mm1
        m3=mm3
        m4=mm4
       m1sq=m1*m1
      m3sq=m3*m3
      m4sq=m4*m4
c      PrintB=.False. 
c      PrintC=.true. 
c      PrintD=.False. 
c      PrintE=.False.
c************************************************************************************
c************************************************************************************
        if (comp.gt.0) then    
c    Calling C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************
       call  XEgetPenOddG_R5_Div(M1,M1,M1,M3,M4,p1sq,p2sq,p3sq,p4sq,p5sq,p
     -   1p2,p1p3,p1p4,p2p3,p2p4,p3p4,p1p5,p4p5,musq,Div)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call FPenOddG1(F(1))
       call FPenOddG2(F(373))
       call FPenOddG3(F(745))
       call FPenOddG4(F(1117))
       call FPenOddG5(F(1489))
       call FPenOddG6(F(1861))
       call FPenOddG7(F(2233))
       call FPenOddG8(F(2605))
       call FPenOddG9(F(2977))
       call FPenOddG10(F(3349))
       call FPenOddG11(F(3721))
       call FPenOddG12(F(4093))
       call FPenOddG13(F(4465))
       call FPenOddG14(F(4837))
       call FPenOddG15(F(5209))
       call FPenOddG16(F(5581))
       call FPenOddG17(F(5953))
       call FPenOddG18(F(6325))
       call FPenOddG19(F(6697))
       call FPenOddG20(F(7069))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup1 = dotrc(p1,mup1)
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrc(p1,mup3)
       p1mup4 = dotrc(p1,mup4)
       p1mup5 = dotrc(p1,mup5)
       p2mup1 = dotrc(p2,mup1)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup4 = dotrc(p2,mup4)
       p2mup5 = dotrc(p2,mup5)
       p3mup1 = dotrc(p3,mup1)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup4 = dotrc(p3,mup4)
       p3mup5 = dotrc(p3,mup5)
       p4mup1 = dotrc(p4,mup1)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup4 = dotrc(p4,mup4)
       p4mup5 = dotrc(p4,mup5)
       p5mup1 = dotrc(p5,mup1)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrc(p5,mup3)
       p5mup4 = dotrc(p5,mup4)
       p5mup5 = dotrc(p5,mup5)
       mup1mup2 = dotcc(mup1,mup2)
       mup1mup3 = dotcc(mup1,mup3)
       mup1mup4 = dotcc(mup1,mup4)
       mup1mup5 = dotcc(mup1,mup5)
       mup2mup3 = dotcc(mup2,mup3)
       mup2mup4 = dotcc(mup2,mup4)
       mup2mup5 = dotcc(mup2,mup5)
       mup3mup4 = dotcc(mup3,mup4)
       mup3mup5 = dotcc(mup3,mup5)
       mup4mup5 = dotcc(mup4,mup5)
c************** Calling the Fa functions**********************************************************************
       EpsF(1)=EPSCCRR(mup3,mup4,p1,p2)
       EpsF(2)=EPSCCRR(mup3,mup4,p2,p3)
       EpsF(3)=EPSCCRR(mup3,mup4,p2,p4)
       EpsF(4)=EPSCCRR(mup2,mup4,p1,p2)
       EpsF(5)=EPSCCRR(mup2,mup4,p2,p3)
       EpsF(6)=EPSCCRR(mup2,mup4,p2,p4)
       EpsF(7)=EPSCCRR(mup2,mup3,p1,p2)
       EpsF(8)=EPSCCRR(mup2,mup3,p2,p3)
       EpsF(9)=EPSCCRR(mup2,mup3,p2,p4)
       EpsF(10)=EPSCCCR(mup2,mup3,mup4,p1)
       EpsF(11)=EPSCCCR(mup2,mup3,mup4,p2)
       EpsF(12)=EPSCCCR(mup2,mup3,mup4,p3)
       EpsF(13)=EPSCCCR(mup2,mup3,mup4,p4)
       EpsF(14)=EPSCCRR(mup1,mup4,p1,p2)
       EpsF(15)=EPSCCRR(mup1,mup4,p2,p3)
       EpsF(16)=EPSCCRR(mup1,mup4,p2,p4)
       EpsF(17)=EPSCCRR(mup1,mup3,p1,p2)
       EpsF(18)=EPSCCRR(mup1,mup3,p2,p3)
       EpsF(19)=EPSCCRR(mup1,mup3,p2,p4)
       EpsF(20)=EPSCCCC(mup1,mup2,mup3,mup4)
       EpsF(21)=EPSCCCR(mup1,mup3,mup4,p1)
       EpsF(22)=EPSCCCR(mup1,mup3,mup4,p2)
       EpsF(23)=EPSCCCR(mup1,mup3,mup4,p3)
       EpsF(24)=EPSCCCR(mup1,mup3,mup4,p4)
       EpsF(25)=EPSCCRR(mup1,mup2,p1,p2)
       EpsF(26)=EPSCCRR(mup1,mup2,p2,p3)
       EpsF(27)=EPSCCRR(mup1,mup2,p2,p4)
       EpsF(28)=EPSCCCR(mup1,mup2,mup4,p1)
       EpsF(29)=EPSCCCR(mup1,mup2,mup4,p2)
       EpsF(30)=EPSCCCR(mup1,mup2,mup4,p3)
       EpsF(31)=EPSCCCR(mup1,mup2,mup4,p4)
       EpsF(32)=EPSCCCR(mup1,mup2,mup3,p1)
       EpsF(33)=EPSCCCR(mup1,mup2,mup3,p2)
       EpsF(34)=EPSCCCR(mup1,mup2,mup3,p3)
       EpsF(35)=EPSCCCR(mup1,mup2,mup3,p4)
       EpsF(36)=EPSCCCC(mup2,mup3,mup4,mup5)
       EpsF(37)=EPSCCCR(mup2,mup3,mup5,p1)
       EpsF(38)=EPSCCCR(mup2,mup3,mup5,p2)
       EpsF(39)=EPSCCCR(mup2,mup3,mup5,p3)
       EpsF(40)=EPSCCCR(mup2,mup3,mup5,p4)
       EpsF(41)=EPSCCRR(mup2,mup3,p1,p4)
       EpsF(42)=EPSCCRR(mup2,mup3,p3,p4)
       EpsF(43)=EPSCCCR(mup2,mup4,mup5,p1)
       EpsF(44)=EPSCCCR(mup2,mup4,mup5,p2)
       EpsF(45)=EPSCCCR(mup2,mup4,mup5,p3)
       EpsF(46)=EPSCCCR(mup2,mup4,mup5,p4)
       EpsF(47)=EPSCCRR(mup2,mup4,p1,p4)
       EpsF(48)=EPSCCRR(mup2,mup4,p3,p4)
       EpsF(49)=EPSCCRR(mup2,mup5,p1,p4)
       EpsF(50)=EPSCCRR(mup2,mup5,p2,p4)
       EpsF(51)=EPSCCRR(mup2,mup5,p3,p4)
       EpsF(52)=EPSCCCR(mup3,mup4,mup5,p1)
       EpsF(53)=EPSCCCR(mup3,mup4,mup5,p2)
       EpsF(54)=EPSCCCR(mup3,mup4,mup5,p3)
       EpsF(55)=EPSCCCR(mup3,mup4,mup5,p4)
       EpsF(56)=EPSCCRR(mup3,mup4,p1,p4)
       EpsF(57)=EPSCCRR(mup3,mup4,p3,p4)
       EpsF(58)=EPSCCRR(mup3,mup5,p1,p4)
       EpsF(59)=EPSCCRR(mup3,mup5,p2,p4)
       EpsF(60)=EPSCCRR(mup3,mup5,p3,p4)
       EpsF(61)=EPSCCRR(mup4,mup5,p1,p4)
       EpsF(62)=EPSCCRR(mup4,mup5,p2,p4)
       EpsF(63)=EPSCCRR(mup4,mup5,p3,p4)
       EpsF(64)=EPSCCCC(mup1,mup2,mup3,mup5)
       EpsF(65)=EPSCCCR(mup1,mup2,mup5,p2)
       EpsF(66)=EPSCCCR(mup1,mup2,mup5,p3)
       EpsF(67)=EPSCCRR(mup1,mup2,p1,p3)
       EpsF(68)=EPSCCRR(mup1,mup2,p3,p4)
       EpsF(69)=EPSCCCC(mup1,mup3,mup4,mup5)
       EpsF(70)=EPSCCCR(mup1,mup3,mup5,p3)
       EpsF(71)=EPSCCRR(mup1,mup3,p1,p3)
       EpsF(72)=EPSCCRR(mup1,mup3,p3,p4)
       EpsF(73)=EPSCCCR(mup1,mup3,mup5,p2)
       EpsF(74)=EPSCCCR(mup1,mup3,mup5,p4)
       EpsF(75)=EPSCCRR(mup1,mup3,p1,p4)
       EpsF(76)=EPSCCCC(mup1,mup2,mup4,mup5)
       EpsF(77)=EPSCCCR(mup1,mup4,mup5,p1)
       EpsF(78)=EPSCCCR(mup1,mup4,mup5,p2)
       EpsF(79)=EPSCCCR(mup1,mup4,mup5,p3)
       EpsF(80)=EPSCCRR(mup1,mup4,p1,p3)
       EpsF(81)=EPSCCCR(mup1,mup4,mup5,p4)
       EpsF(82)=EPSCCRR(mup1,mup4,p1,p4)
       EpsF(83)=EPSCCRR(mup1,mup4,p3,p4)
       EpsF(84)=EPSCCCR(mup1,mup2,mup5,p1)
       EpsF(85)=EPSCCCR(mup1,mup3,mup5,p1)
       EpsF(86)=EPSCCRR(mup1,mup5,p1,p2)
       EpsF(87)=EPSCCRR(mup1,mup5,p1,p3)
       EpsF(88)=EPSCCRR(mup1,mup5,p2,p3)
       EpsF(89)=EPSCCCR(mup1,mup2,mup5,p4)
       EpsF(90)=EPSCCRR(mup1,mup5,p1,p4)
       EpsF(91)=EPSCCRR(mup1,mup5,p2,p4)
       EpsF(92)=EPSCCRR(mup1,mup5,p3,p4)
       EpsF(93)=EPSCCRR(mup2,mup3,p1,p3)
       EpsF(94)=EPSCCRR(mup2,mup5,p1,p2)
       EpsF(95)=EPSCRRR(mup2,p1,p2,p3)
       EpsF(96)=EPSCRRR(mup2,p1,p2,p4)
       EpsF(97)=EPSCCRR(mup3,mup4,p1,p3)
       EpsF(98)=EPSCCRR(mup3,mup5,p1,p3)
       EpsF(99)=EPSCRRR(mup3,p1,p2,p3)
       EpsF(100)=EPSCRRR(mup3,p1,p3,p4)
       EpsF(101)=EPSCCRR(mup4,mup5,p1,p2)
       EpsF(102)=EPSCCRR(mup4,mup5,p1,p3)
       EpsF(103)=EPSCRRR(mup4,p1,p2,p4)
       EpsF(104)=EPSCRRR(mup4,p1,p3,p4)
       EpsF(105)=EPSCCRR(mup2,mup5,p1,p3)
       EpsF(106)=EPSCCRR(mup3,mup5,p1,p2)
       EpsF(107)=EPSCRRR(mup5,p1,p2,p3)
       EpsF(108)=EPSCRRR(mup5,p1,p2,p4)
       EpsF(109)=EPSCRRR(mup5,p1,p3,p4)
       EpsF(110)=EPSCRRR(mup1,p1,p2,p3)
       EpsF(111)=EPSCRRR(mup1,p2,p3,p4)
       EpsF(112)=EPSCCRR(mup2,mup5,p2,p3)
       EpsF(113)=EPSCRRR(mup2,p2,p3,p4)
       EpsF(114)=EPSCCRR(mup3,mup5,p2,p3)
       EpsF(115)=EPSCRRR(mup3,p2,p3,p4)
       EpsF(116)=EPSCCRR(mup2,mup4,p1,p3)
       EpsF(117)=EPSCRRR(mup2,p1,p3,p4)
       EpsF(118)=EPSCRRR(mup3,p1,p2,p4)
       EpsF(119)=EPSCRRR(mup4,p1,p2,p3)
       EpsF(120)=EPSRRRR(p1,p2,p3,p4)
       EpsF(121)=EPSCCRR(mup1,mup2,p1,p4)
       EpsF(122)=EPSCRRR(mup1,p1,p2,p4)
       EpsF(123)=EPSCRRR(mup4,p2,p3,p4)
       EpsF(124)=EPSCCRR(mup4,mup5,p2,p3)
       EpsF(125)=EPSCRRR(mup1,p1,p3,p4)
       EpsF(126)=EPSCRRR(mup5,p2,p3,p4)
c************************************************************************************
c************************************************************************************
       call FaPenOddG1(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mup
     -   2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p4
     -   mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup4
     -   ,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mup
     -   4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1))
       call FaPenOddG2(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mup
     -   2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p4
     -   mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup4
     -   ,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mup
     -   4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(293))
       call FaPenOddG3(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mup
     -   2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p4
     -   mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup4
     -   ,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mup
     -   4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(585))
       call FaPenOddG4(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mup
     -   2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p4
     -   mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup4
     -   ,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mup
     -   4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(877))
       call FaPenOddG5(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mup
     -   2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p4
     -   mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup4
     -   ,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mup
     -   4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1169))
       call FaPenOddG6(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mup
     -   2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p4
     -   mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup4
     -   ,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mup
     -   4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1461))
       call FaPenOddG7(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mup
     -   2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p4
     -   mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup4
     -   ,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mup
     -   4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1753))
       call FaPenOddG8(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mup
     -   2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p4
     -   mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup4
     -   ,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mup
     -   4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2045))
       call FaPenOddG9(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mup
     -   2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p4
     -   mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup4
     -   ,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mup
     -   4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2337))
       call FaPenOddG10(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2629))
       call FaPenOddG11(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(2921))
       call FaPenOddG12(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(3213))
       call FaPenOddG13(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(3505))
       call FaPenOddG14(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(3797))
       call FaPenOddG15(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(4089))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the Matrix Element  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Amplitude                         
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************


c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       result(1,1) = (0,0.5)*(Cpl(3)*(2*M3*Cpl(9)*(M1*Cpl(10)*(EpsF(20
     -   )*Fa(44)+EpsF(28)*Fa(45)-EpsF(29)*Fa(46)+EpsF(30)*Fa(47)+Eps
     -   F(31)*Fa(48)+EpsF(32)*Fa(49)+EpsF(33)*Fa(50)+EpsF(34)*Fa(51)
     -   +EpsF(35)*Fa(52)+EpsF(36)*Fa(53)+EpsF(64)*Fa(54)+EpsF(69)*Fa
     -   (55)+EpsF(76)*Fa(56)+EpsF(21)*Fa(3134)+EpsF(22)*Fa(3135)+Eps
     -   F(23)*Fa(3136)+EpsF(24)*Fa(3137)+EpsF(65)*Fa(3138)+EpsF(66)*
     -   Fa(3139)+EpsF(84)*Fa(3140)+EpsF(89)*Fa(3141)+EpsF(70)*Fa(314
     -   2)+EpsF(73)*Fa(3143)+EpsF(74)*Fa(3144)+EpsF(77)*Fa(3145)-Eps
     -   F(78)*Fa(3146)+EpsF(79)*Fa(3147)+EpsF(81)*Fa(3148)+EpsF(85)*
     -   Fa(3149)+EpsF(10)*Fa(3158)+EpsF(11)*Fa(3159)+EpsF(12)*Fa(316
     -   0)+EpsF(13)*Fa(3161)+EpsF(37)*Fa(3162)+EpsF(38)*Fa(3163)+Eps
     -   F(39)*Fa(3164)+EpsF(40)*Fa(3165)+EpsF(43)*Fa(3169)+EpsF(44)*
     -   Fa(3170)+EpsF(45)*Fa(3171)+EpsF(46)*Fa(3172)+EpsF(52)*Fa(317
     -   4)+EpsF(53)*Fa(3175)+EpsF(54)*Fa(3176)+EpsF(55)*Fa(3177)-2*(
     -   EpsF(26)*Fa(133)-EpsF(25)*Fa(151)+EpsF(27)*Fa(152)-EpsF(18)*
     -   Fa(209)-EpsF(71)*Fa(212)-EpsF(17)*Fa(230)-EpsF(19)*Fa(231)-E
     -   psF(75)*Fa(234)+EpsF(15)*Fa(277)+EpsF(80)*Fa(282)-EpsF(14)*F
     -   a(291)+EpsF(16)*Fa(292)+EpsF(82)*Fa(296)-EpsF(87)*Fa(366)-Ep
     -   sF(88)*Fa(367)-EpsF(86)*Fa(391)-EpsF(90)*Fa(392)-EpsF(91)*Fa
     -   (393)-EpsF(97)*Fa(3069)+EpsF(98)*Fa(3085)-EpsF(102)*Fa(3088)
     -   -EpsF(1)*Fa(3150)-EpsF(2)*Fa(3151)-EpsF(3)*Fa(3152)-EpsF(4)*
     -   Fa(3153)-EpsF(5)*Fa(3154)-EpsF(6)*Fa(3155)-EpsF(7)*Fa(3156)-
     -   EpsF(9)*Fa(3157)-F(1095)*(EpsF(121)*Fa(137)+EpsF(49)*Fa(541)
     -   +EpsF(108)*Fa(735)-EpsF(41)*Fa(3065)-EpsF(96)*Fa(3068)+EpsF(
     -   118)*Fa(3086)+EpsF(122)*Fa(3090)-EpsF(103)*Fa(3094)+EpsF(47)
     -   *Fa(3166))-F(1051)*(EpsF(67)*Fa(137)+EpsF(105)*Fa(541)+EpsF(
     -   107)*Fa(735)-EpsF(93)*Fa(3065)-EpsF(95)*Fa(3068)+EpsF(99)*Fa
     -   (3086)+EpsF(110)*Fa(3090)-EpsF(119)*Fa(3094)+EpsF(116)*Fa(31
     -   66))-EpsF(56)*Fa(3167)-EpsF(8)*Fa(3168)-EpsF(50)*Fa(3173)+Ep
     -   sF(58)*Fa(3178)-EpsF(59)*Fa(3179)-EpsF(61)*Fa(3180)-EpsF(62)
     -   *Fa(3181)-EpsF(94)*Fa(3182)-EpsF(101)*Fa(3183)-EpsF(106)*Fa(
     -   3184)-EpsF(112)*Fa(3185)-EpsF(114)*Fa(3186)-EpsF(124)*Fa(318
     -   7)))+M4*Cpl(5)*(EpsF(31)*Fa(24)+EpsF(35)*Fa(28)+EpsF(20)*Fa(
     -   33)+EpsF(28)*Fa(34)+EpsF(29)*Fa(35)+EpsF(30)*Fa(36)+EpsF(32)
     -   *Fa(37)+EpsF(33)*Fa(38)+EpsF(34)*Fa(39)+EpsF(36)*Fa(40)+EpsF
     -   (64)*Fa(41)+EpsF(69)*Fa(42)+EpsF(76)*Fa(43)+EpsF(24)*Fa(3036
     -   )+EpsF(89)*Fa(3040)+EpsF(74)*Fa(3043)+EpsF(81)*Fa(3047)+EpsF
     -   (13)*Fa(3060)+EpsF(40)*Fa(3064)+EpsF(46)*Fa(3074)+EpsF(55)*F
     -   a(3079)+EpsF(21)*Fa(3095)+EpsF(22)*Fa(3096)+EpsF(23)*Fa(3097
     -   )+EpsF(65)*Fa(3098)+EpsF(66)*Fa(3099)+EpsF(84)*Fa(3100)+EpsF
     -   (70)*Fa(3101)+EpsF(73)*Fa(3102)+EpsF(77)*Fa(3103)+EpsF(78)*F
     -   a(3104)+EpsF(79)*Fa(3105)+EpsF(85)*Fa(3106)+EpsF(10)*Fa(3112
     -   )+EpsF(11)*Fa(3113)+EpsF(12)*Fa(3114)+EpsF(37)*Fa(3115)+EpsF
     -   (38)*Fa(3116)+EpsF(39)*Fa(3117)+EpsF(43)*Fa(3120)+EpsF(44)*F
     -   a(3121)+EpsF(45)*Fa(3122)+EpsF(52)*Fa(3123)+EpsF(53)*Fa(3124
     -   )+EpsF(54)*Fa(3125)+2*(EpsF(27)*Fa(134)-EpsF(25)*Fa(143)+Eps
     -   F(26)*Fa(144)-EpsF(19)*Fa(210)-EpsF(75)*Fa(215)+EpsF(17)*Fa(
     -   221)-EpsF(18)*Fa(222)-EpsF(71)*Fa(224)+EpsF(16)*Fa(278)+EpsF
     -   (82)*Fa(284)+EpsF(14)*Fa(285)+EpsF(15)*Fa(286)+EpsF(80)*Fa(2
     -   90)-EpsF(90)*Fa(368)-EpsF(91)*Fa(369)+EpsF(86)*Fa(380)-EpsF(
     -   87)*Fa(381)-EpsF(88)*Fa(382)-EpsF(3)*Fa(3051)+EpsF(6)*Fa(305
     -   4)-EpsF(9)*Fa(3056)+EpsF(56)*Fa(3067)-EpsF(50)*Fa(3075)+EpsF
     -   (58)*Fa(3080)+EpsF(59)*Fa(3081)+EpsF(61)*Fa(3082)+EpsF(62)*F
     -   a(3083)-EpsF(1)*Fa(3107)-EpsF(2)*Fa(3108)+EpsF(4)*Fa(3109)+E
     -   psF(5)*Fa(3110)-EpsF(7)*Fa(3111)+EpsF(97)*Fa(3118)-EpsF(8)*F
     -   a(3119)-EpsF(94)*Fa(3126)+EpsF(98)*Fa(3127)+EpsF(101)*Fa(312
     -   8)+EpsF(102)*Fa(3129)+EpsF(106)*Fa(3130)-EpsF(112)*Fa(3131)+
     -   EpsF(114)*Fa(3132)+EpsF(124)*Fa(3133)-(EpsF(67)*Fa(137)-EpsF
     -   (105)*Fa(498)+EpsF(107)*Fa(735)-EpsF(93)*Fa(3065)-EpsF(116)*
     -   Fa(3066)-EpsF(95)*Fa(3068)+EpsF(99)*Fa(3086)+EpsF(110)*Fa(30
     -   90)-EpsF(119)*Fa(3094))*F(7440)-(EpsF(121)*Fa(137)-EpsF(49)*
     -   Fa(498)+EpsF(108)*Fa(735)-EpsF(41)*Fa(3065)-EpsF(47)*Fa(3066
     -   )-EpsF(96)*Fa(3068)+EpsF(118)*Fa(3086)+EpsF(122)*Fa(3090)-Ep
     -   sF(103)*Fa(3094))*F(7441))))+Cpl(4)*(Cpl(5)*(EpsF(20)*Fa(1)+
     -   EpsF(36)*Fa(13)+EpsF(64)*Fa(14)+EpsF(69)*Fa(17)+EpsF(76)*Fa(
     -   18)+2*(EpsF(28)*Fa(5)+EpsF(29)*Fa(6)+EpsF(30)*Fa(7)+EpsF(31)
     -   *Fa(8)+EpsF(32)*Fa(9)+EpsF(33)*Fa(10)+EpsF(34)*Fa(11)+EpsF(3
     -   5)*Fa(12)+EpsF(25)*Fa(2920)+EpsF(26)*Fa(2921)+EpsF(27)*Fa(29
     -   22)+EpsF(67)*Fa(2923)+EpsF(68)*Fa(2924)+EpsF(121)*Fa(2925)+E
     -   psF(17)*Fa(2926)+EpsF(18)*Fa(2927)+EpsF(19)*Fa(2928)+EpsF(21
     -   )*Fa(2929)+EpsF(22)*Fa(2930)+EpsF(23)*Fa(2931)+EpsF(24)*Fa(2
     -   932)+EpsF(65)*Fa(2933)+EpsF(66)*Fa(2934)+EpsF(71)*Fa(2935)+E
     -   psF(72)*Fa(2936)+EpsF(75)*Fa(2937)+EpsF(84)*Fa(2938)+EpsF(89
     -   )*Fa(2939)+EpsF(14)*Fa(2940)+EpsF(15)*Fa(2941)+EpsF(16)*Fa(2
     -   942)+EpsF(70)*Fa(2943)+EpsF(73)*Fa(2944)+EpsF(74)*Fa(2945)+E
     -   psF(77)*Fa(2946)+EpsF(78)*Fa(2947)+EpsF(79)*Fa(2948)+EpsF(80
     -   )*Fa(2949)+EpsF(81)*Fa(2950)+EpsF(82)*Fa(2951)+EpsF(83)*Fa(2
     -   952)+EpsF(85)*Fa(2953)+EpsF(86)*Fa(2954)+EpsF(87)*Fa(2955)+E
     -   psF(88)*Fa(2956)+EpsF(90)*Fa(2957)+EpsF(91)*Fa(2958)+EpsF(92
     -   )*Fa(2959)+EpsF(1)*Fa(2960)+EpsF(2)*Fa(2961)+EpsF(3)*Fa(2962
     -   )+EpsF(4)*Fa(2963)+EpsF(5)*Fa(2964)+EpsF(6)*Fa(2965)+EpsF(7)
     -   *Fa(2966)+EpsF(9)*Fa(2967)+EpsF(10)*Fa(2968)+EpsF(11)*Fa(296
     -   9)+EpsF(12)*Fa(2970)+EpsF(13)*Fa(2971)+EpsF(37)*Fa(2972)+Eps
     -   F(38)*Fa(2973)+EpsF(39)*Fa(2974)+EpsF(40)*Fa(2975)+EpsF(41)*
     -   Fa(2976)+EpsF(42)*Fa(2977)+EpsF(47)*Fa(2978)+EpsF(48)*Fa(297
     -   9)+EpsF(56)*Fa(2980)+EpsF(57)*Fa(2981)+EpsF(93)*Fa(2982)+Eps
     -   F(95)*Fa(2983)+EpsF(96)*Fa(2984)+EpsF(97)*Fa(2985)+EpsF(113)
     -   *Fa(2986)+EpsF(116)*Fa(2987)+EpsF(117)*Fa(2988)+EpsF(8)*Fa(2
     -   989)+EpsF(43)*Fa(2990)+EpsF(44)*Fa(2991)+EpsF(45)*Fa(2992)+E
     -   psF(46)*Fa(2993)+EpsF(49)*Fa(2994)+EpsF(50)*Fa(2995)+EpsF(51
     -   )*Fa(2996)+EpsF(52)*Fa(2997)+EpsF(53)*Fa(2998)+EpsF(54)*Fa(2
     -   999)+EpsF(55)*Fa(3000)+EpsF(58)*Fa(3001)+EpsF(59)*Fa(3002)+E
     -   psF(60)*Fa(3003)+EpsF(61)*Fa(3004)+EpsF(62)*Fa(3005)+EpsF(63
     -   )*Fa(3006)+EpsF(94)*Fa(3007)+EpsF(98)*Fa(3008)+EpsF(99)*Fa(3
     -   009)+EpsF(100)*Fa(3010)+EpsF(101)*Fa(3011)+EpsF(102)*Fa(3012
     -   )+EpsF(105)*Fa(3013)+EpsF(106)*Fa(3014)+EpsF(110)*Fa(3015)+E
     -   psF(111)*Fa(3016)+EpsF(112)*Fa(3017)+EpsF(114)*Fa(3018)+EpsF
     -   (115)*Fa(3019)+EpsF(118)*Fa(3020)+EpsF(122)*Fa(3021)+EpsF(12
     -   4)*Fa(3022)+EpsF(125)*Fa(3023)+EpsF(103)*Fa(3024)+EpsF(104)*
     -   Fa(3025)+EpsF(119)*Fa(3026)+EpsF(123)*Fa(3027)+EpsF(107)*Fa(
     -   3028)+EpsF(108)*Fa(3029)+EpsF(109)*Fa(3030)+EpsF(126)*Fa(303
     -   1)+4*EpsF(120)*Fa(3032)))+2*M1*M4*Cpl(10)*(EpsF(20)*Fa(20)+E
     -   psF(28)*Fa(21)+EpsF(29)*Fa(22)+EpsF(30)*Fa(23)+EpsF(31)*Fa(2
     -   4)+EpsF(32)*Fa(25)+EpsF(33)*Fa(26)+EpsF(34)*Fa(27)+EpsF(35)*
     -   Fa(28)+EpsF(36)*Fa(29)+EpsF(64)*Fa(30)+EpsF(69)*Fa(31)+EpsF(
     -   76)*Fa(32)+EpsF(21)*Fa(3033)+EpsF(22)*Fa(3034)+EpsF(23)*Fa(3
     -   035)+EpsF(24)*Fa(3036)+EpsF(65)*Fa(3037)+EpsF(66)*Fa(3038)+E
     -   psF(84)*Fa(3039)+EpsF(89)*Fa(3040)+EpsF(70)*Fa(3041)+EpsF(73
     -   )*Fa(3042)+EpsF(74)*Fa(3043)+EpsF(77)*Fa(3044)+EpsF(78)*Fa(3
     -   045)+EpsF(79)*Fa(3046)+EpsF(81)*Fa(3047)+EpsF(85)*Fa(3048)+E
     -   psF(10)*Fa(3057)+EpsF(11)*Fa(3058)+EpsF(12)*Fa(3059)+EpsF(13
     -   )*Fa(3060)+EpsF(37)*Fa(3061)+EpsF(38)*Fa(3062)+EpsF(39)*Fa(3
     -   063)+EpsF(40)*Fa(3064)+EpsF(43)*Fa(3071)+EpsF(44)*Fa(3072)+E
     -   psF(45)*Fa(3073)+EpsF(46)*Fa(3074)+EpsF(52)*Fa(3076)+EpsF(53
     -   )*Fa(3077)+EpsF(54)*Fa(3078)+EpsF(55)*Fa(3079)-2*(EpsF(25)*F
     -   a(132)-EpsF(26)*Fa(133)-EpsF(27)*Fa(134)-EpsF(17)*Fa(208)+Ep
     -   sF(18)*Fa(209)+EpsF(19)*Fa(210)+EpsF(71)*Fa(212)+EpsF(75)*Fa
     -   (215)-EpsF(14)*Fa(276)-EpsF(15)*Fa(277)-EpsF(16)*Fa(278)-Eps
     -   F(80)*Fa(282)-EpsF(82)*Fa(284)-EpsF(86)*Fa(365)+EpsF(87)*Fa(
     -   366)+EpsF(88)*Fa(367)+EpsF(90)*Fa(368)+EpsF(91)*Fa(369)+EpsF
     -   (1)*Fa(3049)-EpsF(2)*Fa(3050)+EpsF(3)*Fa(3051)-EpsF(4)*Fa(30
     -   52)-EpsF(5)*Fa(3053)-EpsF(6)*Fa(3054)+EpsF(7)*Fa(3055)+EpsF(
     -   9)*Fa(3056)-EpsF(56)*Fa(3067)+EpsF(97)*Fa(3069)-EpsF(8)*Fa(3
     -   070)+EpsF(50)*Fa(3075)-EpsF(58)*Fa(3080)-EpsF(59)*Fa(3081)-E
     -   psF(61)*Fa(3082)-EpsF(62)*Fa(3083)-EpsF(94)*Fa(3084)-EpsF(98
     -   )*Fa(3085)-EpsF(101)*Fa(3087)+EpsF(102)*Fa(3088)-EpsF(106)*F
     -   a(3089)-EpsF(112)*Fa(3091)-EpsF(114)*Fa(3092)-EpsF(124)*Fa(3
     -   093)-F(1051)*(-(EpsF(67)*Fa(137))+EpsF(105)*Fa(498)-EpsF(107
     -   )*Fa(735)+EpsF(93)*Fa(3065)+EpsF(116)*Fa(3066)+EpsF(95)*Fa(3
     -   068)-EpsF(99)*Fa(3086)-EpsF(110)*Fa(3090)+EpsF(119)*Fa(3094)
     -   )+(EpsF(121)*Fa(137)-EpsF(49)*Fa(498)+EpsF(108)*Fa(735)-EpsF
     -   (41)*Fa(3065)-EpsF(47)*Fa(3066)-EpsF(96)*Fa(3068)+EpsF(118)*
     -   Fa(3086)+EpsF(122)*Fa(3090)-EpsF(103)*Fa(3094))*F(7441)))))+
     -   Cpl(8)*(Cpl(9)*(Cpl(10)*(EpsF(20)*Fa(90)+EpsF(36)*Fa(102)+Ep
     -   sF(64)*Fa(103)+EpsF(69)*Fa(106)+EpsF(76)*Fa(107)+2*(EpsF(28)
     -   *Fa(94)+EpsF(29)*Fa(95)+EpsF(30)*Fa(96)+EpsF(31)*Fa(97)+EpsF
     -   (32)*Fa(98)+EpsF(33)*Fa(99)+EpsF(34)*Fa(100)+EpsF(35)*Fa(101
     -   )-4*EpsF(120)*Fa(3032)+EpsF(25)*Fa(3288)+EpsF(26)*Fa(3289)+E
     -   psF(27)*Fa(3290)+EpsF(67)*Fa(3291)+EpsF(68)*Fa(3292)+EpsF(12
     -   1)*Fa(3293)+EpsF(17)*Fa(3294)+EpsF(18)*Fa(3295)+EpsF(19)*Fa(
     -   3296)+EpsF(21)*Fa(3297)+EpsF(22)*Fa(3298)+EpsF(23)*Fa(3299)+
     -   EpsF(24)*Fa(3300)+EpsF(65)*Fa(3301)+EpsF(66)*Fa(3302)+EpsF(7
     -   1)*Fa(3303)+EpsF(72)*Fa(3304)+EpsF(75)*Fa(3305)+EpsF(84)*Fa(
     -   3306)+EpsF(89)*Fa(3307)+EpsF(14)*Fa(3308)+EpsF(15)*Fa(3309)+
     -   EpsF(16)*Fa(3310)+EpsF(70)*Fa(3311)+EpsF(73)*Fa(3312)+EpsF(7
     -   4)*Fa(3313)+EpsF(77)*Fa(3314)+EpsF(78)*Fa(3315)+EpsF(79)*Fa(
     -   3316)+EpsF(80)*Fa(3317)+EpsF(81)*Fa(3318)+EpsF(82)*Fa(3319)+
     -   EpsF(83)*Fa(3320)+EpsF(85)*Fa(3321)+EpsF(86)*Fa(3322)+EpsF(8
     -   7)*Fa(3323)+EpsF(88)*Fa(3324)+EpsF(90)*Fa(3325)+EpsF(91)*Fa(
     -   3326)+EpsF(92)*Fa(3327)+EpsF(1)*Fa(3328)+EpsF(2)*Fa(3329)+Ep
     -   sF(3)*Fa(3330)+EpsF(4)*Fa(3331)+EpsF(5)*Fa(3332)+EpsF(6)*Fa(
     -   3333)+EpsF(7)*Fa(3334)+EpsF(9)*Fa(3335)+EpsF(10)*Fa(3336)+Ep
     -   sF(11)*Fa(3337)+EpsF(12)*Fa(3338)+EpsF(13)*Fa(3339)+EpsF(37)
     -   *Fa(3340)+EpsF(38)*Fa(3341)+EpsF(39)*Fa(3342)+EpsF(40)*Fa(33
     -   43)+EpsF(41)*Fa(3344)+EpsF(42)*Fa(3345)+EpsF(47)*Fa(3346)+Ep
     -   sF(48)*Fa(3347)+EpsF(56)*Fa(3348)+EpsF(57)*Fa(3349)+EpsF(93)
     -   *Fa(3350)+EpsF(95)*Fa(3351)+EpsF(96)*Fa(3352)+EpsF(97)*Fa(33
     -   53)+EpsF(113)*Fa(3354)+EpsF(116)*Fa(3355)+EpsF(117)*Fa(3356)
     -   +EpsF(8)*Fa(3357)+EpsF(43)*Fa(3358)+EpsF(44)*Fa(3359)+EpsF(4
     -   5)*Fa(3360)+EpsF(46)*Fa(3361)+EpsF(49)*Fa(3362)+EpsF(50)*Fa(
     -   3363)+EpsF(51)*Fa(3364)+EpsF(52)*Fa(3365)+EpsF(53)*Fa(3366)+
     -   EpsF(54)*Fa(3367)+EpsF(55)*Fa(3368)+EpsF(58)*Fa(3369)+EpsF(5
     -   9)*Fa(3370)+EpsF(60)*Fa(3371)+EpsF(61)*Fa(3372)+EpsF(62)*Fa(
     -   3373)+EpsF(63)*Fa(3374)+EpsF(94)*Fa(3375)+EpsF(98)*Fa(3376)+
     -   EpsF(99)*Fa(3377)+EpsF(100)*Fa(3378)+EpsF(101)*Fa(3379)+EpsF
     -   (102)*Fa(3380)+EpsF(105)*Fa(3381)+EpsF(106)*Fa(3382)+EpsF(11
     -   0)*Fa(3383)+EpsF(111)*Fa(3384)+EpsF(112)*Fa(3385)+EpsF(114)*
     -   Fa(3386)+EpsF(115)*Fa(3387)+EpsF(118)*Fa(3388)+EpsF(122)*Fa(
     -   3389)+EpsF(124)*Fa(3390)+EpsF(125)*Fa(3391)+EpsF(103)*Fa(339
     -   2)+EpsF(104)*Fa(3393)+EpsF(119)*Fa(3394)+EpsF(123)*Fa(3395)+
     -   EpsF(107)*Fa(3396)+EpsF(108)*Fa(3397)+EpsF(109)*Fa(3398)+Eps
     -   F(126)*Fa(3399)))+2*M1*M4*Cpl(5)*(EpsF(30)*Fa(47)+EpsF(34)*F
     -   a(51)+EpsF(31)*Fa(72)+EpsF(35)*Fa(76)+EpsF(20)*Fa(81)+EpsF(2
     -   8)*Fa(82)-EpsF(29)*Fa(83)+EpsF(32)*Fa(84)+EpsF(33)*Fa(85)+Ep
     -   sF(36)*Fa(86)+EpsF(64)*Fa(87)+EpsF(69)*Fa(88)+EpsF(76)*Fa(89
     -   )+EpsF(23)*Fa(3136)+EpsF(66)*Fa(3139)+EpsF(70)*Fa(3142)+EpsF
     -   (79)*Fa(3147)+EpsF(12)*Fa(3160)+EpsF(39)*Fa(3164)+EpsF(45)*F
     -   a(3171)+EpsF(54)*Fa(3176)+EpsF(24)*Fa(3226)+EpsF(89)*Fa(3230
     -   )+EpsF(74)*Fa(3233)+EpsF(81)*Fa(3237)+EpsF(13)*Fa(3245)+EpsF
     -   (40)*Fa(3249)+EpsF(46)*Fa(3255)+EpsF(55)*Fa(3259)+EpsF(21)*F
     -   a(3266)+EpsF(22)*Fa(3267)+EpsF(65)*Fa(3268)+EpsF(84)*Fa(3269
     -   )+EpsF(73)*Fa(3270)+EpsF(77)*Fa(3271)-EpsF(78)*Fa(3272)+EpsF
     -   (85)*Fa(3273)+EpsF(10)*Fa(3277)-EpsF(11)*Fa(3278)+EpsF(37)*F
     -   a(3279)+EpsF(38)*Fa(3280)+EpsF(43)*Fa(3281)+EpsF(44)*Fa(3282
     -   )+EpsF(52)*Fa(3283)+EpsF(53)*Fa(3284)-2*(EpsF(26)*Fa(133)+Ep
     -   sF(27)*Fa(134)-EpsF(25)*Fa(169)-EpsF(18)*Fa(209)-EpsF(19)*Fa
     -   (210)-EpsF(71)*Fa(212)-EpsF(75)*Fa(215)-EpsF(17)*Fa(241)+Eps
     -   F(15)*Fa(277)+EpsF(16)*Fa(278)+EpsF(80)*Fa(282)+EpsF(82)*Fa(
     -   284)-EpsF(14)*Fa(300)-EpsF(87)*Fa(366)-EpsF(88)*Fa(367)-EpsF
     -   (90)*Fa(368)-EpsF(91)*Fa(369)-EpsF(86)*Fa(403)-EpsF(3)*Fa(30
     -   51)-EpsF(9)*Fa(3056)-EpsF(97)*Fa(3069)-EpsF(50)*Fa(3075)+Eps
     -   F(58)*Fa(3080)+EpsF(59)*Fa(3081)+EpsF(98)*Fa(3085)-EpsF(102)
     -   *Fa(3088)-EpsF(2)*Fa(3151)-EpsF(5)*Fa(3154)-F(1051)*(EpsF(67
     -   )*Fa(137)+EpsF(105)*Fa(541)+EpsF(107)*Fa(735)-EpsF(93)*Fa(30
     -   65)-EpsF(95)*Fa(3068)+EpsF(99)*Fa(3086)+EpsF(110)*Fa(3090)-E
     -   psF(119)*Fa(3094)+EpsF(116)*Fa(3166))-EpsF(8)*Fa(3168)-EpsF(
     -   112)*Fa(3185)-EpsF(114)*Fa(3186)-EpsF(124)*Fa(3187)-EpsF(6)*
     -   Fa(3241)-EpsF(56)*Fa(3250)-EpsF(61)*Fa(3260)-EpsF(62)*Fa(326
     -   1)-EpsF(1)*Fa(3274)+EpsF(4)*Fa(3275)-EpsF(7)*Fa(3276)-EpsF(9
     -   4)*Fa(3285)-EpsF(101)*Fa(3286)-EpsF(106)*Fa(3287)-(EpsF(121)
     -   *Fa(137)+EpsF(49)*Fa(541)+EpsF(108)*Fa(735)-EpsF(41)*Fa(3065
     -   )-EpsF(96)*Fa(3068)+EpsF(118)*Fa(3086)+EpsF(122)*Fa(3090)-Ep
     -   sF(103)*Fa(3094)+EpsF(47)*Fa(3166))*F(7441))))+2*M3*Cpl(4)*(
     -   M1*Cpl(5)*(EpsF(30)*Fa(23)+EpsF(34)*Fa(27)+EpsF(20)*Fa(57)+E
     -   psF(28)*Fa(58)+EpsF(29)*Fa(59)+EpsF(31)*Fa(60)+EpsF(32)*Fa(6
     -   1)-EpsF(33)*Fa(62)+EpsF(35)*Fa(63)+EpsF(36)*Fa(64)+EpsF(64)*
     -   Fa(65)+EpsF(69)*Fa(66)+EpsF(76)*Fa(67)+EpsF(23)*Fa(3035)+Eps
     -   F(66)*Fa(3038)+EpsF(70)*Fa(3041)+EpsF(79)*Fa(3046)+EpsF(12)*
     -   Fa(3059)+EpsF(39)*Fa(3063)+EpsF(45)*Fa(3073)+EpsF(54)*Fa(307
     -   8)+EpsF(21)*Fa(3188)+EpsF(22)*Fa(3189)+EpsF(24)*Fa(3190)+Eps
     -   F(65)*Fa(3191)+EpsF(84)*Fa(3192)+EpsF(89)*Fa(3193)-EpsF(73)*
     -   Fa(3194)+EpsF(74)*Fa(3195)+EpsF(77)*Fa(3196)+EpsF(78)*Fa(319
     -   7)+EpsF(81)*Fa(3198)+EpsF(85)*Fa(3199)+EpsF(10)*Fa(3206)+Eps
     -   F(11)*Fa(3207)+EpsF(13)*Fa(3208)+EpsF(37)*Fa(3209)+EpsF(38)*
     -   Fa(3210)+EpsF(40)*Fa(3211)+EpsF(43)*Fa(3212)+EpsF(44)*Fa(321
     -   3)+EpsF(46)*Fa(3214)+EpsF(52)*Fa(3216)+EpsF(53)*Fa(3217)+Eps
     -   F(55)*Fa(3218)+2*(EpsF(26)*Fa(133)+EpsF(27)*Fa(152)+EpsF(25)
     -   *Fa(159)-EpsF(18)*Fa(209)-EpsF(71)*Fa(212)-EpsF(19)*Fa(231)-
     -   EpsF(75)*Fa(234)+EpsF(17)*Fa(238)+EpsF(15)*Fa(277)+EpsF(80)*
     -   Fa(282)+EpsF(16)*Fa(292)+EpsF(82)*Fa(296)+EpsF(14)*Fa(297)-E
     -   psF(87)*Fa(366)-EpsF(88)*Fa(367)-EpsF(90)*Fa(392)-EpsF(91)*F
     -   a(393)+EpsF(86)*Fa(397)+EpsF(2)*Fa(3050)+EpsF(5)*Fa(3053)-Ep
     -   sF(97)*Fa(3069)+EpsF(8)*Fa(3070)+EpsF(98)*Fa(3085)-EpsF(102)
     -   *Fa(3088)+EpsF(112)*Fa(3091)+EpsF(114)*Fa(3092)+EpsF(124)*Fa
     -   (3093)-F(1095)*(EpsF(121)*Fa(137)-EpsF(49)*Fa(498)+EpsF(108)
     -   *Fa(735)-EpsF(41)*Fa(3065)-EpsF(47)*Fa(3066)-EpsF(96)*Fa(306
     -   8)+EpsF(118)*Fa(3086)+EpsF(122)*Fa(3090)-EpsF(103)*Fa(3094))
     -   -F(1051)*(EpsF(67)*Fa(137)-EpsF(105)*Fa(498)+EpsF(107)*Fa(73
     -   5)-EpsF(93)*Fa(3065)-EpsF(116)*Fa(3066)-EpsF(95)*Fa(3068)+Ep
     -   sF(99)*Fa(3086)+EpsF(110)*Fa(3090)-EpsF(119)*Fa(3094))-EpsF(
     -   56)*Fa(3167)+EpsF(58)*Fa(3178)-EpsF(61)*Fa(3180)-EpsF(94)*Fa
     -   (3182)+EpsF(1)*Fa(3200)+EpsF(3)*Fa(3201)+EpsF(4)*Fa(3202)+Ep
     -   sF(6)*Fa(3203)+EpsF(7)*Fa(3204)+EpsF(9)*Fa(3205)+EpsF(50)*Fa
     -   (3215)+EpsF(59)*Fa(3219)+EpsF(62)*Fa(3220)+EpsF(101)*Fa(3221
     -   )+EpsF(106)*Fa(3222)))+M4*Cpl(10)*(EpsF(20)*Fa(68)+EpsF(28)*
     -   Fa(69)+EpsF(29)*Fa(70)+EpsF(30)*Fa(71)+EpsF(31)*Fa(72)+EpsF(
     -   32)*Fa(73)+EpsF(33)*Fa(74)+EpsF(34)*Fa(75)+EpsF(35)*Fa(76)+E
     -   psF(36)*Fa(77)+EpsF(64)*Fa(78)+EpsF(69)*Fa(79)+EpsF(76)*Fa(8
     -   0)+EpsF(21)*Fa(3223)+EpsF(22)*Fa(3224)+EpsF(23)*Fa(3225)+Eps
     -   F(24)*Fa(3226)+EpsF(65)*Fa(3227)+EpsF(66)*Fa(3228)+EpsF(84)*
     -   Fa(3229)+EpsF(89)*Fa(3230)+EpsF(70)*Fa(3231)+EpsF(73)*Fa(323
     -   2)+EpsF(74)*Fa(3233)+EpsF(77)*Fa(3234)+EpsF(78)*Fa(3235)+Eps
     -   F(79)*Fa(3236)+EpsF(81)*Fa(3237)+EpsF(85)*Fa(3238)+EpsF(10)*
     -   Fa(3242)-EpsF(11)*Fa(3243)+EpsF(12)*Fa(3244)+EpsF(13)*Fa(324
     -   5)+EpsF(37)*Fa(3246)+EpsF(38)*Fa(3247)+EpsF(39)*Fa(3248)+Eps
     -   F(40)*Fa(3249)+EpsF(43)*Fa(3252)+EpsF(44)*Fa(3253)+EpsF(45)*
     -   Fa(3254)+EpsF(46)*Fa(3255)+EpsF(52)*Fa(3256)+EpsF(53)*Fa(325
     -   7)+EpsF(54)*Fa(3258)+EpsF(55)*Fa(3259)-2*(EpsF(27)*Fa(134)-E
     -   psF(25)*Fa(143)+EpsF(26)*Fa(144)-EpsF(19)*Fa(210)-EpsF(75)*F
     -   a(215)-EpsF(18)*Fa(222)-EpsF(71)*Fa(224)-EpsF(17)*Fa(240)+Ep
     -   sF(16)*Fa(278)+EpsF(82)*Fa(284)+EpsF(15)*Fa(286)+EpsF(80)*Fa
     -   (290)-EpsF(14)*Fa(299)-EpsF(90)*Fa(368)-EpsF(91)*Fa(369)-Eps
     -   F(87)*Fa(381)-EpsF(88)*Fa(382)-EpsF(86)*Fa(398)-EpsF(3)*Fa(3
     -   051)-EpsF(9)*Fa(3056)-EpsF(50)*Fa(3075)+EpsF(58)*Fa(3080)+Ep
     -   sF(59)*Fa(3081)-EpsF(2)*Fa(3108)+EpsF(4)*Fa(3109)-EpsF(7)*Fa
     -   (3111)-EpsF(8)*Fa(3119)-EpsF(94)*Fa(3126)+EpsF(98)*Fa(3127)-
     -   EpsF(112)*Fa(3131)+EpsF(114)*Fa(3132)-EpsF(1)*Fa(3239)-EpsF(
     -   5)*Fa(3240)-EpsF(6)*Fa(3241)-EpsF(56)*Fa(3250)-EpsF(97)*Fa(3
     -   251)-EpsF(61)*Fa(3260)-EpsF(62)*Fa(3261)-EpsF(101)*Fa(3262)-
     -   EpsF(102)*Fa(3263)-EpsF(106)*Fa(3264)-EpsF(124)*Fa(3265)-(Ep
     -   sF(67)*Fa(137)+EpsF(105)*Fa(541)+EpsF(107)*Fa(735)-EpsF(93)*
     -   Fa(3065)-EpsF(95)*Fa(3068)+EpsF(99)*Fa(3086)+EpsF(110)*Fa(30
     -   90)-EpsF(119)*Fa(3094)+EpsF(116)*Fa(3166))*F(7440)-(EpsF(121
     -   )*Fa(137)+EpsF(49)*Fa(541)+EpsF(108)*Fa(735)-EpsF(41)*Fa(306
     -   5)-EpsF(96)*Fa(3068)+EpsF(118)*Fa(3086)+EpsF(122)*Fa(3090)-E
     -   psF(103)*Fa(3094)+EpsF(47)*Fa(3166))*F(7441))))))
       result(1,2) = (0,0.5)*(Cpl(3)*(Cpl(4)*(Cpl(5)*(EpsF(36)*F(2892)
     -   +2*(EpsF(4)*Fa(1064)+EpsF(5)*Fa(1065)+EpsF(6)*Fa(1066)+EpsF(
     -   7)*Fa(1067)+EpsF(8)*Fa(1068)+EpsF(9)*Fa(1069)+EpsF(10)*Fa(10
     -   70)+EpsF(11)*Fa(1071)+EpsF(12)*Fa(1072)+EpsF(13)*Fa(1073)+Ep
     -   sF(37)*Fa(1074)+EpsF(38)*Fa(1075)+EpsF(39)*Fa(1076)+EpsF(40)
     -   *Fa(1077)+EpsF(41)*Fa(1078)+EpsF(42)*Fa(1079)+EpsF(43)*Fa(10
     -   80)+EpsF(44)*Fa(1081)+EpsF(45)*Fa(1082)+EpsF(46)*Fa(1083)+Ep
     -   sF(47)*Fa(1084)+EpsF(48)*Fa(1085)+EpsF(52)*Fa(1086)+EpsF(53)
     -   *Fa(1087)+EpsF(54)*Fa(1088)+EpsF(55)*Fa(1089)+EpsF(93)*Fa(10
     -   90)+EpsF(116)*Fa(1094)+EpsF(1)*Fa(3400)+EpsF(2)*Fa(3401)+Eps
     -   F(3)*Fa(3402)+EpsF(56)*Fa(3403)+EpsF(57)*Fa(3404)+EpsF(97)*F
     -   a(3407)+EpsF(117)*Fa(3409)+EpsF(49)*Fa(3410)+EpsF(50)*Fa(341
     -   1)+EpsF(51)*Fa(3412)+EpsF(58)*Fa(3413)+EpsF(59)*Fa(3414)+Eps
     -   F(60)*Fa(3415)+EpsF(61)*Fa(3416)+EpsF(62)*Fa(3417)+EpsF(63)*
     -   Fa(3418)+EpsF(94)*Fa(3419)+EpsF(98)*Fa(3420)+EpsF(101)*Fa(34
     -   23)+EpsF(102)*Fa(3424)+EpsF(105)*Fa(3425)+EpsF(106)*Fa(3426)
     -   +EpsF(112)*Fa(3427)+EpsF(114)*Fa(3428)+EpsF(115)*Fa(3429)+Ep
     -   sF(124)*Fa(3431)+EpsF(119)*Fa(3434)+EpsF(107)*Fa(3436)+EpsF(
     -   109)*Fa(3438)-2*(EpsF(95)*Fa(3405)+EpsF(96)*Fa(3406)-EpsF(11
     -   3)*Fa(3408)-EpsF(99)*Fa(3421)+EpsF(100)*Fa(3422)-EpsF(118)*F
     -   a(3430)-EpsF(103)*Fa(3432)-EpsF(104)*Fa(3433)-EpsF(123)*Fa(3
     -   435)+EpsF(108)*Fa(3437)-EpsF(126)*Fa(3439))+EpsF(120)*Fa(344
     -   0)))+2*M1*M4*Cpl(10)*(-(EpsF(36)*F(3113))+2*(EpsF(4)*Fa(1096
     -   )+EpsF(5)*Fa(1097)+EpsF(6)*Fa(1098)+EpsF(7)*Fa(1099)+EpsF(8)
     -   *Fa(1100)-EpsF(9)*Fa(1101)+EpsF(10)*Fa(1102)+EpsF(11)*Fa(110
     -   3)+EpsF(12)*Fa(1104)+EpsF(13)*Fa(1105)+EpsF(37)*Fa(1106)+Eps
     -   F(38)*Fa(1107)-EpsF(39)*Fa(1108)-EpsF(40)*Fa(1109)-EpsF(41)*
     -   Fa(1110)+EpsF(43)*Fa(1111)+EpsF(44)*Fa(1112)+EpsF(45)*Fa(111
     -   3)+EpsF(46)*Fa(1114)+EpsF(47)*Fa(1115)+EpsF(52)*Fa(1116)+Eps
     -   F(53)*Fa(1117)+EpsF(54)*Fa(1118)+EpsF(55)*Fa(1119)-EpsF(93)*
     -   Fa(1120)+EpsF(116)*Fa(1121)-2*(EpsF(95)*Fa(1271)+EpsF(96)*Fa
     -   (1272)-EpsF(99)*Fa(1345)-EpsF(118)*Fa(1350)+EpsF(103)*Fa(141
     -   1)+EpsF(119)*Fa(1412)-EpsF(107)*Fa(1438)-EpsF(108)*Fa(1439))
     -   +EpsF(1)*Fa(3441)+EpsF(2)*Fa(3442)-EpsF(3)*Fa(3443)+EpsF(56)
     -   *Fa(3444)+EpsF(97)*Fa(3445)+EpsF(49)*Fa(3446)-EpsF(50)*Fa(34
     -   47)+EpsF(58)*Fa(3448)+EpsF(59)*Fa(3449)+EpsF(61)*Fa(3450)+Ep
     -   sF(62)*Fa(3451)+EpsF(94)*Fa(3452)+EpsF(98)*Fa(3453)+EpsF(101
     -   )*Fa(3454)+EpsF(102)*Fa(3455)+EpsF(105)*Fa(3456)+EpsF(106)*F
     -   a(3457)+EpsF(112)*Fa(3458)+EpsF(114)*Fa(3459)+EpsF(124)*Fa(3
     -   460))))+2*M3*Cpl(9)*(M4*Cpl(5)*(EpsF(36)*F(3166)+2*(EpsF(6)*
     -   Fa(1098)-EpsF(9)*Fa(1101)+EpsF(13)*Fa(1105)-EpsF(40)*Fa(1109
     -   )-EpsF(41)*Fa(1110)+EpsF(46)*Fa(1114)+EpsF(47)*Fa(1115)+EpsF
     -   (55)*Fa(1119)+EpsF(4)*Fa(1122)+EpsF(5)*Fa(1123)+EpsF(7)*Fa(1
     -   124)-EpsF(8)*Fa(1125)+EpsF(10)*Fa(1126)+EpsF(11)*Fa(1127)+Ep
     -   sF(12)*Fa(1128)-EpsF(37)*Fa(1129)-EpsF(38)*Fa(1130)-EpsF(39)
     -   *Fa(1131)+EpsF(43)*Fa(1132)+EpsF(44)*Fa(1133)+EpsF(45)*Fa(11
     -   34)+EpsF(52)*Fa(1135)+EpsF(53)*Fa(1136)+EpsF(54)*Fa(1137)-Ep
     -   sF(93)*Fa(1138)+EpsF(116)*Fa(1139)-2*(EpsF(96)*Fa(1272)+EpsF
     -   (95)*Fa(1279)-EpsF(118)*Fa(1350)-EpsF(99)*Fa(1353)+EpsF(103)
     -   *Fa(1411)+EpsF(119)*Fa(1415)-EpsF(108)*Fa(1439)-EpsF(107)*Fa
     -   (1440))-EpsF(3)*Fa(3443)+EpsF(56)*Fa(3444)+EpsF(49)*Fa(3446)
     -   -EpsF(50)*Fa(3447)+EpsF(58)*Fa(3448)+EpsF(59)*Fa(3449)+EpsF(
     -   61)*Fa(3450)+EpsF(62)*Fa(3451)+EpsF(1)*Fa(3461)-EpsF(2)*Fa(3
     -   462)+EpsF(97)*Fa(3463)+EpsF(94)*Fa(3464)+EpsF(98)*Fa(3465)+E
     -   psF(101)*Fa(3466)+EpsF(102)*Fa(3467)+EpsF(105)*Fa(3468)-EpsF
     -   (106)*Fa(3469)-EpsF(112)*Fa(3470)+EpsF(114)*Fa(3471)+EpsF(12
     -   4)*Fa(3472)))+M1*Cpl(10)*(EpsF(36)*F(3200)-2*(EpsF(4)*Fa(114
     -   0)-EpsF(5)*Fa(1141)-EpsF(6)*Fa(1142)-EpsF(7)*Fa(1143)-EpsF(8
     -   )*Fa(1144)-EpsF(9)*Fa(1145)-EpsF(10)*Fa(1146)-EpsF(11)*Fa(11
     -   47)+EpsF(12)*Fa(1148)+EpsF(13)*Fa(1149)-EpsF(37)*Fa(1150)-Ep
     -   sF(38)*Fa(1151)-EpsF(39)*Fa(1152)-EpsF(40)*Fa(1153)-EpsF(41)
     -   *Fa(1154)+EpsF(43)*Fa(1155)+EpsF(44)*Fa(1156)+EpsF(45)*Fa(11
     -   57)+EpsF(46)*Fa(1158)+EpsF(47)*Fa(1159)-EpsF(52)*Fa(1160)-Ep
     -   sF(53)*Fa(1161)+EpsF(54)*Fa(1162)-EpsF(55)*Fa(1163)-EpsF(93)
     -   *Fa(1164)+EpsF(116)*Fa(1165)-2*(EpsF(95)*Fa(1271)+EpsF(96)*F
     -   a(1290)-EpsF(99)*Fa(1345)-EpsF(118)*Fa(1366)+EpsF(119)*Fa(14
     -   12)+EpsF(103)*Fa(1419)-EpsF(107)*Fa(1438)-EpsF(108)*Fa(1441)
     -   )-EpsF(1)*Fa(3473)-EpsF(2)*Fa(3474)-EpsF(3)*Fa(3475)-EpsF(56
     -   )*Fa(3476)-EpsF(97)*Fa(3477)-EpsF(49)*Fa(3478)-EpsF(50)*Fa(3
     -   479)-EpsF(58)*Fa(3480)-EpsF(59)*Fa(3481)-EpsF(61)*Fa(3482)-E
     -   psF(62)*Fa(3483)-EpsF(94)*Fa(3484)-EpsF(98)*Fa(3485)+EpsF(10
     -   1)*Fa(3486)-EpsF(102)*Fa(3487)-EpsF(105)*Fa(3488)-EpsF(106)*
     -   Fa(3489)-EpsF(112)*Fa(3490)-EpsF(114)*Fa(3491)-EpsF(124)*Fa(
     -   3492)))))+Cpl(8)*(2*M3*Cpl(4)*(M1*Cpl(5)*(-(EpsF(36)*F(3200)
     -   )+2*(EpsF(5)*Fa(1097)+EpsF(8)*Fa(1100)+EpsF(12)*Fa(1104)-Eps
     -   F(39)*Fa(1108)+EpsF(45)*Fa(1113)+EpsF(54)*Fa(1118)-EpsF(93)*
     -   Fa(1120)+EpsF(116)*Fa(1121)+EpsF(4)*Fa(1166)+EpsF(6)*Fa(1167
     -   )-EpsF(7)*Fa(1168)+EpsF(9)*Fa(1169)+EpsF(10)*Fa(1170)+EpsF(1
     -   1)*Fa(1171)+EpsF(13)*Fa(1172)-EpsF(37)*Fa(1173)-EpsF(38)*Fa(
     -   1174)-EpsF(40)*Fa(1175)-EpsF(41)*Fa(1176)+EpsF(43)*Fa(1177)+
     -   EpsF(44)*Fa(1178)+EpsF(46)*Fa(1179)+EpsF(47)*Fa(1180)+EpsF(5
     -   2)*Fa(1181)+EpsF(53)*Fa(1182)+EpsF(55)*Fa(1183)-2*(EpsF(95)*
     -   Fa(1271)+EpsF(96)*Fa(1290)-EpsF(99)*Fa(1345)-EpsF(118)*Fa(13
     -   66)+EpsF(119)*Fa(1412)+EpsF(103)*Fa(1419)-EpsF(107)*Fa(1438)
     -   -EpsF(108)*Fa(1441))+EpsF(2)*Fa(3442)+EpsF(97)*Fa(3445)+EpsF
     -   (98)*Fa(3453)+EpsF(102)*Fa(3455)+EpsF(105)*Fa(3456)+EpsF(112
     -   )*Fa(3458)+EpsF(114)*Fa(3459)+EpsF(124)*Fa(3460)+EpsF(1)*Fa(
     -   3493)+EpsF(3)*Fa(3494)+EpsF(56)*Fa(3495)+EpsF(49)*Fa(3496)+E
     -   psF(50)*Fa(3497)+EpsF(58)*Fa(3498)+EpsF(59)*Fa(3499)+EpsF(61
     -   )*Fa(3500)+EpsF(62)*Fa(3501)+EpsF(94)*Fa(3502)+EpsF(101)*Fa(
     -   3503)-EpsF(106)*Fa(3504)))+M4*Cpl(10)*(EpsF(36)*F(3249)-2*(E
     -   psF(6)*Fa(1098)-EpsF(9)*Fa(1101)+EpsF(5)*Fa(1123)-EpsF(8)*Fa
     -   (1125)-EpsF(4)*Fa(1184)-EpsF(7)*Fa(1185)+EpsF(10)*Fa(1186)+E
     -   psF(11)*Fa(1187)+EpsF(12)*Fa(1188)+EpsF(13)*Fa(1189)-EpsF(37
     -   )*Fa(1190)-EpsF(38)*Fa(1191)-EpsF(39)*Fa(1192)-EpsF(40)*Fa(1
     -   193)-EpsF(41)*Fa(1194)+EpsF(43)*Fa(1195)+EpsF(44)*Fa(1196)+E
     -   psF(45)*Fa(1197)+EpsF(46)*Fa(1198)+EpsF(47)*Fa(1199)-EpsF(52
     -   )*Fa(1200)-EpsF(53)*Fa(1201)+EpsF(54)*Fa(1202)+EpsF(55)*Fa(1
     -   203)-EpsF(93)*Fa(1204)+EpsF(116)*Fa(1205)-2*(EpsF(96)*Fa(127
     -   2)+EpsF(95)*Fa(1279)-EpsF(118)*Fa(1350)-EpsF(99)*Fa(1353)+Ep
     -   sF(103)*Fa(1411)+EpsF(119)*Fa(1415)-EpsF(108)*Fa(1439)-EpsF(
     -   107)*Fa(1440))+EpsF(1)*Fa(3505)-EpsF(2)*Fa(3506)-EpsF(3)*Fa(
     -   3507)-EpsF(56)*Fa(3508)-EpsF(97)*Fa(3509)-EpsF(49)*Fa(3510)-
     -   EpsF(50)*Fa(3511)-EpsF(58)*Fa(3512)+EpsF(59)*Fa(3513)-EpsF(6
     -   1)*Fa(3514)-EpsF(62)*Fa(3515)-EpsF(94)*Fa(3516)-EpsF(98)*Fa(
     -   3517)+EpsF(101)*Fa(3518)-EpsF(102)*Fa(3519)-EpsF(105)*Fa(352
     -   0)-EpsF(106)*Fa(3521)-EpsF(112)*Fa(3522)+EpsF(114)*Fa(3523)-
     -   EpsF(124)*Fa(3524))))+Cpl(9)*(2*M1*M4*Cpl(5)*(EpsF(36)*F(311
     -   3)-2*(EpsF(6)*Fa(1098)-EpsF(9)*Fa(1101)-EpsF(5)*Fa(1141)-Eps
     -   F(8)*Fa(1144)+EpsF(12)*Fa(1148)-EpsF(39)*Fa(1152)+EpsF(45)*F
     -   a(1157)+EpsF(54)*Fa(1162)-EpsF(93)*Fa(1164)+EpsF(116)*Fa(116
     -   5)+EpsF(13)*Fa(1189)-EpsF(40)*Fa(1193)-EpsF(41)*Fa(1194)+Eps
     -   F(46)*Fa(1198)+EpsF(47)*Fa(1199)+EpsF(55)*Fa(1203)+EpsF(4)*F
     -   a(1206)+EpsF(7)*Fa(1207)-EpsF(10)*Fa(1208)-EpsF(11)*Fa(1209)
     -   -EpsF(37)*Fa(1210)-EpsF(38)*Fa(1211)+EpsF(43)*Fa(1212)+EpsF(
     -   44)*Fa(1213)-EpsF(52)*Fa(1214)-EpsF(53)*Fa(1215)-2*(EpsF(95)
     -   *Fa(1271)+EpsF(96)*Fa(1272)-EpsF(99)*Fa(1345)-EpsF(118)*Fa(1
     -   350)+EpsF(103)*Fa(1411)+EpsF(119)*Fa(1412)-EpsF(107)*Fa(1438
     -   )-EpsF(108)*Fa(1439))-EpsF(2)*Fa(3474)-EpsF(97)*Fa(3477)-Eps
     -   F(98)*Fa(3485)-EpsF(102)*Fa(3487)-EpsF(105)*Fa(3488)-EpsF(11
     -   2)*Fa(3490)-EpsF(114)*Fa(3491)-EpsF(124)*Fa(3492)-EpsF(3)*Fa
     -   (3507)-EpsF(56)*Fa(3508)-EpsF(49)*Fa(3510)-EpsF(50)*Fa(3511)
     -   -EpsF(58)*Fa(3512)+EpsF(59)*Fa(3513)-EpsF(61)*Fa(3514)-EpsF(
     -   62)*Fa(3515)-EpsF(1)*Fa(3525)-EpsF(94)*Fa(3526)+EpsF(101)*Fa
     -   (3527)-EpsF(106)*Fa(3528)))+Cpl(10)*(EpsF(36)*F(3341)+4*EpsF
     -   (100)*Fa(3422)+2*(EpsF(4)*Fa(1216)+EpsF(5)*Fa(1217)+EpsF(6)*
     -   Fa(1218)+EpsF(7)*Fa(1219)+EpsF(8)*Fa(1220)+EpsF(9)*Fa(1221)+
     -   EpsF(10)*Fa(1222)+EpsF(11)*Fa(1223)+EpsF(12)*Fa(1224)+EpsF(1
     -   3)*Fa(1225)+EpsF(37)*Fa(1226)+EpsF(38)*Fa(1227)+EpsF(39)*Fa(
     -   1228)+EpsF(40)*Fa(1229)+EpsF(41)*Fa(1230)+EpsF(42)*Fa(1231)+
     -   EpsF(43)*Fa(1232)+EpsF(44)*Fa(1233)+EpsF(45)*Fa(1234)+EpsF(4
     -   6)*Fa(1235)+EpsF(47)*Fa(1236)+EpsF(48)*Fa(1237)+EpsF(52)*Fa(
     -   1238)+EpsF(53)*Fa(1239)+EpsF(54)*Fa(1240)+EpsF(55)*Fa(1241)+
     -   EpsF(93)*Fa(1242)+EpsF(116)*Fa(1243)+EpsF(1)*Fa(3529)+EpsF(2
     -   )*Fa(3530)+EpsF(3)*Fa(3531)+EpsF(56)*Fa(3532)+EpsF(57)*Fa(35
     -   33)+EpsF(97)*Fa(3536)+EpsF(113)*Fa(3537)+EpsF(117)*Fa(3538)+
     -   EpsF(49)*Fa(3539)+EpsF(50)*Fa(3540)+EpsF(51)*Fa(3541)+EpsF(5
     -   8)*Fa(3542)+EpsF(59)*Fa(3543)+EpsF(60)*Fa(3544)+EpsF(61)*Fa(
     -   3545)+EpsF(62)*Fa(3546)+EpsF(63)*Fa(3547)+EpsF(94)*Fa(3548)+
     -   EpsF(98)*Fa(3549)+EpsF(101)*Fa(3551)+EpsF(102)*Fa(3552)+EpsF
     -   (105)*Fa(3553)+EpsF(106)*Fa(3554)+EpsF(112)*Fa(3555)+EpsF(11
     -   4)*Fa(3556)+EpsF(124)*Fa(3559)+2*(EpsF(95)*Fa(3534)+EpsF(96)
     -   *Fa(3535)-EpsF(99)*Fa(3550)+EpsF(115)*Fa(3557)-EpsF(118)*Fa(
     -   3558)-EpsF(103)*Fa(3560)-EpsF(104)*Fa(3561)+EpsF(119)*Fa(356
     -   2)+EpsF(123)*Fa(3563)+EpsF(107)*Fa(3564)+EpsF(108)*Fa(3565)+
     -   EpsF(109)*Fa(3566))+EpsF(126)*Fa(3567)+EpsF(120)*Fa(3568))))
     -   ))
       result(1,3) = (0,0.5)*(Cpl(3)*(Cpl(4)*(Cpl(5)*(EpsF(69)*F(4617)
     -   +2*(EpsF(14)*Fa(1456)+EpsF(15)*Fa(1457)+EpsF(16)*Fa(1458)+Ep
     -   sF(17)*Fa(1459)+EpsF(18)*Fa(1460)+EpsF(19)*Fa(1461)+EpsF(21)
     -   *Fa(1462)+EpsF(22)*Fa(1463)+EpsF(23)*Fa(1464)+EpsF(24)*Fa(14
     -   65)+EpsF(52)*Fa(1466)+EpsF(53)*Fa(1467)+EpsF(54)*Fa(1468)+Ep
     -   sF(55)*Fa(1469)+EpsF(70)*Fa(1470)+EpsF(71)*Fa(1471)+EpsF(72)
     -   *Fa(1472)+EpsF(73)*Fa(1473)+EpsF(74)*Fa(1474)+EpsF(75)*Fa(14
     -   75)+EpsF(77)*Fa(1476)+EpsF(78)*Fa(1477)+EpsF(79)*Fa(1478)+Ep
     -   sF(80)*Fa(1479)+EpsF(81)*Fa(1480)+EpsF(82)*Fa(1481)+EpsF(83)
     -   *Fa(1482)+EpsF(85)*Fa(1483)+EpsF(1)*Fa(1611)+EpsF(2)*Fa(1612
     -   )+EpsF(97)*Fa(1618)+EpsF(86)*Fa(3569)+EpsF(87)*Fa(3570)+EpsF
     -   (88)*Fa(3571)+EpsF(90)*Fa(3572)+EpsF(91)*Fa(3573)+EpsF(92)*F
     -   a(3574)+EpsF(3)*Fa(3575)+EpsF(56)*Fa(3576)+EpsF(57)*Fa(3577)
     -   +EpsF(58)*Fa(3578)+EpsF(59)*Fa(3579)+EpsF(60)*Fa(3580)+EpsF(
     -   61)*Fa(3581)+EpsF(62)*Fa(3582)+EpsF(63)*Fa(3583)+EpsF(98)*Fa
     -   (3584)+EpsF(101)*Fa(3586)+EpsF(102)*Fa(3587)+EpsF(106)*Fa(35
     -   88)+EpsF(114)*Fa(3589)+EpsF(124)*Fa(3592)+4*(EpsF(111)*Fa(15
     -   93)-EpsF(122)*Fa(1594)-EpsF(125)*Fa(1595)+EpsF(100)*Fa(3585)
     -   -EpsF(115)*Fa(3590)+EpsF(118)*Fa(3591)-EpsF(103)*Fa(3593)-Ep
     -   sF(104)*Fa(3594)+EpsF(123)*Fa(3595)+EpsF(108)*Fa(3596)+EpsF(
     -   109)*Fa(3597)-EpsF(126)*Fa(3598))))+2*M1*M4*Cpl(10)*(EpsF(69
     -   )*F(4749)+2*(mup3mup5*EpsF(14)*F(4730)-mup4mup5*EpsF(17)*F(4
     -   730)-mup3mup4*EpsF(86)*F(4730)-mup1mup3*EpsF(101)*F(4730)+mu
     -   p1mup4*EpsF(106)*F(4730)-mup3mup5*EpsF(15)*F(4731)+mup4mup5*
     -   EpsF(18)*F(4731)+mup3mup4*EpsF(88)*F(4731)-mup1mup4*EpsF(114
     -   )*F(4731)+mup1mup3*EpsF(124)*F(4731)+mup4mup5*EpsF(71)*F(474
     -   5)+mup3mup4*EpsF(87)*F(4745)+mup1mup3*EpsF(102)*F(4745)+mup1
     -   mup3*EpsF(61)*F(4746)+mup4mup5*EpsF(75)*F(4746)+mup3mup4*Eps
     -   F(90)*F(4746)+mup1mup5*(-(EpsF(1)*F(4730))+EpsF(2)*F(4731)+E
     -   psF(97)*F(4745)+EpsF(56)*F(4746))+mup3mup5*EpsF(80)*F(4747)+
     -   mup1mup4*EpsF(98)*F(4747)+(mup1mup4*EpsF(58)+mup3mup5*EpsF(8
     -   2))*F(4748)+EpsF(21)*Fa(1484)+EpsF(22)*Fa(1485)+EpsF(23)*Fa(
     -   1486)+EpsF(24)*Fa(1487)+EpsF(52)*Fa(1488)+EpsF(53)*Fa(1489)+
     -   EpsF(54)*Fa(1490)+EpsF(55)*Fa(1491)+EpsF(70)*Fa(1492)-EpsF(7
     -   3)*Fa(1493)+EpsF(74)*Fa(1494)+EpsF(77)*Fa(1495)+EpsF(78)*Fa(
     -   1496)+EpsF(79)*Fa(1497)+EpsF(81)*Fa(1498)+EpsF(85)*Fa(1499)+
     -   (mup1mup5*EpsF(3)-mup3mup5*EpsF(16)+mup4mup5*EpsF(19)-mup1mu
     -   p4*EpsF(59)+mup1mup3*EpsF(62)+mup3mup4*EpsF(91))*F(7439))))+
     -   2*M3*Cpl(9)*(M1*Cpl(10)*(EpsF(69)*F(4788)+2*(mup1mup3*EpsF(1
     -   01)*F(4730)-mup1mup4*EpsF(106)*F(4730)+mup1mup4*EpsF(114)*F(
     -   4731)-mup1mup3*EpsF(124)*F(4731)+mup1mup4*EpsF(98)*F(4745)+m
     -   up1mup3*EpsF(102)*F(4747)+mup1mup4*EpsF(59)*F(4779)-mup1mup3
     -   *EpsF(62)*F(4779)+mup1mup3*EpsF(61)*F(4786)+mup1mup5*(EpsF(1
     -   )*F(4730)-EpsF(2)*F(4731)+EpsF(97)*F(4747)-EpsF(3)*F(4779)+E
     -   psF(56)*F(4786))+mup4mup5*(EpsF(17)*F(4730)-EpsF(18)*F(4731)
     -   +EpsF(71)*F(4747)-EpsF(19)*F(4779)+EpsF(75)*F(4786))+mup3mup
     -   4*(EpsF(86)*F(4730)-EpsF(88)*F(4731)+EpsF(87)*F(4747)-EpsF(9
     -   1)*F(4779)+EpsF(90)*F(4786))+mup1mup4*EpsF(58)*F(4787)+mup3m
     -   up5*(-(EpsF(14)*F(4730))+EpsF(15)*F(4731)+EpsF(80)*F(4745)+E
     -   psF(16)*F(4779)+EpsF(82)*F(4787))+EpsF(21)*Fa(1512)-EpsF(22)
     -   *Fa(1513)+EpsF(23)*Fa(1514)+EpsF(24)*Fa(1515)+EpsF(52)*Fa(15
     -   16)-EpsF(53)*Fa(1517)+EpsF(54)*Fa(1518)+EpsF(55)*Fa(1519)+Ep
     -   sF(70)*Fa(1520)-EpsF(73)*Fa(1521)+EpsF(74)*Fa(1522)+EpsF(77)
     -   *Fa(1523)+EpsF(78)*Fa(1524)+EpsF(79)*Fa(1525)+EpsF(81)*Fa(15
     -   26)+EpsF(85)*Fa(1527)))+M4*Cpl(5)*(EpsF(69)*F(4775)+2*((mup1
     -   mup5*EpsF(56)+mup1mup3*EpsF(61)+mup4mup5*EpsF(75)+mup3mup4*E
     -   psF(90))*F(4746)+(mup1mup4*EpsF(58)+mup3mup5*EpsF(82))*F(474
     -   8)+(mup4mup5*EpsF(71)+mup3mup4*EpsF(87)+mup1mup5*EpsF(97)+mu
     -   p1mup3*EpsF(102))*F(4773)+(mup3mup5*EpsF(80)+mup1mup4*EpsF(9
     -   8))*F(4774)+EpsF(24)*Fa(1487)+EpsF(55)*Fa(1491)+EpsF(74)*Fa(
     -   1494)+EpsF(81)*Fa(1498)+EpsF(21)*Fa(1500)-EpsF(22)*Fa(1501)+
     -   EpsF(23)*Fa(1502)+EpsF(52)*Fa(1503)+EpsF(53)*Fa(1504)+EpsF(5
     -   4)*Fa(1505)+EpsF(70)*Fa(1506)+EpsF(73)*Fa(1507)+EpsF(77)*Fa(
     -   1508)-EpsF(78)*Fa(1509)+EpsF(79)*Fa(1510)+EpsF(85)*Fa(1511)-
     -   (mup1mup5*EpsF(1)-mup3mup5*EpsF(14)+mup4mup5*EpsF(17)+mup3mu
     -   p4*EpsF(86)+mup1mup3*EpsF(101)-mup1mup4*EpsF(106))*F(7437)+(
     -   mup1mup5*EpsF(2)-mup3mup5*EpsF(15)+mup4mup5*EpsF(18)+mup3mup
     -   4*EpsF(88)-mup1mup4*EpsF(114)+mup1mup3*EpsF(124))*F(7438)+(m
     -   up1mup5*EpsF(3)-mup3mup5*EpsF(16)+mup4mup5*EpsF(19)-mup1mup4
     -   *EpsF(59)+mup1mup3*EpsF(62)+mup3mup4*EpsF(91))*F(7439)))))+C
     -   pl(8)*(Cpl(9)*(Cpl(10)*(EpsF(69)*F(4829)+2*(EpsF(14)*Fa(1559
     -   )+EpsF(15)*Fa(1560)+EpsF(16)*Fa(1561)+EpsF(17)*Fa(1562)+EpsF
     -   (18)*Fa(1563)+EpsF(19)*Fa(1564)+EpsF(21)*Fa(1565)+EpsF(22)*F
     -   a(1566)+EpsF(23)*Fa(1567)+EpsF(24)*Fa(1568)+EpsF(52)*Fa(1569
     -   )+EpsF(53)*Fa(1570)+EpsF(54)*Fa(1571)+EpsF(55)*Fa(1572)+EpsF
     -   (70)*Fa(1573)+EpsF(71)*Fa(1574)+EpsF(72)*Fa(1575)+EpsF(73)*F
     -   a(1576)+EpsF(74)*Fa(1577)+EpsF(75)*Fa(1578)+EpsF(77)*Fa(1579
     -   )+EpsF(78)*Fa(1580)+EpsF(79)*Fa(1581)+EpsF(80)*Fa(1582)+EpsF
     -   (81)*Fa(1583)+EpsF(82)*Fa(1584)+EpsF(83)*Fa(1585)+EpsF(85)*F
     -   a(1586)+EpsF(1)*Fa(1628)+EpsF(2)*Fa(1629)+EpsF(97)*Fa(1632)+
     -   EpsF(86)*Fa(3599)+EpsF(87)*Fa(3600)+EpsF(88)*Fa(3601)+EpsF(9
     -   0)*Fa(3602)+EpsF(91)*Fa(3603)+EpsF(92)*Fa(3604)+EpsF(3)*Fa(3
     -   605)+EpsF(56)*Fa(3606)+EpsF(57)*Fa(3607)+EpsF(58)*Fa(3608)+E
     -   psF(59)*Fa(3609)+EpsF(60)*Fa(3610)+EpsF(61)*Fa(3611)+EpsF(62
     -   )*Fa(3612)+EpsF(63)*Fa(3613)+EpsF(98)*Fa(3614)+EpsF(101)*Fa(
     -   3615)+EpsF(102)*Fa(3616)+EpsF(106)*Fa(3617)+EpsF(114)*Fa(361
     -   8)+EpsF(124)*Fa(3619)-4*(EpsF(111)*Fa(1593)-EpsF(122)*Fa(159
     -   4)-EpsF(125)*Fa(1595)+EpsF(100)*Fa(3585)-EpsF(115)*Fa(3590)+
     -   EpsF(118)*Fa(3591)-EpsF(104)*Fa(3594)+EpsF(123)*Fa(3595)+Eps
     -   F(108)*Fa(3596)+EpsF(109)*Fa(3597)-EpsF(126)*Fa(3598)-EpsF(1
     -   03)*Fa(3620))))+2*M1*M4*Cpl(5)*(EpsF(69)*F(4800)-2*(-(mup1mu
     -   p5*EpsF(1)*F(4730))+mup3mup5*EpsF(14)*F(4730)-mup4mup5*EpsF(
     -   17)*F(4730)-mup3mup4*EpsF(86)*F(4730)-mup1mup3*EpsF(101)*F(4
     -   730)+mup1mup4*EpsF(106)*F(4730)+mup1mup5*EpsF(2)*F(4731)-mup
     -   3mup5*EpsF(15)*F(4731)+mup4mup5*EpsF(18)*F(4731)+mup3mup4*Ep
     -   sF(88)*F(4731)-mup1mup4*EpsF(114)*F(4731)+mup1mup3*EpsF(124)
     -   *F(4731)-mup3mup5*EpsF(80)*F(4745)-mup1mup4*EpsF(98)*F(4745)
     -   -mup1mup4*EpsF(58)*F(4746)-mup3mup5*EpsF(82)*F(4746)-mup4mup
     -   5*EpsF(71)*F(4747)-mup3mup4*EpsF(87)*F(4747)-mup1mup5*EpsF(9
     -   7)*F(4747)-mup1mup3*EpsF(102)*F(4747)-(mup1mup5*EpsF(56)+mup
     -   1mup3*EpsF(61)+mup4mup5*EpsF(75)+mup3mup4*EpsF(90))*F(4748)+
     -   EpsF(53)*Fa(1489)-EpsF(23)*Fa(1514)-EpsF(54)*Fa(1518)-EpsF(7
     -   0)*Fa(1520)-EpsF(79)*Fa(1525)-EpsF(24)*Fa(1540)-EpsF(55)*Fa(
     -   1543)-EpsF(74)*Fa(1546)-EpsF(81)*Fa(1550)-EpsF(21)*Fa(1552)+
     -   EpsF(22)*Fa(1553)-EpsF(52)*Fa(1554)-EpsF(73)*Fa(1555)-EpsF(7
     -   7)*Fa(1556)-EpsF(78)*Fa(1557)-EpsF(85)*Fa(1558)+(mup1mup5*Ep
     -   sF(3)-mup3mup5*EpsF(16)+mup4mup5*EpsF(19)-mup1mup4*EpsF(59)+
     -   mup1mup3*EpsF(62)+mup3mup4*EpsF(91))*F(7439))))+2*M3*Cpl(4)*
     -   (M1*Cpl(5)*(EpsF(69)*F(4798)+2*(-(mup1mup5*EpsF(1)*F(4730))+
     -   mup3mup5*EpsF(14)*F(4730)-mup4mup5*EpsF(17)*F(4730)-mup3mup4
     -   *EpsF(86)*F(4730)-mup1mup3*EpsF(101)*F(4730)+mup1mup4*EpsF(1
     -   06)*F(4730)+mup1mup5*EpsF(2)*F(4731)-mup3mup5*EpsF(15)*F(473
     -   1)+mup4mup5*EpsF(18)*F(4731)+mup3mup4*EpsF(88)*F(4731)-mup1m
     -   up4*EpsF(114)*F(4731)+mup1mup3*EpsF(124)*F(4731)+mup4mup5*Ep
     -   sF(71)*F(4745)+mup3mup4*EpsF(87)*F(4745)+mup1mup5*EpsF(97)*F
     -   (4745)+mup1mup3*EpsF(102)*F(4745)+mup3mup5*EpsF(80)*F(4747)+
     -   mup1mup4*EpsF(98)*F(4747)+mup1mup5*EpsF(3)*F(4779)-mup3mup5*
     -   EpsF(16)*F(4779)+mup4mup5*EpsF(19)*F(4779)-mup1mup4*EpsF(59)
     -   *F(4779)+mup1mup3*EpsF(62)*F(4779)+mup3mup4*EpsF(91)*F(4779)
     -   +mup1mup4*EpsF(58)*F(4786)+mup3mup5*EpsF(82)*F(4786)+(mup1mu
     -   p5*EpsF(56)+mup1mup3*EpsF(61)+mup4mup5*EpsF(75)+mup3mup4*Eps
     -   F(90))*F(4787)+EpsF(23)*Fa(1486)+EpsF(54)*Fa(1490)+EpsF(70)*
     -   Fa(1492)+EpsF(79)*Fa(1497)+EpsF(53)*Fa(1517)+EpsF(73)*Fa(152
     -   1)-EpsF(78)*Fa(1524)+EpsF(21)*Fa(1528)+EpsF(22)*Fa(1529)+Eps
     -   F(24)*Fa(1530)+EpsF(52)*Fa(1531)+EpsF(55)*Fa(1532)+EpsF(74)*
     -   Fa(1533)+EpsF(77)*Fa(1534)+EpsF(81)*Fa(1535)+EpsF(85)*Fa(153
     -   6)))+M4*Cpl(10)*(EpsF(69)*F(4799)-2*(-((mup1mup4*EpsF(58)+mu
     -   p3mup5*EpsF(82))*F(4746))-(mup1mup5*EpsF(56)+mup1mup3*EpsF(6
     -   1)+mup4mup5*EpsF(75)+mup3mup4*EpsF(90))*F(4748)-(mup3mup5*Ep
     -   sF(80)+mup1mup4*EpsF(98))*F(4773)-(mup4mup5*EpsF(71)+mup3mup
     -   4*EpsF(87)+mup1mup5*EpsF(97)+mup1mup3*EpsF(102))*F(4774)+Eps
     -   F(53)*Fa(1504)-EpsF(21)*Fa(1537)-EpsF(22)*Fa(1538)-EpsF(23)*
     -   Fa(1539)-EpsF(24)*Fa(1540)-EpsF(52)*Fa(1541)-EpsF(54)*Fa(154
     -   2)-EpsF(55)*Fa(1543)-EpsF(70)*Fa(1544)+EpsF(73)*Fa(1545)-Eps
     -   F(74)*Fa(1546)-EpsF(77)*Fa(1547)-EpsF(78)*Fa(1548)-EpsF(79)*
     -   Fa(1549)-EpsF(81)*Fa(1550)-EpsF(85)*Fa(1551)-(mup1mup5*EpsF(
     -   1)-mup3mup5*EpsF(14)+mup4mup5*EpsF(17)+mup3mup4*EpsF(86)+mup
     -   1mup3*EpsF(101)-mup1mup4*EpsF(106))*F(7437)+(mup1mup5*EpsF(2
     -   )-mup3mup5*EpsF(15)+mup4mup5*EpsF(18)+mup3mup4*EpsF(88)-mup1
     -   mup4*EpsF(114)+mup1mup3*EpsF(124))*F(7438)+(mup1mup5*EpsF(3)
     -   -mup3mup5*EpsF(16)+mup4mup5*EpsF(19)-mup1mup4*EpsF(59)+mup1m
     -   up3*EpsF(62)+mup3mup4*EpsF(91))*F(7439))))))
       result(1,4) = (0,0.5)*(Cpl(3)*(Cpl(4)*(Cpl(5)*(EpsF(76)*F(5082)
     -   +4*(EpsF(110)*Fa(1845)-EpsF(111)*Fa(1846)-EpsF(125)*Fa(1848)
     -   +EpsF(120)*F(5708)*Fa(3086)+EpsF(123)*Fa(3658)-EpsF(109)*Fa(
     -   3661))+2*(EpsF(25)*Fa(1649)+EpsF(26)*Fa(1650)+EpsF(27)*Fa(16
     -   51)+EpsF(28)*Fa(1652)+EpsF(29)*Fa(1653)+EpsF(30)*Fa(1654)+Ep
     -   sF(31)*Fa(1655)+EpsF(43)*Fa(1656)+EpsF(44)*Fa(1657)+EpsF(45)
     -   *Fa(1658)+EpsF(46)*Fa(1659)+EpsF(65)*Fa(1660)+EpsF(66)*Fa(16
     -   61)+EpsF(67)*Fa(1662)+EpsF(68)*Fa(1663)+EpsF(77)*Fa(1664)+Ep
     -   sF(78)*Fa(1665)+EpsF(79)*Fa(1666)+EpsF(81)*Fa(1667)+EpsF(84)
     -   *Fa(1668)+EpsF(89)*Fa(1669)+EpsF(121)*Fa(1670)+EpsF(122)*Fa(
     -   1847)+EpsF(14)*Fa(3621)+EpsF(15)*Fa(3622)+EpsF(16)*Fa(3623)+
     -   EpsF(80)*Fa(3624)+EpsF(82)*Fa(3625)+EpsF(83)*Fa(3626)+EpsF(8
     -   6)*Fa(3627)+EpsF(87)*Fa(3628)+EpsF(88)*Fa(3629)+EpsF(90)*Fa(
     -   3630)+EpsF(91)*Fa(3631)+EpsF(92)*Fa(3632)+EpsF(4)*Fa(3633)+E
     -   psF(5)*Fa(3634)+EpsF(6)*Fa(3635)+EpsF(47)*Fa(3636)+EpsF(48)*
     -   Fa(3637)+EpsF(95)*Fa(3638)+EpsF(96)*Fa(3639)+EpsF(113)*Fa(36
     -   40)+EpsF(116)*Fa(3641)+EpsF(117)*Fa(3642)+EpsF(49)*Fa(3643)+
     -   EpsF(50)*Fa(3644)+EpsF(51)*Fa(3645)+EpsF(61)*Fa(3646)+EpsF(6
     -   2)*Fa(3647)+EpsF(63)*Fa(3648)+EpsF(94)*Fa(3649)+EpsF(101)*Fa
     -   (3650)+EpsF(102)*Fa(3651)+EpsF(105)*Fa(3652)+EpsF(112)*Fa(36
     -   53)+EpsF(124)*Fa(3654)+EpsF(103)*Fa(3655)+EpsF(104)*Fa(3656)
     -   +EpsF(119)*Fa(3657)+EpsF(107)*Fa(3659)+EpsF(108)*Fa(3660)+Ep
     -   sF(126)*Fa(3662)))+M1*M4*Cpl(10)*(EpsF(76)*F(5151)+4*(EpsF(1
     -   11)*Fa(210)+EpsF(125)*Fa(215))+2*(EpsF(68)*Fa(28)+EpsF(25)*F
     -   a(1671)+EpsF(26)*Fa(1672)+EpsF(27)*Fa(1673)+EpsF(28)*Fa(1674
     -   )+EpsF(29)*Fa(1675)-EpsF(30)*Fa(1676)+EpsF(31)*Fa(1677)+EpsF
     -   (43)*Fa(1678)+EpsF(44)*Fa(1679)+EpsF(45)*Fa(1680)+EpsF(46)*F
     -   a(1681)+EpsF(65)*Fa(1682)+EpsF(66)*Fa(1683)+EpsF(67)*Fa(1684
     -   )+EpsF(77)*Fa(1685)+EpsF(78)*Fa(1686)+EpsF(79)*Fa(1687)+EpsF
     -   (81)*Fa(1688)+EpsF(84)*Fa(1689)+EpsF(89)*Fa(1690)+EpsF(63)*F
     -   a(3079)+EpsF(83)*Fa(3226)+EpsF(92)*Fa(3233)+EpsF(48)*Fa(3245
     -   )+EpsF(51)*Fa(3249)+EpsF(14)*Fa(3663)+EpsF(15)*Fa(3664)+EpsF
     -   (16)*Fa(3665)+EpsF(80)*Fa(3666)+EpsF(82)*Fa(3667)+EpsF(86)*F
     -   a(3668)+EpsF(87)*Fa(3669)+EpsF(88)*Fa(3670)+EpsF(90)*Fa(3671
     -   )+EpsF(91)*Fa(3672)+EpsF(4)*Fa(3673)+EpsF(5)*Fa(3674)+EpsF(6
     -   )*Fa(3675)+EpsF(116)*Fa(3679)+EpsF(50)*Fa(3681)+EpsF(61)*Fa(
     -   3682)+EpsF(62)*Fa(3683)+EpsF(94)*Fa(3684)+EpsF(101)*Fa(3685)
     -   +EpsF(102)*Fa(3686)+EpsF(105)*Fa(3687)+EpsF(112)*Fa(3688)+Ep
     -   sF(124)*Fa(3689)+2*(EpsF(110)*Fa(1854)-EpsF(123)*Fa(3051)+Ep
     -   sF(113)*Fa(3056)+EpsF(104)*Fa(3067)+EpsF(109)*Fa(3080)+EpsF(
     -   126)*Fa(3081)-EpsF(95)*Fa(3677)+EpsF(119)*Fa(3691)+EpsF(107)
     -   *Fa(3692))+(EpsF(121)*Fa(1691)+EpsF(47)*Fa(3676)+EpsF(49)*Fa
     -   (3680)-2*(EpsF(122)*Fa(1855)+EpsF(108)*Fa(1954)+EpsF(117)*Fa
     -   (3065)+EpsF(120)*Fa(3086)-EpsF(96)*Fa(3678)-EpsF(103)*Fa(369
     -   0)))*F(7441))))+M3*Cpl(9)*(M1*Cpl(10)*(EpsF(76)*F(5258)-4*(E
     -   psF(111)*Fa(231)+EpsF(125)*Fa(234)-EpsF(110)*Fa(1861)-EpsF(1
     -   23)*Fa(3152)-EpsF(104)*Fa(3167)+EpsF(109)*Fa(3178)-EpsF(126)
     -   *Fa(3179)-EpsF(113)*Fa(3205))+2*(EpsF(68)*Fa(52)+EpsF(25)*Fa
     -   (1712)-EpsF(26)*Fa(1713)+EpsF(27)*Fa(1714)+EpsF(28)*Fa(1715)
     -   +EpsF(29)*Fa(1716)+EpsF(30)*Fa(1717)+EpsF(31)*Fa(1718)+EpsF(
     -   43)*Fa(1719)+EpsF(44)*Fa(1720)+EpsF(45)*Fa(1721)+EpsF(46)*Fa
     -   (1722)+EpsF(65)*Fa(1723)+EpsF(66)*Fa(1724)+EpsF(67)*Fa(1725)
     -   +EpsF(77)*Fa(1726)+EpsF(78)*Fa(1727)+EpsF(79)*Fa(1728)+EpsF(
     -   81)*Fa(1729)+EpsF(84)*Fa(1730)+EpsF(89)*Fa(1731)+EpsF(63)*Fa
     -   (3177)+EpsF(83)*Fa(3190)+EpsF(92)*Fa(3195)+EpsF(48)*Fa(3208)
     -   +EpsF(51)*Fa(3211)+EpsF(14)*Fa(3713)+EpsF(15)*Fa(3714)+EpsF(
     -   16)*Fa(3715)+EpsF(80)*Fa(3716)+EpsF(82)*Fa(3717)+EpsF(86)*Fa
     -   (3718)+EpsF(87)*Fa(3719)+EpsF(88)*Fa(3720)+EpsF(90)*Fa(3721)
     -   +EpsF(91)*Fa(3722)+EpsF(4)*Fa(3723)-EpsF(5)*Fa(3724)+EpsF(6)
     -   *Fa(3725)+EpsF(116)*Fa(3728)+F(1095)*(-(EpsF(121)*Fa(1691))+
     -   2*(EpsF(122)*Fa(1855)+EpsF(108)*Fa(1954)+EpsF(117)*Fa(3065)+
     -   EpsF(120)*Fa(3086)-EpsF(96)*Fa(3678)-EpsF(103)*Fa(3690))+Eps
     -   F(47)*Fa(3726)+EpsF(49)*Fa(3729))+EpsF(50)*Fa(3730)+EpsF(61)
     -   *Fa(3731)+EpsF(62)*Fa(3732)+EpsF(94)*Fa(3733)+EpsF(101)*Fa(3
     -   734)+EpsF(102)*Fa(3735)+EpsF(105)*Fa(3736)+EpsF(112)*Fa(3737
     -   )+EpsF(124)*Fa(3738)+2*(EpsF(95)*Fa(3727)-EpsF(119)*Fa(3739)
     -   +EpsF(107)*Fa(3740))))+M4*Cpl(5)*(EpsF(76)*F(5203)+4*(EpsF(1
     -   11)*Fa(210)+EpsF(125)*Fa(215))+2*(EpsF(68)*Fa(28)+EpsF(25)*F
     -   a(1692)+EpsF(26)*Fa(1693)+EpsF(27)*Fa(1694)+EpsF(28)*Fa(1695
     -   )+EpsF(29)*Fa(1696)+EpsF(30)*Fa(1697)+EpsF(31)*Fa(1698)+EpsF
     -   (43)*Fa(1699)+EpsF(44)*Fa(1700)+EpsF(45)*Fa(1701)+EpsF(46)*F
     -   a(1702)+EpsF(65)*Fa(1703)+EpsF(66)*Fa(1704)+EpsF(67)*Fa(1705
     -   )+EpsF(77)*Fa(1706)+EpsF(78)*Fa(1707)+EpsF(79)*Fa(1708)+EpsF
     -   (81)*Fa(1709)+EpsF(84)*Fa(1710)+EpsF(89)*Fa(1711)+EpsF(63)*F
     -   a(3079)+EpsF(83)*Fa(3226)+EpsF(92)*Fa(3233)+EpsF(48)*Fa(3245
     -   )+EpsF(51)*Fa(3249)+EpsF(82)*Fa(3667)+EpsF(90)*Fa(3671)+EpsF
     -   (61)*Fa(3682)+2*(EpsF(110)*Fa(1854)-EpsF(123)*Fa(3051)+EpsF(
     -   113)*Fa(3056)+EpsF(104)*Fa(3067)+EpsF(109)*Fa(3080)+EpsF(126
     -   )*Fa(3081)-EpsF(95)*Fa(3677)+EpsF(119)*Fa(3691)+EpsF(107)*Fa
     -   (3692))+EpsF(14)*Fa(3693)+EpsF(15)*Fa(3694)+EpsF(16)*Fa(3695
     -   )+EpsF(80)*Fa(3696)+EpsF(86)*Fa(3697)+EpsF(87)*Fa(3698)+EpsF
     -   (88)*Fa(3699)+EpsF(91)*Fa(3700)+EpsF(4)*Fa(3701)+EpsF(5)*Fa(
     -   3702)+EpsF(6)*Fa(3703)+EpsF(116)*Fa(3704)+EpsF(50)*Fa(3705)+
     -   EpsF(62)*Fa(3706)+EpsF(94)*Fa(3707)+EpsF(101)*Fa(3708)+EpsF(
     -   102)*Fa(3709)+EpsF(105)*Fa(3710)+EpsF(112)*Fa(3711)+EpsF(124
     -   )*Fa(3712)+(EpsF(121)*Fa(1691)+EpsF(47)*Fa(3676)+EpsF(49)*Fa
     -   (3680)-2*(EpsF(122)*Fa(1855)+EpsF(108)*Fa(1954)+EpsF(117)*Fa
     -   (3065)+EpsF(120)*Fa(3086)-EpsF(96)*Fa(3678)-EpsF(103)*Fa(369
     -   0)))*F(7441)))))+Cpl(8)*(M3*Cpl(4)*(M1*Cpl(5)*(EpsF(76)*F(53
     -   00)+4*(EpsF(111)*Fa(231)+EpsF(125)*Fa(234)+EpsF(110)*Fa(1864
     -   )+EpsF(113)*Fa(3157)-EpsF(104)*Fa(3167)+EpsF(109)*Fa(3178)+E
     -   psF(123)*Fa(3201)+EpsF(126)*Fa(3219))+2*(EpsF(68)*Fa(63)+Eps
     -   F(25)*Fa(1732)+EpsF(26)*Fa(1733)+EpsF(27)*Fa(1734)+EpsF(28)*
     -   Fa(1735)+EpsF(29)*Fa(1736)-EpsF(30)*Fa(1737)+EpsF(31)*Fa(173
     -   8)+EpsF(43)*Fa(1739)+EpsF(44)*Fa(1740)+EpsF(45)*Fa(1741)+Eps
     -   F(46)*Fa(1742)+EpsF(65)*Fa(1743)+EpsF(66)*Fa(1744)+EpsF(67)*
     -   Fa(1745)+EpsF(77)*Fa(1746)+EpsF(78)*Fa(1747)+EpsF(79)*Fa(174
     -   8)+EpsF(81)*Fa(1749)+EpsF(84)*Fa(1750)+EpsF(89)*Fa(1751)+Eps
     -   F(83)*Fa(3137)+EpsF(92)*Fa(3144)+EpsF(48)*Fa(3161)+EpsF(51)*
     -   Fa(3165)+EpsF(63)*Fa(3218)+F(1095)*(EpsF(121)*Fa(1691)+EpsF(
     -   47)*Fa(3676)+EpsF(49)*Fa(3680)-2*(EpsF(122)*Fa(1855)+EpsF(10
     -   8)*Fa(1954)+EpsF(117)*Fa(3065)+EpsF(120)*Fa(3086)-EpsF(96)*F
     -   a(3678)-EpsF(103)*Fa(3690)))+EpsF(14)*Fa(3741)+EpsF(15)*Fa(3
     -   742)+EpsF(16)*Fa(3743)+EpsF(80)*Fa(3744)+EpsF(82)*Fa(3745)+E
     -   psF(86)*Fa(3746)+EpsF(87)*Fa(3747)+EpsF(88)*Fa(3748)+EpsF(90
     -   )*Fa(3749)+EpsF(91)*Fa(3750)+EpsF(4)*Fa(3751)+EpsF(5)*Fa(375
     -   2)+EpsF(6)*Fa(3753)+EpsF(116)*Fa(3755)+EpsF(50)*Fa(3756)+Eps
     -   F(61)*Fa(3757)+EpsF(62)*Fa(3758)+EpsF(94)*Fa(3759)+EpsF(101)
     -   *Fa(3760)+EpsF(102)*Fa(3761)+EpsF(105)*Fa(3762)+EpsF(112)*Fa
     -   (3763)+EpsF(124)*Fa(3764)-2*(EpsF(95)*Fa(3754)-EpsF(119)*Fa(
     -   3765)-EpsF(107)*Fa(3766))))+M4*Cpl(10)*(EpsF(76)*F(5323)-4*(
     -   EpsF(111)*Fa(210)+EpsF(125)*Fa(215))+2*(EpsF(68)*Fa(76)+EpsF
     -   (25)*Fa(1752)+EpsF(26)*Fa(1753)+EpsF(27)*Fa(1754)+EpsF(28)*F
     -   a(1755)+EpsF(29)*Fa(1756)+EpsF(30)*Fa(1757)+EpsF(31)*Fa(1758
     -   )+EpsF(43)*Fa(1759)+EpsF(44)*Fa(1760)+EpsF(45)*Fa(1761)+EpsF
     -   (46)*Fa(1762)+EpsF(65)*Fa(1763)+EpsF(66)*Fa(1764)+EpsF(67)*F
     -   a(1765)+EpsF(77)*Fa(1766)+EpsF(78)*Fa(1767)+EpsF(79)*Fa(1768
     -   )+EpsF(81)*Fa(1769)+EpsF(84)*Fa(1770)+EpsF(89)*Fa(1771)+EpsF
     -   (83)*Fa(3036)+EpsF(92)*Fa(3043)+EpsF(48)*Fa(3060)+EpsF(51)*F
     -   a(3064)+EpsF(63)*Fa(3259)+EpsF(14)*Fa(3767)+EpsF(15)*Fa(3768
     -   )+EpsF(16)*Fa(3769)+EpsF(80)*Fa(3770)+EpsF(82)*Fa(3771)+EpsF
     -   (86)*Fa(3772)+EpsF(87)*Fa(3773)+EpsF(88)*Fa(3774)+EpsF(90)*F
     -   a(3775)+EpsF(91)*Fa(3776)+EpsF(4)*Fa(3777)-EpsF(5)*Fa(3778)+
     -   EpsF(6)*Fa(3779)+EpsF(116)*Fa(3780)+EpsF(50)*Fa(3781)+EpsF(6
     -   1)*Fa(3782)+EpsF(62)*Fa(3783)+EpsF(94)*Fa(3784)+EpsF(101)*Fa
     -   (3785)+EpsF(102)*Fa(3786)+EpsF(105)*Fa(3787)+EpsF(112)*Fa(37
     -   88)+EpsF(124)*Fa(3789)+2*(EpsF(110)*Fa(1866)+EpsF(123)*Fa(30
     -   51)-EpsF(113)*Fa(3056)-EpsF(109)*Fa(3080)-EpsF(126)*Fa(3081)
     -   +EpsF(104)*Fa(3250)+EpsF(95)*Fa(3677)-EpsF(119)*Fa(3790)+Eps
     -   F(107)*Fa(3791))+(-(EpsF(121)*Fa(1691))+2*(EpsF(122)*Fa(1855
     -   )+EpsF(108)*Fa(1954)+EpsF(117)*Fa(3065)+EpsF(120)*Fa(3086)-E
     -   psF(96)*Fa(3678)-EpsF(103)*Fa(3690))+EpsF(47)*Fa(3726)+EpsF(
     -   49)*Fa(3729))*F(7441))))+Cpl(9)*(Cpl(10)*(EpsF(76)*F(5413)+4
     -   *(EpsF(111)*Fa(1846)+EpsF(125)*Fa(1848)-EpsF(120)*F(5708)*Fa
     -   (3086)+EpsF(109)*Fa(3661)-EpsF(104)*Fa(3847))+2*(EpsF(25)*Fa
     -   (1792)+EpsF(26)*Fa(1793)+EpsF(27)*Fa(1794)+EpsF(28)*Fa(1795)
     -   +EpsF(29)*Fa(1796)+EpsF(30)*Fa(1797)+EpsF(31)*Fa(1798)+EpsF(
     -   43)*Fa(1799)+EpsF(44)*Fa(1800)+EpsF(45)*Fa(1801)+EpsF(46)*Fa
     -   (1802)+EpsF(65)*Fa(1803)+EpsF(66)*Fa(1804)+EpsF(67)*Fa(1805)
     -   +EpsF(68)*Fa(1806)+EpsF(77)*Fa(1807)+EpsF(78)*Fa(1808)+EpsF(
     -   79)*Fa(1809)+EpsF(81)*Fa(1810)+EpsF(84)*Fa(1811)+EpsF(89)*Fa
     -   (1812)+EpsF(121)*Fa(1813)+EpsF(110)*Fa(1873)+EpsF(122)*Fa(18
     -   74)+EpsF(14)*Fa(3812)+EpsF(15)*Fa(3813)+EpsF(16)*Fa(3814)+Ep
     -   sF(80)*Fa(3815)+EpsF(82)*Fa(3816)+EpsF(83)*Fa(3817)+EpsF(86)
     -   *Fa(3818)+EpsF(87)*Fa(3819)+EpsF(88)*Fa(3820)+EpsF(90)*Fa(38
     -   21)+EpsF(91)*Fa(3822)+EpsF(92)*Fa(3823)+EpsF(4)*Fa(3824)+Eps
     -   F(5)*Fa(3825)+EpsF(6)*Fa(3826)+EpsF(47)*Fa(3827)+EpsF(48)*Fa
     -   (3828)+EpsF(95)*Fa(3829)+EpsF(96)*Fa(3830)+EpsF(113)*Fa(3831
     -   )+EpsF(116)*Fa(3832)+EpsF(117)*Fa(3833)+EpsF(49)*Fa(3834)+Ep
     -   sF(50)*Fa(3835)+EpsF(51)*Fa(3836)+EpsF(61)*Fa(3837)+EpsF(62)
     -   *Fa(3838)+EpsF(63)*Fa(3839)+EpsF(94)*Fa(3840)+EpsF(101)*Fa(3
     -   841)+EpsF(102)*Fa(3842)+EpsF(105)*Fa(3843)+EpsF(112)*Fa(3844
     -   )+EpsF(124)*Fa(3845)+EpsF(103)*Fa(3846)+EpsF(119)*Fa(3848)+E
     -   psF(123)*Fa(3849)+EpsF(107)*Fa(3850)+EpsF(108)*Fa(3851)+EpsF
     -   (126)*Fa(3852)))+Cpl(5)*(M1*M4*(EpsF(76)*F(5341)-4*(EpsF(111
     -   )*Fa(210)+EpsF(125)*Fa(215))+2*(EpsF(68)*Fa(76)+EpsF(25)*Fa(
     -   1772)+EpsF(26)*Fa(1773)+EpsF(27)*Fa(1774)+EpsF(28)*Fa(1775)+
     -   EpsF(29)*Fa(1776)+EpsF(30)*Fa(1777)+EpsF(31)*Fa(1778)+EpsF(4
     -   3)*Fa(1779)+EpsF(44)*Fa(1780)+EpsF(45)*Fa(1781)+EpsF(46)*Fa(
     -   1782)+EpsF(65)*Fa(1783)-EpsF(66)*Fa(1784)+EpsF(67)*Fa(1785)+
     -   EpsF(77)*Fa(1786)+EpsF(78)*Fa(1787)+EpsF(81)*Fa(1789)+EpsF(8
     -   4)*Fa(1790)+EpsF(89)*Fa(1791)+EpsF(83)*Fa(3036)+EpsF(92)*Fa(
     -   3043)+EpsF(48)*Fa(3060)+EpsF(51)*Fa(3064)+EpsF(63)*Fa(3259)+
     -   EpsF(82)*Fa(3771)+EpsF(90)*Fa(3775)+EpsF(61)*Fa(3782)+2*(Eps
     -   F(110)*Fa(1866)+EpsF(123)*Fa(3051)-EpsF(113)*Fa(3056)-EpsF(1
     -   09)*Fa(3080)-EpsF(126)*Fa(3081)+EpsF(104)*Fa(3250)+EpsF(95)*
     -   Fa(3677)-EpsF(119)*Fa(3790)+EpsF(107)*Fa(3791))+EpsF(14)*Fa(
     -   3792)+EpsF(15)*Fa(3793)+EpsF(16)*Fa(3794)+EpsF(80)*Fa(3795)+
     -   EpsF(86)*Fa(3796)+EpsF(87)*Fa(3797)+EpsF(88)*Fa(3798)+EpsF(9
     -   1)*Fa(3799)+EpsF(4)*Fa(3800)-EpsF(5)*Fa(3801)+EpsF(6)*Fa(380
     -   2)+EpsF(116)*Fa(3803)+EpsF(50)*Fa(3804)+EpsF(62)*Fa(3805)+Ep
     -   sF(94)*Fa(3806)+EpsF(101)*Fa(3807)+EpsF(102)*Fa(3808)+EpsF(1
     -   05)*Fa(3809)+EpsF(112)*Fa(3810)+EpsF(124)*Fa(3811)))+2*(EpsF
     -   (79)*Fa(1788)+M1*M4*(-(EpsF(121)*Fa(1691))+2*(EpsF(122)*Fa(1
     -   855)+EpsF(108)*Fa(1954)+EpsF(117)*Fa(3065)+EpsF(120)*Fa(3086
     -   )-EpsF(96)*Fa(3678)-EpsF(103)*Fa(3690))+EpsF(47)*Fa(3726)+Ep
     -   sF(49)*Fa(3729))*F(7441))))))
       result(1,5) = (0,0.5)*(Cpl(3)*(Cpl(4)*(Cpl(5)*(EpsF(64)*F(5785)
     -   +4*F(6350)*(EpsF(120)*Fa(3094)+EpsF(117)*Fa(3166))+2*(EpsF(2
     -   5)*Fa(2048)+EpsF(26)*Fa(2049)+EpsF(27)*Fa(2050)+EpsF(32)*Fa(
     -   2051)+EpsF(33)*Fa(2052)+EpsF(34)*Fa(2053)+EpsF(35)*Fa(2054)+
     -   EpsF(37)*Fa(2055)+EpsF(38)*Fa(2056)+EpsF(39)*Fa(2057)+EpsF(4
     -   0)*Fa(2058)+EpsF(65)*Fa(2059)+EpsF(66)*Fa(2060)+EpsF(67)*Fa(
     -   2061)+EpsF(68)*Fa(2062)+EpsF(70)*Fa(2063)+EpsF(73)*Fa(2064)+
     -   EpsF(74)*Fa(2065)+EpsF(84)*Fa(2066)+EpsF(85)*Fa(2067)+EpsF(8
     -   9)*Fa(2068)+EpsF(121)*Fa(2069)+EpsF(110)*Fa(2256)+EpsF(111)*
     -   Fa(2257)+EpsF(122)*Fa(2258)+EpsF(125)*Fa(2259)+EpsF(17)*Fa(3
     -   853)+EpsF(18)*Fa(3854)+EpsF(19)*Fa(3855)+EpsF(71)*Fa(3856)+E
     -   psF(72)*Fa(3857)+EpsF(75)*Fa(3858)+EpsF(86)*Fa(3859)+EpsF(87
     -   )*Fa(3860)+EpsF(88)*Fa(3861)+EpsF(90)*Fa(3862)+EpsF(91)*Fa(3
     -   863)+EpsF(92)*Fa(3864)+EpsF(7)*Fa(3865)+EpsF(9)*Fa(3866)+Eps
     -   F(41)*Fa(3867)+EpsF(42)*Fa(3868)+EpsF(93)*Fa(3869)+EpsF(95)*
     -   Fa(3870)+EpsF(96)*Fa(3871)+EpsF(113)*Fa(3872)+EpsF(8)*Fa(387
     -   3)+EpsF(49)*Fa(3874)+EpsF(50)*Fa(3875)+EpsF(51)*Fa(3876)+Eps
     -   F(58)*Fa(3877)+EpsF(59)*Fa(3878)+EpsF(60)*Fa(3879)+EpsF(94)*
     -   Fa(3880)+EpsF(98)*Fa(3881)+EpsF(99)*Fa(3882)+EpsF(100)*Fa(38
     -   83)+EpsF(105)*Fa(3884)+EpsF(106)*Fa(3885)+EpsF(112)*Fa(3886)
     -   +EpsF(114)*Fa(3887)+EpsF(115)*Fa(3888)+EpsF(118)*Fa(3889)+Ep
     -   sF(107)*Fa(3890)+EpsF(108)*Fa(3891)+EpsF(109)*Fa(3892)+EpsF(
     -   126)*Fa(3893)))+M1*M4*Cpl(10)*(EpsF(64)*F(5866)+4*(EpsF(111)
     -   *Fa(277)+EpsF(125)*Fa(282)+EpsF(122)*Fa(2266)+EpsF(115)*Fa(3
     -   050)+EpsF(113)*Fa(3053)-EpsF(100)*Fa(3069)+EpsF(109)*Fa(3088
     -   )+EpsF(126)*Fa(3187))+2*(EpsF(68)*Fa(47)-EpsF(25)*Fa(2070)+E
     -   psF(26)*Fa(2071)-EpsF(27)*Fa(2072)+EpsF(32)*Fa(2073)+EpsF(33
     -   )*Fa(2074)+EpsF(34)*Fa(2075)+EpsF(35)*Fa(2076)+EpsF(37)*Fa(2
     -   077)+EpsF(38)*Fa(2078)+EpsF(39)*Fa(2079)+EpsF(40)*Fa(2080)+E
     -   psF(65)*Fa(2081)+EpsF(66)*Fa(2082)+EpsF(70)*Fa(2084)+EpsF(73
     -   )*Fa(2085)+EpsF(74)*Fa(2086)+EpsF(84)*Fa(2087)+EpsF(85)*Fa(2
     -   088)+EpsF(89)*Fa(2089)+EpsF(121)*Fa(2090)+EpsF(92)*Fa(3046)+
     -   EpsF(51)*Fa(3073)+EpsF(60)*Fa(3078)+EpsF(72)*Fa(3136)+EpsF(4
     -   2)*Fa(3160)+EpsF(17)*Fa(3894)+EpsF(18)*Fa(3895)+EpsF(19)*Fa(
     -   3896)+EpsF(71)*Fa(3897)+EpsF(75)*Fa(3898)+EpsF(86)*Fa(3899)+
     -   EpsF(87)*Fa(3900)+EpsF(88)*Fa(3901)+EpsF(90)*Fa(3902)+EpsF(9
     -   1)*Fa(3903)+EpsF(7)*Fa(3904)+EpsF(9)*Fa(3905)+EpsF(41)*Fa(39
     -   06)+EpsF(8)*Fa(3910)+EpsF(49)*Fa(3911)+EpsF(50)*Fa(3912)+Eps
     -   F(58)*Fa(3913)+EpsF(59)*Fa(3914)+EpsF(94)*Fa(3915)+EpsF(98)*
     -   Fa(3916)+EpsF(106)*Fa(3919)+EpsF(112)*Fa(3920)+EpsF(114)*Fa(
     -   3921)+F(1051)*(EpsF(67)*Fa(2083)+EpsF(93)*Fa(3907)+EpsF(105)
     -   *Fa(3918)+2*(EpsF(110)*Fa(2265)+EpsF(117)*Fa(3066)-EpsF(120)
     -   *Fa(3094)+EpsF(95)*Fa(3908)-EpsF(99)*Fa(3917)+EpsF(107)*Fa(3
     -   923)))+2*(EpsF(96)*Fa(3909)+EpsF(118)*Fa(3922)+EpsF(108)*Fa(
     -   3924)))))+M3*Cpl(9)*(M1*Cpl(10)*(EpsF(64)*F(5971)-4*(EpsF(11
     -   1)*Fa(277)+EpsF(125)*Fa(282)-EpsF(122)*Fa(2277)-EpsF(100)*Fa
     -   (3069)+EpsF(109)*Fa(3088)-EpsF(126)*Fa(3093)-EpsF(115)*Fa(31
     -   51)-EpsF(113)*Fa(3154))+2*(EpsF(68)*Fa(23)+EpsF(25)*Fa(2111)
     -   -EpsF(26)*Fa(2112)+EpsF(27)*Fa(2113)+EpsF(32)*Fa(2114)+EpsF(
     -   33)*Fa(2115)+EpsF(34)*Fa(2116)+EpsF(35)*Fa(2117)+EpsF(37)*Fa
     -   (2118)+EpsF(38)*Fa(2119)+EpsF(39)*Fa(2120)+EpsF(40)*Fa(2121)
     -   +EpsF(65)*Fa(2122)+EpsF(66)*Fa(2123)+EpsF(70)*Fa(2125)+EpsF(
     -   73)*Fa(2126)+EpsF(74)*Fa(2127)+EpsF(84)*Fa(2128)+EpsF(85)*Fa
     -   (2129)+EpsF(89)*Fa(2130)+EpsF(121)*Fa(2131)+EpsF(72)*Fa(3035
     -   )+EpsF(42)*Fa(3059)+EpsF(92)*Fa(3147)+EpsF(51)*Fa(3171)+EpsF
     -   (60)*Fa(3176)+EpsF(17)*Fa(3951)+EpsF(18)*Fa(3952)+EpsF(19)*F
     -   a(3953)+EpsF(71)*Fa(3954)+EpsF(75)*Fa(3955)+EpsF(86)*Fa(3956
     -   )+EpsF(87)*Fa(3957)+EpsF(88)*Fa(3958)+EpsF(90)*Fa(3959)-EpsF
     -   (91)*Fa(3960)+EpsF(7)*Fa(3961)+EpsF(9)*Fa(3962)+EpsF(41)*Fa(
     -   3963)-EpsF(8)*Fa(3966)+EpsF(49)*Fa(3967)+EpsF(50)*Fa(3968)+E
     -   psF(58)*Fa(3969)+EpsF(59)*Fa(3970)+EpsF(94)*Fa(3971)+EpsF(98
     -   )*Fa(3972)+EpsF(106)*Fa(3974)+EpsF(112)*Fa(3975)+EpsF(114)*F
     -   a(3976)+F(1051)*(EpsF(67)*Fa(2124)+EpsF(93)*Fa(3964)+EpsF(10
     -   5)*Fa(3973)+2*(EpsF(110)*Fa(2276)+EpsF(120)*Fa(3094)+EpsF(11
     -   7)*Fa(3166)-EpsF(95)*Fa(3908)+EpsF(99)*Fa(3917)+EpsF(107)*Fa
     -   (3978)))+2*(EpsF(96)*Fa(3965)+EpsF(118)*Fa(3977)+EpsF(108)*F
     -   a(3979))))+M4*Cpl(5)*(EpsF(64)*F(5920)+4*(EpsF(111)*Fa(286)+
     -   EpsF(125)*Fa(290))+2*(EpsF(68)*Fa(71)-EpsF(25)*Fa(2091)+EpsF
     -   (26)*Fa(2092)+EpsF(27)*Fa(2093)+EpsF(32)*Fa(2094)+EpsF(33)*F
     -   a(2095)+EpsF(34)*Fa(2096)+EpsF(35)*Fa(2097)+EpsF(37)*Fa(2098
     -   )+EpsF(38)*Fa(2099)+EpsF(39)*Fa(2100)+EpsF(40)*Fa(2101)+EpsF
     -   (65)*Fa(2102)+EpsF(66)*Fa(2103)+EpsF(70)*Fa(2104)+EpsF(73)*F
     -   a(2105)+EpsF(74)*Fa(2106)+EpsF(84)*Fa(2107)+EpsF(85)*Fa(2108
     -   )+EpsF(89)*Fa(2109)+EpsF(121)*Fa(2110)+EpsF(92)*Fa(3105)+Eps
     -   F(51)*Fa(3122)+EpsF(60)*Fa(3125)+EpsF(72)*Fa(3225)+EpsF(42)*
     -   Fa(3244)+EpsF(17)*Fa(3925)+EpsF(18)*Fa(3926)+EpsF(19)*Fa(392
     -   7)+EpsF(71)*Fa(3928)+EpsF(75)*Fa(3929)+EpsF(86)*Fa(3930)+Eps
     -   F(87)*Fa(3931)+EpsF(88)*Fa(3932)+EpsF(90)*Fa(3933)+EpsF(91)*
     -   Fa(3934)-EpsF(7)*Fa(3935)-EpsF(9)*Fa(3936)+EpsF(41)*Fa(3937)
     -   +EpsF(8)*Fa(3939)+EpsF(49)*Fa(3940)+EpsF(50)*Fa(3941)+EpsF(5
     -   8)*Fa(3942)+EpsF(59)*Fa(3943)+EpsF(94)*Fa(3944)+EpsF(98)*Fa(
     -   3945)+EpsF(106)*Fa(3946)+EpsF(112)*Fa(3947)+EpsF(114)*Fa(394
     -   8)+2*(EpsF(122)*Fa(2272)-EpsF(115)*Fa(3108)+EpsF(113)*Fa(311
     -   0)+EpsF(100)*Fa(3118)+EpsF(109)*Fa(3263)+EpsF(126)*Fa(3265)+
     -   EpsF(96)*Fa(3938)-EpsF(118)*Fa(3949)+EpsF(108)*Fa(3950))+(Ep
     -   sF(67)*Fa(2083)+EpsF(93)*Fa(3907)+EpsF(105)*Fa(3918)+2*(EpsF
     -   (110)*Fa(2265)+EpsF(117)*Fa(3066)-EpsF(120)*Fa(3094)+EpsF(95
     -   )*Fa(3908)-EpsF(99)*Fa(3917)+EpsF(107)*Fa(3923)))*F(7440))))
     -   )+Cpl(8)*(Cpl(9)*(M1*M4*Cpl(5)*(EpsF(64)*F(6066)-4*(EpsF(111
     -   )*Fa(277)+EpsF(125)*Fa(282)-EpsF(122)*Fa(2277)-EpsF(100)*Fa(
     -   3069)+EpsF(109)*Fa(3088)-EpsF(126)*Fa(3093)-EpsF(115)*Fa(315
     -   1)-EpsF(113)*Fa(3154))+2*(EpsF(68)*Fa(23)+EpsF(121)*Fa(2131)
     -   +EpsF(25)*Fa(2171)-EpsF(26)*Fa(2172)+EpsF(27)*Fa(2173)+EpsF(
     -   32)*Fa(2174)+EpsF(33)*Fa(2175)+EpsF(34)*Fa(2176)+EpsF(35)*Fa
     -   (2177)+EpsF(37)*Fa(2178)+EpsF(38)*Fa(2179)+EpsF(39)*Fa(2180)
     -   +EpsF(40)*Fa(2181)+EpsF(65)*Fa(2182)+EpsF(66)*Fa(2183)+EpsF(
     -   70)*Fa(2184)+EpsF(73)*Fa(2185)+EpsF(74)*Fa(2186)+EpsF(84)*Fa
     -   (2187)+EpsF(85)*Fa(2188)+EpsF(89)*Fa(2189)+EpsF(72)*Fa(3035)
     -   +EpsF(42)*Fa(3059)+EpsF(92)*Fa(3147)+EpsF(51)*Fa(3171)+EpsF(
     -   60)*Fa(3176)+EpsF(71)*Fa(3954)+EpsF(75)*Fa(3955)+EpsF(87)*Fa
     -   (3957)+EpsF(90)*Fa(3959)+EpsF(41)*Fa(3963)+EpsF(49)*Fa(3967)
     -   +EpsF(58)*Fa(3969)+EpsF(98)*Fa(3972)+F(1051)*(EpsF(67)*Fa(21
     -   24)+EpsF(93)*Fa(3964)+EpsF(105)*Fa(3973)+2*(EpsF(110)*Fa(227
     -   6)+EpsF(120)*Fa(3094)+EpsF(117)*Fa(3166)-EpsF(95)*Fa(3908)+E
     -   psF(99)*Fa(3917)+EpsF(107)*Fa(3978)))+2*(EpsF(96)*Fa(3965)+E
     -   psF(118)*Fa(3977)+EpsF(108)*Fa(3979))+EpsF(17)*Fa(4021)+EpsF
     -   (18)*Fa(4022)+EpsF(19)*Fa(4023)+EpsF(86)*Fa(4024)+EpsF(88)*F
     -   a(4025)-EpsF(91)*Fa(4026)+EpsF(7)*Fa(4027)+EpsF(9)*Fa(4028)-
     -   EpsF(8)*Fa(4029)+EpsF(50)*Fa(4030)+EpsF(59)*Fa(4031)+EpsF(94
     -   )*Fa(4032)+EpsF(106)*Fa(4033)+EpsF(112)*Fa(4034)+EpsF(114)*F
     -   a(4035)))+Cpl(10)*(EpsF(64)*F(6123)+4*F(6350)*(EpsF(117)*Fa(
     -   3066)-EpsF(120)*Fa(3094))+2*(EpsF(25)*Fa(2190)+EpsF(26)*Fa(2
     -   191)+EpsF(27)*Fa(2192)+EpsF(32)*Fa(2193)+EpsF(33)*Fa(2194)+E
     -   psF(34)*Fa(2195)+EpsF(35)*Fa(2196)+EpsF(37)*Fa(2197)+EpsF(38
     -   )*Fa(2198)+EpsF(39)*Fa(2199)+EpsF(40)*Fa(2200)+EpsF(65)*Fa(2
     -   201)+EpsF(66)*Fa(2202)+EpsF(67)*Fa(2203)+EpsF(68)*Fa(2204)+E
     -   psF(70)*Fa(2205)+EpsF(73)*Fa(2206)+EpsF(74)*Fa(2207)+EpsF(84
     -   )*Fa(2208)+EpsF(85)*Fa(2209)+EpsF(89)*Fa(2210)+EpsF(121)*Fa(
     -   2211)+EpsF(110)*Fa(2287)+EpsF(111)*Fa(2288)+EpsF(122)*Fa(228
     -   9)+EpsF(125)*Fa(2290)+EpsF(17)*Fa(4036)+EpsF(18)*Fa(4037)+Ep
     -   sF(19)*Fa(4038)+EpsF(71)*Fa(4039)+EpsF(72)*Fa(4040)+EpsF(75)
     -   *Fa(4041)+EpsF(86)*Fa(4042)+EpsF(87)*Fa(4043)+EpsF(88)*Fa(40
     -   44)+EpsF(90)*Fa(4045)+EpsF(91)*Fa(4046)+EpsF(92)*Fa(4047)+Ep
     -   sF(7)*Fa(4048)+EpsF(9)*Fa(4049)+EpsF(41)*Fa(4050)+EpsF(42)*F
     -   a(4051)+EpsF(93)*Fa(4052)+EpsF(95)*Fa(4053)+EpsF(96)*Fa(4054
     -   )+EpsF(113)*Fa(4055)+EpsF(8)*Fa(4056)+EpsF(49)*Fa(4057)+EpsF
     -   (50)*Fa(4058)+EpsF(51)*Fa(4059)+EpsF(58)*Fa(4060)+EpsF(59)*F
     -   a(4061)+EpsF(60)*Fa(4062)+EpsF(94)*Fa(4063)+EpsF(98)*Fa(4064
     -   )+EpsF(99)*Fa(4065)+EpsF(100)*Fa(4066)+EpsF(105)*Fa(4067)+Ep
     -   sF(106)*Fa(4068)+EpsF(112)*Fa(4069)+EpsF(114)*Fa(4070)+EpsF(
     -   115)*Fa(4071)+EpsF(118)*Fa(4072)+EpsF(107)*Fa(4073)+EpsF(108
     -   )*Fa(4074)+EpsF(109)*Fa(4075)+EpsF(126)*Fa(4076))))+M3*Cpl(4
     -   )*(M1*Cpl(5)*(EpsF(64)*F(6009)+4*(EpsF(111)*Fa(277)+EpsF(125
     -   )*Fa(282)+EpsF(122)*Fa(2266)+EpsF(115)*Fa(3050)+EpsF(113)*Fa
     -   (3053)-EpsF(100)*Fa(3069)+EpsF(109)*Fa(3088)+EpsF(126)*Fa(31
     -   87))+2*(EpsF(68)*Fa(47)+EpsF(121)*Fa(2090)-EpsF(25)*Fa(2132)
     -   +EpsF(26)*Fa(2133)-EpsF(27)*Fa(2134)+EpsF(32)*Fa(2135)+EpsF(
     -   33)*Fa(2136)+EpsF(34)*Fa(2137)+EpsF(35)*Fa(2138)+EpsF(37)*Fa
     -   (2139)+EpsF(38)*Fa(2140)+EpsF(39)*Fa(2141)+EpsF(40)*Fa(2142)
     -   +EpsF(65)*Fa(2143)+EpsF(66)*Fa(2144)+EpsF(70)*Fa(2145)+EpsF(
     -   73)*Fa(2146)+EpsF(74)*Fa(2147)+EpsF(84)*Fa(2148)+EpsF(85)*Fa
     -   (2149)+EpsF(89)*Fa(2150)+EpsF(92)*Fa(3046)+EpsF(51)*Fa(3073)
     -   +EpsF(60)*Fa(3078)+EpsF(72)*Fa(3136)+EpsF(42)*Fa(3160)+EpsF(
     -   71)*Fa(3897)+EpsF(75)*Fa(3898)+EpsF(87)*Fa(3900)+EpsF(90)*Fa
     -   (3902)+EpsF(41)*Fa(3906)+EpsF(49)*Fa(3911)+EpsF(58)*Fa(3913)
     -   +EpsF(98)*Fa(3916)+F(1051)*(EpsF(67)*Fa(2083)+EpsF(93)*Fa(39
     -   07)+EpsF(105)*Fa(3918)+2*(EpsF(110)*Fa(2265)+EpsF(117)*Fa(30
     -   66)-EpsF(120)*Fa(3094)+EpsF(95)*Fa(3908)-EpsF(99)*Fa(3917)+E
     -   psF(107)*Fa(3923)))+2*(EpsF(96)*Fa(3909)+EpsF(118)*Fa(3922)+
     -   EpsF(108)*Fa(3924))+EpsF(17)*Fa(3980)+EpsF(18)*Fa(3981)+EpsF
     -   (19)*Fa(3982)+EpsF(86)*Fa(3983)+EpsF(88)*Fa(3984)+EpsF(91)*F
     -   a(3985)+EpsF(7)*Fa(3986)+EpsF(9)*Fa(3987)+EpsF(8)*Fa(3988)+E
     -   psF(50)*Fa(3989)+EpsF(59)*Fa(3990)+EpsF(94)*Fa(3991)+EpsF(10
     -   6)*Fa(3992)+EpsF(112)*Fa(3993)+EpsF(114)*Fa(3994)))+M4*Cpl(1
     -   0)*(EpsF(64)*F(6038)-4*(EpsF(111)*Fa(286)+EpsF(125)*Fa(290))
     -   +2*(EpsF(68)*Fa(36)+EpsF(25)*Fa(2151)-EpsF(26)*Fa(2152)-EpsF
     -   (27)*Fa(2153)+EpsF(32)*Fa(2154)+EpsF(33)*Fa(2155)+EpsF(34)*F
     -   a(2156)+EpsF(35)*Fa(2157)+EpsF(37)*Fa(2158)+EpsF(38)*Fa(2159
     -   )+EpsF(39)*Fa(2160)+EpsF(40)*Fa(2161)+EpsF(65)*Fa(2162)+EpsF
     -   (66)*Fa(2163)+EpsF(70)*Fa(2164)+EpsF(73)*Fa(2165)+EpsF(74)*F
     -   a(2166)+EpsF(84)*Fa(2167)+EpsF(85)*Fa(2168)+EpsF(89)*Fa(2169
     -   )+EpsF(121)*Fa(2170)+EpsF(72)*Fa(3097)+EpsF(42)*Fa(3114)+Eps
     -   F(92)*Fa(3236)+EpsF(51)*Fa(3254)+EpsF(60)*Fa(3258)+EpsF(17)*
     -   Fa(3995)+EpsF(18)*Fa(3996)+EpsF(19)*Fa(3997)+EpsF(71)*Fa(399
     -   8)+EpsF(75)*Fa(3999)+EpsF(86)*Fa(4000)+EpsF(87)*Fa(4001)+Eps
     -   F(88)*Fa(4002)+EpsF(90)*Fa(4003)+EpsF(91)*Fa(4004)+EpsF(7)*F
     -   a(4005)+EpsF(9)*Fa(4006)+EpsF(41)*Fa(4007)+EpsF(8)*Fa(4009)+
     -   EpsF(49)*Fa(4010)+EpsF(50)*Fa(4011)+EpsF(58)*Fa(4012)+EpsF(5
     -   9)*Fa(4013)+EpsF(94)*Fa(4014)+EpsF(98)*Fa(4015)+EpsF(106)*Fa
     -   (4016)+EpsF(112)*Fa(4017)+EpsF(114)*Fa(4018)+2*(EpsF(122)*Fa
     -   (2280)+EpsF(115)*Fa(3108)+EpsF(109)*Fa(3129)+EpsF(126)*Fa(31
     -   33)+EpsF(113)*Fa(3240)+EpsF(100)*Fa(3251)+EpsF(96)*Fa(4008)+
     -   EpsF(118)*Fa(4019)+EpsF(108)*Fa(4020))+(EpsF(67)*Fa(2124)+Ep
     -   sF(93)*Fa(3964)+EpsF(105)*Fa(3973)+2*(EpsF(110)*Fa(2276)+Eps
     -   F(120)*Fa(3094)+EpsF(117)*Fa(3166)-EpsF(95)*Fa(3908)+EpsF(99
     -   )*Fa(3917)+EpsF(107)*Fa(3978)))*F(7440))))))
       result(1,6) = (0,0.5)*(Cpl(3)*(Cpl(4)*(Cpl(5)*(EpsF(20)*F(6404)
     -   -4*(-(F(7398)*(EpsF(117)*Fa(498)+EpsF(120)*Fa(735)))+EpsF(11
     -   1)*Fa(2690)+EpsF(125)*Fa(2692)-EpsF(100)*Fa(4117)+EpsF(104)*
     -   Fa(4121)-EpsF(119)*Fa(4122))+2*(EpsF(10)*Fa(2432)+EpsF(11)*F
     -   a(2433)+EpsF(12)*Fa(2434)+EpsF(13)*Fa(2435)+EpsF(21)*Fa(2436
     -   )+EpsF(22)*Fa(2437)+EpsF(23)*Fa(2438)+EpsF(24)*Fa(2439)+EpsF
     -   (28)*Fa(2440)+EpsF(29)*Fa(2441)+EpsF(30)*Fa(2442)+EpsF(31)*F
     -   a(2443)+EpsF(32)*Fa(2444)+EpsF(33)*Fa(2445)+EpsF(34)*Fa(2446
     -   )+EpsF(35)*Fa(2447)+EpsF(110)*Fa(2689)+EpsF(122)*Fa(2691)+Ep
     -   sF(25)*Fa(4077)+EpsF(26)*Fa(4078)+EpsF(27)*Fa(4079)+EpsF(67)
     -   *Fa(4080)+EpsF(68)*Fa(4081)+EpsF(121)*Fa(4082)+EpsF(17)*Fa(4
     -   083)+EpsF(18)*Fa(4084)+EpsF(19)*Fa(4085)+EpsF(71)*Fa(4086)+E
     -   psF(72)*Fa(4087)+EpsF(75)*Fa(4088)+EpsF(14)*Fa(4089)+EpsF(15
     -   )*Fa(4090)+EpsF(16)*Fa(4091)+EpsF(80)*Fa(4092)+EpsF(82)*Fa(4
     -   093)+EpsF(83)*Fa(4094)+EpsF(1)*Fa(4095)+EpsF(2)*Fa(4096)+Eps
     -   F(3)*Fa(4097)+EpsF(4)*Fa(4098)+EpsF(5)*Fa(4099)+EpsF(6)*Fa(4
     -   100)+EpsF(7)*Fa(4101)+EpsF(9)*Fa(4102)+EpsF(41)*Fa(4103)+Eps
     -   F(42)*Fa(4104)+EpsF(47)*Fa(4105)+EpsF(48)*Fa(4106)+EpsF(56)*
     -   Fa(4107)+EpsF(57)*Fa(4108)+EpsF(93)*Fa(4109)+EpsF(95)*Fa(411
     -   0)+EpsF(96)*Fa(4111)+EpsF(97)*Fa(4112)+EpsF(113)*Fa(4113)+Ep
     -   sF(116)*Fa(4114)+EpsF(8)*Fa(4115)+EpsF(99)*Fa(4116)+EpsF(115
     -   )*Fa(4118)+EpsF(118)*Fa(4119)+EpsF(103)*Fa(4120)+EpsF(123)*F
     -   a(4123)))+M1*M4*Cpl(10)*(EpsF(20)*F(6480)+2*(EpsF(10)*Fa(244
     -   8)+EpsF(11)*Fa(2449)+EpsF(12)*Fa(2450)+EpsF(13)*Fa(2451)+Eps
     -   F(21)*Fa(2452)+EpsF(22)*Fa(2453)+EpsF(23)*Fa(2454)+EpsF(24)*
     -   Fa(2455)+EpsF(28)*Fa(2456)+EpsF(29)*Fa(2457)+EpsF(30)*Fa(245
     -   8)+EpsF(31)*Fa(2459)+EpsF(32)*Fa(2460)+EpsF(33)*Fa(2461)+Eps
     -   F(34)*Fa(2462)+EpsF(35)*Fa(2463)+EpsF(25)*Fa(4124)+EpsF(26)*
     -   Fa(4125)+EpsF(27)*Fa(4126)+EpsF(67)*Fa(4127)+EpsF(68)*Fa(412
     -   8)+EpsF(121)*Fa(4129)+EpsF(17)*Fa(4130)+EpsF(18)*Fa(4131)+Ep
     -   sF(19)*Fa(4132)+EpsF(71)*Fa(4133)+EpsF(72)*Fa(4134)+EpsF(75)
     -   *Fa(4135)+EpsF(14)*Fa(4136)+EpsF(15)*Fa(4137)+EpsF(16)*Fa(41
     -   38)+EpsF(80)*Fa(4139)+EpsF(82)*Fa(4140)+EpsF(83)*Fa(4141)+Ep
     -   sF(1)*Fa(4142)+EpsF(2)*Fa(4143)+EpsF(3)*Fa(4144)+EpsF(4)*Fa(
     -   4145)+EpsF(5)*Fa(4146)+EpsF(6)*Fa(4147)+EpsF(7)*Fa(4148)+Eps
     -   F(9)*Fa(4149)+EpsF(41)*Fa(4150)+EpsF(42)*Fa(4151)+EpsF(47)*F
     -   a(4152)+EpsF(48)*Fa(4153)+EpsF(56)*Fa(4154)+EpsF(57)*Fa(4155
     -   )+EpsF(93)*Fa(4156)+EpsF(97)*Fa(4159)+EpsF(116)*Fa(4161)+Eps
     -   F(8)*Fa(4162))+4*(F(7418)*(EpsF(117)*Fa(541)-EpsF(120)*Fa(73
     -   5))+EpsF(110)*Fa(2698)+EpsF(111)*Fa(2699)+EpsF(122)*Fa(2700)
     -   +EpsF(125)*Fa(2701)+EpsF(95)*Fa(4157)+EpsF(96)*Fa(4158)+EpsF
     -   (113)*Fa(4160)+EpsF(99)*Fa(4163)-EpsF(100)*Fa(4164)+EpsF(115
     -   )*Fa(4165)+EpsF(118)*Fa(4166)+EpsF(103)*Fa(4167)+EpsF(104)*F
     -   a(4168)+EpsF(119)*Fa(4169)+EpsF(123)*Fa(4170))))+Cpl(9)*(M3*
     -   M4*Cpl(5)*(EpsF(20)*F(6553)+4*(EpsF(110)*Fa(2707)+EpsF(111)*
     -   Fa(2708)+EpsF(122)*Fa(2709)+EpsF(125)*Fa(2710)+EpsF(117)*Fa(
     -   2796)+EpsF(120)*Fa(2909))+2*(EpsF(10)*Fa(2464)+EpsF(11)*Fa(2
     -   465)+EpsF(12)*Fa(2466)+EpsF(13)*Fa(2467)+EpsF(21)*Fa(2468)+E
     -   psF(22)*Fa(2469)+EpsF(23)*Fa(2470)+EpsF(24)*Fa(2471)+EpsF(28
     -   )*Fa(2472)+EpsF(29)*Fa(2473)+EpsF(30)*Fa(2474)+EpsF(31)*Fa(2
     -   475)+EpsF(32)*Fa(2476)+EpsF(33)*Fa(2477)+EpsF(34)*Fa(2478)+E
     -   psF(35)*Fa(2479)+EpsF(25)*Fa(4171)+EpsF(26)*Fa(4172)+EpsF(27
     -   )*Fa(4173)+EpsF(67)*Fa(4174)+EpsF(68)*Fa(4175)+EpsF(121)*Fa(
     -   4176)+EpsF(17)*Fa(4177)+EpsF(18)*Fa(4178)+EpsF(19)*Fa(4179)+
     -   EpsF(71)*Fa(4180)+EpsF(72)*Fa(4181)+EpsF(75)*Fa(4182)+EpsF(1
     -   4)*Fa(4183)+EpsF(15)*Fa(4184)+EpsF(16)*Fa(4185)+EpsF(80)*Fa(
     -   4186)+EpsF(82)*Fa(4187)+EpsF(83)*Fa(4188)+EpsF(1)*Fa(4189)+E
     -   psF(2)*Fa(4190)+EpsF(3)*Fa(4191)+EpsF(4)*Fa(4192)+EpsF(5)*Fa
     -   (4193)+EpsF(6)*Fa(4194)+EpsF(7)*Fa(4195)+EpsF(9)*Fa(4196)+Ep
     -   sF(41)*Fa(4197)+EpsF(42)*Fa(4198)+EpsF(47)*Fa(4199)+EpsF(48)
     -   *Fa(4200)+EpsF(56)*Fa(4201)+EpsF(57)*Fa(4202)+EpsF(93)*Fa(42
     -   03)+EpsF(97)*Fa(4206)+EpsF(116)*Fa(4208)+EpsF(8)*Fa(4209)+2*
     -   (EpsF(95)*Fa(4204)+EpsF(96)*Fa(4205)+EpsF(113)*Fa(4207)+EpsF
     -   (99)*Fa(4210)+EpsF(100)*Fa(4211)+EpsF(115)*Fa(4212)+EpsF(118
     -   )*Fa(4213)+EpsF(103)*Fa(4214)+EpsF(104)*Fa(4215)+EpsF(119)*F
     -   a(4216)+EpsF(123)*Fa(4217))))+Cpl(10)*(2*EpsF(33)*Fa(2493)+M
     -   1*M3*(EpsF(20)*F(6623)+4*(EpsF(110)*Fa(2716)+EpsF(111)*Fa(27
     -   17)+EpsF(122)*Fa(2718)+EpsF(125)*Fa(2719)+EpsF(117)*Fa(2819)
     -   +EpsF(120)*Fa(2914))+2*(EpsF(10)*Fa(2480)+EpsF(11)*Fa(2481)+
     -   EpsF(12)*Fa(2482)+EpsF(13)*Fa(2483)+EpsF(21)*Fa(2484)+EpsF(2
     -   2)*Fa(2485)+EpsF(23)*Fa(2486)+EpsF(24)*Fa(2487)+EpsF(28)*Fa(
     -   2488)+EpsF(29)*Fa(2489)+EpsF(30)*Fa(2490)+EpsF(31)*Fa(2491)+
     -   EpsF(32)*Fa(2492)+EpsF(34)*Fa(2494)+EpsF(35)*Fa(2495)+EpsF(2
     -   5)*Fa(4218)+EpsF(26)*Fa(4219)+EpsF(27)*Fa(4220)+EpsF(67)*Fa(
     -   4221)+EpsF(68)*Fa(4222)+EpsF(121)*Fa(4223)+EpsF(17)*Fa(4224)
     -   +EpsF(18)*Fa(4225)+EpsF(19)*Fa(4226)+EpsF(71)*Fa(4227)+EpsF(
     -   72)*Fa(4228)+EpsF(75)*Fa(4229)+EpsF(14)*Fa(4230)+EpsF(15)*Fa
     -   (4231)+EpsF(16)*Fa(4232)+EpsF(80)*Fa(4233)+EpsF(82)*Fa(4234)
     -   +EpsF(83)*Fa(4235)+EpsF(1)*Fa(4236)+EpsF(2)*Fa(4237)+EpsF(3)
     -   *Fa(4238)+EpsF(4)*Fa(4239)+EpsF(5)*Fa(4240)+EpsF(6)*Fa(4241)
     -   +EpsF(7)*Fa(4242)+EpsF(9)*Fa(4243)+EpsF(41)*Fa(4244)+EpsF(42
     -   )*Fa(4245)+EpsF(47)*Fa(4246)+EpsF(48)*Fa(4247)+EpsF(56)*Fa(4
     -   248)+EpsF(57)*Fa(4249)+EpsF(93)*Fa(4250)+EpsF(97)*Fa(4253)+E
     -   psF(116)*Fa(4255)+EpsF(8)*Fa(4256)+2*(EpsF(95)*Fa(4251)+EpsF
     -   (96)*Fa(4252)+EpsF(113)*Fa(4254)+EpsF(99)*Fa(4257)+EpsF(100)
     -   *Fa(4258)+EpsF(115)*Fa(4259)+EpsF(118)*Fa(4260)+EpsF(103)*Fa
     -   (4261)+EpsF(104)*Fa(4262)+EpsF(119)*Fa(4263)+EpsF(123)*Fa(42
     -   64)))))))+Cpl(8)*(Cpl(4)*(Cpl(5)*(2*(EpsF(23)*Fa(2502)+EpsF(
     -   24)*Fa(2503)+EpsF(29)*Fa(2505))+M1*M3*(EpsF(20)*F(6657)+4*(E
     -   psF(110)*Fa(2707)+EpsF(111)*Fa(2708)+EpsF(122)*Fa(2709)+EpsF
     -   (125)*Fa(2710)+EpsF(117)*Fa(2796)+EpsF(120)*Fa(2909))+2*(Eps
     -   F(10)*Fa(2496)+EpsF(11)*Fa(2497)+EpsF(12)*Fa(2498)+EpsF(13)*
     -   Fa(2499)+EpsF(21)*Fa(2500)+EpsF(22)*Fa(2501)+EpsF(28)*Fa(250
     -   4)+EpsF(30)*Fa(2506)+EpsF(31)*Fa(2507)+EpsF(32)*Fa(2508)+Eps
     -   F(33)*Fa(2509)+EpsF(34)*Fa(2510)+EpsF(35)*Fa(2511)+EpsF(67)*
     -   Fa(4174)+EpsF(68)*Fa(4175)+EpsF(121)*Fa(4176)+EpsF(71)*Fa(41
     -   80)+EpsF(72)*Fa(4181)+EpsF(75)*Fa(4182)+EpsF(80)*Fa(4186)+Ep
     -   sF(82)*Fa(4187)+EpsF(83)*Fa(4188)+EpsF(41)*Fa(4197)+EpsF(42)
     -   *Fa(4198)+EpsF(47)*Fa(4199)+EpsF(48)*Fa(4200)+EpsF(56)*Fa(42
     -   01)+EpsF(57)*Fa(4202)+EpsF(93)*Fa(4203)+EpsF(97)*Fa(4206)+Ep
     -   sF(116)*Fa(4208)+2*(EpsF(95)*Fa(4204)+EpsF(96)*Fa(4205)+EpsF
     -   (113)*Fa(4207)+EpsF(99)*Fa(4210)+EpsF(100)*Fa(4211)+EpsF(115
     -   )*Fa(4212)+EpsF(118)*Fa(4213)+EpsF(103)*Fa(4214)+EpsF(104)*F
     -   a(4215)+EpsF(119)*Fa(4216)+EpsF(123)*Fa(4217))+EpsF(25)*Fa(4
     -   265)+EpsF(26)*Fa(4266)+EpsF(27)*Fa(4267)+EpsF(17)*Fa(4268)+E
     -   psF(18)*Fa(4269)+EpsF(19)*Fa(4270)+EpsF(14)*Fa(4271)+EpsF(15
     -   )*Fa(4272)+EpsF(16)*Fa(4273)+EpsF(1)*Fa(4274)+EpsF(2)*Fa(427
     -   5)+EpsF(3)*Fa(4276)+EpsF(4)*Fa(4277)+EpsF(5)*Fa(4278)+EpsF(6
     -   )*Fa(4279)+EpsF(7)*Fa(4280)+EpsF(9)*Fa(4281)+EpsF(8)*Fa(4282
     -   ))))+M3*M4*Cpl(10)*(EpsF(20)*F(6674)+4*(EpsF(110)*Fa(2716)+E
     -   psF(111)*Fa(2717)+EpsF(122)*Fa(2718)+EpsF(125)*Fa(2719)+EpsF
     -   (117)*Fa(2819)+EpsF(120)*Fa(2914))+2*(EpsF(10)*Fa(2512)+EpsF
     -   (11)*Fa(2513)+EpsF(12)*Fa(2514)+EpsF(13)*Fa(2515)+EpsF(21)*F
     -   a(2516)+EpsF(22)*Fa(2517)+EpsF(23)*Fa(2518)+EpsF(24)*Fa(2519
     -   )+EpsF(28)*Fa(2520)+EpsF(29)*Fa(2521)+EpsF(30)*Fa(2522)+EpsF
     -   (31)*Fa(2523)+EpsF(32)*Fa(2524)+EpsF(33)*Fa(2525)+EpsF(34)*F
     -   a(2526)+EpsF(35)*Fa(2527)+EpsF(67)*Fa(4221)+EpsF(68)*Fa(4222
     -   )+EpsF(121)*Fa(4223)+EpsF(71)*Fa(4227)+EpsF(72)*Fa(4228)+Eps
     -   F(75)*Fa(4229)+EpsF(80)*Fa(4233)+EpsF(82)*Fa(4234)+EpsF(83)*
     -   Fa(4235)+EpsF(41)*Fa(4244)+EpsF(42)*Fa(4245)+EpsF(47)*Fa(424
     -   6)+EpsF(48)*Fa(4247)+EpsF(56)*Fa(4248)+EpsF(57)*Fa(4249)+Eps
     -   F(93)*Fa(4250)+EpsF(97)*Fa(4253)+EpsF(116)*Fa(4255)+2*(EpsF(
     -   95)*Fa(4251)+EpsF(96)*Fa(4252)+EpsF(113)*Fa(4254)+EpsF(99)*F
     -   a(4257)+EpsF(100)*Fa(4258)+EpsF(115)*Fa(4259)+EpsF(118)*Fa(4
     -   260)+EpsF(103)*Fa(4261)+EpsF(104)*Fa(4262)+EpsF(119)*Fa(4263
     -   )+EpsF(123)*Fa(4264))+EpsF(25)*Fa(4283)+EpsF(26)*Fa(4284)+Ep
     -   sF(27)*Fa(4285)+EpsF(17)*Fa(4286)+EpsF(18)*Fa(4287)+EpsF(19)
     -   *Fa(4288)+EpsF(14)*Fa(4289)+EpsF(15)*Fa(4290)+EpsF(16)*Fa(42
     -   91)+EpsF(1)*Fa(4292)+EpsF(2)*Fa(4293)+EpsF(3)*Fa(4294)+EpsF(
     -   4)*Fa(4295)+EpsF(5)*Fa(4296)+EpsF(6)*Fa(4297)+EpsF(7)*Fa(429
     -   8)+EpsF(9)*Fa(4299)+EpsF(8)*Fa(4300))))+Cpl(9)*(M1*M4*Cpl(5)
     -   *(EpsF(20)*F(6697)+2*(EpsF(10)*Fa(2528)+EpsF(11)*Fa(2529)+Ep
     -   sF(12)*Fa(2530)+EpsF(13)*Fa(2531)+EpsF(21)*Fa(2532)+EpsF(22)
     -   *Fa(2533)+EpsF(23)*Fa(2534)-EpsF(24)*Fa(2535)+EpsF(28)*Fa(25
     -   36)+EpsF(29)*Fa(2537)+EpsF(30)*Fa(2538)+EpsF(31)*Fa(2539)+Ep
     -   sF(32)*Fa(2540)+EpsF(33)*Fa(2541)+EpsF(34)*Fa(2542)+EpsF(35)
     -   *Fa(2543)+EpsF(25)*Fa(4301)+EpsF(26)*Fa(4302)+EpsF(27)*Fa(43
     -   03)+EpsF(67)*Fa(4304)+EpsF(68)*Fa(4305)+EpsF(121)*Fa(4306)+E
     -   psF(17)*Fa(4307)+EpsF(18)*Fa(4308)+EpsF(19)*Fa(4309)+EpsF(71
     -   )*Fa(4310)+EpsF(72)*Fa(4311)+EpsF(75)*Fa(4312)+EpsF(14)*Fa(4
     -   313)+EpsF(15)*Fa(4314)-EpsF(16)*Fa(4315)+EpsF(80)*Fa(4316)+E
     -   psF(82)*Fa(4317)+EpsF(83)*Fa(4318)+EpsF(1)*Fa(4319)+EpsF(2)*
     -   Fa(4320)+EpsF(3)*Fa(4321)+EpsF(4)*Fa(4322)+EpsF(5)*Fa(4323)+
     -   EpsF(6)*Fa(4324)+EpsF(7)*Fa(4325)+EpsF(9)*Fa(4326)+EpsF(41)*
     -   Fa(4327)+EpsF(42)*Fa(4328)+EpsF(47)*Fa(4329)+EpsF(48)*Fa(433
     -   0)+EpsF(56)*Fa(4331)+EpsF(57)*Fa(4332)+EpsF(93)*Fa(4333)+Eps
     -   F(97)*Fa(4336)+EpsF(116)*Fa(4338)+EpsF(8)*Fa(4339))-4*(-(F(7
     -   418)*(EpsF(117)*Fa(498)+EpsF(120)*Fa(735)))+EpsF(111)*Fa(269
     -   9)+EpsF(125)*Fa(2701)-EpsF(110)*Fa(2721)-EpsF(122)*Fa(2722)-
     -   EpsF(100)*Fa(4164)+EpsF(104)*Fa(4168)-EpsF(95)*Fa(4334)-EpsF
     -   (96)*Fa(4335)-EpsF(113)*Fa(4337)-EpsF(99)*Fa(4340)-EpsF(115)
     -   *Fa(4341)-EpsF(118)*Fa(4342)-EpsF(103)*Fa(4343)-EpsF(119)*Fa
     -   (4344)-EpsF(123)*Fa(4345)))+Cpl(10)*(EpsF(20)*F(6756)+4*(F(7
     -   398)*(EpsF(117)*Fa(541)-EpsF(120)*Fa(735))+EpsF(111)*Fa(2690
     -   )+EpsF(125)*Fa(2692)-EpsF(100)*Fa(4117)+EpsF(104)*Fa(4121)+E
     -   psF(113)*Fa(4382)-EpsF(119)*Fa(4389))+2*(EpsF(10)*Fa(2544)+E
     -   psF(11)*Fa(2545)+EpsF(12)*Fa(2546)+EpsF(13)*Fa(2547)+EpsF(21
     -   )*Fa(2548)+EpsF(22)*Fa(2549)+EpsF(23)*Fa(2550)+EpsF(24)*Fa(2
     -   551)+EpsF(28)*Fa(2552)+EpsF(29)*Fa(2553)+EpsF(30)*Fa(2554)+E
     -   psF(31)*Fa(2555)+EpsF(32)*Fa(2556)+EpsF(33)*Fa(2557)+EpsF(34
     -   )*Fa(2558)+EpsF(35)*Fa(2559)+EpsF(110)*Fa(2728)+EpsF(122)*Fa
     -   (2729)+EpsF(25)*Fa(4346)+EpsF(26)*Fa(4347)+EpsF(27)*Fa(4348)
     -   +EpsF(67)*Fa(4349)+EpsF(68)*Fa(4350)+EpsF(121)*Fa(4351)+EpsF
     -   (17)*Fa(4352)+EpsF(18)*Fa(4353)+EpsF(19)*Fa(4354)+EpsF(71)*F
     -   a(4355)+EpsF(72)*Fa(4356)+EpsF(75)*Fa(4357)+EpsF(14)*Fa(4358
     -   )+EpsF(15)*Fa(4359)+EpsF(16)*Fa(4360)+EpsF(80)*Fa(4361)+EpsF
     -   (82)*Fa(4362)+EpsF(83)*Fa(4363)+EpsF(1)*Fa(4364)+EpsF(2)*Fa(
     -   4365)+EpsF(3)*Fa(4366)+EpsF(4)*Fa(4367)+EpsF(5)*Fa(4368)+Eps
     -   F(6)*Fa(4369)+EpsF(7)*Fa(4370)+EpsF(9)*Fa(4371)+EpsF(41)*Fa(
     -   4372)+EpsF(42)*Fa(4373)+EpsF(47)*Fa(4374)+EpsF(48)*Fa(4375)+
     -   EpsF(56)*Fa(4376)+EpsF(57)*Fa(4377)+EpsF(93)*Fa(4378)+EpsF(9
     -   5)*Fa(4379)+EpsF(96)*Fa(4380)+EpsF(97)*Fa(4381)+EpsF(116)*Fa
     -   (4383)+EpsF(8)*Fa(4384)+EpsF(99)*Fa(4385)+EpsF(115)*Fa(4386)
     -   +EpsF(118)*Fa(4387)+EpsF(103)*Fa(4388)+EpsF(123)*Fa(4390))))
     -   ))
       result(1,1) =-result(1,1)
       result(1,2) =-result(1,2)
       result(1,3) =-result(1,3)
       result(1,4) =-result(1,4)
       result(1,5) =-result(1,5)
       result(1,6) =-result(1,6)
c************************************************************************************
c************************************************************************************
       Return
       End
