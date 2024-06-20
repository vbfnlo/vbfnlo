       subroutine PenEvenG_Div(cpl,MM1,MM3,MM4,p1,p2,p3,p4,p5, mup1, mup2,
     -   mup3,mup4,mup5,musq,comp,result,Div)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/02/2012
c Modified:8/2/2012
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
       real*8 cpl(10) 
      Complex*16 F(4300),Fa(150)
      Real*8 P(2563)
      Complex*16 dotrc,dotcc,result(1,6)
       Real*8 musq
      EXTERNAL   dotrc,dotcc

       COMMON/FaPenEvenG/Fa
       COMMON/FPenEvenG/F
       COMMON/PPenEvenG/P
      Save/FPenEvenG/
      Save/PPenEvenG/
       COMMON/InvariantsPenEvenG/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,
     -   s45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       COMMON/EVALFunctionsPenEvenG/ E012345R,E12345R,E012345I,E12345I
       COMMON/DVALFunctionsPenEvenG/ D01234R,D1234R,D01234I,D1234I, D0
     -   1235R,D1235R,D01235I,D1235I, D01245R,D1245R,D01245I,D1245I, 
     -   D01345R,D1345R,D01345I,D1345I, D02345R,D2345R,D02345I,D2345I
       COMMON/CVALFunctionsPenEvenG/ C0123R,C123R,C0123I,C123I, C0124R
     -   ,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0134R,C134R
     -   ,C0134I,C134I, C0135R,C135R,C0135I,C135I, C0145R,C145R,C0145
     -   I,C145I, C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235
     -   I, C0245R,C245R,C0245I,C245I, C0345R,C345R,C0345I,C345I 
       COMMON/BVALFunctionsPenEvenG/ B012R,B012I, B013R,B013I, B014R,B
     -   014I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,B025I, B0
     -   34R,B034I, B035R,B035I, B045R,B045I, B12R,B12I, B13R,B13I, B
     -   14R,B14I, B15R,B15I, B23R,B23I, B24R,B24I, B25R,B25I, B34R,B
     -   34I, B35R,B35I, B45R,B45I 

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
      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False.
c************************************************************************************
c************************************************************************************
        if (comp.gt.0) then    
c    Calling C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************
       call  XEget_R5_Div(M1,M1,M1,M3,M4,p1sq,p2sq,p3sq,p4sq,p5sq,p1
     -   p2,p1p3,p1p4,p2p3,p2p4,p3p4,p1p5,p4p5,musq,Div)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call FPenEvenG1(F(1))
       call FPenEvenG2(F(216))
       call FPenEvenG3(F(431))
       call FPenEvenG4(F(646))
       call FPenEvenG5(F(861))
       call FPenEvenG6(F(1076))
       call FPenEvenG7(F(1291))
       call FPenEvenG8(F(1506))
       call FPenEvenG9(F(1721))
       call FPenEvenG10(F(1936))
       call FPenEvenG11(F(2151))
       call FPenEvenG12(F(2366))
       call FPenEvenG13(F(2581))
       call FPenEvenG14(F(2796))
       call FPenEvenG15(F(3011))
       call FPenEvenG16(F(3226))
       call FPenEvenG17(F(3441))
       call FPenEvenG18(F(3656))
       call FPenEvenG19(F(3871))
       call FPenEvenG20(F(4086))
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
c************************************************************************************
c************************************************************************************
       call FaPenEvenG1(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1))
       call FaPenEvenG2(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(11))
       call FaPenEvenG3(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(21))
       call FaPenEvenG4(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(31))
       call FaPenEvenG5(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(41))
       call FaPenEvenG6(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(51))
       call FaPenEvenG7(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(61))
       call FaPenEvenG8(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(71))
       call FaPenEvenG9(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2mu
     -   p2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,p
     -   4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mup
     -   4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2mu
     -   p4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(81))
       call FaPenEvenG10(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mu
     -   p4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2m
     -   up4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(91))
       call FaPenEvenG11(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mu
     -   p4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2m
     -   up4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(101))
       call FaPenEvenG12(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mu
     -   p4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2m
     -   up4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(111))
       call FaPenEvenG13(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mu
     -   p4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2m
     -   up4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(121))
       call FaPenEvenG14(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mu
     -   p4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2m
     -   up4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(131))
       call FaPenEvenG15(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup5,
     -   p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5mu
     -   p4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup2m
     -   up4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(141))
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
c       print*,("ngluon",cpl(ngluon),ngluon=1,10)

       result(1,1) = Cpl(8)*(Cpl(9)*(Cpl(10)*Fa(127)+M1*M4*Cpl(5)*Fa(1
     -   28))+M3*Cpl(4)*(M4*Cpl(10)*Fa(129)+M1*Cpl(5)*Fa(130)))+Cpl(3
     -   )*(Cpl(4)*(Cpl(5)*Fa(127)+M1*M4*Cpl(10)*Fa(128))+M3*Cpl(9)*(
     -   M4*Cpl(5)*Fa(129)+M1*Cpl(10)*Fa(130)))
       result(1,2) = Cpl(8)*(Cpl(9)*(Cpl(10)*Fa(131)+M1*M4*Cpl(5)*Fa(1
     -   32))+M3*Cpl(4)*(M4*Cpl(10)*Fa(133)+M1*Cpl(5)*Fa(134)))+Cpl(3
     -   )*(Cpl(4)*(Cpl(5)*Fa(131)+M1*M4*Cpl(10)*Fa(132))+M3*Cpl(9)*(
     -   M4*Cpl(5)*Fa(133)+M1*Cpl(10)*Fa(134)))
       result(1,3) = Cpl(8)*(Cpl(9)*(Cpl(10)*Fa(135)+M1*M4*Cpl(5)*Fa(1
     -   36))+M3*Cpl(4)*(M4*Cpl(10)*Fa(137)+M1*Cpl(5)*Fa(138)))+Cpl(3
     -   )*(Cpl(4)*(Cpl(5)*Fa(135)+M1*M4*Cpl(10)*Fa(136))+M3*Cpl(9)*(
     -   M4*Cpl(5)*Fa(137)+M1*Cpl(10)*Fa(138)))
       result(1,4) = (Cpl(8)*(Cpl(9)*(2*Cpl(10)*Fa(139)+M1*M4*Cpl(5)*F
     -   a(140))+M3*Cpl(4)*(M4*Cpl(10)*Fa(141)+M1*Cpl(5)*Fa(142)))+Cp
     -   l(3)*(Cpl(4)*(2*Cpl(5)*Fa(139)+M1*M4*Cpl(10)*Fa(140))+M3*Cpl
     -   (9)*(M4*Cpl(5)*Fa(141)+M1*Cpl(10)*Fa(142))))/2.
       result(1,5) = (Cpl(8)*(Cpl(9)*(2*Cpl(10)*Fa(143)+M1*M4*Cpl(5)*F
     -   a(144))+M3*Cpl(4)*(M4*Cpl(10)*Fa(145)+M1*Cpl(5)*Fa(146)))+Cp
     -   l(3)*(Cpl(4)*(2*Cpl(5)*Fa(143)+M1*M4*Cpl(10)*Fa(144))+M3*Cpl
     -   (9)*(M4*Cpl(5)*Fa(145)+M1*Cpl(10)*Fa(146))))/2.
       result(1,6) = (Cpl(8)*(Cpl(9)*(2*Cpl(10)*Fa(147)+M1*M4*Cpl(5)*F
     -   a(148))+M3*Cpl(4)*(M4*Cpl(10)*Fa(149)+M1*Cpl(5)*Fa(150)))+Cp
     -   l(3)*(Cpl(4)*(2*Cpl(5)*Fa(147)+M1*M4*Cpl(10)*Fa(148))+M3*Cpl
     -   (9)*(M4*Cpl(5)*Fa(149)+M1*Cpl(10)*Fa(150))))/2.
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
