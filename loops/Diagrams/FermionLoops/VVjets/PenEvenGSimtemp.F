       subroutine PenEvenGSim(cpl,MM1,MM3,MM4,p1,p2,p3,p4,p5, mup1, mu
     -   p2,mup3,mup4,mup5,musq,comp,result)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/02/2012
c Modified:17/7/2012
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
       Integer comp 
       Real*8  MM1,MM3,MM4,M1,M3,M4,m1sq,m3sq,m4sq 
       Real*8 cpl(10) 
      Complex*16 F(3806),Fa(123)
      Real*8 P(1320)
      Complex*16 dotrc,dotcc,result(1,6)
       Real*8 musq
      EXTERNAL   dotrc,dotcc

       COMMON/FaPenEvenGSim/Fa
       COMMON/FPenEvenGSim/F
       COMMON/PPenEvenGSim/P
      Save/FPenEvenGSim/
      Save/PPenEvenGSim/
       COMMON/InvariantsPenEvenGSim/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s
     -   34,s45,s234,m1,m1sq,m3,m3sq,m4,m4sq
       COMMON/EVALFunctionsPenEvenGSim/ E012345R,E12345R,E012345I,E123
     -   45I
       COMMON/DVALFunctionsPenEvenGSim/ D01234R,D1234R,D01234I,D1234I,
     -    D01235R,D1235R,D01235I,D1235I, D01245R,D1245R,D01245I,D1245
     -   I, D01345R,D1345R,D01345I,D1345I, D02345R,D2345R,D02345I,D23
     -   45I
       COMMON/CVALFunctionsPenEvenGSim/ C0123R,C123R,C0123I,C123I, C01
     -   24R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0134R,C1
     -   34R,C0134I,C134I, C0135R,C135R,C0135I,C135I, C0145R,C145R,C0
     -   145I,C145I, C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C
     -   235I, C0245R,C245R,C0245I,C245I, C0345R,C345R,C0345I,C345I 
       COMMON/BVALFunctionsPenEvenGSim/ B012R,B012I, B013R,B013I, B014
     -   R,B014I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,B025I,
     -    B034R,B034I, B035R,B035I, B045R,B045I, B12R,B12I, B13R,B13I
     -   , B14R,B14I, B15R,B15I, B23R,B23I, B24R,B24I, B25R,B25I, B34
     -   R,B34I, B35R,B35I, B45R,B45I 


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


c       Write(*,'(a5,F20.10)')," m1sq",m1sq 
c       Write(*,'(a5,F20.10)')," m3sq",m3sq
c       Write(*,'(a5,F20.10)')," m4sq",m4sq
       
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
       call  XEgetSim_R5(M1,M1,M1,M3,M4,p1sq,p2sq,p3sq,p4sq,p5sq,p1
     -   p2,p1p3,p1p4,p2p3,p2p4,p3p4,p1p5,p4p5,musq)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call FPenEvenGSim1(F(1))
       call FPenEvenGSim2(F(191))
       call FPenEvenGSim3(F(381))
       call FPenEvenGSim4(F(571))
       call FPenEvenGSim5(F(761))
       call FPenEvenGSim6(F(951))
       call FPenEvenGSim7(F(1141))
       call FPenEvenGSim8(F(1331))
       call FPenEvenGSim9(F(1521))
       call FPenEvenGSim10(F(1711))
       call FPenEvenGSim11(F(1901))
       call FPenEvenGSim12(F(2091))
       call FPenEvenGSim13(F(2281))
       call FPenEvenGSim14(F(2471))
       call FPenEvenGSim15(F(2661))
       call FPenEvenGSim16(F(2851))
       call FPenEvenGSim17(F(3041))
       call FPenEvenGSim18(F(3231))
       call FPenEvenGSim19(F(3421))
       call FPenEvenGSim20(F(3611))
cfc       do j1=1,1320
cfc          print*,"P",j1,P(j1)
cfc       enddo
cfc       do j1=1,3806
cfc          print*,"F",j1,Dble(F(j1))
cfc          print*,"F",j1,DImag(F(j1))
cfc       enddo


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

cfc       print*, "p1mup1", p1mup1
cfc       print*, "p1mup2", p1mup2
cfc       print*, "p1mup3", p1mup3
cfc       print*, "p1mup4", p1mup4
cfc       print*, "p1mup5", p1mup5
cfc       print*, "p2mup1", p2mup1
cfc       print*, "p2mup2", p2mup2
cfc       print*, "p2mup3", p2mup3
cfc       print*, "p2mup4", p2mup4
cfc       print*, "p2mup5", p2mup5
cfc       print*, "p3mup1", p3mup1
cfc       print*, "p3mup2", p3mup2
cfc       print*, "p3mup3", p3mup3
cfc       print*, "p3mup4", p3mup4
cfc       print*, "p3mup5", p3mup5
cfc       print*, "p4mup1", p4mup1
cfc       print*, "p4mup2", p4mup2
cfc       print*, "p4mup3", p4mup3
cfc       print*, "p4mup4", p4mup4
cfc       print*, "p4mup5", p4mup5
cfc       print*, "p5mup1", p5mup1
cfc       print*, "p5mup2", p5mup2
cfc       print*, "p5mup3", p5mup3
cfc       print*, "p5mup4", p5mup4
cfc       print*, "p5mup5", p5mup5
cfc       print*, "mup1mup2",mup1mup2
cfc       print*, "mup1mup3",mup1mup3
cfc       print*, "mup1mup4",mup1mup4
cfc       print*, "mup1mup5",mup1mup5
cfc       print*, "mup2mup3",mup2mup3
cfc       print*, "mup2mup4",mup2mup4
cfc       print*, "mup2mup5",mup2mup5
cfc       print*, "mup3mup4",mup3mup4
cfc       print*, "mup3mup5",mup3mup5
cfc       print*, "mup4mup5",mup4mup5


c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call FaPenEvenGSim1(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(1))
       call FaPenEvenGSim2(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(9))
       call FaPenEvenGSim3(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(17))
       call FaPenEvenGSim4(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(25))
       call FaPenEvenGSim5(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(33))
       call FaPenEvenGSim6(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(41))
       if(comp.gt.0) then
       call FaPenEvenGSim7(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(49))
       call FaPenEvenGSim8(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(57))
       call FaPenEvenGSim9(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(65))
       call FaPenEvenGSim10(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(73))
       call FaPenEvenGSim11(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(81))
       call FaPenEvenGSim12(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(89))
        call FaPenEvenGSim13(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(97))
         endif
       call FaPenEvenGSim14(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(104))
       if(comp.gt.0) then
       call FaPenEvenGSim15(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,
     -   p2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mu
     -   p5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p
     -   5mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mu
     -   p2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa(113))
       endif


cfc      do j1=1,123
cfc          print*,"Fa",j1,Fa(j1)
cfc       enddo
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
       result(1,1) = Cpl(8)*(Cpl(9)*(Cpl(10)*Fa(104)+M1*M4*Cpl(5)*Fa(1
     -   05))+M3*Cpl(4)*(M4*Cpl(10)*Fa(106)+M1*Cpl(5)*Fa(107)))+Cpl(3
     -   )*(Cpl(4)*(Cpl(5)*Fa(104)+M1*M4*Cpl(10)*Fa(105))+M3*Cpl(9)*(
     -   M4*Cpl(5)*Fa(106)+M1*Cpl(10)*Fa(107)))
       if(comp.gt.0)then
       result(1,2) = Cpl(8)*(Cpl(9)*(Cpl(10)*Fa(108)+M1*M4*Cpl(5)*Fa(1
     -   09))+M3*Cpl(4)*(M4*Cpl(10)*Fa(110)+M1*Cpl(5)*Fa(111)))+Cpl(3
     -   )*(Cpl(4)*(Cpl(5)*Fa(108)+M1*M4*Cpl(10)*Fa(109))+M3*Cpl(9)*(
     -   M4*Cpl(5)*Fa(110)+M1*Cpl(10)*Fa(111)))
       result(1,3) = Cpl(8)*(Cpl(9)*(Cpl(10)*Fa(112)+M1*M4*Cpl(5)*Fa(1
     -   13))+M3*Cpl(4)*(M4*Cpl(10)*Fa(114)+M1*Cpl(5)*Fa(115)))+Cpl(3
     -   )*(Cpl(4)*(Cpl(5)*Fa(112)+M1*M4*Cpl(10)*Fa(113))+M3*Cpl(9)*(
     -   M4*Cpl(5)*Fa(114)+M1*Cpl(10)*Fa(115)))
       result(1,4) = (Cpl(8)*(Cpl(9)*(2*Cpl(10)*Fa(116)+M1*M4*Cpl(5)*F
     -   a(117))+M3*Cpl(4)*(M4*Cpl(10)*Fa(118)+M1*Cpl(5)*Fa(119)))+Cp
     -   l(3)*(Cpl(4)*(2*Cpl(5)*Fa(116)+M1*M4*Cpl(10)*Fa(117))+M3*Cpl
     -   (9)*(M4*Cpl(5)*Fa(118)+M1*Cpl(10)*Fa(119))))/2.
       result(1,5) = (Cpl(8)*(Cpl(9)*(2*Cpl(10)*Fa(120)+M1*M4*Cpl(5)*F
     -   a(121))+M3*Cpl(4)*(M4*Cpl(10)*Fa(122)+M1*Cpl(5)*Fa(123)))+Cp
     -   l(3)*(Cpl(4)*(2*Cpl(5)*Fa(120)+M1*M4*Cpl(10)*Fa(121))+M3*Cpl
     -   (9)*(M4*Cpl(5)*Fa(122)+M1*Cpl(10)*Fa(123))))/2.
       endif
       result(1,6) = 0
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
