       subroutine NoAbe4vertBoxLineFFhl1(F1)
       IMPLICIT NONE
       Real*8 P(3),FRe(1:21),FIm(1:21)
       Complex*16 F1(1:21)
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12, s15 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45


       Real*8  B012R,B013R,B014R,B015R 
       Real*8  B023R,B024R,B025R 
       Real*8  B034R,B035R 
       Real*8  B045R 
       Real*8  B012I,B013I,B014I,B015I 
       Real*8  B023I,B024I,B025I 
       Real*8  B034I,B035I 
       Real*8  B045I


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
       Real*8 Cij123R(4,2),Cij124R(4,2),Cij125R(4,2) 
       Real*8 Cij134R(4,2),Cij135R(4,2) 
       Real*8 Cij145R(4,2) 
       Real*8 Cij234R(4,2),Cij235R(4,2) 
       Real*8 Cij245R(4,2) 
       Real*8 Cij345R(4,2) 
       Real*8 Cij123I(4,2),Cij124I(4,2),Cij125I(4,2) 
       Real*8 Cij134I(4,2),Cij135I(4,2) 
       Real*8 Cij145I(4,2) 
       Real*8 Cij234I(4,2),Cij235I(4,2) 
       Real*8 Cij245I(4,2) 
       Real*8 Cij345I(4,2)


       Real*8 D01234R,D01235R 
       Real*8 D01245R 
       Real*8 D01345R 
       Real*8 D02345R 
       Real*8 D01234I,D01235I 
       Real*8 D01245I 
       Real*8 D01345I 
       Real*8 D02345I 
       Real*8 Dij1234R(13,3),Dij1235R(13,3) 
       Real*8 Dij1245R(13,3) 
       Real*8 Dij1345R(13,3) 
       Real*8 Dij2345R(13,3) 
       Real*8 Dij1234I(13,3),Dij1235I(13,3) 
       Real*8 Dij1245I(13,3) 
       Real*8 Dij1345I(13,3) 
       Real*8 Dij2345I(13,3) 

       Real*8 EE0R 
       Real*8 EE0I 
       Real*8 EijR(46,4) 
       Real*8 EijI(46,4) 

       COMMON/NoAbe4vertBoxLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12
     &   ,s23,s34,s45,s15
       COMMON/NoAbe4vertBoxLinePFunctions/P
       COMMON/NoAbe4vertBoxLineEVALFunctions/ EE0R,EijR,EE0I,EijI 
       COMMON/NoAbe4vertBoxLineDVALFunctions/ D01234R,Dij1234R,D01234I
     &   ,Dij1234I, D01235R,Dij1235R,D01235I,Dij1235I, D01245R,Dij124
     &   5R,D01245I,Dij1245I, D01345R,Dij1345R,D01345I,Dij1345I, D023
     &   45R,Dij2345R,D02345I,Dij2345I 
       COMMON/NoAbe4vertBoxLineCVALFunctions/ C0123R,Cij123R,C0123I,Ci
     &   j123I, C0124R,Cij124R,C0124I,Cij124I, C0125R,Cij125R,C0125I,
     &   Cij125I, C0134R,Cij134R,C0134I,Cij134I, C0135R,Cij135R,C0135
     &   I,Cij135I, C0145R,Cij145R,C0145I,Cij145I, C0234R,Cij234R,C02
     &   34I,Cij234I, C0235R,Cij235R,C0235I,Cij235I, C0245R,Cij245R,C
     &   0245I,Cij245I, C0345R,Cij345R,C0345I,Cij345I 
       COMMON/NoAbe4vertBoxLineBVALFunctions/ B012R,B012I, B013R,B013I
     &   , B014R,B014I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,
     &   B025I, B034R,B034I, B035R,B035I, B045R,B045I
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       FRe(1) = Cij124R(2,1)
       FIm(1) = Cij124I(2,1)
       F1(1)=DCMPLX(FRe(1),FIm(1))
       FRe(2) = C0134R+Cij134R(1,1)
       FIm(2) = C0134I+Cij134I(1,1)
       F1(2)=DCMPLX(FRe(2),FIm(2))
       P(1) = -p2sq+s12
       FRe(3) = 2*(C0134R+p2sq*Dij1234R(2,1)-2*Dij1234R(7,2)+(D01234R+
     &   Dij1234R(1,1))*P(1))
       FIm(3) = 2*(C0134I+p2sq*Dij1234I(2,1)-2*Dij1234I(7,2)+(D01234I+
     &   Dij1234I(1,1))*P(1))
       F1(3)=DCMPLX(FRe(3),FIm(3))
       FRe(4) = D01234R+Dij1234R(1,1)+Dij1234R(2,1)+Dij1234R(4,2)
       FIm(4) = D01234I+Dij1234I(1,1)+Dij1234I(2,1)+Dij1234I(4,2)
       F1(4)=DCMPLX(FRe(4),FIm(4))
       FRe(5) = Dij1234R(2,1)+Dij1234R(2,2)
       FIm(5) = Dij1234I(2,1)+Dij1234I(2,2)
       F1(5)=DCMPLX(FRe(5),FIm(5))
       FRe(6) = Dij1234R(6,2)
       FIm(6) = Dij1234I(6,2)
       F1(6)=DCMPLX(FRe(6),FIm(6))
       FRe(7) = (C0124R+Cij124R(1,1))/s23
       FIm(7) = (C0124I+Cij124I(1,1))/s23
       F1(7)=DCMPLX(FRe(7),FIm(7))
       FRe(8) = Cij134R(2,1)
       FIm(8) = Cij134I(2,1)
       F1(8)=DCMPLX(FRe(8),FIm(8))
       FRe(9) = (C0134R+Cij134R(1,1))/s12
       FIm(9) = (C0134I+Cij134I(1,1))/s12
       F1(9)=DCMPLX(FRe(9),FIm(9))
       FRe(10) = -C0134R-Cij134R(1,1)
       FIm(10) = -C0134I-Cij134I(1,1)
       F1(10)=DCMPLX(FRe(10),FIm(10))
       FRe(11) = 4*(D01234R+2*Dij1234R(1,1)+Dij1234R(1,2))
       FIm(11) = 4*(D01234I+2*Dij1234I(1,1)+Dij1234I(1,2))
       F1(11)=DCMPLX(FRe(11),FIm(11))
       FRe(12) = Dij1234R(2,1)+Dij1234R(4,2)
       FIm(12) = Dij1234I(2,1)+Dij1234I(4,2)
       F1(12)=DCMPLX(FRe(12),FIm(12))
       FRe(13) = Dij1234R(3,1)+Dij1234R(5,2)
       FIm(13) = Dij1234I(3,1)+Dij1234I(5,2)
       F1(13)=DCMPLX(FRe(13),FIm(13))
       FRe(14) = Dij1234R(3,1)+Dij1234R(6,2)
       FIm(14) = Dij1234I(3,1)+Dij1234I(6,2)
       F1(14)=DCMPLX(FRe(14),FIm(14))
       P(2) = -p2sq+s23
       P(3) = p2sq-s23
       FRe(15) = -2*(C0134R+p2sq*Dij1234R(2,1)-2*Dij1234R(7,2)-Dij1234
     &   R(3,1)*P(3))
       FIm(15) = -2*(C0134I+p2sq*Dij1234I(2,1)-2*Dij1234I(7,2)-Dij1234
     &   I(3,1)*P(3))
       F1(15)=DCMPLX(FRe(15),FIm(15))
       FRe(16) = Dij1234R(3,1)
       FIm(16) = Dij1234I(3,1)
       F1(16)=DCMPLX(FRe(16),FIm(16))
       FRe(17) = Dij1234R(3,2)
       FIm(17) = Dij1234I(3,2)
       F1(17)=DCMPLX(FRe(17),FIm(17))
       FRe(18) = C0134R+p2sq*Dij1234R(2,1)-2*Dij1234R(7,2)+Dij1234R(3,
     &   1)*P(2)
       FIm(18) = C0134I+p2sq*Dij1234I(2,1)-2*Dij1234I(7,2)+Dij1234I(3,
     &   1)*P(2)
       F1(18)=DCMPLX(FRe(18),FIm(18))
       FRe(19) = D01234R+Dij1234R(1,1)
       FIm(19) = D01234I+Dij1234I(1,1)
       F1(19)=DCMPLX(FRe(19),FIm(19))
       FRe(20) = 2*(C0134R+p2sq*Dij1234R(2,1)-2*Dij1234R(7,2)+Dij1234R
     &   (3,1)*P(2))
       FIm(20) = 2*(C0134I+p2sq*Dij1234I(2,1)-2*Dij1234I(7,2)+Dij1234I
     &   (3,1)*P(2))
       F1(20)=DCMPLX(FRe(20),FIm(20))
       FRe(21) = -C0134R-p2sq*Dij1234R(2,1)+2*Dij1234R(7,2)+Dij1234R(3
     &   ,1)*P(3)
       FIm(21) = -C0134I-p2sq*Dij1234I(2,1)+2*Dij1234I(7,2)+Dij1234I(3
     &   ,1)*P(3)
       F1(21)=DCMPLX(FRe(21),FIm(21))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
