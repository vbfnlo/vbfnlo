       subroutine PenBoxFFhl1(F1)
       IMPLICIT NONE
       Real*8 P(48),FRe(1:43),FIm(1:43)
       Complex*16 F1(1:43)
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

       COMMON/PenBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s4
     &   5,s15
       COMMON/PenBoxPFunctions/P
       COMMON/PenBoxEVALFunctions/ EE0R,EijR,EE0I,EijI 
       COMMON/PenBoxDVALFunctions/ D01234R,Dij1234R,D01234I,Dij1234I, 
     &   D01235R,Dij1235R,D01235I,Dij1235I, D01245R,Dij1245R,D01245I,
     &   Dij1245I, D01345R,Dij1345R,D01345I,Dij1345I, D02345R,Dij2345
     &   R,D02345I,Dij2345I 
       COMMON/PenBoxCVALFunctions/ C0123R,Cij123R,C0123I,Cij123I, C012
     &   4R,Cij124R,C0124I,Cij124I, C0125R,Cij125R,C0125I,Cij125I, C0
     &   134R,Cij134R,C0134I,Cij134I, C0135R,Cij135R,C0135I,Cij135I, 
     &   C0145R,Cij145R,C0145I,Cij145I, C0234R,Cij234R,C0234I,Cij234I
     &   , C0235R,Cij235R,C0235I,Cij235I, C0245R,Cij245R,C0245I,Cij24
     &   5I, C0345R,Cij345R,C0345I,Cij345I 
       COMMON/PenBoxBVALFunctions/ B012R,B012I, B013R,B013I, B014R,B01
     &   4I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,B025I, B034
     &   R,B034I, B035R,B035I, B045R,B045I
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       P(1) = s12+s34
       P(2) = s12+s15
       P(3) = s15-s34
       FRe(1) = 4*(s12*(Dij1235R(1,1)+Dij1235R(2,1)+2*Dij1235R(4,2))+8
     &   *Dij1235R(7,2)-Dij1235R(3,1)*P(1)+2*(s34*Dij1235R(3,2)-Dij12
     &   35R(5,2)*P(2)+Dij1235R(6,2)*P(3)))
       FIm(1) = 4*(s12*(Dij1235I(1,1)+Dij1235I(2,1)+2*Dij1235I(4,2))+8
     &   *Dij1235I(7,2)-Dij1235I(3,1)*P(1)+2*(s34*Dij1235I(3,2)-Dij12
     &   35I(5,2)*P(2)+Dij1235I(6,2)*P(3)))
       F1(1)=DCMPLX(FRe(1),FIm(1))
       FRe(2) = -4*(s34*(Dij1235R(3,2)-Dij1235R(6,2))+6*Dij1235R(7,2)+
     &   (Dij1235R(2,1)-Dij1235R(3,1)+Dij1235R(4,2)-Dij1235R(5,2))*P(
     &   2))
       FIm(2) = -4*(s34*(Dij1235I(3,2)-Dij1235I(6,2))+6*Dij1235I(7,2)+
     &   (Dij1235I(2,1)-Dij1235I(3,1)+Dij1235I(4,2)-Dij1235I(5,2))*P(
     &   2))
       F1(2)=DCMPLX(FRe(2),FIm(2))
       FRe(3) = 4*s12*(Dij1235R(1,1)+Dij1235R(4,2)+Dij1235R(5,2)-2*Dij
     &   1235R(6,2))
       FIm(3) = 4*s12*(Dij1235I(1,1)+Dij1235I(4,2)+Dij1235I(5,2)-2*Dij
     &   1235I(6,2))
       F1(3)=DCMPLX(FRe(3),FIm(3))
       FRe(4) = 4*(2*(s34*Dij1235R(3,2)-s34*Dij1235R(5,2))+s12*(Dij123
     &   5R(1,1)+Dij1235R(2,1)+2*(Dij1235R(4,2)-Dij1235R(6,2)))+6*Dij
     &   1235R(7,2)-Dij1235R(3,1)*P(1))
       FIm(4) = 4*(2*(s34*Dij1235I(3,2)-s34*Dij1235I(5,2))+s12*(Dij123
     &   5I(1,1)+Dij1235I(2,1)+2*(Dij1235I(4,2)-Dij1235I(6,2)))+6*Dij
     &   1235I(7,2)-Dij1235I(3,1)*P(1))
       F1(4)=DCMPLX(FRe(4),FIm(4))
       FRe(5) = -4*(2*(s34*Dij1235R(3,2)-s34*Dij1235R(5,2))+s12*(Dij12
     &   35R(1,1)+Dij1235R(2,1)+2*(Dij1235R(4,2)-Dij1235R(6,2)))+8*Di
     &   j1235R(7,2)-Dij1235R(3,1)*P(1))
       FIm(5) = -4*(2*(s34*Dij1235I(3,2)-s34*Dij1235I(5,2))+s12*(Dij12
     &   35I(1,1)+Dij1235I(2,1)+2*(Dij1235I(4,2)-Dij1235I(6,2)))+8*Di
     &   j1235I(7,2)-Dij1235I(3,1)*P(1))
       F1(5)=DCMPLX(FRe(5),FIm(5))
       P(4) = 2*s12+s15-s34
       P(5) = 2*s12+s15
       P(6) = 3*s12+s15-s34
       P(7) = s12+s15+s34
       FRe(6) = -2*(s12*Dij1235R(1,1)+2*(s34*Dij1235R(3,2)-s12*Dij1235
     &   R(6,2))+6*Dij1235R(7,2)+Dij1235R(2,1)*P(4)-Dij1235R(3,1)*P(5
     &   )+Dij1235R(4,2)*P(6)-Dij1235R(5,2)*P(7))
       FIm(6) = -2*(s12*Dij1235I(1,1)+2*(s34*Dij1235I(3,2)-s12*Dij1235
     &   I(6,2))+6*Dij1235I(7,2)+Dij1235I(2,1)*P(4)-Dij1235I(3,1)*P(5
     &   )+Dij1235I(4,2)*P(6)-Dij1235I(5,2)*P(7))
       F1(6)=DCMPLX(FRe(6),FIm(6))
       FRe(7) = -4*(s12*(Dij1235R(1,1)-Dij1235R(2,1)+Dij1235R(3,1)+2*(
     &   Dij1235R(5,2)-Dij1235R(6,2)))+2*Dij1235R(7,2))
       FIm(7) = -4*(s12*(Dij1235I(1,1)-Dij1235I(2,1)+Dij1235I(3,1)+2*(
     &   Dij1235I(5,2)-Dij1235I(6,2)))+2*Dij1235I(7,2))
       F1(7)=DCMPLX(FRe(7),FIm(7))
       FRe(8) = -2*s12*(Dij1235R(1,1)+Dij1235R(4,2)+Dij1235R(5,2)-2*Di
     &   j1235R(6,2))+12*Dij1235R(7,2)
       FIm(8) = -2*s12*(Dij1235I(1,1)+Dij1235I(4,2)+Dij1235I(5,2)-2*Di
     &   j1235I(6,2))+12*Dij1235I(7,2)
       F1(8)=DCMPLX(FRe(8),FIm(8))
       FRe(9) = -Dij1235R(5,2)+Dij1235R(6,2)
       FIm(9) = -Dij1235I(5,2)+Dij1235I(6,2)
       F1(9)=DCMPLX(FRe(9),FIm(9))
       FRe(10) = -8*(Dij1235R(2,1)-Dij1235R(3,1)+Dij1235R(4,2)+Dij1235
     &   R(5,2)-2*Dij1235R(6,2))
       FIm(10) = -8*(Dij1235I(2,1)-Dij1235I(3,1)+Dij1235I(4,2)+Dij1235
     &   I(5,2)-2*Dij1235I(6,2))
       F1(10)=DCMPLX(FRe(10),FIm(10))
       FRe(11) = 4*(Dij1235R(2,1)-Dij1235R(3,1)+Dij1235R(4,2)+Dij1235R
     &   (5,2)-2*Dij1235R(6,2))
       FIm(11) = 4*(Dij1235I(2,1)-Dij1235I(3,1)+Dij1235I(4,2)+Dij1235I
     &   (5,2)-2*Dij1235I(6,2))
       F1(11)=DCMPLX(FRe(11),FIm(11))
       FRe(12) = Dij1235R(7,2)
       FIm(12) = Dij1235I(7,2)
       F1(12)=DCMPLX(FRe(12),FIm(12))
       P(8) = -s12-s15
       FRe(13) = s34*Dij1235R(3,2)-s34*Dij1235R(6,2)+6*Dij1235R(7,2)+(
     &   Dij1235R(2,1)+Dij1235R(4,2))*P(2)+(Dij1235R(3,1)+Dij1235R(5,
     &   2))*P(8)
       FIm(13) = s34*Dij1235I(3,2)-s34*Dij1235I(6,2)+6*Dij1235I(7,2)+(
     &   Dij1235I(2,1)+Dij1235I(4,2))*P(2)+(Dij1235I(3,1)+Dij1235I(5,
     &   2))*P(8)
       F1(13)=DCMPLX(FRe(13),FIm(13))
       FRe(14) = Dij1235R(2,1)-Dij1235R(3,1)+Dij1235R(4,2)-Dij1235R(5,
     &   2)
       FIm(14) = Dij1235I(2,1)-Dij1235I(3,1)+Dij1235I(4,2)-Dij1235I(5,
     &   2)
       F1(14)=DCMPLX(FRe(14),FIm(14))
       FRe(15) = -Dij1235R(2,1)+Dij1235R(3,1)-Dij1235R(4,2)+Dij1235R(5
     &   ,2)
       FIm(15) = -Dij1235I(2,1)+Dij1235I(3,1)-Dij1235I(4,2)+Dij1235I(5
     &   ,2)
       F1(15)=DCMPLX(FRe(15),FIm(15))
       FRe(16) = 2*(s34*(Dij1235R(3,2)-Dij1235R(6,2))+6*Dij1235R(7,2)+
     &   (Dij1235R(2,1)-Dij1235R(3,1)+Dij1235R(4,2)-Dij1235R(5,2))*P(
     &   2))
       FIm(16) = 2*(s34*(Dij1235I(3,2)-Dij1235I(6,2))+6*Dij1235I(7,2)+
     &   (Dij1235I(2,1)-Dij1235I(3,1)+Dij1235I(4,2)-Dij1235I(5,2))*P(
     &   2))
       F1(16)=DCMPLX(FRe(16),FIm(16))
       P(9) = s23-s45
       P(10) = s15-s23
       P(11) = s15+s23-s45
       P(12) = s15-s23+s45
       P(13) = -2*s15+s23
       FRe(17) = 4*(-(s23*Dij1245R(2,2))-s15*Dij1245R(3,1)-2*Dij1245R(
     &   7,2)+Dij1245R(1,1)*P(9)+Dij1245R(2,1)*P(10)+Dij1245R(4,2)*P(
     &   11)+Dij1245R(5,2)*P(12)+Dij1245R(6,2)*P(13))
       FIm(17) = 4*(-(s23*Dij1245I(2,2))-s15*Dij1245I(3,1)-2*Dij1245I(
     &   7,2)+Dij1245I(1,1)*P(9)+Dij1245I(2,1)*P(10)+Dij1245I(4,2)*P(
     &   11)+Dij1245I(5,2)*P(12)+Dij1245I(6,2)*P(13))
       F1(17)=DCMPLX(FRe(17),FIm(17))
       FRe(18) = -2*s45*(Dij1245R(1,1)+Dij1245R(2,1)-Dij1245R(3,1)+2*(
     &   Dij1245R(4,2)-Dij1245R(6,2)))
       FIm(18) = -2*s45*(Dij1245I(1,1)+Dij1245I(2,1)-Dij1245I(3,1)+2*(
     &   Dij1245I(4,2)-Dij1245I(6,2)))
       F1(18)=DCMPLX(FRe(18),FIm(18))
       P(14) = s23+s45
       FRe(19) = -4*(s23*Dij1245R(2,2)-s45*Dij1245R(6,2))+2*(s45*Dij12
     &   45R(3,1)+(Dij1245R(1,1)+2*Dij1245R(4,2))*P(9)-Dij1245R(2,1)*
     &   P(14))
       FIm(19) = -4*(s23*Dij1245I(2,2)-s45*Dij1245I(6,2))+2*(s45*Dij12
     &   45I(3,1)+(Dij1245I(1,1)+2*Dij1245I(4,2))*P(9)-Dij1245I(2,1)*
     &   P(14))
       F1(19)=DCMPLX(FRe(19),FIm(19))
       FRe(20) = -4*(Dij1245R(2,1)-Dij1245R(3,1)+Dij1245R(4,2)+Dij1245
     &   R(5,2)-2*Dij1245R(6,2))
       FIm(20) = -4*(Dij1245I(2,1)-Dij1245I(3,1)+Dij1245I(4,2)+Dij1245
     &   I(5,2)-2*Dij1245I(6,2))
       F1(20)=DCMPLX(FRe(20),FIm(20))
       FRe(21) = Dij1245R(7,2)
       FIm(21) = Dij1245I(7,2)
       F1(21)=DCMPLX(FRe(21),FIm(21))
       P(15) = -s15-s45
       P(16) = s15+s45
       P(17) = -s15+s23-s45
       FRe(22) = -(s23*Dij1245R(2,2))+s23*Dij1245R(6,2)-6*Dij1245R(7,2
     &   )+Dij1245R(5,2)*P(12)+Dij1245R(2,1)*P(15)+Dij1245R(3,1)*P(16
     &   )+Dij1245R(4,2)*P(17)
       FIm(22) = -(s23*Dij1245I(2,2))+s23*Dij1245I(6,2)-6*Dij1245I(7,2
     &   )+Dij1245I(5,2)*P(12)+Dij1245I(2,1)*P(15)+Dij1245I(3,1)*P(16
     &   )+Dij1245I(4,2)*P(17)
       F1(22)=DCMPLX(FRe(22),FIm(22))
       FRe(23) = Dij1245R(2,1)-Dij1245R(3,1)+Dij1245R(4,2)-Dij1245R(5,
     &   2)
       FIm(23) = Dij1245I(2,1)-Dij1245I(3,1)+Dij1245I(4,2)-Dij1245I(5,
     &   2)
       F1(23)=DCMPLX(FRe(23),FIm(23))
       FRe(24) = 2*(s23*Dij1245R(2,2)-s23*Dij1245R(6,2)+6*Dij1245R(7,2
     &   )+(Dij1245R(4,2)-Dij1245R(5,2))*P(12)+(Dij1245R(2,1)-Dij1245
     &   R(3,1))*P(16))
       FIm(24) = 2*(s23*Dij1245I(2,2)-s23*Dij1245I(6,2)+6*Dij1245I(7,2
     &   )+(Dij1245I(4,2)-Dij1245I(5,2))*P(12)+(Dij1245I(2,1)-Dij1245
     &   I(3,1))*P(16))
       F1(24)=DCMPLX(FRe(24),FIm(24))
       P(18) = s12-s15
       FRe(25) = -4*(s12*(EijR(1,1)-EijR(2,1)-EijR(7,2)-EijR(8,2)+EijR
     &   (9,2))+s15*(-EijR(3,1)+EijR(4,1)+EijR(10,2))+2*(EijR(11,2)+E
     &   ijR(23,3)-EijR(24,3))+EijR(6,2)*P(18))
       FIm(25) = -4*(s12*(EijI(1,1)-EijI(2,1)-EijI(7,2)-EijI(8,2)+EijI
     &   (9,2))+s15*(-EijI(3,1)+EijI(4,1)+EijI(10,2))+2*(EijI(11,2)+E
     &   ijI(23,3)-EijI(24,3))+EijI(6,2)*P(18))
       F1(25)=DCMPLX(FRe(25),FIm(25))
       P(19) = p3sq+s12+2*s15-s23-s34
       P(20) = p3sq-s23-s34
       P(21) = s12+4*s15+s45+2*P(20)
       P(22) = p3sq+s15-s23-s34
       P(23) = s12+2*P(22)
       P(24) = -p3sq-s15+s23+s34
       P(25) = s12+2*s15
       P(26) = 2*s12+6*s15+s45+4*P(20)
       P(27) = s45+2*P(22)
       P(28) = s45+4*P(22)
       FRe(26) = 4*(s12*(EijR(5,2)-EijR(6,2)+EijR(8,3))+2*(EijR(11,2)+
     &   EijR(22,3)-EijR(23,3)-EijR(17,3)*P(2))+(EijR(3,2)-EijR(8,2))
     &   *P(19)-EijR(3,3)*P(21)+(EijR(9,2)+2*EijR(10,3))*P(22)-EijR(9
     &   ,3)*P(23)+EijR(10,2)*P(24)+EijR(11,3)*P(25)+EijR(12,3)*P(26)
     &   +EijR(13,3)*P(27)-EijR(20,3)*P(28))
       FIm(26) = 4*(s12*(EijI(5,2)-EijI(6,2)+EijI(8,3))+2*(EijI(11,2)+
     &   EijI(22,3)-EijI(23,3)-EijI(17,3)*P(2))+(EijI(3,2)-EijI(8,2))
     &   *P(19)-EijI(3,3)*P(21)+(EijI(9,2)+2*EijI(10,3))*P(22)-EijI(9
     &   ,3)*P(23)+EijI(10,2)*P(24)+EijI(11,3)*P(25)+EijI(12,3)*P(26)
     &   +EijI(13,3)*P(27)-EijI(20,3)*P(28))
       F1(26)=DCMPLX(FRe(26),FIm(26))
       P(29) = p3sq+2*s15-s23-s34+s45
       P(30) = s23+s34
       P(31) = 3*s15+s45-2*P(30)
       P(32) = 4*p3sq+s12+2*P(31)
       P(33) = s12+4*P(22)
       FRe(27) = -4*(-(s15*EijR(6,2))+s15*EijR(7,2)+s12*EijR(17,3)-s12
     &   *EijR(18,3)+EijR(3,3)*P(21)+(-EijR(8,2)+EijR(9,2)-2*EijR(15,
     &   3))*P(22)-EijR(12,3)*P(23)+(-EijR(11,3)+EijR(19,3))*P(25)+Ei
     &   jR(16,3)*P(27)+(EijR(3,2)-EijR(10,2))*P(29)-EijR(13,3)*P(32)
     &   +EijR(20,3)*P(33))
       FIm(27) = -4*(-(s15*EijI(6,2))+s15*EijI(7,2)+s12*EijI(17,3)-s12
     &   *EijI(18,3)+EijI(3,3)*P(21)+(-EijI(8,2)+EijI(9,2)-2*EijI(15,
     &   3))*P(22)-EijI(12,3)*P(23)+(-EijI(11,3)+EijI(19,3))*P(25)+Ei
     &   jI(16,3)*P(27)+(EijI(3,2)-EijI(10,2))*P(29)-EijI(13,3)*P(32)
     &   +EijI(20,3)*P(33))
       F1(27)=DCMPLX(FRe(27),FIm(27))
       P(34) = p3sq+s15-s23-s34+s45
       P(35) = p3sq+s12-s23-s34
       P(36) = 6*s15+s45+2*P(35)
       P(37) = p3sq+s12+s15-s23-s34
       FRe(28) = 4*(D02345R+Dij2345R(2,1)+s12*EijR(5,3)+s45*EijR(7,2)+
     &   EijR(6,2)*P(20)-EijR(3,3)*P(21)+EijR(12,3)*P(23)+(EijR(2,1)+
     &   EijR(5,2)+EijR(8,2))*P(24)-EijR(6,3)*P(25)+(EijR(13,3)-EijR(
     &   19,3))*P(27)+EijR(3,2)*P(29)+(EijR(3,1)-2*EijR(10,2))*P(34)+
     &   EijR(11,3)*P(36)-2*(EijR(11,2)-EijR(21,3)+EijR(23,3)-(EijR(9
     &   ,2)+EijR(18,3)-EijR(20,3))*P(22)+EijR(17,3)*P(37)))
       FIm(28) = 4*(D02345I+Dij2345I(2,1)+s12*EijI(5,3)+s45*EijI(7,2)+
     &   EijI(6,2)*P(20)-EijI(3,3)*P(21)+EijI(12,3)*P(23)+(EijI(2,1)+
     &   EijI(5,2)+EijI(8,2))*P(24)-EijI(6,3)*P(25)+(EijI(13,3)-EijI(
     &   19,3))*P(27)+EijI(3,2)*P(29)+(EijI(3,1)-2*EijI(10,2))*P(34)+
     &   EijI(11,3)*P(36)-2*(EijI(11,2)-EijI(21,3)+EijI(23,3)-(EijI(9
     &   ,2)+EijI(18,3)-EijI(20,3))*P(22)+EijI(17,3)*P(37)))
       F1(28)=DCMPLX(FRe(28),FIm(28))
       FRe(29) = 4*(s15*EijR(6,2)-s12*EijR(17,3)-EijR(3,3)*P(21)+(EijR
     &   (8,2)-2*EijR(20,3))*P(22)+EijR(12,3)*P(23)+EijR(11,3)*P(25)+
     &   EijR(13,3)*P(27)-EijR(3,2)*P(29))
       FIm(29) = 4*(s15*EijI(6,2)-s12*EijI(17,3)-EijI(3,3)*P(21)+(EijI
     &   (8,2)-2*EijI(20,3))*P(22)+EijI(12,3)*P(23)+EijI(11,3)*P(25)+
     &   EijI(13,3)*P(27)-EijI(3,2)*P(29))
       F1(29)=DCMPLX(FRe(29),FIm(29))
       P(38) = -p3sq+s12+s23
       P(39) = -p3sq+s12+s23-2*s45
       P(40) = s12-s45
       P(41) = p3sq-s23
       P(42) = s12+s15-s34+s45+2*P(41)
       P(43) = p3sq+s15-s23-s34+2*s45
       FRe(30) = -2*(D02345R+Dij2345R(2,1)-s12*EijR(1,2)-s45*EijR(4,1)
     &   -2*(EijR(11,2)-EijR(21,3)+EijR(23,3))-(EijR(7,2)+EijR(8,2)-2
     &   *EijR(9,2))*P(22)+(EijR(2,1)+EijR(5,2))*P(38)-EijR(3,1)*P(39
     &   )-EijR(3,2)*P(40)+EijR(6,2)*P(42)-EijR(10,2)*P(43))
       FIm(30) = -2*(D02345I+Dij2345I(2,1)-s12*EijI(1,2)-s45*EijI(4,1)
     &   -2*(EijI(11,2)-EijI(21,3)+EijI(23,3))-(EijI(7,2)+EijI(8,2)-2
     &   *EijI(9,2))*P(22)+(EijI(2,1)+EijI(5,2))*P(38)-EijI(3,1)*P(39
     &   )-EijI(3,2)*P(40)+EijI(6,2)*P(42)-EijI(10,2)*P(43))
       F1(30)=DCMPLX(FRe(30),FIm(30))
       P(44) = -s12+s45
       FRe(31) = 4*(Dij2345R(2,1)+s12*EijR(8,2)-2*EijR(23,3)+EijR(10,2
     &   )*P(22)+EijR(9,2)*P(24)+EijR(3,2)*P(44))
       FIm(31) = 4*(Dij2345I(2,1)+s12*EijI(8,2)-2*EijI(23,3)+EijI(10,2
     &   )*P(22)+EijI(9,2)*P(24)+EijI(3,2)*P(44))
       F1(31)=DCMPLX(FRe(31),FIm(31))
       P(45) = p3sq+2*s12+s15-s23-s34
       FRe(32) = 2*(Dij2345R(2,1)-s12*(EijR(1,1)+EijR(5,2))+8*EijR(11,
     &   2)-2*EijR(23,3)+(-2*EijR(9,2)+EijR(10,2))*P(22)+EijR(4,1)*P(
     &   24)+EijR(3,1)*P(37)+EijR(3,2)*P(44)+EijR(8,2)*P(45))
       FIm(32) = 2*(Dij2345I(2,1)-s12*(EijI(1,1)+EijI(5,2))+8*EijI(11,
     &   2)-2*EijI(23,3)+(-2*EijI(9,2)+EijI(10,2))*P(22)+EijI(4,1)*P(
     &   24)+EijI(3,1)*P(37)+EijI(3,2)*P(44)+EijI(8,2)*P(45))
       F1(32)=DCMPLX(FRe(32),FIm(32))
       FRe(33) = 8*(EijR(3,2)-EijR(3,3)-EijR(8,2)-EijR(9,3)+2*EijR(12,
     &   3))
       FIm(33) = 8*(EijI(3,2)-EijI(3,3)-EijI(8,2)-EijI(9,3)+2*EijI(12,
     &   3))
       F1(33)=DCMPLX(FRe(33),FIm(33))
       FRe(34) = -8*(EijR(3,3)-EijR(12,3))
       FIm(34) = -8*(EijI(3,3)-EijI(12,3))
       F1(34)=DCMPLX(FRe(34),FIm(34))
       FRe(35) = EijR(3,2)-EijR(3,3)-EijR(8,2)+EijR(11,3)+EijR(12,3)-E
     &   ijR(17,3)
       FIm(35) = EijI(3,2)-EijI(3,3)-EijI(8,2)+EijI(11,3)+EijI(12,3)-E
     &   ijI(17,3)
       F1(35)=DCMPLX(FRe(35),FIm(35))
       FRe(36) = EijR(3,3)-EijR(12,3)-EijR(13,3)+EijR(20,3)
       FIm(36) = EijI(3,3)-EijI(12,3)-EijI(13,3)+EijI(20,3)
       F1(36)=DCMPLX(FRe(36),FIm(36))
       FRe(37) = 4*(-EijR(3,2)+EijR(6,2)+2*(EijR(3,3)+EijR(6,3))-4*Eij
     &   R(11,3))
       FIm(37) = 4*(-EijI(3,2)+EijI(6,2)+2*(EijI(3,3)+EijI(6,3))-4*Eij
     &   I(11,3))
       F1(37)=DCMPLX(FRe(37),FIm(37))
       FRe(38) = 4*(EijR(3,2)-EijR(6,2)+2*(EijR(3,3)-EijR(11,3)))
       FIm(38) = 4*(EijI(3,2)-EijI(6,2)+2*(EijI(3,3)-EijI(11,3)))
       F1(38)=DCMPLX(FRe(38),FIm(38))
       FRe(39) = 4*(-EijR(3,2)+EijR(8,2)+2*(EijR(3,3)-EijR(11,3)-EijR(
     &   12,3)+EijR(17,3)))
       FIm(39) = 4*(-EijI(3,2)+EijI(8,2)+2*(EijI(3,3)-EijI(11,3)-EijI(
     &   12,3)+EijI(17,3)))
       F1(39)=DCMPLX(FRe(39),FIm(39))
       FRe(40) = 4*(EijR(3,2)-EijR(6,2)+EijR(7,2)-EijR(10,2)+2*(EijR(3
     &   ,3)-EijR(11,3)-EijR(13,3)+EijR(19,3)))
       FIm(40) = 4*(EijI(3,2)-EijI(6,2)+EijI(7,2)-EijI(10,2)+2*(EijI(3
     &   ,3)-EijI(11,3)-EijI(13,3)+EijI(19,3)))
       F1(40)=DCMPLX(FRe(40),FIm(40))
       FRe(41) = -EijR(7,2)+EijR(10,2)
       FIm(41) = -EijI(7,2)+EijI(10,2)
       F1(41)=DCMPLX(FRe(41),FIm(41))
       FRe(42) = -EijR(9,2)+EijR(10,2)
       FIm(42) = -EijI(9,2)+EijI(10,2)
       F1(42)=DCMPLX(FRe(42),FIm(42))
       FRe(43) = 4*(EijR(3,2)-EijR(8,2)+EijR(9,2)-EijR(10,2)+2*(EijR(3
     &   ,3)-EijR(12,3)-EijR(15,3)+EijR(16,3)-2*(EijR(13,3)-EijR(20,3
     &   ))))
       FIm(43) = 4*(EijI(3,2)-EijI(8,2)+EijI(9,2)-EijI(10,2)+2*(EijI(3
     &   ,3)-EijI(12,3)-EijI(15,3)+EijI(16,3)-2*(EijI(13,3)-EijI(20,3
     &   ))))
       F1(43)=DCMPLX(FRe(43),FIm(43))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
