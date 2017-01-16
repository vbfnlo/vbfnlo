       subroutine NoAbe2gaPenLineFFhl6(F6)
       IMPLICIT NONE
       Real*8 P(767),FRe(211:252),FIm(211:252)
       Complex*16 F6(211:252)
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

       COMMON/NoAbe2gaPenLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s
     &   23,s34,s45,s15
       COMMON/NoAbe2gaPenLinePFunctions/P
       COMMON/NoAbe2gaPenLineEVALFunctions/ EE0R,EijR,EE0I,EijI 
       COMMON/NoAbe2gaPenLineDVALFunctions/ D01234R,Dij1234R,D01234I,D
     &   ij1234I, D01235R,Dij1235R,D01235I,Dij1235I, D01245R,Dij1245R
     &   ,D01245I,Dij1245I, D01345R,Dij1345R,D01345I,Dij1345I, D02345
     &   R,Dij2345R,D02345I,Dij2345I 
       COMMON/NoAbe2gaPenLineCVALFunctions/ C0123R,Cij123R,C0123I,Cij1
     &   23I, C0124R,Cij124R,C0124I,Cij124I, C0125R,Cij125R,C0125I,Ci
     &   j125I, C0134R,Cij134R,C0134I,Cij134I, C0135R,Cij135R,C0135I,
     &   Cij135I, C0145R,Cij145R,C0145I,Cij145I, C0234R,Cij234R,C0234
     &   I,Cij234I, C0235R,Cij235R,C0235I,Cij235I, C0245R,Cij245R,C02
     &   45I,Cij245I, C0345R,Cij345R,C0345I,Cij345I 
       COMMON/NoAbe2gaPenLineBVALFunctions/ B012R,B012I, B013R,B013I, 
     &   B014R,B014I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,B0
     &   25I, B034R,B034I, B035R,B035I, B045R,B045I
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       FRe(211) = 4*(8*EijR(2,2)+2*(EijR(2,1)-EijR(4,1))+7*EijR(4,2)-1
     &   5*EijR(9,2)-4*(EijR(4,4)-EijR(10,4))+6*(EijR(2,3)+EijR(15,3)
     &   )-12*(EijR(10,3)-EijR(15,4)+EijR(21,4)))
       FIm(211) = 4*(8*EijI(2,2)+2*(EijI(2,1)-EijI(4,1))+7*EijI(4,2)-1
     &   5*EijI(9,2)-4*(EijI(4,4)-EijI(10,4))+6*(EijI(2,3)+EijI(15,3)
     &   )-12*(EijI(10,3)-EijI(15,4)+EijI(21,4)))
       F6(211)=DCMPLX(FRe(211),FIm(211))
       FRe(212) = -8*(EijR(4,2)+EijR(8,2)-EijR(9,2)-EijR(10,2)-3*(EijR
     &   (4,3)-EijR(9,3)+EijR(10,3)-EijR(16,3))+2*(EijR(4,4)+EijR(9,4
     &   )-EijR(10,4)-EijR(16,4))+6*(EijR(15,3)-EijR(15,4)-EijR(20,3)
     &   +EijR(21,4)-EijR(28,4)+EijR(34,4)))
       FIm(212) = -8*(EijI(4,2)+EijI(8,2)-EijI(9,2)-EijI(10,2)-3*(EijI
     &   (4,3)-EijI(9,3)+EijI(10,3)-EijI(16,3))+2*(EijI(4,4)+EijI(9,4
     &   )-EijI(10,4)-EijI(16,4))+6*(EijI(15,3)-EijI(15,4)-EijI(20,3)
     &   +EijI(21,4)-EijI(28,4)+EijI(34,4)))
       F6(212)=DCMPLX(FRe(212),FIm(212))
       FRe(213) = 8*(EijR(4,3)+EijR(10,3)-2*(EijR(4,4)-EijR(10,4)+EijR
     &   (15,3))+6*(EijR(15,4)-EijR(21,4)))
       FIm(213) = 8*(EijI(4,3)+EijI(10,3)-2*(EijI(4,4)-EijI(10,4)+EijI
     &   (15,3))+6*(EijI(15,4)-EijI(21,4)))
       F6(213)=DCMPLX(FRe(213),FIm(213))
       FRe(214) = -8*(EijR(2,1)-EijR(4,1)+4*EijR(4,2)-5*EijR(4,3)+EijR
     &   (5,2)-EijR(7,2)-7*EijR(9,2)-9*EijR(10,3)+3*(EijR(2,2)+EijR(8
     &   ,3)+EijR(14,3))+2*(EijR(2,3)+EijR(4,4)+EijR(8,4)-EijR(10,4)-
     &   EijR(14,4))+12*EijR(15,3)-6*(EijR(15,4)+EijR(18,3)-EijR(21,4
     &   )+EijR(27,4)-EijR(32,4)))
       FIm(214) = -8*(EijI(2,1)-EijI(4,1)+4*EijI(4,2)-5*EijI(4,3)+EijI
     &   (5,2)-EijI(7,2)-7*EijI(9,2)-9*EijI(10,3)+3*(EijI(2,2)+EijI(8
     &   ,3)+EijI(14,3))+2*(EijI(2,3)+EijI(4,4)+EijI(8,4)-EijI(10,4)-
     &   EijI(14,4))+12*EijI(15,3)-6*(EijI(15,4)+EijI(18,3)-EijI(21,4
     &   )+EijI(27,4)-EijI(32,4)))
       F6(214)=DCMPLX(FRe(214),FIm(214))
       FRe(215) = 8*(EijR(4,3)-EijR(9,3)+EijR(10,3)+6*EijR(15,4)-EijR(
     &   16,3)-2*(EijR(4,4)+EijR(9,4)-EijR(10,4)+EijR(15,3)-EijR(16,4
     &   )-EijR(20,3)+3*(EijR(21,4)-EijR(28,4)+EijR(34,4))))
       FIm(215) = 8*(EijI(4,3)-EijI(9,3)+EijI(10,3)+6*EijI(15,4)-EijI(
     &   16,3)-2*(EijI(4,4)+EijI(9,4)-EijI(10,4)+EijI(15,3)-EijI(16,4
     &   )-EijI(20,3)+3*(EijI(21,4)-EijI(28,4)+EijI(34,4))))
       F6(215)=DCMPLX(FRe(215),FIm(215))
       FRe(216) = -8*(EijR(2,2)+EijR(2,3)+EijR(4,2)-EijR(4,3)+2*(EijR(
     &   4,4)-EijR(9,2)-EijR(10,4))-3*(EijR(10,3)-EijR(15,3))-6*(EijR
     &   (15,4)-EijR(21,4)))
       FIm(216) = -8*(EijI(2,2)+EijI(2,3)+EijI(4,2)-EijI(4,3)+2*(EijI(
     &   4,4)-EijI(9,2)-EijI(10,4))-3*(EijI(10,3)-EijI(15,3))-6*(EijI
     &   (15,4)-EijI(21,4)))
       F6(216)=DCMPLX(FRe(216),FIm(216))
       FRe(217) = 8*(EijR(2,2)+EijR(2,3)-EijR(4,2)+3*EijR(4,3)+EijR(10
     &   ,3)-5*EijR(15,3)-4*(EijR(9,3)+EijR(16,3))-2*(EijR(4,4)+EijR(
     &   8,2)+EijR(9,4)-EijR(10,2)-EijR(10,4)-EijR(16,4))+8*EijR(20,3
     &   )+6*(EijR(15,4)-EijR(21,4)+EijR(28,4)-EijR(34,4)))
       FIm(217) = 8*(EijI(2,2)+EijI(2,3)-EijI(4,2)+3*EijI(4,3)+EijI(10
     &   ,3)-5*EijI(15,3)-4*(EijI(9,3)+EijI(16,3))-2*(EijI(4,4)+EijI(
     &   8,2)+EijI(9,4)-EijI(10,2)-EijI(10,4)-EijI(16,4))+8*EijI(20,3
     &   )+6*(EijI(15,4)-EijI(21,4)+EijI(28,4)-EijI(34,4)))
       F6(217)=DCMPLX(FRe(217),FIm(217))
       FRe(218) = 8*(EijR(4,3)-EijR(12,3)+EijR(13,3)-EijR(15,3)-2*(Eij
     &   R(4,4)+EijR(16,3)-EijR(20,3)+EijR(20,4)+EijR(21,4)+EijR(22,4
     &   )-2*(EijR(15,4)+EijR(16,4)+EijR(28,4)+EijR(31,4)-2*EijR(34,4
     &   ))))
       FIm(218) = 8*(EijI(4,3)-EijI(12,3)+EijI(13,3)-EijI(15,3)-2*(Eij
     &   I(4,4)+EijI(16,3)-EijI(20,3)+EijI(20,4)+EijI(21,4)+EijI(22,4
     &   )-2*(EijI(15,4)+EijI(16,4)+EijI(28,4)+EijI(31,4)-2*EijI(34,4
     &   ))))
       F6(218)=DCMPLX(FRe(218),FIm(218))
       FRe(219) = 4*(EijR(3,1)-EijR(4,1)+7*(EijR(4,2)-EijR(10,2))+6*(E
     &   ijR(8,2)-EijR(9,2)+EijR(9,3)-EijR(10,3))+10*(EijR(15,3)-EijR
     &   (20,3))-4*(EijR(4,3)+EijR(4,4)-EijR(16,3)-EijR(16,4)+EijR(21
     &   ,4)-EijR(28,4))+8*(EijR(15,4)-EijR(34,4)))
       FIm(219) = 4*(EijI(3,1)-EijI(4,1)+7*(EijI(4,2)-EijI(10,2))+6*(E
     &   ijI(8,2)-EijI(9,2)+EijI(9,3)-EijI(10,3))+10*(EijI(15,3)-EijI
     &   (20,3))-4*(EijI(4,3)+EijI(4,4)-EijI(16,3)-EijI(16,4)+EijI(21
     &   ,4)-EijI(28,4))+8*(EijI(15,4)-EijI(34,4)))
       F6(219)=DCMPLX(FRe(219),FIm(219))
       FRe(220) = -8*(EijR(4,2)-EijR(4,3)+EijR(8,2)-EijR(9,2)+EijR(9,3
     &   )-EijR(10,2)-EijR(10,3)-4*EijR(15,4)+EijR(16,3)+2*(EijR(4,4)
     &   +EijR(15,3)-EijR(16,4)-EijR(20,3)+EijR(21,4)-EijR(28,4)+2*Ei
     &   jR(34,4)))
       FIm(220) = -8*(EijI(4,2)-EijI(4,3)+EijI(8,2)-EijI(9,2)+EijI(9,3
     &   )-EijI(10,2)-EijI(10,3)-4*EijI(15,4)+EijI(16,3)+2*(EijI(4,4)
     &   +EijI(15,3)-EijI(16,4)-EijI(20,3)+EijI(21,4)-EijI(28,4)+2*Ei
     &   jI(34,4)))
       F6(220)=DCMPLX(FRe(220),FIm(220))
       FRe(221) = -8*(EijR(4,2)-EijR(8,2)+EijR(9,2)-EijR(9,3)-3*(EijR(
     &   4,3)+EijR(10,2))+EijR(10,3)+4*(EijR(12,3)-EijR(13,3)-EijR(15
     &   ,4))+7*EijR(16,3)+2*(EijR(3,2)+EijR(4,4)+EijR(15,3)-3*EijR(2
     &   0,3)+EijR(20,4)+EijR(21,4)+EijR(22,4)-2*(EijR(16,4)+EijR(28,
     &   4)+EijR(31,4)-2*EijR(34,4))))
       FIm(221) = -8*(EijI(4,2)-EijI(8,2)+EijI(9,2)-EijI(9,3)-3*(EijI(
     &   4,3)+EijI(10,2))+EijI(10,3)+4*(EijI(12,3)-EijI(13,3)-EijI(15
     &   ,4))+7*EijI(16,3)+2*(EijI(3,2)+EijI(4,4)+EijI(15,3)-3*EijI(2
     &   0,3)+EijI(20,4)+EijI(21,4)+EijI(22,4)-2*(EijI(16,4)+EijI(28,
     &   4)+EijI(31,4)-2*EijI(34,4))))
       F6(221)=DCMPLX(FRe(221),FIm(221))
       FRe(222) = 4*(7*EijR(4,2)-4*EijR(4,4)+8*EijR(8,2)-11*EijR(9,2)+
     &   6*EijR(9,3)+2*(EijR(3,1)-EijR(4,1)-4*(EijR(10,3)-EijR(15,3)-
     &   EijR(15,4))-EijR(16,3)-2*(EijR(10,2)-EijR(16,4)+EijR(20,3)+E
     &   ijR(21,4)-EijR(28,4)+2*EijR(34,4))))
       FIm(222) = 4*(7*EijI(4,2)-4*EijI(4,4)+8*EijI(8,2)-11*EijI(9,2)+
     &   6*EijI(9,3)+2*(EijI(3,1)-EijI(4,1)-4*(EijI(10,3)-EijI(15,3)-
     &   EijI(15,4))-EijI(16,3)-2*(EijI(10,2)-EijI(16,4)+EijI(20,3)+E
     &   ijI(21,4)-EijI(28,4)+2*EijI(34,4))))
       F6(222)=DCMPLX(FRe(222),FIm(222))
       FRe(223) = 8*(EijR(4,3)+EijR(9,3)-EijR(10,3)-3*EijR(16,3)-2*(Ei
     &   jR(4,4)+EijR(12,3)-EijR(13,3)-EijR(20,3)+EijR(20,4)+EijR(21,
     &   4)+EijR(22,4))+4*(EijR(15,4)+EijR(16,4)+EijR(28,4)+EijR(31,4
     &   ))-8*EijR(34,4))
       FIm(223) = 8*(EijI(4,3)+EijI(9,3)-EijI(10,3)-3*EijI(16,3)-2*(Ei
     &   jI(4,4)+EijI(12,3)-EijI(13,3)-EijI(20,3)+EijI(20,4)+EijI(21,
     &   4)+EijI(22,4))+4*(EijI(15,4)+EijI(16,4)+EijI(28,4)+EijI(31,4
     &   ))-8*EijI(34,4))
       F6(223)=DCMPLX(FRe(223),FIm(223))
       FRe(224) = -8*(EijR(4,3)+EijR(9,3)-EijR(10,3)-EijR(16,3)+2*(Eij
     &   R(4,2)+EijR(4,4)+EijR(8,2)-EijR(9,2)-EijR(10,2)-EijR(16,4)+E
     &   ijR(21,4)-EijR(28,4))-4*(EijR(15,4)-EijR(34,4)))
       FIm(224) = -8*(EijI(4,3)+EijI(9,3)-EijI(10,3)-EijI(16,3)+2*(Eij
     &   I(4,2)+EijI(4,4)+EijI(8,2)-EijI(9,2)-EijI(10,2)-EijI(16,4)+E
     &   ijI(21,4)-EijI(28,4))-4*(EijI(15,4)-EijI(34,4)))
       F6(224)=DCMPLX(FRe(224),FIm(224))
       FRe(225) = -8*(EijR(2,1)+3*EijR(2,2)-EijR(4,1)+4*EijR(4,2)+EijR
     &   (5,2)-EijR(7,2)-7*EijR(9,2)-5*(EijR(4,3)-EijR(8,3)+EijR(10,3
     &   )-EijR(14,3))+2*(EijR(4,4)+EijR(17,4)+5*(EijR(15,3)-EijR(18,
     &   3))+EijR(19,4)+EijR(21,4)-2*(EijR(14,4)+EijR(15,4)+EijR(24,4
     &   )+EijR(27,4)-2*EijR(32,4))))
       FIm(225) = -8*(EijI(2,1)+3*EijI(2,2)-EijI(4,1)+4*EijI(4,2)+EijI
     &   (5,2)-EijI(7,2)-7*EijI(9,2)-5*(EijI(4,3)-EijI(8,3)+EijI(10,3
     &   )-EijI(14,3))+2*(EijI(4,4)+EijI(17,4)+5*(EijI(15,3)-EijI(18,
     &   3))+EijI(19,4)+EijI(21,4)-2*(EijI(14,4)+EijI(15,4)+EijI(24,4
     &   )+EijI(27,4)-2*EijI(32,4))))
       F6(225)=DCMPLX(FRe(225),FIm(225))
       FRe(226) = 8*(EijR(4,3)+3*EijR(10,3)-2*(EijR(4,4)-EijR(14,3)-Ei
     &   jR(14,4)+EijR(18,3)+EijR(21,4)-EijR(27,4)+2*(EijR(15,3)-EijR
     &   (15,4)+EijR(32,4))))
       FIm(226) = 8*(EijI(4,3)+3*EijI(10,3)-2*(EijI(4,4)-EijI(14,3)-Ei
     &   jI(14,4)+EijI(18,3)+EijI(21,4)-EijI(27,4)+2*(EijI(15,3)-EijI
     &   (15,4)+EijI(32,4))))
       F6(226)=DCMPLX(FRe(226),FIm(226))
       FRe(227) = 8*(EijR(4,3)-3*(EijR(9,3)-EijR(10,3))-4*(EijR(15,3)-
     &   EijR(15,4))-EijR(16,3)-2*(EijR(4,4)-EijR(14,3)-EijR(14,4)-Ei
     &   jR(16,4)-EijR(17,3)+EijR(18,3)+EijR(19,3)+EijR(21,4)+EijR(26
     &   ,4)-EijR(27,4)-EijR(28,4)+EijR(33,4)-2*(EijR(20,3)-EijR(32,4
     &   )-EijR(34,4)+EijR(35,4))))
       FIm(227) = 8*(EijI(4,3)-3*(EijI(9,3)-EijI(10,3))-4*(EijI(15,3)-
     &   EijI(15,4))-EijI(16,3)-2*(EijI(4,4)-EijI(14,3)-EijI(14,4)-Ei
     &   jI(16,4)-EijI(17,3)+EijI(18,3)+EijI(19,3)+EijI(21,4)+EijI(26
     &   ,4)-EijI(27,4)-EijI(28,4)+EijI(33,4)-2*(EijI(20,3)-EijI(32,4
     &   )-EijI(34,4)+EijI(35,4))))
       F6(227)=DCMPLX(FRe(227),FIm(227))
       FRe(228) = -4*(EijR(1,1)+EijR(1,2)+EijR(2,1)+8*EijR(4,2)-10*Eij
     &   R(4,3)-9*EijR(7,2)+6*(EijR(5,3)-EijR(7,3))+4*(EijR(4,4)+EijR
     &   (8,3))+7*(EijR(5,2)-EijR(9,2))-2*(EijR(4,1)-8*EijR(14,3)+4*E
     &   ijR(14,4)-7*EijR(15,3)+2*(EijR(10,3)-EijR(17,4)+5*EijR(18,3)
     &   -EijR(19,4)-EijR(21,4)+2*(EijR(15,4)+EijR(24,4)+EijR(27,4)-2
     &   *EijR(32,4)))))
       FIm(228) = -4*(EijI(1,1)+EijI(1,2)+EijI(2,1)+8*EijI(4,2)-10*Eij
     &   I(4,3)-9*EijI(7,2)+6*(EijI(5,3)-EijI(7,3))+4*(EijI(4,4)+EijI
     &   (8,3))+7*(EijI(5,2)-EijI(9,2))-2*(EijI(4,1)-8*EijI(14,3)+4*E
     &   ijI(14,4)-7*EijI(15,3)+2*(EijI(10,3)-EijI(17,4)+5*EijI(18,3)
     &   -EijI(19,4)-EijI(21,4)+2*(EijI(15,4)+EijI(24,4)+EijI(27,4)-2
     &   *EijI(32,4)))))
       F6(228)=DCMPLX(FRe(228),FIm(228))
       FRe(229) = 4*(-EijR(7,2)+EijR(9,2)+2*(EijR(4,3)-EijR(14,3)-EijR
     &   (15,3)+EijR(18,3))-4*(EijR(4,4)-EijR(14,4)+EijR(21,4)-EijR(2
     &   7,4))+8*(EijR(15,4)-EijR(32,4)))
       FIm(229) = 4*(-EijI(7,2)+EijI(9,2)+2*(EijI(4,3)-EijI(14,3)-EijI
     &   (15,3)+EijI(18,3))-4*(EijI(4,4)-EijI(14,4)+EijI(21,4)-EijI(2
     &   7,4))+8*(EijI(15,4)-EijI(32,4)))
       F6(229)=DCMPLX(FRe(229),FIm(229))
       FRe(230) = 4*(-4*EijR(4,4)+EijR(6,2)-EijR(7,2)-EijR(8,2)+EijR(9
     &   ,2)+2*(EijR(4,3)-EijR(14,3)-EijR(15,3)+4*EijR(15,4)-EijR(16,
     &   3)-EijR(17,3)+EijR(18,3)+EijR(19,3)+EijR(20,3)+2*(EijR(14,4)
     &   +EijR(16,4)-EijR(21,4)-EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR
     &   (33,4)-2*(EijR(32,4)+EijR(34,4)-EijR(35,4)))))
       FIm(230) = 4*(-4*EijI(4,4)+EijI(6,2)-EijI(7,2)-EijI(8,2)+EijI(9
     &   ,2)+2*(EijI(4,3)-EijI(14,3)-EijI(15,3)+4*EijI(15,4)-EijI(16,
     &   3)-EijI(17,3)+EijI(18,3)+EijI(19,3)+EijI(20,3)+2*(EijI(14,4)
     &   +EijI(16,4)-EijI(21,4)-EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI
     &   (33,4)-2*(EijI(32,4)+EijI(34,4)-EijI(35,4)))))
       F6(230)=DCMPLX(FRe(230),FIm(230))
       FRe(231) = -8*(EijR(2,2)+EijR(4,2)+3*(EijR(2,3)-EijR(4,3))+2*(E
     &   ijR(4,4)+EijR(8,4)-EijR(9,2)-EijR(10,4)-EijR(14,4))-9*(EijR(
     &   10,3)-EijR(15,3))-6*(EijR(15,4)-EijR(21,4)+EijR(27,4)-EijR(3
     &   2,4)))
       FIm(231) = -8*(EijI(2,2)+EijI(4,2)+3*(EijI(2,3)-EijI(4,3))+2*(E
     &   ijI(4,4)+EijI(8,4)-EijI(9,2)-EijI(10,4)-EijI(14,4))-9*(EijI(
     &   10,3)-EijI(15,3))-6*(EijI(15,4)-EijI(21,4)+EijI(27,4)-EijI(3
     &   2,4)))
       F6(231)=DCMPLX(FRe(231),FIm(231))
       FRe(232) = -4*(EijR(2,2)+4*EijR(4,4)+EijR(5,2)-EijR(7,2)-6*(Eij
     &   R(4,3)-EijR(8,3))-3*EijR(9,2)+2*(EijR(4,2)-3*(EijR(10,3)-Eij
     &   R(14,3))+2*(EijR(17,4)+3*(EijR(15,3)-EijR(18,3))+EijR(19,4)+
     &   EijR(21,4)-2*(EijR(14,4)+EijR(15,4)+EijR(24,4)+EijR(27,4)-2*
     &   EijR(32,4)))))
       FIm(232) = -4*(EijI(2,2)+4*EijI(4,4)+EijI(5,2)-EijI(7,2)-6*(Eij
     &   I(4,3)-EijI(8,3))-3*EijI(9,2)+2*(EijI(4,2)-3*(EijI(10,3)-Eij
     &   I(14,3))+2*(EijI(17,4)+3*(EijI(15,3)-EijI(18,3))+EijI(19,4)+
     &   EijI(21,4)-2*(EijI(14,4)+EijI(15,4)+EijI(24,4)+EijI(27,4)-2*
     &   EijI(32,4)))))
       F6(232)=DCMPLX(FRe(232),FIm(232))
       FRe(233) = 4*(2*(EijR(2,1)-EijR(4,1))+7*EijR(4,2)-11*EijR(9,2)+
     &   6*(EijR(8,3)+EijR(14,3))+4*(EijR(2,2)-EijR(4,4)+EijR(5,2)-Ei
     &   jR(7,2)+EijR(14,4)-3*EijR(18,3)-EijR(21,4)+EijR(27,4)+2*(Eij
     &   R(15,4)-EijR(32,4))))
       FIm(233) = 4*(2*(EijI(2,1)-EijI(4,1))+7*EijI(4,2)-11*EijI(9,2)+
     &   6*(EijI(8,3)+EijI(14,3))+4*(EijI(2,2)-EijI(4,4)+EijI(5,2)-Ei
     &   jI(7,2)+EijI(14,4)-3*EijI(18,3)-EijI(21,4)+EijI(27,4)+2*(Eij
     &   I(15,4)-EijI(32,4))))
       F6(233)=DCMPLX(FRe(233),FIm(233))
       FRe(234) = -8*(EijR(4,2)+EijR(8,2)-EijR(9,2)-EijR(10,2)-3*(EijR
     &   (4,3)-EijR(9,3)+EijR(10,3)-EijR(16,3))+6*(EijR(15,3)-EijR(20
     &   ,3))+2*(EijR(4,4)-EijR(14,4)-EijR(16,4)+EijR(21,4)+EijR(26,4
     &   )-EijR(27,4)-EijR(28,4)+EijR(33,4))-4*(EijR(15,4)-EijR(32,4)
     &   -EijR(34,4)+EijR(35,4)))
       FIm(234) = -8*(EijI(4,2)+EijI(8,2)-EijI(9,2)-EijI(10,2)-3*(EijI
     &   (4,3)-EijI(9,3)+EijI(10,3)-EijI(16,3))+6*(EijI(15,3)-EijI(20
     &   ,3))+2*(EijI(4,4)-EijI(14,4)-EijI(16,4)+EijI(21,4)+EijI(26,4
     &   )-EijI(27,4)-EijI(28,4)+EijI(33,4))-4*(EijI(15,4)-EijI(32,4)
     &   -EijI(34,4)+EijI(35,4)))
       F6(234)=DCMPLX(FRe(234),FIm(234))
       FRe(235) = -4*(EijR(2,2)+EijR(5,2)-EijR(7,2)-6*(EijR(4,3)-EijR(
     &   8,3))+4*(EijR(4,4)+EijR(8,4))-3*EijR(9,2)+2*(EijR(4,2)-3*(Ei
     &   jR(10,3)-EijR(14,3))-2*(EijR(10,4)+EijR(14,4)-3*(EijR(15,3)-
     &   EijR(15,4)-EijR(18,3)+EijR(21,4)-EijR(27,4)+EijR(32,4)))))
       FIm(235) = -4*(EijI(2,2)+EijI(5,2)-EijI(7,2)-6*(EijI(4,3)-EijI(
     &   8,3))+4*(EijI(4,4)+EijI(8,4))-3*EijI(9,2)+2*(EijI(4,2)-3*(Ei
     &   jI(10,3)-EijI(14,3))-2*(EijI(10,4)+EijI(14,4)-3*(EijI(15,3)-
     &   EijI(15,4)-EijI(18,3)+EijI(21,4)-EijI(27,4)+EijI(32,4)))))
       F6(235)=DCMPLX(FRe(235),FIm(235))
       FRe(236) = 4*(3*EijR(2,2)-4*EijR(4,4)-5*(EijR(5,2)-EijR(7,2))-E
     &   ijR(9,2)-2*(EijR(4,2)-EijR(4,3)+EijR(8,3)-EijR(10,3)+EijR(14
     &   ,3)-2*(EijR(14,4)-EijR(15,3)+EijR(18,3)-EijR(21,4)+EijR(27,4
     &   )+2*(EijR(15,4)-EijR(32,4)))))
       FIm(236) = 4*(3*EijI(2,2)-4*EijI(4,4)-5*(EijI(5,2)-EijI(7,2))-E
     &   ijI(9,2)-2*(EijI(4,2)-EijI(4,3)+EijI(8,3)-EijI(10,3)+EijI(14
     &   ,3)-2*(EijI(14,4)-EijI(15,3)+EijI(18,3)-EijI(21,4)+EijI(27,4
     &   )+2*(EijI(15,4)-EijI(32,4)))))
       F6(236)=DCMPLX(FRe(236),FIm(236))
       FRe(237) = -4*(EijR(2,2)-6*EijR(4,3)-3*EijR(5,2)+4*(EijR(4,4)+E
     &   ijR(6,2))-EijR(7,2)+EijR(9,2)+2*(EijR(4,2)-EijR(8,3)+EijR(10
     &   ,3)+3*EijR(14,3)-2*(EijR(10,2)+EijR(14,4)-EijR(15,3)+EijR(16
     &   ,4)+EijR(18,3)-EijR(21,4)-EijR(26,4)+EijR(27,4)+EijR(28,4)-E
     &   ijR(33,4)+2*(EijR(15,4)-EijR(16,3)-EijR(17,3)+EijR(19,3)+Eij
     &   R(20,3)-EijR(32,4)-EijR(34,4)+EijR(35,4)))))
       FIm(237) = -4*(EijI(2,2)-6*EijI(4,3)-3*EijI(5,2)+4*(EijI(4,4)+E
     &   ijI(6,2))-EijI(7,2)+EijI(9,2)+2*(EijI(4,2)-EijI(8,3)+EijI(10
     &   ,3)+3*EijI(14,3)-2*(EijI(10,2)+EijI(14,4)-EijI(15,3)+EijI(16
     &   ,4)+EijI(18,3)-EijI(21,4)-EijI(26,4)+EijI(27,4)+EijI(28,4)-E
     &   ijI(33,4)+2*(EijI(15,4)-EijI(16,3)-EijI(17,3)+EijI(19,3)+Eij
     &   I(20,3)-EijI(32,4)-EijI(34,4)+EijI(35,4)))))
       F6(237)=DCMPLX(FRe(237),FIm(237))
       FRe(238) = -8*(-EijR(4,3)+EijR(10,3)+2*(EijR(3,1)-EijR(4,1)+Eij
     &   R(4,4)+EijR(8,2)-EijR(10,2)+EijR(16,3)-EijR(16,4)-EijR(20,3)
     &   +EijR(21,4)-EijR(28,4)-2*(EijR(15,4)-EijR(34,4))))
       FIm(238) = -8*(-EijI(4,3)+EijI(10,3)+2*(EijI(3,1)-EijI(4,1)+Eij
     &   I(4,4)+EijI(8,2)-EijI(10,2)+EijI(16,3)-EijI(16,4)-EijI(20,3)
     &   +EijI(21,4)-EijI(28,4)-2*(EijI(15,4)-EijI(34,4))))
       F6(238)=DCMPLX(FRe(238),FIm(238))
       FRe(239) = -8*(EijR(2,1)-EijR(2,2)-EijR(4,1)-5*EijR(4,3)+EijR(5
     &   ,2)-EijR(7,2)-EijR(8,3)-EijR(10,3)-3*(EijR(9,2)-EijR(14,3))+
     &   6*(EijR(15,3)+EijR(16,3))-8*EijR(20,3)+2*(EijR(4,4)+EijR(9,3
     &   )-EijR(14,4)-EijR(16,4)-EijR(18,3)+EijR(21,4)+EijR(26,4)-Eij
     &   R(27,4)-EijR(28,4)+EijR(33,4))+4*(EijR(4,2)+EijR(8,2)-EijR(1
     &   0,2)-EijR(15,4)+EijR(17,3)-EijR(19,3)+EijR(32,4)+EijR(34,4)-
     &   EijR(35,4)))
       FIm(239) = -8*(EijI(2,1)-EijI(2,2)-EijI(4,1)-5*EijI(4,3)+EijI(5
     &   ,2)-EijI(7,2)-EijI(8,3)-EijI(10,3)-3*(EijI(9,2)-EijI(14,3))+
     &   6*(EijI(15,3)+EijI(16,3))-8*EijI(20,3)+2*(EijI(4,4)+EijI(9,3
     &   )-EijI(14,4)-EijI(16,4)-EijI(18,3)+EijI(21,4)+EijI(26,4)-Eij
     &   I(27,4)-EijI(28,4)+EijI(33,4))+4*(EijI(4,2)+EijI(8,2)-EijI(1
     &   0,2)-EijI(15,4)+EijI(17,3)-EijI(19,3)+EijI(32,4)+EijI(34,4)-
     &   EijI(35,4)))
       F6(239)=DCMPLX(FRe(239),FIm(239))
       FRe(240) = 8*(EijR(3,1)-EijR(4,1)+EijR(4,3)-EijR(15,3)-EijR(16,
     &   3)+EijR(20,3)+2*(EijR(4,2)-EijR(4,4)+EijR(8,2)-EijR(9,2)-Eij
     &   R(10,2)+EijR(16,4)-EijR(21,4)+EijR(28,4))+4*(EijR(15,4)-EijR
     &   (34,4)))
       FIm(240) = 8*(EijI(3,1)-EijI(4,1)+EijI(4,3)-EijI(15,3)-EijI(16,
     &   3)+EijI(20,3)+2*(EijI(4,2)-EijI(4,4)+EijI(8,2)-EijI(9,2)-Eij
     &   I(10,2)+EijI(16,4)-EijI(21,4)+EijI(28,4))+4*(EijI(15,4)-EijI
     &   (34,4)))
       F6(240)=DCMPLX(FRe(240),FIm(240))
       FRe(241) = -8*(EijR(3,1)-EijR(4,1)+EijR(6,2)-EijR(7,2)+4*(EijR(
     &   4,2)-EijR(10,2)-EijR(15,4))-5*(EijR(4,3)-EijR(16,3))+3*(EijR
     &   (8,2)-EijR(9,2)+EijR(14,3)+EijR(17,3)-EijR(18,3)-EijR(19,3))
     &   +7*(EijR(15,3)-EijR(20,3))+2*(EijR(4,4)+EijR(9,3)-EijR(10,3)
     &   -EijR(14,4)-EijR(16,4)+EijR(21,4)+EijR(26,4)-EijR(27,4)-EijR
     &   (28,4)+EijR(33,4)+2*(EijR(32,4)+EijR(34,4)-EijR(35,4))))
       FIm(241) = -8*(EijI(3,1)-EijI(4,1)+EijI(6,2)-EijI(7,2)+4*(EijI(
     &   4,2)-EijI(10,2)-EijI(15,4))-5*(EijI(4,3)-EijI(16,3))+3*(EijI
     &   (8,2)-EijI(9,2)+EijI(14,3)+EijI(17,3)-EijI(18,3)-EijI(19,3))
     &   +7*(EijI(15,3)-EijI(20,3))+2*(EijI(4,4)+EijI(9,3)-EijI(10,3)
     &   -EijI(14,4)-EijI(16,4)+EijI(21,4)+EijI(26,4)-EijI(27,4)-EijI
     &   (28,4)+EijI(33,4)+2*(EijI(32,4)+EijI(34,4)-EijI(35,4))))
       F6(241)=DCMPLX(FRe(241),FIm(241))
       FRe(242) = -8*(EijR(3,1)-EijR(4,1)+EijR(4,3)+EijR(15,3)-EijR(16
     &   ,3)-EijR(20,3)+2*(EijR(4,4)-EijR(15,4)-EijR(16,4)+EijR(34,4)
     &   ))
       FIm(242) = -8*(EijI(3,1)-EijI(4,1)+EijI(4,3)+EijI(15,3)-EijI(16
     &   ,3)-EijI(20,3)+2*(EijI(4,4)-EijI(15,4)-EijI(16,4)+EijI(34,4)
     &   ))
       F6(242)=DCMPLX(FRe(242),FIm(242))
       FRe(243) = 4*(-8*EijR(4,3)+2*(EijR(8,2)-EijR(9,2))+3*(EijR(3,1)
     &   -EijR(4,1)-EijR(4,2)+EijR(10,2))-4*(EijR(4,4)-EijR(15,3)-Eij
     &   R(15,4)-2*EijR(16,3)-EijR(16,4)+EijR(20,3)+EijR(34,4)))
       FIm(243) = 4*(-8*EijI(4,3)+2*(EijI(8,2)-EijI(9,2))+3*(EijI(3,1)
     &   -EijI(4,1)-EijI(4,2)+EijI(10,2))-4*(EijI(4,4)-EijI(15,3)-Eij
     &   I(15,4)-2*EijI(16,3)-EijI(16,4)+EijI(20,3)+EijI(34,4)))
       F6(243)=DCMPLX(FRe(243),FIm(243))
       FRe(244) = -8*(EijR(3,1)-EijR(4,1)+EijR(6,2)-EijR(7,2)+3*(EijR(
     &   8,2)-EijR(9,2))+4*(EijR(4,2)-EijR(10,2)-EijR(14,4))-5*(EijR(
     &   4,3)-EijR(14,3)-EijR(15,3)-EijR(16,3)-EijR(17,3)+EijR(18,3)+
     &   EijR(19,3)+EijR(20,3))+2*(EijR(4,4)-EijR(15,4)-EijR(16,4)+Ei
     &   jR(19,4)+EijR(23,4)-EijR(24,4)-EijR(25,4)+EijR(34,4)+2*(EijR
     &   (32,4)+EijR(33,4)-EijR(35,4))))
       FIm(244) = -8*(EijI(3,1)-EijI(4,1)+EijI(6,2)-EijI(7,2)+3*(EijI(
     &   8,2)-EijI(9,2))+4*(EijI(4,2)-EijI(10,2)-EijI(14,4))-5*(EijI(
     &   4,3)-EijI(14,3)-EijI(15,3)-EijI(16,3)-EijI(17,3)+EijI(18,3)+
     &   EijI(19,3)+EijI(20,3))+2*(EijI(4,4)-EijI(15,4)-EijI(16,4)+Ei
     &   jI(19,4)+EijI(23,4)-EijI(24,4)-EijI(25,4)+EijI(34,4)+2*(EijI
     &   (32,4)+EijI(33,4)-EijI(35,4))))
       F6(244)=DCMPLX(FRe(244),FIm(244))
       FRe(245) = 4*(-4*EijR(4,4)+EijR(6,2)-EijR(7,2)-EijR(8,2)+EijR(9
     &   ,2)+2*(EijR(4,3)-EijR(14,3)+EijR(15,3)+4*EijR(16,4)+EijR(17,
     &   3)-EijR(18,3)-3*(EijR(16,3)-EijR(19,3))-EijR(20,3)-2*(EijR(1
     &   1,3)-EijR(13,3)-EijR(14,4)-EijR(15,4)+EijR(22,4)+EijR(29,4)-
     &   EijR(30,4)-EijR(31,4)+EijR(32,4)+2*(EijR(33,4)+EijR(34,4)-Ei
     &   jR(35,4)))))
       FIm(245) = 4*(-4*EijI(4,4)+EijI(6,2)-EijI(7,2)-EijI(8,2)+EijI(9
     &   ,2)+2*(EijI(4,3)-EijI(14,3)+EijI(15,3)+4*EijI(16,4)+EijI(17,
     &   3)-EijI(18,3)-3*(EijI(16,3)-EijI(19,3))-EijI(20,3)-2*(EijI(1
     &   1,3)-EijI(13,3)-EijI(14,4)-EijI(15,4)+EijI(22,4)+EijI(29,4)-
     &   EijI(30,4)-EijI(31,4)+EijI(32,4)+2*(EijI(33,4)+EijI(34,4)-Ei
     &   jI(35,4)))))
       F6(245)=DCMPLX(FRe(245),FIm(245))
       FRe(246) = -4*(4*(EijR(4,2)+EijR(4,4))+7*(EijR(6,2)-EijR(7,2))-
     &   3*(EijR(8,2)-EijR(9,2))+2*(EijR(4,3)-EijR(14,3)+EijR(15,3)-E
     &   ijR(16,3)+EijR(17,3)-EijR(18,3)+EijR(19,3)-EijR(20,3)-2*(Eij
     &   R(10,2)+EijR(14,4)+EijR(15,4)+EijR(16,4)-EijR(32,4)-EijR(33,
     &   4)-EijR(34,4)+EijR(35,4))))
       FIm(246) = -4*(4*(EijI(4,2)+EijI(4,4))+7*(EijI(6,2)-EijI(7,2))-
     &   3*(EijI(8,2)-EijI(9,2))+2*(EijI(4,3)-EijI(14,3)+EijI(15,3)-E
     &   ijI(16,3)+EijI(17,3)-EijI(18,3)+EijI(19,3)-EijI(20,3)-2*(Eij
     &   I(10,2)+EijI(14,4)+EijI(15,4)+EijI(16,4)-EijI(32,4)-EijI(33,
     &   4)-EijI(34,4)+EijI(35,4))))
       F6(246)=DCMPLX(FRe(246),FIm(246))
       FRe(247) = 8*(EijR(4,3)+EijR(13,3)-3*(EijR(12,3)+EijR(15,3))-2*
     &   (EijR(4,4)-EijR(11,3)-EijR(14,3)-EijR(14,4)-EijR(15,4)+EijR(
     &   16,3)-3*EijR(20,3)+EijR(22,4)+EijR(29,4)-EijR(30,4)-EijR(31,
     &   4)+EijR(32,4)-2*(EijR(16,4)-EijR(19,3)-EijR(33,4)-EijR(34,4)
     &   +EijR(35,4))))
       FIm(247) = 8*(EijI(4,3)+EijI(13,3)-3*(EijI(12,3)+EijI(15,3))-2*
     &   (EijI(4,4)-EijI(11,3)-EijI(14,3)-EijI(14,4)-EijI(15,4)+EijI(
     &   16,3)-3*EijI(20,3)+EijI(22,4)+EijI(29,4)-EijI(30,4)-EijI(31,
     &   4)+EijI(32,4)-2*(EijI(16,4)-EijI(19,3)-EijI(33,4)-EijI(34,4)
     &   +EijI(35,4))))
       F6(247)=DCMPLX(FRe(247),FIm(247))
       FRe(248) = 4*(EijR(6,2)-EijR(7,2)-4*(EijR(4,4)-EijR(6,3)+EijR(7
     &   ,3))-EijR(8,2)+EijR(9,2)+2*(EijR(4,3)+EijR(14,3)+4*EijR(14,4
     &   )-EijR(16,3)-EijR(19,3)-3*(EijR(15,3)+EijR(17,3)-EijR(18,3)-
     &   EijR(20,3))+2*(EijR(15,4)+EijR(16,4)-EijR(19,4)-EijR(23,4)+E
     &   ijR(24,4)+EijR(25,4)-EijR(34,4)-2*(EijR(32,4)+EijR(33,4)-Eij
     &   R(35,4)))))
       FIm(248) = 4*(EijI(6,2)-EijI(7,2)-4*(EijI(4,4)-EijI(6,3)+EijI(7
     &   ,3))-EijI(8,2)+EijI(9,2)+2*(EijI(4,3)+EijI(14,3)+4*EijI(14,4
     &   )-EijI(16,3)-EijI(19,3)-3*(EijI(15,3)+EijI(17,3)-EijI(18,3)-
     &   EijI(20,3))+2*(EijI(15,4)+EijI(16,4)-EijI(19,4)-EijI(23,4)+E
     &   ijI(24,4)+EijI(25,4)-EijI(34,4)-2*(EijI(32,4)+EijI(33,4)-Eij
     &   I(35,4)))))
       F6(248)=DCMPLX(FRe(248),FIm(248))
       FRe(249) = 4*(EijR(3,1)-EijR(4,1)-4*(EijR(4,3)+EijR(4,4)-EijR(8
     &   ,2)+EijR(9,2))+7*(EijR(4,2)-EijR(10,2))+2*(EijR(6,2)-EijR(7,
     &   2)+3*(EijR(17,3)-EijR(18,3))+5*(EijR(14,3)-EijR(19,3))+2*(Ei
     &   jR(14,4)+EijR(15,4)+EijR(16,3)+EijR(16,4)-EijR(32,4)-EijR(33
     &   ,4)-EijR(34,4)+EijR(35,4))))
       FIm(249) = 4*(EijI(3,1)-EijI(4,1)-4*(EijI(4,3)+EijI(4,4)-EijI(8
     &   ,2)+EijI(9,2))+7*(EijI(4,2)-EijI(10,2))+2*(EijI(6,2)-EijI(7,
     &   2)+3*(EijI(17,3)-EijI(18,3))+5*(EijI(14,3)-EijI(19,3))+2*(Ei
     &   jI(14,4)+EijI(15,4)+EijI(16,3)+EijI(16,4)-EijI(32,4)-EijI(33
     &   ,4)-EijI(34,4)+EijI(35,4))))
       F6(249)=DCMPLX(FRe(249),FIm(249))
       FRe(250) = 8*(EijR(3,1)-EijR(4,1)+EijR(4,3)-EijR(16,3)-3*(EijR(
     &   15,3)-EijR(20,3))+2*(EijR(4,2)-EijR(4,4)+EijR(6,2)-EijR(7,2)
     &   -EijR(10,2)+EijR(14,3)+EijR(14,4)+EijR(15,4)+EijR(16,4)-EijR
     &   (19,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4)))
       FIm(250) = 8*(EijI(3,1)-EijI(4,1)+EijI(4,3)-EijI(16,3)-3*(EijI(
     &   15,3)-EijI(20,3))+2*(EijI(4,2)-EijI(4,4)+EijI(6,2)-EijI(7,2)
     &   -EijI(10,2)+EijI(14,3)+EijI(14,4)+EijI(15,4)+EijI(16,4)-EijI
     &   (19,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4)))
       F6(250)=DCMPLX(FRe(250),FIm(250))
       FRe(251) = 8*(EijR(3,1)-EijR(4,1)+EijR(6,2)-EijR(7,2)-EijR(8,2)
     &   +EijR(9,2)-EijR(14,3)-EijR(17,3)+EijR(18,3)+EijR(19,3)+3*(Ei
     &   jR(4,3)-EijR(15,3)-EijR(16,3)+EijR(20,3))-2*(EijR(4,2)+EijR(
     &   4,4)-EijR(10,2)-EijR(14,4)-EijR(15,4)-EijR(16,4)+EijR(32,4)+
     &   EijR(33,4)+EijR(34,4)-EijR(35,4)))
       FIm(251) = 8*(EijI(3,1)-EijI(4,1)+EijI(6,2)-EijI(7,2)-EijI(8,2)
     &   +EijI(9,2)-EijI(14,3)-EijI(17,3)+EijI(18,3)+EijI(19,3)+3*(Ei
     &   jI(4,3)-EijI(15,3)-EijI(16,3)+EijI(20,3))-2*(EijI(4,2)+EijI(
     &   4,4)-EijI(10,2)-EijI(14,4)-EijI(15,4)-EijI(16,4)+EijI(32,4)+
     &   EijI(33,4)+EijI(34,4)-EijI(35,4)))
       F6(251)=DCMPLX(FRe(251),FIm(251))
       FRe(252) = -4*(EijR(1,1)+EijR(1,2)+EijR(2,1)-10*EijR(4,3)+4*Eij
     &   R(4,4)-EijR(5,2)+8*(EijR(4,2)+EijR(6,2)+EijR(6,3))-9*EijR(7,
     &   2)-6*EijR(7,3)+EijR(9,2)-2*(EijR(4,1)+EijR(5,3)-8*EijR(14,3)
     &   +4*(EijR(10,2)+EijR(14,4))-EijR(15,3)+2*(EijR(15,4)-3*EijR(1
     &   6,3)+EijR(16,4)-EijR(17,3)+5*EijR(19,3)-EijR(19,4)+EijR(20,3
     &   )-EijR(23,4)+EijR(24,4)+EijR(25,4)-EijR(34,4)-2*(EijR(32,4)+
     &   EijR(33,4)-EijR(35,4)))))
       FIm(252) = -4*(EijI(1,1)+EijI(1,2)+EijI(2,1)-10*EijI(4,3)+4*Eij
     &   I(4,4)-EijI(5,2)+8*(EijI(4,2)+EijI(6,2)+EijI(6,3))-9*EijI(7,
     &   2)-6*EijI(7,3)+EijI(9,2)-2*(EijI(4,1)+EijI(5,3)-8*EijI(14,3)
     &   +4*(EijI(10,2)+EijI(14,4))-EijI(15,3)+2*(EijI(15,4)-3*EijI(1
     &   6,3)+EijI(16,4)-EijI(17,3)+5*EijI(19,3)-EijI(19,4)+EijI(20,3
     &   )-EijI(23,4)+EijI(24,4)+EijI(25,4)-EijI(34,4)-2*(EijI(32,4)+
     &   EijI(33,4)-EijI(35,4)))))
       F6(252)=DCMPLX(FRe(252),FIm(252))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
