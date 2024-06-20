       subroutine PenBoxFFhl2(F2)
       IMPLICIT NONE
       Real*8 P(48),FRe(44:86),FIm(44:86)
       Complex*16 F2(44:86)
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
       FRe(44) = -4*(EijR(3,2)-EijR(8,2)+EijR(9,2)-EijR(10,2)-2*(EijR(
     &   3,3)+EijR(9,3)-EijR(10,3)-EijR(13,3)-2*(EijR(12,3)-EijR(20,3
     &   ))))
       FIm(44) = -4*(EijI(3,2)-EijI(8,2)+EijI(9,2)-EijI(10,2)-2*(EijI(
     &   3,3)+EijI(9,3)-EijI(10,3)-EijI(13,3)-2*(EijI(12,3)-EijI(20,3
     &   ))))
       F2(44)=DCMPLX(FRe(44),FIm(44))
       FRe(45) = 4*(EijR(3,2)-EijR(8,2)+2*(EijR(3,3)-EijR(12,3)-EijR(1
     &   3,3)+EijR(20,3)))
       FIm(45) = 4*(EijI(3,2)-EijI(8,2)+2*(EijI(3,3)-EijI(12,3)-EijI(1
     &   3,3)+EijI(20,3)))
       F2(45)=DCMPLX(FRe(45),FIm(45))
       FRe(46) = 4*(EijR(2,1)-EijR(3,1)-EijR(3,2)+EijR(5,2)-EijR(6,2)+
     &   EijR(8,2)+2*(EijR(3,3)-EijR(9,2)+EijR(10,2)-EijR(11,3)-EijR(
     &   12,3)-EijR(13,3)+EijR(17,3)-EijR(18,3)+EijR(19,3)+EijR(20,3)
     &   ))
       FIm(46) = 4*(EijI(2,1)-EijI(3,1)-EijI(3,2)+EijI(5,2)-EijI(6,2)+
     &   EijI(8,2)+2*(EijI(3,3)-EijI(9,2)+EijI(10,2)-EijI(11,3)-EijI(
     &   12,3)-EijI(13,3)+EijI(17,3)-EijI(18,3)+EijI(19,3)+EijI(20,3)
     &   ))
       F2(46)=DCMPLX(FRe(46),FIm(46))
       FRe(47) = -2*(2*EijR(11,2)+EijR(23,3))
       FIm(47) = -2*(2*EijI(11,2)+EijI(23,3))
       F2(47)=DCMPLX(FRe(47),FIm(47))
       FRe(48) = EijR(11,2)
       FIm(48) = EijI(11,2)
       F2(48)=DCMPLX(FRe(48),FIm(48))
       P(46) = 2*s12-s45
       FRe(49) = -Dij2345R(1,1)+Dij2345R(2,1)-s12*EijR(2,2)+s15*(EijR(
     &   2,1)-EijR(3,1)+EijR(5,2)-EijR(6,2))+s45*(EijR(9,2)-EijR(10,2
     &   ))+2*(EijR(11,2)+EijR(22,3)-EijR(23,3))+EijR(3,2)*P(44)+EijR
     &   (8,2)*P(46)
       FIm(49) = -Dij2345I(1,1)+Dij2345I(2,1)-s12*EijI(2,2)+s15*(EijI(
     &   2,1)-EijI(3,1)+EijI(5,2)-EijI(6,2))+s45*(EijI(9,2)-EijI(10,2
     &   ))+2*(EijI(11,2)+EijI(22,3)-EijI(23,3))+EijI(3,2)*P(44)+EijI
     &   (8,2)*P(46)
       F2(49)=DCMPLX(FRe(49),FIm(49))
       FRe(50) = -EijR(3,3)+EijR(11,3)+EijR(12,3)-EijR(17,3)
       FIm(50) = -EijI(3,3)+EijI(11,3)+EijI(12,3)-EijI(17,3)
       F2(50)=DCMPLX(FRe(50),FIm(50))
       FRe(51) = EijR(3,3)-EijR(11,3)-EijR(12,3)-EijR(13,3)+EijR(17,3)
     &   -EijR(18,3)+EijR(19,3)+EijR(20,3)
       FIm(51) = EijI(3,3)-EijI(11,3)-EijI(12,3)-EijI(13,3)+EijI(17,3)
     &   -EijI(18,3)+EijI(19,3)+EijI(20,3)
       F2(51)=DCMPLX(FRe(51),FIm(51))
       FRe(52) = 2*(EijR(2,2)+EijR(3,2)-EijR(3,3)+EijR(8,3)-EijR(9,3)+
     &   EijR(11,3)-2*(EijR(8,2)-EijR(12,3)+EijR(17,3)))
       FIm(52) = 2*(EijI(2,2)+EijI(3,2)-EijI(3,3)+EijI(8,3)-EijI(9,3)+
     &   EijI(11,3)-2*(EijI(8,2)-EijI(12,3)+EijI(17,3)))
       F2(52)=DCMPLX(FRe(52),FIm(52))
       FRe(53) = 2*(EijR(3,2)-EijR(3,3)+EijR(5,2)+EijR(5,3)-EijR(6,2)-
     &   EijR(6,3)-EijR(8,2)+EijR(12,3)+2*(EijR(11,3)-EijR(17,3)))
       FIm(53) = 2*(EijI(3,2)-EijI(3,3)+EijI(5,2)+EijI(5,3)-EijI(6,2)-
     &   EijI(6,3)-EijI(8,2)+EijI(12,3)+2*(EijI(11,3)-EijI(17,3)))
       F2(53)=DCMPLX(FRe(53),FIm(53))
       FRe(54) = -2*(EijR(3,3)+EijR(9,3)-2*EijR(12,3))
       FIm(54) = -2*(EijI(3,3)+EijI(9,3)-2*EijI(12,3))
       F2(54)=DCMPLX(FRe(54),FIm(54))
       FRe(55) = -2*(EijR(3,2)+EijR(3,3)-EijR(8,2)-EijR(12,3))
       FIm(55) = -2*(EijI(3,2)+EijI(3,3)-EijI(8,2)-EijI(12,3))
       F2(55)=DCMPLX(FRe(55),FIm(55))
       FRe(56) = EijR(3,2)+EijR(3,3)-EijR(8,2)+EijR(9,2)-EijR(10,2)-Ei
     &   jR(12,3)-EijR(13,3)+EijR(20,3)
       FIm(56) = EijI(3,2)+EijI(3,3)-EijI(8,2)+EijI(9,2)-EijI(10,2)-Ei
     &   jI(12,3)-EijI(13,3)+EijI(20,3)
       F2(56)=DCMPLX(FRe(56),FIm(56))
       FRe(57) = EijR(2,1)-EijR(3,1)+EijR(3,3)+EijR(5,2)-EijR(6,2)-Eij
     &   R(11,3)-EijR(12,3)+EijR(17,3)
       FIm(57) = EijI(2,1)-EijI(3,1)+EijI(3,3)+EijI(5,2)-EijI(6,2)-Eij
     &   I(11,3)-EijI(12,3)+EijI(17,3)
       F2(57)=DCMPLX(FRe(57),FIm(57))
       FRe(58) = -EijR(2,1)+EijR(3,1)-EijR(5,2)+EijR(6,2)
       FIm(58) = -EijI(2,1)+EijI(3,1)-EijI(5,2)+EijI(6,2)
       F2(58)=DCMPLX(FRe(58),FIm(58))
       P(47) = s12-2*s45
       FRe(59) = Dij2345R(2,1)-Dij2345R(3,1)+s45*EijR(4,2)+s15*(-EijR(
     &   3,1)+EijR(4,1)-EijR(6,2)+EijR(7,2))+s12*(EijR(8,2)-EijR(9,2)
     &   )-2*(EijR(11,2)+EijR(23,3)-EijR(24,3))+EijR(3,2)*P(44)+EijR(
     &   10,2)*P(47)
       FIm(59) = Dij2345I(2,1)-Dij2345I(3,1)+s45*EijI(4,2)+s15*(-EijI(
     &   3,1)+EijI(4,1)-EijI(6,2)+EijI(7,2))+s12*(EijI(8,2)-EijI(9,2)
     &   )-2*(EijI(11,2)+EijI(23,3)-EijI(24,3))+EijI(3,2)*P(44)+EijI(
     &   10,2)*P(47)
       F2(59)=DCMPLX(FRe(59),FIm(59))
       FRe(60) = -8*(EijR(3,3)+EijR(6,3)-EijR(7,3)-EijR(13,3)-2*(EijR(
     &   11,3)-EijR(19,3)))
       FIm(60) = -8*(EijI(3,3)+EijI(6,3)-EijI(7,3)-EijI(13,3)-2*(EijI(
     &   11,3)-EijI(19,3)))
       F2(60)=DCMPLX(FRe(60),FIm(60))
       FRe(61) = EijR(3,2)+EijR(3,3)-EijR(6,2)+EijR(7,2)-EijR(10,2)-Ei
     &   jR(11,3)-EijR(14,3)+EijR(16,3)-2*(EijR(13,3)-EijR(19,3))
       FIm(61) = EijI(3,2)+EijI(3,3)-EijI(6,2)+EijI(7,2)-EijI(10,2)-Ei
     &   jI(11,3)-EijI(14,3)+EijI(16,3)-2*(EijI(13,3)-EijI(19,3))
       F2(61)=DCMPLX(FRe(61),FIm(61))
       FRe(62) = EijR(3,2)+EijR(3,3)-EijR(6,2)-EijR(11,3)-EijR(13,3)+E
     &   ijR(19,3)
       FIm(62) = EijI(3,2)+EijI(3,3)-EijI(6,2)-EijI(11,3)-EijI(13,3)+E
     &   ijI(19,3)
       F2(62)=DCMPLX(FRe(62),FIm(62))
       FRe(63) = EijR(2,1)-EijR(3,1)+EijR(3,3)+EijR(5,2)-EijR(6,2)-Eij
     &   R(9,2)+EijR(10,2)-EijR(11,3)-EijR(12,3)-EijR(13,3)+EijR(17,3
     &   )-EijR(18,3)+EijR(19,3)+EijR(20,3)
       FIm(63) = EijI(2,1)-EijI(3,1)+EijI(3,3)+EijI(5,2)-EijI(6,2)-Eij
     &   I(9,2)+EijI(10,2)-EijI(11,3)-EijI(12,3)-EijI(13,3)+EijI(17,3
     &   )-EijI(18,3)+EijI(19,3)+EijI(20,3)
       F2(63)=DCMPLX(FRe(63),FIm(63))
       FRe(64) = 2*(EijR(3,2)-EijR(3,3)-EijR(6,2)-EijR(6,3)+EijR(7,2)+
     &   EijR(7,3)-EijR(10,2)+EijR(13,3)+2*(EijR(11,3)-EijR(19,3)))
       FIm(64) = 2*(EijI(3,2)-EijI(3,3)-EijI(6,2)-EijI(6,3)+EijI(7,2)+
     &   EijI(7,3)-EijI(10,2)+EijI(13,3)+2*(EijI(11,3)-EijI(19,3)))
       F2(64)=DCMPLX(FRe(64),FIm(64))
       FRe(65) = EijR(3,3)-EijR(11,3)-EijR(13,3)+EijR(19,3)
       FIm(65) = EijI(3,3)-EijI(11,3)-EijI(13,3)+EijI(19,3)
       F2(65)=DCMPLX(FRe(65),FIm(65))
       FRe(66) = -EijR(3,3)+EijR(11,3)+EijR(14,3)-EijR(16,3)+2*(EijR(1
     &   3,3)-EijR(19,3))
       FIm(66) = -EijI(3,3)+EijI(11,3)+EijI(14,3)-EijI(16,3)+2*(EijI(1
     &   3,3)-EijI(19,3))
       F2(66)=DCMPLX(FRe(66),FIm(66))
       FRe(67) = EijR(3,2)-EijR(3,3)-EijR(8,2)+EijR(9,2)-EijR(10,2)+Ei
     &   jR(11,3)+EijR(12,3)+EijR(13,3)-EijR(17,3)+EijR(18,3)-EijR(19
     &   ,3)-EijR(20,3)
       FIm(67) = EijI(3,2)-EijI(3,3)-EijI(8,2)+EijI(9,2)-EijI(10,2)+Ei
     &   jI(11,3)+EijI(12,3)+EijI(13,3)-EijI(17,3)+EijI(18,3)-EijI(19
     &   ,3)-EijI(20,3)
       F2(67)=DCMPLX(FRe(67),FIm(67))
       FRe(68) = -2*(EijR(3,2)+EijR(3,3)-EijR(10,2)-EijR(13,3))
       FIm(68) = -2*(EijI(3,2)+EijI(3,3)-EijI(10,2)-EijI(13,3))
       F2(68)=DCMPLX(FRe(68),FIm(68))
       FRe(69) = EijR(3,1)-EijR(3,3)-EijR(4,1)+EijR(6,2)-EijR(7,2)+Eij
     &   R(11,3)+EijR(13,3)-EijR(19,3)
       FIm(69) = EijI(3,1)-EijI(3,3)-EijI(4,1)+EijI(6,2)-EijI(7,2)+Eij
     &   I(11,3)+EijI(13,3)-EijI(19,3)
       F2(69)=DCMPLX(FRe(69),FIm(69))
       FRe(70) = EijR(3,2)+EijR(3,3)+EijR(4,2)-2*(EijR(10,2)+EijR(13,3
     &   ))+EijR(16,3)
       FIm(70) = EijI(3,2)+EijI(3,3)+EijI(4,2)-2*(EijI(10,2)+EijI(13,3
     &   ))+EijI(16,3)
       F2(70)=DCMPLX(FRe(70),FIm(70))
       FRe(71) = EijR(3,1)-EijR(4,1)+EijR(6,2)-EijR(7,2)
       FIm(71) = EijI(3,1)-EijI(4,1)+EijI(6,2)-EijI(7,2)
       F2(71)=DCMPLX(FRe(71),FIm(71))
       FRe(72) = -2*(EijR(3,1)-EijR(4,1)+EijR(8,2)-2*EijR(9,2)+EijR(10
     &   ,2))
       FIm(72) = -2*(EijI(3,1)-EijI(4,1)+EijI(8,2)-2*EijI(9,2)+EijI(10
     &   ,2))
       F2(72)=DCMPLX(FRe(72),FIm(72))
       FRe(73) = EijR(11,2)+EijR(21,3)-EijR(23,3)
       FIm(73) = EijI(11,2)+EijI(21,3)-EijI(23,3)
       F2(73)=DCMPLX(FRe(73),FIm(73))
       FRe(74) = EijR(22,3)-EijR(23,3)
       FIm(74) = EijI(22,3)-EijI(23,3)
       F2(74)=DCMPLX(FRe(74),FIm(74))
       FRe(75) = EijR(23,3)-EijR(24,3)
       FIm(75) = EijI(23,3)-EijI(24,3)
       F2(75)=DCMPLX(FRe(75),FIm(75))
       FRe(76) = EijR(23,3)
       FIm(76) = EijI(23,3)
       F2(76)=DCMPLX(FRe(76),FIm(76))
       FRe(77) = EijR(2,1)-EijR(3,1)+EijR(5,2)-EijR(6,2)
       FIm(77) = EijI(2,1)-EijI(3,1)+EijI(5,2)-EijI(6,2)
       F2(77)=DCMPLX(FRe(77),FIm(77))
       FRe(78) = -EijR(11,2)+EijR(21,3)-EijR(23,3)
       FIm(78) = -EijI(11,2)+EijI(21,3)-EijI(23,3)
       F2(78)=DCMPLX(FRe(78),FIm(78))
       FRe(79) = -2*(EijR(2,1)-EijR(3,1)+EijR(5,2)+EijR(7,2)+EijR(8,2)
     &   -2*(EijR(6,2)+EijR(9,2))+EijR(10,2))
       FIm(79) = -2*(EijI(2,1)-EijI(3,1)+EijI(5,2)+EijI(7,2)+EijI(8,2)
     &   -2*(EijI(6,2)+EijI(9,2))+EijI(10,2))
       F2(79)=DCMPLX(FRe(79),FIm(79))
       FRe(80) = s15*(EijR(2,1)+EijR(5,2))+s45*EijR(9,2)+(EijR(3,1)+Ei
     &   jR(8,2))*P(15)
       FIm(80) = s15*(EijI(2,1)+EijI(5,2))+s45*EijI(9,2)+(EijI(3,1)+Ei
     &   jI(8,2))*P(15)
       F2(80)=DCMPLX(FRe(80),FIm(80))
       P(48) = -p3sq-s15+s23+s34-s45
       FRe(81) = (EijR(2,1)+EijR(5,2))*P(22)+EijR(9,2)*P(24)+EijR(10,2
     &   )*P(34)+(EijR(3,1)+EijR(6,2))*P(48)
       FIm(81) = (EijI(2,1)+EijI(5,2))*P(22)+EijI(9,2)*P(24)+EijI(10,2
     &   )*P(34)+(EijI(3,1)+EijI(6,2))*P(48)
       F2(81)=DCMPLX(FRe(81),FIm(81))
       FRe(82) = EijR(2,1)-EijR(3,1)+EijR(5,2)-EijR(8,2)
       FIm(82) = EijI(2,1)-EijI(3,1)+EijI(5,2)-EijI(8,2)
       F2(82)=DCMPLX(FRe(82),FIm(82))
       FRe(83) = EijR(3,1)-EijR(4,1)+EijR(6,2)-EijR(10,2)
       FIm(83) = EijI(3,1)-EijI(4,1)+EijI(6,2)-EijI(10,2)
       F2(83)=DCMPLX(FRe(83),FIm(83))
       FRe(84) = EijR(11,2)+EijR(22,3)-EijR(23,3)
       FIm(84) = EijI(11,2)+EijI(22,3)-EijI(23,3)
       F2(84)=DCMPLX(FRe(84),FIm(84))
       FRe(85) = EijR(2,1)-EijR(3,1)+EijR(5,2)-EijR(6,2)-EijR(9,2)+Eij
     &   R(10,2)
       FIm(85) = EijI(2,1)-EijI(3,1)+EijI(5,2)-EijI(6,2)-EijI(9,2)+Eij
     &   I(10,2)
       F2(85)=DCMPLX(FRe(85),FIm(85))
       FRe(86) = EijR(11,2)+EijR(23,3)-EijR(24,3)
       FIm(86) = EijI(11,2)+EijI(23,3)-EijI(24,3)
       F2(86)=DCMPLX(FRe(86),FIm(86))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
