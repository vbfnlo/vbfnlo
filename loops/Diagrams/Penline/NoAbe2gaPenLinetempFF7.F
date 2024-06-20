       subroutine NoAbe2gaPenLineFFhl7(F7)
       IMPLICIT NONE
       Real*8 P(767),FRe(253:294),FIm(253:294)
       Complex*16 F7(253:294)
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
       FRe(253) = -4*(4*(EijR(3,1)-EijR(4,1)+EijR(4,4)+EijR(6,2))+EijR
     &   (7,2)-EijR(9,2)-2*(EijR(4,3)-EijR(14,3)+EijR(15,3)-EijR(18,3
     &   )+2*(EijR(10,2)+EijR(14,4)+EijR(15,4)-EijR(16,3)+EijR(16,4)+
     &   EijR(19,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4))))
       FIm(253) = -4*(4*(EijI(3,1)-EijI(4,1)+EijI(4,4)+EijI(6,2))+EijI
     &   (7,2)-EijI(9,2)-2*(EijI(4,3)-EijI(14,3)+EijI(15,3)-EijI(18,3
     &   )+2*(EijI(10,2)+EijI(14,4)+EijI(15,4)-EijI(16,3)+EijI(16,4)+
     &   EijI(19,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4))))
       F7(253)=DCMPLX(FRe(253),FIm(253))
       FRe(254) = 4*(EE0R+EijR(1,1)+4*(EijR(2,1)+EijR(4,3)-EijR(4,4)+E
     &   ijR(5,2))-3*(EijR(4,1)-EijR(4,2)+EijR(7,2))+6*(EijR(2,2)+Eij
     &   R(8,3))-2*(5*EijR(9,2)+EijR(10,3)+EijR(15,3)+3*EijR(18,3)-2*
     &   (EijR(14,4)-EijR(21,4)+EijR(27,4))-4*(EijR(15,4)-EijR(32,4))
     &   ))
       FIm(254) = 4*(EE0I+EijI(1,1)+4*(EijI(2,1)+EijI(4,3)-EijI(4,4)+E
     &   ijI(5,2))-3*(EijI(4,1)-EijI(4,2)+EijI(7,2))+6*(EijI(2,2)+Eij
     &   I(8,3))-2*(5*EijI(9,2)+EijI(10,3)+EijI(15,3)+3*EijI(18,3)-2*
     &   (EijI(14,4)-EijI(21,4)+EijI(27,4))-4*(EijI(15,4)-EijI(32,4))
     &   ))
       F7(254)=DCMPLX(FRe(254),FIm(254))
       FRe(255) = 4*(EE0R+EijR(1,1)+4*(EijR(3,1)+EijR(4,3)-EijR(4,4)+E
     &   ijR(6,2))-3*(EijR(4,1)-EijR(4,2)+EijR(7,2))+2*(-EijR(10,2)+3
     &   *(EijR(8,2)-EijR(16,3)+EijR(17,3))-4*(EijR(9,2)+EijR(18,3))+
     &   EijR(19,3)-EijR(20,3)+2*(EijR(14,4)+EijR(15,3)+EijR(15,4)+Ei
     &   jR(16,4)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4))))
       FIm(255) = 4*(EE0I+EijI(1,1)+4*(EijI(3,1)+EijI(4,3)-EijI(4,4)+E
     &   ijI(6,2))-3*(EijI(4,1)-EijI(4,2)+EijI(7,2))+2*(-EijI(10,2)+3
     &   *(EijI(8,2)-EijI(16,3)+EijI(17,3))-4*(EijI(9,2)+EijI(18,3))+
     &   EijI(19,3)-EijI(20,3)+2*(EijI(14,4)+EijI(15,3)+EijI(15,4)+Ei
     &   jI(16,4)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4))))
       F7(255)=DCMPLX(FRe(255),FIm(255))
       FRe(256) = 4*(EE0R+EijR(1,1)-3*(EijR(4,1)-EijR(4,2))+10*EijR(5,
     &   2)-9*EijR(7,2)+6*(EijR(5,3)-EijR(7,3))+2*(EijR(14,3)-EijR(18
     &   ,3))+4*(EijR(2,1)+EijR(4,3)-EijR(4,4)-EijR(9,2)-EijR(15,3)+E
     &   ijR(15,4)-EijR(19,4)+EijR(24,4))+8*(EijR(14,4)-EijR(32,4)))
       FIm(256) = 4*(EE0I+EijI(1,1)-3*(EijI(4,1)-EijI(4,2))+10*EijI(5,
     &   2)-9*EijI(7,2)+6*(EijI(5,3)-EijI(7,3))+2*(EijI(14,3)-EijI(18
     &   ,3))+4*(EijI(2,1)+EijI(4,3)-EijI(4,4)-EijI(9,2)-EijI(15,3)+E
     &   ijI(15,4)-EijI(19,4)+EijI(24,4))+8*(EijI(14,4)-EijI(32,4)))
       F7(256)=DCMPLX(FRe(256),FIm(256))
       FRe(257) = -4*(EijR(4,1)-7*EijR(4,2)+2*(EijR(7,2)-5*EijR(14,3)+
     &   3*EijR(18,3)+2*(EijR(4,3)+EijR(4,4)+EijR(9,2)-EijR(14,4)-Eij
     &   R(15,4)+EijR(32,4))))
       FIm(257) = -4*(EijI(4,1)-7*EijI(4,2)+2*(EijI(7,2)-5*EijI(14,3)+
     &   3*EijI(18,3)+2*(EijI(4,3)+EijI(4,4)+EijI(9,2)-EijI(14,4)-Eij
     &   I(15,4)+EijI(32,4))))
       F7(257)=DCMPLX(FRe(257),FIm(257))
       FRe(258) = -4*(EE0R+EijR(1,1)+5*EijR(4,1)-EijR(4,2)-2*(EijR(2,1
     &   )+EijR(5,2))+EijR(7,2)+6*EijR(9,2)+4*(EijR(4,4)-EijR(14,3)-E
     &   ijR(14,4)-EijR(15,4)+EijR(18,3)+EijR(32,4)))
       FIm(258) = -4*(EE0I+EijI(1,1)+5*EijI(4,1)-EijI(4,2)-2*(EijI(2,1
     &   )+EijI(5,2))+EijI(7,2)+6*EijI(9,2)+4*(EijI(4,4)-EijI(14,3)-E
     &   ijI(14,4)-EijI(15,4)+EijI(18,3)+EijI(32,4)))
       F7(258)=DCMPLX(FRe(258),FIm(258))
       FRe(259) = -4*(EijR(4,1)-7*EijR(4,2)+6*(EijR(9,2)+EijR(10,3))-1
     &   0*EijR(15,3)-8*EijR(15,4)+4*(EijR(4,3)+EijR(4,4)+EijR(21,4))
     &   )
       FIm(259) = -4*(EijI(4,1)-7*EijI(4,2)+6*(EijI(9,2)+EijI(10,3))-1
     &   0*EijI(15,3)-8*EijI(15,4)+4*(EijI(4,3)+EijI(4,4)+EijI(21,4))
     &   )
       F7(259)=DCMPLX(FRe(259),FIm(259))
       FRe(260) = -4*(EijR(3,1)-7*EijR(4,2)-8*EijR(15,3)-2*(EijR(4,1)-
     &   EijR(10,2)+EijR(16,3))+6*EijR(20,3)+4*(EijR(4,3)+EijR(4,4)-E
     &   ijR(15,4)-EijR(16,4)+EijR(34,4)))
       FIm(260) = -4*(EijI(3,1)-7*EijI(4,2)-8*EijI(15,3)-2*(EijI(4,1)-
     &   EijI(10,2)+EijI(16,3))+6*EijI(20,3)+4*(EijI(4,3)+EijI(4,4)-E
     &   ijI(15,4)-EijI(16,4)+EijI(34,4)))
       F7(260)=DCMPLX(FRe(260),FIm(260))
       FRe(261) = -4*(3*(EijR(4,1)+EijR(4,2))+2*(4*EijR(4,3)+EijR(9,2)
     &   +2*(EijR(4,4)-EijR(15,3)-EijR(15,4))))
       FIm(261) = -4*(3*(EijI(4,1)+EijI(4,2))+2*(4*EijI(4,3)+EijI(9,2)
     &   +2*(EijI(4,4)-EijI(15,3)-EijI(15,4))))
       F7(261)=DCMPLX(FRe(261),FIm(261))
       P(653) = p5sq+s15-s23-s45+2*P(187)
       P(654) = s23+2*s34+s45
       P(655) = p2sq+4*p4sq+p5sq+s12+s15-2*P(654)
       FRe(262) = 4*(Dij1245R(2,2)+Dij1245R(3,2)-2*Dij1245R(6,2)-Dij13
     &   45R(2,1)+Dij1345R(3,1)-Dij2345R(2,1)+Dij2345R(3,1)+4*(EijR(2
     &   3,3)-EijR(24,3)-EijR(39,4)-EijR(43,4)+EijR(44,4)+EijR(45,4))
     &   +EijR(3,2)*P(286)+(-EijR(8,2)+EijR(9,2))*P(292)+EijR(4,2)*P(
     &   653)-EijR(10,2)*P(655))
       FIm(262) = 4*(Dij1245I(2,2)+Dij1245I(3,2)-2*Dij1245I(6,2)-Dij13
     &   45I(2,1)+Dij1345I(3,1)-Dij2345I(2,1)+Dij2345I(3,1)+4*(EijI(2
     &   3,3)-EijI(24,3)-EijI(39,4)-EijI(43,4)+EijI(44,4)+EijI(45,4))
     &   +EijI(3,2)*P(286)+(-EijI(8,2)+EijI(9,2))*P(292)+EijI(4,2)*P(
     &   653)-EijI(10,2)*P(655))
       F7(262)=DCMPLX(FRe(262),FIm(262))
       P(656) = p2sq-p5sq+s12-s15+s23+s45-2*P(187)
       FRe(263) = 4*(Dij1245R(2,2)+Dij1245R(3,2)-2*Dij1245R(6,2)-Dij13
     &   45R(1,1)+Dij1345R(3,1)-Dij2345R(1,1)+Dij2345R(3,1)+4*(EijR(2
     &   2,3)-EijR(24,3)-EijR(37,4)-EijR(39,4))+8*EijR(44,4)+(EijR(8,
     &   2)-EijR(10,2))*P(286)-EijR(2,2)*P(292)+EijR(4,2)*P(653)+EijR
     &   (9,2)*P(656))
       FIm(263) = 4*(Dij1245I(2,2)+Dij1245I(3,2)-2*Dij1245I(6,2)-Dij13
     &   45I(1,1)+Dij1345I(3,1)-Dij2345I(1,1)+Dij2345I(3,1)+4*(EijI(2
     &   2,3)-EijI(24,3)-EijI(37,4)-EijI(39,4))+8*EijI(44,4)+(EijI(8,
     &   2)-EijI(10,2))*P(286)-EijI(2,2)*P(292)+EijI(4,2)*P(653)+EijI
     &   (9,2)*P(656))
       F7(263)=DCMPLX(FRe(263),FIm(263))
       P(657) = -p5sq-s15+s23+s45-2*P(187)
       FRe(264) = 4*(D02345R+Dij1245R(3,2)+Dij1245R(4,2)-Dij1245R(5,2)
     &   -Dij1245R(6,2)-Dij1345R(1,1)+Dij1345R(3,1)+Dij2345R(3,1)+2*(
     &   EijR(11,2)+2*(EijR(21,3)-EijR(24,3)-EijR(39,4)-EijR(40,4)+Ei
     &   jR(42,4)+EijR(44,4)))+EijR(6,2)*P(286)+(-EijR(5,2)+EijR(9,2)
     &   )*P(292)+EijR(10,2)*P(533)+EijR(4,2)*P(653)+EijR(7,2)*P(657)
     &   )
       FIm(264) = 4*(D02345I+Dij1245I(3,2)+Dij1245I(4,2)-Dij1245I(5,2)
     &   -Dij1245I(6,2)-Dij1345I(1,1)+Dij1345I(3,1)+Dij2345I(3,1)+2*(
     &   EijI(11,2)+2*(EijI(21,3)-EijI(24,3)-EijI(39,4)-EijI(40,4)+Ei
     &   jI(42,4)+EijI(44,4)))+EijI(6,2)*P(286)+(-EijI(5,2)+EijI(9,2)
     &   )*P(292)+EijI(10,2)*P(533)+EijI(4,2)*P(653)+EijI(7,2)*P(657)
     &   )
       F7(264)=DCMPLX(FRe(264),FIm(264))
       FRe(265) = 2*(D01345R+D02345R-Dij1245R(2,1)+Dij1245R(3,1)+2*(Di
     &   j1245R(3,2)-Dij1245R(6,2)+Dij1345R(3,1)+Dij2345R(3,1)+6*(Eij
     &   R(11,2)+EijR(22,3))-10*EijR(24,3)-4*(EijR(39,4)-EijR(44,4))-
     &   EijR(10,2)*P(286))+(EijR(2,1)+2*EijR(9,2))*P(292)+EijR(3,1)*
     &   P(533)+(EijR(4,1)+2*EijR(4,2))*P(653))
       FIm(265) = 2*(D01345I+D02345I-Dij1245I(2,1)+Dij1245I(3,1)+2*(Di
     &   j1245I(3,2)-Dij1245I(6,2)+Dij1345I(3,1)+Dij2345I(3,1)+6*(Eij
     &   I(11,2)+EijI(22,3))-10*EijI(24,3)-4*(EijI(39,4)-EijI(44,4))-
     &   EijI(10,2)*P(286))+(EijI(2,1)+2*EijI(9,2))*P(292)+EijI(3,1)*
     &   P(533)+(EijI(4,1)+2*EijI(4,2))*P(653))
       F7(265)=DCMPLX(FRe(265),FIm(265))
       P(658) = p4sq+p5sq-s12+s34-3*s45
       P(659) = p4sq-2*s12+s34-3*s45
       P(660) = 3*p4sq+s15-s23+4*P(19)+2*P(49)
       P(661) = p2sq+s15-7*s34-4*s45+6*P(199)+2*P(314)
       FRe(266) = EijR(3,2)*P(594)+2*(-D02345R-Dij1234R(3,1)+Dij1245R(
     &   2,1)+Dij1245R(2,2)-Dij1245R(3,1)+Dij1245R(3,2)-2*Dij1245R(6,
     &   2)+Dij1345R(2,1)-Dij1345R(3,1)-s12*EijR(2,1)+10*EijR(11,2)+4
     &   *(EijR(22,3)+EijR(24,3))-8*(EijR(23,3)+EijR(39,4)+EijR(43,4)
     &   -EijR(44,4)-EijR(45,4))+(-EijR(6,2)+EijR(7,2))*P(21)+(EijR(8
     &   ,2)-EijR(9,2))*P(270)+EijR(3,1)*P(658)-EijR(4,1)*P(659)-EijR
     &   (4,2)*P(660)+EijR(10,2)*P(661))
       FIm(266) = EijI(3,2)*P(594)+2*(-D02345I-Dij1234I(3,1)+Dij1245I(
     &   2,1)+Dij1245I(2,2)-Dij1245I(3,1)+Dij1245I(3,2)-2*Dij1245I(6,
     &   2)+Dij1345I(2,1)-Dij1345I(3,1)-s12*EijI(2,1)+10*EijI(11,2)+4
     &   *(EijI(22,3)+EijI(24,3))-8*(EijI(23,3)+EijI(39,4)+EijI(43,4)
     &   -EijI(44,4)-EijI(45,4))+(-EijI(6,2)+EijI(7,2))*P(21)+(EijI(8
     &   ,2)-EijI(9,2))*P(270)+EijI(3,1)*P(658)-EijI(4,1)*P(659)-EijI
     &   (4,2)*P(660)+EijI(10,2)*P(661))
       F7(266)=DCMPLX(FRe(266),FIm(266))
       P(662) = p4sq+p5sq+2*s12-s34-s45
       P(663) = p4sq+s12-s34
       P(664) = s15-s23+3*P(49)+4*P(663)
       P(665) = p2sq-5*s12-s15+s23-3*P(49)-4*P(187)
       FRe(267) = -2*(D02345R-Dij1245R(2,1)-Dij1245R(2,2)+Dij1245R(3,1
     &   )-Dij1245R(3,2)-Dij1345R(1,1)+Dij1345R(3,1)+Dij2345R(1,1)-Di
     &   j2345R(3,1)+s12*EijR(2,1)+2*(Dij1245R(6,2)-EijR(11,2))+4*(Ei
     &   jR(22,3)-EijR(24,3))+8*(EijR(37,4)+EijR(39,4))-16*EijR(44,4)
     &   -EijR(3,1)*P(230)-EijR(2,2)*P(266)-EijR(10,2)*P(333)-EijR(8,
     &   2)*P(568)-EijR(4,1)*P(662)+EijR(4,2)*P(664)+EijR(9,2)*P(665)
     &   )
       FIm(267) = -2*(D02345I-Dij1245I(2,1)-Dij1245I(2,2)+Dij1245I(3,1
     &   )-Dij1245I(3,2)-Dij1345I(1,1)+Dij1345I(3,1)+Dij2345I(1,1)-Di
     &   j2345I(3,1)+s12*EijI(2,1)+2*(Dij1245I(6,2)-EijI(11,2))+4*(Ei
     &   jI(22,3)-EijI(24,3))+8*(EijI(37,4)+EijI(39,4))-16*EijI(44,4)
     &   -EijI(3,1)*P(230)-EijI(2,2)*P(266)-EijI(10,2)*P(333)-EijI(8,
     &   2)*P(568)-EijI(4,1)*P(662)+EijI(4,2)*P(664)+EijI(9,2)*P(665)
     &   )
       F7(267)=DCMPLX(FRe(267),FIm(267))
       FRe(268) = 2*(Dij1245R(2,1)-Dij1245R(3,1)+Dij1245R(3,2)+Dij1245
     &   R(4,2)-Dij1245R(5,2)-Dij1245R(6,2)+Dij1345R(1,1)-Dij1345R(3,
     &   1)+Dij2345R(3,1)-s12*EijR(2,1)-12*EijR(22,3)+4*EijR(24,3)+8*
     &   (EijR(21,3)-EijR(39,4)-EijR(40,4)+EijR(42,4)+EijR(44,4))+Eij
     &   R(3,1)*P(230)+EijR(5,2)*P(266)+EijR(10,2)*P(333)+EijR(9,2)*P
     &   (337)+EijR(6,2)*P(568)+EijR(4,1)*P(662)+(-EijR(4,2)+EijR(7,2
     &   ))*P(664))
       FIm(268) = 2*(Dij1245I(2,1)-Dij1245I(3,1)+Dij1245I(3,2)+Dij1245
     &   I(4,2)-Dij1245I(5,2)-Dij1245I(6,2)+Dij1345I(1,1)-Dij1345I(3,
     &   1)+Dij2345I(3,1)-s12*EijI(2,1)-12*EijI(22,3)+4*EijI(24,3)+8*
     &   (EijI(21,3)-EijI(39,4)-EijI(40,4)+EijI(42,4)+EijI(44,4))+Eij
     &   I(3,1)*P(230)+EijI(5,2)*P(266)+EijI(10,2)*P(333)+EijI(9,2)*P
     &   (337)+EijI(6,2)*P(568)+EijI(4,1)*P(662)+(-EijI(4,2)+EijI(7,2
     &   ))*P(664))
       F7(268)=DCMPLX(FRe(268),FIm(268))
       P(666) = p4sq-2*s12+s15-s23-s34
       FRe(269) = 2*(D01345R+2*(D02345R-Dij1245R(2,1)+Dij1245R(3,1))+D
     &   ij1245R(3,2)-Dij1245R(6,2)-Dij1345R(3,1)+Dij2345R(3,1)-22*Ei
     &   jR(11,2)-4*(EijR(22,3)+EijR(24,3))-8*(EijR(39,4)-EijR(44,4))
     &   +EijR(2,1)*P(324)+EijR(10,2)*P(333)+EijR(9,2)*P(337)+EijR(3,
     &   1)*P(631)-EijR(4,2)*P(664)+EijR(4,1)*P(666))
       FIm(269) = 2*(D01345I+2*(D02345I-Dij1245I(2,1)+Dij1245I(3,1))+D
     &   ij1245I(3,2)-Dij1245I(6,2)-Dij1345I(3,1)+Dij2345I(3,1)-22*Ei
     &   jI(11,2)-4*(EijI(22,3)+EijI(24,3))-8*(EijI(39,4)-EijI(44,4))
     &   +EijI(2,1)*P(324)+EijI(10,2)*P(333)+EijI(9,2)*P(337)+EijI(3,
     &   1)*P(631)-EijI(4,2)*P(664)+EijI(4,1)*P(666))
       F7(269)=DCMPLX(FRe(269),FIm(269))
       P(667) = 4*p4sq+10*p5sq-8*s45
       P(668) = -6*p4sq-4*p5sq+10*s45
       P(669) = -4*p4sq-10*s15+8*s45
       P(670) = 4*p4sq-6*s15+8*P(49)
       P(671) = 2*p4sq+p5sq+4*P(67)
       P(672) = 2*p5sq-s12+s34-5*s45+3*P(448)
       P(673) = p4sq+s15-s45
       P(674) = p4sq+p5sq-s12-s15+s34-s45
       P(675) = p4sq-s12+s34-s45
       FRe(270) = EijR(4,1)*P(667)+EijR(4,2)*P(668)+EijR(5,2)*P(669)+E
     &   ijR(7,2)*P(670)+2*(Dij1234R(2,1)+Dij1234R(4,2)-2*(Dij1235R(2
     &   ,1)-Dij1235R(3,1)+Dij1235R(3,2)+Dij1235R(4,2)-Dij1235R(5,2)-
     &   Dij1235R(6,2))-Dij2345R(3,2)+Dij2345R(5,2)-8*(2*(EijR(11,2)+
     &   EijR(21,3))+EijR(22,3)-3*EijR(24,3)+EijR(39,4)+EijR(40,4)-Ei
     &   jR(42,4)-EijR(44,4))+(EijR(2,2)+EijR(8,3)-EijR(10,3))*P(1)-(
     &   4*(EijR(1,1)+EijR(1,2))-EijR(5,3)+EijR(7,3))*P(21)+(-EijR(8,
     &   2)+EijR(10,2)-EijR(16,3)-EijR(17,3)+EijR(19,3)+EijR(20,3))*P
     &   (233)+EijR(4,3)*P(385)-EijR(2,1)*P(671)+EijR(9,2)*P(672)-Eij
     &   R(14,3)*P(673)-EijR(15,3)*P(674)+EijR(18,3)*P(675))
       FIm(270) = EijI(4,1)*P(667)+EijI(4,2)*P(668)+EijI(5,2)*P(669)+E
     &   ijI(7,2)*P(670)+2*(Dij1234I(2,1)+Dij1234I(4,2)-2*(Dij1235I(2
     &   ,1)-Dij1235I(3,1)+Dij1235I(3,2)+Dij1235I(4,2)-Dij1235I(5,2)-
     &   Dij1235I(6,2))-Dij2345I(3,2)+Dij2345I(5,2)-8*(2*(EijI(11,2)+
     &   EijI(21,3))+EijI(22,3)-3*EijI(24,3)+EijI(39,4)+EijI(40,4)-Ei
     &   jI(42,4)-EijI(44,4))+(EijI(2,2)+EijI(8,3)-EijI(10,3))*P(1)-(
     &   4*(EijI(1,1)+EijI(1,2))-EijI(5,3)+EijI(7,3))*P(21)+(-EijI(8,
     &   2)+EijI(10,2)-EijI(16,3)-EijI(17,3)+EijI(19,3)+EijI(20,3))*P
     &   (233)+EijI(4,3)*P(385)-EijI(2,1)*P(671)+EijI(9,2)*P(672)-Eij
     &   I(14,3)*P(673)-EijI(15,3)*P(674)+EijI(18,3)*P(675))
       F7(270)=DCMPLX(FRe(270),FIm(270))
       P(676) = -4*p4sq-10*p5sq+8*s45
       P(677) = 2*p4sq+p5sq-s12+3*s15+s34-4*s45
       P(678) = 5*p4sq-s12+s34-9*s45+3*P(401)
       P(679) = p4sq+p5sq-s45-2*P(1)
       P(680) = s12+s15-s34-2*P(385)
       FRe(271) = EijR(4,1)*P(667)+EijR(4,2)*P(668)+EijR(2,1)*P(676)+2
     &   *(Dij1234R(2,1)+Dij1234R(2,2)+4*Dij1235R(6,2)-Dij2345R(1,1)-
     &   Dij2345R(1,2)+Dij2345R(3,1)-Dij2345R(3,2)-2*(Dij1235R(2,1)+D
     &   ij1235R(2,2)-Dij1235R(3,1)+Dij1235R(3,2)-Dij2345R(5,2)+4*(Ei
     &   jR(11,2)+EijR(37,4)+EijR(39,4)+2*(EijR(22,3)-EijR(24,3)-EijR
     &   (44,4))))+EijR(2,3)*P(1)+(-3*(EijR(5,2)-EijR(7,2))+EijR(8,3)
     &   +EijR(14,3)-2*EijR(18,3))*P(21)-(EijR(8,2)+EijR(9,3)-EijR(10
     &   ,2)+EijR(16,3)-2*EijR(20,3))*P(233)+EijR(4,3)*P(385)-EijR(2,
     &   2)*P(677)+EijR(9,2)*P(678)+EijR(10,3)*P(679)+EijR(15,3)*P(68
     &   0))
       FIm(271) = EijI(4,1)*P(667)+EijI(4,2)*P(668)+EijI(2,1)*P(676)+2
     &   *(Dij1234I(2,1)+Dij1234I(2,2)+4*Dij1235I(6,2)-Dij2345I(1,1)-
     &   Dij2345I(1,2)+Dij2345I(3,1)-Dij2345I(3,2)-2*(Dij1235I(2,1)+D
     &   ij1235I(2,2)-Dij1235I(3,1)+Dij1235I(3,2)-Dij2345I(5,2)+4*(Ei
     &   jI(11,2)+EijI(37,4)+EijI(39,4)+2*(EijI(22,3)-EijI(24,3)-EijI
     &   (44,4))))+EijI(2,3)*P(1)+(-3*(EijI(5,2)-EijI(7,2))+EijI(8,3)
     &   +EijI(14,3)-2*EijI(18,3))*P(21)-(EijI(8,2)+EijI(9,3)-EijI(10
     &   ,2)+EijI(16,3)-2*EijI(20,3))*P(233)+EijI(4,3)*P(385)-EijI(2,
     &   2)*P(677)+EijI(9,2)*P(678)+EijI(10,3)*P(679)+EijI(15,3)*P(68
     &   0))
       F7(271)=DCMPLX(FRe(271),FIm(271))
       P(681) = 4*p4sq+2*p5sq+s12-s34-6*s45
       P(682) = p5sq+s12-s34+2*P(64)
       P(683) = p5sq-s15+2*P(64)
       FRe(272) = EijR(4,1)*P(667)+EijR(4,2)*P(668)+EijR(3,1)*P(676)+2
     &   *(Dij1234R(3,1)+Dij1234R(6,2)-Dij2345R(2,1)+Dij2345R(3,1)-Di
     &   j2345R(3,2)-Dij2345R(4,2)+Dij2345R(5,2)+Dij2345R(6,2)-8*(Eij
     &   R(11,2)+EijR(23,3)-EijR(24,3)+EijR(39,4)+EijR(43,4)-EijR(44,
     &   4)-EijR(45,4))+(EijR(9,3)-EijR(10,3))*P(1)+(-3*(EijR(6,2)-Ei
     &   jR(7,2))+EijR(14,3)+EijR(17,3)-EijR(18,3)-EijR(19,3))*P(21)-
     &   (EijR(3,2)+EijR(12,3)-EijR(13,3))*P(233)+EijR(4,3)*P(385)-Ei
     &   jR(15,3)*P(674)+(-EijR(8,2)+EijR(9,2))*P(677)+EijR(10,2)*P(6
     &   81)-EijR(16,3)*P(682)+EijR(20,3)*P(683))
       FIm(272) = EijI(4,1)*P(667)+EijI(4,2)*P(668)+EijI(3,1)*P(676)+2
     &   *(Dij1234I(3,1)+Dij1234I(6,2)-Dij2345I(2,1)+Dij2345I(3,1)-Di
     &   j2345I(3,2)-Dij2345I(4,2)+Dij2345I(5,2)+Dij2345I(6,2)-8*(Eij
     &   I(11,2)+EijI(23,3)-EijI(24,3)+EijI(39,4)+EijI(43,4)-EijI(44,
     &   4)-EijI(45,4))+(EijI(9,3)-EijI(10,3))*P(1)+(-3*(EijI(6,2)-Ei
     &   jI(7,2))+EijI(14,3)+EijI(17,3)-EijI(18,3)-EijI(19,3))*P(21)-
     &   (EijI(3,2)+EijI(12,3)-EijI(13,3))*P(233)+EijI(4,3)*P(385)-Ei
     &   jI(15,3)*P(674)+(-EijI(8,2)+EijI(9,2))*P(677)+EijI(10,2)*P(6
     &   81)-EijI(16,3)*P(682)+EijI(20,3)*P(683))
       F7(272)=DCMPLX(FRe(272),FIm(272))
       P(684) = p4sq+s15-s23-2*s34+s45
       P(685) = 6*p5sq-4*P(684)
       FRe(273) = 4*(EijR(3,1)*P(286)-EijR(2,1)*P(292))+EijR(4,2)*P(66
     &   8)-2*(-Dij2345R(3,1)+Dij2345R(3,2)-Dij2345R(5,2)+2*(D01345R+
     &   D02345R-Dij1235R(3,1)+Dij1235R(3,2)-Dij1235R(6,2)-Dij1245R(2
     &   ,1)+Dij1245R(3,1)+EijR(11,2)+EijR(22,3)-5*EijR(24,3)+4*(EijR
     &   (39,4)-EijR(44,4)))+EijR(10,3)*P(1)-(3*EijR(7,2)+EijR(14,3)-
     &   EijR(18,3))*P(21)-(EijR(10,2)-EijR(16,3)+EijR(20,3))*P(233)-
     &   EijR(4,3)*P(385)+EijR(15,3)*P(674)-EijR(9,2)*P(677))+EijR(4,
     &   1)*P(685)
       FIm(273) = 4*(EijI(3,1)*P(286)-EijI(2,1)*P(292))+EijI(4,2)*P(66
     &   8)-2*(-Dij2345I(3,1)+Dij2345I(3,2)-Dij2345I(5,2)+2*(D01345I+
     &   D02345I-Dij1235I(3,1)+Dij1235I(3,2)-Dij1235I(6,2)-Dij1245I(2
     &   ,1)+Dij1245I(3,1)+EijI(11,2)+EijI(22,3)-5*EijI(24,3)+4*(EijI
     &   (39,4)-EijI(44,4)))+EijI(10,3)*P(1)-(3*EijI(7,2)+EijI(14,3)-
     &   EijI(18,3))*P(21)-(EijI(10,2)-EijI(16,3)+EijI(20,3))*P(233)-
     &   EijI(4,3)*P(385)+EijI(15,3)*P(674)-EijI(9,2)*P(677))+EijI(4,
     &   1)*P(685)
       F7(273)=DCMPLX(FRe(273),FIm(273))
       FRe(274) = -4*(EijR(2,1)-EijR(4,1)+EijR(5,2)-3*EijR(9,2)-6*(Eij
     &   R(4,3)-EijR(15,3))+2*(EijR(7,2)+EijR(14,3)-EijR(18,3))+4*(Ei
     &   jR(4,4)-EijR(14,4)-EijR(15,4)+EijR(32,4)))
       FIm(274) = -4*(EijI(2,1)-EijI(4,1)+EijI(5,2)-3*EijI(9,2)-6*(Eij
     &   I(4,3)-EijI(15,3))+2*(EijI(7,2)+EijI(14,3)-EijI(18,3))+4*(Ei
     &   jI(4,4)-EijI(14,4)-EijI(15,4)+EijI(32,4)))
       F7(274)=DCMPLX(FRe(274),FIm(274))
       FRe(275) = -8*(EijR(3,2)+EijR(4,2)+EijR(4,3)+EijR(13,3)+2*(EijR
     &   (4,4)-EijR(10,2)-EijR(13,4)-EijR(16,3)-3*(EijR(16,4)-EijR(22
     &   ,4))))
       FIm(275) = -8*(EijI(3,2)+EijI(4,2)+EijI(4,3)+EijI(13,3)+2*(EijI
     &   (4,4)-EijI(10,2)-EijI(13,4)-EijI(16,3)-3*(EijI(16,4)-EijI(22
     &   ,4))))
       F7(275)=DCMPLX(FRe(275),FIm(275))
       FRe(276) = -4*(EijR(3,2)-3*EijR(10,2)+2*(EijR(4,2)-EijR(13,3)+E
     &   ijR(16,3))-8*EijR(16,4)+4*(EijR(4,4)+EijR(22,4)))
       FIm(276) = -4*(EijI(3,2)-3*EijI(10,2)+2*(EijI(4,2)-EijI(13,3)+E
     &   ijI(16,3))-8*EijI(16,4)+4*(EijI(4,4)+EijI(22,4)))
       F7(276)=DCMPLX(FRe(276),FIm(276))
       FRe(277) = -8*(EijR(3,2)+EijR(4,2)+3*(EijR(4,3)+EijR(13,3))-6*E
     &   ijR(16,3)-4*EijR(16,4)+2*(EijR(4,4)-EijR(10,2)+EijR(22,4)))
       FIm(277) = -8*(EijI(3,2)+EijI(4,2)+3*(EijI(4,3)+EijI(13,3))-6*E
     &   ijI(16,3)-4*EijI(16,4)+2*(EijI(4,4)-EijI(10,2)+EijI(22,4)))
       F7(277)=DCMPLX(FRe(277),FIm(277))
       FRe(278) = -8*(EijR(4,2)+EijR(4,3)-EijR(10,2)-EijR(16,3)-4*EijR
     &   (16,4)+2*(EijR(4,4)+EijR(22,4)))
       FIm(278) = -8*(EijI(4,2)+EijI(4,3)-EijI(10,2)-EijI(16,3)-4*EijI
     &   (16,4)+2*(EijI(4,4)+EijI(22,4)))
       F7(278)=DCMPLX(FRe(278),FIm(278))
       FRe(279) = -8*(EijR(3,1)-EijR(4,1)-EijR(4,3)+EijR(6,2)-EijR(7,2
     &   )-EijR(14,3)+3*EijR(16,3)+EijR(19,3)+2*(EijR(4,4)-EijR(13,3)
     &   -EijR(14,4)+EijR(22,4)-EijR(30,4))-4*(EijR(16,4)-EijR(33,4))
     &   )
       FIm(279) = -8*(EijI(3,1)-EijI(4,1)-EijI(4,3)+EijI(6,2)-EijI(7,2
     &   )-EijI(14,3)+3*EijI(16,3)+EijI(19,3)+2*(EijI(4,4)-EijI(13,3)
     &   -EijI(14,4)+EijI(22,4)-EijI(30,4))-4*(EijI(16,4)-EijI(33,4))
     &   )
       F7(279)=DCMPLX(FRe(279),FIm(279))
       FRe(280) = -4*(EijR(4,2)-EijR(10,2)+4*(EijR(4,3)+EijR(4,4)-EijR
     &   (16,3)-EijR(16,4)))
       FIm(280) = -4*(EijI(4,2)-EijI(10,2)+4*(EijI(4,3)+EijI(4,4)-EijI
     &   (16,3)-EijI(16,4)))
       F7(280)=DCMPLX(FRe(280),FIm(280))
       FRe(281) = -8*(EijR(4,2)-EijR(10,2)+3*(EijR(4,3)-EijR(16,3))+2*
     &   (EijR(4,4)-EijR(16,4)))
       FIm(281) = -8*(EijI(4,2)-EijI(10,2)+3*(EijI(4,3)-EijI(16,3))+2*
     &   (EijI(4,4)-EijI(16,4)))
       F7(281)=DCMPLX(FRe(281),FIm(281))
       FRe(282) = 8*(EijR(4,3)+7*(EijR(6,2)-EijR(7,2))-EijR(16,3)+3*(E
     &   ijR(3,1)-EijR(4,1)+EijR(14,3)-EijR(19,3))-2*(EijR(4,4)-EijR(
     &   16,4)+EijR(19,4)-EijR(25,4))+4*(EijR(4,2)+EijR(6,3)-EijR(7,3
     &   )-EijR(10,2)+EijR(14,4)-EijR(33,4)))
       FIm(282) = 8*(EijI(4,3)+7*(EijI(6,2)-EijI(7,2))-EijI(16,3)+3*(E
     &   ijI(3,1)-EijI(4,1)+EijI(14,3)-EijI(19,3))-2*(EijI(4,4)-EijI(
     &   16,4)+EijI(19,4)-EijI(25,4))+4*(EijI(4,2)+EijI(6,3)-EijI(7,3
     &   )-EijI(10,2)+EijI(14,4)-EijI(33,4)))
       F7(282)=DCMPLX(FRe(282),FIm(282))
       FRe(283) = 8*(3*(EijR(3,2)+EijR(4,2))-EijR(4,3)-6*EijR(10,2)+4*
     &   EijR(11,3)-EijR(13,3)-2*(EijR(4,4)-EijR(14,4)-EijR(16,3)+4*E
     &   ijR(19,3)+EijR(22,4)-EijR(30,4)-2*(EijR(14,3)+EijR(16,4)-Eij
     &   R(33,4))))
       FIm(283) = 8*(3*(EijI(3,2)+EijI(4,2))-EijI(4,3)-6*EijI(10,2)+4*
     &   EijI(11,3)-EijI(13,3)-2*(EijI(4,4)-EijI(14,4)-EijI(16,3)+4*E
     &   ijI(19,3)+EijI(22,4)-EijI(30,4)-2*(EijI(14,3)+EijI(16,4)-Eij
     &   I(33,4))))
       F7(283)=DCMPLX(FRe(283),FIm(283))
       FRe(284) = -4*(EijR(3,1)-EijR(4,1)+EijR(6,2)-EijR(7,2)+2*(EijR(
     &   4,2)-EijR(10,2))-6*(EijR(4,3)-EijR(14,3)-EijR(16,3)+EijR(19,
     &   3))+4*(EijR(4,4)-EijR(16,4)+EijR(19,4)-EijR(25,4))-8*(EijR(1
     &   4,4)-EijR(33,4)))
       FIm(284) = -4*(EijI(3,1)-EijI(4,1)+EijI(6,2)-EijI(7,2)+2*(EijI(
     &   4,2)-EijI(10,2))-6*(EijI(4,3)-EijI(14,3)-EijI(16,3)+EijI(19,
     &   3))+4*(EijI(4,4)-EijI(16,4)+EijI(19,4)-EijI(25,4))-8*(EijI(1
     &   4,4)-EijI(33,4)))
       F7(284)=DCMPLX(FRe(284),FIm(284))
       FRe(285) = -4*(2*(EijR(3,1)-EijR(4,1)+EijR(6,2)-EijR(7,2))-5*(E
     &   ijR(4,2)-EijR(10,2))-6*(EijR(14,3)-EijR(19,3))+4*(EijR(4,4)-
     &   EijR(14,4)-EijR(16,4)+EijR(33,4)))
       FIm(285) = -4*(2*(EijI(3,1)-EijI(4,1)+EijI(6,2)-EijI(7,2))-5*(E
     &   ijI(4,2)-EijI(10,2))-6*(EijI(14,3)-EijI(19,3))+4*(EijI(4,4)-
     &   EijI(14,4)-EijI(16,4)+EijI(33,4)))
       F7(285)=DCMPLX(FRe(285),FIm(285))
       FRe(286) = -4*(EijR(3,1)-EijR(4,1)+10*EijR(4,2)+EijR(6,2)-EijR(
     &   7,2)-18*EijR(10,2)-6*(EijR(4,3)-EijR(14,3))+14*(EijR(16,3)-E
     &   ijR(19,3))+4*(EijR(4,4)-EijR(14,4)+EijR(22,4)-EijR(30,4))+8*
     &   (EijR(3,2)+EijR(11,3)-EijR(13,3)-EijR(16,4)+EijR(33,4)))
       FIm(286) = -4*(EijI(3,1)-EijI(4,1)+10*EijI(4,2)+EijI(6,2)-EijI(
     &   7,2)-18*EijI(10,2)-6*(EijI(4,3)-EijI(14,3))+14*(EijI(16,3)-E
     &   ijI(19,3))+4*(EijI(4,4)-EijI(14,4)+EijI(22,4)-EijI(30,4))+8*
     &   (EijI(3,2)+EijI(11,3)-EijI(13,3)-EijI(16,4)+EijI(33,4)))
       F7(286)=DCMPLX(FRe(286),FIm(286))
       FRe(287) = -4*(EijR(3,1)-EijR(4,1)+4*EijR(4,4)+EijR(6,2)-EijR(7
     &   ,2)+2*(EijR(4,2)-EijR(4,3)-EijR(10,2)+EijR(14,3)+EijR(16,3)-
     &   EijR(19,3)-2*(EijR(14,4)+EijR(16,4)-EijR(33,4))))
       FIm(287) = -4*(EijI(3,1)-EijI(4,1)+4*EijI(4,4)+EijI(6,2)-EijI(7
     &   ,2)+2*(EijI(4,2)-EijI(4,3)-EijI(10,2)+EijI(14,3)+EijI(16,3)-
     &   EijI(19,3)-2*(EijI(14,4)+EijI(16,4)-EijI(33,4))))
       F7(287)=DCMPLX(FRe(287),FIm(287))
       FRe(288) = 8*(-EijR(4,3)+3*(EijR(4,2)-EijR(10,2))+EijR(16,3)+4*
     &   (EijR(14,3)-EijR(19,3))-2*(EijR(4,4)-EijR(14,4)-EijR(16,4)+E
     &   ijR(33,4)))
       FIm(288) = 8*(-EijI(4,3)+3*(EijI(4,2)-EijI(10,2))+EijI(16,3)+4*
     &   (EijI(14,3)-EijI(19,3))-2*(EijI(4,4)-EijI(14,4)-EijI(16,4)+E
     &   ijI(33,4)))
       F7(288)=DCMPLX(FRe(288),FIm(288))
       FRe(289) = -8*(EijR(3,1)-EijR(4,1)+EijR(4,3)+EijR(6,2)-EijR(7,2
     &   )-EijR(16,3)-3*(EijR(14,3)-EijR(19,3))-2*(EijR(4,2)-EijR(4,4
     &   )-EijR(10,2)+EijR(14,4)+EijR(16,4)-EijR(33,4)))
       FIm(289) = -8*(EijI(3,1)-EijI(4,1)+EijI(4,3)+EijI(6,2)-EijI(7,2
     &   )-EijI(16,3)-3*(EijI(14,3)-EijI(19,3))-2*(EijI(4,2)-EijI(4,4
     &   )-EijI(10,2)+EijI(14,4)+EijI(16,4)-EijI(33,4)))
       F7(289)=DCMPLX(FRe(289),FIm(289))
       FRe(290) = -4*(EijR(2,1)+EijR(2,2)-EijR(4,1)-6*EijR(4,3)-EijR(9
     &   ,2)-2*EijR(10,3)+8*(EijR(15,3)-EijR(15,4))+4*(EijR(4,4)+EijR
     &   (21,4)))
       FIm(290) = -4*(EijI(2,1)+EijI(2,2)-EijI(4,1)-6*EijI(4,3)-EijI(9
     &   ,2)-2*EijI(10,3)+8*(EijI(15,3)-EijI(15,4))+4*(EijI(4,4)+EijI
     &   (21,4)))
       F7(290)=DCMPLX(FRe(290),FIm(290))
       FRe(291) = -4*(EijR(2,1)+EijR(2,2)-EijR(4,1)-10*EijR(4,3)-13*Ei
     &   jR(9,2)-6*EijR(10,3)+16*EijR(15,3)+12*(EijR(4,2)+EijR(16,3))
     &   -20*EijR(20,3)+4*(EijR(4,4)-EijR(16,4)+EijR(21,4)-EijR(28,4)
     &   )+8*(EijR(8,2)+EijR(9,3)-EijR(10,2)-EijR(15,4)+EijR(34,4)))
       FIm(291) = -4*(EijI(2,1)+EijI(2,2)-EijI(4,1)-10*EijI(4,3)-13*Ei
     &   jI(9,2)-6*EijI(10,3)+16*EijI(15,3)+12*(EijI(4,2)+EijI(16,3))
     &   -20*EijI(20,3)+4*(EijI(4,4)-EijI(16,4)+EijI(21,4)-EijI(28,4)
     &   )+8*(EijI(8,2)+EijI(9,3)-EijI(10,2)-EijI(15,4)+EijI(34,4)))
       F7(291)=DCMPLX(FRe(291),FIm(291))
       FRe(292) = 8*(3*(EijR(3,2)+EijR(4,2))-EijR(4,3)-6*EijR(10,2)+4*
     &   EijR(12,3)-EijR(13,3)-2*(EijR(4,4)-EijR(15,4)-EijR(16,3)+4*E
     &   ijR(20,3)+EijR(22,4)-EijR(31,4)-2*(EijR(15,3)+EijR(16,4)-Eij
     &   R(34,4))))
       FIm(292) = 8*(3*(EijI(3,2)+EijI(4,2))-EijI(4,3)-6*EijI(10,2)+4*
     &   EijI(12,3)-EijI(13,3)-2*(EijI(4,4)-EijI(15,4)-EijI(16,3)+4*E
     &   ijI(20,3)+EijI(22,4)-EijI(31,4)-2*(EijI(15,3)+EijI(16,4)-Eij
     &   I(34,4))))
       F7(292)=DCMPLX(FRe(292),FIm(292))
       FRe(293) = -4*(2*(EijR(3,1)-EijR(4,1)+EijR(8,2)-EijR(9,2))-5*(E
     &   ijR(4,2)-EijR(10,2))-6*(EijR(15,3)-EijR(20,3))+4*(EijR(4,4)-
     &   EijR(15,4)-EijR(16,4)+EijR(34,4)))
       FIm(293) = -4*(2*(EijI(3,1)-EijI(4,1)+EijI(8,2)-EijI(9,2))-5*(E
     &   ijI(4,2)-EijI(10,2))-6*(EijI(15,3)-EijI(20,3))+4*(EijI(4,4)-
     &   EijI(15,4)-EijI(16,4)+EijI(34,4)))
       F7(293)=DCMPLX(FRe(293),FIm(293))
       FRe(294) = -4*(EijR(3,1)-EijR(4,1)+10*EijR(4,2)+EijR(8,2)-EijR(
     &   9,2)-18*EijR(10,2)-6*(EijR(4,3)-EijR(15,3))+14*(EijR(16,3)-E
     &   ijR(20,3))+4*(EijR(4,4)-EijR(15,4)+EijR(22,4)-EijR(31,4))+8*
     &   (EijR(3,2)+EijR(12,3)-EijR(13,3)-EijR(16,4)+EijR(34,4)))
       FIm(294) = -4*(EijI(3,1)-EijI(4,1)+10*EijI(4,2)+EijI(8,2)-EijI(
     &   9,2)-18*EijI(10,2)-6*(EijI(4,3)-EijI(15,3))+14*(EijI(16,3)-E
     &   ijI(20,3))+4*(EijI(4,4)-EijI(15,4)+EijI(22,4)-EijI(31,4))+8*
     &   (EijI(3,2)+EijI(12,3)-EijI(13,3)-EijI(16,4)+EijI(34,4)))
       F7(294)=DCMPLX(FRe(294),FIm(294))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
