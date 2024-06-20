       subroutine NoAbe2gaPenLineFFhl10(F10)
       IMPLICIT NONE
       Real*8 P(767),FRe(379:426),FIm(379:426)
       Complex*16 F10(379:426)
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
       FRe(379) = -Dij2345R(3,1)+Dij2345R(3,2)-Dij2345R(5,2)-2*(Dij123
     &   5R(3,1)-Dij1235R(3,2)+Dij1235R(6,2)-EijR(22,3)-9*(EijR(11,2)
     &   -EijR(24,3)))+EijR(10,3)*P(1)+(EijR(7,2)+EijR(18,3))*P(21)+(
     &   EijR(10,2)+EijR(20,3))*P(230)+EijR(16,3)*P(233)+EijR(14,3)*P
     &   (386)+EijR(4,1)*P(540)+EijR(15,3)*P(674)+EijR(4,3)*P(740)+Ei
     &   jR(4,2)*P(741)+EijR(9,2)*P(750)
       FIm(379) = -Dij2345I(3,1)+Dij2345I(3,2)-Dij2345I(5,2)-2*(Dij123
     &   5I(3,1)-Dij1235I(3,2)+Dij1235I(6,2)-EijI(22,3)-9*(EijI(11,2)
     &   -EijI(24,3)))+EijI(10,3)*P(1)+(EijI(7,2)+EijI(18,3))*P(21)+(
     &   EijI(10,2)+EijI(20,3))*P(230)+EijI(16,3)*P(233)+EijI(14,3)*P
     &   (386)+EijI(4,1)*P(540)+EijI(15,3)*P(674)+EijI(4,3)*P(740)+Ei
     &   jI(4,2)*P(741)+EijI(9,2)*P(750)
       F10(379)=DCMPLX(FRe(379),FIm(379))
       FRe(380) = 2*(Dij1235R(3,2)-Dij1235R(6,2))+Dij2345R(3,2)-Dij234
     &   5R(5,2)+34*(EijR(22,3)-EijR(24,3))+EijR(10,3)*P(1)+EijR(18,3
     &   )*P(21)+EijR(20,3)*P(230)+EijR(16,3)*P(233)+EijR(14,3)*P(386
     &   )+EijR(9,2)*P(540)+EijR(4,2)*P(541)+EijR(15,3)*P(674)+EijR(4
     &   ,3)*P(740)
       FIm(380) = 2*(Dij1235I(3,2)-Dij1235I(6,2))+Dij2345I(3,2)-Dij234
     &   5I(5,2)+34*(EijI(22,3)-EijI(24,3))+EijI(10,3)*P(1)+EijI(18,3
     &   )*P(21)+EijI(20,3)*P(230)+EijI(16,3)*P(233)+EijI(14,3)*P(386
     &   )+EijI(9,2)*P(540)+EijI(4,2)*P(541)+EijI(15,3)*P(674)+EijI(4
     &   ,3)*P(740)
       F10(380)=DCMPLX(FRe(380),FIm(380))
       FRe(381) = Dij2345R(3,2)-Dij2345R(6,2)+34*(EijR(23,3)-EijR(24,3
     &   ))+EijR(20,3)*P(1)+EijR(19,3)*P(21)+EijR(13,3)*P(230)+EijR(1
     &   4,3)*P(386)+EijR(10,2)*P(540)+EijR(4,2)*P(541)+EijR(16,3)*P(
     &   682)+EijR(15,3)*P(721)+EijR(4,3)*P(740)
       FIm(381) = Dij2345I(3,2)-Dij2345I(6,2)+34*(EijI(23,3)-EijI(24,3
     &   ))+EijI(20,3)*P(1)+EijI(19,3)*P(21)+EijI(13,3)*P(230)+EijI(1
     &   4,3)*P(386)+EijI(10,2)*P(540)+EijI(4,2)*P(541)+EijI(16,3)*P(
     &   682)+EijI(15,3)*P(721)+EijI(4,3)*P(740)
       F10(381)=DCMPLX(FRe(381),FIm(381))
       FRe(382) = -2*(Dij1235R(3,1)-Dij1235R(3,2)+Dij1235R(5,2))+Dij23
     &   45R(3,2)+34*(EijR(11,2)+EijR(21,3)-EijR(24,3))+(EijR(9,2)+Ei
     &   jR(18,3))*P(1)+EijR(7,3)*P(21)+(EijR(10,2)+EijR(19,3))*P(230
     &   )+EijR(16,3)*P(233)+EijR(4,1)*P(540)+EijR(14,3)*P(673)+EijR(
     &   15,3)*P(721)+EijR(4,3)*P(740)+EijR(4,2)*P(741)+EijR(7,2)*P(7
     &   43)
       FIm(382) = -2*(Dij1235I(3,1)-Dij1235I(3,2)+Dij1235I(5,2))+Dij23
     &   45I(3,2)+34*(EijI(11,2)+EijI(21,3)-EijI(24,3))+(EijI(9,2)+Ei
     &   jI(18,3))*P(1)+EijI(7,3)*P(21)+(EijI(10,2)+EijI(19,3))*P(230
     &   )+EijI(16,3)*P(233)+EijI(4,1)*P(540)+EijI(14,3)*P(673)+EijI(
     &   15,3)*P(721)+EijI(4,3)*P(740)+EijI(4,2)*P(741)+EijI(7,2)*P(7
     &   43)
       F10(382)=DCMPLX(FRe(382),FIm(382))
       FRe(383) = 2*Dij1235R(3,2)+Dij2345R(3,2)-36*EijR(24,3)+EijR(16,
     &   3)*P(233)+EijR(14,3)*P(386)+EijR(4,2)*P(541)+EijR(15,3)*P(72
     &   1)+EijR(4,3)*P(740)
       FIm(383) = 2*Dij1235I(3,2)+Dij2345I(3,2)-36*EijI(24,3)+EijI(16,
     &   3)*P(233)+EijI(14,3)*P(386)+EijI(4,2)*P(541)+EijI(15,3)*P(72
     &   1)+EijI(4,3)*P(740)
       F10(383)=DCMPLX(FRe(383),FIm(383))
       P(753) = -s12-4*P(21)
       P(754) = p4sq+5*p5sq+2*s12-6*s15-s34-s45
       P(755) = s23-3*P(84)-4*P(663)
       P(756) = -p2sq+5*s12-s23+3*P(84)+4*P(187)
       FRe(384) = -D01345R-D02345R+Dij1245R(2,1)+Dij1245R(2,2)-Dij1245
     &   R(3,1)+Dij1245R(3,2)+4*Dij1345R(5,2)-Dij2345R(1,1)+Dij2345R(
     &   3,1)+6*EijR(11,2)+8*(EijR(22,3)-EijR(24,3))+EijR(3,1)*P(230)
     &   +EijR(2,2)*P(266)+EijR(10,2)*P(333)-2*(Dij1245R(6,2)+Dij1345
     &   R(1,1)+Dij1345R(1,2)-Dij1345R(3,1)+Dij1345R(3,2)+(EijR(1,1)+
     &   2*(EijR(5,2)-EijR(7,2)))*P(21)-EE0R*P(386))+EijR(8,2)*P(568)
     &   +EijR(2,1)*P(753)+EijR(4,1)*P(754)+EijR(4,2)*P(755)+EijR(9,2
     &   )*P(756)
       FIm(384) = -D01345I-D02345I+Dij1245I(2,1)+Dij1245I(2,2)-Dij1245
     &   I(3,1)+Dij1245I(3,2)+4*Dij1345I(5,2)-Dij2345I(1,1)+Dij2345I(
     &   3,1)+6*EijI(11,2)+8*(EijI(22,3)-EijI(24,3))+EijI(3,1)*P(230)
     &   +EijI(2,2)*P(266)+EijI(10,2)*P(333)-2*(Dij1245I(6,2)+Dij1345
     &   I(1,1)+Dij1345I(1,2)-Dij1345I(3,1)+Dij1345I(3,2)+(EijI(1,1)+
     &   2*(EijI(5,2)-EijI(7,2)))*P(21)-EE0I*P(386))+EijI(8,2)*P(568)
     &   +EijI(2,1)*P(753)+EijI(4,1)*P(754)+EijI(4,2)*P(755)+EijI(9,2
     &   )*P(756)
       F10(384)=DCMPLX(FRe(384),FIm(384))
       P(757) = -p4sq+4*p5sq+2*s12-6*s15-s34+3*s45
       P(758) = s23-4*P(19)-3*P(47)-2*P(49)
       P(759) = -p2sq+2*s12+s15-s34
       P(760) = p2sq-3*s15-7*s34-4*s45+6*P(199)+2*P(314)
       FRe(385) = -D01345R-D02345R-Dij1234R(3,1)+Dij1245R(2,1)+Dij1245
     &   R(2,2)-Dij1245R(3,1)+Dij1245R(3,2)-Dij1345R(1,1)-Dij1345R(2,
     &   1)-s12*EijR(2,1)+14*EijR(11,2)+8*(EijR(23,3)-EijR(24,3))+(-2
     &   *EijR(1,1)-5*(EijR(6,2)-EijR(7,2)))*P(21)+EijR(8,2)*P(270)-2
     &   *(Dij1245R(6,2)-Dij1345R(3,1)+Dij1345R(3,2)+Dij1345R(4,2)-Di
     &   j1345R(5,2)-Dij1345R(6,2)-EE0R*P(386))+EijR(3,1)*P(596)+EijR
     &   (3,2)*P(616)+EijR(4,1)*P(757)+EijR(4,2)*P(758)+EijR(9,2)*P(7
     &   59)+EijR(10,2)*P(760)
       FIm(385) = -D01345I-D02345I-Dij1234I(3,1)+Dij1245I(2,1)+Dij1245
     &   I(2,2)-Dij1245I(3,1)+Dij1245I(3,2)-Dij1345I(1,1)-Dij1345I(2,
     &   1)-s12*EijI(2,1)+14*EijI(11,2)+8*(EijI(23,3)-EijI(24,3))+(-2
     &   *EijI(1,1)-5*(EijI(6,2)-EijI(7,2)))*P(21)+EijI(8,2)*P(270)-2
     &   *(Dij1245I(6,2)-Dij1345I(3,1)+Dij1345I(3,2)+Dij1345I(4,2)-Di
     &   j1345I(5,2)-Dij1345I(6,2)-EE0I*P(386))+EijI(3,1)*P(596)+EijI
     &   (3,2)*P(616)+EijI(4,1)*P(757)+EijI(4,2)*P(758)+EijI(9,2)*P(7
     &   59)+EijI(10,2)*P(760)
       F10(385)=DCMPLX(FRe(385),FIm(385))
       P(761) = -s23-3*s45+7*P(21)+4*P(663)
       FRe(386) = -D01345R+Dij1245R(2,1)-Dij1245R(3,1)+Dij1245R(3,2)+D
     &   ij1245R(4,2)-Dij1245R(5,2)-Dij1245R(6,2)+Dij2345R(3,1)-s12*E
     &   ijR(2,1)+4*(Dij1345R(5,2)-EijR(11,2))+8*(EijR(21,3)-EijR(24,
     &   3))+EijR(3,1)*P(230)+EijR(5,2)*P(266)+EijR(10,2)*P(333)+EijR
     &   (9,2)*P(337)-2*(Dij1345R(1,1)+Dij1345R(1,2)-Dij1345R(3,1)+Di
     &   j1345R(3,2)+(3*EijR(1,1)+2*EijR(1,2))*P(21)-EE0R*P(386))+Eij
     &   R(6,2)*P(568)+EijR(4,1)*P(754)+EijR(4,2)*P(755)+EijR(7,2)*P(
     &   761)
       FIm(386) = -D01345I+Dij1245I(2,1)-Dij1245I(3,1)+Dij1245I(3,2)+D
     &   ij1245I(4,2)-Dij1245I(5,2)-Dij1245I(6,2)+Dij2345I(3,1)-s12*E
     &   ijI(2,1)+4*(Dij1345I(5,2)-EijI(11,2))+8*(EijI(21,3)-EijI(24,
     &   3))+EijI(3,1)*P(230)+EijI(5,2)*P(266)+EijI(10,2)*P(333)+EijI
     &   (9,2)*P(337)-2*(Dij1345I(1,1)+Dij1345I(1,2)-Dij1345I(3,1)+Di
     &   j1345I(3,2)+(3*EijI(1,1)+2*EijI(1,2))*P(21)-EE0I*P(386))+Eij
     &   I(6,2)*P(568)+EijI(4,1)*P(754)+EijI(4,2)*P(755)+EijI(7,2)*P(
     &   761)
       F10(386)=DCMPLX(FRe(386),FIm(386))
       P(762) = p4sq+4*p5sq-2*s12-s15-s23-s34
       FRe(387) = Dij1245R(3,2)-Dij1245R(6,2)+Dij1345R(1,1)+Dij2345R(3
     &   ,1)-26*EijR(11,2)-8*EijR(24,3)+2*(D01345R+D02345R-Dij1245R(2
     &   ,1)+Dij1245R(3,1)-Dij1345R(3,2)+Dij1345R(5,2)+(EE0R+EijR(1,1
     &   )+2*EijR(7,2))*P(21))+EijR(2,1)*P(324)+EijR(10,2)*P(333)+Eij
     &   R(9,2)*P(337)+EijR(3,1)*P(631)+EijR(4,2)*P(755)+EijR(4,1)*P(
     &   762)
       FIm(387) = Dij1245I(3,2)-Dij1245I(6,2)+Dij1345I(1,1)+Dij2345I(3
     &   ,1)-26*EijI(11,2)-8*EijI(24,3)+2*(D01345I+D02345I-Dij1245I(2
     &   ,1)+Dij1245I(3,1)-Dij1345I(3,2)+Dij1345I(5,2)+(EE0I+EijI(1,1
     &   )+2*EijI(7,2))*P(21))+EijI(2,1)*P(324)+EijI(10,2)*P(333)+Eij
     &   I(9,2)*P(337)+EijI(3,1)*P(631)+EijI(4,2)*P(755)+EijI(4,1)*P(
     &   762)
       F10(387)=DCMPLX(FRe(387),FIm(387))
       P(763) = 2*s12+s15-3*s23-s34
       FRe(388) = -Dij1245R(3,1)+Dij1245R(3,2)-Dij1245R(5,2)-7*Dij1345
     &   R(3,1)-2*(Dij1345R(3,2)-Dij1345R(5,2))+3*(D01345R+p2sq*EijR(
     &   2,1))+8*(Dij1345R(1,1)+p2sq*EijR(5,2)-p2sq*EijR(9,2))-6*EijR
     &   (11,2)-16*(EijR(21,3)-EijR(24,3))+4*(-EijR(4,2)+EijR(7,2))*P
     &   (4)+(-3*EijR(3,1)-8*(EijR(6,2)-EijR(10,2)))*P(264)+EijR(4,1)
     &   *P(763)
       FIm(388) = -Dij1245I(3,1)+Dij1245I(3,2)-Dij1245I(5,2)-7*Dij1345
     &   I(3,1)-2*(Dij1345I(3,2)-Dij1345I(5,2))+3*(D01345I+p2sq*EijI(
     &   2,1))+8*(Dij1345I(1,1)+p2sq*EijI(5,2)-p2sq*EijI(9,2))-6*EijI
     &   (11,2)-16*(EijI(21,3)-EijI(24,3))+4*(-EijI(4,2)+EijI(7,2))*P
     &   (4)+(-3*EijI(3,1)-8*(EijI(6,2)-EijI(10,2)))*P(264)+EijI(4,1)
     &   *P(763)
       F10(388)=DCMPLX(FRe(388),FIm(388))
       P(764) = 2*p2sq+s12+s15-4*s23-s34
       FRe(389) = -Dij1245R(3,1)+Dij1245R(3,2)-Dij1245R(6,2)-7*Dij1345
     &   R(3,1)-2*(Dij1345R(3,2)-Dij1345R(6,2))+3*(D01345R+p2sq*EijR(
     &   2,1))+8*(Dij1345R(2,1)+p2sq*EijR(8,2)-p2sq*EijR(9,2))-6*EijR
     &   (11,2)-16*(EijR(23,3)-EijR(24,3))-(3*EijR(3,1)+8*EijR(3,2))*
     &   P(264)+EijR(4,1)*P(763)-4*(EijR(4,2)*P(4)-EijR(10,2)*P(764))
       FIm(389) = -Dij1245I(3,1)+Dij1245I(3,2)-Dij1245I(6,2)-7*Dij1345
     &   I(3,1)-2*(Dij1345I(3,2)-Dij1345I(6,2))+3*(D01345I+p2sq*EijI(
     &   2,1))+8*(Dij1345I(2,1)+p2sq*EijI(8,2)-p2sq*EijI(9,2))-6*EijI
     &   (11,2)-16*(EijI(23,3)-EijI(24,3))-(3*EijI(3,1)+8*EijI(3,2))*
     &   P(264)+EijI(4,1)*P(763)-4*(EijI(4,2)*P(4)-EijI(10,2)*P(764))
       F10(389)=DCMPLX(FRe(389),FIm(389))
       P(765) = s12+s15-s34-2*P(488)
       FRe(390) = -Dij1245R(3,1)+Dij1245R(3,2)-Dij1245R(6,2)-7*Dij1345
     &   R(3,1)-2*(Dij1345R(3,2)-Dij1345R(5,2))+3*(D01345R+p2sq*EijR(
     &   2,1))+8*(Dij1345R(1,1)+p2sq*EijR(2,2))+18*EijR(11,2)-16*(Eij
     &   R(22,3)-EijR(24,3))+(-3*EijR(3,1)-8*(EijR(8,2)-EijR(10,2)))*
     &   P(264)+EijR(4,1)*P(763)-4*(EijR(4,2)*P(4)-EijR(9,2)*P(765))
       FIm(390) = -Dij1245I(3,1)+Dij1245I(3,2)-Dij1245I(6,2)-7*Dij1345
     &   I(3,1)-2*(Dij1345I(3,2)-Dij1345I(5,2))+3*(D01345I+p2sq*EijI(
     &   2,1))+8*(Dij1345I(1,1)+p2sq*EijI(2,2))+18*EijI(11,2)-16*(Eij
     &   I(22,3)-EijI(24,3))+(-3*EijI(3,1)-8*(EijI(8,2)-EijI(10,2)))*
     &   P(264)+EijI(4,1)*P(763)-4*(EijI(4,2)*P(4)-EijI(9,2)*P(765))
       F10(390)=DCMPLX(FRe(390),FIm(390))
       P(766) = -2*s12+s15+s23-s34
       FRe(391) = -D01345R+Dij1245R(3,2)-9*Dij1345R(3,1)-p2sq*(EijR(2,
     &   1)+8*EijR(9,2))+2*(Dij1245R(3,1)-Dij1345R(3,2)+EijR(11,2)+8*
     &   EijR(24,3))-4*EijR(4,2)*P(4)+(EijR(3,1)+8*EijR(10,2))*P(264)
     &   +EijR(4,1)*P(766)
       FIm(391) = -D01345I+Dij1245I(3,2)-9*Dij1345I(3,1)-p2sq*(EijI(2,
     &   1)+8*EijI(9,2))+2*(Dij1245I(3,1)-Dij1345I(3,2)+EijI(11,2)+8*
     &   EijI(24,3))-4*EijI(4,2)*P(4)+(EijI(3,1)+8*EijI(10,2))*P(264)
     &   +EijI(4,1)*P(766)
       F10(391)=DCMPLX(FRe(391),FIm(391))
       FRe(392) = 2*(EijR(2,1)-EijR(4,1)+EijR(5,2)-4*(EijR(4,2)+EijR(7
     &   ,2))+7*EijR(9,2)+8*(EijR(4,3)-EijR(14,3)-EijR(15,3)+EijR(18,
     &   3)))
       FIm(392) = 2*(EijI(2,1)-EijI(4,1)+EijI(5,2)-4*(EijI(4,2)+EijI(7
     &   ,2))+7*EijI(9,2)+8*(EijI(4,3)-EijI(14,3)-EijI(15,3)+EijI(18,
     &   3)))
       F10(392)=DCMPLX(FRe(392),FIm(392))
       FRe(393) = -2*(3*(EijR(2,1)-EijR(4,1))+12*(EijR(4,2)-EijR(7,2))
     &   +11*(EijR(5,2)-EijR(9,2))-8*(EijR(4,3)-EijR(5,3)+EijR(7,3)-E
     &   ijR(15,3)-2*(EijR(14,3)-EijR(18,3))))
       FIm(393) = -2*(3*(EijI(2,1)-EijI(4,1))+12*(EijI(4,2)-EijI(7,2))
     &   +11*(EijI(5,2)-EijI(9,2))-8*(EijI(4,3)-EijI(5,3)+EijI(7,3)-E
     &   ijI(15,3)-2*(EijI(14,3)-EijI(18,3))))
       F10(393)=DCMPLX(FRe(393),FIm(393))
       FRe(394) = 2*(EijR(2,1)+EijR(2,2)-EijR(4,1)-4*EijR(4,2)+3*EijR(
     &   9,2)+8*(EijR(4,3)+EijR(10,3))-16*EijR(15,3))
       FIm(394) = 2*(EijI(2,1)+EijI(2,2)-EijI(4,1)-4*EijI(4,2)+3*EijI(
     &   9,2)+8*(EijI(4,3)+EijI(10,3))-16*EijI(15,3))
       F10(394)=DCMPLX(FRe(394),FIm(394))
       FRe(395) = -2*(12*EijR(4,2)+3*(EijR(2,1)-EijR(4,1)+EijR(5,2))-4
     &   *EijR(7,2)-19*EijR(9,2)+8*(EijR(2,2)-EijR(4,3)+EijR(8,3)-Eij
     &   R(10,3)+EijR(14,3))+16*(EijR(15,3)-EijR(18,3)))
       FIm(395) = -2*(12*EijI(4,2)+3*(EijI(2,1)-EijI(4,1)+EijI(5,2))-4
     &   *EijI(7,2)-19*EijI(9,2)+8*(EijI(2,2)-EijI(4,3)+EijI(8,3)-Eij
     &   I(10,3)+EijI(14,3))+16*(EijI(15,3)-EijI(18,3)))
       F10(395)=DCMPLX(FRe(395),FIm(395))
       FRe(396) = -2*(3*(EijR(2,1)+EijR(2,2)-EijR(4,1))+12*EijR(4,2)-1
     &   5*EijR(9,2)-8*(EijR(4,3)-EijR(5,2)+EijR(7,2)-EijR(8,3)+EijR(
     &   10,3)-EijR(14,3))+16*(EijR(15,3)-EijR(18,3)))
       FIm(396) = -2*(3*(EijI(2,1)+EijI(2,2)-EijI(4,1))+12*EijI(4,2)-1
     &   5*EijI(9,2)-8*(EijI(4,3)-EijI(5,2)+EijI(7,2)-EijI(8,3)+EijI(
     &   10,3)-EijI(14,3))+16*(EijI(15,3)-EijI(18,3)))
       F10(396)=DCMPLX(FRe(396),FIm(396))
       FRe(397) = -2*(11*EijR(2,2)+3*(EijR(2,1)-EijR(4,1))+12*EijR(4,2
     &   )+8*(EijR(2,3)-EijR(4,3))-23*EijR(9,2)-24*(EijR(10,3)-EijR(1
     &   5,3)))
       FIm(397) = -2*(11*EijI(2,2)+3*(EijI(2,1)-EijI(4,1))+12*EijI(4,2
     &   )+8*(EijI(2,3)-EijI(4,3))-23*EijI(9,2)-24*(EijI(10,3)-EijI(1
     &   5,3)))
       F10(397)=DCMPLX(FRe(397),FIm(397))
       FRe(398) = 2*(EijR(3,1)-EijR(4,1)+EijR(8,2)-EijR(9,2)-4*(EijR(4
     &   ,2)-EijR(10,2))+8*(EijR(4,3)-EijR(15,3)-EijR(16,3)+EijR(20,3
     &   )))
       FIm(398) = 2*(EijI(3,1)-EijI(4,1)+EijI(8,2)-EijI(9,2)-4*(EijI(4
     &   ,2)-EijI(10,2))+8*(EijI(4,3)-EijI(15,3)-EijI(16,3)+EijI(20,3
     &   )))
       F10(398)=DCMPLX(FRe(398),FIm(398))
       FRe(399) = -2*(3*(EijR(3,1)-EijR(4,1)+EijR(8,2)-EijR(9,2))+12*(
     &   EijR(4,2)-EijR(10,2))-8*(EijR(4,3)-EijR(6,2)+EijR(7,2)-EijR(
     &   14,3)-EijR(15,3)-EijR(16,3)-EijR(17,3)+EijR(18,3)+EijR(19,3)
     &   +EijR(20,3)))
       FIm(399) = -2*(3*(EijI(3,1)-EijI(4,1)+EijI(8,2)-EijI(9,2))+12*(
     &   EijI(4,2)-EijI(10,2))-8*(EijI(4,3)-EijI(6,2)+EijI(7,2)-EijI(
     &   14,3)-EijI(15,3)-EijI(16,3)-EijI(17,3)+EijI(18,3)+EijI(19,3)
     &   +EijI(20,3)))
       F10(399)=DCMPLX(FRe(399),FIm(399))
       FRe(400) = -2*(12*EijR(4,2)+3*(EijR(2,1)-EijR(4,1)+EijR(5,2))-4
     &   *EijR(7,2)-11*EijR(9,2)-8*(EijR(4,3)-EijR(8,2)+EijR(10,2)-Ei
     &   jR(14,3)-EijR(15,3)-EijR(16,3)-EijR(17,3)+EijR(18,3)+EijR(19
     &   ,3)+EijR(20,3)))
       FIm(400) = -2*(12*EijI(4,2)+3*(EijI(2,1)-EijI(4,1)+EijI(5,2))-4
     &   *EijI(7,2)-11*EijI(9,2)-8*(EijI(4,3)-EijI(8,2)+EijI(10,2)-Ei
     &   jI(14,3)-EijI(15,3)-EijI(16,3)-EijI(17,3)+EijI(18,3)+EijI(19
     &   ,3)+EijI(20,3)))
       F10(400)=DCMPLX(FRe(400),FIm(400))
       FRe(401) = -2*(3*(EijR(3,1)-EijR(4,1))+11*(EijR(8,2)-EijR(9,2))
     &   +12*(EijR(4,2)-EijR(10,2))-8*(EijR(4,3)-EijR(9,3)+EijR(10,3)
     &   -EijR(16,3))+16*(EijR(15,3)-EijR(20,3)))
       FIm(401) = -2*(3*(EijI(3,1)-EijI(4,1))+11*(EijI(8,2)-EijI(9,2))
     &   +12*(EijI(4,2)-EijI(10,2))-8*(EijI(4,3)-EijI(9,3)+EijI(10,3)
     &   -EijI(16,3))+16*(EijI(15,3)-EijI(20,3)))
       F10(401)=DCMPLX(FRe(401),FIm(401))
       FRe(402) = -2*(3*(EijR(2,1)+EijR(2,2)-EijR(4,1))+12*EijR(4,2)-1
     &   5*EijR(9,2)-8*(EijR(4,3)-EijR(8,2)-EijR(9,3)+EijR(10,2)+EijR
     &   (10,3)-EijR(16,3)-2*(EijR(15,3)-EijR(20,3))))
       FIm(402) = -2*(3*(EijI(2,1)+EijI(2,2)-EijI(4,1))+12*EijI(4,2)-1
     &   5*EijI(9,2)-8*(EijI(4,3)-EijI(8,2)-EijI(9,3)+EijI(10,2)+EijI
     &   (10,3)-EijI(16,3)-2*(EijI(15,3)-EijI(20,3))))
       F10(402)=DCMPLX(FRe(402),FIm(402))
       FRe(403) = -2*(12*EijR(4,2)+8*(EijR(3,2)-EijR(4,3))+3*(EijR(3,1
     &   )-EijR(4,1)+EijR(8,2)-EijR(9,2))+4*(-5*EijR(10,2)+2*(EijR(12
     &   ,3)-EijR(13,3)+EijR(15,3)+2*(EijR(16,3)-EijR(20,3)))))
       FIm(403) = -2*(12*EijI(4,2)+8*(EijI(3,2)-EijI(4,3))+3*(EijI(3,1
     &   )-EijI(4,1)+EijI(8,2)-EijI(9,2))+4*(-5*EijI(10,2)+2*(EijI(12
     &   ,3)-EijI(13,3)+EijI(15,3)+2*(EijI(16,3)-EijI(20,3)))))
       F10(403)=DCMPLX(FRe(403),FIm(403))
       FRe(404) = 2*(EijR(4,1)+4*EijR(7,2)+3*EijR(9,2)-8*(EijR(4,2)-Ei
     &   jR(4,3)+EijR(14,3)+EijR(15,3)-EijR(18,3)))
       FIm(404) = 2*(EijI(4,1)+4*EijI(7,2)+3*EijI(9,2)-8*(EijI(4,2)-Ei
     &   jI(4,3)+EijI(14,3)+EijI(15,3)-EijI(18,3)))
       F10(404)=DCMPLX(FRe(404),FIm(404))
       FRe(405) = 2*(EijR(4,1)+7*EijR(9,2)-8*(EijR(4,2)-EijR(4,3)-EijR
     &   (10,3))-16*EijR(15,3))
       FIm(405) = 2*(EijI(4,1)+7*EijI(9,2)-8*(EijI(4,2)-EijI(4,3)-EijI
     &   (10,3))-16*EijI(15,3))
       F10(405)=DCMPLX(FRe(405),FIm(405))
       FRe(406) = 2*(EijR(4,1)+3*EijR(9,2)+4*EijR(10,2)-8*(EijR(4,2)-E
     &   ijR(4,3)+EijR(15,3)+EijR(16,3)-EijR(20,3)))
       FIm(406) = 2*(EijI(4,1)+3*EijI(9,2)+4*EijI(10,2)-8*(EijI(4,2)-E
     &   ijI(4,3)+EijI(15,3)+EijI(16,3)-EijI(20,3)))
       F10(406)=DCMPLX(FRe(406),FIm(406))
       FRe(407) = 2*(EijR(4,1)-EijR(9,2)+8*(EijR(4,3)-EijR(15,3)))
       FIm(407) = 2*(EijI(4,1)-EijI(9,2)+8*(EijI(4,3)-EijI(15,3)))
       F10(407)=DCMPLX(FRe(407),FIm(407))
       FRe(408) = -8*(3*EijR(11,2)+2*(EijR(22,3)-EijR(24,3)))
       FIm(408) = -8*(3*EijI(11,2)+2*(EijI(22,3)-EijI(24,3)))
       F10(408)=DCMPLX(FRe(408),FIm(408))
       FRe(409) = -4*(EijR(4,2)-EijR(9,2))-8*(EijR(4,3)-EijR(15,3))
       FIm(409) = -4*(EijI(4,2)-EijI(9,2))-8*(EijI(4,3)-EijI(15,3))
       F10(409)=DCMPLX(FRe(409),FIm(409))
       FRe(410) = 4*(EijR(4,2)-EijR(9,2)+2*(EijR(2,1)-EijR(4,1)-EijR(4
     &   ,3)+EijR(5,2)-EijR(7,2)+EijR(14,3)+EijR(15,3)-EijR(18,3)))
       FIm(410) = 4*(EijI(4,2)-EijI(9,2)+2*(EijI(2,1)-EijI(4,1)-EijI(4
     &   ,3)+EijI(5,2)-EijI(7,2)+EijI(14,3)+EijI(15,3)-EijI(18,3)))
       F10(410)=DCMPLX(FRe(410),FIm(410))
       FRe(411) = 4*(EijR(4,2)-EijR(9,2)-2*(EijR(4,3)+EijR(10,3))+4*Ei
     &   jR(15,3))
       FIm(411) = 4*(EijI(4,2)-EijI(9,2)-2*(EijI(4,3)+EijI(10,3))+4*Ei
     &   jI(15,3))
       F10(411)=DCMPLX(FRe(411),FIm(411))
       FRe(412) = -4*(EijR(3,1)-EijR(4,1)-EijR(4,2)+EijR(10,2)+2*(EijR
     &   (4,3)-EijR(16,3)))
       FIm(412) = -4*(EijI(3,1)-EijI(4,1)-EijI(4,2)+EijI(10,2)+2*(EijI
     &   (4,3)-EijI(16,3)))
       F10(412)=DCMPLX(FRe(412),FIm(412))
       FRe(413) = 4*(4*(EijR(6,2)-EijR(7,2))+3*(EijR(3,1)-EijR(4,1)+Ei
     &   jR(4,2)-EijR(10,2))-2*(EijR(4,3)-EijR(14,3)-EijR(16,3)+EijR(
     &   19,3)))
       FIm(413) = 4*(4*(EijI(6,2)-EijI(7,2))+3*(EijI(3,1)-EijI(4,1)+Ei
     &   jI(4,2)-EijI(10,2))-2*(EijI(4,3)-EijI(14,3)-EijI(16,3)+EijI(
     &   19,3)))
       F10(413)=DCMPLX(FRe(413),FIm(413))
       FRe(414) = 4*(EijR(4,2)-EijR(9,2)-2*(EijR(3,1)-EijR(4,1)+EijR(4
     &   ,3)+EijR(8,2)-EijR(10,2)-EijR(15,3)-EijR(16,3)+EijR(20,3)))
       FIm(414) = 4*(EijI(4,2)-EijI(9,2)-2*(EijI(3,1)-EijI(4,1)+EijI(4
     &   ,3)+EijI(8,2)-EijI(10,2)-EijI(15,3)-EijI(16,3)+EijI(20,3)))
       F10(414)=DCMPLX(FRe(414),FIm(414))
       FRe(415) = 4*(EijR(3,1)-EijR(4,1)+3*(EijR(4,2)-EijR(10,2))-2*(E
     &   ijR(4,3)-EijR(8,2)+EijR(9,2)-EijR(15,3)-EijR(16,3)+EijR(20,3
     &   )))
       FIm(415) = 4*(EijI(3,1)-EijI(4,1)+3*(EijI(4,2)-EijI(10,2))-2*(E
     &   ijI(4,3)-EijI(8,2)+EijI(9,2)-EijI(15,3)-EijI(16,3)+EijI(20,3
     &   )))
       F10(415)=DCMPLX(FRe(415),FIm(415))
       FRe(416) = -4*(EijR(3,1)-EijR(4,1)-3*(EijR(4,2)-EijR(10,2))+2*(
     &   EijR(4,3)+EijR(13,3))-4*EijR(16,3))
       FIm(416) = -4*(EijI(3,1)-EijI(4,1)-3*(EijI(4,2)-EijI(10,2))+2*(
     &   EijI(4,3)+EijI(13,3))-4*EijI(16,3))
       F10(416)=DCMPLX(FRe(416),FIm(416))
       FRe(417) = -2*(EE0R+EijR(1,1)-2*EijR(4,1)+6*(EijR(4,2)-EijR(7,2
     &   ))+4*(EijR(4,3)-EijR(14,3)))
       FIm(417) = -2*(EE0I+EijI(1,1)-2*EijI(4,1)+6*(EijI(4,2)-EijI(7,2
     &   ))+4*(EijI(4,3)-EijI(14,3)))
       F10(417)=DCMPLX(FRe(417),FIm(417))
       FRe(418) = -2*(EE0R+EijR(1,1)-2*(EijR(4,1)-EijR(4,2)+EijR(7,2)-
     &   2*(EijR(4,3)+EijR(7,3))+4*EijR(14,3)))
       FIm(418) = -2*(EE0I+EijI(1,1)-2*(EijI(4,1)-EijI(4,2)+EijI(7,2)-
     &   2*(EijI(4,3)+EijI(7,3))+4*EijI(14,3)))
       F10(418)=DCMPLX(FRe(418),FIm(418))
       FRe(419) = -2*(EE0R+EijR(1,1)+2*(-3*EijR(4,1)+EijR(4,2)-EijR(7,
     &   2)+2*(EijR(2,1)+EijR(4,3)+EijR(5,2)-EijR(9,2)-EijR(14,3)-Eij
     &   R(15,3)+EijR(18,3))))
       FIm(419) = -2*(EE0I+EijI(1,1)+2*(-3*EijI(4,1)+EijI(4,2)-EijI(7,
     &   2)+2*(EijI(2,1)+EijI(4,3)+EijI(5,2)-EijI(9,2)-EijI(14,3)-Eij
     &   I(15,3)+EijI(18,3))))
       F10(419)=DCMPLX(FRe(419),FIm(419))
       FRe(420) = -2*(EE0R+EijR(1,1)+2*(-5*EijR(4,1)+EijR(4,2)+4*(EijR
     &   (3,1)+EijR(6,2))-EijR(7,2)+2*(EijR(4,3)-2*EijR(10,2)-EijR(14
     &   ,3)-EijR(16,3)+EijR(19,3))))
       FIm(420) = -2*(EE0I+EijI(1,1)+2*(-5*EijI(4,1)+EijI(4,2)+4*(EijI
     &   (3,1)+EijI(6,2))-EijI(7,2)+2*(EijI(4,3)-2*EijI(10,2)-EijI(14
     &   ,3)-EijI(16,3)+EijI(19,3))))
       F10(420)=DCMPLX(FRe(420),FIm(420))
       FRe(421) = 2*(EijR(4,1)+2*EijR(9,2)-4*(EijR(4,3)-EijR(15,3)))
       FIm(421) = 2*(EijI(4,1)+2*EijI(9,2)-4*(EijI(4,3)-EijI(15,3)))
       F10(421)=DCMPLX(FRe(421),FIm(421))
       FRe(422) = -2*EijR(3,1)+12*EijR(10,2)+8*(EijR(4,1)-EijR(4,3)+Ei
     &   jR(16,3))
       FIm(422) = -2*EijI(3,1)+12*EijI(10,2)+8*(EijI(4,1)-EijI(4,3)+Ei
     &   jI(16,3))
       F10(422)=DCMPLX(FRe(422),FIm(422))
       FRe(423) = -2*(3*EijR(4,1)+2*(EE0R+EijR(1,1)+EijR(7,2))+4*(EijR
     &   (4,3)-EijR(14,3)))
       FIm(423) = -2*(3*EijI(4,1)+2*(EE0I+EijI(1,1)+EijI(7,2))+4*(EijI
     &   (4,3)-EijI(14,3)))
       F10(423)=DCMPLX(FRe(423),FIm(423))
       FRe(424) = -2*(EijR(4,1)+4*(EijR(4,2)+EijR(4,3)))
       FIm(424) = -2*(EijI(4,1)+4*(EijI(4,2)+EijI(4,3)))
       F10(424)=DCMPLX(FRe(424),FIm(424))
       P(767) = -p5sq+4*s34-2*P(702)
       FRe(425) = Dij2345R(3,1)-10*EijR(11,2)-8*EijR(24,3)-4*(D01345R+
     &   EijR(1,1)*P(21))+EijR(10,2)*P(233)-2*(D02345R-Dij1235R(3,1)-
     &   Dij1245R(2,1)+Dij1245R(3,1)+Dij1345R(1,1)-EijR(3,1)*P(286)+E
     &   ijR(2,1)*P(292))+(4*EE0R+EijR(7,2))*P(386)+EijR(9,2)*P(721)+
     &   EijR(4,2)*P(740)+EijR(4,1)*P(767)
       FIm(425) = Dij2345I(3,1)-10*EijI(11,2)-8*EijI(24,3)-4*(D01345I+
     &   EijI(1,1)*P(21))+EijI(10,2)*P(233)-2*(D02345I-Dij1235I(3,1)-
     &   Dij1245I(2,1)+Dij1245I(3,1)+Dij1345I(1,1)-EijI(3,1)*P(286)+E
     &   ijI(2,1)*P(292))+(4*EE0I+EijI(7,2))*P(386)+EijI(9,2)*P(721)+
     &   EijI(4,2)*P(740)+EijI(4,1)*P(767)
       F10(425)=DCMPLX(FRe(425),FIm(425))
       FRe(426) = (-2*Dij1235R(3,1)-Dij2345R(3,1)+34*EijR(11,2)+EijR(9
     &   ,2)*P(1)+EijR(7,2)*P(21)+EijR(10,2)*P(230)+EijR(4,2)*P(385)+
     &   EijR(4,1)*P(540))/2.
       FIm(426) = (-2*Dij1235I(3,1)-Dij2345I(3,1)+34*EijI(11,2)+EijI(9
     &   ,2)*P(1)+EijI(7,2)*P(21)+EijI(10,2)*P(230)+EijI(4,2)*P(385)+
     &   EijI(4,1)*P(540))/2.
       F10(426)=DCMPLX(FRe(426),FIm(426))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
