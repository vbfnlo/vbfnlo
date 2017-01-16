       subroutine NoAbe2gaPenLineFFhl9(F9)
       IMPLICIT NONE
       Real*8 P(767),FRe(337:378),FIm(337:378)
       Complex*16 F9(337:378)
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
       FRe(337) = -4*(EijR(4,2)+4*(EijR(4,3)+EijR(4,4)))
       FIm(337) = -4*(EijI(4,2)+4*(EijI(4,3)+EijI(4,4)))
       F9(337)=DCMPLX(FRe(337),FIm(337))
       FRe(338) = -4*(-Dij1245R(3,2)+Dij1245R(5,2)+2*(Dij1345R(1,1)-Di
     &   j1345R(3,1)+p2sq*(EijR(5,2)-EijR(9,2))+2*(EijR(11,2)+EijR(39
     &   ,4)-EijR(42,4))+(-EijR(6,2)+EijR(10,2))*P(264)-EijR(7,2)*P(3
     &   29)-EijR(4,2)*P(352)))
       FIm(338) = -4*(-Dij1245I(3,2)+Dij1245I(5,2)+2*(Dij1345I(1,1)-Di
     &   j1345I(3,1)+p2sq*(EijI(5,2)-EijI(9,2))+2*(EijI(11,2)+EijI(39
     &   ,4)-EijI(42,4))+(-EijI(6,2)+EijI(10,2))*P(264)-EijI(7,2)*P(3
     &   29)-EijI(4,2)*P(352)))
       F9(338)=DCMPLX(FRe(338),FIm(338))
       P(686) = p2sq+s15-2*s23-s34
       FRe(339) = 4*(Dij1245R(3,2)-Dij1245R(6,2)+8*(EijR(23,3)-EijR(24
     &   ,3))-4*(EijR(39,4)-EijR(45,4))-2*(Dij1345R(2,1)-Dij1345R(3,1
     &   )+p2sq*EijR(8,2)-p2sq*EijR(9,2)-EijR(3,2)*P(264)-EijR(4,2)*P
     &   (352)+EijR(10,2)*P(686)))
       FIm(339) = 4*(Dij1245I(3,2)-Dij1245I(6,2)+8*(EijI(23,3)-EijI(24
     &   ,3))-4*(EijI(39,4)-EijI(45,4))-2*(Dij1345I(2,1)-Dij1345I(3,1
     &   )+p2sq*EijI(8,2)-p2sq*EijI(9,2)-EijI(3,2)*P(264)-EijI(4,2)*P
     &   (352)+EijI(10,2)*P(686)))
       F9(339)=DCMPLX(FRe(339),FIm(339))
       P(687) = p2sq-s15+s23+s34
       FRe(340) = -4*(-Dij1245R(3,2)+Dij1245R(6,2)+8*EijR(11,2)+4*(Eij
     &   R(39,4)-EijR(44,4))+2*(Dij1345R(1,1)-Dij1345R(3,1)+p2sq*EijR
     &   (2,2)+(-EijR(8,2)+EijR(10,2))*P(264)-EijR(4,2)*P(352)-EijR(9
     &   ,2)*P(687)))
       FIm(340) = -4*(-Dij1245I(3,2)+Dij1245I(6,2)+8*EijI(11,2)+4*(Eij
     &   I(39,4)-EijI(44,4))+2*(Dij1345I(1,1)-Dij1345I(3,1)+p2sq*EijI
     &   (2,2)+(-EijI(8,2)+EijI(10,2))*P(264)-EijI(4,2)*P(352)-EijI(9
     &   ,2)*P(687)))
       F9(340)=DCMPLX(FRe(340),FIm(340))
       FRe(341) = 2*(Dij1245R(3,1)+2*(D01345R+Dij1245R(3,2)+p2sq*(EijR
     &   (2,1)+2*EijR(9,2))+2*(Dij1345R(3,1)-EijR(11,2)-3*EijR(24,3)-
     &   2*EijR(39,4))-(EijR(3,1)+2*EijR(10,2))*P(264)+(EijR(4,1)+2*E
     &   ijR(4,2))*P(352)))
       FIm(341) = 2*(Dij1245I(3,1)+2*(D01345I+Dij1245I(3,2)+p2sq*(EijI
     &   (2,1)+2*EijI(9,2))+2*(Dij1345I(3,1)-EijI(11,2)-3*EijI(24,3)-
     &   2*EijI(39,4))-(EijI(3,1)+2*EijI(10,2))*P(264)+(EijI(4,1)+2*E
     &   ijI(4,2))*P(352)))
       F9(341)=DCMPLX(FRe(341),FIm(341))
       P(688) = 2*s12+s15-s23-s34
       FRe(342) = 2*(D01345R-Dij1245R(3,1)+Dij1245R(3,2)-Dij1245R(5,2)
     &   +p2sq*(EijR(2,1)+4*(EijR(5,2)-EijR(9,2)))+6*EijR(11,2)+12*Ei
     &   jR(24,3)-8*(EijR(39,4)-EijR(42,4))+4*(Dij1345R(1,1)-Dij1345R
     &   (3,1)+(-EijR(4,2)+EijR(7,2))*P(3)+(-EijR(6,2)+EijR(10,2))*P(
     &   264))+EijR(3,1)*P(368)+EijR(4,1)*P(688))
       FIm(342) = 2*(D01345I-Dij1245I(3,1)+Dij1245I(3,2)-Dij1245I(5,2)
     &   +p2sq*(EijI(2,1)+4*(EijI(5,2)-EijI(9,2)))+6*EijI(11,2)+12*Ei
     &   jI(24,3)-8*(EijI(39,4)-EijI(42,4))+4*(Dij1345I(1,1)-Dij1345I
     &   (3,1)+(-EijI(4,2)+EijI(7,2))*P(3)+(-EijI(6,2)+EijI(10,2))*P(
     &   264))+EijI(3,1)*P(368)+EijI(4,1)*P(688))
       F9(342)=DCMPLX(FRe(342),FIm(342))
       P(689) = p2sq+s12+s15-2*s23-s34
       FRe(343) = 2*(D01345R-Dij1245R(3,1)+Dij1245R(3,2)-Dij1245R(6,2)
     &   +p2sq*(EijR(2,1)+4*(EijR(8,2)-EijR(9,2)))-2*(EijR(11,2)+8*Ei
     &   jR(23,3)-6*EijR(24,3)+4*(EijR(39,4)-EijR(45,4)))+EijR(3,1)*P
     &   (368)+EijR(4,1)*P(688)+4*(Dij1345R(2,1)-Dij1345R(3,1)-EijR(4
     &   ,2)*P(3)-EijR(3,2)*P(264)+EijR(10,2)*P(689)))
       FIm(343) = 2*(D01345I-Dij1245I(3,1)+Dij1245I(3,2)-Dij1245I(6,2)
     &   +p2sq*(EijI(2,1)+4*(EijI(8,2)-EijI(9,2)))-2*(EijI(11,2)+8*Ei
     &   jI(23,3)-6*EijI(24,3)+4*(EijI(39,4)-EijI(45,4)))+EijI(3,1)*P
     &   (368)+EijI(4,1)*P(688)+4*(Dij1345I(2,1)-Dij1345I(3,1)-EijI(4
     &   ,2)*P(3)-EijI(3,2)*P(264)+EijI(10,2)*P(689)))
       F9(343)=DCMPLX(FRe(343),FIm(343))
       P(690) = p2sq-s12-s15+s23+s34
       FRe(344) = 2*(D01345R-Dij1245R(3,1)+Dij1245R(3,2)-Dij1245R(6,2)
     &   +p2sq*EijR(2,1)+14*EijR(11,2)+12*EijR(24,3)-8*(EijR(39,4)-Ei
     &   jR(44,4))+EijR(3,1)*P(368)+EijR(4,1)*P(688)+4*(Dij1345R(1,1)
     &   -Dij1345R(3,1)+p2sq*EijR(2,2)-EijR(4,2)*P(3)+(-EijR(8,2)+Eij
     &   R(10,2))*P(264)-EijR(9,2)*P(690)))
       FIm(344) = 2*(D01345I-Dij1245I(3,1)+Dij1245I(3,2)-Dij1245I(6,2)
     &   +p2sq*EijI(2,1)+14*EijI(11,2)+12*EijI(24,3)-8*(EijI(39,4)-Ei
     &   jI(44,4))+EijI(3,1)*P(368)+EijI(4,1)*P(688)+4*(Dij1345I(1,1)
     &   -Dij1345I(3,1)+p2sq*EijI(2,2)-EijI(4,2)*P(3)+(-EijI(8,2)+Eij
     &   I(10,2))*P(264)-EijI(9,2)*P(690)))
       F9(344)=DCMPLX(FRe(344),FIm(344))
       P(691) = 2*s12-s15+s23+s34
       FRe(345) = 2*(D01345R+Dij1245R(3,2)+p2sq*EijR(2,1)+2*(Dij1245R(
     &   3,1)-EijR(11,2))-8*EijR(39,4)-4*(Dij1345R(3,1)+p2sq*EijR(9,2
     &   )-EijR(24,3)+EijR(4,2)*P(3)-EijR(10,2)*P(264))+EijR(3,1)*P(3
     &   68)-EijR(4,1)*P(691))
       FIm(345) = 2*(D01345I+Dij1245I(3,2)+p2sq*EijI(2,1)+2*(Dij1245I(
     &   3,1)-EijI(11,2))-8*EijI(39,4)-4*(Dij1345I(3,1)+p2sq*EijI(9,2
     &   )-EijI(24,3)+EijI(4,2)*P(3)-EijI(10,2)*P(264))+EijI(3,1)*P(3
     &   68)-EijI(4,1)*P(691))
       F9(345)=DCMPLX(FRe(345),FIm(345))
       P(692) = p2sq-s12-2*P(273)
       P(693) = p4sq-s15+s23+s45
       P(694) = p5sq+2*P(693)
       P(695) = p2sq+p4sq-s12-3*s15+s23+2*s34+s45
       P(696) = p5sq+2*P(695)
       P(697) = -p4sq-p5sq+s12+s15-s34+s45
       FRe(346) = 2*(-Dij2345R(3,2)+Dij2345R(5,2)-18*EijR(22,3)+26*Eij
     &   R(24,3)-8*(EijR(39,4)-EijR(44,4))-EijR(10,3)*P(1)+(-4*(EijR(
     &   5,2)-EijR(7,2))+EijR(14,3))*P(21)+4*(-EijR(2,1)+EijR(4,1))*P
     &   (181)+EijR(16,3)*P(230)+EijR(20,3)*P(233)+EijR(4,3)*P(385)+E
     &   ijR(18,3)*P(386)-2*(Dij1235R(3,2)-Dij1235R(6,2)+Dij1245R(2,1
     &   )-Dij1245R(3,1)-Dij1345R(1,1)+Dij1345R(3,1)+Dij2345R(1,1)-Di
     &   j2345R(3,1)-(-EijR(8,2)+EijR(10,2))*P(347)-EijR(2,2)*P(692))
     &   +EijR(4,2)*P(694)-EijR(9,2)*P(696)+EijR(15,3)*P(697))
       FIm(346) = 2*(-Dij2345I(3,2)+Dij2345I(5,2)-18*EijI(22,3)+26*Eij
     &   I(24,3)-8*(EijI(39,4)-EijI(44,4))-EijI(10,3)*P(1)+(-4*(EijI(
     &   5,2)-EijI(7,2))+EijI(14,3))*P(21)+4*(-EijI(2,1)+EijI(4,1))*P
     &   (181)+EijI(16,3)*P(230)+EijI(20,3)*P(233)+EijI(4,3)*P(385)+E
     &   ijI(18,3)*P(386)-2*(Dij1235I(3,2)-Dij1235I(6,2)+Dij1245I(2,1
     &   )-Dij1245I(3,1)-Dij1345I(1,1)+Dij1345I(3,1)+Dij2345I(1,1)-Di
     &   j2345I(3,1)-(-EijI(8,2)+EijI(10,2))*P(347)-EijI(2,2)*P(692))
     &   +EijI(4,2)*P(694)-EijI(9,2)*P(696)+EijI(15,3)*P(697))
       F9(346)=DCMPLX(FRe(346),FIm(346))
       P(698) = -4*P(266)+8*P(273)
       P(699) = s12-s15
       P(700) = 6*p4sq+p5sq+4*P(423)+2*P(699)
       P(701) = 4*p2sq-2*P(700)
       FRe(347) = -36*(EijR(23,3)-EijR(24,3))-16*(EijR(39,4)-EijR(45,4
     &   ))+8*(-EijR(3,1)+EijR(4,1))*P(181)-4*(Dij1245R(2,1)-Dij1245R
     &   (3,1)-Dij1345R(2,1)+Dij1345R(3,1)+Dij2345R(2,1)-Dij2345R(3,1
     &   )+EijR(3,2)*P(347)-EijR(8,2)*P(692))-2*(Dij2345R(3,2)-Dij234
     &   5R(6,2)-(EijR(15,3)-EijR(20,3))*P(1)-(-4*(EijR(6,2)-EijR(7,2
     &   ))+EijR(14,3))*P(21)-EijR(13,3)*P(233)-EijR(4,3)*P(385)-EijR
     &   (19,3)*P(386)+EijR(16,3)*P(682)-EijR(4,2)*P(694))+EijR(9,2)*
     &   P(698)+EijR(10,2)*P(701)
       FIm(347) = -36*(EijI(23,3)-EijI(24,3))-16*(EijI(39,4)-EijI(45,4
     &   ))+8*(-EijI(3,1)+EijI(4,1))*P(181)-4*(Dij1245I(2,1)-Dij1245I
     &   (3,1)-Dij1345I(2,1)+Dij1345I(3,1)+Dij2345I(2,1)-Dij2345I(3,1
     &   )+EijI(3,2)*P(347)-EijI(8,2)*P(692))-2*(Dij2345I(3,2)-Dij234
     &   5I(6,2)-(EijI(15,3)-EijI(20,3))*P(1)-(-4*(EijI(6,2)-EijI(7,2
     &   ))+EijI(14,3))*P(21)-EijI(13,3)*P(233)-EijI(4,3)*P(385)-EijI
     &   (19,3)*P(386)+EijI(16,3)*P(682)-EijI(4,2)*P(694))+EijI(9,2)*
     &   P(698)+EijI(10,2)*P(701)
       F9(347)=DCMPLX(FRe(347),FIm(347))
       P(702) = p4sq+s15-s23+s45
       P(703) = 6*p5sq-4*P(702)
       P(704) = -4*p2sq+2*s12+6*P(273)
       P(705) = p4sq+3*s45-2*P(126)
       P(706) = p4sq-p5sq+s23+s45
       P(707) = s15+2*P(706)
       P(708) = -s12+s45-3*P(187)+2*P(264)
       FRe(348) = 8*EE0R*P(386)+4*(EijR(3,1)*P(286)-EijR(6,2)*P(347)-E
     &   ijR(2,1)*P(630)+EijR(5,2)*P(692))+EijR(4,1)*P(703)+EijR(9,2)
     &   *P(704)+2*(-4*D01345R+2*(Dij1235R(3,1)-Dij1235R(3,2)+Dij1235
     &   R(5,2)-Dij1345R(3,1)+Dij2345R(3,1))-Dij2345R(3,2)-14*EijR(11
     &   ,2)-18*EijR(21,3)+26*EijR(24,3)-8*(EijR(39,4)-EijR(42,4))+(E
     &   ijR(15,3)-EijR(18,3))*P(1)-(8*EijR(1,1)+4*EijR(1,2)+EijR(7,3
     &   ))*P(21)+(-EijR(16,3)+EijR(19,3))*P(233)+EijR(4,3)*P(385)-Ei
     &   jR(14,3)*P(673)+EijR(4,2)*P(705)-EijR(7,2)*P(707)+EijR(10,2)
     &   *P(708))
       FIm(348) = 8*EE0I*P(386)+4*(EijI(3,1)*P(286)-EijI(6,2)*P(347)-E
     &   ijI(2,1)*P(630)+EijI(5,2)*P(692))+EijI(4,1)*P(703)+EijI(9,2)
     &   *P(704)+2*(-4*D01345I+2*(Dij1235I(3,1)-Dij1235I(3,2)+Dij1235
     &   I(5,2)-Dij1345I(3,1)+Dij2345I(3,1))-Dij2345I(3,2)-14*EijI(11
     &   ,2)-18*EijI(21,3)+26*EijI(24,3)-8*(EijI(39,4)-EijI(42,4))+(E
     &   ijI(15,3)-EijI(18,3))*P(1)-(8*EijI(1,1)+4*EijI(1,2)+EijI(7,3
     &   ))*P(21)+(-EijI(16,3)+EijI(19,3))*P(233)+EijI(4,3)*P(385)-Ei
     &   jI(14,3)*P(673)+EijI(4,2)*P(705)-EijI(7,2)*P(707)+EijI(10,2)
     &   *P(708))
       F9(348)=DCMPLX(FRe(348),FIm(348))
       FRe(349) = 4*EijR(10,2)*P(347)+8*(EijR(3,1)*P(264)+EijR(4,1)*P(
     &   507))-2*(2*(Dij1235R(3,2)+Dij1345R(3,1)-Dij2345R(3,1))+Dij23
     &   45R(3,2)+4*(D01345R+p2sq*EijR(2,1))-20*EijR(24,3)-8*(EijR(11
     &   ,2)-EijR(39,4))-EijR(15,3)*P(1)-(4*EijR(7,2)+EijR(14,3))*P(2
     &   1)+EijR(16,3)*P(233)-EijR(4,3)*P(385)-EijR(4,2)*P(694))+EijR
     &   (9,2)*P(698)
       FIm(349) = 4*EijI(10,2)*P(347)+8*(EijI(3,1)*P(264)+EijI(4,1)*P(
     &   507))-2*(2*(Dij1235I(3,2)+Dij1345I(3,1)-Dij2345I(3,1))+Dij23
     &   45I(3,2)+4*(D01345I+p2sq*EijI(2,1))-20*EijI(24,3)-8*(EijI(11
     &   ,2)-EijI(39,4))-EijI(15,3)*P(1)-(4*EijI(7,2)+EijI(14,3))*P(2
     &   1)+EijI(16,3)*P(233)-EijI(4,3)*P(385)-EijI(4,2)*P(694))+EijI
     &   (9,2)*P(698)
       F9(349)=DCMPLX(FRe(349),FIm(349))
       P(709) = -2*s45+4*P(507)
       P(710) = p5sq-2*s15+s23
       P(711) = 2*s45-4*P(710)
       P(712) = -4*P(21)-2*P(78)
       P(713) = 2*p2sq+s12-s23
       P(714) = p4sq+p5sq-2*s34-s45
       P(715) = -p4sq-p5sq+s12-s23+s45+2*P(310)
       P(716) = -2*P(184)-4*P(340)
       FRe(350) = -Dij1345R(1,2)-3*(Dij1345R(1,1)-Dij1345R(3,1))-Dij13
     &   45R(3,2)-4*(Dij1245R(6,2)-(EE0R+EijR(1,1))*P(21))+EijR(2,1)*
     &   P(709)+EijR(4,1)*P(711)+EijR(7,2)*P(712)+2*(D01345R-Dij1235R
     &   (2,1)+Dij1235R(3,1)+Dij1245R(2,2)+Dij1245R(3,2)+Dij1345R(5,2
     &   )-Dij2345R(1,1)+Dij2345R(3,1)-2*(EijR(11,2)-6*(EijR(22,3)-Ei
     &   jR(24,3)))+EijR(5,2)*P(145)+EijR(8,2)*P(341)-EijR(2,2)*P(713
     &   )+EijR(4,2)*P(714)+EijR(9,2)*P(715))+EijR(10,2)*P(716)
       FIm(350) = -Dij1345I(1,2)-3*(Dij1345I(1,1)-Dij1345I(3,1))-Dij13
     &   45I(3,2)-4*(Dij1245I(6,2)-(EE0I+EijI(1,1))*P(21))+EijI(2,1)*
     &   P(709)+EijI(4,1)*P(711)+EijI(7,2)*P(712)+2*(D01345I-Dij1235I
     &   (2,1)+Dij1235I(3,1)+Dij1245I(2,2)+Dij1245I(3,2)+Dij1345I(5,2
     &   )-Dij2345I(1,1)+Dij2345I(3,1)-2*(EijI(11,2)-6*(EijI(22,3)-Ei
     &   jI(24,3)))+EijI(5,2)*P(145)+EijI(8,2)*P(341)-EijI(2,2)*P(713
     &   )+EijI(4,2)*P(714)+EijI(9,2)*P(715))+EijI(10,2)*P(716)
       F9(350)=DCMPLX(FRe(350),FIm(350))
       P(717) = p2sq-s23-2*s34-s45
       P(718) = 3*p4sq+p5sq+s12+2*P(717)
       FRe(351) = -5*Dij1345R(2,1)+3*Dij1345R(3,1)-Dij1345R(3,2)-Dij13
     &   45R(4,2)+Dij1345R(5,2)+Dij1345R(6,2)+24*(EijR(23,3)-EijR(24,
     &   3))-4*(Dij1245R(6,2)+EijR(11,2)-(EE0R+EijR(1,1))*P(21))+EijR
     &   (3,1)*P(709)+EijR(4,1)*P(711)+EijR(7,2)*P(712)+2*(D01345R+Di
     &   j1245R(2,2)+Dij1245R(3,2)+Dij1345R(1,1)-Dij2345R(2,1)+Dij234
     &   5R(3,1)+EijR(6,2)*P(145)+EijR(3,2)*P(341)+(-EijR(8,2)+EijR(9
     &   ,2))*P(713)+EijR(4,2)*P(714)-EijR(10,2)*P(718))
       FIm(351) = -5*Dij1345I(2,1)+3*Dij1345I(3,1)-Dij1345I(3,2)-Dij13
     &   45I(4,2)+Dij1345I(5,2)+Dij1345I(6,2)+24*(EijI(23,3)-EijI(24,
     &   3))-4*(Dij1245I(6,2)+EijI(11,2)-(EE0I+EijI(1,1))*P(21))+EijI
     &   (3,1)*P(709)+EijI(4,1)*P(711)+EijI(7,2)*P(712)+2*(D01345I+Di
     &   j1245I(2,2)+Dij1245I(3,2)+Dij1345I(1,1)-Dij2345I(2,1)+Dij234
     &   5I(3,1)+EijI(6,2)*P(145)+EijI(3,2)*P(341)+(-EijI(8,2)+EijI(9
     &   ,2))*P(713)+EijI(4,2)*P(714)-EijI(10,2)*P(718))
       F9(351)=DCMPLX(FRe(351),FIm(351))
       P(719) = 2*p4sq+8*p5sq-7*s15+4*s23-6*s45
       P(720) = 3*P(21)+2*P(78)
       P(721) = -s12-s15+s34
       P(722) = -3*p4sq-6*p5sq+8*s15-4*s23+7*s45
       P(723) = -4*s34+3*P(385)
       P(724) = -4*p2sq-3*s12-s15+2*s23+s34
       P(725) = 3*P(184)+5*P(187)+4*P(264)
       P(726) = -6*p5sq-2*s23+4*s34+5*s45-3*P(47)
       P(727) = 4*p2sq+3*s12+s15-2*s23-s34
       P(728) = -3*P(184)-5*P(187)-4*P(264)
       FRe(352) = -D01234R-Dij1234R(1,1)-Dij1345R(1,2)-3*(Dij1345R(1,1
     &   )-Dij1345R(3,1))-Dij1345R(3,2)+2*(D01235R+D01345R+D02345R+Di
     &   j1245R(3,2)+Dij1245R(4,2)-Dij1245R(5,2)-Dij1245R(6,2)+Dij134
     &   5R(5,2))+Dij2345R(3,1)+4*EijR(11,2)+24*(EijR(21,3)-EijR(24,3
     &   ))+EijR(3,1)*P(233)+EE0R*P(492)+EijR(1,1)*P(719)+EijR(1,2)*P
     &   (720)+EijR(2,1)*P(721)+EijR(4,1)*P(722)+EijR(4,2)*P(723)+Eij
     &   R(5,2)*P(724)+EijR(6,2)*P(725)+EijR(7,2)*P(726)+EijR(9,2)*P(
     &   727)+EijR(10,2)*P(728)
       FIm(352) = -D01234I-Dij1234I(1,1)-Dij1345I(1,2)-3*(Dij1345I(1,1
     &   )-Dij1345I(3,1))-Dij1345I(3,2)+2*(D01235I+D01345I+D02345I+Di
     &   j1245I(3,2)+Dij1245I(4,2)-Dij1245I(5,2)-Dij1245I(6,2)+Dij134
     &   5I(5,2))+Dij2345I(3,1)+4*EijI(11,2)+24*(EijI(21,3)-EijI(24,3
     &   ))+EijI(3,1)*P(233)+EE0I*P(492)+EijI(1,1)*P(719)+EijI(1,2)*P
     &   (720)+EijI(2,1)*P(721)+EijI(4,1)*P(722)+EijI(4,2)*P(723)+Eij
     &   I(5,2)*P(724)+EijI(6,2)*P(725)+EijI(7,2)*P(726)+EijI(9,2)*P(
     &   727)+EijI(10,2)*P(728)
       F9(352)=DCMPLX(FRe(352),FIm(352))
       P(729) = p5sq-s15-2*s23+s45
       P(730) = -s12+s45-2*P(340)
       P(731) = -p4sq+7*s15-2*s34+5*s45-4*P(404)
       P(732) = -3*P(21)-2*P(78)
       FRe(353) = D01235R+4*D01345R+D02345R-Dij1245R(2,1)+Dij1245R(3,1
     &   )+2*(Dij1245R(3,2)-Dij1245R(6,2)+Dij1345R(1,1))+3*Dij1345R(3
     &   ,1)-Dij1345R(3,2)+Dij1345R(5,2)+Dij2345R(3,1)+6*EijR(11,2)-2
     &   4*EijR(24,3)+EijR(1,1)*P(43)+EijR(2,1)*P(713)+EijR(4,2)*P(72
     &   3)+EijR(9,2)*P(727)+EijR(10,2)*P(728)+EE0R*P(729)+EijR(3,1)*
     &   P(730)+EijR(4,1)*P(731)+EijR(7,2)*P(732)
       FIm(353) = D01235I+4*D01345I+D02345I-Dij1245I(2,1)+Dij1245I(3,1
     &   )+2*(Dij1245I(3,2)-Dij1245I(6,2)+Dij1345I(1,1))+3*Dij1345I(3
     &   ,1)-Dij1345I(3,2)+Dij1345I(5,2)+Dij2345I(3,1)+6*EijI(11,2)-2
     &   4*EijI(24,3)+EijI(1,1)*P(43)+EijI(2,1)*P(713)+EijI(4,2)*P(72
     &   3)+EijI(9,2)*P(727)+EijI(10,2)*P(728)+EE0I*P(729)+EijI(3,1)*
     &   P(730)+EijI(4,1)*P(731)+EijI(7,2)*P(732)
       F9(353)=DCMPLX(FRe(353),FIm(353))
       FRe(354) = 2*(Dij1245R(3,2)-Dij1245R(5,2))-8*Dij1345R(1,1)+6*Di
     &   j1345R(3,1)-Dij1345R(3,2)+Dij1345R(5,2)-4*(D01345R-s23*EijR(
     &   4,1)+p2sq*(EijR(2,1)+2*(EijR(5,2)-EijR(9,2)))-2*EijR(11,2)-4
     &   *(EijR(21,3)-EijR(24,3))-(EijR(3,1)+2*(EijR(6,2)-EijR(10,2))
     &   )*P(264)-EijR(7,2)*P(329)-EijR(4,2)*P(352))
       FIm(354) = 2*(Dij1245I(3,2)-Dij1245I(5,2))-8*Dij1345I(1,1)+6*Di
     &   j1345I(3,1)-Dij1345I(3,2)+Dij1345I(5,2)-4*(D01345I-s23*EijI(
     &   4,1)+p2sq*(EijI(2,1)+2*(EijI(5,2)-EijI(9,2)))-2*EijI(11,2)-4
     &   *(EijI(21,3)-EijI(24,3))-(EijI(3,1)+2*(EijI(6,2)-EijI(10,2))
     &   )*P(264)-EijI(7,2)*P(329)-EijI(4,2)*P(352))
       F9(354)=DCMPLX(FRe(354),FIm(354))
       P(733) = -2*p2sq-s15+3*s23+s34
       FRe(355) = 2*(Dij1245R(3,2)-Dij1245R(6,2))-8*Dij1345R(2,1)+6*Di
     &   j1345R(3,1)-Dij1345R(3,2)+Dij1345R(6,2)-4*(D01345R-s23*EijR(
     &   4,1)+p2sq*(EijR(2,1)+2*(EijR(8,2)-EijR(9,2)))-2*EijR(11,2)-4
     &   *(EijR(23,3)-EijR(24,3))-(EijR(3,1)+2*EijR(3,2))*P(264)-EijR
     &   (4,2)*P(352)-EijR(10,2)*P(733))
       FIm(355) = 2*(Dij1245I(3,2)-Dij1245I(6,2))-8*Dij1345I(2,1)+6*Di
     &   j1345I(3,1)-Dij1345I(3,2)+Dij1345I(6,2)-4*(D01345I-s23*EijI(
     &   4,1)+p2sq*(EijI(2,1)+2*(EijI(8,2)-EijI(9,2)))-2*EijI(11,2)-4
     &   *(EijI(23,3)-EijI(24,3))-(EijI(3,1)+2*EijI(3,2))*P(264)-EijI
     &   (4,2)*P(352)-EijI(10,2)*P(733))
       F9(355)=DCMPLX(FRe(355),FIm(355))
       P(734) = s23+2*P(273)
       FRe(356) = D01345R+Dij1245R(3,1)+6*Dij1345R(3,1)-Dij1345R(3,2)+
     &   p2sq*(EijR(2,1)+8*EijR(9,2))+2*(Dij1245R(3,2)-EijR(11,2)-8*E
     &   ijR(24,3))-8*EijR(10,2)*P(264)+4*EijR(4,2)*P(352)+EijR(3,1)*
     &   P(368)+EijR(4,1)*P(734)
       FIm(356) = D01345I+Dij1245I(3,1)+6*Dij1345I(3,1)-Dij1345I(3,2)+
     &   p2sq*(EijI(2,1)+8*EijI(9,2))+2*(Dij1245I(3,2)-EijI(11,2)-8*E
     &   ijI(24,3))-8*EijI(10,2)*P(264)+4*EijI(4,2)*P(352)+EijI(3,1)*
     &   P(368)+EijI(4,1)*P(734)
       F9(356)=DCMPLX(FRe(356),FIm(356))
       P(735) = 2*p4sq+p5sq+4*P(78)
       P(736) = -p4sq-p5sq+s45+4*P(352)
       P(737) = 8*p2sq-s12+4*s23-5*P(273)
       P(738) = p4sq+s12-s34-s45-8*P(264)
       FRe(357) = -8*Dij1345R(1,1)+6*Dij1345R(3,1)-Dij1345R(3,2)+Dij13
     &   45R(5,2)+Dij2345R(3,1)+2*(Dij1235R(3,1)+Dij1245R(3,2)-Dij124
     &   5R(6,2)-2*p2sq*EijR(2,1)-4*p2sq*EijR(2,2)-9*EijR(11,2)+8*(Ei
     &   jR(22,3)-EijR(24,3)))-4*(D01345R-(EijR(3,1)+2*EijR(8,2))*P(2
     &   64))+EijR(7,2)*P(386)+EijR(4,1)*P(735)+EijR(4,2)*P(736)+EijR
     &   (9,2)*P(737)+EijR(10,2)*P(738)
       FIm(357) = -8*Dij1345I(1,1)+6*Dij1345I(3,1)-Dij1345I(3,2)+Dij13
     &   45I(5,2)+Dij2345I(3,1)+2*(Dij1235I(3,1)+Dij1245I(3,2)-Dij124
     &   5I(6,2)-2*p2sq*EijI(2,1)-4*p2sq*EijI(2,2)-9*EijI(11,2)+8*(Ei
     &   jI(22,3)-EijI(24,3)))-4*(D01345I-(EijI(3,1)+2*EijI(8,2))*P(2
     &   64))+EijI(7,2)*P(386)+EijI(4,1)*P(735)+EijI(4,2)*P(736)+EijI
     &   (9,2)*P(737)+EijI(10,2)*P(738)
       F9(357)=DCMPLX(FRe(357),FIm(357))
       FRe(358) = 8*(-EijR(7,2)-EijR(9,2)+2*(EijR(4,2)+EijR(15,3)-EijR
     &   (18,3)))
       FIm(358) = 8*(-EijI(7,2)-EijI(9,2)+2*(EijI(4,2)+EijI(15,3)-EijI
     &   (18,3)))
       F9(358)=DCMPLX(FRe(358),FIm(358))
       FRe(359) = 8*(EijR(2,1)-EijR(4,1)+3*(EijR(5,2)-EijR(9,2))+2*(Ei
     &   jR(4,2)+EijR(5,3)-EijR(7,2)+EijR(15,3))-4*EijR(18,3))
       FIm(359) = 8*(EijI(2,1)-EijI(4,1)+3*(EijI(5,2)-EijI(9,2))+2*(Ei
     &   jI(4,2)+EijI(5,3)-EijI(7,2)+EijI(15,3))-4*EijI(18,3))
       F9(359)=DCMPLX(FRe(359),FIm(359))
       FRe(360) = EijR(4,2)-EijR(9,2)-EijR(10,3)+EijR(15,3)
       FIm(360) = EijI(4,2)-EijI(9,2)-EijI(10,3)+EijI(15,3)
       F9(360)=DCMPLX(FRe(360),FIm(360))
       FRe(361) = 8*(EijR(2,1)-EijR(4,1)+EijR(5,2)-EijR(7,2)+2*(EijR(2
     &   ,2)+EijR(4,2)+EijR(8,3)-2*EijR(9,2)-EijR(10,3)+EijR(15,3)-Ei
     &   jR(18,3)))
       FIm(361) = 8*(EijI(2,1)-EijI(4,1)+EijI(5,2)-EijI(7,2)+2*(EijI(2
     &   ,2)+EijI(4,2)+EijI(8,3)-2*EijI(9,2)-EijI(10,3)+EijI(15,3)-Ei
     &   jI(18,3)))
       F9(361)=DCMPLX(FRe(361),FIm(361))
       FRe(362) = 8*(EijR(2,1)+EijR(2,2)-EijR(4,1)-EijR(7,2)-4*EijR(9,
     &   2)+2*(EijR(4,2)+EijR(5,2)+EijR(8,3)-EijR(10,3)+EijR(15,3)-Ei
     &   jR(18,3)))
       FIm(362) = 8*(EijI(2,1)+EijI(2,2)-EijI(4,1)-EijI(7,2)-4*EijI(9,
     &   2)+2*(EijI(4,2)+EijI(5,2)+EijI(8,3)-EijI(10,3)+EijI(15,3)-Ei
     &   jI(18,3)))
       F9(362)=DCMPLX(FRe(362),FIm(362))
       FRe(363) = 8*(EijR(2,1)+3*EijR(2,2)-EijR(4,1)-5*EijR(9,2)-4*Eij
     &   R(10,3)+2*(EijR(2,3)+EijR(4,2)+EijR(15,3)))
       FIm(363) = 8*(EijI(2,1)+3*EijI(2,2)-EijI(4,1)-5*EijI(9,2)-4*Eij
     &   I(10,3)+2*(EijI(2,3)+EijI(4,2)+EijI(15,3)))
       F9(363)=DCMPLX(FRe(363),FIm(363))
       FRe(364) = 8*(-EijR(9,2)-EijR(10,2)+2*(EijR(4,2)+EijR(15,3)-Eij
     &   R(20,3)))
       FIm(364) = 8*(-EijI(9,2)-EijI(10,2)+2*(EijI(4,2)+EijI(15,3)-Eij
     &   I(20,3)))
       F9(364)=DCMPLX(FRe(364),FIm(364))
       FRe(365) = 8*(EijR(3,1)-EijR(4,1)-EijR(7,2)+EijR(8,2)-EijR(9,2)
     &   -3*EijR(10,2)+2*(EijR(4,2)+EijR(6,2)+EijR(15,3)+EijR(17,3)-E
     &   ijR(18,3)-EijR(20,3)))
       FIm(365) = 8*(EijI(3,1)-EijI(4,1)-EijI(7,2)+EijI(8,2)-EijI(9,2)
     &   -3*EijI(10,2)+2*(EijI(4,2)+EijI(6,2)+EijI(15,3)+EijI(17,3)-E
     &   ijI(18,3)-EijI(20,3)))
       F9(365)=DCMPLX(FRe(365),FIm(365))
       FRe(366) = 8*(EijR(2,1)-EijR(4,1)+EijR(5,2)-EijR(7,2)-3*EijR(9,
     &   2)-EijR(10,2)+2*(EijR(4,2)+EijR(8,2)+EijR(15,3)+EijR(17,3)-E
     &   ijR(18,3)-EijR(20,3)))
       FIm(366) = 8*(EijI(2,1)-EijI(4,1)+EijI(5,2)-EijI(7,2)-3*EijI(9,
     &   2)-EijI(10,2)+2*(EijI(4,2)+EijI(8,2)+EijI(15,3)+EijI(17,3)-E
     &   ijI(18,3)-EijI(20,3)))
       F9(366)=DCMPLX(FRe(366),FIm(366))
       FRe(367) = 8*(EijR(3,1)-EijR(4,1)+3*(EijR(8,2)-EijR(10,2))+2*(E
     &   ijR(4,2)-EijR(9,2)+EijR(9,3)-EijR(10,3)+EijR(15,3)-EijR(20,3
     &   )))
       FIm(367) = 8*(EijI(3,1)-EijI(4,1)+3*(EijI(8,2)-EijI(10,2))+2*(E
     &   ijI(4,2)-EijI(9,2)+EijI(9,3)-EijI(10,3)+EijI(15,3)-EijI(20,3
     &   )))
       F9(367)=DCMPLX(FRe(367),FIm(367))
       FRe(368) = 8*(EijR(2,1)+EijR(2,2)-EijR(4,1)-4*EijR(9,2)-EijR(10
     &   ,2)+2*(EijR(4,2)+EijR(8,2)+EijR(9,3)-EijR(10,3)+EijR(15,3)-E
     &   ijR(20,3)))
       FIm(368) = 8*(EijI(2,1)+EijI(2,2)-EijI(4,1)-4*EijI(9,2)-EijI(10
     &   ,2)+2*(EijI(4,2)+EijI(8,2)+EijI(9,3)-EijI(10,3)+EijI(15,3)-E
     &   ijI(20,3)))
       F9(368)=DCMPLX(FRe(368),FIm(368))
       FRe(369) = 8*(EijR(3,1)-EijR(4,1)+EijR(8,2)-EijR(9,2)+2*(EijR(3
     &   ,2)+EijR(4,2)+EijR(12,3)+EijR(15,3)-2*(EijR(10,2)+EijR(20,3)
     &   )))
       FIm(369) = 8*(EijI(3,1)-EijI(4,1)+EijI(8,2)-EijI(9,2)+2*(EijI(3
     &   ,2)+EijI(4,2)+EijI(12,3)+EijI(15,3)-2*(EijI(10,2)+EijI(20,3)
     &   )))
       F9(369)=DCMPLX(FRe(369),FIm(369))
       FRe(370) = -2*(EijR(2,1)+EijR(4,1)+EijR(5,2)+3*EijR(7,2)-7*(Eij
     &   R(4,2)-EijR(9,2))-8*(EijR(15,3)-EijR(18,3)))
       FIm(370) = -2*(EijI(2,1)+EijI(4,1)+EijI(5,2)+3*EijI(7,2)-7*(Eij
     &   I(4,2)-EijI(9,2))-8*(EijI(15,3)-EijI(18,3)))
       F9(370)=DCMPLX(FRe(370),FIm(370))
       FRe(371) = -2*(EijR(2,1)+EijR(2,2)+EijR(4,1)-7*EijR(4,2)+10*Eij
     &   R(9,2)+8*(EijR(10,3)-EijR(15,3)))
       FIm(371) = -2*(EijI(2,1)+EijI(2,2)+EijI(4,1)-7*EijI(4,2)+10*Eij
     &   I(9,2)+8*(EijI(10,3)-EijI(15,3)))
       F9(371)=DCMPLX(FRe(371),FIm(371))
       FRe(372) = -2*(EijR(3,1)+EijR(4,1)-7*EijR(4,2)+EijR(8,2)-EijR(9
     &   ,2)+11*EijR(10,2)-8*(EijR(15,3)-EijR(20,3)))
       FIm(372) = -2*(EijI(3,1)+EijI(4,1)-7*EijI(4,2)+EijI(8,2)-EijI(9
     &   ,2)+11*EijI(10,2)-8*(EijI(15,3)-EijI(20,3)))
       F9(372)=DCMPLX(FRe(372),FIm(372))
       FRe(373) = 2*(2*EijR(4,1)+7*EijR(4,2)+EijR(9,2)+8*EijR(15,3))
       FIm(373) = 2*(2*EijI(4,1)+7*EijI(4,2)+EijI(9,2)+8*EijI(15,3))
       F9(373)=DCMPLX(FRe(373),FIm(373))
       FRe(374) = -2*Dij1235R(3,1)-Dij2345R(3,1)+26*EijR(11,2)+16*EijR
     &   (22,3)+EijR(9,2)*P(1)+EijR(7,2)*P(21)+EijR(10,2)*P(230)+EijR
     &   (4,2)*P(385)+EijR(4,1)*P(540)
       FIm(374) = -2*Dij1235I(3,1)-Dij2345I(3,1)+26*EijI(11,2)+16*EijI
     &   (22,3)+EijI(9,2)*P(1)+EijI(7,2)*P(21)+EijI(10,2)*P(230)+EijI
     &   (4,2)*P(385)+EijI(4,1)*P(540)
       F9(374)=DCMPLX(FRe(374),FIm(374))
       P(739) = 2*p4sq+p5sq+13*s12-4*P(51)+9*P(273)
       P(740) = -p4sq-p5sq+s45
       FRe(375) = (C0134R-C0345R+p2sq*(Dij1234R(2,1)-2*Dij1235R(2,1)-D
     &   ij2345R(1,1))+p2sq*EijR(2,2)*P(1)+Dij2345R(3,1)*P(144)+Dij13
     &   45R(1,1)*P(190)+Dij1345R(2,1)*P(230)+Dij2345R(2,1)*P(264)+Ei
     &   jR(3,2)*P(233)*P(264)+P(21)*(p2sq*EijR(5,2)+EijR(7,2)*P(3)-E
     &   ijR(6,2)*P(264))+Dij1234R(3,1)*P(368)+Dij1345R(3,1)*P(385)+E
     &   ijR(4,2)*P(3)*P(385)+D01345R*P(540)-(p2sq*EijR(2,1)+EijR(4,1
     &   )*P(3)-EijR(3,1)*P(264))*P(541)+EijR(8,2)*P(542)+EijR(9,2)*P
     &   (543)+EijR(10,2)*P(544)-2*(C0135R+Dij1234R(7,2)-2*Dij1235R(7
     &   ,2)+Dij1245R(7,2)+Dij1345R(7,2)-Dij2345R(7,2)-EijR(23,3)*P(2
     &   33)-EijR(21,3)*P(386)-Dij1235R(3,1)*P(392)-EijR(22,3)*P(721)
     &   -EijR(11,2)*P(739)-EijR(24,3)*P(740)))/2.
       FIm(375) = (C0134I-C0345I+p2sq*(Dij1234I(2,1)-2*Dij1235I(2,1)-D
     &   ij2345I(1,1))+p2sq*EijI(2,2)*P(1)+Dij2345I(3,1)*P(144)+Dij13
     &   45I(1,1)*P(190)+Dij1345I(2,1)*P(230)+Dij2345I(2,1)*P(264)+Ei
     &   jI(3,2)*P(233)*P(264)+P(21)*(p2sq*EijI(5,2)+EijI(7,2)*P(3)-E
     &   ijI(6,2)*P(264))+Dij1234I(3,1)*P(368)+Dij1345I(3,1)*P(385)+E
     &   ijI(4,2)*P(3)*P(385)+D01345I*P(540)-(p2sq*EijI(2,1)+EijI(4,1
     &   )*P(3)-EijI(3,1)*P(264))*P(541)+EijI(8,2)*P(542)+EijI(9,2)*P
     &   (543)+EijI(10,2)*P(544)-2*(C0135I+Dij1234I(7,2)-2*Dij1235I(7
     &   ,2)+Dij1245I(7,2)+Dij1345I(7,2)-Dij2345I(7,2)-EijI(23,3)*P(2
     &   33)-EijI(21,3)*P(386)-Dij1235I(3,1)*P(392)-EijI(22,3)*P(721)
     &   -EijI(11,2)*P(739)-EijI(24,3)*P(740)))/2.
       F9(375)=DCMPLX(FRe(375),FIm(375))
       P(741) = 3*p4sq+2*p5sq-5*s45
       P(742) = 2*p4sq+s15-4*s45
       P(743) = -2*p4sq-s15+4*s45
       P(744) = -3*p4sq-2*p5sq+s12+s15-s34+5*s45
       P(745) = -p4sq+s12-s34+s45
       FRe(376) = -Dij1234R(2,1)-Dij1234R(4,2)+2*(Dij1235R(2,1)-Dij123
     &   5R(3,1)+Dij1235R(3,2)+Dij1235R(4,2)-Dij1235R(5,2)-Dij1235R(6
     &   ,2))+Dij2345R(3,2)-Dij2345R(5,2)+8*EijR(11,2)+16*(EijR(21,3)
     &   -EijR(24,3))+EijR(10,3)*P(1)+EijR(7,3)*P(21)+(EijR(10,2)+Eij
     &   R(19,3)+EijR(20,3))*P(230)+(EijR(8,2)+EijR(16,3)+EijR(17,3))
     &   *P(233)+EijR(5,3)*P(386)+EijR(4,1)*P(540)+EijR(2,1)*P(541)+E
     &   ijR(14,3)*P(673)+EijR(15,3)*P(674)+(EijR(2,2)+EijR(8,3))*P(7
     &   21)+EijR(4,3)*P(740)+EijR(4,2)*P(741)+EijR(5,2)*P(742)+EijR(
     &   7,2)*P(743)+EijR(9,2)*P(744)+EijR(18,3)*P(745)
       FIm(376) = -Dij1234I(2,1)-Dij1234I(4,2)+2*(Dij1235I(2,1)-Dij123
     &   5I(3,1)+Dij1235I(3,2)+Dij1235I(4,2)-Dij1235I(5,2)-Dij1235I(6
     &   ,2))+Dij2345I(3,2)-Dij2345I(5,2)+8*EijI(11,2)+16*(EijI(21,3)
     &   -EijI(24,3))+EijI(10,3)*P(1)+EijI(7,3)*P(21)+(EijI(10,2)+Eij
     &   I(19,3)+EijI(20,3))*P(230)+(EijI(8,2)+EijI(16,3)+EijI(17,3))
     &   *P(233)+EijI(5,3)*P(386)+EijI(4,1)*P(540)+EijI(2,1)*P(541)+E
     &   ijI(14,3)*P(673)+EijI(15,3)*P(674)+(EijI(2,2)+EijI(8,3))*P(7
     &   21)+EijI(4,3)*P(740)+EijI(4,2)*P(741)+EijI(5,2)*P(742)+EijI(
     &   7,2)*P(743)+EijI(9,2)*P(744)+EijI(18,3)*P(745)
       F9(376)=DCMPLX(FRe(376),FIm(376))
       P(746) = 2*p4sq+p5sq-s12-s15+s34-4*s45
       P(747) = -5*p4sq-3*p5sq+s12+s15-s34+9*s45
       P(748) = -p4sq-p5sq+s45+2*P(1)
       P(749) = -s12-s15+s34+2*P(385)
       FRe(377) = -Dij1234R(2,1)-Dij1234R(2,2)-4*Dij1235R(6,2)+Dij2345
     &   R(1,1)+Dij2345R(1,2)-Dij2345R(3,1)+Dij2345R(3,2)+2*(Dij1235R
     &   (2,1)+Dij1235R(2,2)-Dij1235R(3,1)+Dij1235R(3,2)-Dij2345R(5,2
     &   ))+8*EijR(11,2)+16*(EijR(22,3)-EijR(24,3))+(EijR(7,2)+2*EijR
     &   (18,3))*P(21)+EijR(10,2)*P(230)+(EijR(8,2)+EijR(9,3)+EijR(16
     &   ,3)-2*EijR(20,3))*P(233)+(EijR(5,2)+EijR(8,3)+EijR(14,3))*P(
     &   386)+EijR(4,1)*P(540)+EijR(2,1)*P(541)+EijR(2,3)*P(721)+EijR
     &   (4,3)*P(740)+EijR(4,2)*P(741)+EijR(2,2)*P(746)+EijR(9,2)*P(7
     &   47)+EijR(10,3)*P(748)+EijR(15,3)*P(749)
       FIm(377) = -Dij1234I(2,1)-Dij1234I(2,2)-4*Dij1235I(6,2)+Dij2345
     &   I(1,1)+Dij2345I(1,2)-Dij2345I(3,1)+Dij2345I(3,2)+2*(Dij1235I
     &   (2,1)+Dij1235I(2,2)-Dij1235I(3,1)+Dij1235I(3,2)-Dij2345I(5,2
     &   ))+8*EijI(11,2)+16*(EijI(22,3)-EijI(24,3))+(EijI(7,2)+2*EijI
     &   (18,3))*P(21)+EijI(10,2)*P(230)+(EijI(8,2)+EijI(9,3)+EijI(16
     &   ,3)-2*EijI(20,3))*P(233)+(EijI(5,2)+EijI(8,3)+EijI(14,3))*P(
     &   386)+EijI(4,1)*P(540)+EijI(2,1)*P(541)+EijI(2,3)*P(721)+EijI
     &   (4,3)*P(740)+EijI(4,2)*P(741)+EijI(2,2)*P(746)+EijI(9,2)*P(7
     &   47)+EijI(10,3)*P(748)+EijI(15,3)*P(749)
       F9(377)=DCMPLX(FRe(377),FIm(377))
       P(750) = -2*p4sq-p5sq+s12+s15-s34+4*s45
       P(751) = -4*p4sq-2*p5sq-s12+s34+6*s45
       P(752) = -p5sq+s15-2*P(64)
       FRe(378) = -Dij1234R(3,1)-Dij1234R(6,2)+Dij2345R(2,1)-Dij2345R(
     &   3,1)+Dij2345R(3,2)+Dij2345R(4,2)-Dij2345R(5,2)-Dij2345R(6,2)
     &   +8*EijR(11,2)+16*(EijR(23,3)-EijR(24,3))+EijR(10,3)*P(1)+(Ei
     &   jR(7,2)+EijR(18,3)+EijR(19,3))*P(21)+EijR(13,3)*P(230)+(EijR
     &   (3,2)+EijR(12,3))*P(233)+(EijR(6,2)+EijR(14,3)+EijR(17,3))*P
     &   (386)+EijR(4,1)*P(540)+EijR(3,1)*P(541)+EijR(15,3)*P(674)+Ei
     &   jR(16,3)*P(682)+EijR(9,3)*P(721)+EijR(4,3)*P(740)+EijR(4,2)*
     &   P(741)+EijR(8,2)*P(746)+EijR(9,2)*P(750)+EijR(10,2)*P(751)+E
     &   ijR(20,3)*P(752)
       FIm(378) = -Dij1234I(3,1)-Dij1234I(6,2)+Dij2345I(2,1)-Dij2345I(
     &   3,1)+Dij2345I(3,2)+Dij2345I(4,2)-Dij2345I(5,2)-Dij2345I(6,2)
     &   +8*EijI(11,2)+16*(EijI(23,3)-EijI(24,3))+EijI(10,3)*P(1)+(Ei
     &   jI(7,2)+EijI(18,3)+EijI(19,3))*P(21)+EijI(13,3)*P(230)+(EijI
     &   (3,2)+EijI(12,3))*P(233)+(EijI(6,2)+EijI(14,3)+EijI(17,3))*P
     &   (386)+EijI(4,1)*P(540)+EijI(3,1)*P(541)+EijI(15,3)*P(674)+Ei
     &   jI(16,3)*P(682)+EijI(9,3)*P(721)+EijI(4,3)*P(740)+EijI(4,2)*
     &   P(741)+EijI(8,2)*P(746)+EijI(9,2)*P(750)+EijI(10,2)*P(751)+E
     &   ijI(20,3)*P(752)
       F9(378)=DCMPLX(FRe(378),FIm(378))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
