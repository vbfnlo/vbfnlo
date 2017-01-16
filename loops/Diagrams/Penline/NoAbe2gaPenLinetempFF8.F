       subroutine NoAbe2gaPenLineFFhl8(F8)
       IMPLICIT NONE
       Real*8 P(767),FRe(295:336),FIm(295:336)
       Complex*16 F8(295:336)
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
       FRe(295) = -4*(EijR(3,1)-EijR(4,1)+4*EijR(4,4)+EijR(8,2)-EijR(9
     &   ,2)+2*(EijR(4,2)-EijR(4,3)-EijR(10,2)+EijR(15,3)+EijR(16,3)-
     &   EijR(20,3)-2*(EijR(15,4)+EijR(16,4)-EijR(34,4))))
       FIm(295) = -4*(EijI(3,1)-EijI(4,1)+4*EijI(4,4)+EijI(8,2)-EijI(9
     &   ,2)+2*(EijI(4,2)-EijI(4,3)-EijI(10,2)+EijI(15,3)+EijI(16,3)-
     &   EijI(20,3)-2*(EijI(15,4)+EijI(16,4)-EijI(34,4))))
       F8(295)=DCMPLX(FRe(295),FIm(295))
       FRe(296) = -8*(EijR(4,3)+EijR(8,2)-EijR(9,2)-EijR(16,3)-3*(EijR
     &   (15,3)-EijR(20,3))+2*(EijR(4,4)-EijR(15,4)-EijR(16,4)+EijR(3
     &   4,4)))
       FIm(296) = -8*(EijI(4,3)+EijI(8,2)-EijI(9,2)-EijI(16,3)-3*(EijI
     &   (15,3)-EijI(20,3))+2*(EijI(4,4)-EijI(15,4)-EijI(16,4)+EijI(3
     &   4,4)))
       F8(296)=DCMPLX(FRe(296),FIm(296))
       FRe(297) = -4*(EijR(2,1)-EijR(4,1)-3*EijR(5,2)+EijR(7,2)-6*EijR
     &   (9,2)-10*(EijR(4,3)-EijR(14,3)-EijR(15,3)+EijR(18,3))+4*(Eij
     &   R(2,2)+EijR(4,2)+EijR(4,4)-EijR(14,4)+EijR(21,4)-EijR(27,4))
     &   -8*(EijR(15,4)-EijR(32,4)))
       FIm(297) = -4*(EijI(2,1)-EijI(4,1)-3*EijI(5,2)+EijI(7,2)-6*EijI
     &   (9,2)-10*(EijI(4,3)-EijI(14,3)-EijI(15,3)+EijI(18,3))+4*(Eij
     &   I(2,2)+EijI(4,2)+EijI(4,4)-EijI(14,4)+EijI(21,4)-EijI(27,4))
     &   -8*(EijI(15,4)-EijI(32,4)))
       F8(297)=DCMPLX(FRe(297),FIm(297))
       FRe(298) = -4*(EijR(2,1)-EijR(4,1)-6*EijR(4,3)+4*EijR(4,4)+EijR
     &   (5,2)+EijR(7,2)-2*(EijR(9,2)-3*EijR(14,3)-EijR(15,3)+EijR(18
     &   ,3)+2*(EijR(14,4)+EijR(15,4)-EijR(32,4))))
       FIm(298) = -4*(EijI(2,1)-EijI(4,1)-6*EijI(4,3)+4*EijI(4,4)+EijI
     &   (5,2)+EijI(7,2)-2*(EijI(9,2)-3*EijI(14,3)-EijI(15,3)+EijI(18
     &   ,3)+2*(EijI(14,4)+EijI(15,4)-EijI(32,4))))
       F8(298)=DCMPLX(FRe(298),FIm(298))
       FRe(299) = -4*(EijR(2,1)-EijR(4,1)+12*EijR(4,2)-10*EijR(4,3)+4*
     &   EijR(4,4)+EijR(5,2)-3*EijR(7,2)+2*(-5*(EijR(9,2)-EijR(14,3))
     &   +6*EijR(16,3)+4*(EijR(8,2)-EijR(10,2)+EijR(17,3))+3*(EijR(15
     &   ,3)-EijR(18,3))-2*(EijR(14,4)+EijR(15,4)+EijR(16,4)+3*EijR(1
     &   9,3)+2*EijR(20,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4
     &   ))))
       FIm(299) = -4*(EijI(2,1)-EijI(4,1)+12*EijI(4,2)-10*EijI(4,3)+4*
     &   EijI(4,4)+EijI(5,2)-3*EijI(7,2)+2*(-5*(EijI(9,2)-EijI(14,3))
     &   +6*EijI(16,3)+4*(EijI(8,2)-EijI(10,2)+EijI(17,3))+3*(EijI(15
     &   ,3)-EijI(18,3))-2*(EijI(14,4)+EijI(15,4)+EijI(16,4)+3*EijI(1
     &   9,3)+2*EijI(20,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4
     &   ))))
       F8(299)=DCMPLX(FRe(299),FIm(299))
       FRe(300) = 8*(-EijR(4,3)+3*(EijR(4,2)-EijR(10,2))+EijR(16,3)+4*
     &   (EijR(15,3)-EijR(20,3))-2*(EijR(4,4)-EijR(15,4)-EijR(16,4)+E
     &   ijR(34,4)))
       FIm(300) = 8*(-EijI(4,3)+3*(EijI(4,2)-EijI(10,2))+EijI(16,3)+4*
     &   (EijI(15,3)-EijI(20,3))-2*(EijI(4,4)-EijI(15,4)-EijI(16,4)+E
     &   ijI(34,4)))
       F8(300)=DCMPLX(FRe(300),FIm(300))
       FRe(301) = 8*(EijR(3,1)-EijR(4,1)+EijR(4,3)+EijR(6,2)-EijR(7,2)
     &   +EijR(14,3)-EijR(16,3)+4*(EijR(4,2)+EijR(8,2)-EijR(9,2)-EijR
     &   (10,2)+EijR(17,3)-EijR(18,3))-EijR(19,3)-2*(EijR(4,4)-EijR(1
     &   4,4)-EijR(15,3)-EijR(15,4)-EijR(16,4)+EijR(20,3)+EijR(32,4)+
     &   EijR(33,4)+EijR(34,4)-EijR(35,4)))
       FIm(301) = 8*(EijI(3,1)-EijI(4,1)+EijI(4,3)+EijI(6,2)-EijI(7,2)
     &   +EijI(14,3)-EijI(16,3)+4*(EijI(4,2)+EijI(8,2)-EijI(9,2)-EijI
     &   (10,2)+EijI(17,3)-EijI(18,3))-EijI(19,3)-2*(EijI(4,4)-EijI(1
     &   4,4)-EijI(15,3)-EijI(15,4)-EijI(16,4)+EijI(20,3)+EijI(32,4)+
     &   EijI(33,4)+EijI(34,4)-EijI(35,4)))
       F8(301)=DCMPLX(FRe(301),FIm(301))
       FRe(302) = -8*(-EijR(4,3)+EijR(8,2)-EijR(9,2)-EijR(15,3)+3*EijR
     &   (16,3)+EijR(20,3)+2*(EijR(4,2)+EijR(4,4)-EijR(10,2)-EijR(13,
     &   3)-EijR(15,4)+EijR(22,4)-EijR(31,4))-4*(EijR(16,4)-EijR(34,4
     &   )))
       FIm(302) = -8*(-EijI(4,3)+EijI(8,2)-EijI(9,2)-EijI(15,3)+3*EijI
     &   (16,3)+EijI(20,3)+2*(EijI(4,2)+EijI(4,4)-EijI(10,2)-EijI(13,
     &   3)-EijI(15,4)+EijI(22,4)-EijI(31,4))-4*(EijI(16,4)-EijI(34,4
     &   )))
       F8(302)=DCMPLX(FRe(302),FIm(302))
       FRe(303) = -4*(EijR(4,1)+EijR(8,2)+EijR(10,2)+6*EijR(16,3)+2*(E
     &   ijR(4,2)-EijR(9,2)-EijR(20,3))-4*(EijR(4,3)-EijR(4,4)+EijR(1
     &   5,4)+EijR(16,4)-EijR(34,4)))
       FIm(303) = -4*(EijI(4,1)+EijI(8,2)+EijI(10,2)+6*EijI(16,3)+2*(E
     &   ijI(4,2)-EijI(9,2)-EijI(20,3))-4*(EijI(4,3)-EijI(4,4)+EijI(1
     &   5,4)+EijI(16,4)-EijI(34,4)))
       F8(303)=DCMPLX(FRe(303),FIm(303))
       FRe(304) = 4*(-2*EijR(7,2)+3*(EijR(2,1)-EijR(4,1)+EijR(5,2)+Eij
     &   R(9,2))+10*(EijR(4,3)-EijR(14,3)-EijR(15,3)+EijR(18,3))-4*(E
     &   ijR(4,2)+EijR(4,4)-EijR(15,4)+EijR(19,4)-EijR(24,4))+8*(EijR
     &   (14,4)-EijR(32,4)))
       FIm(304) = 4*(-2*EijI(7,2)+3*(EijI(2,1)-EijI(4,1)+EijI(5,2)+Eij
     &   I(9,2))+10*(EijI(4,3)-EijI(14,3)-EijI(15,3)+EijI(18,3))-4*(E
     &   ijI(4,2)+EijI(4,4)-EijI(15,4)+EijI(19,4)-EijI(24,4))+8*(EijI
     &   (14,4)-EijI(32,4)))
       F8(304)=DCMPLX(FRe(304),FIm(304))
       FRe(305) = -4*(EijR(3,1)-EijR(4,1)-6*EijR(4,3)+4*(EijR(4,4)+Eij
     &   R(6,2)-EijR(7,2))-3*(EijR(8,2)-EijR(9,2))+2*(EijR(4,2)-EijR(
     &   10,2)+3*(EijR(15,3)+EijR(16,3)-EijR(20,3))-2*(EijR(14,4)+Eij
     &   R(15,4)+EijR(16,4)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,
     &   4))))
       FIm(305) = -4*(EijI(3,1)-EijI(4,1)-6*EijI(4,3)+4*(EijI(4,4)+Eij
     &   I(6,2)-EijI(7,2))-3*(EijI(8,2)-EijI(9,2))+2*(EijI(4,2)-EijI(
     &   10,2)+3*(EijI(15,3)+EijI(16,3)-EijI(20,3))-2*(EijI(14,4)+Eij
     &   I(15,4)+EijI(16,4)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,
     &   4))))
       F8(305)=DCMPLX(FRe(305),FIm(305))
       FRe(306) = 4*(EijR(4,1)-2*EijR(4,2)-3*EijR(9,2)+6*(EijR(4,3)-Ei
     &   jR(15,3))-4*(EijR(4,4)-EijR(7,2)-EijR(14,4)-EijR(15,4)+EijR(
     &   32,4)))
       FIm(306) = 4*(EijI(4,1)-2*EijI(4,2)-3*EijI(9,2)+6*(EijI(4,3)-Ei
     &   jI(15,3))-4*(EijI(4,4)-EijI(7,2)-EijI(14,4)-EijI(15,4)+EijI(
     &   32,4)))
       F8(306)=DCMPLX(FRe(306),FIm(306))
       FRe(307) = -4*(EijR(3,1)-EijR(4,1)-6*EijR(4,3)-3*(EijR(6,2)-Eij
     &   R(7,2))+4*(EijR(4,4)+EijR(8,2))+2*(EijR(4,2)-EijR(10,2)+3*(E
     &   ijR(14,3)+EijR(16,3)-EijR(19,3))-2*(EijR(9,2)+EijR(14,4)+Eij
     &   R(15,4)+EijR(16,4)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,
     &   4))))
       FIm(307) = -4*(EijI(3,1)-EijI(4,1)-6*EijI(4,3)-3*(EijI(6,2)-Eij
     &   I(7,2))+4*(EijI(4,4)+EijI(8,2))+2*(EijI(4,2)-EijI(10,2)+3*(E
     &   ijI(14,3)+EijI(16,3)-EijI(19,3))-2*(EijI(9,2)+EijI(14,4)+Eij
     &   I(15,4)+EijI(16,4)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,
     &   4))))
       F8(307)=DCMPLX(FRe(307),FIm(307))
       FRe(308) = -4*(EE0R+EijR(1,2)-7*EijR(4,1)+8*EijR(4,2)-10*EijR(4
     &   ,3)+4*(EijR(2,1)+EijR(4,4)+EijR(5,2))-9*EijR(7,2)+2*(EijR(1,
     &   1)-3*EijR(7,3)-2*(EijR(9,2)-4*EijR(14,3)-EijR(15,3)+EijR(15,
     &   4)+EijR(18,3)-EijR(19,4)+EijR(24,4)+2*(EijR(14,4)-EijR(32,4)
     &   ))))
       FIm(308) = -4*(EE0I+EijI(1,2)-7*EijI(4,1)+8*EijI(4,2)-10*EijI(4
     &   ,3)+4*(EijI(2,1)+EijI(4,4)+EijI(5,2))-9*EijI(7,2)+2*(EijI(1,
     &   1)-3*EijI(7,3)-2*(EijI(9,2)-4*EijI(14,3)-EijI(15,3)+EijI(15,
     &   4)+EijI(18,3)-EijI(19,4)+EijI(24,4)+2*(EijI(14,4)-EijI(32,4)
     &   ))))
       F8(308)=DCMPLX(FRe(308),FIm(308))
       FRe(309) = 4*(EijR(4,1)-2*EijR(4,2)-3*EijR(7,2)+6*(EijR(4,3)-Ei
     &   jR(14,3))-4*(EijR(4,4)-EijR(9,2)-EijR(14,4)-EijR(15,4)+EijR(
     &   32,4)))
       FIm(309) = 4*(EijI(4,1)-2*EijI(4,2)-3*EijI(7,2)+6*(EijI(4,3)-Ei
     &   jI(14,3))-4*(EijI(4,4)-EijI(9,2)-EijI(14,4)-EijI(15,4)+EijI(
     &   32,4)))
       F8(309)=DCMPLX(FRe(309),FIm(309))
       FRe(310) = -4*(EijR(3,1)+EijR(6,2)-3*EijR(10,2)+6*EijR(16,3)-2*
     &   (EijR(4,1)-EijR(4,2)+EijR(7,2)+EijR(19,3))-4*(EijR(4,3)-EijR
     &   (4,4)+EijR(14,4)+EijR(16,4)-EijR(33,4)))
       FIm(310) = -4*(EijI(3,1)+EijI(6,2)-3*EijI(10,2)+6*EijI(16,3)-2*
     &   (EijI(4,1)-EijI(4,2)+EijI(7,2)+EijI(19,3))-4*(EijI(4,3)-EijI
     &   (4,4)+EijI(14,4)+EijI(16,4)-EijI(33,4)))
       F8(310)=DCMPLX(FRe(310),FIm(310))
       FRe(311) = -4*(-EijR(4,1)-5*EijR(4,2)+EijR(7,2)+2*(EijR(2,1)+Ei
     &   jR(5,2)-EijR(15,3))+6*EijR(18,3)-4*(EijR(4,3)-EijR(4,4)-EijR
     &   (9,2)+EijR(14,4)+EijR(15,4)-EijR(32,4)))
       FIm(311) = -4*(-EijI(4,1)-5*EijI(4,2)+EijI(7,2)+2*(EijI(2,1)+Ei
     &   jI(5,2)-EijI(15,3))+6*EijI(18,3)-4*(EijI(4,3)-EijI(4,4)-EijI
     &   (9,2)+EijI(14,4)+EijI(15,4)-EijI(32,4)))
       F8(311)=DCMPLX(FRe(311),FIm(311))
       FRe(312) = 4*(EijR(4,1)+EijR(9,2)-2*(EijR(4,2)-EijR(4,3)+EijR(1
     &   5,3))-4*(EijR(4,4)-EijR(15,4)))
       FIm(312) = 4*(EijI(4,1)+EijI(9,2)-2*(EijI(4,2)-EijI(4,3)+EijI(1
     &   5,3))-4*(EijI(4,4)-EijI(15,4)))
       F8(312)=DCMPLX(FRe(312),FIm(312))
       FRe(313) = -4*(EijR(2,1)+EijR(2,2)-EijR(4,1)-5*EijR(9,2)-10*(Ei
     &   jR(4,3)+EijR(10,3))+4*(EijR(4,2)+EijR(4,4)-EijR(10,4))+20*Ei
     &   jR(15,3)-12*(EijR(15,4)-EijR(21,4)))
       FIm(313) = -4*(EijI(2,1)+EijI(2,2)-EijI(4,1)-5*EijI(9,2)-10*(Ei
     &   jI(4,3)+EijI(10,3))+4*(EijI(4,2)+EijI(4,4)-EijI(10,4))+20*Ei
     &   jI(15,3)-12*(EijI(15,4)-EijI(21,4)))
       F8(313)=DCMPLX(FRe(313),FIm(313))
       FRe(314) = 8*(EijR(4,3)-EijR(16,3)+3*(EijR(8,2)-EijR(9,2)+EijR(
     &   15,3)-EijR(20,3))+2*(EijR(4,2)-EijR(4,4)-EijR(10,2)+EijR(16,
     &   4)-EijR(21,4)+EijR(28,4))+4*(EijR(9,3)-EijR(10,3)+EijR(15,4)
     &   -EijR(34,4)))
       FIm(314) = 8*(EijI(4,3)-EijI(16,3)+3*(EijI(8,2)-EijI(9,2)+EijI(
     &   15,3)-EijI(20,3))+2*(EijI(4,2)-EijI(4,4)-EijI(10,2)+EijI(16,
     &   4)-EijI(21,4)+EijI(28,4))+4*(EijI(9,3)-EijI(10,3)+EijI(15,4)
     &   -EijI(34,4)))
       F8(314)=DCMPLX(FRe(314),FIm(314))
       FRe(315) = -4*(-EijR(4,1)-5*(EijR(4,2)-EijR(9,2))+6*EijR(10,3)+
     &   2*(EijR(2,1)+EijR(2,2)-EijR(15,3))-8*EijR(15,4)-4*(EijR(4,3)
     &   -EijR(4,4)-EijR(21,4)))
       FIm(315) = -4*(-EijI(4,1)-5*(EijI(4,2)-EijI(9,2))+6*EijI(10,3)+
     &   2*(EijI(2,1)+EijI(2,2)-EijI(15,3))-8*EijI(15,4)-4*(EijI(4,3)
     &   -EijI(4,4)-EijI(21,4)))
       F8(315)=DCMPLX(FRe(315),FIm(315))
       FRe(316) = -4*(EijR(2,1)-3*EijR(2,2)-EijR(4,1)-10*EijR(4,3)-Eij
     &   R(9,2)-6*EijR(10,3)+16*EijR(15,3)+4*(EijR(4,2)+EijR(4,4)+Eij
     &   R(5,2)-EijR(7,2)+EijR(14,3)-EijR(14,4)-EijR(18,3)+EijR(21,4)
     &   -EijR(27,4))-8*(EijR(15,4)-EijR(32,4)))
       FIm(316) = -4*(EijI(2,1)-3*EijI(2,2)-EijI(4,1)-10*EijI(4,3)-Eij
     &   I(9,2)-6*EijI(10,3)+16*EijI(15,3)+4*(EijI(4,2)+EijI(4,4)+Eij
     &   I(5,2)-EijI(7,2)+EijI(14,3)-EijI(14,4)-EijI(18,3)+EijI(21,4)
     &   -EijI(27,4))-8*(EijI(15,4)-EijI(32,4)))
       F8(316)=DCMPLX(FRe(316),FIm(316))
       FRe(317) = -4*(EijR(2,1)-EijR(4,1)-10*EijR(4,3)+4*(EijR(4,2)+Ei
     &   jR(4,4))+EijR(5,2)-3*EijR(7,2)-2*(EijR(9,2)-7*EijR(14,3)-3*(
     &   EijR(15,3)-EijR(18,3))+2*(EijR(7,3)+EijR(15,4)-EijR(19,4)+Ei
     &   jR(24,4))+4*(EijR(14,4)-EijR(32,4))))
       FIm(317) = -4*(EijI(2,1)-EijI(4,1)-10*EijI(4,3)+4*(EijI(4,2)+Ei
     &   jI(4,4))+EijI(5,2)-3*EijI(7,2)-2*(EijI(9,2)-7*EijI(14,3)-3*(
     &   EijI(15,3)-EijI(18,3))+2*(EijI(7,3)+EijI(15,4)-EijI(19,4)+Ei
     &   jI(24,4))+4*(EijI(14,4)-EijI(32,4))))
       F8(317)=DCMPLX(FRe(317),FIm(317))
       FRe(318) = 8*(EijR(4,3)+3*(EijR(8,2)-EijR(9,2))+EijR(15,3)-EijR
     &   (16,3)+4*(EijR(17,3)-EijR(18,3))-EijR(20,3)+2*(EijR(4,2)-Eij
     &   R(4,4)-EijR(10,2)+EijR(14,3)+EijR(14,4)+EijR(15,4)+EijR(16,4
     &   )-EijR(19,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4)))
       FIm(318) = 8*(EijI(4,3)+3*(EijI(8,2)-EijI(9,2))+EijI(15,3)-EijI
     &   (16,3)+4*(EijI(17,3)-EijI(18,3))-EijI(20,3)+2*(EijI(4,2)-Eij
     &   I(4,4)-EijI(10,2)+EijI(14,3)+EijI(14,4)+EijI(15,4)+EijI(16,4
     &   )-EijI(19,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4)))
       F8(318)=DCMPLX(FRe(318),FIm(318))
       FRe(319) = -4*(-EijR(4,1)-5*(EijR(4,2)-EijR(9,2))+2*(EijR(2,1)+
     &   EijR(5,2)-EijR(14,3))+6*EijR(18,3)-4*(EijR(4,3)-EijR(4,4)+Ei
     &   jR(14,4)+EijR(15,4)-EijR(32,4)))
       FIm(319) = -4*(-EijI(4,1)-5*(EijI(4,2)-EijI(9,2))+2*(EijI(2,1)+
     &   EijI(5,2)-EijI(14,3))+6*EijI(18,3)-4*(EijI(4,3)-EijI(4,4)+Ei
     &   jI(14,4)+EijI(15,4)-EijI(32,4)))
       F8(319)=DCMPLX(FRe(319),FIm(319))
       FRe(320) = -4*(EE0R+2*EijR(1,1)+EijR(1,2)-3*EijR(4,1)+8*EijR(4,
     &   2)-9*EijR(7,2)-10*(EijR(4,3)+EijR(7,3))+4*(EijR(4,4)-EijR(7,
     &   4))+20*EijR(14,3)-12*(EijR(14,4)-EijR(19,4)))
       FIm(320) = -4*(EE0I+2*EijI(1,1)+EijI(1,2)-3*EijI(4,1)+8*EijI(4,
     &   2)-9*EijI(7,2)-10*(EijI(4,3)+EijI(7,3))+4*(EijI(4,4)-EijI(7,
     &   4))+20*EijI(14,3)-12*(EijI(14,4)-EijI(19,4)))
       F8(320)=DCMPLX(FRe(320),FIm(320))
       FRe(321) = 4*(EijR(4,1)-2*EijR(4,2)+EijR(7,2)+6*(EijR(4,3)-EijR
     &   (14,3))+8*EijR(14,4)-4*(EijR(4,4)+EijR(19,4)))
       FIm(321) = 4*(EijI(4,1)-2*EijI(4,2)+EijI(7,2)+6*(EijI(4,3)-EijI
     &   (14,3))+8*EijI(14,4)-4*(EijI(4,4)+EijI(19,4)))
       F8(321)=DCMPLX(FRe(321),FIm(321))
       FRe(322) = -4*(EE0R+EijR(1,2)-3*EijR(4,1)-6*EijR(4,3)-5*EijR(7,
     &   2)+2*(EijR(1,1)-EijR(7,3))+8*(EijR(14,3)-EijR(14,4))+4*(EijR
     &   (4,2)+EijR(4,4)+EijR(19,4)))
       FIm(322) = -4*(EE0I+EijI(1,2)-3*EijI(4,1)-6*EijI(4,3)-5*EijI(7,
     &   2)+2*(EijI(1,1)-EijI(7,3))+8*(EijI(14,3)-EijI(14,4))+4*(EijI
     &   (4,2)+EijI(4,4)+EijI(19,4)))
       F8(322)=DCMPLX(FRe(322),FIm(322))
       FRe(323) = 4*(EijR(4,1)+EijR(7,2)-2*(EijR(4,2)-EijR(4,3)+EijR(1
     &   4,3))-4*(EijR(4,4)-EijR(14,4)))
       FIm(323) = 4*(EijI(4,1)+EijI(7,2)-2*(EijI(4,2)-EijI(4,3)+EijI(1
     &   4,3))-4*(EijI(4,4)-EijI(14,4)))
       F8(323)=DCMPLX(FRe(323),FIm(323))
       FRe(324) = -4*(EE0R+EijR(1,2)-11*EijR(4,1)-10*EijR(4,3)+4*EijR(
     &   4,4)+16*(EijR(4,2)+EijR(6,2))+8*(EijR(3,1)+EijR(6,3))-17*Eij
     &   R(7,2)+2*(EijR(1,1)-3*EijR(7,3)-2*(4*(EijR(10,2)-EijR(14,3))
     &   -3*EijR(16,3)+EijR(16,4)+5*EijR(19,3)-EijR(19,4)+EijR(25,4)+
     &   2*(EijR(14,4)-EijR(33,4)))))
       FIm(324) = -4*(EE0I+EijI(1,2)-11*EijI(4,1)-10*EijI(4,3)+4*EijI(
     &   4,4)+16*(EijI(4,2)+EijI(6,2))+8*(EijI(3,1)+EijI(6,3))-17*Eij
     &   I(7,2)+2*(EijI(1,1)-3*EijI(7,3)-2*(4*(EijI(10,2)-EijI(14,3))
     &   -3*EijI(16,3)+EijI(16,4)+5*EijI(19,3)-EijI(19,4)+EijI(25,4)+
     &   2*(EijI(14,4)-EijI(33,4)))))
       F8(324)=DCMPLX(FRe(324),FIm(324))
       FRe(325) = 4*(EijR(4,1)-10*EijR(4,2)+EijR(7,2)+8*EijR(10,2)+6*(
     &   EijR(4,3)-EijR(14,3))-4*(EijR(4,4)-EijR(14,4)-EijR(16,4)+2*(
     &   EijR(16,3)-EijR(19,3))+EijR(33,4)))
       FIm(325) = 4*(EijI(4,1)-10*EijI(4,2)+EijI(7,2)+8*EijI(10,2)+6*(
     &   EijI(4,3)-EijI(14,3))-4*(EijI(4,4)-EijI(14,4)-EijI(16,4)+2*(
     &   EijI(16,3)-EijI(19,3))+EijI(33,4)))
       F8(325)=DCMPLX(FRe(325),FIm(325))
       FRe(326) = -4*(EijR(2,1)-EijR(4,1)-10*EijR(4,3)+EijR(5,2)-2*Eij
     &   R(7,2)-3*EijR(9,2)+14*EijR(15,3)+6*(EijR(14,3)-EijR(18,3))+4
     &   *(EijR(4,2)+EijR(4,4)-EijR(10,3)-EijR(14,4)+EijR(21,4)-EijR(
     &   27,4))-8*(EijR(15,4)-EijR(32,4)))
       FIm(326) = -4*(EijI(2,1)-EijI(4,1)-10*EijI(4,3)+EijI(5,2)-2*Eij
     &   I(7,2)-3*EijI(9,2)+14*EijI(15,3)+6*(EijI(14,3)-EijI(18,3))+4
     &   *(EijI(4,2)+EijI(4,4)-EijI(10,3)-EijI(14,4)+EijI(21,4)-EijI(
     &   27,4))-8*(EijI(15,4)-EijI(32,4)))
       F8(326)=DCMPLX(FRe(326),FIm(326))
       FRe(327) = -4*(EijR(3,1)-EijR(4,1)+EijR(8,2)-EijR(9,2)+2*(EijR(
     &   4,2)-EijR(10,2))-6*(EijR(4,3)-EijR(15,3)-EijR(16,3)+EijR(20,
     &   3))+4*(EijR(4,4)-EijR(16,4)+EijR(21,4)-EijR(28,4))-8*(EijR(1
     &   5,4)-EijR(34,4)))
       FIm(327) = -4*(EijI(3,1)-EijI(4,1)+EijI(8,2)-EijI(9,2)+2*(EijI(
     &   4,2)-EijI(10,2))-6*(EijI(4,3)-EijI(15,3)-EijI(16,3)+EijI(20,
     &   3))+4*(EijI(4,4)-EijI(16,4)+EijI(21,4)-EijI(28,4))-8*(EijI(1
     &   5,4)-EijI(34,4)))
       F8(327)=DCMPLX(FRe(327),FIm(327))
       FRe(328) = 4*(EijR(4,1)-2*EijR(4,2)+EijR(9,2)+6*(EijR(4,3)-EijR
     &   (15,3))+8*EijR(15,4)-4*(EijR(4,4)+EijR(21,4)))
       FIm(328) = 4*(EijI(4,1)-2*EijI(4,2)+EijI(9,2)+6*(EijI(4,3)-EijI
     &   (15,3))+8*EijI(15,4)-4*(EijI(4,4)+EijI(21,4)))
       F8(328)=DCMPLX(FRe(328),FIm(328))
       FRe(329) = -4*(EijR(2,1)-5*EijR(4,1)+12*EijR(4,2)-10*EijR(4,3)+
     &   EijR(5,2)+4*(EijR(3,1)+EijR(4,4)+EijR(6,2))-6*EijR(7,2)+8*Ei
     &   jR(8,2)-7*EijR(9,2)-2*(-5*EijR(15,3)-3*(EijR(14,3)-EijR(18,3
     &   ))-4*(EijR(17,3)-EijR(19,3))+6*(EijR(10,2)-EijR(16,3)+EijR(2
     &   0,3))+2*(EijR(14,4)+EijR(15,4)+EijR(16,4)-EijR(32,4)-EijR(33
     &   ,4)-EijR(34,4)+EijR(35,4))))
       FIm(329) = -4*(EijI(2,1)-5*EijI(4,1)+12*EijI(4,2)-10*EijI(4,3)+
     &   EijI(5,2)+4*(EijI(3,1)+EijI(4,4)+EijI(6,2))-6*EijI(7,2)+8*Ei
     &   jI(8,2)-7*EijI(9,2)-2*(-5*EijI(15,3)-3*(EijI(14,3)-EijI(18,3
     &   ))-4*(EijI(17,3)-EijI(19,3))+6*(EijI(10,2)-EijI(16,3)+EijI(2
     &   0,3))+2*(EijI(14,4)+EijI(15,4)+EijI(16,4)-EijI(32,4)-EijI(33
     &   ,4)-EijI(34,4)+EijI(35,4))))
       F8(329)=DCMPLX(FRe(329),FIm(329))
       FRe(330) = 4*(EijR(4,1)-10*EijR(4,2)+EijR(9,2)+8*EijR(10,2)+6*(
     &   EijR(4,3)-EijR(15,3))-4*(EijR(4,4)-EijR(15,4)-EijR(16,4)+2*(
     &   EijR(16,3)-EijR(20,3))+EijR(34,4)))
       FIm(330) = 4*(EijI(4,1)-10*EijI(4,2)+EijI(9,2)+8*EijI(10,2)+6*(
     &   EijI(4,3)-EijI(15,3))-4*(EijI(4,4)-EijI(15,4)-EijI(16,4)+2*(
     &   EijI(16,3)-EijI(20,3))+EijI(34,4)))
       F8(330)=DCMPLX(FRe(330),FIm(330))
       FRe(331) = 4*(5*EijR(4,2)+2*(EijR(4,1)+EijR(9,2))+6*EijR(15,3)-
     &   4*(EijR(4,4)-EijR(15,4)))
       FIm(331) = 4*(5*EijI(4,2)+2*(EijI(4,1)+EijI(9,2))+6*EijI(15,3)-
     &   4*(EijI(4,4)-EijI(15,4)))
       F8(331)=DCMPLX(FRe(331),FIm(331))
       FRe(332) = 4*(EijR(10,2)-2*(EijR(4,2)+EijR(16,3))-4*(EijR(4,4)-
     &   EijR(16,4)))
       FIm(332) = -4*(-EijI(10,2)+2*(EijI(4,2)+EijI(16,3))+4*(EijI(4,4
     &   )-EijI(16,4)))
       F8(332)=DCMPLX(FRe(332),FIm(332))
       FRe(333) = 4*(EijR(4,1)+3*EijR(4,2)+EijR(7,2)-4*(EijR(4,4)-EijR
     &   (14,3)-EijR(14,4)))
       FIm(333) = 4*(EijI(4,1)+3*EijI(4,2)+EijI(7,2)-4*(EijI(4,4)-EijI
     &   (14,3)-EijI(14,4)))
       F8(333)=DCMPLX(FRe(333),FIm(333))
       FRe(334) = 4*(EijR(4,1)+3*EijR(4,2)+EijR(9,2)-4*(EijR(4,4)-EijR
     &   (15,3)-EijR(15,4)))
       FIm(334) = 4*(EijI(4,1)+3*EijI(4,2)+EijI(9,2)-4*(EijI(4,4)-EijI
     &   (15,3)-EijI(15,4)))
       F8(334)=DCMPLX(FRe(334),FIm(334))
       FRe(335) = -4*(3*EijR(4,1)-5*EijR(4,2)+9*EijR(7,2)+6*EijR(7,3)+
     &   2*(EE0R+EijR(1,2)-EijR(14,3))-8*EijR(14,4)+4*(EijR(1,1)-EijR
     &   (4,3)+EijR(4,4)+EijR(19,4)))
       FIm(335) = -4*(3*EijI(4,1)-5*EijI(4,2)+9*EijI(7,2)+6*EijI(7,3)+
     &   2*(EE0I+EijI(1,2)-EijI(14,3))-8*EijI(14,4)+4*(EijI(1,1)-EijI
     &   (4,3)+EijI(4,4)+EijI(19,4)))
       F8(335)=DCMPLX(FRe(335),FIm(335))
       FRe(336) = 4*(5*EijR(4,2)+2*(EijR(4,1)+EijR(7,2))+6*EijR(14,3)-
     &   4*(EijR(4,4)-EijR(14,4)))
       FIm(336) = 4*(5*EijI(4,2)+2*(EijI(4,1)+EijI(7,2))+6*EijI(14,3)-
     &   4*(EijI(4,4)-EijI(14,4)))
       F8(336)=DCMPLX(FRe(336),FIm(336))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
