       subroutine HexBoxFFhl1(F1)
       Real*8 P(1362),FRe(1:65),FIm(1:65)
       Complex*16 F1(1:65)
       Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   p6sq 
       Real*8   s12, s16 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s56 
       Real*8   s123,s234,s345         
       Real*8  B012R,B013R,B014R,B015R,B016R 
       Real*8  B023R,B024R,B025R,B026R 
       Real*8  B034R,B035R,B036R 
       Real*8  B045R,B046R 
       Real*8  B056R  
       Real*8  B012I,B013I,B014I,B015I,B016I 
       Real*8  B023I,B024I,B025I,B026I 
       Real*8  B034I,B035I,B036I 
       Real*8  B045I,B046I 
       Real*8  B056I          
       Real*8 C0123R,C0124R,C0125R,C0126R 
       Real*8 C0134R,C0135R,C0136R 
       Real*8 C0145R,C0146R 
       Real*8 C0156R 
       Real*8 C0234R,C0235R,C0236R 
       Real*8 C0245R,C0246R 
       Real*8 C0256R 
       Real*8 C0345R,C0346R 
       Real*8 C0356R 
       Real*8 C0456R    
       Real*8 C0123I,C0124I,C0125I,C0126I 
       Real*8 C0134I,C0135I,C0136I 
       Real*8 C0145I,C0146I 
       Real*8 C0156I 
       Real*8 C0234I,C0235I,C0236I 
       Real*8 C0245I,C0246I 
       Real*8 C0256I 
       Real*8 C0345I,C0346I 
       Real*8 C0356I 
       Real*8 C0456I      
       Real*8 C123R(4,2),C124R(4,2),C125R(4,2),C126R(4,2) 
       Real*8 C134R(4,2),C135R(4,2),C136R(4,2) 
       Real*8 C145R(4,2),C146R(4,2) 
       Real*8 C156R(4,2) 
       Real*8 C234R(4,2),C235R(4,2),C236R(4,2) 
       Real*8 C245R(4,2),C246R(4,2) 
       Real*8 C256R(4,2) 
       Real*8 C345R(4,2),C346R(4,2) 
       Real*8 C356R(4,2) 
       Real*8 C456R(4,2)  
       Real*8 C123I(4,2),C124I(4,2),C125I(4,2),C126I(4,2) 
       Real*8 C134I(4,2),C135I(4,2),C136I(4,2) 
       Real*8 C145I(4,2),C146I(4,2) 
       Real*8 C156I(4,2) 
       Real*8 C234I(4,2),C235I(4,2),C236I(4,2) 
       Real*8 C245I(4,2),C246I(4,2) 
       Real*8 C256I(4,2) 
       Real*8 C345I(4,2),C346I(4,2) 
       Real*8 C356I(4,2) 
       Real*8 C456I(4,2)    
       Real*8 D01234R,D01235R,D01236R 
       Real*8 D01245R,D01246R 
       Real*8 D01256R 
       Real*8 D01345R,D01346R 
       Real*8 D01356R 
       Real*8 D01456R 
       Real*8 D02345R,D02346R 
       Real*8 D02356R 
       Real*8 D02456R 
       Real*8 D03456R  
       Real*8 D01234I,D01235I,D01236I 
       Real*8 D01245I,D01246I 
       Real*8 D01256I 
       Real*8 D01345I,D01346I 
       Real*8 D01356I 
       Real*8 D01456I 
       Real*8 D02345I,D02346I 
       Real*8 D02356I 
       Real*8 D02456I 
       Real*8 D03456I   
       Real*8 D1234R(13,3),D1235R(13,3),D1236R(13,3) 
       Real*8 D1245R(13,3),D1246R(13,3) 
       Real*8 D1256R(13,3) 
       Real*8 D1345R(13,3),D1346R(13,3) 
       Real*8 D1356R(13,3) 
       Real*8 D1456R(13,3) 
       Real*8 D2345R(13,3),D2346R(13,3) 
       Real*8 D2356R(13,3) 
       Real*8 D2456R(13,3) 
       Real*8 D3456R(13,3)   
       Real*8 D1234I(13,3),D1235I(13,3),D1236I(13,3) 
       Real*8 D1245I(13,3),D1246I(13,3) 
       Real*8 D1256I(13,3) 
       Real*8 D1345I(13,3),D1346I(13,3) 
       Real*8 D1356I(13,3) 
       Real*8 D1456I(13,3) 
       Real*8 D2345I(13,3),D2346I(13,3) 
       Real*8 D2356I(13,3) 
       Real*8 D2456I(13,3) 
       Real*8 D3456I(13,3)    
       Real*8 E012345R,E012346R,E012356R,E012456R,E013456R,E023456R 
       Real*8 E012345I,E012346I,E012356I,E012456I,E013456I,E023456I  
       Real*8 E12345R(46,4),E12346R(46,4),E12356R(46,4) 
       Real*8 E12456R(46,4),E13456R(46,4),E23456R(46,4)  
       Real*8 E12345I(46,4),E12346I(46,4),E12356I(46,4) 
       Real*8 E12456I(46,4),E13456I(46,4),E23456I(46,4)   
       Real*8 F0R 
       Real*8 F0I 
       Real*8 FR(166,5) 
       Real*8 FI(166,5) 

       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
       COMMON/HexBoxPFunctions/P
       COMMON/HexBoxFVALFunctions/F0R,FR,F0I,FI
       COMMON/HexBoxEVALFunctions/ E012345R,E12345R,E012345I,E12345I, 
     &   E012346R,E12346R,E012346I,E12346I, E012356R,E12356R,E012356I
     &   ,E12356I, E012456R,E12456R,E012456I,E12456I, E013456R,E13456
     &   R,E013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/HexBoxDVALFunctions/ D01234R,D1234R,D01234I,D1234I, D012
     &   35R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D1236I, D0
     &   1245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1246I, 
     &   D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,D1345I
     &   , D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356I,D135
     &   6I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D02345I,D2
     &   345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D02356I,
     &   D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,D03456
     &   I,D3456I
       COMMON/HexBoxCVALFunctions/ C0123R,C123R,C0123I,C123I, C0124R,C
     &   124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C126R,C
     &   0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C0135I,
     &   C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,C145I,
     &    C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I, C0234
     &   R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236R,C236
     &   R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246R,C024
     &   6I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C0345I,C34
     &   5I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C356I, C0
     &   456R,C456R,C0456I,C456I
       COMMON/HexBoxBVALFunctions/ B012R,B012I, B013R,B013I, B014R,B01
     &   4I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,B024I, B025
     &   R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,B036I,
     &    B045R,B045I, B046R,B046I, B056R,B056I
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       P(1) = s234-s56
       P(2) = s45+s56
       P(3) = s16*s234+s123*P(1)-P(1)*P(2)
       P(4) = s16+s56
       P(5) = s234+s56
       P(6) = s16*s234-s16*s45-s234*s45-2*s45*s56-s56**2+p4sq*P(4)+s12
     &   3*P(5)
       P(7) = -s123+2*s45+s56
       P(8) = s16*s45-p4sq*P(4)+s56*P(7)
       P(9) = s123+s16
       P(10) = p4sq-3*s45-s56+2*P(9)
       P(11) = s123+s234
       P(12) = p4sq-3*s45+2*P(11)
       P(13) = s123*s234+s16*s234+s56**2
       P(14) = p4sq*s16-p4sq*s234-s16*s45+3*s234*s45+s56*P(12)-2*P(13)
       P(15) = p4sq-s45
       P(16) = s16-s234+s56
       P(17) = s45-s56
       P(18) = p4sq*s234-s234*P(17)
       P(19) = -s123+s45+s56
       P(20) = -4*P(18)+8*s56*P(19)
       FRe(1) = -4*D1256R(1,1)*P(3)+4*D1256R(2,1)*P(6)+4*D1256R(3,1)*P
     &   (8)+4*s234*D1256R(2,2)*P(10)+4*D1256R(4,2)*P(14)+24*D1256R(7
     &   ,2)*P(15)-4*D1256R(5,2)*P(15)*P(16)+D1256R(6,2)*P(20)
       FIm(1) = -4*D1256I(1,1)*P(3)+4*D1256I(2,1)*P(6)+4*D1256I(3,1)*P
     &   (8)+4*s234*D1256I(2,2)*P(10)+4*D1256I(4,2)*P(14)+24*D1256I(7
     &   ,2)*P(15)-4*D1256I(5,2)*P(15)*P(16)+D1256I(6,2)*P(20)
       F1(1)=DCMPLX(FRe(1),FIm(1))
       P(21) = s123-2*s45-s56
       P(22) = -(s16*s45)+p4sq*P(4)+s56*P(21)
       P(23) = -(p4sq*s16)+p4sq*s234+s16*s45-3*s234*s45-s56*P(12)+2*P(
     &   13)
       P(24) = -s45+s56
       P(25) = s123*s56-s56*P(2)
       P(26) = p4sq*s234+s234*P(24)+2*P(25)
       FRe(2) = D1256R(1,1)*P(3)-D1256R(2,1)*P(6)-s234*D1256R(2,2)*P(1
     &   0)-6*D1256R(7,2)*P(15)+D1256R(5,2)*P(15)*P(16)+D1256R(3,1)*P
     &   (22)+D1256R(4,2)*P(23)+D1256R(6,2)*P(26)
       FIm(2) = D1256I(1,1)*P(3)-D1256I(2,1)*P(6)-s234*D1256I(2,2)*P(1
     &   0)-6*D1256I(7,2)*P(15)+D1256I(5,2)*P(15)*P(16)+D1256I(3,1)*P
     &   (22)+D1256I(4,2)*P(23)+D1256I(6,2)*P(26)
       F1(2)=DCMPLX(FRe(2),FIm(2))
       P(27) = s16-s234
       P(28) = s16+s234-s56
       P(29) = -2*s16+s234
       FRe(3) = -(s234*D1256R(2,2))-s16*D1256R(3,1)-2*D1256R(7,2)+D125
     &   6R(1,1)*P(1)+D1256R(5,2)*P(16)+D1256R(2,1)*P(27)+D1256R(4,2)
     &   *P(28)+D1256R(6,2)*P(29)
       FIm(3) = -(s234*D1256I(2,2))-s16*D1256I(3,1)-2*D1256I(7,2)+D125
     &   6I(1,1)*P(1)+D1256I(5,2)*P(16)+D1256I(2,1)*P(27)+D1256I(4,2)
     &   *P(28)+D1256I(6,2)*P(29)
       F1(3)=DCMPLX(FRe(3),FIm(3))
       P(30) = s123+s234-s45-s56
       P(31) = -s123+s234+s45+s56
       P(32) = s123-s45-s56
       P(33) = s234+2*P(32)
       FRe(4) = -(s234*D1256R(1,1))+s234*D1256R(2,2)+D1256R(3,1)*P(19)
     &   +D1256R(2,1)*P(30)+D1256R(5,2)*P(30)-D1256R(4,2)*P(31)-D1256
     &   R(6,2)*P(33)
       FIm(4) = -(s234*D1256I(1,1))+s234*D1256I(2,2)+D1256I(3,1)*P(19)
     &   +D1256I(2,1)*P(30)+D1256I(5,2)*P(30)-D1256I(4,2)*P(31)-D1256
     &   I(6,2)*P(33)
       F1(4)=DCMPLX(FRe(4),FIm(4))
       P(34) = 16*s16-8*s234
       FRe(5) = 8*s234*D1256R(2,2)+8*s16*D1256R(3,1)+16*D1256R(7,2)-8*
     &   D1256R(1,1)*P(1)-8*D1256R(5,2)*P(16)-8*D1256R(2,1)*P(27)-8*D
     &   1256R(4,2)*P(28)+D1256R(6,2)*P(34)
       FIm(5) = 8*s234*D1256I(2,2)+8*s16*D1256I(3,1)+16*D1256I(7,2)-8*
     &   D1256I(1,1)*P(1)-8*D1256I(5,2)*P(16)-8*D1256I(2,1)*P(27)-8*D
     &   1256I(4,2)*P(28)+D1256I(6,2)*P(34)
       F1(5)=DCMPLX(FRe(5),FIm(5))
       FRe(6) = s234*D1256R(1,1)-s234*D1256R(2,2)-D1256R(2,1)*P(30)-D1
     &   256R(5,2)*P(30)+D1256R(4,2)*P(31)+D1256R(3,1)*P(32)+D1256R(6
     &   ,2)*P(33)
       FIm(6) = s234*D1256I(1,1)-s234*D1256I(2,2)-D1256I(2,1)*P(30)-D1
     &   256I(5,2)*P(30)+D1256I(4,2)*P(31)+D1256I(3,1)*P(32)+D1256I(6
     &   ,2)*P(33)
       F1(6)=DCMPLX(FRe(6),FIm(6))
       FRe(7) = s234*D1256R(2,2)-s234*D1256R(6,2)+6*D1256R(7,2)+D1256R
     &   (2,1)*P(4)-D1256R(3,1)*P(4)+D1256R(4,2)*P(16)-D1256R(5,2)*P(
     &   16)
       FIm(7) = s234*D1256I(2,2)-s234*D1256I(6,2)+6*D1256I(7,2)+D1256I
     &   (2,1)*P(4)-D1256I(3,1)*P(4)+D1256I(4,2)*P(16)-D1256I(5,2)*P(
     &   16)
       F1(7)=DCMPLX(FRe(7),FIm(7))
       FRe(8) = D1256R(1,1)+D1256R(2,1)-D1256R(3,1)+2*D1256R(4,2)-2*D1
     &   256R(6,2)
       FIm(8) = D1256I(1,1)+D1256I(2,1)-D1256I(3,1)+2*D1256I(4,2)-2*D1
     &   256I(6,2)
       F1(8)=DCMPLX(FRe(8),FIm(8))
       FRe(9) = -2*s234*D1256R(2,2)+s56*D1256R(3,1)+2*s56*D1256R(6,2)+
     &   D1256R(1,1)*P(1)+2*D1256R(4,2)*P(1)-D1256R(2,1)*P(5)
       FIm(9) = -2*s234*D1256I(2,2)+s56*D1256I(3,1)+2*s56*D1256I(6,2)+
     &   D1256I(1,1)*P(1)+2*D1256I(4,2)*P(1)-D1256I(2,1)*P(5)
       F1(9)=DCMPLX(FRe(9),FIm(9))
       FRe(10) = D1256R(1,1)-D1256R(2,1)-2*D1256R(2,2)+2*D1256R(4,2)
       FIm(10) = D1256I(1,1)-D1256I(2,1)-2*D1256I(2,2)+2*D1256I(4,2)
       F1(10)=DCMPLX(FRe(10),FIm(10))
       FRe(11) = s234*D1256R(2,2)-s234*D1256R(6,2)+D1256R(2,1)*P(19)-D
     &   1256R(4,2)*P(30)+D1256R(5,2)*P(30)+D1256R(3,1)*P(32)
       FIm(11) = s234*D1256I(2,2)-s234*D1256I(6,2)+D1256I(2,1)*P(19)-D
     &   1256I(4,2)*P(30)+D1256I(5,2)*P(30)+D1256I(3,1)*P(32)
       F1(11)=DCMPLX(FRe(11),FIm(11))
       FRe(12) = D1256R(2,1)+D1256R(2,2)-D1256R(3,1)-D1256R(6,2)
       FIm(12) = D1256I(2,1)+D1256I(2,2)-D1256I(3,1)-D1256I(6,2)
       F1(12)=DCMPLX(FRe(12),FIm(12))
       P(35) = -s123-s234+s45+s56
       FRe(13) = -(s234*D1256R(2,2))+s234*D1256R(6,2)-D1256R(2,1)*P(19
     &   )+D1256R(3,1)*P(19)+D1256R(4,2)*P(30)+D1256R(5,2)*P(35)
       FIm(13) = -(s234*D1256I(2,2))+s234*D1256I(6,2)-D1256I(2,1)*P(19
     &   )+D1256I(3,1)*P(19)+D1256I(4,2)*P(30)+D1256I(5,2)*P(35)
       F1(13)=DCMPLX(FRe(13),FIm(13))
       FRe(14) = -(s234*D1256R(2,2))+s234*D1256R(6,2)-6*D1256R(7,2)-D1
     &   256R(2,1)*P(4)+D1256R(3,1)*P(4)-D1256R(4,2)*P(16)+D1256R(5,2
     &   )*P(16)
       FIm(14) = -(s234*D1256I(2,2))+s234*D1256I(6,2)-6*D1256I(7,2)-D1
     &   256I(2,1)*P(4)+D1256I(3,1)*P(4)-D1256I(4,2)*P(16)+D1256I(5,2
     &   )*P(16)
       F1(14)=DCMPLX(FRe(14),FIm(14))
       FRe(15) = -D1256R(2,1)-D1256R(2,2)+D1256R(3,1)+D1256R(6,2)
       FIm(15) = -D1256I(2,1)-D1256I(2,2)+D1256I(3,1)+D1256I(6,2)
       F1(15)=DCMPLX(FRe(15),FIm(15))
       FRe(16) = D1256R(2,1)-D1256R(3,1)+D1256R(4,2)+D1256R(5,2)-2*D12
     &   56R(6,2)
       FIm(16) = D1256I(2,1)-D1256I(3,1)+D1256I(4,2)+D1256I(5,2)-2*D12
     &   56I(6,2)
       F1(16)=DCMPLX(FRe(16),FIm(16))
       FRe(17) = D1256R(2,1)-D1256R(3,1)+D1256R(4,2)-D1256R(5,2)
       FIm(17) = D1256I(2,1)-D1256I(3,1)+D1256I(4,2)-D1256I(5,2)
       F1(17)=DCMPLX(FRe(17),FIm(17))
       FRe(18) = -D1256R(2,1)+D1256R(3,1)-D1256R(4,2)+D1256R(5,2)
       FIm(18) = -D1256I(2,1)+D1256I(3,1)-D1256I(4,2)+D1256I(5,2)
       F1(18)=DCMPLX(FRe(18),FIm(18))
       P(36) = s23-s234+s45+s56
       P(37) = s123**2+s23*s45-s123*P(36)
       P(38) = -2*s16+s23
       P(39) = s23-s234-s45+s56
       P(40) = -s123**2+2*p4sq*P(9)+s45*P(38)+s123*P(39)
       P(41) = -s16+s23-s234+s56
       P(42) = -s123**2+2*p4sq*P(9)+s123*P(39)+s45*P(41)
       P(43) = p4sq-s123+s23-s234+s56
       P(44) = s123+s16-s23
       P(45) = -s123+s23-s234+s56
       P(46) = -(s16*s45)+p4sq*P(44)+s123*P(45)
       P(47) = p4sq*P(44)+s45*P(45)
       P(48) = s23+s45
       P(49) = -(s16*s45)+s123*P(45)+p4sq*P(48)
       P(50) = s123-s23+s234-s56
       P(51) = 3*p4sq-s45-2*P(50)
       FRe(19) = 2*p4sq*s23*D1246R(2,2)-D1246R(1,1)*P(37)+D1246R(2,1)*
     &   P(40)-D1246R(3,1)*P(42)+2*s45*D1246R(3,2)*P(43)+2*D1246R(4,2
     &   )*P(46)-2*D1246R(5,2)*P(47)-2*D1246R(6,2)*P(49)+4*D1246R(7,2
     &   )*P(51)
       FIm(19) = 2*p4sq*s23*D1246I(2,2)-D1246I(1,1)*P(37)+D1246I(2,1)*
     &   P(40)-D1246I(3,1)*P(42)+2*s45*D1246I(3,2)*P(43)+2*D1246I(4,2
     &   )*P(46)-2*D1246I(5,2)*P(47)-2*D1246I(6,2)*P(49)+4*D1246I(7,2
     &   )*P(51)
       F1(19)=DCMPLX(FRe(19),FIm(19))
       P(52) = -s16+s23
       P(53) = p4sq*P(9)+s45*P(52)
       P(54) = s23-s234+s56
       P(55) = p4sq*P(9)+s45*P(54)
       P(56) = p4sq+s45
       P(57) = p4sq+s23-s234+s56
       P(58) = s16+s23
       P(59) = s123**2-s123*P(36)-p4sq*P(44)+s45*P(58)
       P(60) = s16+s23-2*P(1)
       P(61) = s123**2-s123*P(36)+p4sq*P(44)+s45*P(60)
       P(62) = -2*s16+s234-s56
       P(63) = s123**2-s123*P(36)
       P(64) = p4sq*P(48)+s45*P(62)-2*P(63)
       P(65) = p4sq+s23-s234
       FRe(20) = D1246R(1,1)*P(37)-D1246R(2,1)*P(53)+D1246R(3,1)*P(55)
     &   -s23*D1246R(2,2)*P(56)-s45*D1246R(3,2)*P(57)+D1246R(4,2)*P(5
     &   9)+D1246R(5,2)*P(61)+D1246R(6,2)*P(64)-6*D1246R(7,2)*P(65)
       FIm(20) = D1246I(1,1)*P(37)-D1246I(2,1)*P(53)+D1246I(3,1)*P(55)
     &   -s23*D1246I(2,2)*P(56)-s45*D1246I(3,2)*P(57)+D1246I(4,2)*P(5
     &   9)+D1246I(5,2)*P(61)+D1246I(6,2)*P(64)-6*D1246I(7,2)*P(65)
       F1(20)=DCMPLX(FRe(20),FIm(20))
       P(66) = -s123**2+s123*P(36)+p4sq*P(44)-s45*P(58)
       P(67) = p4sq+s123-s45-s56
       FRe(21) = -(D1246R(1,1)*P(37))+D1246R(2,1)*P(53)-D1246R(3,1)*P(
     &   55)+s23*D1246R(2,2)*P(56)+s45*D1246R(3,2)*P(57)-D1246R(5,2)*
     &   P(61)-D1246R(6,2)*P(64)+D1246R(4,2)*P(66)+6*D1246R(7,2)*P(67
     &   )
       FIm(21) = -(D1246I(1,1)*P(37))+D1246I(2,1)*P(53)-D1246I(3,1)*P(
     &   55)+s23*D1246I(2,2)*P(56)+s45*D1246I(3,2)*P(57)-D1246I(5,2)*
     &   P(61)-D1246I(6,2)*P(64)+D1246I(4,2)*P(66)+6*D1246I(7,2)*P(67
     &   )
       F1(21)=DCMPLX(FRe(21),FIm(21))
       P(68) = s123-s23
       P(69) = s16-s23
       P(70) = s16-s45
       P(71) = -s123+s16+s23
       P(72) = -2*s16+s23+s45
       FRe(22) = s23*D1246R(2,2)+s45*D1246R(3,2)+2*D1246R(7,2)-D1246R(
     &   5,2)*P(44)+D1246R(1,1)*P(68)-D1246R(2,1)*P(69)+D1246R(3,1)*P
     &   (70)-D1246R(4,2)*P(71)-D1246R(6,2)*P(72)
       FIm(22) = s23*D1246I(2,2)+s45*D1246I(3,2)+2*D1246I(7,2)-D1246I(
     &   5,2)*P(44)+D1246I(1,1)*P(68)-D1246I(2,1)*P(69)+D1246I(3,1)*P
     &   (70)-D1246I(4,2)*P(71)-D1246I(6,2)*P(72)
       F1(22)=DCMPLX(FRe(22),FIm(22))
       FRe(23) = D1246R(2,1)-D1246R(3,1)+D1246R(3,2)+D1246R(4,2)-D1246
     &   R(5,2)-D1246R(6,2)
       FIm(23) = D1246I(2,1)-D1246I(3,1)+D1246I(3,2)+D1246I(4,2)-D1246
     &   I(5,2)-D1246I(6,2)
       F1(23)=DCMPLX(FRe(23),FIm(23))
       P(73) = s123-s23+s234+s45-s56
       P(74) = s45-2*P(50)
       FRe(24) = -(s45*D1246R(3,2))+D1246R(2,1)*P(45)+D1246R(4,2)*P(45
     &   )+3*D1246R(5,2)*P(50)+D1246R(3,1)*P(73)+D1246R(6,2)*P(74)
       FIm(24) = -(s45*D1246I(3,2))+D1246I(2,1)*P(45)+D1246I(4,2)*P(45
     &   )+3*D1246I(5,2)*P(50)+D1246I(3,1)*P(73)+D1246I(6,2)*P(74)
       F1(24)=DCMPLX(FRe(24),FIm(24))
       P(75) = 2*s123+s16+s23
       P(76) = 2*s123+s16+s45
       P(77) = s16+3*P(68)
       P(78) = 2*s16-3*P(48)
       FRe(25) = 3*s23*D1246R(2,2)+3*s45*D1246R(3,2)+14*D1246R(7,2)-3*
     &   D1246R(5,2)*P(44)+D1246R(1,1)*P(68)+D1246R(2,1)*P(75)-D1246R
     &   (3,1)*P(76)+D1246R(4,2)*P(77)+D1246R(6,2)*P(78)
       FIm(25) = 3*s23*D1246I(2,2)+3*s45*D1246I(3,2)+14*D1246I(7,2)-3*
     &   D1246I(5,2)*P(44)+D1246I(1,1)*P(68)+D1246I(2,1)*P(75)-D1246I
     &   (3,1)*P(76)+D1246I(4,2)*P(77)+D1246I(6,2)*P(78)
       F1(25)=DCMPLX(FRe(25),FIm(25))
       P(79) = -8*s16+12*P(48)
       FRe(26) = 56*D1246R(7,2)+12*(s23*D1246R(2,2)+s45*D1246R(3,2)-D1
     &   246R(5,2)*P(44))+4*(D1246R(1,1)*P(68)+D1246R(2,1)*P(75)-D124
     &   6R(3,1)*P(76)+D1246R(4,2)*P(77))-D1246R(6,2)*P(79)
       FIm(26) = 56*D1246I(7,2)+12*(s23*D1246I(2,2)+s45*D1246I(3,2)-D1
     &   246I(5,2)*P(44))+4*(D1246I(1,1)*P(68)+D1246I(2,1)*P(75)-D124
     &   6I(3,1)*P(76)+D1246I(4,2)*P(77))-D1246I(6,2)*P(79)
       F1(26)=DCMPLX(FRe(26),FIm(26))
       FRe(27) = s23*D1246R(2,2)+s45*D1246R(3,2)+6*D1246R(7,2)+D1246R(
     &   2,1)*P(9)-D1246R(3,1)*P(9)+D1246R(4,2)*P(44)-D1246R(5,2)*P(4
     &   4)-D1246R(6,2)*P(48)
       FIm(27) = s23*D1246I(2,2)+s45*D1246I(3,2)+6*D1246I(7,2)+D1246I(
     &   2,1)*P(9)-D1246I(3,1)*P(9)+D1246I(4,2)*P(44)-D1246I(5,2)*P(4
     &   4)-D1246I(6,2)*P(48)
       F1(27)=DCMPLX(FRe(27),FIm(27))
       P(80) = s123+s45
       FRe(28) = s123*D1246R(1,1)+s123*D1246R(2,1)+2*s45*D1246R(3,2)+2
     &   *s123*D1246R(4,2)-2*s45*D1246R(5,2)-2*s123*D1246R(6,2)+8*D12
     &   46R(7,2)-D1246R(3,1)*P(80)
       FIm(28) = s123*D1246I(1,1)+s123*D1246I(2,1)+2*s45*D1246I(3,2)+2
     &   *s123*D1246I(4,2)-2*s45*D1246I(5,2)-2*s123*D1246I(6,2)+8*D12
     &   46I(7,2)-D1246I(3,1)*P(80)
       F1(28)=DCMPLX(FRe(28),FIm(28))
       P(81) = s123+s23
       FRe(29) = 2*s23*D1246R(2,2)+2*s45*D1246R(3,2)+8*D1246R(7,2)-2*D
     &   1246R(6,2)*P(36)-2*D1246R(5,2)*P(50)+D1246R(1,1)*P(68)+2*D12
     &   46R(4,2)*P(68)-D1246R(3,1)*P(80)+D1246R(2,1)*P(81)
       FIm(29) = 2*s23*D1246I(2,2)+2*s45*D1246I(3,2)+8*D1246I(7,2)-2*D
     &   1246I(6,2)*P(36)-2*D1246I(5,2)*P(50)+D1246I(1,1)*P(68)+2*D12
     &   46I(4,2)*P(68)-D1246I(3,1)*P(80)+D1246I(2,1)*P(81)
       F1(29)=DCMPLX(FRe(29),FIm(29))
       P(82) = -s16+s23+s45
       FRe(30) = 2*s23*D1246R(2,2)+2*s45*D1246R(3,2)+8*D1246R(7,2)-2*D
     &   1246R(5,2)*P(44)+D1246R(1,1)*P(68)+2*D1246R(4,2)*P(68)-D1246
     &   R(3,1)*P(80)+D1246R(2,1)*P(81)-2*D1246R(6,2)*P(82)
       FIm(30) = 2*s23*D1246I(2,2)+2*s45*D1246I(3,2)+8*D1246I(7,2)-2*D
     &   1246I(5,2)*P(44)+D1246I(1,1)*P(68)+2*D1246I(4,2)*P(68)-D1246
     &   I(3,1)*P(80)+D1246I(2,1)*P(81)-2*D1246I(6,2)*P(82)
       F1(30)=DCMPLX(FRe(30),FIm(30))
       P(83) = 2*s123+s234-s56
       P(84) = 2*s123-s23+s234+s45-s56
       P(85) = s234-s56+3*P(68)
       P(86) = 3*s45+2*P(54)
       FRe(31) = 2*s23*D1246R(2,2)+3*s45*D1246R(3,2)+14*D1246R(7,2)+3*
     &   D1246R(5,2)*P(45)+D1246R(1,1)*P(68)+D1246R(2,1)*P(83)-D1246R
     &   (3,1)*P(84)+D1246R(4,2)*P(85)-D1246R(6,2)*P(86)
       FIm(31) = 2*s23*D1246I(2,2)+3*s45*D1246I(3,2)+14*D1246I(7,2)+3*
     &   D1246I(5,2)*P(45)+D1246I(1,1)*P(68)+D1246I(2,1)*P(83)-D1246I
     &   (3,1)*P(84)+D1246I(4,2)*P(85)-D1246I(6,2)*P(86)
       F1(31)=DCMPLX(FRe(31),FIm(31))
       P(87) = -s234+s56-3*P(68)
       P(88) = 12*s45+8*P(54)
       FRe(32) = -8*s23*D1246R(2,2)-12*s45*D1246R(3,2)-56*D1246R(7,2)+
     &   12*D1246R(5,2)*P(50)-4*D1246R(1,1)*P(68)-4*D1246R(2,1)*P(83)
     &   +4*D1246R(3,1)*P(84)+4*D1246R(4,2)*P(87)+D1246R(6,2)*P(88)
       FIm(32) = -8*s23*D1246I(2,2)-12*s45*D1246I(3,2)-56*D1246I(7,2)+
     &   12*D1246I(5,2)*P(50)-4*D1246I(1,1)*P(68)-4*D1246I(2,1)*P(83)
     &   +4*D1246I(3,1)*P(84)+4*D1246I(4,2)*P(87)+D1246I(6,2)*P(88)
       F1(32)=DCMPLX(FRe(32),FIm(32))
       P(89) = s45+2*P(50)
       P(90) = -s123+s23-s234+s45+s56
       FRe(33) = 2*s45*D1246R(3,2)+4*D1246R(5,2)*P(45)+2*D1246R(2,1)*P
     &   (50)+2*D1246R(4,2)*P(50)-D1246R(3,1)*P(89)-2*D1246R(6,2)*P(9
     &   0)
       FIm(33) = 2*s45*D1246I(3,2)+4*D1246I(5,2)*P(45)+2*D1246I(2,1)*P
     &   (50)+2*D1246I(4,2)*P(50)-D1246I(3,1)*P(89)-2*D1246I(6,2)*P(9
     &   0)
       F1(33)=DCMPLX(FRe(33),FIm(33))
       FRe(34) = 4*(-2*s45*D1246R(3,2)+2*D1246R(2,1)*P(45)+2*D1246R(4,
     &   2)*P(45)+4*D1246R(5,2)*P(50)+D1246R(3,1)*P(89)+2*D1246R(6,2)
     &   *P(90))
       FIm(34) = 4*(-2*s45*D1246I(3,2)+2*D1246I(2,1)*P(45)+2*D1246I(4,
     &   2)*P(45)+4*D1246I(5,2)*P(50)+D1246I(3,1)*P(89)+2*D1246I(6,2)
     &   *P(90))
       F1(34)=DCMPLX(FRe(34),FIm(34))
       P(91) = 8*p4sq-4*s45
       FRe(35) = -8*(-(p4sq*D1246R(5,2))+s45*D1246R(6,2)+(D1246R(2,1)+
     &   D1246R(4,2))*P(15))+D1246R(3,1)*P(91)
       FIm(35) = -8*(-(p4sq*D1246I(5,2))+s45*D1246I(6,2)+(D1246I(2,1)+
     &   D1246I(4,2))*P(15))+D1246I(3,1)*P(91)
       F1(35)=DCMPLX(FRe(35),FIm(35))
       FRe(36) = 2*s23*D1246R(2,2)-s123*D1246R(3,1)+2*s45*D1246R(3,2)+
     &   8*D1246R(7,2)+D1246R(1,1)*P(68)+2*D1246R(4,2)*P(68)-2*D1246R
     &   (6,2)*P(80)+D1246R(2,1)*P(81)
       FIm(36) = 2*s23*D1246I(2,2)-s123*D1246I(3,1)+2*s45*D1246I(3,2)+
     &   8*D1246I(7,2)+D1246I(1,1)*P(68)+2*D1246I(4,2)*P(68)-2*D1246I
     &   (6,2)*P(80)+D1246I(2,1)*P(81)
       F1(36)=DCMPLX(FRe(36),FIm(36))
       P(92) = s123-2*s45
       P(93) = s123-s45
       FRe(37) = s123*D1246R(1,1)-6*D1246R(7,2)+D1246R(2,1)*P(92)-(D12
     &   46R(3,1)-2*D1246R(4,2)+2*D1246R(6,2))*P(93)
       FIm(37) = s123*D1246I(1,1)-6*D1246I(7,2)+D1246I(2,1)*P(92)-(D12
     &   46I(3,1)-2*D1246I(4,2)+2*D1246I(6,2))*P(93)
       F1(37)=DCMPLX(FRe(37),FIm(37))
       FRe(38) = s123*D1246R(1,1)+D1246R(2,1)*P(92)-(D1246R(3,1)-2*D12
     &   46R(4,2)+2*D1246R(6,2))*P(93)
       FIm(38) = s123*D1246I(1,1)+D1246I(2,1)*P(92)-(D1246I(3,1)-2*D12
     &   46I(4,2)+2*D1246I(6,2))*P(93)
       F1(38)=DCMPLX(FRe(38),FIm(38))
       P(94) = p4sq-2*s45
       FRe(39) = p4sq*D1246R(5,2)-2*s45*D1246R(6,2)+D1246R(3,1)*P(15)-
     &   D1246R(2,1)*P(94)-D1246R(4,2)*P(94)
       FIm(39) = p4sq*D1246I(5,2)-2*s45*D1246I(6,2)+D1246I(3,1)*P(15)-
     &   D1246I(2,1)*P(94)-D1246I(4,2)*P(94)
       F1(39)=DCMPLX(FRe(39),FIm(39))
       P(95) = 8*s123-4*s45
       P(96) = -2*s123+s45
       FRe(40) = D1246R(4,2)*P(95)+4*(s123*D1246R(1,1)-s123*D1246R(3,1
     &   )+s45*D1246R(3,2)-s45*D1246R(5,2)+D1246R(2,1)*P(93)+D1246R(6
     &   ,2)*P(96))
       FIm(40) = D1246I(4,2)*P(95)+4*(s123*D1246I(1,1)-s123*D1246I(3,1
     &   )+s45*D1246I(3,2)-s45*D1246I(5,2)+D1246I(2,1)*P(93)+D1246I(6
     &   ,2)*P(96))
       F1(40)=DCMPLX(FRe(40),FIm(40))
       P(97) = 2*s123+s45
       FRe(41) = 2*s23*D1246R(2,2)-s123*D1246R(3,1)+s45*D1246R(3,2)+D1
     &   246R(1,1)*P(68)+2*D1246R(4,2)*P(68)+D1246R(2,1)*P(81)-D1246R
     &   (6,2)*P(97)
       FIm(41) = 2*s23*D1246I(2,2)-s123*D1246I(3,1)+s45*D1246I(3,2)+D1
     &   246I(1,1)*P(68)+2*D1246I(4,2)*P(68)+D1246I(2,1)*P(81)-D1246I
     &   (6,2)*P(97)
       F1(41)=DCMPLX(FRe(41),FIm(41))
       P(98) = 2*s23-s234+s56
       P(99) = -s123+s23+s234-s56
       FRe(42) = -2*s23*D1246R(2,2)+2*s123*D1246R(6,2)+D1246R(5,2)*P(4
     &   5)+D1246R(3,1)*P(54)-D1246R(1,1)*P(68)-D1246R(2,1)*P(98)+D12
     &   46R(4,2)*P(99)
       FIm(42) = -2*s23*D1246I(2,2)+2*s123*D1246I(6,2)+D1246I(5,2)*P(4
     &   5)+D1246I(3,1)*P(54)-D1246I(1,1)*P(68)-D1246I(2,1)*P(98)+D12
     &   46I(4,2)*P(99)
       F1(42)=DCMPLX(FRe(42),FIm(42))
       FRe(43) = -(p4sq*D1246R(3,1))+2*s45*D1246R(6,2)+D1246R(2,1)*P(1
     &   5)+D1246R(4,2)*P(15)-D1246R(5,2)*P(56)
       FIm(43) = -(p4sq*D1246I(3,1))+2*s45*D1246I(6,2)+D1246I(2,1)*P(1
     &   5)+D1246I(4,2)*P(15)-D1246I(5,2)*P(56)
       F1(43)=DCMPLX(FRe(43),FIm(43))
       FRe(44) = -(p4sq*D1246R(5,2))+2*s45*D1246R(6,2)-D1246R(3,1)*P(1
     &   5)+D1246R(2,1)*P(94)+D1246R(4,2)*P(94)
       FIm(44) = -(p4sq*D1246I(5,2))+2*s45*D1246I(6,2)-D1246I(3,1)*P(1
     &   5)+D1246I(2,1)*P(94)+D1246I(4,2)*P(94)
       F1(44)=DCMPLX(FRe(44),FIm(44))
       FRe(45) = s45*D1246R(2,1)-s45*D1246R(3,1)+s45*D1246R(3,2)+s45*D
     &   1246R(4,2)-s45*D1246R(5,2)-s45*D1246R(6,2)+6*D1246R(7,2)
       FIm(45) = s45*D1246I(2,1)-s45*D1246I(3,1)+s45*D1246I(3,2)+s45*D
     &   1246I(4,2)-s45*D1246I(5,2)-s45*D1246I(6,2)+6*D1246I(7,2)
       F1(45)=DCMPLX(FRe(45),FIm(45))
       FRe(46) = s45*D1246R(3,2)-s45*D1246R(6,2)+6*D1246R(7,2)+D1246R(
     &   3,1)*P(45)+D1246R(5,2)*P(45)+D1246R(2,1)*P(50)+D1246R(4,2)*P
     &   (50)
       FIm(46) = s45*D1246I(3,2)-s45*D1246I(6,2)+6*D1246I(7,2)+D1246I(
     &   3,1)*P(45)+D1246I(5,2)*P(45)+D1246I(2,1)*P(50)+D1246I(4,2)*P
     &   (50)
       F1(46)=DCMPLX(FRe(46),FIm(46))
       FRe(47) = s45*D1246R(3,2)-s45*D1246R(6,2)+D1246R(3,1)*P(45)+D12
     &   46R(5,2)*P(45)+D1246R(2,1)*P(50)+D1246R(4,2)*P(50)
       FIm(47) = s45*D1246I(3,2)-s45*D1246I(6,2)+D1246I(3,1)*P(45)+D12
     &   46I(5,2)*P(45)+D1246I(2,1)*P(50)+D1246I(4,2)*P(50)
       F1(47)=DCMPLX(FRe(47),FIm(47))
       FRe(48) = -(s45*D1246R(3,2))+s45*D1246R(6,2)-3*D1246R(7,2)+D124
     &   6R(2,1)*P(45)+D1246R(4,2)*P(45)+D1246R(3,1)*P(50)+D1246R(5,2
     &   )*P(50)
       FIm(48) = -(s45*D1246I(3,2))+s45*D1246I(6,2)-3*D1246I(7,2)+D124
     &   6I(2,1)*P(45)+D1246I(4,2)*P(45)+D1246I(3,1)*P(50)+D1246I(5,2
     &   )*P(50)
       F1(48)=DCMPLX(FRe(48),FIm(48))
       FRe(49) = -(s45*D1246R(3,2))+s45*D1246R(6,2)+D1246R(2,1)*P(45)+
     &   D1246R(4,2)*P(45)+D1246R(3,1)*P(50)+D1246R(5,2)*P(50)
       FIm(49) = -(s45*D1246I(3,2))+s45*D1246I(6,2)+D1246I(2,1)*P(45)+
     &   D1246I(4,2)*P(45)+D1246I(3,1)*P(50)+D1246I(5,2)*P(50)
       F1(49)=DCMPLX(FRe(49),FIm(49))
       FRe(50) = -(s45*D1246R(3,2))+s45*D1246R(6,2)-6*D1246R(7,2)+D124
     &   6R(2,1)*P(45)+D1246R(4,2)*P(45)+D1246R(3,1)*P(50)+D1246R(5,2
     &   )*P(50)
       FIm(50) = -(s45*D1246I(3,2))+s45*D1246I(6,2)-6*D1246I(7,2)+D124
     &   6I(2,1)*P(45)+D1246I(4,2)*P(45)+D1246I(3,1)*P(50)+D1246I(5,2
     &   )*P(50)
       F1(50)=DCMPLX(FRe(50),FIm(50))
       P(100) = 4*s123-2*s45
       FRe(51) = -2*s123*D1246R(1,1)+2*s123*D1246R(3,1)-2*s45*D1246R(3
     &   ,2)+2*s45*D1246R(5,2)-2*D1246R(2,1)*P(93)+2*D1246R(4,2)*P(96
     &   )+D1246R(6,2)*P(100)
       FIm(51) = -2*s123*D1246I(1,1)+2*s123*D1246I(3,1)-2*s45*D1246I(3
     &   ,2)+2*s45*D1246I(5,2)-2*D1246I(2,1)*P(93)+2*D1246I(4,2)*P(96
     &   )+D1246I(6,2)*P(100)
       F1(51)=DCMPLX(FRe(51),FIm(51))
       FRe(52) = -(s45*D1246R(2,1))+s45*D1246R(3,1)-s45*D1246R(3,2)-s4
     &   5*D1246R(4,2)+s45*D1246R(5,2)+s45*D1246R(6,2)-6*D1246R(7,2)
       FIm(52) = -(s45*D1246I(2,1))+s45*D1246I(3,1)-s45*D1246I(3,2)-s4
     &   5*D1246I(4,2)+s45*D1246I(5,2)+s45*D1246I(6,2)-6*D1246I(7,2)
       F1(52)=DCMPLX(FRe(52),FIm(52))
       P(101) = 1/s23+1/s45
       FRe(53) = s23*D1246R(2,2)+s45*D1246R(3,2)+6*D1246R(7,2)+D1246R(
     &   2,1)*P(9)-D1246R(3,1)*P(9)+D1246R(4,2)*P(44)-D1246R(5,2)*P(4
     &   4)-s23*s45*D1246R(6,2)*P(101)
       FIm(53) = s23*D1246I(2,2)+s45*D1246I(3,2)+6*D1246I(7,2)+D1246I(
     &   2,1)*P(9)-D1246I(3,1)*P(9)+D1246I(4,2)*P(44)-D1246I(5,2)*P(4
     &   4)-s23*s45*D1246I(6,2)*P(101)
       F1(53)=DCMPLX(FRe(53),FIm(53))
       FRe(54) = -D1246R(2,1)+D1246R(3,1)-D1246R(3,2)-D1246R(4,2)+D124
     &   6R(5,2)+D1246R(6,2)
       FIm(54) = -D1246I(2,1)+D1246I(3,1)-D1246I(3,2)-D1246I(4,2)+D124
     &   6I(5,2)+D1246I(6,2)
       F1(54)=DCMPLX(FRe(54),FIm(54))
       FRe(55) = D1246R(2,1)-D1246R(3,1)+D1246R(4,2)-3*D1246R(5,2)+2*D
     &   1246R(6,2)
       FIm(55) = D1246I(2,1)-D1246I(3,1)+D1246I(4,2)-3*D1246I(5,2)+2*D
     &   1246I(6,2)
       F1(55)=DCMPLX(FRe(55),FIm(55))
       FRe(56) = D1246R(2,1)-D1246R(3,1)+D1246R(4,2)+D1246R(5,2)-2*D12
     &   46R(6,2)
       FIm(56) = D1246I(2,1)-D1246I(3,1)+D1246I(4,2)+D1246I(5,2)-2*D12
     &   46I(6,2)
       F1(56)=DCMPLX(FRe(56),FIm(56))
       FRe(57) = D1246R(2,1)-D1246R(3,1)+D1246R(4,2)-D1246R(5,2)
       FIm(57) = D1246I(2,1)-D1246I(3,1)+D1246I(4,2)-D1246I(5,2)
       F1(57)=DCMPLX(FRe(57),FIm(57))
       FRe(58) = D1246R(5,2)-D1246R(6,2)
       FIm(58) = D1246I(5,2)-D1246I(6,2)
       F1(58)=DCMPLX(FRe(58),FIm(58))
       FRe(59) = -D1246R(2,1)+D1246R(3,1)-D1246R(4,2)+D1246R(5,2)
       FIm(59) = -D1246I(2,1)+D1246I(3,1)-D1246I(4,2)+D1246I(5,2)
       F1(59)=DCMPLX(FRe(59),FIm(59))
       P(102) = s12+s16
       P(103) = s23-s234+s34
       P(104) = p3sq*P(102)-s16*P(103)
       P(105) = s23-s234
       P(106) = s123-2*s34-s56-3*P(105)
       P(107) = 2*P(104)+s12*P(106)
       P(108) = s123-2*s34-3*P(105)
       P(109) = -s23+s234
       P(110) = s12+s345
       P(111) = s123*s345+2*P(104)+s12*P(108)+s345*P(109)-s56*P(110)
       P(112) = p3sq+s123-s34-s56-2*P(105)
       P(113) = s123-s34-s56-2*P(105)
       P(114) = p3sq*P(102)-s16*P(103)+s12*P(113)
       P(115) = s16*s23-s16*s234+s16*s34-s123*s345+s23*s345-s234*s345+
     &   s345*s56-p3sq*P(102)+s12*P(103)
       P(116) = p3sq-s23+s234-s34
       P(117) = s12*P(50)+s345*P(116)
       P(118) = p3sq-s34
       P(119) = s123-s56
       P(120) = -8*P(105)+3*P(118)+5*P(119)
       FRe(60) = s12*D1236R(1,1)*P(50)+D1236R(2,1)*P(107)-D1236R(3,1)*
     &   P(111)+2*s345*D1236R(3,2)*P(112)+2*D1236R(4,2)*P(114)+2*D123
     &   6R(5,2)*P(115)-2*D1236R(6,2)*P(117)+4*D1236R(7,2)*P(120)
       FIm(60) = s12*D1236I(1,1)*P(50)+D1236I(2,1)*P(107)-D1236I(3,1)*
     &   P(111)+2*s345*D1236I(3,2)*P(112)+2*D1236I(4,2)*P(114)+2*D123
     &   6I(5,2)*P(115)-2*D1236I(6,2)*P(117)+4*D1236I(7,2)*P(120)
       F1(60)=DCMPLX(FRe(60),FIm(60))
       P(121) = -(s16*s23)+s16*s234-s16*s34+p3sq*P(102)
       P(122) = s123*s345-s23*s345+s234*s345-s345*s56+s12*P(106)+2*P(1
     &   21)
       P(123) = -(s16*s23)+s16*s234-s16*s34+s123*s345-s23*s345+s234*s3
     &   45-s345*s56+p3sq*P(102)-s12*P(103)
       FRe(61) = s12*D1236R(1,1)*P(45)-D1236R(2,1)*P(107)-2*s345*D1236
     &   R(3,2)*P(112)-2*D1236R(4,2)*P(114)+2*D1236R(6,2)*P(117)-4*D1
     &   236R(7,2)*P(120)+D1236R(3,1)*P(122)+2*D1236R(5,2)*P(123)
       FIm(61) = s12*D1236I(1,1)*P(45)-D1236I(2,1)*P(107)-2*s345*D1236
     &   I(3,2)*P(112)-2*D1236I(4,2)*P(114)+2*D1236I(6,2)*P(117)-4*D1
     &   236I(7,2)*P(120)+D1236I(3,1)*P(122)+2*D1236I(5,2)*P(123)
       F1(61)=DCMPLX(FRe(61),FIm(61))
       P(124) = -s123+s34+s56+2*P(105)
       P(125) = s16*s23-s16*s234+s16*s34-s123*s345+s23*s345-s234*s345+
     &   s345*s56-p3sq*P(102)+s12*P(124)
       P(126) = p3sq-s34-3*P(105)+2*P(119)
       P(127) = -s34-3*P(105)+2*P(119)
       P(128) = p3sq*P(102)-s16*P(103)+s12*P(127)
       P(129) = s123*s345-s23*s345+s234*s345-s345*s56
       P(130) = s16*s23-s16*s234+s16*s34-p3sq*P(102)+s12*P(103)-2*P(12
     &   9)
       P(131) = -p3sq+s23-s234+s34
       P(132) = -8*s12*P(50)+4*s345*P(131)
       P(133) = p3sq-s34+s45-2*P(105)
       FRe(62) = 4*(s12*D1236R(1,1)*P(50)+D1236R(2,1)*P(114)+D1236R(3,
     &   1)*P(125)+s345*D1236R(3,2)*P(126)+D1236R(4,2)*P(128)+D1236R(
     &   5,2)*P(130))+D1236R(6,2)*P(132)+24*D1236R(7,2)*P(133)
       FIm(62) = 4*(s12*D1236I(1,1)*P(50)+D1236I(2,1)*P(114)+D1236I(3,
     &   1)*P(125)+s345*D1236I(3,2)*P(126)+D1236I(4,2)*P(128)+D1236I(
     &   5,2)*P(130))+D1236I(6,2)*P(132)+24*D1236I(7,2)*P(133)
       F1(62)=DCMPLX(FRe(62),FIm(62))
       FRe(63) = s12*D1236R(1,1)+s12*D1236R(2,1)+2*s345*D1236R(3,2)+2*
     &   s12*D1236R(4,2)-2*s345*D1236R(5,2)-2*s12*D1236R(6,2)+6*D1236
     &   R(7,2)-D1236R(3,1)*P(110)
       FIm(63) = s12*D1236I(1,1)+s12*D1236I(2,1)+2*s345*D1236I(3,2)+2*
     &   s12*D1236I(4,2)-2*s345*D1236I(5,2)-2*s12*D1236I(6,2)+6*D1236
     &   I(7,2)-D1236I(3,1)*P(110)
       F1(63)=DCMPLX(FRe(63),FIm(63))
       FRe(64) = D1236R(1,1)+D1236R(4,2)+D1236R(5,2)-2*D1236R(6,2)
       FIm(64) = D1236I(1,1)+D1236I(4,2)+D1236I(5,2)-2*D1236I(6,2)
       F1(64)=DCMPLX(FRe(64),FIm(64))
       FRe(65) = D1236R(2,1)-D1236R(3,1)+D1236R(4,2)-D1236R(5,2)
       FIm(65) = D1236I(2,1)-D1236I(3,1)+D1236I(4,2)-D1236I(5,2)
       F1(65)=DCMPLX(FRe(65),FIm(65))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
