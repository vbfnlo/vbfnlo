c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:8/3/2012
c************************************************************************************
c************************************************************************************
c hel = 66,67,76,77 for rr,rl,...
c or 2*hel1 + hel2 -> 3,1,-1,-3

       subroutine HexBoxCross2T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div,hel)
       implicit none
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16 barpsi_p1(4),psi_p6(4),barpsi_p2(4),psi_p5(4)
       Complex*16   mup3(0:3),mup4(0:3)
       Real*8 musq
       Integer comp,ngluon,posgluon,Div,hel
       Complex*16 result(9),resultn,resultborn

       if (hel.eq.77.or.hel.eq.-3) then
         call HexBoxCross277T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       elseif (hel.eq.67.or.hel.eq.1) then
         call HexBoxCross267T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       elseif (hel.eq.76.or.hel.eq.-1) then
         call HexBoxCross276T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       elseif (hel.eq.66.or.hel.eq.3) then
         call HexBoxCross266T(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_
     &   p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn,Div)
       else
	 print*,'wrong helicity selected in HexBoxtempT.F'
       endif
       Return
       End

       subroutine HexBoxCross266T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn,Div)
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE
      Complex*16 Fa(392),F(1853)

     &


       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Real*8 P(1355)
       Complex*16  SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc
       EXTERNAL    SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc 
       Real*8 delta
       External  delta
       Complex*16 ten2cc
       External  ten2cc

       Integer comp
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup3(0:3),mup6(0:3)


       Real*8 dotrr
       EXTERNAL dotrr
      Complex*16 B0fin,C0fin,D0fin,E0fin,F0fin
       EXTERNAL B0fin,C0fin,D0fin,E0fin,F0fin
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

       Real*8 F0123456R 
       Real*8 F0123456I 
       Real*8 F123456R(166,5) 
       Real*8 F123456I(166,5) 


      Complex*16 dotrc,dotcc,result(9),resultn,resultborn
       Real*8 musq
      EXTERNAL   dotrc,dotcc

       COMMON/HexBoxCross2FaFunctions/Fa
       COMMON/HexBoxCross2FhlFunctions/F
      Save/HexBoxCross2FhlFunctions/
       COMMON/HexBoxCross2PFunctions/P
      Save/HexBoxCross2PFunctions/
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
       COMMON/HexBoxCross2FVALFunctions/F0123456R,F123456R,F0123456I,F
     &   123456I
       COMMON/HexBoxCross2EVALFunctions/ E012345R,E12345R,E012345I,E12
     &   345I, E012346R,E12346R,E012346I,E12346I, E012356R,E12356R,E0
     &   12356I,E12356I, E012456R,E12456R,E012456I,E12456I, E013456R,
     &   E13456R,E013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/HexBoxCross2DVALFunctions/ D01234R,D1234R,D01234I,D1234I
     &   , D01235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D123
     &   6I, D01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1
     &   246I, D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,
     &   D1345I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356
     &   I,D1356I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D023
     &   45I,D2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D0
     &   2356I,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,
     &   D03456I,D3456I
       COMMON/HexBoxCross2CVALFunctions/ C0123R,C123R,C0123I,C123I, C0
     &   124R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C
     &   126R,C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C
     &   0135I,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,
     &   C145I, C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I,
     &    C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236
     &   R,C236R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246
     &   R,C0246I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C034
     &   5I,C345I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C35
     &   6I, C0456R,C456R,C0456I,C456I
       COMMON/HexBoxCross2BVALFunctions/ B012R,B012I, B013R,B013I, B01
     &   4R,B014I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,B024I
     &   , B025R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,
     &   B036I, B045R,B045I, B046R,B046I, B056R,B056I
       Integer ngluon, posgluon,Div

       if(div.eq.0) then
       call HexBoxCross266(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn)
       else
       call HexBoxCross266Div(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn,Div)

       endif


       RETURN
       END
      
       subroutine HexBoxCross277T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn,Div)
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE
      Complex*16 Fa(392),F(1853)

     &


       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Real*8 P(1355)
       Complex*16  SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc
       EXTERNAL    SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc 
       Real*8 delta
       External  delta
       Complex*16 ten2cc
       External  ten2cc

       Integer comp
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup3(0:3),mup6(0:3)

       Real*8 dotrr
       EXTERNAL dotrr
      Complex*16 B0fin,C0fin,D0fin,E0fin,F0fin
       EXTERNAL B0fin,C0fin,D0fin,E0fin,F0fin
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

       Real*8 F0123456R 
       Real*8 F0123456I 
       Real*8 F123456R(166,5) 
       Real*8 F123456I(166,5) 


      Complex*16 dotrc,dotcc,result(9),resultn,resultborn
       Real*8 musq
      EXTERNAL   dotrc,dotcc

       COMMON/HexBoxCross2FaFunctions/Fa
       COMMON/HexBoxCross2FhlFunctions/F
      Save/HexBoxCross2FhlFunctions/
       COMMON/HexBoxCross2PFunctions/P
      Save/HexBoxCross2PFunctions/
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
       COMMON/HexBoxCross2FVALFunctions/F0123456R,F123456R,F0123456I,F
     &   123456I
       COMMON/HexBoxCross2EVALFunctions/ E012345R,E12345R,E012345I,E12
     &   345I, E012346R,E12346R,E012346I,E12346I, E012356R,E12356R,E0
     &   12356I,E12356I, E012456R,E12456R,E012456I,E12456I, E013456R,
     &   E13456R,E013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/HexBoxCross2DVALFunctions/ D01234R,D1234R,D01234I,D1234I
     &   , D01235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D123
     &   6I, D01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1
     &   246I, D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,
     &   D1345I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356
     &   I,D1356I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D023
     &   45I,D2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D0
     &   2356I,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,
     &   D03456I,D3456I
       COMMON/HexBoxCross2CVALFunctions/ C0123R,C123R,C0123I,C123I, C0
     &   124R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C
     &   126R,C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C
     &   0135I,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,
     &   C145I, C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I,
     &    C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236
     &   R,C236R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246
     &   R,C0246I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C034
     &   5I,C345I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C35
     &   6I, C0456R,C456R,C0456I,C456I
       COMMON/HexBoxCross2BVALFunctions/ B012R,B012I, B013R,B013I, B01
     &   4R,B014I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,B024I
     &   , B025R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,
     &   B036I, B045R,B045I, B046R,B046I, B056R,B056I
       Integer ngluon, posgluon,Div

       if(div.eq.0) then
       call HexBoxCross277(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn)
       else
       call HexBoxCross277Div(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn,Div)

       endif


       RETURN
       END


       subroutine HexBoxCross267T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn,Div)
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE
      Complex*16 Fa(392),F(1853)

     &


       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Real*8 P(1355)
       Complex*16  SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc
       EXTERNAL    SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc 
       Real*8 delta
       External  delta
       Complex*16 ten2cc
       External  ten2cc

       Integer comp
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup3(0:3),mup6(0:3)

       Real*8 dotrr
       EXTERNAL dotrr
      Complex*16 B0fin,C0fin,D0fin,E0fin,F0fin
       EXTERNAL B0fin,C0fin,D0fin,E0fin,F0fin
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

       Real*8 F0123456R 
       Real*8 F0123456I 
       Real*8 F123456R(166,5) 
       Real*8 F123456I(166,5) 


      Complex*16 dotrc,dotcc,result(9),resultn,resultborn
       Real*8 musq
      EXTERNAL   dotrc,dotcc

       COMMON/HexBoxCross2FaFunctions/Fa
       COMMON/HexBoxCross2FhlFunctions/F
      Save/HexBoxCross2FhlFunctions/
       COMMON/HexBoxCross2PFunctions/P
      Save/HexBoxCross2PFunctions/
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
       COMMON/HexBoxCross2FVALFunctions/F0123456R,F123456R,F0123456I,F
     &   123456I
       COMMON/HexBoxCross2EVALFunctions/ E012345R,E12345R,E012345I,E12
     &   345I, E012346R,E12346R,E012346I,E12346I, E012356R,E12356R,E0
     &   12356I,E12356I, E012456R,E12456R,E012456I,E12456I, E013456R,
     &   E13456R,E013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/HexBoxCross2DVALFunctions/ D01234R,D1234R,D01234I,D1234I
     &   , D01235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D123
     &   6I, D01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1
     &   246I, D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,
     &   D1345I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356
     &   I,D1356I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D023
     &   45I,D2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D0
     &   2356I,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,
     &   D03456I,D3456I
       COMMON/HexBoxCross2CVALFunctions/ C0123R,C123R,C0123I,C123I, C0
     &   124R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C
     &   126R,C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C
     &   0135I,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,
     &   C145I, C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I,
     &    C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236
     &   R,C236R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246
     &   R,C0246I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C034
     &   5I,C345I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C35
     &   6I, C0456R,C456R,C0456I,C456I
       COMMON/HexBoxCross2BVALFunctions/ B012R,B012I, B013R,B013I, B01
     &   4R,B014I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,B024I
     &   , B025R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,
     &   B036I, B045R,B045I, B046R,B046I, B056R,B056I
       Integer ngluon, posgluon,Div

       if(div.eq.0) then
       call HexBoxCross267(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn)
       else
       call HexBoxCross267Div(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn,Div)

       endif


       RETURN
       END


       subroutine HexBoxCross276T(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn,Div)
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE
      Complex*16 Fa(392),F(1853)

     &


       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Real*8 P(1355)
       Complex*16  SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc
       EXTERNAL    SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc 
       Real*8 delta
       External  delta
       Complex*16 ten2cc
       External  ten2cc

       Integer comp
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Complex*16   mup3(0:3),mup6(0:3)

       Real*8 dotrr
       EXTERNAL dotrr
      Complex*16 B0fin,C0fin,D0fin,E0fin,F0fin
       EXTERNAL B0fin,C0fin,D0fin,E0fin,F0fin
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

       Real*8 F0123456R 
       Real*8 F0123456I 
       Real*8 F123456R(166,5) 
       Real*8 F123456I(166,5) 


      Complex*16 dotrc,dotcc,result(9),resultn,resultborn
       Real*8 musq
      EXTERNAL   dotrc,dotcc

       COMMON/HexBoxCross2FaFunctions/Fa
       COMMON/HexBoxCross2FhlFunctions/F
      Save/HexBoxCross2FhlFunctions/
       COMMON/HexBoxCross2PFunctions/P
      Save/HexBoxCross2PFunctions/
       COMMON/HexBoxCross2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
       COMMON/HexBoxCross2FVALFunctions/F0123456R,F123456R,F0123456I,F
     &   123456I
       COMMON/HexBoxCross2EVALFunctions/ E012345R,E12345R,E012345I,E12
     &   345I, E012346R,E12346R,E012346I,E12346I, E012356R,E12356R,E0
     &   12356I,E12356I, E012456R,E12456R,E012456I,E12456I, E013456R,
     &   E13456R,E013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/HexBoxCross2DVALFunctions/ D01234R,D1234R,D01234I,D1234I
     &   , D01235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D123
     &   6I, D01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1
     &   246I, D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,
     &   D1345I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356
     &   I,D1356I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D023
     &   45I,D2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D0
     &   2356I,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,
     &   D03456I,D3456I
       COMMON/HexBoxCross2CVALFunctions/ C0123R,C123R,C0123I,C123I, C0
     &   124R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C
     &   126R,C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C
     &   0135I,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,
     &   C145I, C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I,
     &    C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236
     &   R,C236R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246
     &   R,C0246I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C034
     &   5I,C345I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C35
     &   6I, C0456R,C456R,C0456I,C456I
       COMMON/HexBoxCross2BVALFunctions/ B012R,B012I, B013R,B013I, B01
     &   4R,B014I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,B024I
     &   , B025R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,
     &   B036I, B045R,B045I, B046R,B046I, B056R,B056I
       Integer ngluon, posgluon,Div

       if(div.eq.0) then
       call HexBoxCross276(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn)
       else
       call HexBoxCross276Div(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5,ba
     &   rpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,re
     &   sultn,resultborn,Div)

       endif


       RETURN
       END
