c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:22/3/2012
c ************************************************************************************
c determine the  finite virtual corrections of 
c psi(p6) ---->--$$$$photon$$$$$$--->---   psi(p5)
c                   |              |              
c                   |              |$$$$$$$$$$$$  p4,mu_p4 
c                   |              |             
c                   |              |$$$$$$$$$$$$  p3, mu_p3 
c                   |              |             
c                   |              |              
c                   |              |             
c barpsi(p1)-->-$$$photon$$$$$$------>---   bar_psi(p2)
c                                               
c Note: To make it shorter in the promgram: mu_p3,...->mup3,... 
c Notation of External momenta: p1+p2+p3+p4+p5+p6=0 
c mu_p3, should be think as external current 
c alpha is the helicity of the initial spinor 
c musq is the renormalization scale energy  
c comp: integer value.The first time called with p1...p5, comp=1
c ATTENTION: ONLY!!!If you have to call the subroutine consecutively with the same arguments
c(p1,p2,p3,p4,p5). Then, comp=-1 (it safes 4000 lines of code) 
c This applies when you have for examples the same diagram for an off-shell photon
c and a Z boson. The differences are in the coupling and  the part that depends on the
c polarization vector that are calculated at the end of this program.
c************************************************************************************
c************************************************************************************
       subroutine HexBoxCross266Div(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_p5
     &   ,barpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result
     &   ,resultn,resultborn,Div)
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
      Complex*16 SMB(303),Fa(392),F(1853)
       Complex*16 SMB0(32),SMB11(0:3),SMB12(0:3),SMB13(0:3),SMB14(0:3)
     &   ,SMB15(0:3),SMB16(0:3),SMB17(0:3),SMB18(0:3)
      Complex*16 SMB1(133)

      Complex*16 psi_p1_P(2),barpsi_p5_P(2),barpsi_p2_P(2),psi_p4_P(2) 
       Complex*16 psi_p1_M(2),barpsi_p5_M(2),barpsi_p2_M(2),psi_p4_M(2)
       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Real*8 P(1355)
       Complex*16  SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc
       EXTERNAL    SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5ccrrc 
       Real*8 delta
       External  delta
       Complex*16 ten2cc
       External  ten2cc
       Complex*16 v1(0:3) 
       Integer i,comp
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
       Complex*16   p1mup3, p1mup6, p2mup3, p2mup6, p3mup3, p3mup6, 
     -          p4mup3, p4mup6, p5mup3, p5mup6, p6mup3, p6mup6
       Complex*16   mup3mup6
       Real*8 dotrr
       EXTERNAL dotrr
      Complex*16 B0finDiv,C0finDiv,D0finDiv,E0fin,F0fin
       EXTERNAL B0finDiv,C0finDiv,D0finDiv,E0fin,F0fin
        Real*8   p1sq, p1p2, p1p3, p1p4, p1p5,p1p6 
       Real*8   p2sq, p2p3, p2p4, p2p5 
       Real*8   p3sq, p3p4, p3p5 
       Real*8   p4sq, p4p5 
       Real*8   p5sq, p5p6 
       Real*8   p6sq 
       Real*8   s12, s13, s14, s15, s16 
       Real*8   s23, s24, s25,s26 
       Real*8   s34, s35,s36 
       Real*8   s45,s46 
       Real*8   s56 
       Real*8   s123,s234,s345   
       Real*8   p2p6,p3p6,p4p6  
       Complex*16  B012,B013,B014,B015,B016 
       Complex*16  B023,B024,B025,B026 
       Complex*16  B034,B035,B036 
       Complex*16  B045,B046 
       Complex*16  B056   
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
       Complex*16 C0123,C0124,C0125,C0126 
       Complex*16 C0134,C0135,C0136 
       Complex*16 C0145,C0146 
       Complex*16 C0156 
       Complex*16 C0234,C0235,C0236 
       Complex*16 C0245,C0246 
       Complex*16 C0256 
       Complex*16 C0345,C0346 
       Complex*16 C0356 
       Complex*16 C0456   
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
       Complex*16  D01234,D01235,D01236 
       Complex*16 D01245,D01246 
       Complex*16 D01256 
       Complex*16 D01345,D01346 
       Complex*16 D01356 
       Complex*16 D01456 
       Complex*16 D02345,D02346 
       Complex*16 D02356 
       Complex*16 D02456 
       Complex*16 D03456    
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
       Complex*16 E012345,E012346,E012356,E012456,E013456,E023456 
       Real*8 E012345R,E012346R,E012356R,E012456R,E013456R,E023456R 
       Real*8 E012345I,E012346I,E012356I,E012456I,E013456I,E023456I  
       Real*8 E12345R(46,4),E12346R(46,4),E12356R(46,4) 
       Real*8 E12456R(46,4),E13456R(46,4),E23456R(46,4)  
       Real*8 E12345I(46,4),E12346I(46,4),E12356I(46,4) 
       Real*8 E12456I(46,4),E13456I(46,4),E23456I(46,4)  
       Complex*16 F0123456 
       Real*8 F0123456R 
       Real*8 F0123456I 
       Real*8 F123456R(166,5) 
       Real*8 F123456I(166,5) 
       Complex*16 F123456(166,5) 
       Logical PrintB,PrintC,PrintD,PrintE,PrintF 
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
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the scalar products. Not inlcueded the contraction of the
c       moments with the external currents  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p1p5 = dotrr(p1,p5)
       p1p6 = dotrr(p1,p6)
       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)
       p2p6 = dotrr(p2,p6)
       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
       p3p6 = dotrr(p3,p6)
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)
       p4p6 = dotrr(p4,p6)
       p5sq = dotrr(p5,p5)
       p5p6 = dotrr(p5,p6)
       p6sq = dotrr(p6,p6)
       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s15 = (p1sq +p5sq+ 2*p1p5) 
       s16 = (p1sq +p6sq+ 2*p1p6) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s25 = (p2sq +p5sq+ 2*p2p5) 
       s26 = (p2sq +p6sq+ 2*p2p6) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s35 = (p3sq +p5sq+ 2*p3p5) 
       s36 = (p3sq +p6sq+ 2*p3p6) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
       s46 = (p4sq +p6sq+ 2*p4p6) 
       s56 = (p5sq +p6sq+ 2*p5p6) 
c       Write(*,'(a5,F20.10)')," p1sq ", p1sq 
c       Write(*,'(a5,F20.10)')," p1p2 ", p1p2
c       Write(*,'(a5,F20.10)')," p1p3 ", p1p3
c       Write(*,'(a5,F20.10)')," p1p4 ", p1p4
c       Write(*,'(a5,F20.10)')," p1p5 ", p1p5
c       Write(*,'(a5,F20.10)')," p1p6 ", p1p6
c       Write(*,'(a5,F20.10)')," p2sq ", p2sq 
c       Write(*,'(a5,F20.10)')," p2p3 ", p2p3
c       Write(*,'(a5,F20.10)')," p2p4 ", p2p4
c       Write(*,'(a5,F20.10)')," p2p5 ", p2p5
c       Write(*,'(a5,F20.10)')," p2p6 ", p2p6
c       Write(*,'(a5,F20.10)')," p3sq ", p3sq 
c       Write(*,'(a5,F20.10)')," p3p4 ", p3p4
c       Write(*,'(a5,F20.10)')," p3p5 ", p3p5
c       Write(*,'(a5,F20.10)')," p3p6 ", p3p6
c       Write(*,'(a5,F20.10)')," p4sq ", p4sq 
c       Write(*,'(a5,F20.10)')," p4p5 ", p4p5
c       Write(*,'(a5,F20.10)')," p4p6 ", p4p6
c       Write(*,'(a5,F20.10)')," p5sq ", p5sq 
c       Write(*,'(a5,F20.10)')," p5p6 ", p5p6
c       Write(*,'(a5,F20.10)')," p6sq ", p6sq 
      s123=p1sq+p2sq+p3sq+2*(p1p2+p1p3+p2p3) 
      s234=p2sq+p3sq+p4sq+2*(p2p3+p2p4+p3p4) 
      s345=p3sq+p4sq+p5sq+2*(p3p4+p3p5+p4p5)
      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False. 
      PrintF=.False.
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
        if (comp.gt.0) then    
c    Calling C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************
       B012=B0finDiv(p1sq,musq,Div)
       B012R=Dble(B012)
       B012I=DImag(B012)
       B023=B0finDiv(p2sq,musq,Div)
       B023R=Dble(B023)
       B023I=DImag(B023)
       B034=B0finDiv(p3sq,musq,Div)
       B034R=Dble(B034)
       B034I=DImag(B034)
       B045=B0finDiv(p4sq,musq,Div)
       B045R=Dble(B045)
       B045I=DImag(B045)
       B056=B0finDiv(p5sq,musq,Div)
       B056R=Dble(B056)
       B056I=DImag(B056)
       B016=B0finDiv(p6sq,musq,Div)
       B016R=Dble(B016)
       B016I=DImag(B016)
       B013=B0finDiv(s12,musq,Div)
       B013R=Dble(B013)
       B013I=DImag(B013)
       B014=B0finDiv(s123,musq,Div)
       B014R=Dble(B014)
       B014I=DImag(B014)
       B026=B0finDiv(s16,musq,Div)
       B026R=Dble(B026)
       B026I=DImag(B026)
       B024=B0finDiv(s23,musq,Div)
       B024R=Dble(B024)
       B024I=DImag(B024)
       B025=B0finDiv(s234,musq,Div)
       B025R=Dble(B025)
       B025I=DImag(B025)
       B035=B0finDiv(s34,musq,Div)
       B035R=Dble(B035)
       B035I=DImag(B035)
       B036=B0finDiv(s345,musq,Div)
       B036R=Dble(B036)
       B036I=DImag(B036)
       B046=B0finDiv(s45,musq,Div)
       B046R=Dble(B046)
       B046I=DImag(B046)
       B015=B0finDiv(s56,musq,Div)
       B015R=Dble(B015)
       B015I=DImag(B015)
       C0123=C0finDiv(p1sq,p2sq,s12,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,p2sq,s12,B023,B013,B012,C012
     &   3,C0123R,C0123I,C123R,C123I)
       C0126=C0finDiv(p1sq,s16,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,s16,p6sq,B026,B016,B012,C012
     &   6,C0126R,C0126I,C126R,C126I)
       C0124=C0finDiv(p1sq,s23,s123,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,s23,s123,B024,B014,B012,C012
     &   4,C0124R,C0124I,C124R,C124I)
       C0125=C0finDiv(p1sq,s234,s56,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,s234,s56,B025,B015,B012,C012
     &   5,C0125R,C0125I,C125R,C125I)
       C0234=C0finDiv(p2sq,p3sq,s23,musq,Div)
       call tens_red3_new_Re_Com_Div(p2sq,p3sq,s23,B034,B024,B023,C023
     &   4,C0234R,C0234I,C234R,C234I)
       C0235=C0finDiv(p2sq,s34,s234,musq,Div)
       call tens_red3_new_Re_Com_Div(p2sq,s34,s234,B035,B025,B023,C023
     &   5,C0235R,C0235I,C235R,C235I)
       C0236=C0finDiv(p2sq,s345,s16,musq,Div)
       call tens_red3_new_Re_Com_Div(p2sq,s345,s16,B036,B026,B023,C023
     &   6,C0236R,C0236I,C236R,C236I)
       C0345=C0finDiv(p3sq,p4sq,s34,musq,Div)
       call tens_red3_new_Re_Com_Div(p3sq,p4sq,s34,B045,B035,B034,C034
     &   5,C0345R,C0345I,C345R,C345I)
       C0346=C0finDiv(p3sq,s45,s345,musq,Div)
       call tens_red3_new_Re_Com_Div(p3sq,s45,s345,B046,B036,B034,C034
     &   6,C0346R,C0346I,C346R,C346I)
       C0456=C0finDiv(p4sq,p5sq,s45,musq,Div)
       call tens_red3_new_Re_Com_Div(p4sq,p5sq,s45,B056,B046,B045,C045
     &   6,C0456R,C0456I,C456R,C456I)
       C0134=C0finDiv(s12,p3sq,s123,musq,Div)
       call tens_red3_new_Re_Com_Div(s12,p3sq,s123,B034,B014,B013,C013
     &   4,C0134R,C0134I,C134R,C134I)
       C0135=C0finDiv(s12,s34,s56,musq,Div)
       call tens_red3_new_Re_Com_Div(s12,s34,s56,B035,B015,B013,C0135,
     &   C0135R,C0135I,C135R,C135I)
       C0136=C0finDiv(s12,s345,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(s12,s345,p6sq,B036,B016,B013,C013
     &   6,C0136R,C0136I,C136R,C136I)
       C0145=C0finDiv(s123,p4sq,s56,musq,Div)
       call tens_red3_new_Re_Com_Div(s123,p4sq,s56,B045,B015,B014,C014
     &   5,C0145R,C0145I,C145R,C145I)
       C0146=C0finDiv(s123,s45,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(s123,s45,p6sq,B046,B016,B014,C014
     &   6,C0146R,C0146I,C146R,C146I)
       C0245=C0finDiv(s23,p4sq,s234,musq,Div)
       call tens_red3_new_Re_Com_Div(s23,p4sq,s234,B045,B025,B024,C024
     &   5,C0245R,C0245I,C245R,C245I)
       C0246=C0finDiv(s23,s45,s16,musq,Div)
       call tens_red3_new_Re_Com_Div(s23,s45,s16,B046,B026,B024,C0246,
     &   C0246R,C0246I,C246R,C246I)
       C0256=C0finDiv(s234,p5sq,s16,musq,Div)
       call tens_red3_new_Re_Com_Div(s234,p5sq,s16,B056,B026,B025,C025
     &   6,C0256R,C0256I,C256R,C256I)
       C0356=C0finDiv(s34,p5sq,s345,musq,Div)
       call tens_red3_new_Re_Com_Div(s34,p5sq,s345,B056,B036,B035,C035
     &   6,C0356R,C0356I,C356R,C356I)
       C0156=C0finDiv(s56,p5sq,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(s56,p5sq,p6sq,B056,B016,B015,C015
     &   6,C0156R,C0156I,C156R,C156I)
       D01234=D0finDiv(s12,s23,p1sq,p2sq,p3sq,s123,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C02
     &   34R,C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C013
     &   4I,C0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D012
     &   34I,D1234R,D1234I)
       D01235=D0finDiv(s12,s234,p1sq,p2sq,s34,s56,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3
     &   +p2p4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0
     &   235I,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235,D01
     &   235R,D01235I,D1235R,D1235I)
       D01236=D0finDiv(s12,s16,p1sq,p2sq,s345,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p
     &   5,p2p3+p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C12
     &   6R,C123R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123I
     &   ,D01236,D01236R,D01236I,D1236R,D1236I)
       D01245=D0finDiv(s123,s234,p1sq,s23,p4sq,s56,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4
     &   +p3p4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0
     &   245I,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245,D01
     &   245R,D01245I,D1245R,D1245I)
       D01246=D0finDiv(s123,s16,p1sq,s23,s45,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,
     &   p2p4+p2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R,
     &   C126R,C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C1
     &   24I,D01246,D01246R,D01246I,D1246R,D1246I)
       D01256=D0finDiv(s56,s16,p1sq,s234,p5sq,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p
     &   5,p2p5+p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C12
     &   6R,C125R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125I
     &   ,D01256,D01256R,D01256I,D1256R,D1256I)
       D02345=D0finDiv(s23,s34,p2sq,p3sq,p4sq,s234,musq,Div)
       call tens_red4_new_Re_Com_Div(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C03
     &   45R,C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C024
     &   5I,C0235I,C0234I,C345I,C245I,C235I,C234I,D02345,D02345R,D023
     &   45I,D2345R,D2345I)
       D02346=D0finDiv(s23,s345,p2sq,p3sq,s45,s16,musq,Div)
       call tens_red4_new_Re_Com_Div(p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4
     &   +p3p5,C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C0
     &   346I,C0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346,D02
     &   346R,D02346I,D2346R,D2346I)
       D02356=D0finDiv(s234,s345,p2sq,s34,p5sq,s16,musq,Div)
       call tens_red4_new_Re_Com_Div(p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5
     &   +p4p5,C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C0
     &   356I,C0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356,D02
     &   356R,D02356I,D2356R,D2356I)
       D03456=D0finDiv(s34,s45,p3sq,p4sq,p5sq,s345,musq,Div)
       call tens_red4_new_Re_Com_Div(p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C04
     &   56R,C0356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C035
     &   6I,C0346I,C0345I,C456I,C356I,C346I,C345I,D03456,D03456R,D034
     &   56I,D3456R,D3456I)
       D01345=D0finDiv(s123,s34,s12,p3sq,p4sq,s56,musq,Div)
       call tens_red4_new_Re_Com_Div(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4
     &   ,p3p4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0
     &   345I,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345,D01
     &   345R,D01345I,D1345R,D1345I)
       D01346=D0finDiv(s123,s345,s12,p3sq,s45,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+
     &   p2p4+p2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R,
     &   C136R,C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C1
     &   34I,D01346,D01346R,D01346I,D1346R,D1346I)
       D01356=D0finDiv(s56,s345,s12,s34,p5sq,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,
     &   p1p5+p2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R,
     &   C136R,C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C1
     &   35I,D01356,D01356R,D01356I,D1356R,D1356I)
       D01456=D0finDiv(s56,s45,s123,p4sq,p5sq,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p
     &   5+p2p5+p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C14
     &   6R,C145R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145I
     &   ,D01456,D01456R,D01456I,D1456R,D1456I)
       D02456=D0finDiv(s234,s45,s23,p4sq,p5sq,s16,musq,Div)
       call tens_red4_new_Re_Com_Div(s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5
     &   ,p4p5,C0456R,C0256R,C0246R,C0245R,C456R,C256R,C246R,C245R,C0
     &   456I,C0256I,C0246I,C0245I,C456I,C256I,C246I,C245I,D02456,D02
     &   456R,D02456I,D2456R,D2456I)
       E012345=E0fin(p1sq,p2sq,p3sq,p4sq,s56,s12,s23,s34,s123,s234,D02
     &   345,D01345,D01245,D01235,D01234)
       E012345R=Dble(E012345)
       E012345I=DImag(E012345)
       call tens_red5_new_Re_Com(p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
     &   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2345R,
     &   D1345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01235I,
     &   D01234I,D2345I,D1345I,D1245I,D1235I,D1234I,E12345R,E12345I)
       E012346=E0fin(p1sq,p2sq,p3sq,s45,p6sq,s12,s23,s345,s123,s16,D02
     &   346,D01346,D01246,D01236,D01234)
       E012346R=Dble(E012346)
       E012346I=DImag(E012346)
       call tens_red5_new_Re_Com(p1sq,p2sq,p3sq,s45,p1p2,p1p3,p1p4+p1p
     &   5,p2p3,p2p4+p2p5,p3p4+p3p5,D02346R,D01346R,D01246R,D01236R,D
     &   01234R,D2346R,D1346R,D1246R,D1236R,D1234R,D02346I,D01346I,D0
     &   1246I,D01236I,D01234I,D2346I,D1346I,D1246I,D1236I,D1234I,E12
     &   346R,E12346I)
       E012356=E0fin(p1sq,p2sq,s34,p5sq,p6sq,s12,s234,s345,s56,s16,D02
     &   356,D01356,D01256,D01236,D01235)
       E012356R=Dble(E012356)
       E012356I=DImag(E012356)
       call tens_red5_new_Re_Com(p1sq,p2sq,s34,p5sq,p1p2,p1p3+p1p4,p1p
     &   5,p2p3+p2p4,p2p5,p3p5+p4p5,D02356R,D01356R,D01256R,D01236R,D
     &   01235R,D2356R,D1356R,D1256R,D1236R,D1235R,D02356I,D01356I,D0
     &   1256I,D01236I,D01235I,D2356I,D1356I,D1256I,D1236I,D1235I,E12
     &   356R,E12356I)
       E012456=E0fin(p1sq,s23,p4sq,p5sq,p6sq,s123,s234,s45,s56,s16,D02
     &   456,D01456,D01256,D01246,D01245)
       E012456R=Dble(E012456)
       E012456I=DImag(E012456)
       call tens_red5_new_Re_Com(p1sq,s23,p4sq,p5sq,p1p2+p1p3,p1p4,p1p
     &   5,p2p4+p3p4,p2p5+p3p5,p4p5,D02456R,D01456R,D01256R,D01246R,D
     &   01245R,D2456R,D1456R,D1256R,D1246R,D1245R,D02456I,D01456I,D0
     &   1256I,D01246I,D01245I,D2456I,D1456I,D1256I,D1246I,D1245I,E12
     &   456R,E12456I)
       E023456=E0fin(p2sq,p3sq,p4sq,p5sq,s16,s23,s34,s45,s234,s345,D03
     &   456,D02456,D02356,D02346,D02345)
       E023456R=Dble(E023456)
       E023456I=DImag(E023456)
       call tens_red5_new_Re_Com(p2sq,p3sq,p4sq,p5sq,p2p3,p2p4,p2p5,p3
     &   p4,p3p5,p4p5,D03456R,D02456R,D02356R,D02346R,D02345R,D3456R,
     &   D2456R,D2356R,D2346R,D2345R,D03456I,D02456I,D02356I,D02346I,
     &   D02345I,D3456I,D2456I,D2356I,D2346I,D2345I,E23456R,E23456I)
       E013456=E0fin(s12,p3sq,p4sq,p5sq,p6sq,s123,s34,s45,s56,s345,D03
     &   456,D01456,D01356,D01346,D01345)
       E013456R=Dble(E013456)
       E013456I=DImag(E013456)
       call tens_red5_new_Re_Com(s12,p3sq,p4sq,p5sq,p1p3+p2p3,p1p4+p2p
     &   4,p1p5+p2p5,p3p4,p3p5,p4p5,D03456R,D01456R,D01356R,D01346R,D
     &   01345R,D3456R,D1456R,D1356R,D1346R,D1345R,D03456I,D01456I,D0
     &   1356I,D01346I,D01345I,D3456I,D1456I,D1356I,D1346I,D1345I,E13
     &   456R,E13456I)
       F0123456=F0fin(p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45,s5
     &   6,s16,s123,s234,s345,E023456,E013456,E012456,E012356,E012346
     &   ,E012345)
       F0123456R=Dble(F0123456)
       F0123456I=DImag(F0123456)
       call tens_red6_new_Re_Com(p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1
     &   p4,p1p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,E023456R,E013456R,E012
     &   456R,E012356R,E012346R,E012345R,E23456R,E13456R,E12456R,E123
     &   56R,E12346R,E12345R,E023456I,E013456I,E012456I,E012356I,E012
     &   346I,E012345I,E23456I,E13456I,E12456I,E12356I,E12346I,E12345
     &   I,F123456R,F123456I,F123456)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call HexBoxCross2FFhl1(F(1))
       call HexBoxCross2FFhl2(F(93))
       call HexBoxCross2FFhl3(F(185))
       call HexBoxCross2FFhl4(F(277))
       call HexBoxCross2FFhl5(F(369))
       call HexBoxCross2FFhl6(F(461))
       call HexBoxCross2FFhl7(F(553))
       call HexBoxCross2FFhl8(F(645))
       call HexBoxCross2FFhl9(F(737))
       call HexBoxCross2FFhl10(F(829))
       call HexBoxCross2FFhl11(F(921))
       call HexBoxCross2FFhl12(F(1013))
       call HexBoxCross2FFhl13(F(1105))
       call HexBoxCross2FFhl14(F(1197))
       call HexBoxCross2FFhl15(F(1289))
       call HexBoxCross2FFhl16(F(1381))
       call HexBoxCross2FFhl17(F(1473))
       call HexBoxCross2FFhl18(F(1565))
       call HexBoxCross2FFhl19(F(1657))
       call HexBoxCross2FFhl20(F(1749))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup3 = dotrc(p1,mup3)
       p1mup6 = dotrc(p1,mup6)
       p2mup3 = dotrc(p2,mup3)
       p2mup6 = dotrc(p2,mup6)
       p3mup3 = dotrc(p3,mup3)
       p3mup6 = dotrc(p3,mup6)
       p4mup3 = dotrc(p4,mup3)
       p4mup6 = dotrc(p4,mup6)
       p5mup3 = dotrc(p5,mup3)
       p5mup6 = dotrc(p5,mup6)
       p6mup3 = dotrc(p6,mup3)
       p6mup6 = dotrc(p6,mup6)
       mup3mup6 = dotcc(mup3,mup6)
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call HexBoxCross2Fa1(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,
     &   p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(1))
       call HexBoxCross2Fa2(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,
     &   p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(20))
       call HexBoxCross2Fa3(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,
     &   p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(39))
       call HexBoxCross2Fa4(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,
     &   p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(58))
       call HexBoxCross2Fa5(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,
     &   p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(77))
       call HexBoxCross2Fa6(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,
     &   p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(96))
       call HexBoxCross2Fa7(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,
     &   p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(115))
       call HexBoxCross2Fa8(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,
     &   p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(134))
       call HexBoxCross2Fa9(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,
     &   p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(153))
       call HexBoxCross2Fa10(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(172))
       call HexBoxCross2Fa11(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(191))
       call HexBoxCross2Fa12(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(210))
       call HexBoxCross2Fa13(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(229))
       call HexBoxCross2Fa14(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(248))
       call HexBoxCross2Fa15(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(267))
       call HexBoxCross2Fa16(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(286))
       call HexBoxCross2Fa17(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(305))
       call HexBoxCross2Fa18(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(324))
       call HexBoxCross2Fa19(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(343))
       call HexBoxCross2Fa20(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(362))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the Matrix Element  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       barpsi_p5_P(1)=barpsi_p5(1)
       barpsi_p5_P(2)=barpsi_p5(2)
       barpsi_p5_M(1)=barpsi_p5(3)
       barpsi_p5_M(2)=barpsi_p5(4)
       psi_p4_M(1)=psi_p4(1)
       psi_p4_M(2)=psi_p4(2)
       psi_p4_P(1)=psi_p4(3)
       psi_p4_P(2)=psi_p4(4)
       psi_p1_M(1)=psi_p1(1)
       psi_p1_M(2)=psi_p1(2)
       psi_p1_P(1)=psi_p1(3)
       psi_p1_P(2)=psi_p1(4)
       barpsi_p2_P(1)=barpsi_p2(1)
       barpsi_p2_P(2)=barpsi_p2(2)
       barpsi_p2_M(1)=barpsi_p2(3)
       barpsi_p2_M(2)=barpsi_p2(4)
c************************************************************************************
c************************************************************************************
       SMB0(1) = SC1c(barpsi_p2_P,mup3,psi_p4_P,+1)
       SMB0(2) = SC1c(barpsi_p5_P,mup6,psi_p1_P,+1)
       SMB0(3) = SC1r(barpsi_p5_P,p2,psi_p1_P,+1)
       SMB0(4) = SC1r(barpsi_p5_P,p4,psi_p1_P,+1)
       SMB0(5) = SC1r(barpsi_p5_P,p6,psi_p1_P,+1)
       SMB0(6) = SC3crr(barpsi_p5_P,mup6,p4,p2,psi_p1_P,+1)
       SMB0(7) = SC3crr(barpsi_p5_P,mup6,p6,p2,psi_p1_P,+1)
       SMB0(8) = SC3crr(barpsi_p5_P,mup6,p6,p4,psi_p1_P,+1)
       SMB0(9) = SC3rrr(barpsi_p5_P,p6,p4,p2,psi_p1_P,+1)
       SMB0(10) = SC1c(barpsi_p2_P,mup6,psi_p4_P,+1)
       SMB0(11) = SC1c(barpsi_p5_P,mup3,psi_p1_P,+1)
       SMB0(12) = SC3crr(barpsi_p5_P,mup3,p4,p2,psi_p1_P,+1)
       SMB0(13) = SC3crr(barpsi_p5_P,mup3,p6,p2,psi_p1_P,+1)
       SMB0(14) = SC3crr(barpsi_p5_P,mup3,p6,p4,psi_p1_P,+1)
       SMB0(15) = SC1r(barpsi_p2_P,p6,psi_p4_P,+1)
       SMB0(16) = SC3ccr(barpsi_p5_P,mup6,mup3,p4,psi_p1_P,+1)
       SMB0(17) = SC1r(barpsi_p2_P,p1,psi_p4_P,+1)
       SMB0(18) = SC3ccr(barpsi_p5_P,mup6,mup3,p2,psi_p1_P,+1)
       SMB0(19) = SC3ccr(barpsi_p5_P,mup6,mup3,p6,psi_p1_P,+1)
       SMB0(20) = SC5ccrrr(barpsi_p5_P,mup6,mup3,p6,p4,p2,psi_p1_P,+1)
       SMB0(21) = SC1r(barpsi_p2_P,p5,psi_p4_P,+1)
       SMB0(22) = SC3crr(barpsi_p2_P,mup6,p5,p1,psi_p4_P,+1)
       SMB0(23) = SC3crr(barpsi_p2_P,mup6,p6,p1,psi_p4_P,+1)
       SMB0(24) = SC3crr(barpsi_p2_P,mup6,p6,p5,psi_p4_P,+1)
       SMB0(25) = SC3rrr(barpsi_p2_P,p6,p5,p1,psi_p4_P,+1)
       SMB0(26) = SC3crr(barpsi_p2_P,mup3,p5,p1,psi_p4_P,+1)
       SMB0(27) = SC3crr(barpsi_p2_P,mup3,p6,p1,psi_p4_P,+1)
       SMB0(28) = SC3crr(barpsi_p2_P,mup3,p6,p5,psi_p4_P,+1)
       SMB0(29) = SC3ccr(barpsi_p2_P,mup6,mup3,p1,psi_p4_P,+1)
       SMB0(30) = SC3ccr(barpsi_p2_P,mup6,mup3,p5,psi_p4_P,+1)
       SMB0(31) = SC3ccr(barpsi_p2_P,mup6,mup3,p6,psi_p4_P,+1)
       SMB0(32) = SC5ccrrr(barpsi_p2_P,mup6,mup3,p6,p5,p1,psi_p4_P,+1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB11(i) = SC1c(barpsi_p2_P,v1,psi_p4_P,+1)
       SMB12(i) = SC1c(barpsi_p5_P,v1,psi_p1_P,+1)
       SMB13(i) = SC3crc(barpsi_p2_P,mup3,p6,v1,psi_p4_P,+1)
       SMB14(i) = SC3ccc(barpsi_p2_P,mup6,mup3,v1,psi_p4_P,+1)
       SMB15(i) = SC3crc(barpsi_p2_P,mup6,p6,v1,psi_p4_P,+1)
       SMB16(i) = SC3ccc(barpsi_p5_P,mup6,mup3,v1,psi_p1_P,+1)
       SMB17(i) = SC3crc(barpsi_p5_P,mup3,p6,v1,psi_p1_P,+1)
       SMB18(i) = SC3crc(barpsi_p5_P,mup6,p6,v1,psi_p1_P,+1)
      enddo
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       SMB1(1) = SMB0(1)*SMB0(2)
       SMB1(2) = SMB0(1)*SMB0(3)
       SMB1(3) = SMB0(1)*SMB0(4)
       SMB1(4) = SMB0(1)*SMB0(5)
       SMB1(5) = SMB0(1)*SMB0(6)
       SMB1(6) = SMB0(1)*SMB0(7)
       SMB1(7) = SMB0(1)*SMB0(8)
       SMB1(8) = SMB0(1)*SMB0(9)
       SMB1(9) = SMB0(10)*SMB0(11)
       SMB1(10) = SMB0(3)*SMB0(10)
       SMB1(11) = SMB0(4)*SMB0(10)
       SMB1(12) = SMB0(5)*SMB0(10)
       SMB1(13) = SMB0(10)*SMB0(12)
       SMB1(14) = SMB0(10)*SMB0(13)
       SMB1(15) = SMB0(10)*SMB0(14)
       SMB1(16) = SMB0(9)*SMB0(10)
       SMB1(21) = SMB0(4)*SMB0(15)
       SMB1(22) = SMB0(15)*SMB0(16)
       SMB1(24) = SMB0(11)*SMB0(17)
       SMB1(25) = SMB0(2)*SMB0(17)
       SMB1(26) = SMB0(3)*SMB0(17)
       SMB1(27) = SMB0(4)*SMB0(17)
       SMB1(28) = SMB0(5)*SMB0(17)
       SMB1(29) = SMB0(12)*SMB0(17)
       SMB1(30) = SMB0(13)*SMB0(17)
       SMB1(31) = SMB0(14)*SMB0(17)
       SMB1(32) = SMB0(17)*SMB0(18)
       SMB1(33) = SMB0(16)*SMB0(17)
       SMB1(34) = SMB0(17)*SMB0(19)
       SMB1(35) = SMB0(6)*SMB0(17)
       SMB1(36) = SMB0(7)*SMB0(17)
       SMB1(37) = SMB0(8)*SMB0(17)
       SMB1(38) = SMB0(9)*SMB0(17)
       SMB1(39) = SMB0(17)*SMB0(20)
       SMB1(40) = SMB0(11)*SMB0(21)
       SMB1(41) = SMB0(2)*SMB0(21)
       SMB1(42) = SMB0(3)*SMB0(21)
       SMB1(43) = SMB0(4)*SMB0(21)
       SMB1(44) = SMB0(5)*SMB0(21)
       SMB1(45) = SMB0(12)*SMB0(21)
       SMB1(46) = SMB0(13)*SMB0(21)
       SMB1(47) = SMB0(14)*SMB0(21)
       SMB1(48) = SMB0(18)*SMB0(21)
       SMB1(49) = SMB0(16)*SMB0(21)
       SMB1(50) = SMB0(19)*SMB0(21)
       SMB1(51) = SMB0(6)*SMB0(21)
       SMB1(52) = SMB0(7)*SMB0(21)
       SMB1(53) = SMB0(8)*SMB0(21)
       SMB1(54) = SMB0(9)*SMB0(21)
       SMB1(55) = SMB0(20)*SMB0(21)
       SMB1(56) = SMB0(11)*SMB0(15)
       SMB1(57) = SMB0(2)*SMB0(15)
       SMB1(58) = SMB0(3)*SMB0(15)
       SMB1(59) = SMB0(5)*SMB0(15)
       SMB1(60) = SMB0(12)*SMB0(15)
       SMB1(61) = SMB0(13)*SMB0(15)
       SMB1(62) = SMB0(14)*SMB0(15)
       SMB1(63) = SMB0(15)*SMB0(18)
       SMB1(64) = SMB0(15)*SMB0(19)
       SMB1(65) = SMB0(6)*SMB0(15)
       SMB1(66) = SMB0(7)*SMB0(15)
       SMB1(67) = SMB0(8)*SMB0(15)
       SMB1(68) = SMB0(9)*SMB0(15)
       SMB1(69) = SMB0(15)*SMB0(20)
       SMB1(70) = SMB0(11)*SMB0(22)
       SMB1(71) = SMB0(11)*SMB0(23)
       SMB1(72) = SMB0(11)*SMB0(24)
       SMB1(73) = SMB0(11)*SMB0(25)
       SMB1(74) = SMB0(2)*SMB0(26)
       SMB1(75) = SMB0(2)*SMB0(27)
       SMB1(76) = SMB0(2)*SMB0(28)
       SMB1(77) = SMB0(2)*SMB0(25)
       SMB1(78) = SMB0(5)*SMB0(29)
       SMB1(79) = SMB0(3)*SMB0(26)
       SMB1(80) = SMB0(3)*SMB0(27)
       SMB1(81) = SMB0(3)*SMB0(28)
       SMB1(82) = SMB0(3)*SMB0(29)
       SMB1(83) = SMB0(3)*SMB0(30)
       SMB1(84) = SMB0(3)*SMB0(31)
       SMB1(85) = SMB0(3)*SMB0(22)
       SMB1(86) = SMB0(3)*SMB0(23)
       SMB1(87) = SMB0(3)*SMB0(24)
       SMB1(88) = SMB0(3)*SMB0(25)
       SMB1(89) = SMB0(3)*SMB0(32)
       SMB1(90) = SMB0(4)*SMB0(26)
       SMB1(91) = SMB0(4)*SMB0(27)
       SMB1(92) = SMB0(4)*SMB0(28)
       SMB1(93) = SMB0(4)*SMB0(29)
       SMB1(94) = SMB0(4)*SMB0(30)
       SMB1(95) = SMB0(4)*SMB0(31)
       SMB1(96) = SMB0(4)*SMB0(22)
       SMB1(97) = SMB0(4)*SMB0(23)
       SMB1(98) = SMB0(4)*SMB0(24)
       SMB1(99) = SMB0(4)*SMB0(25)
       SMB1(100) = SMB0(4)*SMB0(32)
       SMB1(101) = SMB0(5)*SMB0(26)
       SMB1(102) = SMB0(5)*SMB0(27)
       SMB1(103) = SMB0(5)*SMB0(28)
       SMB1(104) = SMB0(5)*SMB0(30)
       SMB1(105) = SMB0(5)*SMB0(31)
       SMB1(106) = SMB0(5)*SMB0(22)
       SMB1(107) = SMB0(5)*SMB0(23)
       SMB1(108) = SMB0(5)*SMB0(24)
       SMB1(109) = SMB0(5)*SMB0(25)
       SMB1(110) = SMB0(5)*SMB0(32)
       SMB1(112) = SMB0(7)*SMB0(26)
       SMB1(113) = SMB0(8)*SMB0(26)
       SMB1(115) = SMB0(6)*SMB0(27)
       SMB1(116) = SMB0(7)*SMB0(27)
       SMB1(117) = SMB0(8)*SMB0(27)
       SMB1(118) = SMB0(6)*SMB0(28)
       SMB1(119) = SMB0(7)*SMB0(28)
       SMB1(120) = SMB0(8)*SMB0(28)
       SMB1(121) = SMB0(9)*SMB0(29)
       SMB1(122) = SMB0(9)*SMB0(30)
       SMB1(123) = SMB0(9)*SMB0(31)
       SMB1(124) = SMB0(13)*SMB0(22)
       SMB1(125) = SMB0(14)*SMB0(22)
       SMB1(126) = SMB0(13)*SMB0(23)
       SMB1(127) = SMB0(14)*SMB0(23)
       SMB1(128) = SMB0(13)*SMB0(24)
       SMB1(129) = SMB0(14)*SMB0(24)
       SMB1(130) = SMB0(16)*SMB0(25)
       SMB1(131) = SMB0(18)*SMB0(25)
       SMB1(132) = SMB0(19)*SMB0(25)
c************************************************************************************
       SMB1(17) = dotcc(SMB11,SMB12)
       SMB1(18) = dotcc(SMB12,SMB13)
       SMB1(19) = dotcc(SMB12,SMB14)
       SMB1(20) = dotcc(SMB12,SMB15)
       SMB1(23) = dotcc(SMB11,SMB16)
       SMB1(111) = dotcc(SMB11,SMB17)
       SMB1(114) = dotcc(SMB11,SMB18)
       SMB1(133) = dotcc(SMB13,SMB18)
c************************************************************************************
c************************************************************************************
       SMB(1) = SMB1(3)
       SMB(2) = SMB1(4)
       SMB(3) = SMB1(27)
       SMB(4) = SMB1(28)
       SMB(5) = SMB1(93)
       SMB(6) = SMB1(78)
       SMB(7) = SMB1(2)
       SMB(8) = SMB1(10)
       SMB(9) = SMB1(12)
       SMB(10) = SMB1(26)
       SMB(11) = SMB1(82)
       SMB(12) = SMB1(11)
       SMB(13) = SMB1(43)
       SMB(14) = SMB1(44)
       SMB(15) = SMB1(21)
       SMB(16) = SMB1(59)
       SMB(17) = SMB1(94)
       SMB(18) = SMB1(95)
       SMB(19) = SMB1(104)
       SMB(20) = SMB1(105)
       SMB(21) = SMB1(42)
       SMB(22) = SMB1(58)
       SMB(23) = SMB1(83)
       SMB(24) = SMB1(84)
       SMB(25) = SMB1(1)
       SMB(26) = SMB1(9)
       SMB(27) = SMB1(17)
       SMB(28) = 2*mup3mup6*SMB1(17)-SMB1(23)
       SMB(29) = SMB1(25)
       SMB(30) = SMB1(41)
       SMB(31) = SMB1(57)
       SMB(32) = SMB1(74)
       SMB(33) = SMB1(76)
       SMB(34) = 2*mup3mup6*SMB1(17)-SMB1(19)
       SMB(35) = SMB1(24)
       SMB(36) = SMB1(33)
       SMB(37) = SMB1(40)
       SMB(38) = SMB1(49)
       SMB(39) = SMB1(56)
       SMB(40) = SMB1(22)
       SMB(41) = SMB1(70)
       SMB(42) = SMB1(72)
       SMB(43) = (p6sq-s16+s234-s56)*SMB1(17)
       SMB(44) = 0
       SMB(45) = 2*mup3mup6*(p6sq-s16+s234-s56)*SMB1(17)+(-p6sq+s16-s2
     &   34+s56)*SMB1(19)+4*p5mup6*SMB1(24)-4*p5mup3*SMB1(25)
       SMB(46) = -2*p5mup6*SMB1(111)+2*p5mup3*SMB1(114)
       SMB(47) = -2*SMB1(11)+2*p4mup6*SMB1(17)
       SMB(48) = 2*p1mup3*SMB1(17)-2*SMB1(24)
       SMB(49) = 0
       SMB(50) = SMB1(18)
       SMB(51) = -2*(2*p1mup6*SMB1(3)-2*p1mup3*p4mup6*SMB1(17)+2*p4mup
     &   6*SMB1(24)-2*mup3mup6*SMB1(27)+SMB1(93))
       SMB(52) = -2*(2*p5mup6*SMB1(3)-2*mup3mup6*SMB1(43)+SMB1(94))
       SMB(53) = -2*(2*p6mup6*SMB1(3)-2*mup3mup6*SMB1(21)+SMB1(95)+p4m
     &   up6*(-2*p6mup3*SMB1(17)+SMB1(111)))
       SMB(54) = SMB1(34)
       SMB(55) = (-4*p1mup6*p6mup3+4*p1mup3*p6mup6+2*mup3mup6*(p6sq-s1
     &   6))*SMB1(17)+(-p6sq+s16)*SMB1(19)+(-p6sq+s16)*SMB1(23)-4*p6m
     &   up6*SMB1(24)+4*mup3mup6*SMB1(28)-2*SMB1(78)+2*p1mup6*SMB1(11
     &   1)
       SMB(56) = SMB1(32)
       SMB(57) = s12*SMB1(23)-4*p2mup3*SMB1(25)
       SMB(58) = 0
       SMB(59) = 2*p6mup6*SMB1(17)-SMB1(20)
       SMB(60) = SMB1(48)
       SMB(61) = SMB1(50)
       SMB(62) = SMB1(63)
       SMB(63) = SMB1(64)
       SMB(64) = SMB1(75)
       SMB(65) = 4*mup3mup6*SMB1(17)
       SMB(66) = -4*p2mup3*p5mup6*SMB1(17)+(s16-s234+s34-s345)*SMB1(23
     &   )
       SMB(67) = -4*p5mup6*p6mup3*SMB1(17)+(-p6sq+s56)*SMB1(23)+2*p5mu
     &   p6*SMB1(111)
       SMB(68) = -((s12+s16-s345)*SMB1(23))-2*p2mup3*SMB1(114)
       SMB(69) = p6sq*SMB1(23)
       SMB(70) = SMB1(7)
       SMB(71) = SMB1(5)
       SMB(72) = SMB1(71)
       SMB(73) = SMB1(73)
       SMB(74) = SMB1(77)
       SMB(75) = (-p6sq+s16)*SMB1(17)-2*SMB1(28)
       SMB(76) = SMB1(99)
       SMB(77) = SMB1(109)
       SMB(78) = (-p6sq+s16)*SMB1(23)-2*SMB1(34)
       SMB(79) = SMB1(130)
       SMB(80) = SMB1(132)
       SMB(81) = SMB1(88)
       SMB(82) = SMB1(131)
       SMB(83) = SMB1(14)
       SMB(84) = 0
       SMB(85) = SMB1(30)
       SMB(86) = SMB1(46)
       SMB(87) = SMB1(61)
       SMB(88) = 2*p1mup6*SMB1(17)-2*SMB1(25)
       SMB(89) = 0
       SMB(90) = SMB1(20)
       SMB(91) = SMB1(85)
       SMB(92) = SMB1(86)
       SMB(93) = SMB1(87)
       SMB(94) = SMB1(106)
       SMB(95) = SMB1(107)
       SMB(96) = SMB1(108)
       SMB(97) = 2*mup3mup6*s12*SMB1(17)-s12*SMB1(23)-4*p2mup6*SMB1(24
     &   )
       SMB(98) = -4*p2mup6*p5mup3*SMB1(17)+(s16-s234+s34-s345)*(2*mup3
     &   mup6*SMB1(17)-SMB1(23))
       SMB(99) = SMB1(124)
       SMB(100) = -2*mup3mup6*(s12+s16-s345)*SMB1(17)+(s12+s16-s345)*S
     &   MB1(23)-2*p2mup6*SMB1(111)
       SMB(101) = SMB1(126)
       SMB(102) = SMB1(128)
       SMB(103) = SMB1(6)
       SMB(104) = 0
       SMB(105) = 0
       SMB(106) = -4*p2mup6*SMB1(3)+4*p2mup3*SMB1(11)
       SMB(107) = (-4*p2mup6*p6mup3+4*p2mup3*p6mup6)*SMB1(17)+2*p2mup6
     &   *SMB1(111)-2*p2mup3*SMB1(114)
       SMB(108) = SMB1(15)
       SMB(109) = (s123-s45-s56)*SMB1(17)-2*SMB1(21)
       SMB(110) = 2*mup3mup6*(s123-s45-s56)*SMB1(17)+(-s123+s45+s56)*S
     &   MB1(23)-2*SMB1(95)
       SMB(111) = (p6sq-s16+s234-s56)*SMB1(19)
       SMB(112) = (-p6sq+s16)*SMB1(19)-2*SMB1(78)
       SMB(113) = 4*p1mup6*SMB1(17)
       SMB(114) = 4*p5mup6*SMB1(17)
       SMB(115) = 2*(SMB1(20)+SMB1(114))
       SMB(116) = 2*(4*p1mup3*p5mup6+2*mup3mup6*(-p6sq+s16-s234+s56))*
     &   SMB1(17)+2*(p6sq-s16+s234-s56)*SMB1(19)-8*p5mup6*SMB1(24)+4*
     &   SMB1(74)
       SMB(117) = (-8*p1mup6*p6mup3+8*p1mup3*p6mup6+4*mup3mup6*(p6sq-s
     &   16))*SMB1(17)+4*p1mup6*SMB1(18)+2*(-p6sq+s16)*SMB1(19)-8*p6m
     &   up6*SMB1(24)+8*mup3mup6*SMB1(28)+4*SMB1(75)-4*SMB1(78)+4*p1m
     &   up6*SMB1(111)
       SMB(118) = 0
       SMB(119) = 4*p5mup6*SMB1(18)
       SMB(120) = -2*SMB1(3)+2*p4mup3*SMB1(17)
       SMB(121) = 2*p6mup3*SMB1(17)-SMB1(18)
       SMB(122) = 0
       SMB(123) = -2*SMB1(22)+(s123-s45-s56)*SMB1(23)
       SMB(124) = 4*p2mup3*SMB1(17)
       SMB(125) = 4*p4mup3*SMB1(17)
       SMB(126) = 8*p6mup3*SMB1(17)-2*(SMB1(18)+SMB1(111))
       SMB(127) = p2mup3*(8*p6mup6*SMB1(17)-4*SMB1(20))
       SMB(128) = -8*p6mup6*SMB1(3)+4*SMB1(7)+8*p4mup6*p6mup3*SMB1(17)
     &   +8*mup3mup6*SMB1(21)+2*(-s123+s45+s56)*SMB1(23)-4*SMB1(95)-4
     &   *p4mup6*SMB1(111)+4*p4mup3*SMB1(114)
       SMB(129) = SMB1(35)
       SMB(130) = SMB1(36)
       SMB(131) = SMB1(51)
       SMB(132) = SMB1(52)
       SMB(133) = SMB1(65)
       SMB(134) = SMB1(66)
       SMB(135) = SMB1(79)
       SMB(136) = SMB1(80)
       SMB(137) = SMB1(81)
       SMB(138) = SMB1(90)
       SMB(139) = SMB1(91)
       SMB(140) = SMB1(92)
       SMB(141) = SMB1(101)
       SMB(142) = SMB1(102)
       SMB(143) = SMB1(103)
       SMB(144) = SMB1(112)
       SMB(145) = SMB1(115)
       SMB(146) = SMB1(116)
       SMB(147) = SMB1(118)
       SMB(148) = SMB1(119)
       SMB(149) = SMB1(8)
       SMB(150) = SMB1(13)
       SMB(151) = SMB1(16)
       SMB(152) = SMB1(38)
       SMB(153) = SMB1(54)
       SMB(154) = SMB1(68)
       SMB(155) = SMB1(121)
       SMB(156) = SMB1(122)
       SMB(157) = SMB1(123)
       SMB(158) = SMB1(29)
       SMB(159) = SMB1(31)
       SMB(160) = SMB1(45)
       SMB(161) = SMB1(47)
       SMB(162) = SMB1(60)
       SMB(163) = SMB1(62)
       SMB(164) = SMB1(96)
       SMB(165) = SMB1(97)
       SMB(166) = SMB1(98)
       SMB(167) = -4*p1mup3*SMB1(11)+4*p1mup6*p4mup3*SMB1(17)-4*p4mup3
     &   *SMB1(25)+2*SMB1(93)
       SMB(168) = 2*(2*p1mup6*p6mup3-2*p1mup3*p6mup6+mup3mup6*(-p6sq+s
     &   16))*SMB1(17)+(p6sq-s16)*SMB1(19)+(p6sq-s16)*SMB1(23)-4*p6mu
     &   p3*SMB1(25)+2*SMB1(78)+2*p1mup3*SMB1(114)
       SMB(169) = -4*p5mup3*SMB1(11)+2*SMB1(94)
       SMB(170) = -2*(2*p5mup3*p6mup6+mup3mup6*(p6sq-s56))*SMB1(17)+(p
     &   6sq-s56)*SMB1(23)+2*p5mup3*SMB1(114)
       SMB(171) = SMB1(125)
       SMB(172) = 2*(-2*p6mup3*SMB1(11)+2*p4mup3*p6mup6*SMB1(17)+SMB1(
     &   95)-p4mup3*SMB1(114))
       SMB(173) = p6sq*(2*mup3mup6*SMB1(17)-SMB1(23))
       SMB(174) = SMB1(127)
       SMB(175) = SMB1(129)
       SMB(176) = SMB1(37)
       SMB(177) = SMB1(113)
       SMB(178) = SMB1(117)
       SMB(179) = SMB1(53)
       SMB(180) = SMB1(67)
       SMB(181) = SMB1(120)
       SMB(182) = 2*p6mup3*s12*SMB1(17)-4*p2mup3*SMB1(28)-s12*SMB1(111
     &   )
       SMB(183) = 2*((p6sq-s16)*SMB1(3)+p1mup3*(s123-s45-s56)*SMB1(17)
     &   +(-s123+s45+s56)*SMB1(24)+SMB1(91))
       SMB(184) = 2*p2mup3*(p6sq-s56)*SMB1(17)+(s16-s234+s34-s345)*(2*
     &   p6mup3*SMB1(17)-SMB1(111))
       SMB(185) = 2*((p6sq-s56)*SMB1(3)+SMB1(92))
       SMB(186) = -2*(p2mup3*p6sq+p6mup3*(s12+s16-s345))*SMB1(17)+(s12
     &   +s16-s345)*SMB1(111)
       SMB(187) = -2*p6sq*SMB1(3)+(s123-s45-s56)*(2*p6mup3*SMB1(17)-SM
     &   B1(111))
       SMB(188) = 8*p6mup6*SMB1(17)-2*(SMB1(20)+SMB1(114))
       SMB(189) = 2*p6mup6*s12*SMB1(17)-4*p2mup6*SMB1(28)-s12*SMB1(114
     &   )
       SMB(190) = 2*((p6sq-s16)*SMB1(11)+p1mup6*(s123-s45-s56)*SMB1(17
     &   )+(-s123+s45+s56)*SMB1(25)+SMB1(97))
       SMB(191) = 2*p2mup6*(p6sq-s56)*SMB1(17)+(s16-s234+s34-s345)*(2*
     &   p6mup6*SMB1(17)-SMB1(114))
       SMB(192) = 2*((p6sq-s56)*SMB1(11)+SMB1(98))
       SMB(193) = -2*(p2mup6*p6sq+p6mup6*(s12+s16-s345))*SMB1(17)+(s12
     &   +s16-s345)*SMB1(114)
       SMB(194) = -2*p6sq*SMB1(11)+(s123-s45-s56)*(2*p6mup6*SMB1(17)-S
     &   MB1(114))
       SMB(195) = 2*(-p6sq+s16)*SMB1(17)
       SMB(196) = s12*(-p6sq+s56)*SMB1(17)+2*(-s16+s234-s34+s345)*SMB1
     &   (28)
       SMB(197) = -((p6sq-s16+s234-s56)*(-s123+s45+s56)*SMB1(17))-2*(p
     &   6sq-s56)*SMB1(27)+2*(p6sq-s16)*SMB1(43)-2*SMB1(99)
       SMB(198) = p6sq*s12*SMB1(17)+2*(s12+s16-s345)*SMB1(28)
       SMB(199) = -((p6sq-s16)*(s123-s45-s56)*SMB1(17))+2*(p6sq-s16)*S
     &   MB1(21)+2*p6sq*SMB1(27)+2*(-s123+s45+s56)*SMB1(28)
       SMB(200) = 2*(p6sq-s16)*(-2*mup3mup6*SMB1(17)+SMB1(23))
       SMB(201) = 2*(p2mup6*(-p6sq+s16-s234+s56)*SMB1(18)+p2mup3*(p6sq
     &   -s16+s234-s56)*SMB1(20)+(-s16+s234-s34+s345)*SMB1(78)+(-p6sq
     &   +s56)*SMB1(82))
       SMB(202) = 2*(p2mup6*(p6sq-s16)*SMB1(18)+p2mup3*(-p6sq+s16)*SMB
     &   1(20)+(s12+s16-s345)*SMB1(78)+p6sq*SMB1(82))
       SMB(203) = -((p6sq-s16+s234-s56)*(-s123+s45+s56)*SMB1(19))-2*s4
     &   5*SMB1(78)-2*(p6sq-s56)*SMB1(93)-2*(p6sq-s16+s234-s56)*SMB1(
     &   95)
       SMB(204) = -((p6sq-s16)*(s123-s45-s56)*SMB1(19))+2*(-s123+s45+s
     &   56)*SMB1(78)+2*p6sq*SMB1(93)+2*(p6sq-s16)*SMB1(95)
       SMB(205) = (p6sq-s16+s234-s56)*SMB1(18)
       SMB(206) = SMB1(89)
       SMB(207) = SMB1(110)
       SMB(208) = 2*(-s16+s234-s34+s345)*SMB1(25)+4*p5mup6*SMB1(26)
       SMB(209) = 2*p6mup6*(p6sq-s16+s234-s56)*SMB1(17)+(-p6sq+s16-s23
     &   4+s56)*SMB1(20)+2*(p6sq-s56)*SMB1(25)+4*p5mup6*SMB1(28)
       SMB(210) = 2*((s12+s16-s345)*SMB1(25)+SMB1(36))
       SMB(211) = 2*p6mup6*(-p6sq+s16)*SMB1(17)+(p6sq-s16)*SMB1(20)-2*
     &   p6sq*SMB1(25)
       SMB(212) = 2*p5mup6*(s12+s16-s345)*SMB1(17)+(s16-s234+s34-s345)
     &   *(2*p6mup6*SMB1(17)-SMB1(20))
       SMB(213) = -2*(p5mup6*p6sq+p6mup6*(p6sq-s56))*SMB1(17)+(p6sq-s5
     &   6)*SMB1(20)
       SMB(214) = 2*(2*p2mup3*p6mup6+mup3mup6*(s12+s16-s345))*(-p6sq+s
     &   16-s234+s56)*SMB1(17)+(s12+s16-s345)*(p6sq-s16+s234-s56)*SMB
     &   1(19)+2*p2mup3*(p6sq-s16+s234-s56)*SMB1(20)-2*(s16-s234+s34-
     &   s345)*SMB1(75)+4*p5mup6*SMB1(80)
       SMB(215) = 2*(p6sq-s56)*SMB1(75)+4*p5mup6*SMB1(102)+(p6sq-s16+s
     &   234-s56)*SMB1(133)
       SMB(216) = SMB1(100)
       SMB(217) = 2*(-p6sq+s16-s234+s56)*SMB1(11)+2*p4mup6*(p6sq-s16+s
     &   234-s56)*SMB1(17)-2*s45*SMB1(25)+4*p5mup6*SMB1(27)
       SMB(218) = 2*((p6sq-s16)*SMB1(11)+p4mup6*(-p6sq+s16)*SMB1(17)+(
     &   -s123+s45+s56)*SMB1(25)+SMB1(37))
       SMB(219) = 2*p6mup6*s45*SMB1(17)-s45*SMB1(20)-4*p5mup6*SMB1(21)
     &   
       SMB(220) = -4*p6mup6*(p6sq-s16+s234-s56)*SMB1(3)+2*p4mup6*(p6sq
     &   -s16+s234-s56)*SMB1(18)+4*mup3mup6*(p6sq-s16+s234-s56)*SMB1(
     &   21)-2*s45*SMB1(75)+4*p5mup6*SMB1(91)+2*(-p6sq+s16-s234+s56)*
     &   SMB1(95)
       SMB(221) = 2*(-p6sq+s56)*SMB1(17)
       SMB(222) = 4*p6sq*SMB1(17)
       SMB(223) = 2*(p6sq-s56)*(-2*mup3mup6*SMB1(17)+SMB1(23))
       SMB(224) = 8*mup3mup6*p6sq*SMB1(17)-4*p6sq*SMB1(23)
       SMB(225) = 4*p2mup6*SMB1(17)
       SMB(226) = 2*s12*SMB1(17)
       SMB(227) = 2*(s16-s234+s34-s345)*SMB1(17)
       SMB(228) = -2*(s12+s16-s345)*SMB1(17)
       SMB(229) = (-(p6sq*(s12+s234-s34))+(s12+s16-s345)*s56)*SMB1(17)
     &   
       SMB(230) = 4*mup3mup6*s12*SMB1(17)-2*s12*SMB1(23)
       SMB(231) = 2*(s16-s234+s34-s345)*(2*mup3mup6*SMB1(17)-SMB1(23))
     &   
       SMB(232) = 2*(s12+s16-s345)*(-2*mup3mup6*SMB1(17)+SMB1(23))
       SMB(233) = -2*mup3mup6*(p6sq*(s12+s234-s34)-(s12+s16-s345)*s56)
     &   *SMB1(17)+2*p2mup6*(p6sq-s56)*SMB1(18)+2*p2mup3*(-p6sq+s56)*
     &   SMB1(20)+(p6sq*(s12+s234-s34)-(s12+s16-s345)*s56)*SMB1(23)
       SMB(234) = 2*s12*SMB1(3)-4*p2mup3*SMB1(27)
       SMB(235) = 2*(s16-s234+s34-s345)*SMB1(3)-4*p2mup3*SMB1(43)
       SMB(236) = -2*(s12+s16-s345)*SMB1(3)-4*p2mup3*SMB1(21)
       SMB(237) = 4*p4mup6*SMB1(17)
       SMB(238) = 2*(-s123+s23-s234+s56)*SMB1(17)
       SMB(239) = 2*s45*SMB1(17)
       SMB(240) = 2*(s123-s45-s56)*SMB1(17)
       SMB(241) = -2*s56*SMB1(21)+2*p6sq*(SMB1(21)+SMB1(43))
       SMB(242) = 2*(s123-s23+s234-s56)*(-2*mup3mup6*SMB1(17)+SMB1(23)
     &   )
       SMB(243) = 4*mup3mup6*s45*SMB1(17)-2*s45*SMB1(23)
       SMB(244) = 2*(s123-s45-s56)*(2*mup3mup6*SMB1(17)-SMB1(23))
       SMB(245) = 2*mup3mup6*(p6sq*(s123-s56)+s56*(-s123+s45+s56))*SMB
     &   1(17)+(p6sq-s56)*(-s123+s45+s56)*SMB1(19)+(p6sq*(-s123+s56)-
     &   s56*(-s123+s45+s56))*SMB1(23)+2*(p6sq-s56)*SMB1(95)
       SMB(246) = SMB1(39)
       SMB(247) = SMB1(55)
       SMB(248) = SMB1(69)
       SMB(249) = 4*p1mup3*SMB1(17)
       SMB(250) = 4*p5mup3*SMB1(17)
       SMB(251) = 2*(-s16+s234-s34+s345)*SMB1(24)+4*p5mup3*SMB1(26)
       SMB(252) = 2*p6mup3*(p6sq-s16+s234-s56)*SMB1(17)+(-p6sq+s16-s23
     &   4+s56)*SMB1(18)+2*(p6sq-s56)*SMB1(24)+4*p5mup3*SMB1(28)
       SMB(253) = -2*mup3mup6*(s16**2-s16*(s234-s34+s345)+p6sq*(s12-s1
     &   6+s234-s34+s345)-s12*s56)*SMB1(17)+(s16**2-s16*(s234-s34+s34
     &   5)+p6sq*(s12-s16+s234-s34+s345)-s12*s56)*SMB1(19)-(p6sq-s16)
     &   *(s16-s234+s34-s345)*SMB1(23)-2*(s16-s234+s34-s345)*SMB1(34)
     &   +2*p5mup6*s12*SMB1(111)-2*p5mup3*s12*SMB1(114)
       SMB(254) = 2*(SMB1(18)+SMB1(111))
       SMB(255) = 2*((s12+s16-s345)*SMB1(24)+SMB1(30))
       SMB(256) = 2*p6mup3*(-p6sq+s16)*SMB1(17)+(p6sq-s16)*SMB1(18)-2*
     &   p6sq*SMB1(24)
       SMB(257) = 2*mup3mup6*(s16*(s12+s16-s345)+p6sq*(-s16+s345))*SMB
     &   1(17)+(p6sq*(s16-s345)-s16*(s12+s16-s345))*SMB1(19)+(p6sq-s1
     &   6)*(s12+s16-s345)*SMB1(23)+2*(s12+s16-s345)*SMB1(34)
       SMB(258) = 2*p5mup3*(s12+s16-s345)*SMB1(17)+(s16-s234+s34-s345)
     &   *(2*p6mup3*SMB1(17)-SMB1(18))
       SMB(259) = -2*(p5mup3*p6sq+p6mup3*(p6sq-s56))*SMB1(17)+(p6sq-s5
     &   6)*SMB1(18)
       SMB(260) = -2*mup3mup6*(p6sq*(s12+s234-s34)-(s12+s16-s345)*s56)
     &   *SMB1(17)+p6sq*(s12+s234-s34)*SMB1(19)-(s12+s16-s345)*(s56*S
     &   MB1(19)-2*p5mup6*SMB1(111)+2*p5mup3*SMB1(114))
       SMB(261) = 0
       SMB(262) = 2*s12*SMB1(7)-4*p2mup3*SMB1(37)
       SMB(263) = 2*(s16-s234+s34-s345)*SMB1(7)-4*p2mup3*SMB1(53)
       SMB(264) = -2*(s12+s16-s345)*SMB1(7)-4*p2mup3*SMB1(67)
       SMB(265) = 16*SMB1(17)
       SMB(266) = 4*p6mup6*s12*SMB1(17)-2*s12*SMB1(20)
       SMB(267) = 2*(s16-s234+s34-s345)*(2*p6mup6*SMB1(17)-SMB1(20))
       SMB(268) = 2*(s12+s16-s345)*(-2*p6mup6*SMB1(17)+SMB1(20))
       SMB(269) = 32*(p1mup3*SMB1(17)-SMB1(24))
       SMB(270) = 0
       SMB(271) = 32*p6mup3*SMB1(17)-16*SMB1(111)
       SMB(272) = 0
       SMB(273) = 4*p6mup3*s12*SMB1(17)-2*s12*SMB1(111)
       SMB(274) = 2*(s16-s234+s34-s345)*(2*p6mup3*SMB1(17)-SMB1(111))
       SMB(275) = 0
       SMB(276) = 2*(p6sq-s16)*(-2*p6mup3*SMB1(17)+SMB1(111))
       SMB(277) = 2*(p6sq-s56)*(-2*p6mup3*SMB1(17)+SMB1(111))
       SMB(278) = 2*s12*(-4*p5mup6*p6mup3*SMB1(17)+(-p6sq+s56)*SMB1(23
     &   )+2*p5mup6*SMB1(111))
       SMB(279) = 2*p6sq*s12*(2*mup3mup6*SMB1(17)-SMB1(19))
       SMB(280) = 2*p6sq*(s16-s234+s34-s345)*(2*mup3mup6*SMB1(17)-SMB1
     &   (19))
       SMB(281) = 2*(-p6sq+s16-s234+s56)*SMB1(3)+2*p4mup3*(p6sq-s16+s2
     &   34-s56)*SMB1(17)-2*s45*SMB1(24)+4*p5mup3*SMB1(27)
       SMB(282) = 2*((p6sq-s16)*SMB1(3)+p4mup3*(-p6sq+s16)*SMB1(17)+(-
     &   s123+s45+s56)*SMB1(24)+SMB1(31))
       SMB(283) = 2*p6mup3*s45*SMB1(17)-s45*SMB1(18)-4*p5mup3*SMB1(21)
     &   
       SMB(284) = 0
       SMB(285) = 2*(s123-s23+s234-s56)*(-2*p6mup3*SMB1(17)+SMB1(111))
     &   
       SMB(286) = 4*p6mup3*s45*SMB1(17)-2*s45*SMB1(111)
       SMB(287) = 2*(p6sq-s16)*(-2*mup3mup6*SMB1(17)+SMB1(19))
       SMB(288) = 2*(p6sq-s56)*(-2*mup3mup6*SMB1(17)+SMB1(19))
       SMB(289) = -((p6sq-s16+s234-s56)*(-s123+s45+s56)*SMB1(23))+4*p5
     &   mup6*(s123-s45-s56)*SMB1(24)+4*p5mup3*(-s123+s45+s56)*SMB1(2
     &   5)-2*(p6sq-s56)*SMB1(33)+2*(p6sq-s16)*SMB1(49)-2*SMB1(130)
       SMB(290) = 8*mup3mup6*p6sq*SMB1(17)-4*p6sq*SMB1(19)
       SMB(291) = 2*(p6sq-s16)*SMB1(22)-(p6sq-s16)*(s123-s45-s56)*SMB1
     &   (23)+2*p6sq*SMB1(33)+2*(-s123+s45+s56)*SMB1(34)
       SMB(292) = 2*((p6sq-s56)*SMB1(22)+p6sq*SMB1(49)+p5mup6*(-s123+s
     &   45+s56)*SMB1(111)+p5mup3*(s123-s45-s56)*SMB1(114))
       SMB(293) = 32*p6mup6*SMB1(17)-16*SMB1(20)
       SMB(294) = 2*(s123-s23+s234-s56)*(-2*p6mup6*SMB1(17)+SMB1(20))
       SMB(295) = 4*p6mup6*s45*SMB1(17)-2*s45*SMB1(20)
       SMB(296) = 2*(s123-s45-s56)*(2*p6mup6*SMB1(17)-SMB1(20))
       SMB(297) = 16*(2*(-2*p1mup6*p6mup3+2*p1mup3*p6mup6+mup3mup6*(p6
     &   sq-s16))*SMB1(17)+(-p6sq+s16)*SMB1(19)+(-p6sq+s16)*SMB1(23)-
     &   4*p6mup6*SMB1(24)+4*mup3mup6*SMB1(28)-2*SMB1(78)+2*p1mup6*SM
     &   B1(111))
       SMB(298) = 16*(-p6sq+s56)*SMB1(23)+32*p5mup6*(-2*p6mup3*SMB1(17
     &   )+SMB1(111))
       SMB(299) = 2*(s123-s23+s234-s56)*(4*p5mup6*p6mup3*SMB1(17)+(p6s
     &   q-s56)*SMB1(23)-2*p5mup6*SMB1(111))
       SMB(300) = 16*p6sq*(2*mup3mup6*SMB1(17)-SMB1(19))
       SMB(301) = 2*p6sq*(s123-s23+s234-s56)*(-2*mup3mup6*SMB1(17)+SMB
     &   1(19))
       SMB(302) = 2*p6sq*s45*(2*mup3mup6*SMB1(17)-SMB1(19))
       SMB(303) = 4*mup3mup6*SMB1(17)
c       Print*," SMB(1) ", SMB(1)
c       Print*," SMB(2) ", SMB(2)
c       Print*," SMB(3) ", SMB(3)
c       Print*," SMB(4) ", SMB(4)
c       Print*," SMB(5) ", SMB(5)
c       Print*," SMB(6) ", SMB(6)
c       Print*," SMB(7) ", SMB(7)
c       Print*," SMB(8) ", SMB(8)
c       Print*," SMB(9) ", SMB(9)
c       Print*," SMB(10) ", SMB(10)
c       Print*," SMB(11) ", SMB(11)
c       Print*," SMB(12) ", SMB(12)
c       Print*," SMB(13) ", SMB(13)
c       Print*," SMB(14) ", SMB(14)
c       Print*," SMB(15) ", SMB(15)
c       Print*," SMB(16) ", SMB(16)
c       Print*," SMB(17) ", SMB(17)
c       Print*," SMB(18) ", SMB(18)
c       Print*," SMB(19) ", SMB(19)
c       Print*," SMB(20) ", SMB(20)
c       Print*," SMB(21) ", SMB(21)
c       Print*," SMB(22) ", SMB(22)
c       Print*," SMB(23) ", SMB(23)
c       Print*," SMB(24) ", SMB(24)
c       Print*," SMB(25) ", SMB(25)
c       Print*," SMB(26) ", SMB(26)
c       Print*," SMB(27) ", SMB(27)
c       Print*," SMB(28) ", SMB(28)
c       Print*," SMB(29) ", SMB(29)
c       Print*," SMB(30) ", SMB(30)
c       Print*," SMB(31) ", SMB(31)
c       Print*," SMB(32) ", SMB(32)
c       Print*," SMB(33) ", SMB(33)
c       Print*," SMB(34) ", SMB(34)
c       Print*," SMB(35) ", SMB(35)
c       Print*," SMB(36) ", SMB(36)
c       Print*," SMB(37) ", SMB(37)
c       Print*," SMB(38) ", SMB(38)
c       Print*," SMB(39) ", SMB(39)
c       Print*," SMB(40) ", SMB(40)
c       Print*," SMB(41) ", SMB(41)
c       Print*," SMB(42) ", SMB(42)
c       Print*," SMB(43) ", SMB(43)
c       Print*," SMB(44) ", SMB(44)
c       Print*," SMB(45) ", SMB(45)
c       Print*," SMB(46) ", SMB(46)
c       Print*," SMB(47) ", SMB(47)
c       Print*," SMB(48) ", SMB(48)
c       Print*," SMB(49) ", SMB(49)
c       Print*," SMB(50) ", SMB(50)
c       Print*," SMB(51) ", SMB(51)
c       Print*," SMB(52) ", SMB(52)
c       Print*," SMB(53) ", SMB(53)
c       Print*," SMB(54) ", SMB(54)
c       Print*," SMB(55) ", SMB(55)
c       Print*," SMB(56) ", SMB(56)
c       Print*," SMB(57) ", SMB(57)
c       Print*," SMB(58) ", SMB(58)
c       Print*," SMB(59) ", SMB(59)
c       Print*," SMB(60) ", SMB(60)
c       Print*," SMB(61) ", SMB(61)
c       Print*," SMB(62) ", SMB(62)
c       Print*," SMB(63) ", SMB(63)
c       Print*," SMB(64) ", SMB(64)
c       Print*," SMB(65) ", SMB(65)
c       Print*," SMB(66) ", SMB(66)
c       Print*," SMB(67) ", SMB(67)
c       Print*," SMB(68) ", SMB(68)
c       Print*," SMB(69) ", SMB(69)
c       Print*," SMB(70) ", SMB(70)
c       Print*," SMB(71) ", SMB(71)
c       Print*," SMB(72) ", SMB(72)
c       Print*," SMB(73) ", SMB(73)
c       Print*," SMB(74) ", SMB(74)
c       Print*," SMB(75) ", SMB(75)
c       Print*," SMB(76) ", SMB(76)
c       Print*," SMB(77) ", SMB(77)
c       Print*," SMB(78) ", SMB(78)
c       Print*," SMB(79) ", SMB(79)
c       Print*," SMB(80) ", SMB(80)
c       Print*," SMB(81) ", SMB(81)
c       Print*," SMB(82) ", SMB(82)
c       Print*," SMB(83) ", SMB(83)
c       Print*," SMB(84) ", SMB(84)
c       Print*," SMB(85) ", SMB(85)
c       Print*," SMB(86) ", SMB(86)
c       Print*," SMB(87) ", SMB(87)
c       Print*," SMB(88) ", SMB(88)
c       Print*," SMB(89) ", SMB(89)
c       Print*," SMB(90) ", SMB(90)
c       Print*," SMB(91) ", SMB(91)
c       Print*," SMB(92) ", SMB(92)
c       Print*," SMB(93) ", SMB(93)
c       Print*," SMB(94) ", SMB(94)
c       Print*," SMB(95) ", SMB(95)
c       Print*," SMB(96) ", SMB(96)
c       Print*," SMB(97) ", SMB(97)
c       Print*," SMB(98) ", SMB(98)
c       Print*," SMB(99) ", SMB(99)
c       Print*," SMB(100) ", SMB(100)
c       Print*," SMB(101) ", SMB(101)
c       Print*," SMB(102) ", SMB(102)
c       Print*," SMB(103) ", SMB(103)
c       Print*," SMB(104) ", SMB(104)
c       Print*," SMB(105) ", SMB(105)
c       Print*," SMB(106) ", SMB(106)
c       Print*," SMB(107) ", SMB(107)
c       Print*," SMB(108) ", SMB(108)
c       Print*," SMB(109) ", SMB(109)
c       Print*," SMB(110) ", SMB(110)
c       Print*," SMB(111) ", SMB(111)
c       Print*," SMB(112) ", SMB(112)
c       Print*," SMB(113) ", SMB(113)
c       Print*," SMB(114) ", SMB(114)
c       Print*," SMB(115) ", SMB(115)
c       Print*," SMB(116) ", SMB(116)
c       Print*," SMB(117) ", SMB(117)
c       Print*," SMB(118) ", SMB(118)
c       Print*," SMB(119) ", SMB(119)
c       Print*," SMB(120) ", SMB(120)
c       Print*," SMB(121) ", SMB(121)
c       Print*," SMB(122) ", SMB(122)
c       Print*," SMB(123) ", SMB(123)
c       Print*," SMB(124) ", SMB(124)
c       Print*," SMB(125) ", SMB(125)
c       Print*," SMB(126) ", SMB(126)
c       Print*," SMB(127) ", SMB(127)
c       Print*," SMB(128) ", SMB(128)
c       Print*," SMB(129) ", SMB(129)
c       Print*," SMB(130) ", SMB(130)
c       Print*," SMB(131) ", SMB(131)
c       Print*," SMB(132) ", SMB(132)
c       Print*," SMB(133) ", SMB(133)
c       Print*," SMB(134) ", SMB(134)
c       Print*," SMB(135) ", SMB(135)
c       Print*," SMB(136) ", SMB(136)
c       Print*," SMB(137) ", SMB(137)
c       Print*," SMB(138) ", SMB(138)
c       Print*," SMB(139) ", SMB(139)
c       Print*," SMB(140) ", SMB(140)
c       Print*," SMB(141) ", SMB(141)
c       Print*," SMB(142) ", SMB(142)
c       Print*," SMB(143) ", SMB(143)
c       Print*," SMB(144) ", SMB(144)
c       Print*," SMB(145) ", SMB(145)
c       Print*," SMB(146) ", SMB(146)
c       Print*," SMB(147) ", SMB(147)
c       Print*," SMB(148) ", SMB(148)
c       Print*," SMB(149) ", SMB(149)
c       Print*," SMB(150) ", SMB(150)
c       Print*," SMB(151) ", SMB(151)
c       Print*," SMB(152) ", SMB(152)
c       Print*," SMB(153) ", SMB(153)
c       Print*," SMB(154) ", SMB(154)
c       Print*," SMB(155) ", SMB(155)
c       Print*," SMB(156) ", SMB(156)
c       Print*," SMB(157) ", SMB(157)
c       Print*," SMB(158) ", SMB(158)
c       Print*," SMB(159) ", SMB(159)
c       Print*," SMB(160) ", SMB(160)
c       Print*," SMB(161) ", SMB(161)
c       Print*," SMB(162) ", SMB(162)
c       Print*," SMB(163) ", SMB(163)
c       Print*," SMB(164) ", SMB(164)
c       Print*," SMB(165) ", SMB(165)
c       Print*," SMB(166) ", SMB(166)
c       Print*," SMB(167) ", SMB(167)
c       Print*," SMB(168) ", SMB(168)
c       Print*," SMB(169) ", SMB(169)
c       Print*," SMB(170) ", SMB(170)
c       Print*," SMB(171) ", SMB(171)
c       Print*," SMB(172) ", SMB(172)
c       Print*," SMB(173) ", SMB(173)
c       Print*," SMB(174) ", SMB(174)
c       Print*," SMB(175) ", SMB(175)
c       Print*," SMB(176) ", SMB(176)
c       Print*," SMB(177) ", SMB(177)
c       Print*," SMB(178) ", SMB(178)
c       Print*," SMB(179) ", SMB(179)
c       Print*," SMB(180) ", SMB(180)
c       Print*," SMB(181) ", SMB(181)
c       Print*," SMB(182) ", SMB(182)
c       Print*," SMB(183) ", SMB(183)
c       Print*," SMB(184) ", SMB(184)
c       Print*," SMB(185) ", SMB(185)
c       Print*," SMB(186) ", SMB(186)
c       Print*," SMB(187) ", SMB(187)
c       Print*," SMB(188) ", SMB(188)
c       Print*," SMB(189) ", SMB(189)
c       Print*," SMB(190) ", SMB(190)
c       Print*," SMB(191) ", SMB(191)
c       Print*," SMB(192) ", SMB(192)
c       Print*," SMB(193) ", SMB(193)
c       Print*," SMB(194) ", SMB(194)
c       Print*," SMB(195) ", SMB(195)
c       Print*," SMB(196) ", SMB(196)
c       Print*," SMB(197) ", SMB(197)
c       Print*," SMB(198) ", SMB(198)
c       Print*," SMB(199) ", SMB(199)
c       Print*," SMB(200) ", SMB(200)
c       Print*," SMB(201) ", SMB(201)
c       Print*," SMB(202) ", SMB(202)
c       Print*," SMB(203) ", SMB(203)
c       Print*," SMB(204) ", SMB(204)
c       Print*," SMB(205) ", SMB(205)
c       Print*," SMB(206) ", SMB(206)
c       Print*," SMB(207) ", SMB(207)
c       Print*," SMB(208) ", SMB(208)
c       Print*," SMB(209) ", SMB(209)
c       Print*," SMB(210) ", SMB(210)
c       Print*," SMB(211) ", SMB(211)
c       Print*," SMB(212) ", SMB(212)
c       Print*," SMB(213) ", SMB(213)
c       Print*," SMB(214) ", SMB(214)
c       Print*," SMB(215) ", SMB(215)
c       Print*," SMB(216) ", SMB(216)
c       Print*," SMB(217) ", SMB(217)
c       Print*," SMB(218) ", SMB(218)
c       Print*," SMB(219) ", SMB(219)
c       Print*," SMB(220) ", SMB(220)
c       Print*," SMB(221) ", SMB(221)
c       Print*," SMB(222) ", SMB(222)
c       Print*," SMB(223) ", SMB(223)
c       Print*," SMB(224) ", SMB(224)
c       Print*," SMB(225) ", SMB(225)
c       Print*," SMB(226) ", SMB(226)
c       Print*," SMB(227) ", SMB(227)
c       Print*," SMB(228) ", SMB(228)
c       Print*," SMB(229) ", SMB(229)
c       Print*," SMB(230) ", SMB(230)
c       Print*," SMB(231) ", SMB(231)
c       Print*," SMB(232) ", SMB(232)
c       Print*," SMB(233) ", SMB(233)
c       Print*," SMB(234) ", SMB(234)
c       Print*," SMB(235) ", SMB(235)
c       Print*," SMB(236) ", SMB(236)
c       Print*," SMB(237) ", SMB(237)
c       Print*," SMB(238) ", SMB(238)
c       Print*," SMB(239) ", SMB(239)
c       Print*," SMB(240) ", SMB(240)
c       Print*," SMB(241) ", SMB(241)
c       Print*," SMB(242) ", SMB(242)
c       Print*," SMB(243) ", SMB(243)
c       Print*," SMB(244) ", SMB(244)
c       Print*," SMB(245) ", SMB(245)
c       Print*," SMB(246) ", SMB(246)
c       Print*," SMB(247) ", SMB(247)
c       Print*," SMB(248) ", SMB(248)
c       Print*," SMB(249) ", SMB(249)
c       Print*," SMB(250) ", SMB(250)
c       Print*," SMB(251) ", SMB(251)
c       Print*," SMB(252) ", SMB(252)
c       Print*," SMB(253) ", SMB(253)
c       Print*," SMB(254) ", SMB(254)
c       Print*," SMB(255) ", SMB(255)
c       Print*," SMB(256) ", SMB(256)
c       Print*," SMB(257) ", SMB(257)
c       Print*," SMB(258) ", SMB(258)
c       Print*," SMB(259) ", SMB(259)
c       Print*," SMB(260) ", SMB(260)
c       Print*," SMB(261) ", SMB(261)
c       Print*," SMB(262) ", SMB(262)
c       Print*," SMB(263) ", SMB(263)
c       Print*," SMB(264) ", SMB(264)
c       Print*," SMB(265) ", SMB(265)
c       Print*," SMB(266) ", SMB(266)
c       Print*," SMB(267) ", SMB(267)
c       Print*," SMB(268) ", SMB(268)
c       Print*," SMB(269) ", SMB(269)
c       Print*," SMB(270) ", SMB(270)
c       Print*," SMB(271) ", SMB(271)
c       Print*," SMB(272) ", SMB(272)
c       Print*," SMB(273) ", SMB(273)
c       Print*," SMB(274) ", SMB(274)
c       Print*," SMB(275) ", SMB(275)
c       Print*," SMB(276) ", SMB(276)
c       Print*," SMB(277) ", SMB(277)
c       Print*," SMB(278) ", SMB(278)
c       Print*," SMB(279) ", SMB(279)
c       Print*," SMB(280) ", SMB(280)
c       Print*," SMB(281) ", SMB(281)
c       Print*," SMB(282) ", SMB(282)
c       Print*," SMB(283) ", SMB(283)
c       Print*," SMB(284) ", SMB(284)
c       Print*," SMB(285) ", SMB(285)
c       Print*," SMB(286) ", SMB(286)
c       Print*," SMB(287) ", SMB(287)
c       Print*," SMB(288) ", SMB(288)
c       Print*," SMB(289) ", SMB(289)
c       Print*," SMB(290) ", SMB(290)
c       Print*," SMB(291) ", SMB(291)
c       Print*," SMB(292) ", SMB(292)
c       Print*," SMB(293) ", SMB(293)
c       Print*," SMB(294) ", SMB(294)
c       Print*," SMB(295) ", SMB(295)
c       Print*," SMB(296) ", SMB(296)
c       Print*," SMB(297) ", SMB(297)
c       Print*," SMB(298) ", SMB(298)
c       Print*," SMB(299) ", SMB(299)
c       Print*," SMB(300) ", SMB(300)
c       Print*," SMB(301) ", SMB(301)
c       Print*," SMB(302) ", SMB(302)
c       Print*," SMB(303) ", SMB(303)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Amplitude                         
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************


c The Finite virtual piece should !!NOT!!!! be multiplied for (-1)
c. The Born has a factor 1 and the Virtual (-I). The factorization
c from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, (-I)*I=(1)!!!
       result(1) = (8*Fa(1)*SMB(1)+Fa(2)*SMB(2)+F(11)*SMB(18)+F(12)*SM
     &   B(19)+F(18)*SMB(33)+Fa(3)*SMB(48)+Fa(4)*SMB(49)+Fa(5)*SMB(50
     &   )+F(30)*SMB(53)+F(32)*SMB(65)+F(33)*SMB(67)-mup3mup6*(-(F(9)
     &   *SMB(14))-4*(F(20)*SMB(43)-F(21)*SMB(44))+16*F(35)*SMB(76)+8
     &   *(F(6)*SMB(3)+F(7)*SMB(4)+s16*F(2)*SMB(13)+F(10)*SMB(15)+F(3
     &   )*SMB(16)-F(15)*SMB(27)+F(22)*SMB(109)))-F(36)*SMB(116)+F(37
     &   )*SMB(119)+F(22)*(4*SMB(110)+2*SMB(123)+Fa(6)*(4*SMB(120)-2*
     &   SMB(125))+SMB(128))+4*(F(13)*SMB(20)+F(14)*SMB(25)+F(17)*SMB
     &   (32)-p6mup3*F(22)*SMB(47)+F(8)*(SMB(6)+SMB(64))+F(34)*SMB(70
     &   )+F(40)*SMB(177))-2*(F(16)*SMB(28)-F(19)*SMB(34)-F(6)*(2*SMB
     &   (5)+SMB(51))-s16*F(2)*(2*SMB(17)+SMB(52))-F(31)*(SMB(55)+SMB
     &   (69))+F(21)*SMB(111)-F(20)*SMB(118)-Fa(7)*(-2*SMB(121)+SMB(1
     &   26))-8*p4mup6*F(39)*SMB(141)-4*(p4mup6*F(38)*SMB(143)-F(39)*
     &   SMB(181)+Fa(6)*(F(38)*SMB(138)+2*F(39)*SMB(140)-F(35)*SMB(18
     &   5)))-2*p5mup6*F(35)*(-4*SMB(139)+2*SMB(183)+SMB(187)))+F(35)
     &   *(8*(SMB(216)+mup3mup6*(-SMB(197)+SMB(241)))+4*(SMB(203)+SMB
     &   (220)-SMB(245))-2*(p4mup6*(2*SMB(205)-SMB(275)+SMB(277))-Fa(
     &   6)*(-SMB(284)+SMB(286)))+SMB(299))+F(1842)*(-12*mup3mup6*(SM
     &   B(195)+SMB(221)+SMB(222))+6*(SMB(117)+SMB(200)+SMB(223)+SMB(
     &   224))-2*Fa(6)*(SMB(269)+SMB(270)+SMB(271))+SMB(297)+SMB(298)
     &   +SMB(300))+F(41)*SMB(302))/(s16*s23)
       result(2) = (Fa(8)*SMB(2)+32*Fa(9)*SMB(4)+Fa(12)*SMB(9)+F(55)*S
     &   MB(19)-16*(Fa(14)*SMB(14)-Fa(15)*SMB(16)+Fa(17)*SMB(22))+Fa(
     &   21)*SMB(31)+F(69)*SMB(33)+Fa(23)*SMB(37)+Fa(28)*SMB(49)+Fa(2
     &   9)*SMB(50)+Fa(30)*SMB(58)+Fa(31)*SMB(59)+F(83)*SMB(61)+F(84)
     &   *SMB(65)+F(85)*SMB(67)+8*(Fa(10)*SMB(7)+Fa(11)*SMB(8)+F(59)*
     &   SMB(26)-Fa(18)*SMB(27)+Fa(19)*SMB(29)-Fa(20)*SMB(30)-F(66)*S
     &   MB(32)+Fa(22)*SMB(35)+Fa(24)*SMB(39)-F(74)*SMB(42)+Fa(32)*SM
     &   B(75))+F(87)*SMB(103)+F(63)*(8*SMB(41)+4*SMB(45)-SMB(116)-2*
     &   (SMB(111)-p5mup3*(2*(SMB(88)+SMB(90))+SMB(113))+Fa(33)*(2*SM
     &   B(89)+SMB(114))-SMB(118)))+F(88)*SMB(119)+F(64)*(4*SMB(107)+
     &   SMB(127))+2*(F(65)*SMB(28)+F(70)*SMB(34)-F(50)*(2*(SMB(11)-S
     &   MB(56))+SMB(57))-s16*F(49)*(2*(SMB(23)-SMB(60))+SMB(66))-F(8
     &   6)*SMB(68)-F(81)*(SMB(55)+SMB(69))-8*p2mup3*F(46)*SMB(74)+4*
     &   F(46)*(SMB(82)-2*(p5mup6*SMB(85)-Fa(6)*SMB(86)+p5mup6*SMB(87
     &   )))+Fa(35)*SMB(115)+F(52)*(4*SMB(83)+2*Fa(6)*SMB(84)+SMB(122
     &   ))+Fa(36)*SMB(126)+Fa(6)*(F(52)*SMB(124)+8*F(53)*SMB(135)))+
     &   4*(F(58)*SMB(25)-Fa(25)*SMB(43)+Fa(26)*SMB(44)-F(75)*SMB(46)
     &   +Fa(27)*SMB(48)+F(57)*(-SMB(24)+SMB(62))+F(42)*(-SMB(20)+SMB
     &   (63))-F(47)*(SMB(6)-SMB(54)+SMB(64))+p2mup6*(4*F(46)*SMB(73)
     &   +F(64)*SMB(121))+F(54)*(-2*(Fa(6)*SMB(137)-p2mup6*SMB(143))+
     &   SMB(148))-2*(F(45)*(2*p2mup6*SMB(141)+SMB(144))-Fa(33)*(-(F(
     &   54)*SMB(132))+F(46)*(2*SMB(93)-SMB(191)+SMB(212))))+F(46)*(4
     &   *p5mup6*SMB(136)-p5mup6*SMB(186)-SMB(201)+p2mup6*SMB(205)+2*
     &   (SMB(99)-SMB(102)+Fa(34)*SMB(105)-p5mup6*SMB(182)+Fa(6)*SMB(
     &   184)-SMB(206))+p2mup3*(4*(SMB(94)-SMB(96))+SMB(209)-SMB(213)
     &   )-SMB(214)-Fa(38)*SMB(221))+p5mup3*(4*F(45)*SMB(130)+2*F(54)
     &   *SMB(134)+F(46)*(-4*SMB(92)+2*(SMB(189)+SMB(193)-SMB(210))+p
     &   2mup6*(SMB(195)+SMB(222))))+Fa(16)*(4*F(56)*SMB(21)+F(46)*SM
     &   B(227))-Fa(13)*(8*F(53)*SMB(10)+F(46)*(SMB(226)+SMB(228)))+F
     &   (46)*SMB(233))+F(89)*SMB(278)+F(46)*(4*(SMB(253)-SMB(260))-2
     &   *(Fa(33)*SMB(267)-p5mup3*(SMB(266)+SMB(268))-Fa(6)*(4*(SMB(2
     &   51)-SMB(258))+SMB(272)-SMB(274))-p2mup6*(-4*(SMB(252)-SMB(25
     &   9))-SMB(275)+SMB(277)))+SMB(280))-F(1841)*(48*SMB(72)+12*SMB
     &   (78)+24*(SMB(168)+SMB(170)+SMB(173))-4*Fa(39)*SMB(265)-6*(SM
     &   B(117)+SMB(200)+SMB(223)+SMB(224)-SMB(287)-SMB(288)-SMB(290)
     &   )+2*(Fa(6)*(6*(SMB(249)+SMB(250)+SMB(254))+SMB(269)+SMB(270)
     &   +SMB(271))+Fa(37)*(6*SMB(188)+SMB(293)))-SMB(297)-SMB(298)-S
     &   MB(300)))/(s16*s34)
       result(3) = (Fa(40)*SMB(48)-s56*(-(F(97)*SMB(51))+2*(F(91)*(2*S
     &   MB(25)+SMB(28))+F(98)*(SMB(52)+SMB(53)))-F(92)*(2*SMB(34)+SM
     &   B(65))+F(93)*(4*mup3mup6*(SMB(43)+SMB(75))-2*(SMB(111)+SMB(1
     &   12))+SMB(116)+SMB(117)))+F(100)*(2*SMB(123)+SMB(128))+4*(p6m
     &   up3*F(94)*SMB(47)-F(90)*SMB(70)-F(101)*(SMB(177)+SMB(178)))-
     &   2*(F(99)*SMB(55)-F(96)*(SMB(67)+SMB(69))+p5mup6*(-2*(F(96)*(
     &   SMB(49)+SMB(50))+F(94)*SMB(120))-F(100)*SMB(125)+4*(F(90)*SM
     &   B(1)+F(101)*(SMB(138)+SMB(139)))+F(102)*(SMB(284)+SMB(285)))
     &   )+F(1840)*(2*p5mup6*(SMB(269)+SMB(270)+SMB(271))+SMB(297)+SM
     &   B(298)+SMB(300))+F(103)*(SMB(299)+SMB(301)))/(s23*s56)
       result(4) = (s56*(F(107)*SMB(25)+2*F(111)*SMB(28))+F(121)*SMB(5
     &   4)+Fa(47)*SMB(58)+Fa(51)*(SMB(114)+SMB(115))+s56*(F(123)*SMB
     &   (56)+F(124)*SMB(57)-4*F(1837)*(SMB(32)+SMB(64))+F(127)*SMB(6
     &   5)+2*F(115)*(-2*(SMB(60)+SMB(62))+SMB(66)+SMB(68))+F(119)*(4
     &   *SMB(45)+2*(SMB(111)+SMB(112))-SMB(116)-SMB(117)))+F(131)*SM
     &   B(127)+8*(Fa(42)*SMB(27)+Fa(43)*SMB(29)+Fa(44)*(SMB(30)+SMB(
     &   31))+p1mup3*F(106)*(SMB(132)+SMB(134))+F(105)*(-2*Fa(52)*SMB
     &   (130)+SMB(144)+SMB(146)))+4*(Fa(46)*SMB(48)+Fa(48)*SMB(59)+F
     &   (114)*(SMB(61)+SMB(63))-Fa(45)*(SMB(43)+SMB(75))+F(129)*SMB(
     &   78)+Fa(50)*(SMB(89)+SMB(90))+F(104)*SMB(103)+F(118)*(2*p6mup
     &   3*(SMB(208)+SMB(210))-p2mup3*(SMB(209)+SMB(211)))-Fa(41)*(8*
     &   F(105)*SMB(10)+F(118)*SMB(226)))+2*(F(116)*SMB(34)-F(122)*SM
     &   B(55)-F(128)*(SMB(67)+SMB(69))+Fa(49)*(-2*SMB(88)+SMB(113))+
     &   F(110)*SMB(122)+p5mup6*(4*(F(104)*SMB(7)-F(117)*SMB(35)-F(11
     &   4)*(SMB(37)+SMB(39)))-F(110)*(2*SMB(84)+SMB(124))+8*F(105)*(
     &   SMB(135)+SMB(136))-2*(F(109)*(SMB(49)+SMB(50))-p1mup3*(4*F(1
     &   06)*(SMB(21)+SMB(22))+F(118)*(SMB(227)+SMB(228))))))+F(118)*
     &   (-4*(SMB(253)+SMB(257))-2*(Fa(52)*SMB(266)-p1mup3*(SMB(267)+
     &   SMB(268))-p5mup6*(4*(SMB(251)+SMB(255))+SMB(272)+SMB(273)))+
     &   SMB(278)+SMB(279))+F(1839)*(12*SMB(46)-4*Fa(53)*SMB(265)-6*(
     &   SMB(287)+SMB(288)+SMB(290))+2*(p5mup6*(6*(SMB(249)+SMB(250)+
     &   SMB(254))+SMB(269)+SMB(270)+SMB(271))-Fa(37)*SMB(293))+SMB(2
     &   97)+SMB(298)+SMB(300)))/(s34*s56)
       result(5) = (Fa(54)*SMB(1)+Fa(55)*SMB(2)+F(147)*SMB(19)+F(148)*
     &   SMB(20)+F(152)*SMB(32)+F(153)*SMB(33)+Fa(56)*SMB(48)+Fa(57)*
     &   SMB(49)+Fa(58)*SMB(50)+F(154)*(2*SMB(34)+SMB(65))+F(175)*SMB
     &   (67)+F(176)*SMB(69)+F(177)*SMB(70)+F(182)*SMB(116)+F(183)*SM
     &   B(117)+F(184)*SMB(119)+F(193)*(2*SMB(123)+SMB(128))+8*(Fa(63
     &   )*SMB(138)+Fa(64)*SMB(139)+Fa(65)*SMB(140)-Fa(66)*SMB(141)-F
     &   a(67)*SMB(142)+Fa(68)*SMB(143))-4*(F(140)*SMB(5)+F(142)*SMB(
     &   6)+F(143)*SMB(17)+F(145)*SMB(18)-F(149)*SMB(25)-F(151)*SMB(2
     &   8)+p6mup3*F(157)*SMB(47)-F(174)*SMB(64)-F(155)*SMB(111)+F(17
     &   8)*SMB(112)+F(156)*SMB(118)-Fa(59)*SMB(120)-Fa(60)*SMB(121)-
     &   F(222)*SMB(177)-F(223)*SMB(178)+F(224)*SMB(181)+Fa(69)*SMB(1
     &   83)-Fa(70)*SMB(185)-Fa(71)*SMB(187)-Fa(72)*SMB(205)-F(179)*S
     &   MB(207))+mup3mup6*(F(144)*SMB(14)+8*(F(140)*SMB(3)+F(141)*SM
     &   B(4)+F(143)*SMB(13)+F(145)*SMB(15))+F(146)*SMB(16)-4*(F(150)
     &   *SMB(27)+2*(F(155)*SMB(43)-F(156)*SMB(44)-F(178)*SMB(75)+F(1
     &   79)*SMB(77))-F(180)*SMB(109)-F(238)*SMB(197)-F(239)*SMB(199)
     &   +F(237)*(SMB(195)+SMB(221)+SMB(222))-F(244)*(SMB(238)+SMB(23
     &   9)+SMB(240))+F(212)*SMB(241)))-2*(F(170)*SMB(51)+F(171)*SMB(
     &   52)+F(172)*SMB(53)-F(173)*SMB(55)-F(181)*SMB(110)+Fa(61)*SMB
     &   (125)-Fa(62)*SMB(126)+F(239)*SMB(204)+F(217)*SMB(215)+F(238)
     &   *(SMB(203)+SMB(220))-F(237)*(SMB(200)+SMB(223)+SMB(224))+F(2
     &   44)*(SMB(242)+SMB(243)+SMB(244))+Fa(73)*(SMB(269)+SMB(270)+S
     &   MB(271))-Fa(74)*SMB(275)-Fa(75)*SMB(276)-Fa(76)*SMB(277)-Fa(
     &   77)*SMB(284)-Fa(78)*SMB(285)-Fa(79)*SMB(286))+F(271)*SMB(299
     &   )+F(1846)*(SMB(297)+SMB(298)+SMB(300))+F(272)*SMB(301)+F(212
     &   )*(2*SMB(245)+SMB(302)))/s23
       result(6) = (Fa(80)*SMB(2)+F(288)*SMB(6)+Fa(82)*SMB(7)+Fa(83)*S
     &   MB(8)+Fa(84)*SMB(9)+F(328)*SMB(19)+F(329)*SMB(20)-16*(Fa(81)
     &   *SMB(4)-Fa(85)*SMB(10)-Fa(86)*SMB(14)+Fa(87)*SMB(16)-Fa(88)*
     &   SMB(21)-Fa(89)*SMB(22))+F(347)*SMB(25)+F(348)*SMB(26)+Fa(90)
     &   *SMB(27)+F(365)*SMB(28)+Fa(91)*SMB(29)+Fa(92)*SMB(30)+Fa(93)
     &   *SMB(31)+F(376)*SMB(32)+F(377)*SMB(33)+Fa(94)*SMB(35)+Fa(95)
     &   *SMB(37)+Fa(96)*SMB(39)+Fa(97)*SMB(43)+Fa(99)*SMB(48)+Fa(100
     &   )*SMB(49)+Fa(101)*SMB(50)+F(427)*SMB(54)+F(428)*SMB(55)+Fa(1
     &   03)*SMB(59)+F(438)*SMB(61)+F(439)*SMB(63)+F(440)*SMB(64)+F(4
     &   41)*SMB(65)+F(442)*SMB(67)+F(444)*SMB(69)+F(486)*SMB(116)+F(
     &   487)*SMB(117)+F(488)*SMB(119)+F(497)*SMB(127)-8*(Fa(104)*SMB
     &   (73)-Fa(105)*SMB(74)+Fa(110)*SMB(87)+Fa(111)*SMB(88)-Fa(112)
     &   *SMB(89)+Fa(113)*SMB(90)-p6mup3*(F(300)*SMB(92)-F(471)*SMB(9
     &   3))-Fa(114)*SMB(94)+Fa(115)*SMB(95)-Fa(116)*SMB(96)-Fa(117)*
     &   SMB(105)+Fa(123)*SMB(130)-Fa(124)*SMB(132)-Fa(125)*SMB(134)-
     &   Fa(126)*SMB(135)+Fa(127)*SMB(136)-Fa(128)*SMB(137)-Fa(129)*S
     &   MB(141)-Fa(130)*SMB(142)+Fa(131)*SMB(143)-F(498)*SMB(144))+F
     &   a(139)*SMB(210)+4*(F(307)*SMB(11)+F(294)*SMB(23)+F(346)*SMB(
     &   24)+F(391)*SMB(41)+F(392)*SMB(42)-Fa(98)*SMB(44)+F(413)*SMB(
     &   45)+F(414)*SMB(46)-F(429)*SMB(56)-Fa(102)*SMB(58)-F(437)*SMB
     &   (60)-F(375)*SMB(62)+F(445)*SMB(72)-Fa(106)*SMB(75)+F(460)*SM
     &   B(78)+F(299)*SMB(83)+F(473)*SMB(103)+F(465)*SMB(111)-F(468)*
     &   SMB(112)-F(469)*SMB(118)+Fa(121)*SMB(121)+F(499)*SMB(146)+F(
     &   501)*SMB(148)+Fa(108)*(2*SMB(85)+SMB(182))-Fa(109)*(2*SMB(86
     &   )+SMB(184))-Fa(132)*SMB(186)-Fa(136)*SMB(195)+Fa(137)*SMB(20
     &   5)-F(472)*SMB(207)+Fa(138)*SMB(209)-Fa(140)*SMB(211)-Fa(141)
     &   *SMB(212))+2*(F(378)*SMB(34)+F(367)*SMB(57)+F(370)*SMB(66)+F
     &   (443)*SMB(68)+F(483)*SMB(107)+Fa(118)*SMB(113)+Fa(119)*SMB(1
     &   14)-Fa(120)*SMB(115)+F(436)*SMB(122)-Fa(107)*(2*SMB(84)+SMB(
     &   124))+Fa(122)*SMB(126)+F(451)*(2*SMB(80)-SMB(215))+F(507)*(-
     &   4*(SMB(168)+SMB(170)+SMB(173))-2*Fa(37)*SMB(188)+SMB(200)+SM
     &   B(223)+SMB(224))+F(466)*(-4*(SMB(97)+SMB(98)+SMB(100))-2*Fa(
     &   37)*SMB(225)+SMB(230)+SMB(231)+SMB(232))+F(471)*(-SMB(233)+2
     &   *(SMB(102)+SMB(260)))+2*(Fa(142)*SMB(213)+Fa(143)*SMB(221)+F
     &   a(144)*SMB(222)+Fa(145)*SMB(226)-Fa(146)*SMB(227)+Fa(147)*SM
     &   B(228)-Fa(155)*SMB(265))-Fa(133)*(2*SMB(189)+SMB(266))+Fa(13
     &   4)*(2*SMB(191)+SMB(267))+Fa(135)*(2*SMB(193)+SMB(268))+Fa(14
     &   8)*(6*(SMB(249)+SMB(250)+SMB(254))+SMB(269)+SMB(270)+SMB(271
     &   ))+Fa(149)*(4*SMB(251)+SMB(272))+Fa(151)*(4*SMB(255)+SMB(273
     &   ))+Fa(153)*(4*SMB(258)+SMB(274))+Fa(150)*(4*SMB(252)+SMB(275
     &   ))+Fa(152)*(4*SMB(256)+SMB(276))+Fa(154)*(4*SMB(259)+SMB(277
     &   )))+F(450)*(-4*(SMB(82)+SMB(99)-p6mup3*(2*SMB(91)+SMB(208)))
     &   +2*(SMB(201)+SMB(214)-2*SMB(253))+SMB(278))+F(300)*(2*SMB(20
     &   2)-4*(SMB(101)+SMB(257))+SMB(279))+F(516)*SMB(280)+F(1845)*(
     &   -6*(SMB(287)+SMB(288)+SMB(290))-2*Fa(37)*SMB(293)+SMB(297)+S
     &   MB(298)+SMB(300)))/s34
       result(7) = (Fa(156)*SMB(1)+Fa(157)*SMB(2)+Fa(159)*SMB(4)+Fa(16
     &   2)*SMB(9)+Fa(166)*SMB(14)+Fa(168)*SMB(16)+F(584)*SMB(19)+F(5
     &   85)*SMB(20)+Fa(171)*SMB(27)+Fa(174)*SMB(31)+F(630)*SMB(33)+F
     &   a(180)*SMB(47)+Fa(181)*SMB(48)+Fa(182)*SMB(49)+Fa(183)*SMB(5
     &   0)+Fa(184)*SMB(58)+Fa(185)*SMB(59)+F(542)*(SMB(6)+SMB(64))+F
     &   (675)*SMB(65)+F(676)*SMB(67)+F(664)*(SMB(55)+SMB(69))+F(678)
     &   *SMB(70)-Fa(196)*SMB(105)+F(700)*SMB(108)-Fa(197)*SMB(109)+F
     &   (702)*SMB(110)+F(643)*(4*SMB(41)-SMB(116)+2*(SMB(45)+SMB(118
     &   )))+F(707)*SMB(119)+Fa(200)*SMB(126)+F(713)*SMB(127)+F(714)*
     &   SMB(128)+Fa(205)*SMB(136)+8*(Fa(158)*SMB(3)+Fa(160)*SMB(7)+F
     &   a(161)*SMB(8)+Fa(163)*SMB(10)-Fa(164)*SMB(12)-Fa(165)*SMB(13
     &   )+Fa(167)*SMB(15)+Fa(169)*SMB(21)-Fa(170)*SMB(22)-Fa(172)*SM
     &   B(29)+Fa(173)*SMB(30)+Fa(175)*SMB(35)+Fa(176)*SMB(37)+Fa(177
     &   )*SMB(39)-Fa(178)*SMB(43)+Fa(179)*SMB(44)+Fa(188)*SMB(75)+Fa
     &   (190)*SMB(86)-Fa(189)*(SMB(85)+SMB(87))-Fa(187)*(SMB(74)-SMB
     &   (94)+SMB(96))+Fa(195)*SMB(104)-Fa(199)*SMB(121)+Fa(201)*SMB(
     &   130)+Fa(202)*SMB(132)+Fa(204)*SMB(134)-Fa(206)*SMB(137)+Fa(2
     &   08)*SMB(140)+Fa(209)*SMB(141)+Fa(210)*SMB(143)+Fa(211)*SMB(1
     &   49)+Fa(213)*SMB(161)+Fa(212)*(SMB(159)+SMB(163))-Fa(216)*SMB
     &   (176)+Fa(217)*SMB(179)-Fa(218)*SMB(180))+Fa(220)*SMB(184)+Fa
     &   (221)*SMB(185)+Fa(6)*(2*F(699)*SMB(124)+F(711)*SMB(125)-8*(F
     &   (688)*SMB(84)+F(708)*SMB(120)-F(715)*SMB(135))-4*(2*(F(718)*
     &   SMB(138)-F(689)*SMB(160)-F(554)*(SMB(158)+SMB(162)))+F(562)*
     &   SMB(235)+F(716)*(SMB(234)+SMB(236))+F(645)*SMB(250)))-4*(F(5
     &   83)*SMB(18)+F(596)*SMB(24)-F(597)*SMB(25)+F(598)*SMB(26)+F(6
     &   29)*SMB(32)-F(631)*SMB(34)+F(633)*SMB(36)-F(637)*SMB(38)-F(6
     &   41)*SMB(40)+F(642)*SMB(42)-F(663)*SMB(54)-F(635)*SMB(56)-F(6
     &   38)*SMB(60)-F(674)*SMB(61)-F(640)*SMB(62)-F(639)*SMB(63)-F(5
     &   18)*SMB(71)-F(687)*SMB(83)+F(523)*SMB(103)+F(699)*SMB(107)+F
     &   (688)*SMB(122)+F(708)*SMB(123)+F(720)*SMB(144)-F(717)*SMB(14
     &   8)-F(726)*SMB(177)+F(719)*SMB(181)-Fa(219)*SMB(182)-Fa(207)*
     &   (2*SMB(139)-SMB(183))-Fa(222)*SMB(186)+Fa(223)*SMB(187)-Fa(2
     &   26)*SMB(209)+Fa(227)*SMB(213)-Fa(228)*SMB(221)-Fa(225)*(SMB(
     &   195)+SMB(222))+Fa(230)*SMB(227)+Fa(229)*(SMB(226)+SMB(228))+
     &   Fa(232)*SMB(239)+Fa(231)*(SMB(238)+SMB(240))-mup3mup6*(-2*(F
     &   (561)*SMB(153)+F(721)*(SMB(152)+SMB(154)))+F(555)*(2*SMB(81)
     &   +SMB(196)-SMB(229))+F(680)*(2*SMB(76)+SMB(197)-SMB(241)))-F(
     &   680)*(SMB(79)+SMB(171)-SMB(175)-SMB(203)-2*SMB(216)-SMB(220)
     &   +SMB(245))-F(689)*SMB(247)-F(554)*(SMB(246)+SMB(248))-Fa(203
     &   )*(2*(SMB(133)-SMB(151))-SMB(261))+Fa(233)*SMB(265))+F(733)*
     &   SMB(278)+F(555)*SMB(280)+2*(F(616)*SMB(28)-F(646)*SMB(46)-F(
     &   541)*(2*SMB(5)+SMB(51))-F(519)*(2*SMB(17)+SMB(52))-F(662)*SM
     &   B(53)-F(556)*(2*SMB(11)+SMB(57))-F(544)*(2*SMB(23)+SMB(66))-
     &   F(677)*SMB(68)-F(634)*SMB(111)+Fa(191)*(2*(SMB(88)+SMB(90))+
     &   SMB(113))+Fa(192)*(2*SMB(89)+SMB(114))+Fa(198)*SMB(115)+F(55
     &   5)*(-4*SMB(206)+SMB(253)+2*(SMB(82)+SMB(99)-SMB(102)-SMB(201
     &   )-SMB(214)+SMB(233)+Fa(6)*(SMB(251)-SMB(258)))-SMB(260))+F(5
     &   62)*(2*(SMB(147)+SMB(156))+SMB(263))+F(716)*(2*(SMB(145)+SMB
     &   (155)+SMB(157))+SMB(262)+SMB(264))+Fa(194)*(-8*SMB(93)+4*(SM
     &   B(191)-SMB(212))+SMB(267))+Fa(193)*(-8*SMB(92)+4*(SMB(189)+S
     &   MB(193)-SMB(210))+SMB(266)+SMB(268))+Fa(186)*(4*SMB(73)+2*(S
     &   MB(205)-SMB(252)+SMB(259))-SMB(275)+SMB(277))+Fa(6)*(F(555)*
     &   (SMB(272)-SMB(274))+F(680)*(2*(SMB(281)-SMB(283))+SMB(284)-S
     &   MB(286)))-F(645)*(4*SMB(170)+SMB(288))-F(679)*(8*SMB(72)+4*(
     &   SMB(168)+SMB(173))+2*(SMB(78)+Fa(6)*(SMB(249)+SMB(254)))+SMB
     &   (287)+SMB(290))+Fa(224)*(6*SMB(188)+SMB(293))+Fa(215)*(-8*SM
     &   B(166)+4*(SMB(192)-SMB(219))+SMB(295))+Fa(214)*(-8*SMB(165)+
     &   4*(SMB(190)+SMB(194)-SMB(218))+SMB(294)+SMB(296)))+F(734)*SM
     &   B(299)+F(1847)*(6*(SMB(117)+SMB(200)+SMB(223)+SMB(224))-2*Fa
     &   (6)*(SMB(269)+SMB(270)+SMB(271))+SMB(297)+SMB(298)+SMB(300))
     &   +F(680)*(2*(SMB(289)-SMB(292))+SMB(302)))/s16
       result(8) = (16*(Fa(234)*SMB(3)-Fa(235)*SMB(10)+Fa(236)*(SMB(13
     &   )+SMB(15))+Fa(237)*(SMB(21)+SMB(22)))+Fa(238)*SMB(27)+s56*(2
     &   *F(754)*SMB(28)+Fa(239)*SMB(29)+Fa(240)*(SMB(30)+SMB(31))+F(
     &   764)*SMB(34)+4*(F(749)*SMB(25)-F(758)*SMB(32)-F(757)*SMB(36)
     &   ))+F(779)*SMB(54)+F(780)*SMB(55)+Fa(245)*SMB(59)+F(787)*(SMB
     &   (61)+SMB(63))+F(790)*(SMB(67)+SMB(69))+F(791)*SMB(70)+Fa(241
     &   )*(SMB(43)+SMB(75))+F(801)*SMB(103)+s56*(Fa(243)*SMB(48)+F(7
     &   78)*SMB(51)+F(1843)*(4*(SMB(38)+SMB(40))-2*(SMB(52)+SMB(53))
     &   )+Fa(244)*SMB(58)+F(788)*SMB(65)+4*(Fa(242)*SMB(47)-F(792)*S
     &   MB(71)+Fa(246)*SMB(88))-2*(F(756)*SMB(57)-F(789)*(SMB(66)+SM
     &   B(68))-F(758)*(SMB(45)-2*SMB(64)+SMB(111)+SMB(112)))+Fa(249)
     &   *SMB(113)+Fa(247)*(4*(SMB(89)+SMB(90))-2*(SMB(114)+SMB(115))
     &   )+F(808)*(SMB(116)+SMB(117)))+F(812)*SMB(127)+F(813)*SMB(128
     &   )+F(817)*(SMB(144)+SMB(146))-8*(Fa(248)*SMB(104)+Fa(250)*SMB
     &   (130)-Fa(251)*(SMB(132)+SMB(134))-Fa(252)*SMB(176)-Fa(253)*(
     &   SMB(179)+SMB(180)))+4*(F(771)*SMB(46)-F(794)*SMB(122)-F(809)
     &   *SMB(123)+F(816)*(SMB(177)+SMB(178))+Fa(254)*(SMB(209)+SMB(2
     &   11))+p6mup3*(2*(F(742)*SMB(129)+F(814)*(SMB(131)+SMB(133)))+
     &   F(741)*(SMB(208)+SMB(210))-F(819)*(SMB(217)+SMB(218)))-Fa(25
     &   5)*SMB(226)+Fa(256)*(SMB(227)+SMB(228))-Fa(257)*SMB(238)-Fa(
     &   258)*(SMB(239)+SMB(240))-F(742)*SMB(246)-F(814)*(SMB(247)+SM
     &   B(248))+Fa(259)*SMB(261)+Fa(260)*SMB(265))+F(741)*(-2*(SMB(2
     &   53)+SMB(257))+SMB(278)+SMB(279))+p5mup6*(F(735)*SMB(1)+F(740
     &   )*SMB(7)+F(765)*SMB(35)+F(766)*(SMB(37)+SMB(39))+F(777)*(SMB
     &   (49)+SMB(50))+8*(F(794)*SMB(84)+F(809)*SMB(120))+F(810)*SMB(
     &   124)+F(811)*SMB(125)+2*(4*(F(815)*(SMB(135)+SMB(136))+F(816)
     &   *(SMB(138)+SMB(139))+F(742)*SMB(158)+F(814)*(SMB(160)+SMB(16
     &   2)))+2*(F(828)*SMB(234)+F(743)*(SMB(235)+SMB(236))+F(834)*SM
     &   B(249)+F(771)*(SMB(250)+SMB(254)))+F(741)*(2*(SMB(251)+SMB(2
     &   55))+SMB(272)+SMB(273))-F(819)*(2*(SMB(281)+SMB(282))+SMB(28
     &   4)+SMB(285))))+2*(F(793)*SMB(78)+F(828)*SMB(262)+F(743)*(SMB
     &   (263)+SMB(264))-Fa(261)*SMB(266)+Fa(262)*(SMB(267)+SMB(268))
     &   -F(834)*SMB(287)-F(771)*(SMB(288)+SMB(290))+F(819)*(SMB(289)
     &   +SMB(291))+Fa(263)*SMB(293)-Fa(264)*SMB(294)-Fa(265)*(SMB(29
     &   5)+SMB(296)))+F(1844)*(2*p5mup6*(SMB(269)+SMB(270)+SMB(271))
     &   +SMB(297)+SMB(298)+SMB(300))+F(839)*(SMB(299)+SMB(301)))/s56
       result(9) = Fa(281)*SMB(27)+F(1134)*SMB(32)+Fa(285)*SMB(35)+Fa(
     &   290)*SMB(47)+Fa(291)*SMB(48)+Fa(292)*SMB(49)+Fa(293)*SMB(50)
     &   +Fa(295)*SMB(59)+F(1246)*SMB(62)+F(1248)*SMB(64)+F(1249)*SMB
     &   (65)+F(1251)*SMB(67)+F(1253)*SMB(69)+F(1297)*SMB(83)+Fa(310)
     &   *SMB(105)+F(1376)*SMB(107)+F(1377)*SMB(108)+Fa(311)*SMB(109)
     &   +F(1402)*SMB(110)+Fa(312)*SMB(113)+Fa(313)*SMB(114)+F(1420)*
     &   SMB(116)+F(1421)*SMB(117)+F(1423)*SMB(119)+Fa(317)*SMB(124)+
     &   Fa(318)*SMB(125)+Fa(319)*SMB(126)+F(1450)*SMB(127)+F(1451)*S
     &   MB(128)+F(1544)*SMB(170)+8*(Fa(266)*SMB(1)+Fa(267)*SMB(2)-Fa
     &   (268)*SMB(3)+Fa(269)*SMB(4)+Fa(270)*SMB(7)+Fa(271)*SMB(8)+Fa
     &   (272)*SMB(9)-Fa(273)*SMB(10)-Fa(274)*SMB(12)+Fa(275)*SMB(13)
     &   +Fa(276)*SMB(14)+Fa(277)*SMB(15)+Fa(278)*SMB(16)-Fa(279)*SMB
     &   (21)-Fa(280)*SMB(22)+Fa(282)*SMB(29)-Fa(283)*SMB(30)+Fa(284)
     &   *SMB(31)+Fa(286)*SMB(37)+Fa(287)*SMB(39)+Fa(288)*SMB(43)-Fa(
     &   289)*SMB(44)-Fa(296)*SMB(73)+Fa(297)*SMB(74)+Fa(298)*SMB(75)
     &   -Fa(299)*SMB(84)+Fa(300)*SMB(85)-Fa(301)*SMB(86)-Fa(302)*SMB
     &   (87)-Fa(303)*SMB(88)+Fa(304)*SMB(89)+Fa(305)*SMB(90)-Fa(306)
     &   *SMB(94)+Fa(307)*SMB(95)+Fa(308)*SMB(96)+Fa(309)*SMB(104)-Fa
     &   (315)*SMB(120)+Fa(316)*SMB(121)+Fa(320)*SMB(130)+Fa(321)*SMB
     &   (132)+Fa(322)*SMB(133)-Fa(323)*SMB(134)-Fa(324)*SMB(135)-Fa(
     &   325)*SMB(136)+Fa(326)*SMB(137)-Fa(327)*SMB(138)-Fa(328)*SMB(
     &   139)+Fa(329)*SMB(140)+Fa(330)*SMB(141)+Fa(331)*SMB(142)+Fa(3
     &   32)*SMB(143)+Fa(333)*SMB(149)+F(1319)*SMB(150)+Fa(334)*SMB(1
     &   51)+mup3mup6*(F(1292)*SMB(77)-F(1511)*SMB(152)-F(1512)*SMB(1
     &   53)-F(1520)*SMB(154))-Fa(335)*SMB(158)-Fa(336)*SMB(159)+Fa(3
     &   37)*SMB(160)+Fa(338)*SMB(161)+Fa(339)*SMB(162)+Fa(340)*SMB(1
     &   63)+p6mup3*(F(1296)*SMB(91)+F(1304)*SMB(92)+F(1306)*SMB(93)+
     &   F(1303)*SMB(129)-F(1305)*SMB(131)-F(1294)*SMB(164)-F(1542)*S
     &   MB(165)-F(1328)*SMB(166))-F(1851)*SMB(173)+Fa(341)*SMB(176)-
     &   Fa(342)*SMB(179)-Fa(343)*SMB(180))+12*(F(1304)*SMB(146)+F(15
     &   60)*SMB(181))+Fa(361)*SMB(213)+Fa(365)*SMB(222)+Fa(376)*SMB(
     &   240)+Fa(381)*SMB(254)+F(1296)*SMB(278)+F(1304)*SMB(279)+F(13
     &   06)*SMB(280)+2*(F(1205)*SMB(45)+F(1206)*SMB(46)-F(1227)*SMB(
     &   51)-F(1228)*SMB(52)-F(1229)*SMB(53)+F(1231)*SMB(55)+F(1233)*
     &   SMB(57)+F(1250)*SMB(66)+F(1252)*SMB(68)+F(1293)*SMB(78)+Fa(3
     &   14)*SMB(115)-F(1643)*SMB(215)+F(1628)*(SMB(200)+SMB(223)+SMB
     &   (224))+F(1319)*(SMB(230)+SMB(231)+SMB(232))-F(1320)*(SMB(242
     &   )+SMB(243)+SMB(244))+F(1296)*(SMB(201)+SMB(214)-SMB(253))+F(
     &   1304)*(SMB(202)-SMB(257))+F(1306)*(SMB(233)-SMB(260))+F(1808
     &   )*SMB(262)+F(1809)*SMB(263)+F(1810)*SMB(264)+Fa(350)*SMB(266
     &   )+Fa(352)*SMB(267)-Fa(354)*SMB(268)+Fa(388)*(SMB(269)+SMB(27
     &   0)+SMB(271))-Fa(379)*SMB(272)-Fa(382)*SMB(273)-Fa(384)*SMB(2
     &   74)+Fa(380)*SMB(275)+Fa(383)*SMB(276)+Fa(385)*SMB(277)-Fa(38
     &   9)*SMB(284)-Fa(390)*SMB(285)+Fa(391)*SMB(286)-F(1831)*SMB(28
     &   7)-F(1832)*SMB(288)-F(1294)*(SMB(203)+SMB(220)-SMB(289))-F(1
     &   833)*SMB(290)+F(1542)*(-SMB(204)+SMB(291))+F(1328)*(-SMB(245
     &   )+SMB(292))+Fa(392)*SMB(293)+Fa(351)*SMB(294)-Fa(353)*SMB(29
     &   5)-Fa(355)*SMB(296))+F(1851)*SMB(297)+F(1851)*SMB(298)+F(183
     &   5)*SMB(299)+F(1851)*SMB(300)+F(1836)*SMB(301)+F(1560)*SMB(30
     &   2)+4*(F(901)*SMB(5)-F(902)*SMB(6)+F(935)*SMB(11)+F(1033)*SMB
     &   (17)+F(1034)*SMB(18)-F(1035)*SMB(19)+F(1036)*SMB(20)+F(1088)
     &   *SMB(23)+F(1089)*SMB(24)+F(1090)*SMB(25)+F(1091)*SMB(26)+F(1
     &   118)*SMB(28)+F(1135)*SMB(33)+F(1136)*SMB(34)+F(1142)*SMB(36)
     &   +F(1148)*SMB(38)+F(1154)*SMB(40)+F(1155)*SMB(41)+F(1156)*SMB
     &   (42)+F(1230)*SMB(54)+F(1232)*SMB(56)-Fa(294)*SMB(58)+F(1244)
     &   *SMB(60)+F(1245)*SMB(61)+F(1247)*SMB(63)-F(1254)*SMB(70)-F(1
     &   255)*SMB(71)+F(1256)*SMB(72)+F(1295)*SMB(80)-F(1319)*(SMB(97
     &   )+SMB(98)+SMB(100))-F(1304)*SMB(101)-F(1306)*SMB(102)+F(1329
     &   )*SMB(103)-F(1375)*SMB(106)+F(1403)*SMB(111)-F(1404)*SMB(112
     &   )-F(1422)*SMB(118)-F(1433)*SMB(122)-F(1434)*SMB(123)+F(1510)
     &   *SMB(144)+F(1513)*SMB(148)+F(1511)*(SMB(145)+SMB(155))+F(151
     &   2)*(SMB(147)+SMB(156))+F(1520)*SMB(157)-F(1543)*SMB(168)+F(1
     &   320)*(SMB(167)+SMB(169)+SMB(172))+F(1542)*SMB(174)+F(1328)*S
     &   MB(175)-F(1550)*SMB(177)+F(1551)*SMB(178)+Fa(344)*SMB(182)+F
     &   a(345)*SMB(183)-Fa(346)*SMB(184)+Fa(347)*SMB(185)-Fa(348)*SM
     &   B(186)+Fa(340)*SMB(187)+Fa(349)*SMB(188)+Fa(350)*SMB(189)+Fa
     &   (351)*SMB(190)+Fa(352)*SMB(191)-Fa(353)*SMB(192)-Fa(354)*SMB
     &   (193)-Fa(355)*SMB(194)-Fa(356)*SMB(195)-Fa(357)*SMB(205)-F(1
     &   292)*SMB(207)-F(1296)*(SMB(82)+SMB(99)+mup3mup6*SMB(196)-p6m
     &   up3*SMB(208))+Fa(358)*SMB(209)-Fa(359)*SMB(210)+Fa(307)*SMB(
     &   211)+Fa(360)*SMB(212)+F(1294)*(SMB(79)+SMB(171)+mup3mup6*SMB
     &   (197)-p6mup3*SMB(217))-Fa(362)*SMB(218)+Fa(363)*SMB(219)+Fa(
     &   364)*SMB(221)-Fa(366)*SMB(225)-Fa(367)*SMB(226)-Fa(368)*SMB(
     &   227)-Fa(369)*SMB(228)-Fa(370)*SMB(234)+Fa(371)*SMB(235)+Fa(3
     &   72)*SMB(236)-Fa(373)*SMB(237)-Fa(374)*SMB(238)+Fa(375)*SMB(2
     &   39)+mup3mup6*(-(F(1304)*SMB(198))+F(1542)*SMB(199)-F(1306)*S
     &   MB(229)+F(1328)*SMB(241))-F(1303)*SMB(246)+F(1305)*SMB(247)+
     &   F(1307)*SMB(248)+Fa(377)*SMB(249)-Fa(378)*SMB(250)-Fa(379)*S
     &   MB(251)+Fa(380)*SMB(252)-Fa(382)*SMB(255)+Fa(383)*SMB(256)-F
     &   a(384)*SMB(258)+Fa(385)*SMB(259)+Fa(386)*SMB(261)-Fa(387)*SM
     &   B(265)-Fa(389)*SMB(281)-Fa(390)*SMB(282)+Fa(391)*SMB(283)-F(
     &   1853)*SMB(303))
       If (ngluon.eq.0) then
       resultn=(result(1)+result(2)+result(3)+result(4)+result(5)+resu
     &   lt(6)+result(7)+result(8)+result(9))
       elseif (ngluon.eq.1) then
       If (posgluon.eq.3) then
       result(1)=(result(1)+result(2)+result(3)+result(5)+result(6))
       result(2)=(result(4))
       elseIf (posgluon.eq.4) then
       result(1)=(result(1)+result(2)+result(3)+result(4)+result(6))
       result(2)=(result(5))
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     &   . Look to the heading for explanation" 
       endif
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     &   . Look to the heading for explanation" 
       endif
       resultborn = -(((2*p5mup6*(2*SMB(1)+SMB(48)+SMB(49)+SMB(50))+SM
     &   B(55)+SMB(67)+SMB(69)+2*SMB(70))/(s23*s56)+(4*(p4mup6+p5mup6
     &   )*SMB(2)+2*SMB(6)-4*mup3mup6*(SMB(4)+SMB(14)+SMB(16))+2*SMB(
     &   19)+2*SMB(20)+2*(p6sq-s123-s16+s45)*SMB(25)+2*SMB(33)-2*p1mu
     &   p6*(2*SMB(1)-2*SMB(2)+SMB(48)+SMB(49)+SMB(50))-2*p6mup6*(2*S
     &   MB(1)-2*SMB(2)+SMB(48)+SMB(49)+SMB(50))+SMB(55)+2*SMB(64)+SM
     &   B(67)+SMB(69)+2*SMB(70))/(s16*s23)+(-2*s56*SMB(25)-2*SMB(54)
     &   +SMB(55)+2*(2*p6mup3*(SMB(29)+SMB(30)+SMB(31))+p5mup6*(-2*SM
     &   B(7)-2*(p1mup3+p5mup3+p6mup3)*SMB(27)+2*SMB(35)+2*SMB(37)+2*
     &   SMB(39)+SMB(48)+SMB(49)+SMB(50))-(p1mup3+p5mup3+p6mup3)*SMB(
     &   59))-2*SMB(61)-2*SMB(63)+SMB(67)+SMB(69)-2*SMB(103))/(s34*s5
     &   6)+(-4*p2mup6*SMB(2)+2*SMB(6)-4*(p1mup3+p5mup3+p6mup3)*SMB(9
     &   )+2*SMB(19)+2*SMB(20)+2*(p6sq-s12-s16+s345-s56)*SMB(25)+2*SM
     &   B(33)+p1mup6*(4*SMB(2)+4*SMB(7)+4*(p1mup3+p5mup3+p6mup3)*SMB
     &   (27)-2*(2*SMB(35)+2*SMB(37)+2*SMB(39)+SMB(48)+SMB(49)+SMB(50
     &   )))-2*SMB(54)+SMB(55)+2*(2*p5mup6*SMB(2)+p6mup6*(2*SMB(2)+2*
     &   SMB(7)+2*(p1mup3+p5mup3+p6mup3)*SMB(27)-2*SMB(35)-2*SMB(37)-
     &   2*SMB(39)-SMB(48)-SMB(49)-SMB(50))+(p1mup3+p5mup3+p6mup3)*(2
     &   *SMB(31)-SMB(59)))-2*SMB(61)-2*SMB(63)+2*SMB(64)+SMB(67)+SMB
     &   (69)-2*SMB(103))/(s16*s34))/s234)
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
