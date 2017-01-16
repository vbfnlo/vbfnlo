c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:29/5/2012
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
       subroutine HexBox66(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,barpsi_p
     &   2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,resultn,
     &   resultborn)
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
      Complex*16 SMB(291),Fa(259),F(1314)
       Complex*16 SMB0(31),SMB11(0:3),SMB12(0:3),SMB13(0:3),SMB14(0:3)
     &   ,SMB15(0:3),SMB16(0:3),SMB17(0:3),SMB18(0:3)
      Complex*16 SMB1(127)

      Complex*16 barpsi_p1_P(2),psi_p6_P(2),barpsi_p2_P(2),psi_p5_P(2) 
       Complex*16 barpsi_p1_M(2),psi_p6_M(2),barpsi_p2_M(2),psi_p5_M(2)
       Complex*16 barpsi_p1(4),psi_p6(4),barpsi_p2(4),psi_p5(4)
       Real*8 P(1362)
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
       Complex*16   mup3(0:3),mup4(0:3)
       Complex*16   p1mup3, p1mup4, p2mup3, p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4
       Real*8 dotrr
       Complex*16 B0fin,C0fin,D0fin,E0fin,F0fin
       EXTERNAL dotrr,B0fin,C0fin,D0fin,E0fin,F0fin
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
      Complex*16 dotrc,dotcc,result(6),resultn,resultborn
       Real*8 musq
      EXTERNAL   dotrc,dotcc

       COMMON/HexBoxFaFunctions/Fa
       COMMON/HexBoxFhlFunctions/F
      Save/HexBoxFhlFunctions/
       COMMON/HexBoxPFunctions/P
      Save/HexBoxPFunctions/
       COMMON/HexBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s
     &   34,s45,s56,s16,s123,s234,s345
       COMMON/HexBoxFVALFunctions/F0123456R,F123456R,F0123456I,F123456
     &   I
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
       Integer ngluon, posgluon
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
       B012=B0fin(p1sq,musq)
       B012R=Dble(B012)
       B012I=DImag(B012)
       B023=B0fin(p2sq,musq)
       B023R=Dble(B023)
       B023I=DImag(B023)
       B034=B0fin(p3sq,musq)
       B034R=Dble(B034)
       B034I=DImag(B034)
       B045=B0fin(p4sq,musq)
       B045R=Dble(B045)
       B045I=DImag(B045)
       B056=B0fin(p5sq,musq)
       B056R=Dble(B056)
       B056I=DImag(B056)
       B016=B0fin(p6sq,musq)
       B016R=Dble(B016)
       B016I=DImag(B016)
       B013=B0fin(s12,musq)
       B013R=Dble(B013)
       B013I=DImag(B013)
       B014=B0fin(s123,musq)
       B014R=Dble(B014)
       B014I=DImag(B014)
       B026=B0fin(s16,musq)
       B026R=Dble(B026)
       B026I=DImag(B026)
       B024=B0fin(s23,musq)
       B024R=Dble(B024)
       B024I=DImag(B024)
       B025=B0fin(s234,musq)
       B025R=Dble(B025)
       B025I=DImag(B025)
       B035=B0fin(s34,musq)
       B035R=Dble(B035)
       B035I=DImag(B035)
       B036=B0fin(s345,musq)
       B036R=Dble(B036)
       B036I=DImag(B036)
       B046=B0fin(s45,musq)
       B046R=Dble(B046)
       B046I=DImag(B046)
       B015=B0fin(s56,musq)
       B015R=Dble(B015)
       B015I=DImag(B015)
       C0123=C0fin(p1sq,p2sq,s12,musq)
       call tens_red3_new_Re_Com(p1sq,p2sq,s12,B023,B013,B012,C0123,C0
     &   123R,C0123I,C123R,C123I)
       C0126=C0fin(p1sq,s16,p6sq,musq)
       call tens_red3_new_Re_Com(p1sq,s16,p6sq,B026,B016,B012,C0126,C0
     &   126R,C0126I,C126R,C126I)
       C0124=C0fin(p1sq,s23,s123,musq)
       call tens_red3_new_Re_Com(p1sq,s23,s123,B024,B014,B012,C0124,C0
     &   124R,C0124I,C124R,C124I)
       C0125=C0fin(p1sq,s234,s56,musq)
       call tens_red3_new_Re_Com(p1sq,s234,s56,B025,B015,B012,C0125,C0
     &   125R,C0125I,C125R,C125I)
       C0234=C0fin(p2sq,p3sq,s23,musq)
       call tens_red3_new_Re_Com(p2sq,p3sq,s23,B034,B024,B023,C0234,C0
     &   234R,C0234I,C234R,C234I)
       C0235=C0fin(p2sq,s34,s234,musq)
       call tens_red3_new_Re_Com(p2sq,s34,s234,B035,B025,B023,C0235,C0
     &   235R,C0235I,C235R,C235I)
       C0236=C0fin(p2sq,s345,s16,musq)
       call tens_red3_new_Re_Com(p2sq,s345,s16,B036,B026,B023,C0236,C0
     &   236R,C0236I,C236R,C236I)
       C0345=C0fin(p3sq,p4sq,s34,musq)
       call tens_red3_new_Re_Com(p3sq,p4sq,s34,B045,B035,B034,C0345,C0
     &   345R,C0345I,C345R,C345I)
       C0346=C0fin(p3sq,s45,s345,musq)
       call tens_red3_new_Re_Com(p3sq,s45,s345,B046,B036,B034,C0346,C0
     &   346R,C0346I,C346R,C346I)
       C0456=C0fin(p4sq,p5sq,s45,musq)
       call tens_red3_new_Re_Com(p4sq,p5sq,s45,B056,B046,B045,C0456,C0
     &   456R,C0456I,C456R,C456I)
       C0134=C0fin(s12,p3sq,s123,musq)
       call tens_red3_new_Re_Com(s12,p3sq,s123,B034,B014,B013,C0134,C0
     &   134R,C0134I,C134R,C134I)
       C0135=C0fin(s12,s34,s56,musq)
       call tens_red3_new_Re_Com(s12,s34,s56,B035,B015,B013,C0135,C013
     &   5R,C0135I,C135R,C135I)
       C0136=C0fin(s12,s345,p6sq,musq)
       call tens_red3_new_Re_Com(s12,s345,p6sq,B036,B016,B013,C0136,C0
     &   136R,C0136I,C136R,C136I)
       C0145=C0fin(s123,p4sq,s56,musq)
       call tens_red3_new_Re_Com(s123,p4sq,s56,B045,B015,B014,C0145,C0
     &   145R,C0145I,C145R,C145I)
       C0146=C0fin(s123,s45,p6sq,musq)
       call tens_red3_new_Re_Com(s123,s45,p6sq,B046,B016,B014,C0146,C0
     &   146R,C0146I,C146R,C146I)
       C0245=C0fin(s23,p4sq,s234,musq)
       call tens_red3_new_Re_Com(s23,p4sq,s234,B045,B025,B024,C0245,C0
     &   245R,C0245I,C245R,C245I)
       C0246=C0fin(s23,s45,s16,musq)
       call tens_red3_new_Re_Com(s23,s45,s16,B046,B026,B024,C0246,C024
     &   6R,C0246I,C246R,C246I)
       C0256=C0fin(s234,p5sq,s16,musq)
       call tens_red3_new_Re_Com(s234,p5sq,s16,B056,B026,B025,C0256,C0
     &   256R,C0256I,C256R,C256I)
       C0356=C0fin(s34,p5sq,s345,musq)
       call tens_red3_new_Re_Com(s34,p5sq,s345,B056,B036,B035,C0356,C0
     &   356R,C0356I,C356R,C356I)
       C0156=C0fin(s56,p5sq,p6sq,musq)
       call tens_red3_new_Re_Com(s56,p5sq,p6sq,B056,B016,B015,C0156,C0
     &   156R,C0156I,C156R,C156I)
       D01234=D0fin(s12,s23,p1sq,p2sq,p3sq,s123,musq)
       call tens_red4_new_Re_Com(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
     &   C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C0134I,C
     &   0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D01234I,
     &   D1234R,D1234I)
       D01235=D0fin(s12,s234,p1sq,p2sq,s34,s56,musq)
       call tens_red4_new_Re_Com(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
     &   4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0235I
     &   ,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235,D01235R
     &   ,D01235I,D1235R,D1235I)
       D01236=D0fin(s12,s16,p1sq,p2sq,s345,p6sq,musq)
       call tens_red4_new_Re_Com(p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p5,p2
     &   p3+p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C126R,C
     &   123R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123I,D01
     &   236,D01236R,D01236I,D1236R,D1236I)
       D01245=D0fin(s123,s234,p1sq,s23,p4sq,s56,musq)
       call tens_red4_new_Re_Com(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
     &   4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0245I
     &   ,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245,D01245R
     &   ,D01245I,D1245R,D1245I)
       D01246=D0fin(s123,s16,p1sq,s23,s45,p6sq,musq)
       call tens_red4_new_Re_Com(p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,p2p4
     &   +p2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R,C126
     &   R,C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C124I,
     &   D01246,D01246R,D01246I,D1246R,D1246I)
       D01256=D0fin(s56,s16,p1sq,s234,p5sq,p6sq,musq)
       call tens_red4_new_Re_Com(p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p5,p2
     &   p5+p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C126R,C
     &   125R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125I,D01
     &   256,D01256R,D01256I,D1256R,D1256I)
       D02345=D0fin(s23,s34,p2sq,p3sq,p4sq,s234,musq)
       call tens_red4_new_Re_Com(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
     &   C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C0245I,C
     &   0235I,C0234I,C345I,C245I,C235I,C234I,D02345,D02345R,D02345I,
     &   D2345R,D2345I)
       D02346=D0fin(s23,s345,p2sq,p3sq,s45,s16,musq)
       call tens_red4_new_Re_Com(p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4+p3p
     &   5,C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C0346I
     &   ,C0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346,D02346R
     &   ,D02346I,D2346R,D2346I)
       D02356=D0fin(s234,s345,p2sq,s34,p5sq,s16,musq)
       call tens_red4_new_Re_Com(p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5+p4p
     &   5,C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C0356I
     &   ,C0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356,D02356R
     &   ,D02356I,D2356R,D2356I)
       D03456=D0fin(s34,s45,p3sq,p4sq,p5sq,s345,musq)
       call tens_red4_new_Re_Com(p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C0456R,
     &   C0356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C0356I,C
     &   0346I,C0345I,C456I,C356I,C346I,C345I,D03456,D03456R,D03456I,
     &   D3456R,D3456I)
       D01345=D0fin(s123,s34,s12,p3sq,p4sq,s56,musq)
       call tens_red4_new_Re_Com(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
     &   4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0345I
     &   ,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345,D01345R
     &   ,D01345I,D1345R,D1345I)
       D01346=D0fin(s123,s345,s12,p3sq,s45,p6sq,musq)
       call tens_red4_new_Re_Com(s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+p2p4
     &   +p2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R,C136
     &   R,C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C134I,
     &   D01346,D01346R,D01346I,D1346R,D1346I)
       D01356=D0fin(s56,s345,s12,s34,p5sq,p6sq,musq)
       call tens_red4_new_Re_Com(s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,p1p5
     &   +p2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R,C136
     &   R,C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C135I,
     &   D01356,D01356R,D01356I,D1356R,D1356I)
       D01456=D0fin(s56,s45,s123,p4sq,p5sq,p6sq,musq)
       call tens_red4_new_Re_Com(s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p5+p2
     &   p5+p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C146R,C
     &   145R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145I,D01
     &   456,D01456R,D01456I,D1456R,D1456I)
       D02456=D0fin(s234,s45,s23,p4sq,p5sq,s16,musq)
       call tens_red4_new_Re_Com(s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5,p4p
     &   5,C0456R,C0256R,C0246R,C0245R,C456R,C256R,C246R,C245R,C0456I
     &   ,C0256I,C0246I,C0245I,C456I,C256I,C246I,C245I,D02456,D02456R
     &   ,D02456I,D2456R,D2456I)
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
       call HexBoxFFhl1(F(1))
       call HexBoxFFhl2(F(66))
       call HexBoxFFhl3(F(131))
       call HexBoxFFhl4(F(196))
       call HexBoxFFhl5(F(261))
       call HexBoxFFhl6(F(326))
       call HexBoxFFhl7(F(391))
       call HexBoxFFhl8(F(456))
       call HexBoxFFhl9(F(521))
       call HexBoxFFhl10(F(586))
       call HexBoxFFhl11(F(651))
       call HexBoxFFhl12(F(716))
       call HexBoxFFhl13(F(781))
       call HexBoxFFhl14(F(846))
       call HexBoxFFhl15(F(911))
       call HexBoxFFhl16(F(976))
       call HexBoxFFhl17(F(1041))
       call HexBoxFFhl18(F(1106))
       call HexBoxFFhl19(F(1171))
       call HexBoxFFhl20(F(1236))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup3 = dotrc(p1,mup3)
       p1mup4 = dotrc(p1,mup4)
       p2mup3 = dotrc(p2,mup3)
       p2mup4 = dotrc(p2,mup4)
       p3mup3 = dotrc(p3,mup3)
       p3mup4 = dotrc(p3,mup4)
       p4mup3 = dotrc(p4,mup3)
       p4mup4 = dotrc(p4,mup4)
       p5mup3 = dotrc(p5,mup3)
       p5mup4 = dotrc(p5,mup4)
       p6mup3 = dotrc(p6,mup3)
       p6mup4 = dotrc(p6,mup4)
       mup3mup4 = dotcc(mup3,mup4)
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call HexBoxFa1(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup3
     &   ,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(1))
       call HexBoxFa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup3
     &   ,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(13))
       call HexBoxFa3(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup3
     &   ,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(25))
       call HexBoxFa4(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup3
     &   ,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(37))
       call HexBoxFa5(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup3
     &   ,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(49))
       call HexBoxFa6(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup3
     &   ,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(61))
       call HexBoxFa7(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup3
     &   ,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(73))
       call HexBoxFa8(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup3
     &   ,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(85))
       call HexBoxFa9(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup3
     &   ,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(97))
       call HexBoxFa10(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(109))
       call HexBoxFa11(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(121))
       call HexBoxFa12(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(133))
       call HexBoxFa13(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(145))
       call HexBoxFa14(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(157))
       call HexBoxFa15(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(169))
       call HexBoxFa16(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(181))
       call HexBoxFa17(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(193))
       call HexBoxFa18(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(205))
       call HexBoxFa19(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(217))
       call HexBoxFa20(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(229))
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
       psi_p6_M(1)=psi_p6(1)
       psi_p6_M(2)=psi_p6(2)
       psi_p6_P(1)=psi_p6(3)
       psi_p6_P(2)=psi_p6(4)
       psi_p5_M(1)=psi_p5(1)
       psi_p5_M(2)=psi_p5(2)
       psi_p5_P(1)=psi_p5(3)
       psi_p5_P(2)=psi_p5(4)
       barpsi_p1_P(1)=barpsi_p1(1)
       barpsi_p1_P(2)=barpsi_p1(2)
       barpsi_p1_M(1)=barpsi_p1(3)
       barpsi_p1_M(2)=barpsi_p1(4)
       barpsi_p2_P(1)=barpsi_p2(1)
       barpsi_p2_P(2)=barpsi_p2(2)
       barpsi_p2_M(1)=barpsi_p2(3)
       barpsi_p2_M(2)=barpsi_p2(4)
c************************************************************************************
c************************************************************************************
       SMB0(1) = SC1c(barpsi_p1_P,mup3,psi_p6_P,+1)
       SMB0(2) = SC1c(barpsi_p2_P,mup4,psi_p5_P,+1)
       SMB0(3) = SC1r(barpsi_p2_P,p1,psi_p5_P,+1)
       SMB0(4) = SC1r(barpsi_p2_P,p4,psi_p5_P,+1)
       SMB0(5) = SC1r(barpsi_p2_P,p6,psi_p5_P,+1)
       SMB0(6) = SC3crr(barpsi_p2_P,mup4,p4,p1,psi_p5_P,+1)
       SMB0(7) = SC3crr(barpsi_p2_P,mup4,p6,p1,psi_p5_P,+1)
       SMB0(8) = SC3crr(barpsi_p2_P,mup4,p6,p4,psi_p5_P,+1)
       SMB0(9) = SC3rrr(barpsi_p2_P,p6,p4,p1,psi_p5_P,+1)
       SMB0(10) = SC1c(barpsi_p1_P,mup4,psi_p6_P,+1)
       SMB0(11) = SC1c(barpsi_p2_P,mup3,psi_p5_P,+1)
       SMB0(12) = SC3crr(barpsi_p2_P,mup3,p4,p1,psi_p5_P,+1)
       SMB0(13) = SC3crr(barpsi_p2_P,mup3,p6,p1,psi_p5_P,+1)
       SMB0(14) = SC3crr(barpsi_p2_P,mup3,p6,p4,psi_p5_P,+1)
       SMB0(15) = SC1r(barpsi_p1_P,p4,psi_p6_P,+1)
       SMB0(16) = SC3ccr(barpsi_p2_P,mup4,mup3,p6,psi_p5_P,+1)
       SMB0(17) = SC1r(barpsi_p1_P,p2,psi_p6_P,+1)
       SMB0(18) = SC3ccr(barpsi_p2_P,mup4,mup3,p1,psi_p5_P,+1)
       SMB0(19) = SC3ccr(barpsi_p2_P,mup4,mup3,p4,psi_p5_P,+1)
       SMB0(20) = SC5ccrrr(barpsi_p2_P,mup4,mup3,p6,p4,p1,psi_p5_P,+1)
       SMB0(21) = SC1r(barpsi_p1_P,p5,psi_p6_P,+1)
       SMB0(22) = SC3crr(barpsi_p1_P,mup4,p4,p2,psi_p6_P,+1)
       SMB0(23) = SC3crr(barpsi_p1_P,mup4,p5,p2,psi_p6_P,+1)
       SMB0(24) = SC3crr(barpsi_p1_P,mup4,p5,p4,psi_p6_P,+1)
       SMB0(25) = SC3rrr(barpsi_p1_P,p5,p4,p2,psi_p6_P,+1)
       SMB0(26) = SC3crr(barpsi_p1_P,mup3,p4,p2,psi_p6_P,+1)
       SMB0(27) = SC3crr(barpsi_p1_P,mup3,p5,p2,psi_p6_P,+1)
       SMB0(28) = SC3crr(barpsi_p1_P,mup3,p5,p4,psi_p6_P,+1)
       SMB0(29) = SC3ccr(barpsi_p1_P,mup4,mup3,p2,psi_p6_P,+1)
       SMB0(30) = SC3ccr(barpsi_p1_P,mup4,mup3,p4,psi_p6_P,+1)
       SMB0(31) = SC3ccr(barpsi_p1_P,mup4,mup3,p5,psi_p6_P,+1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB11(i) = SC1c(barpsi_p1_P,v1,psi_p6_P,+1)
       SMB12(i) = SC1c(barpsi_p2_P,v1,psi_p5_P,+1)
       SMB13(i) = SC3crc(barpsi_p1_P,mup3,p4,v1,psi_p6_P,+1)
       SMB14(i) = SC3ccc(barpsi_p1_P,mup4,mup3,v1,psi_p6_P,+1)
       SMB15(i) = SC3crc(barpsi_p1_P,mup4,p4,v1,psi_p6_P,+1)
       SMB16(i) = SC3ccc(barpsi_p2_P,mup4,mup3,v1,psi_p5_P,+1)
       SMB17(i) = SC3crc(barpsi_p2_P,mup3,p4,v1,psi_p5_P,+1)
       SMB18(i) = SC3crc(barpsi_p2_P,mup4,p4,v1,psi_p5_P,+1)
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
       SMB1(21) = SMB0(5)*SMB0(15)
       SMB1(22) = SMB0(15)*SMB0(16)
       SMB1(23) = SMB0(11)*SMB0(17)
       SMB1(24) = SMB0(2)*SMB0(17)
       SMB1(25) = SMB0(3)*SMB0(17)
       SMB1(26) = SMB0(4)*SMB0(17)
       SMB1(27) = SMB0(5)*SMB0(17)
       SMB1(28) = SMB0(12)*SMB0(17)
       SMB1(29) = SMB0(13)*SMB0(17)
       SMB1(30) = SMB0(14)*SMB0(17)
       SMB1(31) = SMB0(17)*SMB0(18)
       SMB1(32) = SMB0(17)*SMB0(19)
       SMB1(33) = SMB0(16)*SMB0(17)
       SMB1(34) = SMB0(6)*SMB0(17)
       SMB1(35) = SMB0(8)*SMB0(17)
       SMB1(36) = SMB0(9)*SMB0(17)
       SMB1(37) = SMB0(17)*SMB0(20)
       SMB1(38) = SMB0(11)*SMB0(15)
       SMB1(39) = SMB0(2)*SMB0(15)
       SMB1(40) = SMB0(3)*SMB0(15)
       SMB1(41) = SMB0(4)*SMB0(15)
       SMB1(42) = SMB0(12)*SMB0(15)
       SMB1(43) = SMB0(13)*SMB0(15)
       SMB1(44) = SMB0(14)*SMB0(15)
       SMB1(45) = SMB0(15)*SMB0(18)
       SMB1(46) = SMB0(15)*SMB0(19)
       SMB1(47) = SMB0(6)*SMB0(15)
       SMB1(48) = SMB0(7)*SMB0(15)
       SMB1(49) = SMB0(8)*SMB0(15)
       SMB1(50) = SMB0(9)*SMB0(15)
       SMB1(51) = SMB0(15)*SMB0(20)
       SMB1(52) = SMB0(11)*SMB0(21)
       SMB1(53) = SMB0(2)*SMB0(21)
       SMB1(54) = SMB0(3)*SMB0(21)
       SMB1(55) = SMB0(4)*SMB0(21)
       SMB1(56) = SMB0(5)*SMB0(21)
       SMB1(57) = SMB0(12)*SMB0(21)
       SMB1(58) = SMB0(13)*SMB0(21)
       SMB1(59) = SMB0(14)*SMB0(21)
       SMB1(60) = SMB0(18)*SMB0(21)
       SMB1(61) = SMB0(19)*SMB0(21)
       SMB1(62) = SMB0(16)*SMB0(21)
       SMB1(63) = SMB0(6)*SMB0(21)
       SMB1(64) = SMB0(8)*SMB0(21)
       SMB1(65) = SMB0(9)*SMB0(21)
       SMB1(66) = SMB0(20)*SMB0(21)
       SMB1(67) = SMB0(11)*SMB0(22)
       SMB1(68) = SMB0(11)*SMB0(23)
       SMB1(69) = SMB0(11)*SMB0(24)
       SMB1(70) = SMB0(11)*SMB0(25)
       SMB1(71) = SMB0(2)*SMB0(26)
       SMB1(72) = SMB0(2)*SMB0(27)
       SMB1(73) = SMB0(2)*SMB0(28)
       SMB1(74) = SMB0(2)*SMB0(25)
       SMB1(75) = SMB0(3)*SMB0(26)
       SMB1(76) = SMB0(3)*SMB0(27)
       SMB1(77) = SMB0(3)*SMB0(28)
       SMB1(78) = SMB0(3)*SMB0(29)
       SMB1(79) = SMB0(3)*SMB0(30)
       SMB1(80) = SMB0(3)*SMB0(31)
       SMB1(81) = SMB0(3)*SMB0(22)
       SMB1(82) = SMB0(3)*SMB0(24)
       SMB1(83) = SMB0(3)*SMB0(25)
       SMB1(84) = SMB0(4)*SMB0(26)
       SMB1(85) = SMB0(4)*SMB0(27)
       SMB1(86) = SMB0(4)*SMB0(28)
       SMB1(87) = SMB0(4)*SMB0(29)
       SMB1(88) = SMB0(4)*SMB0(30)
       SMB1(89) = SMB0(4)*SMB0(31)
       SMB1(90) = SMB0(4)*SMB0(22)
       SMB1(91) = SMB0(4)*SMB0(23)
       SMB1(92) = SMB0(4)*SMB0(24)
       SMB1(93) = SMB0(4)*SMB0(25)
       SMB1(94) = SMB0(5)*SMB0(26)
       SMB1(95) = SMB0(5)*SMB0(27)
       SMB1(96) = SMB0(5)*SMB0(28)
       SMB1(97) = SMB0(5)*SMB0(29)
       SMB1(98) = SMB0(5)*SMB0(30)
       SMB1(99) = SMB0(5)*SMB0(31)
       SMB1(100) = SMB0(5)*SMB0(22)
       SMB1(101) = SMB0(5)*SMB0(24)
       SMB1(102) = SMB0(5)*SMB0(25)
       SMB1(106) = SMB0(6)*SMB0(26)
       SMB1(107) = SMB0(7)*SMB0(26)
       SMB1(108) = SMB0(8)*SMB0(26)
       SMB1(109) = SMB0(6)*SMB0(27)
       SMB1(110) = SMB0(8)*SMB0(27)
       SMB1(111) = SMB0(6)*SMB0(28)
       SMB1(112) = SMB0(7)*SMB0(28)
       SMB1(113) = SMB0(8)*SMB0(28)
       SMB1(114) = SMB0(9)*SMB0(29)
       SMB1(115) = SMB0(9)*SMB0(30)
       SMB1(116) = SMB0(9)*SMB0(31)
       SMB1(117) = SMB0(12)*SMB0(22)
       SMB1(118) = SMB0(13)*SMB0(22)
       SMB1(119) = SMB0(14)*SMB0(22)
       SMB1(120) = SMB0(12)*SMB0(23)
       SMB1(121) = SMB0(14)*SMB0(23)
       SMB1(122) = SMB0(12)*SMB0(24)
       SMB1(123) = SMB0(13)*SMB0(24)
       SMB1(124) = SMB0(14)*SMB0(24)
       SMB1(125) = SMB0(16)*SMB0(25)
       SMB1(126) = SMB0(18)*SMB0(25)
       SMB1(127) = SMB0(19)*SMB0(25)
c************************************************************************************
       SMB1(17) = dotcc(SMB11,SMB12)
       SMB1(18) = dotcc(SMB12,SMB13)
       SMB1(19) = dotcc(SMB12,SMB14)
       SMB1(20) = dotcc(SMB12,SMB15)
       SMB1(103) = dotcc(SMB11,SMB16)
       SMB1(104) = dotcc(SMB11,SMB17)
       SMB1(105) = dotcc(SMB11,SMB18)
c************************************************************************************
c************************************************************************************
       SMB(1) = SMB1(1)
       SMB(2) = SMB1(9)
       SMB(3) = SMB1(17)
       SMB(4) = 2*mup3mup4*SMB1(17)-SMB1(19)
       SMB(5) = SMB1(52)
       SMB(6) = SMB1(53)
       SMB(7) = SMB1(54)
       SMB(8) = SMB1(55)
       SMB(9) = SMB1(60)
       SMB(10) = SMB1(61)
       SMB(11) = SMB1(68)
       SMB(12) = SMB1(69)
       SMB(13) = SMB1(72)
       SMB(14) = SMB1(73)
       SMB(15) = 0
       SMB(16) = 2*SMB1(55)
       SMB(17) = -4*p2mup4*SMB1(52)+4*p2mup3*SMB1(53)
       SMB(18) = 2*SMB1(61)
       SMB(19) = SMB1(23)
       SMB(20) = SMB1(24)
       SMB(21) = SMB1(27)
       SMB(22) = SMB1(33)
       SMB(23) = SMB1(25)
       SMB(24) = SMB1(26)
       SMB(25) = SMB1(31)
       SMB(26) = SMB1(32)
       SMB(27) = SMB1(38)
       SMB(28) = SMB1(39)
       SMB(29) = SMB1(40)
       SMB(30) = SMB1(41)
       SMB(31) = SMB1(21)
       SMB(32) = SMB1(45)
       SMB(33) = SMB1(46)
       SMB(34) = SMB1(22)
       SMB(35) = SMB1(56)
       SMB(36) = SMB1(62)
       SMB(37) = SMB1(70)
       SMB(38) = SMB1(74)
       SMB(39) = SMB1(83)
       SMB(40) = SMB1(93)
       SMB(41) = SMB1(102)
       SMB(42) = SMB1(126)
       SMB(43) = SMB1(127)
       SMB(44) = SMB1(125)
       SMB(45) = SMB1(4)
       SMB(46) = SMB1(6)
       SMB(47) = SMB1(7)
       SMB(48) = SMB1(12)
       SMB(49) = SMB1(14)
       SMB(50) = SMB1(15)
       SMB(51) = 0
       SMB(52) = 2*p4mup3*SMB1(17)-SMB1(18)
       SMB(53) = -2*SMB1(4)+2*p6mup3*SMB1(17)
       SMB(54) = 0
       SMB(55) = 2*p4mup4*SMB1(17)-SMB1(20)
       SMB(56) = -2*SMB1(12)+2*p6mup4*SMB1(17)
       SMB(57) = 0
       SMB(58) = 2*SMB1(21)
       SMB(59) = 0
       SMB(60) = 2*SMB1(22)
       SMB(61) = 0
       SMB(62) = SMB1(18)
       SMB(63) = 2*p5mup3*SMB1(17)-2*SMB1(52)
       SMB(64) = 0
       SMB(65) = SMB1(20)
       SMB(66) = 2*p5mup4*SMB1(17)-2*SMB1(53)
       SMB(67) = SMB1(19)
       SMB(68) = 4*mup3mup4*SMB1(17)
       SMB(69) = -4*p1mup3*p2mup4*SMB1(17)+s12*SMB1(103)
       SMB(70) = -4*p2mup4*p4mup3*SMB1(17)+(p3sq-s23+s234-s34)*SMB1(10
     &   3)+2*p2mup4*SMB1(104)
       SMB(71) = -2*(2*p2mup4*SMB1(4)-2*mup3mup4*SMB1(27)+SMB1(97))
       SMB(72) = (-s123+s23-s234+s56)*SMB1(103)-2*p1mup3*SMB1(105)
       SMB(73) = p4sq*SMB1(103)
       SMB(74) = -4*p4mup4*SMB1(4)+2*p6mup4*SMB1(18)+4*mup3mup4*SMB1(2
     &   1)-2*SMB1(98)
       SMB(75) = -4*p1mup3*SMB1(53)-(s16-s234+s56)*SMB1(103)
       SMB(76) = 4*(p4mup4*p5mup3-p4mup3*p5mup4)*SMB1(17)-4*p4mup4*SMB
     &   1(52)+4*mup3mup4*SMB1(55)-2*SMB1(89)+2*p5mup4*SMB1(104)
       SMB(77) = -2*(2*p5mup4*SMB1(4)-2*p5mup3*p6mup4*SMB1(17)+2*p6mup
     &   4*SMB1(52)-2*mup3mup4*SMB1(56)+SMB1(99))
       SMB(78) = 4*mup3mup4*SMB1(17)
       SMB(79) = (-4*p1mup4*p2mup3+2*mup3mup4*s12)*SMB1(17)-s12*SMB1(1
     &   03)
       SMB(80) = (-4*p2mup3*p4mup4+2*mup3mup4*(p3sq-s23+s234-s34))*SMB
     &   1(17)+(-p3sq+s23-s234+s34)*SMB1(103)+2*p2mup3*SMB1(105)
       SMB(81) = -4*p2mup3*SMB1(12)+2*SMB1(97)
       SMB(82) = -2*mup3mup4*(s123-s23+s234-s56)*SMB1(17)+(s123-s23+s2
     &   34-s56)*SMB1(103)-2*p1mup4*SMB1(104)
       SMB(83) = p4sq*(2*mup3mup4*SMB1(17)-SMB1(103))
       SMB(84) = 2*(-2*p4mup3*SMB1(12)+p6mup3*SMB1(20)+SMB1(98))
       SMB(85) = -2*mup3mup4*(s16-s234+s56)*SMB1(17)-4*p1mup4*SMB1(52)
     &   +(s16-s234+s56)*SMB1(103)
       SMB(86) = (-4*p4mup4*p5mup3+4*p4mup3*p5mup4)*SMB1(17)-4*p4mup3*
     &   SMB1(53)+2*SMB1(89)+2*p5mup3*SMB1(105)
       SMB(87) = -4*p5mup3*SMB1(12)+4*p5mup4*p6mup3*SMB1(17)-4*p6mup3*
     &   SMB1(53)+2*SMB1(99)
       SMB(88) = 16*SMB1(17)
       SMB(89) = 32*mup3mup4*SMB1(17)-16*SMB1(19)
       SMB(90) = 4*p2mup3*SMB1(17)
       SMB(91) = 4*p2mup4*SMB1(17)
       SMB(92) = 2*s12*SMB1(17)
       SMB(93) = 0
       SMB(94) = 2*(p3sq-s23+s234-s34)*SMB1(17)
       SMB(95) = -2*(s12+s16-s345)*SMB1(17)
       SMB(96) = 4*mup3mup4*s12*SMB1(17)-2*s12*SMB1(19)
       SMB(97) = (-4*p2mup4*p4mup3+4*p2mup3*p4mup4)*SMB1(17)+2*p2mup4*
     &   SMB1(18)-2*p2mup3*SMB1(20)
       SMB(98) = 2*(p3sq-s23+s234-s34)*(2*mup3mup4*SMB1(17)-SMB1(19))
       SMB(99) = 2*(s12+s16-s345)*(-2*mup3mup4*SMB1(17)+SMB1(19))
       SMB(100) = 2*(SMB1(18)+SMB1(104))
       SMB(101) = 2*(SMB1(20)+SMB1(105))
       SMB(102) = 0
       SMB(103) = 2*(-s123+s23-s234+s56)*SMB1(17)
       SMB(104) = 4*p4sq*SMB1(17)
       SMB(105) = 2*(s123-s45-s56)*SMB1(17)
       SMB(106) = 0
       SMB(107) = 2*(s123-s23+s234-s56)*(-2*mup3mup4*SMB1(17)+SMB1(19)
     &   )
       SMB(108) = 8*mup3mup4*p4sq*SMB1(17)-4*p4sq*SMB1(19)
       SMB(109) = 2*(s123-s45-s56)*(2*mup3mup4*SMB1(17)-SMB1(19))
       SMB(110) = 4*p5mup3*SMB1(17)
       SMB(111) = 4*p5mup4*SMB1(17)
       SMB(112) = -2*(s16-s234+s56)*SMB1(17)
       SMB(113) = 2*(-p4sq+s45)*SMB1(17)
       SMB(114) = 2*s56*SMB1(17)
       SMB(115) = 2*(s16-s234+s56)*(-2*mup3mup4*SMB1(17)+SMB1(19))
       SMB(116) = 2*(p4sq-s45)*(-2*mup3mup4*SMB1(17)+SMB1(19))
       SMB(117) = 4*mup3mup4*s56*SMB1(17)-2*s56*SMB1(19)
       SMB(118) = SMB1(3)
       SMB(119) = SMB1(11)
       SMB(120) = SMB1(13)
       SMB(121) = SMB1(67)
       SMB(122) = SMB1(71)
       SMB(123) = SMB1(87)
       SMB(124) = SMB1(88)
       SMB(125) = SMB1(89)
       SMB(126) = SMB1(2)
       SMB(127) = SMB1(5)
       SMB(128) = SMB1(10)
       SMB(129) = SMB1(57)
       SMB(130) = SMB1(58)
       SMB(131) = SMB1(59)
       SMB(132) = SMB1(63)
       SMB(133) = SMB1(65)
       SMB(134) = SMB1(66)
       SMB(135) = SMB1(76)
       SMB(136) = SMB1(77)
       SMB(137) = SMB1(80)
       SMB(138) = SMB1(82)
       SMB(139) = SMB1(85)
       SMB(140) = SMB1(86)
       SMB(141) = SMB1(95)
       SMB(142) = SMB1(96)
       SMB(143) = SMB1(109)
       SMB(144) = SMB1(111)
       SMB(145) = SMB1(112)
       SMB(146) = SMB1(120)
       SMB(147) = SMB1(122)
       SMB(148) = SMB1(123)
       SMB(149) = 2*s12*SMB1(52)-4*p2mup3*SMB1(54)
       SMB(150) = 2*(p3sq-s23+s234-s34)*SMB1(52)-4*p2mup3*SMB1(55)
       SMB(151) = -2*(s12+s16-s345)*SMB1(52)-4*p2mup3*SMB1(56)
       SMB(152) = (-p3sq+s23-s234+s34)*(-s16+s234-s56)*SMB1(17)+2*s12*
     &   SMB1(55)
       SMB(153) = 2*mup3mup4*(p3sq-s23+s234-s34)*(s16-s234+s56)*SMB1(1
     &   7)-(p3sq-s23+s234-s34)*(s16-s234+s56)*SMB1(19)+2*s12*SMB1(61
     &   )+2*p2mup4*(s16-s234+s56)*SMB1(104)-2*p2mup3*(s16-s234+s56)*
     &   SMB1(105)
       SMB(154) = -2*(s123-s23+s234-s56)*SMB1(52)-2*SMB1(57)
       SMB(155) = 2*p4sq*SMB1(52)
       SMB(156) = 2*SMB1(59)
       SMB(157) = -2*(s123-s23+s234-s56)*SMB1(53)-2*SMB1(63)
       SMB(158) = p4sq*(s16-s234+s56)*SMB1(17)+2*(-s123+s23-s234+s56)*
     &   SMB1(55)
       SMB(159) = 2*(s16-s234+s56)*SMB1(21)+2*(-s123+s23-s234+s56)*SMB
     &   1(56)
       SMB(160) = 2*mup3mup4*p4sq*(s16-s234+s56)*SMB1(17)-p4sq*(s16-s2
     &   34+s56)*SMB1(19)+2*(-s123+s23-s234+s56)*SMB1(61)
       SMB(161) = 2*(s16-s234+s56)*SMB1(22)+2*(-s123+s23-s234+s56)*SMB
     &   1(62)
       SMB(162) = SMB1(64)
       SMB(163) = SMB1(91)
       SMB(164) = SMB1(92)
       SMB(165) = SMB1(99)
       SMB(166) = SMB1(101)
       SMB(167) = SMB1(110)
       SMB(168) = SMB1(113)
       SMB(169) = SMB1(121)
       SMB(170) = SMB1(124)
       SMB(171) = 2*(p3sq-s23+s234-s34)*SMB1(53)-4*p2mup4*SMB1(55)
       SMB(172) = 2*(p3sq-s23+s234-s34)*SMB1(56)-2*SMB1(102)
       SMB(173) = 2*(p3sq-s23+s234-s34)*SMB1(62)-2*SMB1(125)
       SMB(174) = 2*p4sq*SMB1(53)
       SMB(175) = 2*SMB1(64)
       SMB(176) = 2*p4sq*SMB1(56)
       SMB(177) = 2*p4sq*SMB1(62)
       SMB(178) = SMB1(78)
       SMB(179) = SMB1(97)
       SMB(180) = -2*p4mup3*(s12+s16-s345)*SMB1(17)+(s12+s16-s345)*SMB
     &   1(18)-4*p2mup3*SMB1(21)
       SMB(181) = -2*p4mup4*(s12+s16-s345)*SMB1(17)+(s12+s16-s345)*SMB
     &   1(20)-4*p2mup4*SMB1(21)
       SMB(182) = 2*s12*SMB1(21)+2*(s123-s23+s234-s56)*SMB1(27)
       SMB(183) = 2*(p3sq-s23+s234-s34)*SMB1(21)-2*p4sq*SMB1(27)
       SMB(184) = 2*s12*SMB1(22)+2*(s123-s23+s234-s56)*SMB1(33)
       SMB(185) = 2*(p3sq-s23+s234-s34)*SMB1(22)-2*p4sq*SMB1(33)
       SMB(186) = SMB1(28)
       SMB(187) = SMB1(29)
       SMB(188) = SMB1(30)
       SMB(189) = SMB1(35)
       SMB(190) = SMB1(36)
       SMB(191) = SMB1(37)
       SMB(192) = SMB1(42)
       SMB(193) = SMB1(44)
       SMB(194) = SMB1(90)
       SMB(195) = SMB1(117)
       SMB(196) = SMB1(119)
       SMB(197) = (2*p4mup3*s12+2*p2mup3*(s123-s23+s234-s56))*SMB1(17)
     &   -s12*SMB1(18)
       SMB(198) = (-2*p2mup3*p4sq+2*p4mup3*(p3sq-s23+s234-s34))*SMB1(1
     &   7)+(-p3sq+s23-s234+s34)*SMB1(18)
       SMB(199) = (-2*p2mup4*p4sq+2*p4mup4*(p3sq-s23+s234-s34))*SMB1(1
     &   7)+(-p3sq+s23-s234+s34)*SMB1(20)
       SMB(200) = (p4sq*s12+(p3sq-s23+s234-s34)*(s123-s23+s234-s56))*S
     &   MB1(17)
       SMB(201) = 2*mup3mup4*(p4sq*s12+(p3sq-s23+s234-s34)*(s123-s23+s
     &   234-s56))*SMB1(17)+(-(p4sq*s12)+(p3sq-s23+s234-s34)*(-s123+s
     &   23-s234+s56))*SMB1(19)+2*p2mup4*(s123-s23+s234-s56)*SMB1(104
     &   )+2*p2mup3*(-s123+s23-s234+s56)*SMB1(105)
       SMB(202) = (2*p4mup4*s12+2*p2mup4*(s123-s23+s234-s56))*SMB1(17)
     &   -s12*SMB1(20)
       SMB(203) = SMB1(43)
       SMB(204) = SMB1(48)
       SMB(205) = SMB1(49)
       SMB(206) = SMB1(50)
       SMB(207) = SMB1(51)
       SMB(208) = SMB1(34)
       SMB(209) = 0
       SMB(210) = 0
       SMB(211) = 4*p6mup3*SMB1(17)
       SMB(212) = 2*s12*SMB1(4)-4*p1mup3*SMB1(27)
       SMB(213) = 2*(p3sq-s23+s234-s34)*SMB1(4)-2*SMB1(94)
       SMB(214) = -4*p4mup4*s12*SMB1(4)+4*p1mup3*(p3sq-s23+s234-s34)*S
     &   MB1(12)-2*(-2*p4mup3*p6mup4*s12+2*p1mup3*p6mup4*(p3sq-s23+s2
     &   34-s34)+2*p1mup3*p4mup4*(s12+s16-s345))*SMB1(17)+2*p1mup3*(s
     &   12+s16-s345)*SMB1(20)+(-8*p1mup3*p2mup4+4*mup3mup4*s12)*SMB1
     &   (21)-2*s12*SMB1(98)-2*p6mup4*s12*SMB1(104)
       SMB(215) = -2*(s123-s23+s234-s56)*SMB1(4)-4*p1mup3*SMB1(21)
       SMB(216) = 2*p4sq*SMB1(4)
       SMB(217) = 2*(2*p4mup4*(s123-s23+s234-s56)*SMB1(4)+2*p1mup3*p4s
     &   q*(SMB1(12)-p6mup4*SMB1(17))-(s123-s23+s234-s56)*(p6mup4*SMB
     &   1(18)+2*mup3mup4*SMB1(21)-SMB1(98)))
       SMB(218) = -2*(s16-s234+s56)*SMB1(4)-4*p1mup3*SMB1(56)
       SMB(219) = 2*SMB1(96)
       SMB(220) = 4*p4mup4*(s16-s234+s56)*SMB1(4)-2*p6mup4*(s16-s234+s
     &   56)*SMB1(18)-4*mup3mup4*(s16-s234+s56)*SMB1(21)-4*p1mup3*SMB
     &   1(64)+2*(s16-s234+s56)*SMB1(98)
       SMB(221) = -4*p1mup4*SMB1(4)+4*p1mup3*SMB1(12)
       SMB(222) = -8*p1mup4*(SMB1(4)-p6mup3*SMB1(17))
       SMB(223) = 2*(2*p1mup4*(p3sq-s23+s234-s34)*SMB1(4)-2*p4mup3*s12
     &   *SMB1(12)+(2*p4mup4*p6mup3*s12-2*p1mup4*p6mup3*(p3sq-s23+s23
     &   4-s34)-2*p1mup4*p4mup3*(s12+s16-s345))*SMB1(17)+p1mup4*(s12+
     &   s16-s345)*SMB1(18)-4*p1mup4*p2mup3*SMB1(21)+s12*SMB1(98)-p6m
     &   up3*s12*SMB1(105))
       SMB(224) = 4*p1mup4*p4sq*(SMB1(4)-p6mup3*SMB1(17))+2*(s123-s23+
     &   s234-s56)*(2*p4mup3*SMB1(12)-p6mup3*SMB1(20)-SMB1(98))
       SMB(225) = 4*p4mup3*(s16-s234+s56)*SMB1(12)-2*p6mup3*(s16-s234+
     &   s56)*SMB1(20)-4*p1mup4*SMB1(59)-2*(s16-s234+s56)*SMB1(98)
       SMB(226) = 32*(-SMB1(4)+p6mup3*SMB1(17))
       SMB(227) = -4*s12*(SMB1(4)-p6mup3*SMB1(17))
       SMB(228) = -4*(p3sq-s23+s234-s34)*(SMB1(4)-p6mup3*SMB1(17))
       SMB(229) = 4*s12*SMB1(21)
       SMB(230) = 4*s12*SMB1(22)
       SMB(231) = 4*(s123-s23+s234-s56)*(SMB1(4)-p6mup3*SMB1(17))
       SMB(232) = -8*p4sq*(SMB1(4)-p6mup3*SMB1(17))
       SMB(233) = 4*(-s123+s23-s234+s56)*SMB1(21)
       SMB(234) = 4*(-s123+s23-s234+s56)*SMB1(22)
       SMB(235) = 4*(s16-s234+s56)*(SMB1(4)-p6mup3*SMB1(17))
       SMB(236) = 4*(p4sq-s45)*(SMB1(4)-p6mup3*SMB1(17))
       SMB(237) = -4*(s16-s234+s56)*SMB1(21)
       SMB(238) = -4*(s16-s234+s56)*SMB1(22)
       SMB(239) = SMB1(8)
       SMB(240) = SMB1(16)
       SMB(241) = SMB1(47)
       SMB(242) = SMB1(75)
       SMB(243) = SMB1(79)
       SMB(244) = SMB1(81)
       SMB(245) = SMB1(84)
       SMB(246) = SMB1(94)
       SMB(247) = SMB1(98)
       SMB(248) = 4*p1mup3*SMB1(17)
       SMB(249) = 8*p4mup3*SMB1(17)-2*(SMB1(18)+SMB1(104))
       SMB(250) = p1mup3*(8*p4mup4*SMB1(17)-4*SMB1(20))
       SMB(251) = (2*p4mup3*s12-2*p1mup3*(p3sq-s23+s234-s34))*SMB1(17)
     &   -s12*SMB1(104)
       SMB(252) = -2*p1mup3*p4sq*SMB1(17)-(s123-s23+s234-s56)*(2*p4mup
     &   3*SMB1(17)-SMB1(104))
       SMB(253) = SMB1(106)
       SMB(254) = SMB1(107)
       SMB(255) = -2*p4mup3*(s16-s234+s56)*SMB1(17)-4*p1mup3*SMB1(55)+
     &   (s16-s234+s56)*SMB1(104)
       SMB(256) = (-4*p1mup4*p4mup3+4*p1mup3*p4mup4)*SMB1(17)+2*p1mup4
     &   *SMB1(104)-2*p1mup3*SMB1(105)
       SMB(257) = SMB1(114)
       SMB(258) = SMB1(115)
       SMB(259) = SMB1(116)
       SMB(260) = p1mup4*(8*p4mup3*SMB1(17)-4*SMB1(18))
       SMB(261) = (2*p4mup4*s12-2*p1mup4*(p3sq-s23+s234-s34))*SMB1(17)
     &   -s12*SMB1(105)
       SMB(262) = -2*p1mup4*p4sq*SMB1(17)-(s123-s23+s234-s56)*(2*p4mup
     &   4*SMB1(17)-SMB1(105))
       SMB(263) = SMB1(118)
       SMB(264) = -2*p4mup4*(s16-s234+s56)*SMB1(17)-4*p1mup4*SMB1(55)+
     &   (s16-s234+s56)*SMB1(105)
       SMB(265) = 0
       SMB(266) = 32*p4mup3*SMB1(17)-16*SMB1(18)
       SMB(267) = 0
       SMB(268) = 0
       SMB(269) = 4*p4mup3*s12*SMB1(17)-2*s12*SMB1(18)
       SMB(270) = 4*p4mup4*s12*SMB1(17)-2*s12*SMB1(20)
       SMB(271) = 2*(s123-s23+s234-s56)*(-2*p4mup3*SMB1(17)+SMB1(18))
       SMB(272) = 2*(s123-s23+s234-s56)*(-2*p4mup4*SMB1(17)+SMB1(20))
       SMB(273) = 2*(s16-s234+s56)*(-2*p4mup3*SMB1(17)+SMB1(18))
       SMB(274) = 2*(s16-s234+s56)*(-2*p4mup4*SMB1(17)+SMB1(20))
       SMB(275) = -8*p4mup3*SMB1(12)+4*p6mup4*(SMB1(18)+SMB1(104))
       SMB(276) = 2*SMB1(98)
       SMB(277) = 4*p6mup4*SMB1(17)
       SMB(278) = -8*p4mup4*SMB1(4)+4*p6mup3*(SMB1(20)+SMB1(105))
       SMB(279) = 2*(p3sq-s23+s234-s34)*SMB1(12)-2*SMB1(100)
       SMB(280) = 2*p4sq*SMB1(12)
       SMB(281) = 2*SMB1(101)
       SMB(282) = 32*SMB1(21)
       SMB(283) = 32*SMB1(22)
       SMB(284) = -4*(p3sq-s23+s234-s34)*(SMB1(12)-p6mup4*SMB1(17))
       SMB(285) = -8*p4sq*(SMB1(12)-p6mup4*SMB1(17))
       SMB(286) = 4*(p4sq-s45)*(SMB1(12)-p6mup4*SMB1(17))
       SMB(287) = SMB1(100)
       SMB(288) = SMB1(108)
       SMB(289) = 8*p4mup4*SMB1(17)-2*(SMB1(20)+SMB1(105))
       SMB(290) = 32*p4mup4*SMB1(17)-16*SMB1(20)
       SMB(291) = 4*p1mup4*SMB1(17)
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
       result(1) = (mup3mup4*(F(1)*SMB(3)+8*F(4)*SMB(7)+F(5)*SMB(8)+4*
     &   (F(7)*SMB(16)+s234*F(10)*SMB(57)+F(9)*SMB(58)))+F(13)*SMB(11
     &   5)+F(14)*SMB(116)+s234*F(15)*SMB(117)+4*(F(6)*SMB(9)+F(3)*SM
     &   B(10)+Fa(1)*(2*F(3)*SMB(5)-s56*F(8)*SMB(51)+F(9)*SMB(53))+F(
     &   16)*SMB(134))+2*(F(2)*SMB(4)-s234*F(10)*SMB(59)-F(9)*SMB(60)
     &   -s56*F(8)*SMB(106)-F(7)*(SMB(18)+Fa(1)*(-2*SMB(63)+SMB(110))
     &   )-4*F(16)*Fa(1)*SMB(130)+F(17)*Fa(1)*SMB(235)+mup3mup4*(2*s5
     &   6*F(8)*SMB(102)+F(11)*SMB(112)+F(7)*SMB(113)+s234*F(12)*SMB(
     &   114)-4*F(16)*SMB(133)+F(17)*SMB(237)))+F(18)*SMB(238)+F(1307
     &   )*(P(35)*SMB(89)-SMB(268)-2*(Fa(1)*(SMB(226)+SMB(265))-mup3m
     &   up4*(P(30)*SMB(88)+SMB(267)-SMB(282)))+SMB(283)))/(s23*s234)
       result(2) = (-(F(26)*SMB(10))+F(32)*SMB(34)+F(34)*SMB(36)+F(35)
     &   *SMB(49)+Fa(5)*SMB(52)+F(51)*SMB(106)+F(52)*SMB(107)+F(50)*S
     &   MB(109)+s45*F(54)*SMB(115)-F(53)*(SMB(108)+SMB(116))+F(49)*(
     &   4*SMB(87)+SMB(117))+4*(F(19)*SMB(2)+F(36)*SMB(50)-Fa(4)*SMB(
     &   51)+Fa(6)*SMB(53)+Fa(7)*SMB(62)+Fa(8)*SMB(63)-F(45)*SMB(82)+
     &   F(48)*SMB(84)-s45*F(23)*(SMB(9)+SMB(85))+p1mup4*(-2*F(55)*SM
     &   B(131)+F(57)*SMB(156)+4*F(58)*SMB(193)))-8*(Fa(2)*SMB(5)-F(2
     &   7)*SMB(12)+Fa(3)*SMB(27)+F(57)*SMB(148)+p6mup4*(F(55)*SMB(12
     &   9)+F(57)*SMB(154)-2*F(58)*SMB(192))+Fa(11)*(F(56)*SMB(130)+2
     &   *F(58)*SMB(203)))+mup3mup4*(8*(F(25)*SMB(8)+F(28)*SMB(29)+F(
     &   30)*SMB(30)+F(31)*SMB(31)+F(33)*SMB(35))+4*(F(20)*SMB(3)+F(4
     &   4)*SMB(57)-F(41)*SMB(58))+F(40)*SMB(102)+2*(F(45)*SMB(103)+F
     &   (46)*SMB(105)+s45*F(23)*(4*SMB(7)+SMB(112))+F(27)*(-2*SMB(16
     &   )+SMB(104)+SMB(113))+F(47)*SMB(114)+4*F(55)*SMB(133)-8*F(58)
     &   *SMB(206)+F(57)*(4*SMB(159)+SMB(233)+SMB(237))))+F(59)*(SMB(
     &   234)+SMB(238))+2*(F(21)*SMB(4)+F(39)*SMB(59)+F(41)*SMB(60)-F
     &   (27)*(-SMB(18)+2*(SMB(83)+SMB(86)))-Fa(9)*SMB(100)-Fa(10)*SM
     &   B(110)-2*(F(30)*SMB(33)+F(28)*(SMB(32)+SMB(120))+F(55)*SMB(1
     &   34))+4*F(58)*SMB(207)-F(57)*(2*(SMB(161)+SMB(224)+SMB(225))-
     &   Fa(11)*(2*SMB(209)+SMB(231)+SMB(235))+p1mup4*(SMB(232)+SMB(2
     &   36))+p6mup4*(SMB(271)+SMB(273))))+F(1306)*(P(90)*(12*SMB(67)
     &   +6*SMB(78)-2*mup3mup4*SMB(88)+SMB(89))-SMB(268)-6*(SMB(260)-
     &   SMB(278))-2*(Fa(11)*(SMB(226)+SMB(265))-Fa(12)*SMB(266)-mup3
     &   mup4*SMB(267)+mup3mup4*SMB(282))+SMB(283)))/(s23*s45)
       result(3) = (Fa(13)*SMB(3)+F(74)*SMB(26)+Fa(29)*SMB(57)+Fa(30)*
     &   SMB(58)+F(84)*SMB(99)+F(66)*(-4*(SMB(70)-SMB(80))+2*SMB(97)+
     &   SMB(98)+8*(SMB(121)-SMB(122)-SMB(123)))+4*(F(60)*SMB(1)+F(61
     &   )*SMB(2)+F(72)*P(50)*SMB(22)-F(79)*SMB(34)+F(76)*(-SMB(47)+S
     &   MB(50))-Fa(24)*SMB(51)+Fa(25)*SMB(52)+Fa(26)*SMB(53)+Fa(27)*
     &   SMB(54)-Fa(28)*SMB(55)+F(67)*P(50)*SMB(71)+F(65)*P(50)*SMB(8
     &   1)+F(65)*Fa(16)*SMB(92)+Fa(36)*SMB(102)-F(75)*(SMB(32)+SMB(1
     &   20)-SMB(127))+p1mup4*F(65)*SMB(180)-p1mup3*F(65)*SMB(181)-F(
     &   65)*SMB(184)-p2mup4*F(65)*SMB(209)+p2mup3*F(65)*SMB(210)-F(6
     &   5)*SMB(214)+F(65)*SMB(223))+F(65)*SMB(230)+8*(Fa(14)*SMB(19)
     &   -F(71)*Fa(15)*SMB(21)-Fa(17)*SMB(24)+Fa(20)*SMB(29)+Fa(21)*S
     &   MB(31)-Fa(22)*SMB(45)+F(67)*P(116)*SMB(46)+Fa(18)*(-SMB(27)+
     &   SMB(118))+Fa(19)*(SMB(28)-SMB(119))-Fa(37)*SMB(126)+F(67)*P(
     &   50)*SMB(179)+mup3mup4*F(65)*SMB(182)+F(70)*SMB(191)-p6mup4*F
     &   (65)*SMB(197)-F(65)*Fa(11)*SMB(212)+p1mup4*F(65)*SMB(213)+p6
     &   mup4*F(65)*SMB(251)+F(65)*SMB(254)-F(65)*SMB(257)-F(65)*SMB(
     &   263)-p1mup3*F(65)*SMB(279))+2*(F(68)*SMB(4)+F(65)*(2*Fa(23)*
     &   (4*SMB(48)+SMB(56))+P(116)*(4*SMB(49)+SMB(59)))+F(82)*SMB(60
     &   )+Fa(31)*(2*SMB(61)+SMB(90))-Fa(33)*(2*SMB(93)+SMB(94))+Fa(3
     &   4)*SMB(95)-F(81)*SMB(106)-F(65)*Fa(11)*SMB(227)+p1mup4*F(65)
     &   *SMB(228)-mup3mup4*F(65)*SMB(229)+p6mup4*F(65)*SMB(269)-p6mu
     &   p3*(-(F(67)*P(50)*(2*SMB(64)+SMB(91)))+8*F(70)*(P(50)*SMB(20
     &   )+SMB(208))+F(65)*(8*(P(116)*SMB(128)-SMB(244))-4*(SMB(202)-
     &   SMB(261))+SMB(270)))-p1mup3*F(65)*SMB(284))+16*(p6mup4*F(70)
     &   *SMB(186)+F(73)*(-2*Fa(16)*SMB(23)+Fa(11)*SMB(187))+p1mup4*F
     &   (70)*SMB(188)-p1mup3*F(70)*SMB(189)+mup3mup4*F(73)*SMB(190)-
     &   p2mup4*F(65)*SMB(203)+p2mup3*F(65)*SMB(204)-p2mup4*F(65)*SMB
     &   (239)+p2mup3*F(65)*SMB(240)-p6mup4*F(65)*SMB(242)+p1mup4*F(6
     &   5)*SMB(246)-p1mup3*F(65)*SMB(287))+F(1305)*(12*(SMB(74)-SMB(
     &   84))+P(45)*SMB(89)+48*(SMB(243)+SMB(247))-SMB(268)+24*(SMB(7
     &   2)-SMB(82)+SMB(256)-SMB(276))+6*(P(50)*(SMB(68)-SMB(78))-SMB
     &   (107)-SMB(109)+SMB(250)-SMB(260)-SMB(275)+SMB(278))+SMB(283)
     &   +2*(Fa(32)*SMB(88)-Fa(11)*(SMB(226)+6*(SMB(211)+SMB(248))+SM
     &   B(265))+Fa(12)*(-6*(SMB(100)-SMB(249))+SMB(266))+mup3mup4*(6
     &   *(SMB(103)+SMB(105))+SMB(267)-SMB(282))+Fa(35)*(6*(SMB(101)-
     &   SMB(289))-SMB(290)))))/(s345*s45)
       result(4) = (F(87)*SMB(4)+Fa(38)*SMB(5)+F(94)*SMB(9)+F(95)*SMB(
     &   10)+F(107)*SMB(49)+Fa(40)*SMB(51)+Fa(41)*SMB(52)+Fa(42)*SMB(
     &   53)+F(123)*SMB(59)+F(124)*SMB(60)+Fa(43)*SMB(62)+F(135)*SMB(
     &   83)+Fa(45)*SMB(100)+F(148)*SMB(106)+F(149)*SMB(107)+F(150)*S
     &   MB(108)+F(151)*SMB(109)+F(159)*SMB(115)+F(160)*SMB(116)+F(16
     &   1)*SMB(117)+Fa(51)*SMB(155)+8*(Fa(39)*SMB(27)+Fa(47)*SMB(129
     &   )-Fa(48)*SMB(130)+Fa(49)*SMB(131)-Fa(53)*SMB(193)-Fa(54)*SMB
     &   (203))+4*(F(85)*SMB(2)-F(96)*SMB(12)+F(104)*SMB(32)+F(99)*SM
     &   B(33)+F(103)*SMB(34)+F(106)*SMB(36)-F(108)*SMB(50)+Fa(44)*SM
     &   B(63)+F(162)*SMB(120)+F(175)*SMB(134)+F(176)*SMB(147)+Fa(50)
     &   *SMB(154)+Fa(52)*SMB(156)-F(191)*SMB(170)+F(195)*SMB(207)+Fa
     &   (55)*SMB(209))+F(216)*SMB(234)+mup3mup4*(F(86)*SMB(3)+F(92)*
     &   SMB(7)+F(93)*SMB(8)+4*F(97)*SMB(16)+8*(F(102)*SMB(29)-F(99)*
     &   SMB(30)-F(103)*SMB(31)+F(105)*SMB(35))+F(121)*SMB(57)+F(122)
     &   *SMB(58)+F(144)*SMB(102)+F(146)*SMB(104)-2*(F(145)*SMB(103)+
     &   F(147)*SMB(105)+F(156)*SMB(112)+F(157)*SMB(113)+F(158)*SMB(1
     &   14)+2*(F(189)*SMB(158)-F(192)*SMB(176))+4*(F(175)*SMB(133)+F
     &   (195)*SMB(206))-F(193)*SMB(233)-F(190)*(2*SMB(159)+SMB(237))
     &   ))+F(225)*SMB(238)-2*(F(97)*SMB(18)+F(133)*(2*SMB(67)+SMB(78
     &   ))-F(134)*SMB(82)-F(136)*SMB(84)-F(137)*SMB(85)-F(138)*SMB(8
     &   6)-F(139)*SMB(87)-Fa(46)*SMB(110)-F(189)*SMB(160)+F(192)*SMB
     &   (177)-F(203)*SMB(222)+F(193)*SMB(224)+F(190)*(SMB(161)+SMB(2
     &   25))-Fa(57)*SMB(231)+Fa(58)*SMB(232)-Fa(59)*SMB(235)+Fa(60)*
     &   SMB(236)+F(226)*SMB(260)+Fa(56)*(SMB(226)+SMB(265))+Fa(61)*S
     &   MB(266)-Fa(62)*SMB(271)-Fa(63)*SMB(273)-F(238)*SMB(278))+F(1
     &   309)*(P(90)*SMB(89)-SMB(268)-2*mup3mup4*(P(90)*SMB(88)-SMB(2
     &   67)+SMB(282))+SMB(283)))/s23
       result(5) = (F(240)*SMB(2)+Fa(64)*SMB(3)-Fa(65)*SMB(5)+Fa(66)*S
     &   MB(6)+F(285)*SMB(9)+F(286)*SMB(10)+Fa(69)*SMB(15)+Fa(70)*SMB
     &   (16)+Fa(71)*SMB(19)+Fa(73)*SMB(21)+F(317)*SMB(26)+Fa(77)*SMB
     &   (29)+Fa(83)*SMB(45)+F(347)*SMB(46)+Fa(84)*SMB(48)+F(353)*SMB
     &   (49)+Fa(85)*SMB(51)+Fa(86)*SMB(52)+Fa(87)*SMB(53)+Fa(88)*SMB
     &   (54)+Fa(89)*SMB(55)+Fa(90)*SMB(56)-Fa(91)*SMB(57)+Fa(92)*SMB
     &   (58)+F(381)*SMB(59)+Fa(95)*(SMB(65)+SMB(66))+F(393)*SMB(70)+
     &   F(398)*SMB(80)+Fa(97)*SMB(90)+Fa(98)*SMB(92)+Fa(99)*SMB(93)-
     &   Fa(100)*SMB(94)+Fa(101)*SMB(95)-F(408)*SMB(97)+F(409)*SMB(98
     &   )+F(410)*SMB(99)+Fa(103)*SMB(101)+Fa(104)*SMB(102)+Fa(107)*S
     &   MB(105)+F(419)*SMB(106)+F(421)*SMB(109)+Fa(109)*SMB(111)+Fa(
     &   111)*SMB(114)+F(420)*(SMB(108)+SMB(116))+F(424)*SMB(117)+Fa(
     &   113)*SMB(119)+F(433)*SMB(122)+F(434)*SMB(123)+F(289)*(4*SMB(
     &   14)-2*(SMB(73)+SMB(76)+2*(SMB(124)+SMB(125))))+Fa(114)*SMB(1
     &   26)+Fa(115)*SMB(128)+mup3mup4*(s45*F(307)*SMB(23)-8*(F(267)*
     &   SMB(39)+F(343)*SMB(41)-F(335)*SMB(133)))+s45*(F(399)*(4*SMB(
     &   82)+SMB(107))+F(423)*SMB(115)+2*F(396)*(SMB(75)-2*(SMB(85)-S
     &   MB(137))))+Fa(126)*SMB(181)+F(445)*SMB(230)+F(451)*(SMB(234)
     &   +SMB(238))+F(288)*(s45*(-2*SMB(69)+SMB(96)+4*(SMB(79)-SMB(17
     &   8)))+4*(P(50)*SMB(13)-2*Fa(11)*SMB(135)+SMB(143)+SMB(253)))+
     &   4*(F(239)*SMB(1)+F(306)*SMB(22)+F(318)*SMB(32)+F(321)*SMB(33
     &   )+F(341)*SMB(34)-F(342)*SMB(36)+F(343)*SMB(44)+F(348)*SMB(47
     &   )+F(354)*SMB(50)+F(391)*SMB(67)+F(400)*SMB(84)+F(432)*SMB(12
     &   0)-F(319)*SMB(121)-F(427)*SMB(127)-F(441)*SMB(134)+F(267)*(p
     &   4sq*SMB(25)+SMB(42)-SMB(146))+F(443)*(SMB(145)-2*(SMB(148)+p
     &   6mup4*(SMB(136)+SMB(154)))+p1mup4*(2*SMB(142)+SMB(156)))-Fa(
     &   124)*SMB(175)-Fa(125)*SMB(180)+F(446)*SMB(191)-Fa(132)*SMB(1
     &   97)+Fa(81)*(2*SMB(37)+SMB(150)+SMB(198))+Fa(123)*(SMB(171)+S
     &   MB(199))-F(335)*SMB(207)-Fa(138)*SMB(210)+Fa(137)*(SMB(209)+
     &   2*SMB(239))-F(452)*SMB(247)+Fa(142)*SMB(249)+p6mup4*F(445)*S
     &   MB(251)-p6mup4*F(443)*SMB(252)+F(445)*SMB(254)-p6mup4*F(443)
     &   *SMB(255)-F(445)*SMB(257)+F(443)*SMB(258)+F(443)*SMB(259))+F
     &   a(143)*SMB(279)+8*(Fa(67)*SMB(7)+Fa(68)*SMB(8)-Fa(72)*SMB(20
     &   )+Fa(74)*SMB(24)+Fa(75)*SMB(27)+Fa(76)*SMB(28)+Fa(78)*SMB(30
     &   )+Fa(79)*SMB(31)+Fa(80)*SMB(35)+Fa(93)*SMB(62)-Fa(112)*SMB(1
     &   18)-Fa(116)*SMB(129)-Fa(11)*(F(382)*SMB(61)+F(385)*SMB(63)+F
     &   (442)*SMB(130))+Fa(117)*SMB(131)+Fa(118)*SMB(132)+Fa(121)*SM
     &   B(162)+Fa(127)*SMB(186)-Fa(128)*SMB(188)+Fa(129)*SMB(189)-Fa
     &   (130)*SMB(192)-Fa(131)*SMB(193)-Fa(82)*(SMB(38)+SMB(163)+SMB
     &   (194))+Fa(134)*SMB(203)+Fa(136)*SMB(205)-Fa(135)*(SMB(204)+S
     &   MB(240))+Fa(139)*SMB(241)-Fa(140)*SMB(242)+Fa(120)*(SMB(139)
     &   +SMB(245))-Fa(141)*SMB(246)-F(445)*SMB(263)+Fa(122)*(-2*SMB(
     &   166)+SMB(280)+SMB(281)))+F(1308)*(P(90)*(6*SMB(78)+SMB(89))-
     &   SMB(268)-6*(SMB(260)-SMB(278))-2*(Fa(11)*(SMB(226)+SMB(265))
     &   -Fa(12)*SMB(266)-mup3mup4*SMB(267)+mup3mup4*SMB(282))+SMB(28
     &   3))+Fa(144)*SMB(284)+Fa(145)*(SMB(285)+SMB(286))+2*(F(258)*S
     &   MB(4)-F(261)*(2*SMB(11)+SMB(17))+F(360)*SMB(60)+F(392)*SMB(6
     &   8)+F(395)*SMB(74)-F(287)*(4*SMB(12)+SMB(18)-2*(SMB(83)+SMB(8
     &   6)))-Fa(96)*SMB(88)+Fa(94)*(2*SMB(64)+SMB(91))+Fa(102)*SMB(1
     &   00)+Fa(105)*SMB(103)+Fa(108)*SMB(110)+Fa(110)*SMB(112)-Fa(10
     &   6)*(SMB(104)+SMB(113))-F(397)*(SMB(77)-2*(SMB(87)-SMB(165)))
     &   -F(346)*(SMB(71)-2*(SMB(81)-SMB(179)))+F(343)*(SMB(173)+SMB(
     &   185)+2*(SMB(169)+SMB(196)))-F(267)*(SMB(153)+2*SMB(195)+SMB(
     &   201))-F(445)*(SMB(214)+2*(SMB(184)-p1mup4*SMB(213)-SMB(223))
     &   )-F(443)*(-SMB(217)-SMB(220)+2*(SMB(161)+p1mup4*(SMB(216)+SM
     &   B(219))+SMB(224)+SMB(225)))+p1mup4*F(445)*SMB(228)-mup3mup4*
     &   F(445)*SMB(229)-p1mup4*F(443)*SMB(232)+mup3mup4*F(443)*SMB(2
     &   33)-p1mup4*F(443)*SMB(236)+mup3mup4*F(443)*SMB(237)+F(394)*(
     &   SMB(72)+2*SMB(243))+Fa(11)*(-(F(445)*SMB(227))+F(443)*(SMB(2
     &   31)+SMB(235))-2*F(453)*SMB(248))+F(453)*(SMB(250)+4*SMB(256)
     &   )+p6mup4*F(445)*SMB(269)+Fa(133)*(-8*(SMB(208)-SMB(244))+4*(
     &   SMB(202)-SMB(261))-SMB(270))-p6mup4*F(443)*SMB(271)-p6mup4*F
     &   (443)*SMB(273)+Fa(119)*(8*SMB(138)+4*(SMB(157)+SMB(262)+SMB(
     &   264))+SMB(272)+SMB(274))-F(450)*(SMB(275)+4*SMB(276))-8*Fa(1
     &   46)*SMB(287)+2*(mup3mup4*(-(F(343)*(SMB(172)+SMB(183)))+F(26
     &   7)*(SMB(152)+SMB(200))+2*(F(443)*SMB(159)+F(445)*SMB(182)+F(
     &   448)*SMB(190)-F(442)*SMB(206)))-Fa(11)*(F(267)*SMB(149)+F(34
     &   3)*SMB(151)+2*(F(334)*SMB(141)-F(447)*SMB(187))+F(450)*SMB(2
     &   11)+F(445)*SMB(212)-F(443)*(SMB(215)+SMB(218)))-F(334)*(SMB(
     &   167)+SMB(288)))+Fa(147)*(6*SMB(289)+SMB(290))))/s45
       result(6) = Fa(148)*SMB(3)+F(549)*SMB(11)+Fa(153)*SMB(15)+F(594
     &   )*SMB(17)+Fa(155)*SMB(19)+Fa(164)*SMB(45)+Fa(166)*SMB(51)+Fa
     &   (167)*SMB(52)+Fa(168)*SMB(53)+Fa(170)*SMB(55)+Fa(171)*SMB(56
     &   )+Fa(172)*SMB(57)+Fa(173)*SMB(58)+16*(Fa(175)*SMB(62)+Fa(178
     &   )*SMB(65))+F(805)*SMB(69)+F(806)*SMB(70)+F(807)*SMB(71)+F(80
     &   8)*SMB(72)+F(811)*SMB(75)+F(814)*SMB(78)+F(815)*SMB(79)+F(81
     &   8)*SMB(82)+F(821)*SMB(85)+Fa(180)*SMB(88)+F(841)*SMB(89)+Fa(
     &   181)*SMB(90)+Fa(182)*SMB(91)+Fa(183)*SMB(92)+Fa(184)*SMB(93)
     &   +Fa(185)*SMB(94)+Fa(186)*SMB(95)+F(910)*SMB(96)+F(911)*SMB(9
     &   7)+F(912)*SMB(98)+F(913)*SMB(99)+Fa(187)*SMB(100)+Fa(188)*SM
     &   B(101)+Fa(190)*SMB(103)+Fa(191)*SMB(104)+Fa(192)*SMB(105)+F(
     &   973)*SMB(107)+F(974)*SMB(108)+F(975)*SMB(109)+Fa(193)*SMB(11
     &   0)+Fa(194)*SMB(111)+Fa(195)*SMB(112)+Fa(196)*SMB(113)+Fa(197
     &   )*SMB(114)+F(1030)*SMB(115)+F(1031)*SMB(116)+F(1032)*SMB(117
     &   )+F(1150)*SMB(230)+F(1249)*SMB(234)+F(1250)*SMB(238)+8*(Fa(1
     &   49)*SMB(5)+Fa(150)*SMB(6)+Fa(151)*SMB(7)+Fa(152)*SMB(8)+Fa(1
     &   56)*SMB(20)-Fa(157)*SMB(21)+Fa(158)*SMB(27)+Fa(159)*SMB(28)-
     &   Fa(160)*SMB(31)-Fa(161)*SMB(35)+Fa(162)*SMB(37)+Fa(165)*SMB(
     &   48)-Fa(174)*SMB(61)-Fa(176)*SMB(63)+Fa(177)*SMB(64)-Fa(179)*
     &   SMB(66)+F(803)*SMB(67)+Fa(198)*SMB(118)-Fa(199)*SMB(119)+F(1
     &   044)*SMB(125)+Fa(200)*SMB(126)+Fa(201)*SMB(128)+Fa(202)*SMB(
     &   129)-Fa(203)*SMB(130)+Fa(204)*SMB(131)-Fa(205)*SMB(132)-Fa(2
     &   06)*SMB(135)+Fa(207)*SMB(136)+F(1085)*SMB(137)+Fa(208)*SMB(1
     &   39)+Fa(209)*SMB(140)+Fa(210)*SMB(141)-Fa(211)*SMB(142)-Fa(21
     &   7)*SMB(162)+Fa(163)*(SMB(38)+SMB(163))+Fa(218)*SMB(164)-F(13
     &   14)*SMB(165)+Fa(221)*SMB(187)+Fa(222)*SMB(188)+Fa(223)*SMB(1
     &   89)+Fa(224)*SMB(193)-Fa(229)*SMB(203)-Fa(231)*SMB(205)+mup3m
     &   up4*(F(632)*SMB(23)+F(633)*SMB(24)+F(646)*SMB(29)-F(647)*SMB
     &   (30)-F(692)*SMB(39)+F(693)*SMB(40)+F(694)*SMB(41)-F(1074)*SM
     &   B(133)+F(1166)*SMB(190)-F(1204)*SMB(206))-Fa(242)*SMB(239)-F
     &   a(230)*(SMB(204)+SMB(240)))-F(1312)*SMB(268)+Fa(253)*SMB(280
     &   )+Fa(254)*SMB(281)+F(1312)*SMB(283)+Fa(255)*SMB(285)+Fa(256)
     &   *SMB(286)+2*(F(485)*SMB(4)+F(595)*SMB(18)+F(776)*SMB(59)+F(7
     &   77)*SMB(60)+F(804)*SMB(68)+F(809)*SMB(73)+F(810)*SMB(74)-F(8
     &   12)*SMB(76)+F(813)*SMB(77)+F(816)*SMB(80)-F(817)*SMB(81)+F(8
     &   19)*SMB(83)-F(820)*SMB(84)+F(822)*SMB(86)+F(823)*SMB(87)+F(9
     &   72)*SMB(106)-F(1116)*SMB(153)+F(1100)*SMB(160)+F(1144)*SMB(1
     &   73)+F(1138)*SMB(177)+F(1151)*SMB(185)-F(1186)*SMB(201)+F(124
     &   4)*SMB(222)-Fa(235)*SMB(227)-Fa(236)*SMB(228)-F(1150)*(SMB(1
     &   84)+SMB(214)-SMB(223)+mup3mup4*SMB(229))+Fa(237)*SMB(231)+Fa
     &   (238)*SMB(232)+F(1204)*(SMB(217)-SMB(224)+mup3mup4*SMB(233))
     &   +Fa(239)*SMB(235)+Fa(240)*SMB(236)+F(1127)*(-SMB(161)+SMB(22
     &   0)-SMB(225)+mup3mup4*SMB(237))+F(1266)*SMB(250)-F(1282)*SMB(
     &   260)+Fa(241)*(SMB(226)+SMB(265))+Fa(250)*SMB(266)-Fa(245)*SM
     &   B(269)+Fa(248)*SMB(270)-Fa(246)*SMB(271)-Fa(249)*SMB(272)-Fa
     &   (247)*SMB(273)+Fa(216)*SMB(274)-F(1286)*SMB(275)+F(1288)*SMB
     &   (278)+mup3mup4*F(1312)*(SMB(267)-SMB(282))+Fa(252)*SMB(284)+
     &   Fa(258)*SMB(290))+4*(F(457)*SMB(1)+F(458)*SMB(2)+F(547)*SMB(
     &   9)-F(548)*SMB(10)+F(550)*SMB(12)-F(551)*SMB(13)-F(552)*SMB(1
     &   4)-Fa(154)*SMB(16)+F(631)*SMB(22)+F(634)*SMB(25)+F(635)*SMB(
     &   26)+F(669)*SMB(32)+F(647)*SMB(33)+F(670)*SMB(34)+F(681)*SMB(
     &   36)-F(693)*SMB(43)+F(700)*SMB(46)+F(701)*SMB(47)+F(707)*SMB(
     &   49)-F(708)*SMB(50)+Fa(169)*SMB(54)-Fa(189)*SMB(102)+F(1043)*
     &   (SMB(120)-SMB(127))+F(1074)*SMB(134)+F(692)*(SMB(42)+SMB(143
     &   )-SMB(146))-Fa(212)*SMB(149)+Fa(213)*SMB(150)-Fa(214)*SMB(15
     &   1)+Fa(215)*SMB(154)+Fa(209)*SMB(155)-Fa(211)*SMB(156)-F(1100
     &   )*(SMB(144)-SMB(147)+mup3mup4*SMB(158))-F(694)*(SMB(44)-SMB(
     &   167)+SMB(169))+F(1138)*(-SMB(168)+SMB(170))+Fa(219)*SMB(171)
     &   +Fa(218)*SMB(174)-Fa(220)*SMB(180)-F(1166)*SMB(191)+Fa(225)*
     &   SMB(197)+Fa(226)*SMB(198)+Fa(227)*SMB(199)+mup3mup4*(F(1116)
     &   *SMB(152)+F(1127)*SMB(159)-F(1144)*SMB(172)-F(1138)*SMB(176)
     &   +F(1150)*SMB(182)-F(1151)*SMB(183)+F(1186)*SMB(200))-Fa(228)
     &   *SMB(202)+F(1204)*SMB(207)-Fa(232)*SMB(209)+Fa(233)*SMB(210)
     &   -Fa(234)*SMB(211)-Fa(235)*SMB(212)-Fa(236)*SMB(213)+Fa(237)*
     &   SMB(215)+Fa(238)*SMB(216)+Fa(239)*SMB(218)+Fa(240)*SMB(219)-
     &   F(1244)*SMB(221)-Fa(243)*SMB(248)+Fa(244)*SMB(249)-Fa(245)*S
     &   MB(251)-Fa(246)*SMB(252)-Fa(247)*SMB(255)+F(1281)*SMB(256)+F
     &   a(248)*SMB(261)-Fa(249)*SMB(262)+Fa(216)*(SMB(157)+SMB(264))
     &   -F(1287)*SMB(276)-Fa(251)*SMB(277)+Fa(252)*SMB(279)+Fa(257)*
     &   SMB(289)-Fa(259)*SMB(291))
       If (ngluon.eq.0) then
       resultn=(result(1)+result(2)+result(3)+result(4)+result(5)+res
     &   ult(6))
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
       resultborn = ((s123+s234-s45-s56)*SMB(4)+2*SMB(10)+2*(p1mup4+p4
     &   mup4+p5mup4+p6mup4)*(2*SMB(5)+SMB(51)+SMB(53))-SMB(60)-2*mup
     &   3mup4*((s123+s234-s45-s56)*SMB(3)+2*SMB(8)-SMB(58)+SMB(102))
     &   +SMB(106))/(s16*s23*s234)+(2*s123*SMB(2)-2*s23*SMB(2)+2*s234
     &   *SMB(2)-2*s45*SMB(2)-2*s56*SMB(2)+s123*SMB(4)-s23*SMB(4)+s23
     &   4*SMB(4)-s45*SMB(4)-s56*SMB(4)+4*p5mup4*SMB(5)+2*SMB(10)+4*p
     &   1mup4*SMB(27)+4*p5mup4*SMB(27)+4*p6mup4*SMB(27)+2*SMB(32)+2*
     &   SMB(33)+2*SMB(34)-2*SMB(50)+2*p5mup4*SMB(51)-2*p1mup4*SMB(52
     &   )-2*p6mup4*SMB(52)+2*p5mup4*SMB(53)+2*p4mup4*(2*SMB(5)+2*SMB
     &   (27)+SMB(51)+SMB(53))-SMB(60)-2*mup3mup4*(s123*SMB(3)-(s23-s
     &   234+s45+s56)*SMB(3)+2*SMB(8)+2*SMB(29)+2*SMB(30)+2*SMB(31)-S
     &   MB(58)+SMB(102))+SMB(106)+2*SMB(120))/(s16*s23*s45)+(-2*s234
     &   *SMB(1)+2*s56*SMB(1)+2*s234*SMB(2)-2*s56*SMB(2)-2*mup3mup4*s
     &   234*SMB(3)+2*mup3mup4*s56*SMB(3)+s23*(2*SMB(1)-2*SMB(2)+2*mu
     &   p3mup4*SMB(3)-SMB(4))+s234*SMB(4)-s56*SMB(4)+s123*(-2*SMB(1)
     &   +2*SMB(2)-2*mup3mup4*SMB(3)+SMB(4))+4*mup3mup4*SMB(24)-2*SMB
     &   (26)-4*mup3mup4*SMB(29)-4*mup3mup4*SMB(31)+2*SMB(32)+2*SMB(3
     &   4)+2*SMB(47)-2*SMB(50)+2*mup3mup4*SMB(58)-SMB(60)-2*mup3mup4
     &   *SMB(102)+SMB(106)+2*(p1mup4+p6mup4)*(2*SMB(27)-SMB(52)-2*SM
     &   B(118))+p1mup3*(-4*(p4mup4+p5mup4)*SMB(3)-4*SMB(28)+2*SMB(55
     &   )+4*SMB(119))+2*SMB(120)+2*p4mup4*(-2*p6mup3*SMB(3)-2*SMB(19
     &   )+2*SMB(45)+SMB(51)+SMB(53)+2*SMB(126))+2*(p6mup3*(-2*SMB(28
     &   )+SMB(55)+2*SMB(119))+p5mup4*(-2*p6mup3*SMB(3)-2*SMB(19)+2*S
     &   MB(45)+SMB(51)+SMB(53)+2*SMB(126)))-2*SMB(127))/(s16*s345*s4
     &   5)
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
