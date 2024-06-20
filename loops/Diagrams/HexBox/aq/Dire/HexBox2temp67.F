c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:9/3/2012
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
       subroutine HexBox267(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p5,barpsi_
     &   p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,result,resultn
     &   ,resultborn)
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
      Complex*16 SMB(300),Fa(385),F(1646)
       Complex*16 SMB0(32),SMB11(0:3),SMB12(0:3),SMB13(0:3),SMB14(0:3)
     &   ,SMB15(0:3),SMB16(0:3),SMB17(0:3),SMB18(0:3),SMB19(0:3),SMB1
     &   10(0:3)
      Complex*16 SMB1(133)

      Complex*16 barpsi_p1_P(2),psi_p5_P(2),barpsi_p2_P(2),psi_p4_P(2) 
       Complex*16 barpsi_p1_M(2),psi_p5_M(2),barpsi_p2_M(2),psi_p4_M(2)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p4(4)
       Real*8 P(1185)
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
       Complex*16 B0fin,C0fin,D0fin,E0fin,F0fin
       EXTERNAL B0fin,C0fin,D0fin,E0fin,F0fin
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

       COMMON/HexBox2FaFunctions/Fa
       COMMON/HexBox2FhlFunctions/F
      Save/HexBox2FhlFunctions/
       COMMON/HexBox2PFunctions/P
      Save/HexBox2PFunctions/
       COMMON/HexBox2Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,
     &   s34,s45,s56,s16,s123,s234,s345
       COMMON/HexBox2FVALFunctions/F0123456R,F123456R,F0123456I,F12345
     &   6I
       COMMON/HexBox2EVALFunctions/ E012345R,E12345R,E012345I,E12345I,
     &    E012346R,E12346R,E012346I,E12346I, E012356R,E12356R,E012356
     &   I,E12356I, E012456R,E12456R,E012456I,E12456I, E013456R,E1345
     &   6R,E013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/HexBox2DVALFunctions/ D01234R,D1234R,D01234I,D1234I, D01
     &   235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D1236I, D
     &   01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1246I,
     &    D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,D1345
     &   I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356I,D13
     &   56I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D02345I,D
     &   2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D02356I
     &   ,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,D0345
     &   6I,D3456I
       COMMON/HexBox2CVALFunctions/ C0123R,C123R,C0123I,C123I, C0124R,
     &   C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C126R,
     &   C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C0135I
     &   ,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,C145I
     &   , C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I, C023
     &   4R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236R,C23
     &   6R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246R,C02
     &   46I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C0345I,C3
     &   45I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C356I, C
     &   0456R,C456R,C0456I,C456I
       COMMON/HexBox2BVALFunctions/ B012R,B012I, B013R,B013I, B014R,B0
     &   14I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,B024I, B02
     &   5R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,B036I
     &   , B045R,B045I, B046R,B046I, B056R,B056I
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
       call HexBox2FFhl1(F(1))
       call HexBox2FFhl2(F(83))
       call HexBox2FFhl3(F(165))
       call HexBox2FFhl4(F(247))
       call HexBox2FFhl5(F(329))
       call HexBox2FFhl6(F(411))
       call HexBox2FFhl7(F(493))
       call HexBox2FFhl8(F(575))
       call HexBox2FFhl9(F(657))
       call HexBox2FFhl10(F(739))
       call HexBox2FFhl11(F(821))
       call HexBox2FFhl12(F(903))
       call HexBox2FFhl13(F(985))
       call HexBox2FFhl14(F(1067))
       call HexBox2FFhl15(F(1149))
       call HexBox2FFhl16(F(1231))
       call HexBox2FFhl17(F(1313))
       call HexBox2FFhl18(F(1395))
       call HexBox2FFhl19(F(1477))
       call HexBox2FFhl20(F(1559))
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
       call HexBox2Fa1(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mup
     &   3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(1))
       call HexBox2Fa2(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mup
     &   3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(20))
       call HexBox2Fa3(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mup
     &   3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(39))
       call HexBox2Fa4(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mup
     &   3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(58))
       call HexBox2Fa5(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mup
     &   3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(77))
       call HexBox2Fa6(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mup
     &   3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(96))
       call HexBox2Fa7(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mup
     &   3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(115))
       call HexBox2Fa8(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mup
     &   3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(134))
       call HexBox2Fa9(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mup
     &   3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(153))
       call HexBox2Fa10(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(172))
       call HexBox2Fa11(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(191))
       call HexBox2Fa12(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(210))
       call HexBox2Fa13(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(229))
       call HexBox2Fa14(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(248))
       call HexBox2Fa15(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(267))
       call HexBox2Fa16(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(286))
       call HexBox2Fa17(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(305))
       call HexBox2Fa18(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(324))
       call HexBox2Fa19(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(343))
       call HexBox2Fa20(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6,p4mu
     &   p3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(362))
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
       psi_p5_M(1)=psi_p5(1)
       psi_p5_M(2)=psi_p5(2)
       psi_p5_P(1)=psi_p5(3)
       psi_p5_P(2)=psi_p5(4)
       psi_p4_M(1)=psi_p4(1)
       psi_p4_M(2)=psi_p4(2)
       psi_p4_P(1)=psi_p4(3)
       psi_p4_P(2)=psi_p4(4)
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
       SMB0(1) = SC1c(barpsi_p1_P,mup3,psi_p5_P,+1)
       SMB0(2) = SC1c(barpsi_p2_M,mup6,psi_p4_M,-1)
       SMB0(3) = SC1r(barpsi_p2_M,p1,psi_p4_M,-1)
       SMB0(4) = SC1r(barpsi_p2_M,p5,psi_p4_M,-1)
       SMB0(5) = SC1r(barpsi_p2_M,p6,psi_p4_M,-1)
       SMB0(6) = SC3crr(barpsi_p2_M,mup6,p5,p1,psi_p4_M,-1)
       SMB0(7) = SC3crr(barpsi_p2_M,mup6,p6,p1,psi_p4_M,-1)
       SMB0(8) = SC3crr(barpsi_p2_M,mup6,p6,p5,psi_p4_M,-1)
       SMB0(9) = SC3rrr(barpsi_p2_M,p6,p5,p1,psi_p4_M,-1)
       SMB0(10) = SC1c(barpsi_p1_P,mup6,psi_p5_P,+1)
       SMB0(11) = SC1c(barpsi_p2_M,mup3,psi_p4_M,-1)
       SMB0(12) = SC3crr(barpsi_p2_M,mup3,p5,p1,psi_p4_M,-1)
       SMB0(13) = SC3crr(barpsi_p2_M,mup3,p6,p1,psi_p4_M,-1)
       SMB0(14) = SC3crr(barpsi_p2_M,mup3,p6,p5,psi_p4_M,-1)
       SMB0(15) = SC1r(barpsi_p1_P,p6,psi_p5_P,+1)
       SMB0(16) = SC3ccr(barpsi_p2_M,mup6,mup3,p1,psi_p4_M,-1)
       SMB0(17) = SC1r(barpsi_p1_P,p2,psi_p5_P,+1)
       SMB0(18) = SC3ccr(barpsi_p2_M,mup6,mup3,p5,psi_p4_M,-1)
       SMB0(19) = SC3ccr(barpsi_p2_M,mup6,mup3,p6,psi_p4_M,-1)
       SMB0(20) = SC5ccrrr(barpsi_p2_M,mup6,mup3,p6,p5,p1,psi_p4_M,-1)
       SMB0(21) = SC1r(barpsi_p1_P,p4,psi_p5_P,+1)
       SMB0(22) = SC3crr(barpsi_p1_P,mup6,p4,p2,psi_p5_P,+1)
       SMB0(23) = SC3crr(barpsi_p1_P,mup6,p6,p2,psi_p5_P,+1)
       SMB0(24) = SC3crr(barpsi_p1_P,mup6,p6,p4,psi_p5_P,+1)
       SMB0(25) = SC3rrr(barpsi_p1_P,p6,p4,p2,psi_p5_P,+1)
       SMB0(26) = SC3crr(barpsi_p1_P,mup3,p4,p2,psi_p5_P,+1)
       SMB0(27) = SC3crr(barpsi_p1_P,mup3,p6,p2,psi_p5_P,+1)
       SMB0(28) = SC3crr(barpsi_p1_P,mup3,p6,p4,psi_p5_P,+1)
       SMB0(29) = SC3ccr(barpsi_p1_P,mup6,mup3,p2,psi_p5_P,+1)
       SMB0(30) = SC3ccr(barpsi_p1_P,mup6,mup3,p4,psi_p5_P,+1)
       SMB0(31) = SC3ccr(barpsi_p1_P,mup6,mup3,p6,psi_p5_P,+1)
       SMB0(32) = SC5ccrrr(barpsi_p1_P,mup6,mup3,p6,p4,p2,psi_p5_P,+1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB11(i) = SC1c(barpsi_p1_P,v1,psi_p5_P,+1)
       SMB12(i) = SC1c(barpsi_p2_M,v1,psi_p4_M,-1)
       SMB13(i) = SC3crc(barpsi_p1_P,mup3,p6,v1,psi_p5_P,+1)
       SMB14(i) = SC3ccc(barpsi_p1_P,mup6,mup3,v1,psi_p5_P,+1)
       SMB15(i) = SC3crc(barpsi_p1_P,mup6,p6,v1,psi_p5_P,+1)
       SMB16(i) = SC3ccc(barpsi_p2_M,mup6,mup3,v1,psi_p4_M,-1)
       SMB17(i) = SC3crc(barpsi_p2_M,mup3,p6,v1,psi_p4_M,-1)
       SMB18(i) = SC3crc(barpsi_p2_M,mup6,p6,v1,psi_p4_M,-1)
       SMB19(i) = SC5ccrrc(barpsi_p1_P,mup6,mup3,p6,p2,v1,psi_p5_P,+1)
       SMB110(i) = SC5ccrrc(barpsi_p1_P,mup6,mup3,p6,p4,v1,psi_p5_P,+1)
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
       SMB1(18) = SMB0(11)*SMB0(15)
       SMB1(21) = SMB0(2)*SMB0(15)
       SMB1(23) = SMB0(3)*SMB0(15)
       SMB1(24) = SMB0(15)*SMB0(16)
       SMB1(25) = SMB0(11)*SMB0(17)
       SMB1(26) = SMB0(2)*SMB0(17)
       SMB1(27) = SMB0(3)*SMB0(17)
       SMB1(28) = SMB0(4)*SMB0(17)
       SMB1(29) = SMB0(5)*SMB0(17)
       SMB1(30) = SMB0(12)*SMB0(17)
       SMB1(31) = SMB0(13)*SMB0(17)
       SMB1(32) = SMB0(14)*SMB0(17)
       SMB1(33) = SMB0(16)*SMB0(17)
       SMB1(34) = SMB0(17)*SMB0(18)
       SMB1(35) = SMB0(17)*SMB0(19)
       SMB1(36) = SMB0(6)*SMB0(17)
       SMB1(37) = SMB0(7)*SMB0(17)
       SMB1(38) = SMB0(8)*SMB0(17)
       SMB1(39) = SMB0(9)*SMB0(17)
       SMB1(40) = SMB0(17)*SMB0(20)
       SMB1(41) = SMB0(11)*SMB0(21)
       SMB1(42) = SMB0(2)*SMB0(21)
       SMB1(43) = SMB0(3)*SMB0(21)
       SMB1(44) = SMB0(4)*SMB0(21)
       SMB1(45) = SMB0(5)*SMB0(21)
       SMB1(46) = SMB0(12)*SMB0(21)
       SMB1(47) = SMB0(13)*SMB0(21)
       SMB1(48) = SMB0(14)*SMB0(21)
       SMB1(49) = SMB0(16)*SMB0(21)
       SMB1(50) = SMB0(18)*SMB0(21)
       SMB1(51) = SMB0(19)*SMB0(21)
       SMB1(52) = SMB0(7)*SMB0(21)
       SMB1(53) = SMB0(8)*SMB0(21)
       SMB1(54) = SMB0(9)*SMB0(21)
       SMB1(55) = SMB0(20)*SMB0(21)
       SMB1(56) = SMB0(4)*SMB0(15)
       SMB1(57) = SMB0(5)*SMB0(15)
       SMB1(58) = SMB0(12)*SMB0(15)
       SMB1(59) = SMB0(13)*SMB0(15)
       SMB1(60) = SMB0(14)*SMB0(15)
       SMB1(61) = SMB0(15)*SMB0(18)
       SMB1(62) = SMB0(15)*SMB0(19)
       SMB1(63) = SMB0(6)*SMB0(15)
       SMB1(64) = SMB0(7)*SMB0(15)
       SMB1(65) = SMB0(8)*SMB0(15)
       SMB1(66) = SMB0(9)*SMB0(15)
       SMB1(67) = SMB0(15)*SMB0(20)
       SMB1(68) = SMB0(11)*SMB0(22)
       SMB1(69) = SMB0(11)*SMB0(23)
       SMB1(70) = SMB0(11)*SMB0(24)
       SMB1(71) = SMB0(11)*SMB0(25)
       SMB1(72) = SMB0(2)*SMB0(26)
       SMB1(73) = SMB0(2)*SMB0(27)
       SMB1(74) = SMB0(2)*SMB0(28)
       SMB1(75) = SMB0(2)*SMB0(25)
       SMB1(76) = SMB0(5)*SMB0(29)
       SMB1(78) = SMB0(3)*SMB0(26)
       SMB1(79) = SMB0(3)*SMB0(27)
       SMB1(80) = SMB0(3)*SMB0(28)
       SMB1(81) = SMB0(3)*SMB0(29)
       SMB1(82) = SMB0(3)*SMB0(30)
       SMB1(83) = SMB0(3)*SMB0(31)
       SMB1(84) = SMB0(3)*SMB0(22)
       SMB1(85) = SMB0(3)*SMB0(23)
       SMB1(86) = SMB0(3)*SMB0(24)
       SMB1(87) = SMB0(3)*SMB0(25)
       SMB1(88) = SMB0(3)*SMB0(32)
       SMB1(89) = SMB0(4)*SMB0(26)
       SMB1(90) = SMB0(4)*SMB0(27)
       SMB1(91) = SMB0(4)*SMB0(28)
       SMB1(92) = SMB0(4)*SMB0(29)
       SMB1(93) = SMB0(4)*SMB0(30)
       SMB1(94) = SMB0(4)*SMB0(31)
       SMB1(95) = SMB0(4)*SMB0(22)
       SMB1(96) = SMB0(4)*SMB0(23)
       SMB1(97) = SMB0(4)*SMB0(24)
       SMB1(98) = SMB0(4)*SMB0(25)
       SMB1(99) = SMB0(4)*SMB0(32)
       SMB1(100) = SMB0(5)*SMB0(26)
       SMB1(101) = SMB0(5)*SMB0(27)
       SMB1(102) = SMB0(5)*SMB0(28)
       SMB1(103) = SMB0(5)*SMB0(30)
       SMB1(104) = SMB0(5)*SMB0(31)
       SMB1(105) = SMB0(5)*SMB0(22)
       SMB1(106) = SMB0(5)*SMB0(23)
       SMB1(107) = SMB0(5)*SMB0(24)
       SMB1(108) = SMB0(5)*SMB0(25)
       SMB1(109) = SMB0(5)*SMB0(32)
       SMB1(112) = SMB0(6)*SMB0(27)
       SMB1(113) = SMB0(7)*SMB0(27)
       SMB1(114) = SMB0(8)*SMB0(27)
       SMB1(115) = SMB0(6)*SMB0(28)
       SMB1(116) = SMB0(7)*SMB0(28)
       SMB1(117) = SMB0(9)*SMB0(29)
       SMB1(118) = SMB0(9)*SMB0(30)
       SMB1(119) = SMB0(9)*SMB0(31)
       SMB1(122) = SMB0(13)*SMB0(22)
       SMB1(123) = SMB0(14)*SMB0(22)
       SMB1(125) = SMB0(12)*SMB0(23)
       SMB1(126) = SMB0(13)*SMB0(23)
       SMB1(127) = SMB0(14)*SMB0(23)
       SMB1(128) = SMB0(12)*SMB0(24)
       SMB1(129) = SMB0(13)*SMB0(24)
       SMB1(130) = SMB0(14)*SMB0(24)
       SMB1(131) = SMB0(16)*SMB0(25)
       SMB1(132) = SMB0(18)*SMB0(25)
       SMB1(133) = SMB0(19)*SMB0(25)
c************************************************************************************
       SMB1(17) = dotcc(SMB11,SMB12)
       SMB1(19) = dotcc(SMB12,SMB13)
       SMB1(20) = dotcc(SMB12,SMB14)
       SMB1(22) = dotcc(SMB12,SMB15)
       SMB1(77) = dotcc(SMB11,SMB16)
       SMB1(110) = dotcc(SMB11,SMB17)
       SMB1(111) = dotcc(SMB11,SMB18)
       SMB1(120) = dotcc(SMB12,SMB19)
       SMB1(121) = dotcc(SMB12,SMB110)
       SMB1(124) = dotcc(SMB15,SMB17)
c************************************************************************************
c************************************************************************************
       SMB(1) = 2*SMB1(1)-2*SMB1(9)+SMB1(20)
       SMB(2) = SMB1(20)
       SMB(3) = 4*SMB1(1)
       SMB(4) = SMB1(9)
       SMB(5) = SMB1(13)
       SMB(6) = SMB1(15)
       SMB(7) = 0
       SMB(8) = SMB1(44)
       SMB(9) = SMB1(50)
       SMB(10) = SMB1(56)
       SMB(11) = SMB1(61)
       SMB(12) = 0
       SMB(13) = 0
       SMB(14) = SMB1(11)
       SMB(15) = SMB1(25)
       SMB(16) = SMB1(28)
       SMB(17) = SMB1(30)
       SMB(18) = SMB1(32)
       SMB(19) = SMB1(41)
       SMB(20) = SMB1(46)
       SMB(21) = SMB1(48)
       SMB(22) = SMB1(18)
       SMB(23) = SMB1(58)
       SMB(24) = SMB1(60)
       SMB(25) = SMB1(68)
       SMB(26) = SMB1(69)
       SMB(27) = SMB1(70)
       SMB(28) = SMB1(95)
       SMB(29) = SMB1(96)
       SMB(30) = SMB1(97)
       SMB(31) = SMB1(123)
       SMB(32) = SMB1(125)
       SMB(33) = SMB1(127)
       SMB(34) = SMB1(128)
       SMB(35) = SMB1(130)
       SMB(36) = SMB1(10)
       SMB(37) = SMB1(12)
       SMB(38) = SMB1(14)
       SMB(39) = SMB1(17)
       SMB(40) = -2*SMB1(2)+2*p1mup3*SMB1(17)
       SMB(41) = -2*SMB1(4)+2*SMB1(18)+SMB1(19)
       SMB(42) = 0
       SMB(43) = (p6sq-s16+s234-s56)*SMB1(17)
       SMB(44) = 0
       SMB(45) = (p6sq-s16+s234-s56)*(-2*SMB1(4)+2*SMB1(18)+SMB1(19))
       SMB(46) = (p6sq-s16+s234-s56)*(2*SMB1(1)-2*SMB1(9)+SMB1(20))
       SMB(47) = 0
       SMB(48) = SMB1(42)
       SMB(49) = SMB1(43)
       SMB(50) = SMB1(45)
       SMB(51) = SMB1(49)
       SMB(52) = SMB1(51)
       SMB(53) = SMB1(53)
       SMB(54) = SMB1(54)
       SMB(55) = SMB1(55)
       SMB(56) = SMB1(21)
       SMB(57) = SMB1(23)
       SMB(58) = SMB1(57)
       SMB(59) = SMB1(24)
       SMB(60) = SMB1(62)
       SMB(61) = SMB1(63)
       SMB(62) = SMB1(65)
       SMB(63) = SMB1(66)
       SMB(64) = SMB1(67)
       SMB(65) = 0
       SMB(66) = SMB1(22)
       SMB(67) = 0
       SMB(68) = 0
       SMB(69) = 4*(p6sq-s56)*SMB1(1)+8*p5mup6*SMB1(4)+4*SMB1(15)-8*p5
     &   mup6*SMB1(18)-4*p5mup6*SMB1(19)+8*p5mup3*SMB1(21)-2*(p6sq-s5
     &   6)*SMB1(77)-4*SMB1(94)+4*p5mup6*SMB1(110)-4*p5mup3*SMB1(111)
     &   
       SMB(70) = -2*(s123-s23+s234-s56)*SMB1(1)-4*p4mup3*SMB1(10)+4*p1
     &   mup6*p4mup3*SMB1(17)-4*p1mup6*SMB1(41)+(s123-s23+s234-s56)*S
     &   MB1(77)+2*SMB1(82)
       SMB(71) = 2*(s123-s45-s56)*SMB1(1)-4*p4mup3*SMB1(21)-2*SMB1(70)
     &   +(-s123+s45+s56)*SMB1(77)+2*SMB1(103)+2*p4mup3*SMB1(111)
       SMB(72) = 0
       SMB(73) = 0
       SMB(74) = -((-p6sq+s16-s234+s56)*(2*s123*SMB1(9)-2*s45*SMB1(9)-
     &   2*s56*SMB1(9)-4*p4mup3*SMB1(12)-s123*SMB1(20)+s45*SMB1(20)+s
     &   56*SMB1(20)+2*p4mup3*SMB1(22)-2*SMB1(70)+2*SMB1(103)+SMB1(12
     &   1)))
       SMB(75) = -2*(p6sq-s16)*SMB1(1)-4*p1mup6*SMB1(4)-4*p6mup3*SMB1(
     &   10)+4*p1mup3*SMB1(12)+4*p1mup6*p6mup3*SMB1(17)-4*p1mup3*SMB1
     &   (21)+(p6sq-s16)*SMB1(77)+2*SMB1(83)-2*p1mup6*SMB1(110)+2*p1m
     &   up3*SMB1(111)
       SMB(76) = 2*p6sq*SMB1(1)-4*p6mup3*SMB1(21)-p6sq*SMB1(77)+2*SMB1
     &   (104)+2*p6mup3*SMB1(111)
       SMB(77) = (p6sq-s16+s234-s56)*(2*SMB1(12)-2*SMB1(21)+SMB1(111))
     &   
       SMB(78) = 0
       SMB(79) = (p6sq-s16+s234-s56)*SMB1(124)
       SMB(80) = SMB1(16)
       SMB(81) = -2*SMB1(10)+2*p1mup6*SMB1(17)
       SMB(82) = (-p6sq+s16)*SMB1(17)-2*SMB1(23)
       SMB(83) = -2*(p6sq-s16)*SMB1(1)+2*(p6sq-s16)*SMB1(9)+(-p6sq+s16
     &   )*SMB1(20)-2*SMB1(24)
       SMB(84) = SMB1(47)
       SMB(85) = SMB1(52)
       SMB(86) = SMB1(59)
       SMB(87) = SMB1(64)
       SMB(88) = 0
       SMB(89) = SMB1(86)
       SMB(90) = 0
       SMB(91) = 0
       SMB(92) = 2*(-s123+s45+s56)*SMB1(10)+2*p1mup6*(s123-s45-s56)*SM
     &   B1(17)-2*(s123-s23+s234-s56)*SMB1(21)-4*p1mup6*SMB1(45)+2*SM
     &   B1(86)+(s123-s23+s234-s56)*SMB1(111)
       SMB(93) = 2*s16*SMB1(12)-2*p6sq*(SMB1(10)+SMB1(12)-p1mup6*SMB1(
     &   17))-4*p1mup6*SMB1(57)
       SMB(94) = SMB1(129)
       SMB(95) = 0
       SMB(96) = 0
       SMB(97) = 0
       SMB(98) = 0
       SMB(99) = 0
       SMB(100) = SMB1(19)
       SMB(101) = 2*(2*SMB1(18)+SMB1(19)-SMB1(110))
       SMB(102) = 0
       SMB(103) = 0
       SMB(104) = 0
       SMB(105) = 0
       SMB(106) = 0
       SMB(107) = 0
       SMB(108) = 0
       SMB(109) = 0
       SMB(110) = 0
       SMB(111) = 0
       SMB(112) = 0
       SMB(113) = 0
       SMB(114) = 0
       SMB(115) = 4*p4mup3*(-2*SMB1(12)+2*SMB1(21)+SMB1(22)-SMB1(111))
     &   
       SMB(116) = 0
       SMB(117) = 0
       SMB(118) = 0
       SMB(119) = SMB1(26)
       SMB(120) = SMB1(27)
       SMB(121) = SMB1(29)
       SMB(122) = SMB1(33)
       SMB(123) = SMB1(34)
       SMB(124) = SMB1(35)
       SMB(125) = SMB1(36)
       SMB(126) = SMB1(39)
       SMB(127) = SMB1(40)
       SMB(128) = SMB1(71)
       SMB(129) = 2*p2mup6*SMB1(17)-2*SMB1(26)
       SMB(130) = (p3sq-s23+s234-s34)*SMB1(17)
       SMB(131) = -((s12+s16-s345)*SMB1(17))-2*SMB1(29)
       SMB(132) = SMB1(84)
       SMB(133) = SMB1(87)
       SMB(134) = 2*s12*SMB1(1)-4*p2mup3*SMB1(10)+4*p1mup3*p2mup6*SMB1
     &   (17)-4*p1mup3*SMB1(26)-s12*SMB1(77)+2*SMB1(81)
       SMB(135) = 0
       SMB(136) = -2*(s12+s16-s345)*SMB1(1)-(s12+s16-s345)*SMB1(20)-4*
     &   p2mup3*SMB1(21)-2*SMB1(69)+(s12+s16-s345)*SMB1(77)+2*p2mup3*
     &   SMB1(111)
       SMB(137) = 2*(p6sq-s16+s234-s56)*(p2mup6*SMB1(17)-SMB1(26))
       SMB(138) = 2*(s12+s16-s345)*SMB1(10)+2*p2mup6*(-p6sq+s16)*SMB1(
     &   17)+2*s12*SMB1(21)+2*(p6sq-s16)*SMB1(26)+2*SMB1(85)-s12*SMB1
     &   (111)
       SMB(139) = (-p6sq+s16-s234+s56)*(2*s12*SMB1(9)+2*s16*SMB1(9)-2*
     &   s345*SMB1(9)+4*p2mup3*SMB1(12)-s12*SMB1(20)-s16*SMB1(20)+s34
     &   5*SMB1(20)-2*p2mup3*SMB1(22)+2*SMB1(69)-2*SMB1(76)-SMB1(120)
     &   )
       SMB(140) = SMB1(122)
       SMB(141) = 2*p2mup3*SMB1(17)-2*SMB1(25)
       SMB(142) = 0
       SMB(143) = 0
       SMB(144) = 0
       SMB(145) = 0
       SMB(146) = -2*(p3sq-s23+s234-s34)*(SMB1(2)-p1mup3*SMB1(17))
       SMB(147) = 0
       SMB(148) = (p3sq-s23+s234-s34)*(-2*SMB1(4)+2*SMB1(18)+SMB1(19))
     &   
       SMB(149) = (p3sq-s23+s234-s34)*(2*SMB1(1)-2*SMB1(9)+SMB1(20))
       SMB(150) = 2*(s12+s16-s345)*SMB1(2)+2*s12*SMB1(4)-2*p1mup3*(s12
     &   +s16-s345)*SMB1(17)-s12*SMB1(19)+2*(p6sq-s16)*SMB1(25)+2*SMB
     &   1(31)
       SMB(151) = 2*(s16-s234+s34-s345)*SMB1(4)+2*p2mup3*(-p6sq+s56)*S
     &   MB1(17)+(-s16+s234-s34+s345)*SMB1(19)+2*(p6sq-s56)*SMB1(25)+
     &   2*SMB1(32)-4*p2mup3*SMB1(56)
       SMB(152) = 2*p2mup3*p6sq*SMB1(17)-2*(s12+s16-s345)*SMB1(18)-2*p
     &   6sq*SMB1(25)-4*p2mup3*SMB1(57)
       SMB(153) = -2*(s12+s16-s345)*SMB1(1)-4*p2mup6*SMB1(4)+2*(s12+s1
     &   6-s345)*SMB1(9)+4*p2mup3*SMB1(12)+2*p2mup6*SMB1(19)-(s12+s16
     &   -s345)*SMB1(20)-2*p2mup3*SMB1(22)-2*SMB1(35)
       SMB(154) = 2*(s12+s16-s345)*SMB1(10)+2*s12*SMB1(12)-2*p1mup6*(s
     &   12+s16-s345)*SMB1(17)-s12*SMB1(22)+2*(p6sq-s16)*SMB1(26)+2*S
     &   MB1(37)
       SMB(155) = (s12*(s16-s234)+p6sq*(-s16+s345)+(s16-s345)*(s16-s23
     &   4+s56))*SMB1(17)-2*(p6sq-s56)*SMB1(27)+2*(p6sq-s16)*SMB1(28)
     &   -2*SMB1(39)+2*s12*SMB1(56)
       SMB(156) = (p6sq*(s16-s345)-s16*(s12+s16-s345))*SMB1(17)+2*(s12
     &   +s16-s345)*SMB1(23)+2*p6sq*SMB1(27)+2*(p6sq-s16)*SMB1(29)+2*
     &   s12*SMB1(57)
       SMB(157) = 2*(s16**2-s16*(s234-s34+s345)+p6sq*(s12-s16+s234-s34
     &   +s345)-s12*s56)*SMB1(1)+4*p2mup6*(-p6sq+s56)*SMB1(2)-2*(s16*
     &   *2-s16*(s234-s34+s345)+p6sq*(s12-s16+s234-s34+s345)-s12*s56)
     &   *SMB1(9)+4*p2mup3*(p6sq-s56)*SMB1(10)+4*(p1mup6*p2mup3-p1mup
     &   3*p2mup6)*(-p6sq+s56)*SMB1(17)+(s16**2-s16*(s234-s34+s345)+p
     &   6sq*(s12-s16+s234-s34+s345)-s12*s56)*SMB1(20)-2*(p6sq-s56)*S
     &   MB1(33)+2*(p6sq-s16)*SMB1(34)-2*SMB1(40)-4*p2mup6*SMB1(58)+2
     &   *s12*SMB1(61)+4*p2mup3*SMB1(63)+(-(p6sq*(s12+s234-s34))-s16*
     &   s34+s234*s345+s16*s56-s345*s56+s12*(s16-s234+s56))*SMB1(77)+
     &   2*p2mup6*(p6sq-s16+s234-s56)*SMB1(110)+2*p2mup3*(-p6sq+s16-s
     &   234+s56)*SMB1(111)
       SMB(158) = 2*(p6sq*(s16-s345)-s16*(s12+s16-s345))*SMB1(1)+4*p2m
     &   up6*p6sq*SMB1(2)+4*p2mup6*(p6sq-s16)*SMB1(4)+2*(s16*(s12+s16
     &   -s345)+p6sq*(-s16+s345))*SMB1(9)-4*p2mup3*p6sq*SMB1(10)+4*p2
     &   mup3*(-p6sq+s16)*SMB1(12)+4*(p1mup6*p2mup3-p1mup3*p2mup6)*p6
     &   sq*SMB1(17)+4*p2mup6*(-p6sq+s16)*SMB1(18)+2*p2mup6*(-p6sq+s1
     &   6)*SMB1(19)+(p6sq*(s16-s345)-s16*(s12+s16-s345))*SMB1(20)+4*
     &   p2mup3*(p6sq-s16)*SMB1(21)+2*p2mup3*(p6sq-s16)*SMB1(22)+2*(s
     &   12+s16-s345)*SMB1(24)+2*p6sq*SMB1(33)+2*(p6sq-s16)*SMB1(35)-
     &   4*p2mup6*SMB1(59)+2*s12*SMB1(62)+4*p2mup3*SMB1(64)
       SMB(159) = SMB1(131)
       SMB(160) = 2*(s12+s16-s345)*SMB1(1)-2*(s12+s16-s345)*SMB1(9)-2*
     &   SMB1(76)-(s12+s16-s345)*SMB1(77)
       SMB(161) = -8*p2mup6*p6mup3*SMB1(17)-2*(s12+s16-s345)*SMB1(20)+
     &   4*p2mup3*SMB1(22)+8*p6mup3*SMB1(26)-4*SMB1(69)-4*SMB1(76)
       SMB(162) = (p3sq-s23+s234-s34)*(-2*(p6sq-s16)*SMB1(9)-4*p6mup3*
     &   SMB1(10)+4*p1mup6*p6mup3*SMB1(17)-4*p1mup6*SMB1(18)-2*p1mup6
     &   *SMB1(19)+(p6sq-s16)*SMB1(20)+2*p1mup3*SMB1(22)+2*SMB1(83))
       SMB(163) = 0
       SMB(164) = -4*(p6sq-s16)*(s16-s234+s34-s345)*SMB1(9)+8*p2mup3*(
     &   p6sq-s16)*SMB1(11)-8*p1mup6*(s16-s234+s34-s345)*SMB1(18)+16*
     &   p1mup6*p2mup3*SMB1(56)
       SMB(165) = 4*(p6sq-s16)*(s12+s16-s345)*SMB1(9)+8*p2mup3*(p6sq-s
     &   16)*SMB1(12)+8*p1mup6*(s12+s16-s345)*SMB1(18)+16*p1mup6*p2mu
     &   p3*SMB1(57)
       SMB(166) = SMB1(2)
       SMB(167) = SMB1(3)
       SMB(168) = SMB1(4)
       SMB(169) = SMB1(81)
       SMB(170) = SMB1(83)
       SMB(171) = SMB1(92)
       SMB(172) = SMB1(94)
       SMB(173) = SMB1(76)
       SMB(174) = SMB1(104)
       SMB(175) = SMB1(107)
       SMB(176) = SMB1(85)
       SMB(177) = SMB1(106)
       SMB(178) = SMB1(119)
       SMB(179) = SMB1(126)
       SMB(180) = SMB1(1)
       SMB(181) = SMB1(5)
       SMB(182) = SMB1(6)
       SMB(183) = SMB1(8)
       SMB(184) = -2*SMB1(12)+2*SMB1(21)+SMB1(22)
       SMB(185) = SMB1(31)
       SMB(186) = (p6sq-s16+s234-s56)*(-2*SMB1(1)+2*SMB1(9)+SMB1(77))
       SMB(187) = 2*(p6sq-s16)*SMB1(1)-2*(p6sq-s16)*SMB1(9)-4*p1mup6*S
     &   MB1(18)+4*p1mup3*SMB1(21)+(-p6sq+s16)*SMB1(77)-2*SMB1(83)+2*
     &   p1mup6*SMB1(110)-2*p1mup3*SMB1(111)
       SMB(188) = SMB1(117)
       SMB(189) = 2*(2*SMB1(12)-SMB1(22)+SMB1(111))
       SMB(190) = 0
       SMB(191) = SMB1(73)
       SMB(192) = SMB1(74)
       SMB(193) = SMB1(103)
       SMB(194) = (p3sq-s23+s234-s34)*SMB1(22)
       SMB(195) = 0
       SMB(196) = (p3sq-s23+s234-s34)*SMB1(124)
       SMB(197) = SMB1(79)
       SMB(198) = SMB1(90)
       SMB(199) = SMB1(101)
       SMB(200) = SMB1(75)
       SMB(201) = SMB1(98)
       SMB(202) = SMB1(105)
       SMB(203) = SMB1(108)
       SMB(204) = SMB1(132)
       SMB(205) = SMB1(133)
       SMB(206) = SMB1(82)
       SMB(207) = SMB1(93)
       SMB(208) = 0
       SMB(209) = SMB1(7)
       SMB(210) = SMB1(78)
       SMB(211) = SMB1(80)
       SMB(212) = SMB1(88)
       SMB(213) = SMB1(89)
       SMB(214) = SMB1(91)
       SMB(215) = SMB1(99)
       SMB(216) = SMB1(100)
       SMB(217) = SMB1(102)
       SMB(218) = SMB1(109)
       SMB(219) = SMB1(37)
       SMB(220) = 0
       SMB(221) = 2*(p6sq-s16+s234-s56)*(p2mup3*SMB1(17)-SMB1(25))
       SMB(222) = 2*(s12+s16-s345)*SMB1(2)+2*p2mup3*(-p6sq+s16)*SMB1(1
     &   7)+2*s12*SMB1(18)+2*(p6sq-s16)*SMB1(25)+2*SMB1(79)-s12*SMB1(
     &   110)
       SMB(223) = 0
       SMB(224) = 2*(-s123+s45+s56)*SMB1(2)+2*p1mup3*(s123-s45-s56)*SM
     &   B1(17)-2*(s123-s23+s234-s56)*SMB1(18)-4*p1mup3*SMB1(45)+2*SM
     &   B1(80)+(s123-s23+s234-s56)*SMB1(110)
       SMB(225) = -4*p6mup6*SMB1(2)+4*p1mup6*SMB1(4)-2*(p6sq-s16)*SMB1
     &   (9)-4*p1mup3*SMB1(12)+2*(2*p1mup3*p6mup6+mup3mup6*(p6sq-s16)
     &   )*SMB1(17)-4*p1mup6*SMB1(18)+4*mup3mup6*SMB1(23)+(-p6sq+s16)
     &   *SMB1(77)-2*SMB1(83)+2*p1mup6*SMB1(110)-2*p1mup3*SMB1(111)
       SMB(226) = (p6sq-s16+s234-s56)*(2*SMB1(4)-2*SMB1(18)+SMB1(110))
     &   
       SMB(227) = -2*(-(s16*SMB1(4))+p6sq*(SMB1(2)+SMB1(4)-p1mup3*SMB1
     &   (17))+2*p1mup3*SMB1(57))
       SMB(228) = SMB1(112)
       SMB(229) = SMB1(113)
       SMB(230) = SMB1(115)
       SMB(231) = SMB1(116)
       SMB(232) = SMB1(118)
       SMB(233) = -8*SMB1(2)+8*p1mup3*SMB1(17)
       SMB(234) = 0
       SMB(235) = -8*p1mup6*SMB1(4)+8*p1mup3*SMB1(12)
       SMB(236) = -2*(s16**2-s16*(s234-s34+s345)+p6sq*(s12-s16+s234-s3
     &   4+s345)-s12*s56)*SMB1(1)+4*p1mup6*(s16-s234+s34-s345)*SMB1(4
     &   )+2*(s16**2-s16*(s234-s34+s345)+p6sq*(s12-s16+s234-s34+s345)
     &   -s12*s56)*SMB1(9)-4*p1mup3*(s16-s234+s34-s345)*SMB1(12)-2*p1
     &   mup6*(s16-s234+s34-s345)*SMB1(19)+(-(p6sq*(s12+s234-s34))-s1
     &   6*s34+s234*s345+s16*s56-s345*s56+s12*(s16-s234+s56))*SMB1(20
     &   )+2*p1mup3*(s16-s234+s34-s345)*SMB1(22)+(s16**2-s16*(s234-s3
     &   4+s345)+p6sq*(s12-s16+s234-s34+s345)-s12*s56)*SMB1(77)-2*(p6
     &   sq-s56)*SMB1(81)+2*(p6sq-s16)*SMB1(92)+2*s12*SMB1(94)-2*SMB1
     &   (117)
       SMB(237) = 0
       SMB(238) = (p6sq*(s16-s345)-s16*(s12+s16-s345))*SMB1(20)+2*(p6s
     &   q-s16)*SMB1(76)+2*p6sq*SMB1(81)+2*(s12+s16-s345)*SMB1(83)+2*
     &   s12*SMB1(104)
       SMB(239) = 0
       SMB(240) = 4*(-2*(p6sq-s16)*SMB1(9)-4*p6mup3*SMB1(10)+4*p1mup6*
     &   p6mup3*SMB1(17)-4*p1mup6*SMB1(18)-2*p1mup6*SMB1(19)+(p6sq-s1
     &   6)*SMB1(20)+2*p1mup3*SMB1(22)+2*SMB1(83))
       SMB(241) = SMB1(38)
       SMB(242) = 0
       SMB(243) = 2*(2*SMB1(4)-SMB1(19)+SMB1(110))
       SMB(244) = 0
       SMB(245) = -4*p6mup6*SMB1(18)+4*mup3mup6*SMB1(57)+p6sq*(2*SMB1(
     &   9)-2*mup3mup6*SMB1(17)+SMB1(77))-2*SMB1(104)+2*p6mup6*SMB1(1
     &   10)
       SMB(246) = 4*SMB1(17)
       SMB(247) = 0
       SMB(248) = -8*SMB1(4)+8*SMB1(18)+4*SMB1(19)
       SMB(249) = 0
       SMB(250) = 0
       SMB(251) = 0
       SMB(252) = 0
       SMB(253) = 2*(2*SMB1(21)+SMB1(22)-SMB1(111))
       SMB(254) = 0
       SMB(255) = 4*SMB1(57)
       SMB(256) = 0
       SMB(257) = 4*SMB1(62)
       SMB(258) = 2*(s16-s234+s34-s345)*SMB1(12)+2*p2mup6*(-p6sq+s56)*
     &   SMB1(17)+(-s16+s234-s34+s345)*SMB1(22)+2*(p6sq-s56)*SMB1(26)
     &   +2*SMB1(38)-4*p2mup6*SMB1(56)
       SMB(259) = 2*p2mup6*p6sq*SMB1(17)-2*(s12+s16-s345)*SMB1(21)-2*p
     &   6sq*SMB1(26)-4*p2mup6*SMB1(57)
       SMB(260) = 0
       SMB(261) = 0
       SMB(262) = 4*(2*SMB1(12)-2*SMB1(21)+SMB1(111))
       SMB(263) = 0
       SMB(264) = 0
       SMB(265) = 0
       SMB(266) = 4*SMB1(104)
       SMB(267) = 0
       SMB(268) = 0
       SMB(269) = 0
       SMB(270) = p6sq*(8*SMB1(9)-4*SMB1(20))+8*(p6mup3*(-2*SMB1(12)+S
     &   MB1(22))+SMB1(104))
       SMB(271) = 0
       SMB(272) = 0
       SMB(273) = 2*(p6sq*(s12+s234-s34)-(s12+s16-s345)*s56)*SMB1(1)+4
     &   *p2mup6*(p6sq-s56)*SMB1(4)+2*(-(p6sq*(s12+s234-s34))+(s12+s1
     &   6-s345)*s56)*SMB1(9)+4*p2mup3*(-p6sq+s56)*SMB1(12)+4*p2mup6*
     &   (-p6sq+s56)*SMB1(18)+2*p2mup6*(-p6sq+s56)*SMB1(19)+(p6sq*(s1
     &   2+s234-s34)-(s12+s16-s345)*s56)*SMB1(20)+4*p2mup3*(p6sq-s56)
     &   *SMB1(21)+2*p2mup3*(p6sq-s56)*SMB1(22)+2*p6sq*SMB1(34)+2*(p6
     &   sq-s56)*SMB1(35)-4*p2mup6*SMB1(60)+2*(s12+s16-s345)*SMB1(61)
     &   +2*(s16-s234+s34-s345)*SMB1(62)+4*p2mup3*SMB1(65)
       SMB(274) = 0
       SMB(275) = 2*(-(p6sq*(s12+s234-s34))+(s12+s16-s345)*s56)*SMB1(1
     &   )+2*(p6sq*(s12+s234-s34)-(s12+s16-s345)*s56)*SMB1(9)+2*(p6sq
     &   -s56)*SMB1(76)+(p6sq*(s12+s234-s34)-(s12+s16-s345)*s56)*SMB1
     &   (77)+2*p6sq*SMB1(92)+2*(s12+s16-s345)*SMB1(94)+2*(s16-s234+s
     &   34-s345)*SMB1(104)
       SMB(276) = 0
       SMB(277) = 0
       SMB(278) = 0
       SMB(279) = -4*p2mup6*SMB1(2)+4*p1mup6*p2mup3*SMB1(17)-4*p1mup6*
     &   SMB1(25)+4*mup3mup6*SMB1(27)+s12*(2*SMB1(9)-2*mup3mup6*SMB1(
     &   17)+SMB1(77))-2*SMB1(81)
       SMB(280) = 0
       SMB(281) = -2*(s12+s16-s345)*SMB1(9)-4*p2mup6*SMB1(18)+(s12+s16
     &   -s345)*SMB1(20)-2*SMB1(73)-(s12+s16-s345)*SMB1(77)+2*p2mup6*
     &   SMB1(110)
       SMB(282) = 0
       SMB(283) = SMB1(114)
       SMB(284) = 0
       SMB(285) = 0
       SMB(286) = 0
       SMB(287) = 0
       SMB(288) = -8*p2mup6*SMB1(18)+8*p2mup3*SMB1(21)
       SMB(289) = SMB1(72)
       SMB(290) = -4*p4mup6*SMB1(2)-2*(s123-s23+s234-s56)*SMB1(9)+2*(2
     &   *p1mup3*p4mup6+mup3mup6*(s123-s23+s234-s56))*SMB1(17)-4*p1mu
     &   p3*SMB1(42)+4*mup3mup6*SMB1(43)+(-s123+s23-s234+s56)*SMB1(77
     &   )-2*SMB1(82)
       SMB(291) = 0
       SMB(292) = 2*(s123-s45-s56)*SMB1(9)+2*mup3mup6*(-s123+s45+s56)*
     &   SMB1(17)-4*p4mup6*SMB1(18)+4*mup3mup6*SMB1(45)-2*SMB1(74)+(s
     &   123-s45-s56)*SMB1(77)-2*SMB1(103)+2*p4mup6*SMB1(110)
       SMB(293) = 0
       SMB(294) = 0
       SMB(295) = 0
       SMB(296) = 0
       SMB(297) = 0
       SMB(298) = 0
       SMB(299) = 0
       SMB(300) = 4*SMB1(9)
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
       result(1) = (Fa(1)*SMB(7)+F(8)*SMB(13)+F(12)*SMB(27)+8*F(11)*SM
     &   B(35)+4*(F(10)*SMB(34)-Fa(2)*SMB(40))+s16*(s23*F(2)*SMB(3)+4
     &   *(F(3)*(-(P(1)*SMB(4))+SMB(5))-F(4)*SMB(6))+F(6)*(4*mup3mup6
     &   *(SMB(43)-SMB(44))-2*(SMB(46)-SMB(47))+SMB(68)-SMB(69))+2*(F
     &   (1)*(-SMB(1)+SMB(2))-F(16)*(SMB(70)+SMB(71))))-2*(F(7)*SMB(1
     &   2)+s45*(F(13)*SMB(75)+F(14)*SMB(76)))+F(15)*(-4*p6mup3*SMB(6
     &   5)+2*SMB(114)+SMB(115))+F(17)*SMB(117)+F(1634)*(SMB(240)+2*p
     &   1mup6*(SMB(233)+SMB(247)+SMB(248))+SMB(269)+SMB(270))+F(18)*
     &   SMB(278)+p1mup6*(F(9)*SMB(19)+2*(4*F(10)*SMB(20)+8*F(11)*SMB
     &   (21)-2*s45*F(14)*SMB(41)+F(15)*(-2*SMB(95)+SMB(96))+F(17)*(S
     &   MB(98)+2*SMB(106)-SMB(295)))))/(s16*s23)
       result(2) = (Fa(3)*SMB(7)+F(25)*SMB(13)+Fa(4)*SMB(14)+Fa(11)*SM
     &   B(43)+Fa(13)*SMB(66)+Fa(16)*SMB(129)+F(45)*SMB(161)+F(46)*SM
     &   B(164)+p1mup6*(4*F(33)*SMB(41)+16*F(30)*(SMB(17)-2*p5mup3*(S
     &   MB(120)+SMB(121)))+2*(F(45)*(2*SMB(141)+SMB(142))+F(38)*(-SM
     &   B(144)+2*(p5mup3*SMB(143)-SMB(151))))+F(47)*SMB(167)+8*(F(28
     &   )*SMB(15)-F(29)*SMB(18)+F(32)*(SMB(166)+SMB(168))))+F(48)*SM
     &   B(172)+8*(Fa(7)*(SMB(36)+SMB(37))-Fa(8)*SMB(39)+F(30)*(SMB(3
     &   2)-2*p5mup3*(SMB(176)+SMB(177))))+s16*(F(21)*SMB(3)+F(22)*SM
     &   B(4)+2*(F(19)*SMB(1)-2*F(26)*SMB(171)+F(31)*(SMB(134)+SMB(13
     &   6)-2*(SMB(169)+SMB(173))))-F(23)*(2*(SMB(46)-SMB(47))-SMB(68
     &   )+SMB(69)+4*(SMB(6)+SMB(186))))+Fa(15)*(SMB(90)+SMB(189))+2*
     &   (F(20)*SMB(2)+Fa(10)*(-2*SMB(42)+SMB(67))+F(40)*SMB(75)+F(41
     &   )*SMB(76)+F(44)*SMB(135)+F(35)*SMB(160)+F(38)*(2*p2mup3*(SMB
     &   (77)-SMB(78))+p5mup3*SMB(163)-4*(p6mup3*SMB(137)-p6mup3*SMB(
     &   190))))+4*(F(28)*SMB(26)+F(29)*(4*Fa(5)*SMB(16)+2*Fa(6)*SMB(
     &   29)-SMB(33))-Fa(9)*SMB(40)-Fa(12)*SMB(44)-F(32)*(SMB(170)+SM
     &   B(174))+Fa(14)*(SMB(81)+SMB(184))+F(49)*SMB(208))+F(1633)*(1
     &   2*SMB(187)+SMB(240)-4*Fa(17)*SMB(246)+2*(p1mup6*(SMB(233)+6*
     &   (SMB(220)+SMB(242)+SMB(243))+SMB(247)+SMB(248))-Fa(18)*SMB(2
     &   62))-6*(SMB(235)+SMB(263)+SMB(266))+SMB(269)+SMB(270))+F(38)
     &   *(4*(SMB(236)-Fa(5)*SMB(249)-SMB(275))+SMB(277)-2*(Fa(6)*SMB
     &   (264)-p5mup3*SMB(267)-p1mup6*(2*p5mup3*SMB(250)-4*(SMB(221)-
     &   SMB(282))+SMB(285)))))/(s16*s34)
       result(3) = (F(50)*SMB(1)+F(51)*SMB(2)+F(52)*SMB(3)+F(61)*SMB(9
     &   )-16*p1mup6*F(67)*SMB(21)+8*(F(54)*SMB(5)-mup3mup6*F(60)*SMB
     &   (8)-Fa(20)*SMB(19)+Fa(22)*SMB(22))+F(72)*SMB(27)+Fa(24)*SMB(
     &   41)+F(82)*SMB(52)+F(80)*(SMB(68)+2*(SMB(46)+SMB(83)))+mup3mu
     &   p6*(F(81)*SMB(50)+16*F(67)*SMB(54)+4*F(59)*(SMB(43)+SMB(82))
     &   +8*(F(75)*SMB(39)+s56*F(65)*SMB(49)-F(70)*SMB(57)-F(69)*SMB(
     &   58)+F(83)*SMB(88)))+F(86)*SMB(91)+Fa(21)*(-8*(F(66)*SMB(20)-
     &   F(73)*SMB(84))+4*F(87)*SMB(95))+2*(F(62)*SMB(12)+2*(F(70)*SM
     &   B(59)+F(69)*SMB(60))-s56*F(65)*(2*SMB(51)+SMB(70))+F(84)*SMB
     &   (71)+F(85)*SMB(75)+F(63)*(SMB(13)+SMB(76))+F(83)*Fa(21)*SMB(
     &   96)+Fa(26)*SMB(101)+F(67)*(-4*SMB(55)-p4mup6*(SMB(103)+SMB(1
     &   04))+Fa(21)*(SMB(98)+SMB(99)-4*SMB(105))-2*(SMB(74)-p1mup6*(
     &   SMB(106)+SMB(107)))))+4*(F(53)*SMB(4)-Fa(19)*SMB(7)+F(55)*(S
     &   MB(6)-2*mup3mup6*SMB(10)+SMB(11))+F(74)*SMB(38)-Fa(23)*SMB(4
     &   0)+p6mup3*F(83)*SMB(65)+p4mup6*(F(67)*SMB(45)+4*F(71)*(SMB(2
     &   3)+SMB(86)))-F(73)*(SMB(34)+SMB(94))+Fa(25)*SMB(100)+F(67)*(
     &   2*mup3mup6*(SMB(110)+SMB(111))-SMB(112)-SMB(113)))+F(87)*(4*
     &   SMB(108)+2*SMB(114)+SMB(115))+F(88)*(SMB(117)+SMB(118))+F(16
     &   32)*(SMB(240)-2*Fa(21)*(SMB(233)+SMB(247)+SMB(248))-12*mup3m
     &   up6*(SMB(102)+SMB(254)+SMB(255))+6*(SMB(69)+SMB(234)+SMB(256
     &   )+SMB(257))+SMB(269)+SMB(270)))/(s23*s56)
       result(4) = (F(89)*SMB(1)+F(90)*SMB(2)+F(91)*SMB(3)+F(93)*SMB(5
     &   )+Fa(28)*SMB(10)+Fa(33)*SMB(36)+Fa(34)*SMB(37)+F(115)*SMB(38
     &   )+Fa(37)*SMB(41)+Fa(38)*SMB(43)+Fa(43)*SMB(66)+F(126)*SMB(68
     &   )+Fa(45)*SMB(82)+F(129)*SMB(91)+32*(Fa(41)*SMB(57)-Fa(42)*SM
     &   B(58)+Fa(50)*SMB(121))-F(131)*SMB(124)+Fa(51)*SMB(129)+F(134
     &   )*SMB(161)+F(99)*SMB(164)+F(99)*SMB(165)+Fa(55)*SMB(168)+F(1
     &   38)*SMB(171)+F(139)*SMB(173)+F(140)*SMB(180)+4*(F(92)*SMB(4)
     &   -Fa(27)*SMB(7)-F(109)*SMB(26)-Fa(36)*SMB(40)-F(130)*SMB(123)
     &   -F(99)*(p2mup3*(SMB(77)+4*(SMB(61)-SMB(80)+SMB(87))+SMB(93))
     &   +Fa(47)*SMB(102)-2*(SMB(127)+Fa(52)*SMB(131)+Fa(44)*SMB(138)
     &   )-SMB(139))+Fa(49)*(-8*F(107)*SMB(120)+F(99)*SMB(143))-p1mup
     &   6*F(99)*SMB(151)-p1mup6*F(99)*SMB(152)+F(99)*SMB(157)+F(99)*
     &   SMB(158)+F(110)*(SMB(59)-SMB(170))+F(94)*(SMB(6)+SMB(11)-SMB
     &   (172))+F(108)*(SMB(60)-SMB(174))+p1mup3*F(119)*SMB(184)-F(14
     &   1)*SMB(187))-8*(Fa(29)*SMB(14)-Fa(30)*SMB(15)-Fa(32)*SMB(22)
     &   -F(111)*SMB(32)-Fa(35)*SMB(39)-Fa(39)*SMB(44)+Fa(40)*SMB(56)
     &   +Fa(48)*SMB(119)+F(99)*Fa(44)*SMB(154)+Fa(53)*SMB(166)+Fa(54
     &   )*SMB(167)-F(111)*SMB(179)+F(128)*SMB(182)+F(99)*SMB(188)+p1
     &   mup3*F(99)*SMB(190))-16*(F(106)*Fa(31)*SMB(16)-p1mup6*F(99)*
     &   SMB(18)+F(111)*Fa(44)*SMB(176)-p1mup3*F(111)*SMB(177)+p2mup6
     &   *F(99)*SMB(183)+p1mup6*F(99)*SMB(198)+p1mup6*F(99)*SMB(199)-
     &   F(99)*Fa(44)*SMB(219))+2*(p1mup3*(4*F(100)*SMB(29)+F(119)*(2
     &   *SMB(42)+SMB(67)))+F(127)*SMB(75)+F(101)*(SMB(13)+SMB(76))+F
     &   (119)*(SMB(46)+SMB(83))+Fa(46)*SMB(101)+p2mup6*(8*F(111)*(SM
     &   B(23)+SMB(86))-2*F(117)*SMB(100)+F(99)*(-2*SMB(45)+SMB(103)+
     &   SMB(104)))-F(118)*SMB(135)-F(133)*SMB(136)-F(99)*Fa(44)*SMB(
     &   163)+s56*F(105)*(SMB(134)+2*(SMB(122)-SMB(169)))+F(98)*(Fa(4
     &   4)*(2*SMB(81)+SMB(90))+2*(2*SMB(181)+SMB(186)))+Fa(56)*SMB(1
     &   89)+F(117)*(-2*SMB(153)+SMB(160)+4*SMB(191))+Fa(21)*(F(117)*
     &   (2*SMB(141)+SMB(142))+8*(F(107)*SMB(17)-F(111)*SMB(185))-F(9
     &   9)*(SMB(144)+SMB(145)-4*SMB(150)-8*SMB(197)))+F(99)*(-4*Fa(2
     &   1)*(SMB(221)+SMB(222))-2*(2*(SMB(228)+SMB(229))+SMB(236)+SMB
     &   (238)+4*p1mup3*SMB(241)+Fa(31)*(SMB(249)+SMB(250))-p2mup6*(2
     &   *(SMB(226)+SMB(227))+p1mup3*(SMB(254)+SMB(255))))+p1mup3*(4*
     &   (SMB(258)+SMB(259))+SMB(264)+SMB(267))))+F(1631)*(-12*SMB(20
     &   8)-48*SMB(209)+SMB(240)-24*(SMB(225)+SMB(244)+SMB(245))+4*Fa
     &   (57)*SMB(246)-2*(Fa(21)*(SMB(233)+6*(SMB(220)+SMB(242)+SMB(2
     &   43))+SMB(247)+SMB(248))+Fa(18)*(6*SMB(253)+SMB(262)))+6*(SMB
     &   (69)+SMB(234)-SMB(235)+SMB(256)+SMB(257)-SMB(263)-SMB(266))+
     &   SMB(269)+SMB(270)))/(s34*s56)
       result(5) = (F(144)*SMB(3)+F(147)*SMB(6)+Fa(58)*SMB(7)+F(157)*S
     &   MB(13)+Fa(59)*SMB(19)+F(182)*SMB(27)+Fa(65)*SMB(40)+Fa(66)*S
     &   MB(41)+F(206)*SMB(68)+F(207)*SMB(69)+F(212)*SMB(76)+8*(Fa(60
     &   )*SMB(20)+Fa(61)*SMB(21)-Fa(62)*SMB(22)+Fa(63)*SMB(23)+Fa(64
     &   )*SMB(24)-Fa(68)*SMB(84))+F(220)*SMB(91)+F(270)*(2*SMB(114)+
     &   SMB(115))+F(271)*SMB(117)+F(272)*SMB(118)+F(1637)*(SMB(240)+
     &   SMB(269)+SMB(270))+F(278)*(2*SMB(274)+SMB(278))+2*(F(143)*SM
     &   B(2)+F(156)*SMB(12)+F(208)*SMB(70)+F(209)*SMB(71)+F(211)*SMB
     &   (75)+F(213)*SMB(79)+Fa(70)*SMB(96)-Fa(71)*SMB(98)-Fa(72)*SMB
     &   (99)+Fa(74)*SMB(101)-Fa(75)*SMB(103)-Fa(76)*SMB(104)+F(268)*
     &   SMB(108)-F(210)*(SMB(74)+SMB(112))-F(269)*SMB(113)+Fa(80)*(S
     &   MB(233)+SMB(247)+SMB(248))+F(247)*(SMB(234)+SMB(256)+SMB(257
     &   ))+Fa(81)*SMB(272)-Fa(82)*SMB(295)+F(230)*(SMB(296)+SMB(297)
     &   +SMB(298)))+4*(F(142)*SMB(1)+F(145)*SMB(4)-F(146)*SMB(5)+F(1
     &   53)*SMB(9)+F(155)*SMB(11)+F(183)*SMB(34)+F(184)*SMB(35)+Fa(6
     &   7)*SMB(45)+F(195)*SMB(46)-F(196)*SMB(47)+F(203)*SMB(51)+F(20
     &   4)*SMB(52)+s45*F(185)*(SMB(38)+SMB(59))+F(171)*SMB(60)+F(175
     &   )*SMB(64)-p6mup3*F(205)*SMB(65)-F(214)*SMB(83)+F(221)*SMB(94
     &   )+Fa(69)*SMB(95)-Fa(73)*SMB(100)+Fa(77)*SMB(105)-Fa(78)*SMB(
     &   106)+Fa(79)*SMB(107)-mup3mup6*(-(F(186)*SMB(39))+2*(F(152)*S
     &   MB(8)+F(154)*SMB(10)+F(195)*SMB(43)-F(196)*SMB(44)+F(201)*SM
     &   B(49)+F(202)*SMB(50)+s45*F(185)*SMB(57)+F(171)*SMB(58)+F(175
     &   )*SMB(63)-F(214)*SMB(82))+F(219)*SMB(88)-F(210)*SMB(110)-F(2
     &   69)*SMB(111)+F(230)*(SMB(97)+SMB(251)+SMB(252))+F(247)*(SMB(
     &   102)+SMB(254)+SMB(255))+F(278)*SMB(299))))/s23
       result(6) = (F(283)*SMB(1)+F(285)*SMB(3)+F(288)*SMB(6)+Fa(83)*S
     &   MB(7)+Fa(86)*SMB(15)+Fa(90)*SMB(22)+F(338)*SMB(26)+Fa(94)*SM
     &   B(36)+Fa(95)*SMB(37)+F(349)*SMB(38)+Fa(96)*SMB(39)+Fa(97)*SM
     &   B(40)+Fa(100)*SMB(43)+Fa(101)*SMB(44)+Fa(108)*SMB(66)+F(425)
     &   *SMB(68)+F(426)*SMB(69)+F(427)*SMB(75)+F(442)*SMB(91)-16*(Fa
     &   (84)*SMB(10)+Fa(87)*SMB(16)+Fa(88)*SMB(17)+Fa(104)*SMB(57)-F
     &   a(105)*SMB(58)-Fa(124)*SMB(120)-Fa(125)*SMB(121))+F(485)*SMB
     &   (134)+Fa(128)*SMB(138)+Fa(135)*SMB(152)+F(495)*SMB(161)+Fa(1
     &   37)*SMB(166)+Fa(139)*SMB(168)+Fa(140)*SMB(176)+F(418)*(2*(SM
     &   B(139)+SMB(157))+SMB(164)-4*(SMB(188)+SMB(228)+SMB(236)))+F(
     &   440)*(2*SMB(158)+SMB(165)-4*(SMB(229)+SMB(238)))+4*(F(286)*S
     &   MB(4)+F(287)*SMB(5)+F(305)*SMB(11)-F(340)*SMB(33)+Fa(98)*SMB
     &   (41)+Fa(102)*SMB(45)+F(375)*SMB(46)-F(399)*SMB(47)+F(401)*SM
     &   B(59)-F(400)*SMB(60)+F(416)*SMB(64)+Fa(111)*SMB(78)-Fa(110)*
     &   (SMB(77)-2*SMB(80))+Fa(113)*SMB(82)-F(439)*SMB(83)+Fa(117)*S
     &   MB(93)+Fa(118)*SMB(100)+Fa(120)*SMB(102)-F(470)*SMB(122)-F(4
     &   71)*SMB(123)-F(472)*SMB(124)+Fa(126)*SMB(129)+Fa(130)*SMB(14
     &   3)+Fa(134)*SMB(151)+F(508)*SMB(169)+F(509)*SMB(170)+F(510)*S
     &   MB(171)+F(511)*SMB(172)+F(512)*SMB(173)+F(346)*SMB(174)+F(51
     &   3)*SMB(179)-F(514)*SMB(180)+F(515)*SMB(181)-F(516)*SMB(182)+
     &   F(377)*SMB(186)+F(522)*SMB(187)-Fa(146)*SMB(190)+F(403)*SMB(
     &   191)+Fa(133)*(SMB(150)+2*SMB(197))-F(523)*SMB(208)-F(524)*SM
     &   B(209)+p6mup3*(F(418)*(2*SMB(125)+SMB(137))+2*(F(440)*SMB(21
     &   9)-F(419)*SMB(241)))+Fa(150)*SMB(246)+Fa(151)*SMB(249)-Fa(15
     &   2)*SMB(250)+Fa(153)*SMB(254)-Fa(154)*SMB(255))+F(1636)*(SMB(
     &   240)-2*Fa(18)*SMB(262)-6*(SMB(235)+SMB(263)+SMB(266))+SMB(26
     &   9)+SMB(270))+F(528)*SMB(277)-8*(Fa(85)*SMB(14)-Fa(89)*SMB(18
     &   )-Fa(91)*SMB(23)+Fa(92)*SMB(24)-Fa(93)*SMB(29)+F(341)*SMB(32
     &   )-Fa(99)*SMB(42)-Fa(103)*SMB(56)+Fa(106)*SMB(61)-Fa(107)*SMB
     &   (62)+Fa(112)*SMB(81)+Fa(114)*SMB(86)-Fa(115)*SMB(87)-Fa(123)
     &   *SMB(119)+Fa(127)*SMB(131)-Fa(138)*SMB(167)-Fa(141)*SMB(177)
     &   -Fa(142)*SMB(183)+Fa(143)*SMB(184)-Fa(144)*SMB(185)-Fa(147)*
     &   SMB(198)-Fa(148)*SMB(199)-F(376)*(SMB(279)+SMB(280)+SMB(281)
     &   ))+2*(F(284)*SMB(2)+F(306)*SMB(13)-Fa(109)*SMB(67)+F(428)*SM
     &   B(76)+Fa(116)*SMB(90)+Fa(119)*SMB(101)-F(486)*SMB(135)-F(487
     &   )*SMB(136)+Fa(129)*(2*SMB(141)+SMB(142))-F(494)*SMB(153)+F(4
     &   23)*SMB(160)-Fa(136)*(2*SMB(154)+SMB(163))+F(429)*(SMB(79)-2
     &   *SMB(178))+Fa(145)*SMB(189)+Fa(131)*(SMB(144)+4*SMB(221))+Fa
     &   (132)*(SMB(145)+4*SMB(222))-Fa(121)*(SMB(103)+4*SMB(226))-Fa
     &   (122)*(SMB(104)+4*SMB(227))+Fa(149)*(SMB(233)+6*(SMB(220)+SM
     &   B(242)+SMB(243))+SMB(247)+SMB(248))-F(527)*(-SMB(234)+4*(SMB
     &   (225)+SMB(244)+SMB(245))+2*Fa(18)*SMB(253)-SMB(256)-SMB(257)
     &   )+Fa(155)*(2*SMB(258)+SMB(264))+Fa(156)*(2*SMB(259)+SMB(267)
     &   )-Fa(157)*(4*SMB(271)+SMB(272))-F(419)*(SMB(273)-2*(SMB(275)
     &   +SMB(283)))-Fa(158)*(4*SMB(282)+SMB(285))-F(376)*(-2*Fa(18)*
     &   SMB(284)+SMB(286)+SMB(287)+SMB(288))))/s34
       result(7) = (F(538)*SMB(13)+F(553)*SMB(27)+Fa(165)*SMB(39)+Fa(1
     &   68)*SMB(43)+Fa(169)*SMB(44)-16*(Fa(160)*SMB(8)+Fa(162)*SMB(1
     &   6)-Fa(170)*(SMB(49)+SMB(50)))+Fa(172)*SMB(66)+F(582)*SMB(75)
     &   +F(583)*SMB(76)+F(593)*SMB(115)+F(594)*SMB(117)+F(610)*SMB(1
     &   61)+F(611)*SMB(164)-8*(Fa(163)*SMB(29)+Fa(164)*SMB(30)-Fa(18
     &   1)*SMB(130)-Fa(176)*(SMB(89)+SMB(175)))+s16*(F(530)*SMB(2)+F
     &   (531)*SMB(3)+Fa(159)*SMB(7)+Fa(161)*SMB(14)+Fa(166)*SMB(40)+
     &   F(579)*SMB(69)+F(603)*SMB(135)+F(533)*(-4*SMB(5)+SMB(68)-2*(
     &   SMB(46)-SMB(47)+SMB(186)))+2*(F(529)*SMB(1)+Fa(173)*SMB(67)+
     &   F(602)*(SMB(134)+SMB(136))+Fa(177)*(SMB(90)+SMB(189))-F(580)
     &   *(SMB(70)+SMB(71)-2*(SMB(193)+SMB(206))))+4*(F(532)*SMB(4)-F
     &   (552)*SMB(25)+Fa(171)*SMB(65)+F(535)*SMB(207)))+p1mup6*(F(55
     &   0)*SMB(19)+F(560)*SMB(41)-8*(F(544)*SMB(15)-F(548)*SMB(17)+F
     &   (549)*SMB(18)+F(551)*SMB(20)-F(535)*SMB(21)-F(588)*SMB(95))+
     &   F(589)*SMB(96)-2*(-4*F(605)*SMB(141)+F(606)*SMB(142)-2*(F(60
     &   9)*SMB(147)+F(608)*(SMB(146)+SMB(148))+2*(-(F(612)*SMB(167))
     &   +F(554)*(SMB(166)+SMB(168))-F(536)*SMB(213)-F(601)*(SMB(210)
     &   +SMB(216)))+F(614)*SMB(220))+F(604)*(SMB(144)+2*(SMB(151)+SM
     &   B(221)))+F(581)*(SMB(98)+2*(SMB(106)+SMB(223)))))-4*(F(544)*
     &   SMB(26)-F(548)*SMB(32)+F(549)*SMB(33)+F(551)*SMB(34)-F(535)*
     &   SMB(35)+Fa(167)*SMB(42)-Fa(174)*(-SMB(77)+SMB(78))+F(588)*SM
     &   B(114)+Fa(180)*SMB(129)+F(605)*SMB(160)+F(554)*(SMB(170)+SMB
     &   (174))+Fa(175)*(SMB(81)+SMB(184))+p6mup3*(F(581)*(SMB(72)-SM
     &   B(73))+F(604)*(SMB(137)-SMB(190))+2*(F(536)*SMB(28)-F(554)*(
     &   SMB(36)+SMB(37))+F(601)*(SMB(132)+SMB(202))))-F(536)*SMB(215
     &   )-F(601)*(SMB(212)+SMB(218))-Fa(182)*(SMB(143)+SMB(250))-Fa(
     &   178)*(SMB(97)+SMB(252)))+F(1638)*(SMB(240)+2*p1mup6*(SMB(233
     &   )+SMB(247)+SMB(248))+SMB(269)+SMB(270))+F(604)*(2*(SMB(236)-
     &   SMB(275))+SMB(277))+F(581)*(2*(SMB(237)-SMB(276))+SMB(278))-
     &   2*(F(537)*SMB(12)-F(609)*SMB(195)-F(608)*(SMB(162)+SMB(196))
     &   -F(615)*SMB(208)-F(614)*(2*SMB(187)-SMB(235))-2*(F(613)*SMB(
     &   172)+Fa(184)*SMB(194)+Fa(185)*SMB(246)+Fa(186)*SMB(249)+Fa(1
     &   87)*SMB(251))-Fa(188)*SMB(262)+F(616)*SMB(263)-Fa(189)*SMB(2
     &   64)-Fa(190)*SMB(265)+F(614)*SMB(266)-Fa(183)*(SMB(163)+SMB(2
     &   67))-Fa(179)*(SMB(116)+SMB(268))-p1mup6*(2*(F(616)*SMB(242)+
     &   F(614)*SMB(243))+F(604)*(2*SMB(282)+SMB(285))+F(581)*(2*SMB(
     &   293)+SMB(295)))))/s16
       result(8) = (F(627)*SMB(1)+F(629)*SMB(3)+F(631)*SMB(5)+Fa(191)*
     &   SMB(7)+F(641)*SMB(9)+F(655)*SMB(12)+Fa(194)*SMB(14)+Fa(195)*
     &   SMB(15)+Fa(197)*SMB(18)+Fa(198)*SMB(19)+Fa(200)*SMB(22)+F(68
     &   2)*SMB(25)+F(683)*SMB(26)+F(684)*SMB(27)+Fa(204)*SMB(36)+Fa(
     &   205)*SMB(37)+F(696)*SMB(38)+Fa(206)*SMB(39)+Fa(207)*SMB(40)+
     &   Fa(208)*SMB(41)+Fa(215)*SMB(50)+F(751)*SMB(51)+F(752)*SMB(52
     &   )+Fa(217)*SMB(56)+Fa(218)*SMB(57)+F(769)*SMB(59)+Fa(221)*SMB
     &   (65)+Fa(222)*SMB(66)+Fa(223)*SMB(67)+F(782)*SMB(70)+F(783)*S
     &   MB(71)+F(784)*SMB(75)+Fa(225)*SMB(81)+F(732)*(SMB(68)+2*(SMB
     &   (46)+SMB(83)))+Fa(229)*SMB(88)+Fa(231)*SMB(90)+F(800)*SMB(91
     &   )+Fa(235)*SMB(101)+F(808)*SMB(108)+F(809)*SMB(115)+F(810)*(S
     &   MB(117)+SMB(118))+Fa(240)*SMB(119)+Fa(241)*SMB(121)+F(817)*S
     &   MB(122)+F(818)*SMB(123)+F(819)*SMB(124)+Fa(243)*SMB(129)+Fa(
     &   245)*SMB(131)+F(828)*SMB(134)+F(829)*SMB(135)+F(830)*SMB(136
     &   )+Fa(248)*SMB(150)+F(832)*SMB(153)+F(833)*SMB(161)+Fa(250)*S
     &   MB(166)+Fa(251)*SMB(167)+Fa(252)*SMB(168)+F(845)*SMB(169)+F(
     &   846)*SMB(170)+F(847)*SMB(171)+F(849)*SMB(173)+F(850)*SMB(180
     &   )+F(636)*SMB(181)+F(851)*SMB(182)+Fa(209)*(SMB(42)+SMB(184))
     &   +F(852)*SMB(186)+F(853)*SMB(187)+Fa(256)*SMB(189)+Fa(257)*SM
     &   B(190)+F(858)*SMB(191)+F(859)*SMB(192)+F(860)*SMB(193)+mup3m
     &   up6*(F(815)*SMB(120)+4*(F(676)*(2*SMB(54)+SMB(110)+SMB(111))
     &   -F(670)*(2*SMB(126)+SMB(155)+SMB(156))+2*(F(827)*SMB(133)-F(
     &   863)*(SMB(201)+SMB(203)))))+F(864)*SMB(206)+F(865)*SMB(207)-
     &   8*(Fa(192)*SMB(8)+Fa(193)*SMB(10)-Fa(196)*SMB(16)-Fa(199)*SM
     &   B(21)-Fa(202)*SMB(29)+Fa(203)*SMB(30)+Fa(210)*SMB(43)-Fa(211
     &   )*SMB(44)-Fa(213)*SMB(48)-Fa(214)*SMB(49)+Fa(219)*SMB(58)-Fa
     &   (216)*(2*SMB(53)+SMB(73))+Fa(226)*SMB(82)-Fa(227)*SMB(84)+Fa
     &   (201)*(SMB(23)+SMB(86))+Fa(220)*(SMB(61)-SMB(80)+SMB(87))-Fa
     &   (230)*SMB(89)-Fa(232)*SMB(92)-Fa(234)*SMB(100)-Fa(242)*SMB(1
     &   28)-Fa(244)*SMB(130)+Fa(253)*SMB(175)-Fa(254)*SMB(177)+Fa(25
     &   5)*SMB(185)-Fa(259)*SMB(197)-Fa(260)*(SMB(198)+SMB(199))-Fa(
     &   261)*(-SMB(200)+SMB(202))+Fa(262)*SMB(211)+Fa(263)*(SMB(214)
     &   +SMB(217)))+Fa(21)*(F(802)*SMB(96)+8*(F(669)*SMB(17)-F(675)*
     &   SMB(20)-F(801)*SMB(95)-F(804)*SMB(141))+F(831)*SMB(142)-2*(2
     &   *(F(790)*SMB(146)+F(671)*(SMB(147)+SMB(148))+2*(F(862)*SMB(2
     &   10)+F(639)*(SMB(213)+SMB(216)))+F(789)*SMB(220))+F(670)*(SMB
     &   (144)+SMB(145)+2*(SMB(221)+SMB(222)))-F(676)*(SMB(98)+SMB(99
     &   )+2*(SMB(223)+SMB(224)))))+F(670)*(8*SMB(127)+4*(SMB(139)+SM
     &   B(157)+SMB(158))+SMB(164)+SMB(165)-2*(2*(SMB(188)+SMB(228)+S
     &   MB(229))+SMB(236)+SMB(238)))-16*Fa(264)*SMB(241)+4*(F(628)*S
     &   MB(2)+F(630)*SMB(4)+F(632)*(SMB(6)+SMB(11))+F(680)*SMB(60)+F
     &   a(224)*(SMB(77)+SMB(93))-F(687)*(SMB(34)+SMB(94))+Fa(233)*SM
     &   B(97)-Fa(236)*SMB(102)-Fa(238)*SMB(105)-Fa(239)*(SMB(106)+SM
     &   B(107))-F(801)*SMB(114)+Fa(247)*SMB(143)+Fa(249)*(SMB(151)+S
     &   MB(152))-F(804)*SMB(160)+F(848)*SMB(172)+F(843)*SMB(174)+F(6
     &   81)*(SMB(32)+SMB(179))+Fa(212)*(SMB(45)+2*SMB(183))+Fa(258)*
     &   SMB(194)-F(729)*(SMB(208)+4*SMB(209))-F(862)*SMB(212)-F(639)
     &   *(SMB(215)+SMB(218))+Fa(265)*SMB(246))+Fa(266)*(SMB(249)+SMB
     &   (250))+Fa(267)*(SMB(251)+SMB(252))+2*(F(656)*(SMB(13)+SMB(76
     &   ))+Fa(228)*(-8*SMB(85)+4*SMB(109)+SMB(116))+F(790)*(2*(SMB(1
     &   40)+SMB(159))+SMB(162))+F(671)*(SMB(195)+SMB(196)+2*(SMB(31)
     &   +SMB(204)+SMB(205)))+Fa(246)*(4*(SMB(138)-SMB(154))-SMB(163)
     &   -8*(SMB(176)-SMB(219)))+Fa(237)*(SMB(103)+SMB(104)+2*(SMB(22
     &   6)+SMB(227)))-F(789)*(4*SMB(225)+SMB(235))+F(676)*(-4*SMB(55
     &   )-2*(SMB(74)+SMB(112)+SMB(113)-SMB(230)-SMB(231)-SMB(232))+S
     &   MB(237)+SMB(239))-2*Fa(269)*(SMB(254)+SMB(255))+Fa(268)*(6*S
     &   MB(253)+SMB(262))-F(729)*(2*Fa(21)*(SMB(242)+SMB(243))+4*(SM
     &   B(244)+SMB(245))+SMB(263)+SMB(266))+Fa(270)*(4*(SMB(258)+SMB
     &   (259))+SMB(264)+SMB(267))-Fa(271)*(4*(SMB(260)+SMB(261))+SMB
     &   (265)+SMB(268)))+F(1635)*(SMB(240)-2*Fa(21)*(SMB(233)+SMB(24
     &   7)+SMB(248))+6*(SMB(69)+SMB(234)+SMB(256)+SMB(257))+SMB(269)
     &   +SMB(270)))/s56
       result(9) = F(874)*SMB(3)+Fa(272)*SMB(7)+F(918)*SMB(13)+Fa(289)
     &   *SMB(39)+Fa(290)*SMB(40)+Fa(293)*SMB(43)+Fa(294)*SMB(44)+16*
     &   (Fa(273)*SMB(8)+Fa(274)*SMB(10)+Fa(277)*SMB(16)+Fa(297)*SMB(
     &   49)-Fa(298)*SMB(50))+Fa(301)*SMB(65)+Fa(302)*SMB(66)+F(1155)
     &   *SMB(68)+F(1156)*SMB(69)+F(1158)*SMB(71)+Fa(309)*SMB(82)+F(1
     &   222)*SMB(91)+Fa(315)*SMB(93)+Fa(317)*SMB(96)+Fa(320)*SMB(99)
     &   +Fa(322)*SMB(101)+F(1313)*SMB(108)+F(1316)*SMB(115)+F(1317)*
     &   SMB(117)+F(1318)*SMB(118)+Fa(330)*SMB(129)+F(1364)*SMB(136)+
     &   Fa(335)*SMB(142)+Fa(338)*SMB(145)+F(1424)*SMB(153)+F(1428)*S
     &   MB(161)+F(1369)*SMB(164)+F(1426)*SMB(165)+Fa(363)*SMB(222)+F
     &   a(364)*SMB(224)+F(1643)*SMB(240)+Fa(367)*SMB(243)+F(1519)*SM
     &   B(244)+Fa(373)*SMB(253)+F(1643)*SMB(269)+F(1643)*SMB(270)+F(
     &   1619)*SMB(277)+F(1621)*SMB(278)+8*(Fa(275)*SMB(14)+Fa(276)*S
     &   MB(15)-Fa(278)*SMB(17)+Fa(279)*SMB(18)+Fa(280)*SMB(19)-Fa(28
     &   1)*SMB(20)+Fa(283)*SMB(22)+Fa(284)*SMB(23)-Fa(285)*SMB(24)+F
     &   a(286)*SMB(29)-Fa(287)*SMB(30)+Fa(292)*SMB(42)-Fa(296)*SMB(4
     &   8)+Fa(288)*(SMB(37)-SMB(56))+Fa(299)*SMB(61)-Fa(300)*SMB(62)
     &   +Fa(307)*SMB(80)-Fa(308)*SMB(81)-Fa(310)*SMB(84)-Fa(311)*SMB
     &   (88)+Fa(312)*SMB(89)+Fa(316)*SMB(95)+Fa(321)*SMB(100)-Fa(329
     &   )*SMB(128)-Fa(331)*SMB(130)+Fa(332)*SMB(131)+Fa(334)*SMB(141
     &   )-Fa(346)*SMB(166)-Fa(347)*SMB(167)-Fa(348)*SMB(168)-Fa(349)
     &   *SMB(175)+F(1447)*SMB(180)-Fa(350)*SMB(183)+Fa(351)*SMB(184)
     &   +F(1465)*SMB(192)+Fa(355)*SMB(200)+Fa(356)*SMB(202)+F(1271)*
     &   (SMB(123)+SMB(181)-SMB(209))+Fa(357)*SMB(210)+Fa(358)*SMB(21
     &   1)-Fa(359)*SMB(213)+Fa(282)*(SMB(21)-SMB(214))-Fa(360)*SMB(2
     &   16)-Fa(361)*SMB(217)-F(1520)*SMB(245)-F(1506)*SMB(289))+2*(F
     &   (917)*SMB(12)-Fa(303)*SMB(67)+F(1157)*SMB(70)+F(1164)*SMB(75
     &   )+F(1165)*SMB(76)+F(1172)*SMB(79)-Fa(313)*SMB(90)+Fa(319)*SM
     &   B(98)-Fa(324)*SMB(103)-Fa(325)*SMB(104)-Fa(314)*SMB(116)-F(1
     &   362)*SMB(134)-F(1363)*SMB(135)+Fa(337)*SMB(144)+F(1429)*SMB(
     &   162)+Fa(345)*SMB(163)+F(1453)*SMB(186)+F(1454)*SMB(187)-Fa(3
     &   52)*SMB(189)+F(1472)*SMB(195)+F(1473)*SMB(196)-F(1480)*SMB(2
     &   08)-F(1508)*SMB(235)+F(1369)*(SMB(139)+SMB(157)-SMB(236))-F(
     &   1163)*(SMB(74)+SMB(112)-SMB(237))+F(1426)*(SMB(158)-SMB(238)
     &   )+F(1314)*(-SMB(113)+SMB(239))+Fa(365)*(SMB(233)+SMB(247)+SM
     &   B(248))+Fa(380)*SMB(262)-F(1613)*SMB(263)+Fa(376)*SMB(264)-F
     &   a(378)*SMB(265)+F(1507)*(SMB(234)+SMB(256)+SMB(257)-SMB(266)
     &   )+Fa(377)*SMB(267)-Fa(379)*SMB(268)-Fa(381)*SMB(272)+F(1619)
     &   *(SMB(273)-SMB(275))+F(1620)*(-SMB(274)+SMB(276))+Fa(382)*SM
     &   B(285)-F(1506)*(SMB(286)+SMB(287)+SMB(288))-Fa(384)*SMB(295)
     &   +F(1270)*(SMB(296)+SMB(297)+SMB(298)))+4*(-(F(872)*SMB(1))+F
     &   (873)*SMB(2)-F(875)*SMB(4)+F(876)*SMB(5)-F(877)*SMB(6)+F(892
     &   )*SMB(9)-F(916)*SMB(11)+F(985)*SMB(25)+F(986)*SMB(26)+F(987)
     &   *SMB(27)-F(995)*SMB(32)-F(983)*SMB(33)+F(996)*SMB(34)+Fa(291
     &   )*SMB(41)-Fa(295)*SMB(45)+F(1079)*SMB(46)-F(1080)*SMB(47)+F(
     &   1129)*SMB(51)+F(1130)*SMB(52)+Fa(304)*SMB(73)+Fa(305)*SMB(77
     &   )-Fa(306)*SMB(78)-F(1194)*SMB(83)-F(1232)*SMB(94)+Fa(318)*SM
     &   B(97)-Fa(323)*SMB(102)+Fa(326)*SMB(105)-Fa(327)*SMB(106)-Fa(
     &   328)*SMB(107)+Fa(314)*(SMB(92)-SMB(109))-F(1315)*SMB(114)-Fa
     &   (333)*SMB(138)+Fa(336)*SMB(143)+Fa(339)*SMB(146)-Fa(340)*SMB
     &   (147)-Fa(341)*SMB(148)+F(1410)*SMB(149)+Fa(342)*SMB(150)+Fa(
     &   343)*SMB(151)+Fa(344)*SMB(152)+Fa(345)*SMB(154)-F(1427)*SMB(
     &   160)-F(1002)*(SMB(38)+SMB(59)-SMB(170))+F(1445)*SMB(172)+F(1
     &   131)*(SMB(60)-SMB(174))+F(984)*(SMB(64)-SMB(178))+Fa(353)*SM
     &   B(190)-F(1466)*SMB(193)+Fa(354)*SMB(194)-F(1478)*SMB(206)-F(
     &   1479)*SMB(207)+F(1370)*(SMB(140)+SMB(159)-SMB(212))-F(993)*(
     &   SMB(31)+SMB(204)-SMB(215))+F(1446)*(-SMB(205)+SMB(218))+Fa(3
     &   62)*SMB(220)+Fa(337)*SMB(221)+Fa(319)*SMB(223)-F(1504)*SMB(2
     &   25)-Fa(324)*SMB(226)-Fa(325)*SMB(227)+Fa(366)*SMB(242)+Fa(36
     &   8)*SMB(246)+Fa(369)*SMB(249)+Fa(370)*SMB(250)+Fa(371)*SMB(25
     &   1)+Fa(372)*SMB(252)+Fa(374)*SMB(254)-Fa(375)*SMB(255)+Fa(376
     &   )*SMB(258)+Fa(377)*SMB(259)-Fa(378)*SMB(260)-Fa(379)*SMB(261
     &   )-Fa(381)*SMB(271)+F(1506)*(SMB(279)+SMB(280)+SMB(281))+Fa(3
     &   82)*SMB(282)+Fa(383)*SMB(284)-F(1270)*(SMB(290)+SMB(291)+SMB
     &   (292))-Fa(384)*SMB(293)+Fa(385)*SMB(294)+F(1646)*SMB(300))
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
       resultborn = -((SMB(13)+2*SMB(27)+2*p1mup6*(SMB(7)+2*SMB(19)+SM
     &   B(40)+SMB(41))+SMB(75)+SMB(76))/(s16*s23*s234))-(2*(p6sq-s12
     &   3-s16+s45)*SMB(4)+2*SMB(6)+2*SMB(11)+SMB(13)+4*(p1mup6+p4mup
     &   6)*SMB(22)+2*SMB(27)+2*SMB(38)-2*(p5mup6+p6mup6)*(SMB(7)+2*S
     &   MB(19)-2*SMB(22)+SMB(40)+SMB(41))-4*mup3mup6*(SMB(10)+SMB(57
     &   )+SMB(58))+2*SMB(59)+2*SMB(60)+SMB(75)+SMB(76))/(s23*s234*s5
     &   6)-(2*(p6sq-s12-2*s16+s345)*SMB(4)+2*SMB(6)+2*SMB(11)+SMB(13
     &   )+4*(p1mup6-p2mup6)*SMB(22)-2*SMB(26)+4*(p1mup3+p5mup3+p6mup
     &   3)*SMB(37)+2*SMB(38)-2*(p5mup6+p6mup6)*(SMB(7)-2*SMB(15)-2*S
     &   MB(22)-2*(p1mup3+p5mup3+p6mup3)*SMB(39)+SMB(40)+SMB(41))-4*(
     &   p1mup3+p5mup3+p6mup3)*SMB(56)+2*SMB(59)+2*SMB(60)-2*(p1mup3+
     &   p5mup3+p6mup3)*SMB(66)+SMB(75)+SMB(76)-4*(p5mup6+p6mup6)*SMB
     &   (166)-4*(p5mup6+p6mup6)*SMB(167)-4*(p5mup6+p6mup6)*SMB(168)-
     &   2*(SMB(170)+SMB(172)+SMB(174)))/(s234*s34*s56)-(-2*s16*SMB(4
     &   )+SMB(13)-2*SMB(26)+4*p6mup3*(SMB(14)+SMB(36)+SMB(37))+SMB(7
     &   5)+SMB(76)-2*((p1mup3+p5mup3+p6mup3)*SMB(66)-p1mup6*(SMB(7)-
     &   2*SMB(15)-2*(p1mup3+p5mup3+p6mup3)*SMB(39)+SMB(40)+SMB(41)+2
     &   *(SMB(166)+SMB(167)+SMB(168))))-2*(SMB(170)+SMB(172)+SMB(174
     &   )))/(s16*s234*s34)
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
