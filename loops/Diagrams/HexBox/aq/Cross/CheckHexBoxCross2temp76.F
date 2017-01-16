c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:9/2/2010
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
       subroutine CheckHexBoxCross276(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_
     &   p5,barpsi_p2,psi_p4,mup3,mup6,musq,comp,ngluon,posgluon,resu
     &   lt,resultn,resultborn)
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

       Complex*16 SMB0P(32),SMB1P1(0:3),SMB1P2(0:3),SMB1P3(0:3),SMB1P4
     &   (0:3),SMB1P5(0:3),SMB1P6(0:3),SMB1P7(0:3),SMB1P8(0:3),SMB1P9
     &   (0:3),SMB1P10(0:3),SMB1P11(0:3),SMB1P12(0:3),SMB1P13(0:3),SM
     &   B1P14(0:3),SMB1P15(0:3),SMB1P16(0:3),SMB1P17(0:3),SMB1P18(0:
     &   3),SMB1P19(0:3),SMB1P20(0:3),SMB1P21(0:3),SMB1P22(0:3),SMB1P
     &   23(0:3),SMB1P24(0:3),SMB1P25(0:3),SMB1P26(0:3),SMB1P27(0:3),
     &   SMB1P28(0:3),SMB1P29(0:3),SMB1P30(0:3),SMB2P1(0:3,0:3),SMB2P
     &   2(0:3,0:3),SMB2P3(0:3,0:3),SMB2P4(0:3,0:3),SMB2P5(0:3,0:3),S
     &   MB2P6(0:3,0:3),SMB2P7(0:3,0:3),SMB2P8(0:3,0:3),SMB2P9(0:3,0:
     &   3),SMB2P10(0:3,0:3),SMB2P11(0:3,0:3),SMB2P12(0:3,0:3),SMB2P1
     &   3(0:3,0:3),SMB2P14(0:3,0:3),SMB2P15(0:3,0:3),SMB2P16(0:3,0:3
     &   ),SMB2P17(0:3,0:3),SMB2P18(0:3,0:3),SMB2P19(0:3,0:3),SMB2P20
     &   (0:3,0:3),SMB2P21(0:3,0:3),SMB2P22(0:3,0:3),SMB3P1(0:3,0:3,0
     &   :3),SMB3P2(0:3,0:3,0:3),SMB3P3(0:3,0:3,0:3),SMB3P4(0:3,0:3,0
     &   :3),SMB3P5(0:3,0:3,0:3),SMB3P6(0:3,0:3,0:3),SMB3P7(0:3,0:3,0
     &   :3),SMB3P8(0:3,0:3,0:3)
       Complex*16 SMB0M(32),SMB1M1(0:3),SMB1M2(0:3),SMB1M3(0:3),SMB1M4
     &   (0:3),SMB1M5(0:3),SMB1M6(0:3),SMB1M7(0:3),SMB1M8(0:3),SMB1M9
     &   (0:3),SMB1M10(0:3),SMB1M11(0:3),SMB1M12(0:3),SMB1M13(0:3),SM
     &   B1M14(0:3),SMB1M15(0:3),SMB1M16(0:3),SMB1M17(0:3),SMB1M18(0:
     &   3),SMB1M19(0:3),SMB1M20(0:3),SMB1M21(0:3),SMB1M22(0:3),SMB1M
     &   23(0:3),SMB1M24(0:3),SMB1M25(0:3),SMB1M26(0:3),SMB1M27(0:3),
     &   SMB1M28(0:3),SMB1M29(0:3),SMB1M30(0:3),SMB2M1(0:3,0:3),SMB2M
     &   2(0:3,0:3),SMB2M3(0:3,0:3),SMB2M4(0:3,0:3),SMB2M5(0:3,0:3),S
     &   MB2M6(0:3,0:3),SMB2M7(0:3,0:3),SMB2M8(0:3,0:3),SMB2M9(0:3,0:
     &   3),SMB2M10(0:3,0:3),SMB2M11(0:3,0:3),SMB2M12(0:3,0:3),SMB2M1
     &   3(0:3,0:3),SMB2M14(0:3,0:3),SMB2M15(0:3,0:3),SMB2M16(0:3,0:3
     &   ),SMB2M17(0:3,0:3),SMB2M18(0:3,0:3),SMB2M19(0:3,0:3),SMB2M20
     &   (0:3,0:3),SMB2M21(0:3,0:3),SMB2M22(0:3,0:3),SMB3M1(0:3,0:3,0
     &   :3),SMB3M2(0:3,0:3,0:3),SMB3M3(0:3,0:3,0:3),SMB3M4(0:3,0:3,0
     &   :3),SMB3M5(0:3,0:3,0:3),SMB3M6(0:3,0:3,0:3),SMB3M7(0:3,0:3,0
     &   :3),SMB3M8(0:3,0:3,0:3)

      Complex*16 psi_p1_P(2),barpsi_p5_P(2),barpsi_p2_P(2),psi_p4_P(2) 
       Complex*16 psi_p1_M(2),barpsi_p5_M(2),barpsi_p2_M(2),psi_p4_M(2)
       Complex*16 psi_p1(4),barpsi_p5(4),barpsi_p2(4),psi_p4(4)
       Real*8 P(1355)
       Complex*16 SC1c,SC1r, SC3ccc,SC3ccr,SC3crc,SC3crr,SC3rrc,SC3rrr
     &   , SC5ccccc,SC5cccrc,SC5ccrrc,SC5ccrrr,SC5crccc,SC5crcrc,SC5c
     &   rrrc,SC5rrccc,SC5rrcrc, SC7ccrrccc,SC7ccrrcrc
       EXTERNAL SC1c,SC1r, SC3ccc,SC3ccr,SC3crc,SC3crr,SC3rrc,SC3rrr, 
     &   SC5ccccc,SC5cccrc,SC5ccrrc,SC5ccrrr,SC5crccc,SC5crcrc,SC5crr
     &   rc,SC5rrccc,SC5rrcrc, SC7ccrrccc,SC7ccrrcrc
      Real*8 delta
       External  delta
       Complex*16 ten2cc,ten3cc
       External  ten2cc,ten3cc
       Complex*16 v1(0:3),v2(0:3),v3(0:3) 
       Integer i,j,k,comp
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
       call  HexBoxCross2Fa1(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(1))
       call  HexBoxCross2Fa2(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(20))
       call  HexBoxCross2Fa3(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(39))
       call  HexBoxCross2Fa4(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(58))
       call  HexBoxCross2Fa5(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(77))
       call  HexBoxCross2Fa6(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(96))
       call  HexBoxCross2Fa7(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(115))
       call  HexBoxCross2Fa8(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(134))
       call  HexBoxCross2Fa9(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup6
     &   ,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(153))
       call  HexBoxCross2Fa10(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(172)
     &   )
       call  HexBoxCross2Fa11(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(191)
     &   )
       call  HexBoxCross2Fa12(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(210)
     &   )
       call  HexBoxCross2Fa13(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(229)
     &   )
       call  HexBoxCross2Fa14(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(248)
     &   )
       call  HexBoxCross2Fa15(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(267)
     &   )
       call  HexBoxCross2Fa16(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(286)
     &   )
       call  HexBoxCross2Fa17(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(305)
     &   )
       call  HexBoxCross2Fa18(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(324)
     &   )
       call  HexBoxCross2Fa19(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(343)
     &   )
       call  HexBoxCross2Fa20(p1mup3,p1mup6,p2mup3,p2mup6,p3mup3,p3mup
     &   6,p4mup3,p4mup6,p5mup3,p5mup6,p6mup3,p6mup6,mup3mup6,Fa(362)
     &   )
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
       SMB0P(1) = SC1c(barpsi_p2_P,mup3,psi_p4_P,+1)
       SMB0P(2) = SC1r(barpsi_p5_P,p4,psi_p1_P,+1)
       SMB0P(3) = SC1r(barpsi_p5_P,p6,psi_p1_P,+1)
       SMB0P(4) = SC1r(barpsi_p2_P,p1,psi_p4_P,+1)
       SMB0P(5) = SC3ccr(barpsi_p2_P,mup6,mup3,p1,psi_p4_P,+1)
       SMB0P(6) = SC1r(barpsi_p5_P,p2,psi_p1_P,+1)
       SMB0P(7) = SC1c(barpsi_p2_P,mup6,psi_p4_P,+1)
       SMB0P(8) = SC1r(barpsi_p2_P,p5,psi_p4_P,+1)
       SMB0P(9) = SC1r(barpsi_p2_P,p6,psi_p4_P,+1)
       SMB0P(10) = SC3ccr(barpsi_p2_P,mup6,mup3,p5,psi_p4_P,+1)
       SMB0P(11) = SC3ccr(barpsi_p2_P,mup6,mup3,p6,psi_p4_P,+1)
       SMB0P(12) = SC1c(barpsi_p5_P,mup6,psi_p1_P,+1)
       SMB0P(13) = SC1c(barpsi_p5_P,mup3,psi_p1_P,+1)
       SMB0P(14) = SC3crr(barpsi_p2_P,mup3,p5,p1,psi_p4_P,+1)
       SMB0P(15) = SC3crr(barpsi_p2_P,mup3,p6,p5,psi_p4_P,+1)
       SMB0P(16) = SC3ccr(barpsi_p5_P,mup6,mup3,p4,psi_p1_P,+1)
       SMB0P(17) = SC3crr(barpsi_p2_P,mup6,p5,p1,psi_p4_P,+1)
       SMB0P(18) = SC3crr(barpsi_p2_P,mup6,p6,p5,psi_p4_P,+1)
       SMB0P(19) = SC3ccr(barpsi_p5_P,mup6,mup3,p6,psi_p1_P,+1)
       SMB0P(20) = SC3ccr(barpsi_p5_P,mup6,mup3,p2,psi_p1_P,+1)
       SMB0P(21) = SC3crr(barpsi_p2_P,mup3,p6,p1,psi_p4_P,+1)
       SMB0P(22) = SC3crr(barpsi_p5_P,mup6,p6,p4,psi_p1_P,+1)
       SMB0P(23) = SC3crr(barpsi_p5_P,mup6,p4,p2,psi_p1_P,+1)
       SMB0P(24) = SC3crr(barpsi_p2_P,mup6,p6,p1,psi_p4_P,+1)
       SMB0P(25) = SC3rrr(barpsi_p2_P,p6,p5,p1,psi_p4_P,+1)
       SMB0P(26) = SC3crr(barpsi_p5_P,mup3,p6,p2,psi_p1_P,+1)
       SMB0P(27) = SC3crr(barpsi_p5_P,mup6,p6,p2,psi_p1_P,+1)
       SMB0P(28) = SC3crr(barpsi_p5_P,mup3,p6,p4,psi_p1_P,+1)
       SMB0P(29) = SC3rrr(barpsi_p5_P,p6,p4,p2,psi_p1_P,+1)
       SMB0P(30) = SC3crr(barpsi_p5_P,mup3,p4,p2,psi_p1_P,+1)
       SMB0P(31) = SC5ccrrr(barpsi_p2_P,mup6,mup3,p6,p5,p1,psi_p4_P,+1)
       SMB0P(32) = SC5ccrrr(barpsi_p5_P,mup6,mup3,p6,p4,p2,psi_p1_P,+1)
       SMB0M(1) = SC1c(barpsi_p2_M,mup3,psi_p4_M,-1)
       SMB0M(2) = SC1r(barpsi_p5_M,p4,psi_p1_M,-1)
       SMB0M(3) = SC1r(barpsi_p5_M,p6,psi_p1_M,-1)
       SMB0M(4) = SC1r(barpsi_p2_M,p1,psi_p4_M,-1)
       SMB0M(5) = SC3ccr(barpsi_p2_M,mup6,mup3,p1,psi_p4_M,-1)
       SMB0M(6) = SC1r(barpsi_p5_M,p2,psi_p1_M,-1)
       SMB0M(7) = SC1c(barpsi_p2_M,mup6,psi_p4_M,-1)
       SMB0M(8) = SC1r(barpsi_p2_M,p5,psi_p4_M,-1)
       SMB0M(9) = SC1r(barpsi_p2_M,p6,psi_p4_M,-1)
       SMB0M(10) = SC3ccr(barpsi_p2_M,mup6,mup3,p5,psi_p4_M,-1)
       SMB0M(11) = SC3ccr(barpsi_p2_M,mup6,mup3,p6,psi_p4_M,-1)
       SMB0M(12) = SC1c(barpsi_p5_M,mup6,psi_p1_M,-1)
       SMB0M(13) = SC1c(barpsi_p5_M,mup3,psi_p1_M,-1)
       SMB0M(14) = SC3crr(barpsi_p2_M,mup3,p5,p1,psi_p4_M,-1)
       SMB0M(15) = SC3crr(barpsi_p2_M,mup3,p6,p5,psi_p4_M,-1)
       SMB0M(16) = SC3ccr(barpsi_p5_M,mup6,mup3,p4,psi_p1_M,-1)
       SMB0M(17) = SC3crr(barpsi_p2_M,mup6,p5,p1,psi_p4_M,-1)
       SMB0M(18) = SC3crr(barpsi_p2_M,mup6,p6,p5,psi_p4_M,-1)
       SMB0M(19) = SC3ccr(barpsi_p5_M,mup6,mup3,p6,psi_p1_M,-1)
       SMB0M(20) = SC3ccr(barpsi_p5_M,mup6,mup3,p2,psi_p1_M,-1)
       SMB0M(21) = SC3crr(barpsi_p2_M,mup3,p6,p1,psi_p4_M,-1)
       SMB0M(22) = SC3crr(barpsi_p5_M,mup6,p6,p4,psi_p1_M,-1)
       SMB0M(23) = SC3crr(barpsi_p5_M,mup6,p4,p2,psi_p1_M,-1)
       SMB0M(24) = SC3crr(barpsi_p2_M,mup6,p6,p1,psi_p4_M,-1)
       SMB0M(25) = SC3rrr(barpsi_p2_M,p6,p5,p1,psi_p4_M,-1)
       SMB0M(26) = SC3crr(barpsi_p5_M,mup3,p6,p2,psi_p1_M,-1)
       SMB0M(27) = SC3crr(barpsi_p5_M,mup6,p6,p2,psi_p1_M,-1)
       SMB0M(28) = SC3crr(barpsi_p5_M,mup3,p6,p4,psi_p1_M,-1)
       SMB0M(29) = SC3rrr(barpsi_p5_M,p6,p4,p2,psi_p1_M,-1)
       SMB0M(30) = SC3crr(barpsi_p5_M,mup3,p4,p2,psi_p1_M,-1)
       SMB0M(31) = SC5ccrrr(barpsi_p2_M,mup6,mup3,p6,p5,p1,psi_p4_M,-1)
       SMB0M(32) = SC5ccrrr(barpsi_p5_M,mup6,mup3,p6,p4,p2,psi_p1_M,-1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB1P1(i) = SC1c(barpsi_p2_P,v1,psi_p4_P,+1)
       SMB1P2(i) = SC1c(barpsi_p5_P,v1,psi_p1_P,+1)
       SMB1P3(i) = SC3ccc(barpsi_p2_P,mup6,mup3,v1,psi_p4_P,+1)
       SMB1P4(i) = SC3ccc(barpsi_p5_P,mup6,mup3,v1,psi_p1_P,+1)
       SMB1P5(i) = SC3rrc(barpsi_p2_P,p5,p1,v1,psi_p4_P,+1)
       SMB1P6(i) = SC3rrc(barpsi_p2_P,p6,p5,v1,psi_p4_P,+1)
       SMB1P7(i) = SC3crc(barpsi_p5_P,mup6,p4,v1,psi_p1_P,+1)
       SMB1P8(i) = SC3crc(barpsi_p2_P,mup3,p1,v1,psi_p4_P,+1)
       SMB1P9(i) = SC3crc(barpsi_p2_P,mup3,p5,v1,psi_p4_P,+1)
       SMB1P10(i) = SC3crc(barpsi_p2_P,mup3,p6,v1,psi_p4_P,+1)
       SMB1P11(i) = SC3crc(barpsi_p5_P,mup6,p6,v1,psi_p1_P,+1)
       SMB1P12(i) = SC3crc(barpsi_p5_P,mup6,p2,v1,psi_p1_P,+1)
       SMB1P13(i) = SC3rrc(barpsi_p2_P,p6,p1,v1,psi_p4_P,+1)
       SMB1P14(i) = SC3crc(barpsi_p5_P,mup3,p2,v1,psi_p1_P,+1)
       SMB1P15(i) = SC3crc(barpsi_p2_P,mup6,p1,v1,psi_p4_P,+1)
       SMB1P16(i) = SC3crc(barpsi_p2_P,mup6,p5,v1,psi_p4_P,+1)
       SMB1P17(i) = SC3crc(barpsi_p2_P,mup6,p6,v1,psi_p4_P,+1)
       SMB1P18(i) = SC3rrc(barpsi_p5_P,p4,p2,v1,psi_p1_P,+1)
       SMB1P19(i) = SC3rrc(barpsi_p5_P,p6,p2,v1,psi_p1_P,+1)
       SMB1P20(i) = SC3rrc(barpsi_p5_P,p6,p4,v1,psi_p1_P,+1)
       SMB1P21(i) = SC5ccrrc(barpsi_p2_P,mup6,mup3,p5,p1,v1,psi_p4_P,+1)
       SMB1P22(i) = SC5ccrrc(barpsi_p2_P,mup6,mup3,p6,p1,v1,psi_p4_P,+1)
       SMB1P23(i) = SC5ccrrc(barpsi_p2_P,mup6,mup3,p6,p5,v1,psi_p4_P,+1)
       SMB1P24(i) = SC3crc(barpsi_p5_P,mup3,p4,v1,psi_p1_P,+1)
       SMB1P25(i) = SC3crc(barpsi_p5_P,mup3,p6,v1,psi_p1_P,+1)
       SMB1P26(i) = SC5ccrrc(barpsi_p5_P,mup6,mup3,p6,p2,v1,psi_p1_P,+1)
       SMB1P27(i) = SC5ccrrc(barpsi_p5_P,mup6,mup3,p6,p4,v1,psi_p1_P,+1)
       SMB1P28(i) = SC5crrrc(barpsi_p2_P,mup3,p6,p5,p1,v1,psi_p4_P,+1)
       SMB1P29(i) = SC5ccrrc(barpsi_p5_P,mup6,mup3,p4,p2,v1,psi_p1_P,+1)
       SMB1P30(i) = SC5crrrc(barpsi_p5_P,mup6,p6,p4,p2,v1,psi_p1_P,+1)
       SMB1M1(i) = SC1c(barpsi_p2_M,v1,psi_p4_M,-1)
       SMB1M2(i) = SC1c(barpsi_p5_M,v1,psi_p1_M,-1)
       SMB1M3(i) = SC3ccc(barpsi_p2_M,mup6,mup3,v1,psi_p4_M,-1)
       SMB1M4(i) = SC3ccc(barpsi_p5_M,mup6,mup3,v1,psi_p1_M,-1)
       SMB1M5(i) = SC3rrc(barpsi_p2_M,p5,p1,v1,psi_p4_M,-1)
       SMB1M6(i) = SC3rrc(barpsi_p2_M,p6,p5,v1,psi_p4_M,-1)
       SMB1M7(i) = SC3crc(barpsi_p5_M,mup6,p4,v1,psi_p1_M,-1)
       SMB1M8(i) = SC3crc(barpsi_p2_M,mup3,p1,v1,psi_p4_M,-1)
       SMB1M9(i) = SC3crc(barpsi_p2_M,mup3,p5,v1,psi_p4_M,-1)
       SMB1M10(i) = SC3crc(barpsi_p2_M,mup3,p6,v1,psi_p4_M,-1)
       SMB1M11(i) = SC3crc(barpsi_p5_M,mup6,p6,v1,psi_p1_M,-1)
       SMB1M12(i) = SC3crc(barpsi_p5_M,mup6,p2,v1,psi_p1_M,-1)
       SMB1M13(i) = SC3rrc(barpsi_p2_M,p6,p1,v1,psi_p4_M,-1)
       SMB1M14(i) = SC3crc(barpsi_p5_M,mup3,p2,v1,psi_p1_M,-1)
       SMB1M15(i) = SC3crc(barpsi_p2_M,mup6,p1,v1,psi_p4_M,-1)
       SMB1M16(i) = SC3crc(barpsi_p2_M,mup6,p5,v1,psi_p4_M,-1)
       SMB1M17(i) = SC3crc(barpsi_p2_M,mup6,p6,v1,psi_p4_M,-1)
       SMB1M18(i) = SC3rrc(barpsi_p5_M,p4,p2,v1,psi_p1_M,-1)
       SMB1M19(i) = SC3rrc(barpsi_p5_M,p6,p2,v1,psi_p1_M,-1)
       SMB1M20(i) = SC3rrc(barpsi_p5_M,p6,p4,v1,psi_p1_M,-1)
       SMB1M21(i) = SC5ccrrc(barpsi_p2_M,mup6,mup3,p5,p1,v1,psi_p4_M,-1)
       SMB1M22(i) = SC5ccrrc(barpsi_p2_M,mup6,mup3,p6,p1,v1,psi_p4_M,-1)
       SMB1M23(i) = SC5ccrrc(barpsi_p2_M,mup6,mup3,p6,p5,v1,psi_p4_M,-1)
       SMB1M24(i) = SC3crc(barpsi_p5_M,mup3,p4,v1,psi_p1_M,-1)
       SMB1M25(i) = SC3crc(barpsi_p5_M,mup3,p6,v1,psi_p1_M,-1)
       SMB1M26(i) = SC5ccrrc(barpsi_p5_M,mup6,mup3,p6,p2,v1,psi_p1_M,-1)
       SMB1M27(i) = SC5ccrrc(barpsi_p5_M,mup6,mup3,p6,p4,v1,psi_p1_M,-1)
       SMB1M28(i) = SC5crrrc(barpsi_p2_M,mup3,p6,p5,p1,v1,psi_p4_M,-1)
       SMB1M29(i) = SC5ccrrc(barpsi_p5_M,mup6,mup3,p4,p2,v1,psi_p1_M,-1)
       SMB1M30(i) = SC5crrrc(barpsi_p5_M,mup6,p6,p4,p2,v1,psi_p1_M,-1)
      enddo
c************************************************************************************
c************************************************************************************
      do j=0,3
        do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
        v2(0)=delta(j,0)
        v2(1)=delta(j,1)
        v2(2)=delta(j,2)
        v2(3)=delta(j,3)
       SMB2P1(i,j) = SC3ccc(barpsi_p2_P,v1,mup3,v2,psi_p4_P,+1)
       SMB2P2(i,j) = SC3ccc(barpsi_p5_P,v1,mup6,v2,psi_p1_P,+1)
       SMB2P3(i,j) = SC3crc(barpsi_p2_P,v1,p1,v2,psi_p4_P,+1)
       SMB2P4(i,j) = SC3crc(barpsi_p2_P,v1,p5,v2,psi_p4_P,+1)
       SMB2P5(i,j) = SC3crc(barpsi_p2_P,v1,p6,v2,psi_p4_P,+1)
       SMB2P6(i,j) = SC5crcrc(barpsi_p2_P,mup3,p5,v1,p1,v2,psi_p4_P,+1)
       SMB2P7(i,j) = SC5crcrc(barpsi_p2_P,mup3,p6,v1,p1,v2,psi_p4_P,+1)
       SMB2P8(i,j) = SC5crcrc(barpsi_p2_P,mup3,p6,v1,p5,v2,psi_p4_P,+1)
       SMB2P9(i,j) = SC3crc(barpsi_p5_P,v1,p2,v2,psi_p1_P,+1)
       SMB2P10(i,j) = SC3crc(barpsi_p5_P,v1,p4,v2,psi_p1_P,+1)
       SMB2P11(i,j) = SC3crc(barpsi_p5_P,v1,p6,v2,psi_p1_P,+1)
       SMB2P12(i,j) = SC5crcrc(barpsi_p5_P,mup6,p6,v1,p2,v2,psi_p1_P,+1)
       SMB2P13(i,j) = SC5crcrc(barpsi_p5_P,mup6,p6,v1,p4,v2,psi_p1_P,+1)
       SMB2P14(i,j) = SC3ccc(barpsi_p2_P,v1,mup6,v2,psi_p4_P,+1)
       SMB2P15(i,j) = SC5cccrc(barpsi_p2_P,mup6,mup3,v1,p1,v2,psi_p4_P,+1)
       SMB2P16(i,j) = SC5cccrc(barpsi_p2_P,mup6,mup3,v1,p5,v2,psi_p4_P,+1)
       SMB2P17(i,j) = SC5cccrc(barpsi_p2_P,mup6,mup3,v1,p6,v2,psi_p4_P,+1)
       SMB2P18(i,j) = SC3ccc(barpsi_p5_P,v1,mup3,v2,psi_p1_P,+1)
       SMB2P19(i,j) = SC5cccrc(barpsi_p5_P,mup6,mup3,v1,p2,v2,psi_p1_P,+1)
       SMB2P20(i,j) = SC5crcrc(barpsi_p5_P,mup6,p4,v1,p2,v2,psi_p1_P,+1)
       SMB2P21(i,j) = SC5cccrc(barpsi_p5_P,mup6,mup3,v1,p4,v2,psi_p1_P,+1)
       SMB2P22(i,j) = SC5cccrc(barpsi_p5_P,mup6,mup3,v1,p6,v2,psi_p1_P,+1)
       SMB2M1(i,j) = SC3ccc(barpsi_p2_M,v1,mup3,v2,psi_p4_M,-1)
       SMB2M2(i,j) = SC3ccc(barpsi_p5_M,v1,mup6,v2,psi_p1_M,-1)
       SMB2M3(i,j) = SC3crc(barpsi_p2_M,v1,p1,v2,psi_p4_M,-1)
       SMB2M4(i,j) = SC3crc(barpsi_p2_M,v1,p5,v2,psi_p4_M,-1)
       SMB2M5(i,j) = SC3crc(barpsi_p2_M,v1,p6,v2,psi_p4_M,-1)
       SMB2M6(i,j) = SC5crcrc(barpsi_p2_M,mup3,p5,v1,p1,v2,psi_p4_M,-1)
       SMB2M7(i,j) = SC5crcrc(barpsi_p2_M,mup3,p6,v1,p1,v2,psi_p4_M,-1)
       SMB2M8(i,j) = SC5crcrc(barpsi_p2_M,mup3,p6,v1,p5,v2,psi_p4_M,-1)
       SMB2M9(i,j) = SC3crc(barpsi_p5_M,v1,p2,v2,psi_p1_M,-1)
       SMB2M10(i,j) = SC3crc(barpsi_p5_M,v1,p4,v2,psi_p1_M,-1)
       SMB2M11(i,j) = SC3crc(barpsi_p5_M,v1,p6,v2,psi_p1_M,-1)
       SMB2M12(i,j) = SC5crcrc(barpsi_p5_M,mup6,p6,v1,p2,v2,psi_p1_M,-1)
       SMB2M13(i,j) = SC5crcrc(barpsi_p5_M,mup6,p6,v1,p4,v2,psi_p1_M,-1)
       SMB2M14(i,j) = SC3ccc(barpsi_p2_M,v1,mup6,v2,psi_p4_M,-1)
       SMB2M15(i,j) = SC5cccrc(barpsi_p2_M,mup6,mup3,v1,p1,v2,psi_p4_M,-1)
       SMB2M16(i,j) = SC5cccrc(barpsi_p2_M,mup6,mup3,v1,p5,v2,psi_p4_M,-1)
       SMB2M17(i,j) = SC5cccrc(barpsi_p2_M,mup6,mup3,v1,p6,v2,psi_p4_M,-1)
       SMB2M18(i,j) = SC3ccc(barpsi_p5_M,v1,mup3,v2,psi_p1_M,-1)
       SMB2M19(i,j) = SC5cccrc(barpsi_p5_M,mup6,mup3,v1,p2,v2,psi_p1_M,-1)
       SMB2M20(i,j) = SC5crcrc(barpsi_p5_M,mup6,p4,v1,p2,v2,psi_p1_M,-1)
       SMB2M21(i,j) = SC5cccrc(barpsi_p5_M,mup6,mup3,v1,p4,v2,psi_p1_M,-1)
       SMB2M22(i,j) = SC5cccrc(barpsi_p5_M,mup6,mup3,v1,p6,v2,psi_p1_M,-1)
      enddo
       enddo
c************************************************************************************
c************************************************************************************
      do k=0,3
        do j=0,3
        do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
        v2(0)=delta(j,0)
        v2(1)=delta(j,1)
        v2(2)=delta(j,2)
        v2(3)=delta(j,3)
        v3(0)=delta(k,0)
        v3(1)=delta(k,1)
        v3(2)=delta(k,2)
        v3(3)=delta(k,3)
       SMB3P1(i,j,k) = SC3ccc(barpsi_p2_P,v1,v2,v3,psi_p4_P,+1)
       SMB3P2(i,j,k) = SC3ccc(barpsi_p5_P,v1,v2,v3,psi_p1_P,+1)
       SMB3P3(i,j,k) = SC5crccc(barpsi_p5_P,mup6,p2,v1,v2,v3,psi_p1_P,+1)
       SMB3P4(i,j,k) = SC5crccc(barpsi_p2_P,mup3,p1,v1,v2,v3,psi_p4_P,+1)
       SMB3P5(i,j,k) = SC5crccc(barpsi_p2_P,mup3,p5,v1,v2,v3,psi_p4_P,+1)
       SMB3P6(i,j,k) = SC5crccc(barpsi_p2_P,mup3,p6,v1,v2,v3,psi_p4_P,+1)
       SMB3P7(i,j,k) = SC5crccc(barpsi_p5_P,mup6,p4,v1,v2,v3,psi_p1_P,+1)
       SMB3P8(i,j,k) = SC5crccc(barpsi_p5_P,mup6,p6,v1,v2,v3,psi_p1_P,+1)
       SMB3M1(i,j,k) = SC3ccc(barpsi_p2_M,v1,v2,v3,psi_p4_M,-1)
       SMB3M2(i,j,k) = SC3ccc(barpsi_p5_M,v1,v2,v3,psi_p1_M,-1)
       SMB3M3(i,j,k) = SC5crccc(barpsi_p5_M,mup6,p2,v1,v2,v3,psi_p1_M,-1)
       SMB3M4(i,j,k) = SC5crccc(barpsi_p2_M,mup3,p1,v1,v2,v3,psi_p4_M,-1)
       SMB3M5(i,j,k) = SC5crccc(barpsi_p2_M,mup3,p5,v1,v2,v3,psi_p4_M,-1)
       SMB3M6(i,j,k) = SC5crccc(barpsi_p2_M,mup3,p6,v1,v2,v3,psi_p4_M,-1)
       SMB3M7(i,j,k) = SC5crccc(barpsi_p5_M,mup6,p4,v1,v2,v3,psi_p1_M,-1)
       SMB3M8(i,j,k) = SC5crccc(barpsi_p5_M,mup6,p6,v1,v2,v3,psi_p1_M,-1)
      enddo
       enddo
      enddo
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c   heP-1 heM-1 hMP1 hMM1
c************************************************************************************
       SMB(1) = SMB0M(2)*SMB0P(1)
       SMB(2) = SMB0M(3)*SMB0P(1)
       SMB(3) = SMB0M(2)*SMB0P(4)
       SMB(4) = SMB0M(3)*SMB0P(4)
       SMB(5) = SMB0M(2)*SMB0P(5)
       SMB(6) = SMB0M(3)*SMB0P(5)
       SMB(7) = SMB0M(6)*SMB0P(1)
       SMB(8) = SMB0M(6)*SMB0P(7)
       SMB(9) = SMB0M(3)*SMB0P(7)
       SMB(10) = SMB0M(6)*SMB0P(4)
       SMB(11) = SMB0M(6)*SMB0P(5)
       SMB(12) = SMB0M(2)*SMB0P(7)
       SMB(13) = SMB0M(2)*SMB0P(8)
       SMB(14) = SMB0M(3)*SMB0P(8)
       SMB(15) = SMB0M(2)*SMB0P(9)
       SMB(16) = SMB0M(3)*SMB0P(9)
       SMB(17) = SMB0M(2)*SMB0P(10)
       SMB(18) = SMB0M(2)*SMB0P(11)
       SMB(19) = SMB0M(3)*SMB0P(10)
       SMB(20) = SMB0M(3)*SMB0P(11)
       SMB(21) = SMB0M(6)*SMB0P(8)
       SMB(22) = SMB0M(6)*SMB0P(9)
       SMB(23) = SMB0M(6)*SMB0P(10)
       SMB(24) = SMB0M(6)*SMB0P(11)
       SMB(25) = SMB0M(12)*SMB0P(1)
       SMB(26) = SMB0M(13)*SMB0P(7)
       SMB(29) = SMB0M(12)*SMB0P(4)
       SMB(30) = SMB0M(12)*SMB0P(8)
       SMB(31) = SMB0M(12)*SMB0P(9)
       SMB(32) = SMB0M(12)*SMB0P(14)
       SMB(33) = SMB0M(12)*SMB0P(15)
       SMB(35) = SMB0M(13)*SMB0P(4)
       SMB(36) = SMB0M(16)*SMB0P(4)
       SMB(37) = SMB0M(13)*SMB0P(8)
       SMB(38) = SMB0M(16)*SMB0P(8)
       SMB(39) = SMB0M(13)*SMB0P(9)
       SMB(40) = SMB0M(16)*SMB0P(9)
       SMB(41) = SMB0M(13)*SMB0P(17)
       SMB(42) = SMB0M(13)*SMB0P(18)
       SMB(54) = SMB0M(19)*SMB0P(4)
       SMB(56) = SMB0M(20)*SMB0P(4)
       SMB(60) = SMB0M(20)*SMB0P(8)
       SMB(61) = SMB0M(19)*SMB0P(8)
       SMB(62) = SMB0M(20)*SMB0P(9)
       SMB(63) = SMB0M(19)*SMB0P(9)
       SMB(64) = SMB0M(12)*SMB0P(21)
       SMB(70) = SMB0M(22)*SMB0P(1)
       SMB(71) = SMB0M(23)*SMB0P(1)
       SMB(72) = SMB0M(13)*SMB0P(24)
       SMB(73) = SMB0M(13)*SMB0P(25)
       SMB(74) = SMB0M(12)*SMB0P(25)
       SMB(76) = SMB0M(2)*SMB0P(25)
       SMB(77) = SMB0M(3)*SMB0P(25)
       SMB(79) = SMB0M(16)*SMB0P(25)
       SMB(80) = SMB0M(19)*SMB0P(25)
       SMB(81) = SMB0M(6)*SMB0P(25)
       SMB(82) = SMB0M(20)*SMB0P(25)
       SMB(83) = SMB0M(26)*SMB0P(7)
       SMB(85) = SMB0M(26)*SMB0P(4)
       SMB(86) = SMB0M(26)*SMB0P(8)
       SMB(87) = SMB0M(26)*SMB0P(9)
       SMB(91) = SMB0M(6)*SMB0P(17)
       SMB(92) = SMB0M(6)*SMB0P(24)
       SMB(93) = SMB0M(6)*SMB0P(18)
       SMB(94) = SMB0M(3)*SMB0P(17)
       SMB(95) = SMB0M(3)*SMB0P(24)
       SMB(96) = SMB0M(3)*SMB0P(18)
       SMB(99) = SMB0M(26)*SMB0P(17)
       SMB(101) = SMB0M(26)*SMB0P(24)
       SMB(102) = SMB0M(26)*SMB0P(18)
       SMB(103) = SMB0M(27)*SMB0P(1)
       SMB(108) = SMB0M(28)*SMB0P(7)
       SMB(129) = SMB0M(23)*SMB0P(4)
       SMB(130) = SMB0M(27)*SMB0P(4)
       SMB(131) = SMB0M(23)*SMB0P(8)
       SMB(132) = SMB0M(27)*SMB0P(8)
       SMB(133) = SMB0M(23)*SMB0P(9)
       SMB(134) = SMB0M(27)*SMB0P(9)
       SMB(135) = SMB0M(6)*SMB0P(14)
       SMB(136) = SMB0M(6)*SMB0P(21)
       SMB(137) = SMB0M(6)*SMB0P(15)
       SMB(138) = SMB0M(2)*SMB0P(14)
       SMB(139) = SMB0M(2)*SMB0P(21)
       SMB(140) = SMB0M(2)*SMB0P(15)
       SMB(141) = SMB0M(3)*SMB0P(14)
       SMB(142) = SMB0M(3)*SMB0P(21)
       SMB(143) = SMB0M(3)*SMB0P(15)
       SMB(144) = SMB0M(27)*SMB0P(14)
       SMB(145) = SMB0M(23)*SMB0P(21)
       SMB(146) = SMB0M(27)*SMB0P(21)
       SMB(147) = SMB0M(23)*SMB0P(15)
       SMB(148) = SMB0M(27)*SMB0P(15)
       SMB(149) = SMB0M(29)*SMB0P(1)
       SMB(150) = SMB0M(30)*SMB0P(7)
       SMB(151) = SMB0M(29)*SMB0P(7)
       SMB(152) = SMB0M(29)*SMB0P(4)
       SMB(153) = SMB0M(29)*SMB0P(8)
       SMB(154) = SMB0M(29)*SMB0P(9)
       SMB(155) = SMB0M(29)*SMB0P(5)
       SMB(156) = SMB0M(29)*SMB0P(10)
       SMB(157) = SMB0M(29)*SMB0P(11)
       SMB(158) = SMB0M(30)*SMB0P(4)
       SMB(159) = SMB0M(28)*SMB0P(4)
       SMB(160) = SMB0M(30)*SMB0P(8)
       SMB(161) = SMB0M(28)*SMB0P(8)
       SMB(162) = SMB0M(30)*SMB0P(9)
       SMB(163) = SMB0M(28)*SMB0P(9)
       SMB(164) = SMB0M(2)*SMB0P(17)
       SMB(165) = SMB0M(2)*SMB0P(24)
       SMB(166) = SMB0M(2)*SMB0P(18)
       SMB(171) = SMB0M(28)*SMB0P(17)
       SMB(174) = SMB0M(28)*SMB0P(24)
       SMB(175) = SMB0M(28)*SMB0P(18)
       SMB(176) = SMB0M(22)*SMB0P(4)
       SMB(177) = SMB0M(22)*SMB0P(14)
       SMB(178) = SMB0M(22)*SMB0P(21)
       SMB(179) = SMB0M(22)*SMB0P(8)
       SMB(180) = SMB0M(22)*SMB0P(9)
       SMB(181) = SMB0M(22)*SMB0P(15)
       SMB(206) = SMB0M(6)*SMB0P(31)
       SMB(207) = SMB0M(3)*SMB0P(31)
       SMB(216) = SMB0M(2)*SMB0P(31)
       SMB(246) = SMB0M(32)*SMB0P(4)
       SMB(247) = SMB0M(32)*SMB0P(8)
       SMB(248) = SMB0M(32)*SMB0P(9)
c************************************************************************************
       SMB(27) = dotcc(SMB1M2,SMB1P1)
       SMB(28) = dotcc(SMB1M2,SMB1P3)
       SMB(34) = dotcc(SMB1M4,SMB1P1)
       SMB(43) = dotcc(SMB1M2,SMB1P5)
       SMB(44) = dotcc(SMB1M2,SMB1P6)
       SMB(45) = dotcc(SMB1M4,SMB1P5)
       SMB(46) = dotcc(SMB1M4,SMB1P6)
       SMB(47) = dotcc(SMB1M7,SMB1P1)
       SMB(48) = dotcc(SMB1M2,SMB1P8)
       SMB(49) = dotcc(SMB1M2,SMB1P9)
       SMB(50) = dotcc(SMB1M2,SMB1P10)
       SMB(51) = dotcc(SMB1M7,SMB1P8)
       SMB(52) = dotcc(SMB1M7,SMB1P9)
       SMB(53) = dotcc(SMB1M7,SMB1P10)
       SMB(55) = dotcc(SMB1M11,SMB1P8)
       SMB(57) = dotcc(SMB1M12,SMB1P8)
       SMB(58) = dotcc(SMB1M12,SMB1P1)
       SMB(59) = dotcc(SMB1M11,SMB1P1)
       SMB(66) = dotcc(SMB1M12,SMB1P9)
       SMB(67) = dotcc(SMB1M11,SMB1P9)
       SMB(68) = dotcc(SMB1M12,SMB1P10)
       SMB(69) = dotcc(SMB1M11,SMB1P10)
       SMB(75) = dotcc(SMB1M2,SMB1P13)
       SMB(78) = dotcc(SMB1M4,SMB1P13)
       SMB(84) = dotcc(SMB1M14,SMB1P1)
       SMB(88) = dotcc(SMB1M2,SMB1P15)
       SMB(89) = dotcc(SMB1M2,SMB1P16)
       SMB(90) = dotcc(SMB1M2,SMB1P17)
       SMB(97) = dotcc(SMB1M14,SMB1P15)
       SMB(98) = dotcc(SMB1M14,SMB1P16)
       SMB(100) = dotcc(SMB1M14,SMB1P17)
       SMB(104) = dotcc(SMB1M18,SMB1P1)
       SMB(105) = dotcc(SMB1M19,SMB1P1)
       SMB(106) = dotcc(SMB1M18,SMB1P3)
       SMB(107) = dotcc(SMB1M19,SMB1P3)
       SMB(109) = dotcc(SMB1M20,SMB1P1)
       SMB(110) = dotcc(SMB1M20,SMB1P3)
       SMB(111) = dotcc(SMB1M2,SMB1P21)
       SMB(112) = dotcc(SMB1M2,SMB1P22)
       SMB(118) = dotcc(SMB1M2,SMB1P23)
       SMB(120) = dotcc(SMB1M24,SMB1P1)
       SMB(121) = dotcc(SMB1M25,SMB1P1)
       SMB(122) = dotcc(SMB1M26,SMB1P1)
       SMB(123) = dotcc(SMB1M27,SMB1P1)
       SMB(167) = dotcc(SMB1M24,SMB1P15)
       SMB(168) = dotcc(SMB1M25,SMB1P15)
       SMB(169) = dotcc(SMB1M24,SMB1P16)
       SMB(170) = dotcc(SMB1M25,SMB1P16)
       SMB(172) = dotcc(SMB1M24,SMB1P17)
       SMB(173) = dotcc(SMB1M25,SMB1P17)
       SMB(182) = dotcc(SMB1M19,SMB1P8)
       SMB(183) = dotcc(SMB1M20,SMB1P8)
       SMB(184) = dotcc(SMB1M19,SMB1P9)
       SMB(185) = dotcc(SMB1M20,SMB1P9)
       SMB(186) = dotcc(SMB1M19,SMB1P10)
       SMB(187) = dotcc(SMB1M20,SMB1P10)
       SMB(189) = dotcc(SMB1M19,SMB1P15)
       SMB(190) = dotcc(SMB1M20,SMB1P15)
       SMB(191) = dotcc(SMB1M19,SMB1P16)
       SMB(192) = dotcc(SMB1M20,SMB1P16)
       SMB(193) = dotcc(SMB1M19,SMB1P17)
       SMB(194) = dotcc(SMB1M20,SMB1P17)
       SMB(196) = dotcc(SMB1M19,SMB1P5)
       SMB(197) = dotcc(SMB1M20,SMB1P5)
       SMB(198) = dotcc(SMB1M19,SMB1P13)
       SMB(199) = dotcc(SMB1M20,SMB1P13)
       SMB(201) = dotcc(SMB1M19,SMB1P21)
       SMB(202) = dotcc(SMB1M19,SMB1P22)
       SMB(203) = dotcc(SMB1M20,SMB1P21)
       SMB(204) = dotcc(SMB1M20,SMB1P22)
       SMB(205) = dotcc(SMB1M2,SMB1P28)
       SMB(208) = dotcc(SMB1M12,SMB1P5)
       SMB(209) = dotcc(SMB1M11,SMB1P5)
       SMB(210) = dotcc(SMB1M12,SMB1P13)
       SMB(211) = dotcc(SMB1M11,SMB1P13)
       SMB(212) = dotcc(SMB1M12,SMB1P6)
       SMB(213) = dotcc(SMB1M11,SMB1P6)
       SMB(214) = dotcc(SMB1M12,SMB1P28)
       SMB(215) = dotcc(SMB1M11,SMB1P28)
       SMB(217) = dotcc(SMB1M7,SMB1P5)
       SMB(218) = dotcc(SMB1M7,SMB1P13)
       SMB(219) = dotcc(SMB1M7,SMB1P6)
       SMB(220) = dotcc(SMB1M7,SMB1P28)
       SMB(229) = dotcc(SMB1M19,SMB1P6)
       SMB(233) = dotcc(SMB1M19,SMB1P23)
       SMB(234) = dotcc(SMB1M18,SMB1P8)
       SMB(235) = dotcc(SMB1M18,SMB1P9)
       SMB(236) = dotcc(SMB1M18,SMB1P10)
       SMB(241) = dotcc(SMB1M20,SMB1P6)
       SMB(245) = dotcc(SMB1M20,SMB1P23)
       SMB(251) = dotcc(SMB1M14,SMB1P5)
       SMB(252) = dotcc(SMB1M25,SMB1P5)
       SMB(253) = dotcc(SMB1M26,SMB1P5)
       SMB(255) = dotcc(SMB1M14,SMB1P13)
       SMB(256) = dotcc(SMB1M25,SMB1P13)
       SMB(257) = dotcc(SMB1M26,SMB1P13)
       SMB(258) = dotcc(SMB1M14,SMB1P6)
       SMB(259) = dotcc(SMB1M25,SMB1P6)
       SMB(260) = dotcc(SMB1M26,SMB1P6)
       SMB(261) = dotcc(SMB1M30,SMB1P1)
       SMB(262) = dotcc(SMB1M30,SMB1P8)
       SMB(263) = dotcc(SMB1M30,SMB1P9)
       SMB(264) = dotcc(SMB1M30,SMB1P10)
       SMB(281) = dotcc(SMB1M24,SMB1P5)
       SMB(282) = dotcc(SMB1M24,SMB1P13)
       SMB(283) = dotcc(SMB1M24,SMB1P6)
       SMB(289) = dotcc(SMB1M27,SMB1P5)
       SMB(291) = dotcc(SMB1M27,SMB1P13)
       SMB(292) = dotcc(SMB1M27,SMB1P6)
c************************************************************************************
       SMB(65) = ten2cc(SMB2M2,SMB2P1)
       SMB(113) = ten2cc(SMB2M2,SMB2P3)
       SMB(114) = ten2cc(SMB2M2,SMB2P4)
       SMB(115) = ten2cc(SMB2M2,SMB2P5)
       SMB(116) = ten2cc(SMB2M2,SMB2P6)
       SMB(117) = ten2cc(SMB2M2,SMB2P7)
       SMB(119) = ten2cc(SMB2M2,SMB2P8)
       SMB(124) = ten2cc(SMB2M9,SMB2P1)
       SMB(125) = ten2cc(SMB2M10,SMB2P1)
       SMB(126) = ten2cc(SMB2M11,SMB2P1)
       SMB(127) = ten2cc(SMB2M12,SMB2P1)
       SMB(128) = ten2cc(SMB2M13,SMB2P1)
       SMB(188) = ten2cc(SMB2M11,SMB2P14)
       SMB(195) = ten2cc(SMB2M11,SMB2P3)
       SMB(200) = ten2cc(SMB2M11,SMB2P15)
       SMB(221) = ten2cc(SMB2M11,SMB2P4)
       SMB(222) = ten2cc(SMB2M11,SMB2P5)
       SMB(223) = ten2cc(SMB2M11,SMB2P16)
       SMB(224) = ten2cc(SMB2M11,SMB2P17)
       SMB(225) = ten2cc(SMB2M9,SMB2P14)
       SMB(226) = ten2cc(SMB2M9,SMB2P3)
       SMB(227) = ten2cc(SMB2M9,SMB2P4)
       SMB(228) = ten2cc(SMB2M9,SMB2P5)
       SMB(230) = ten2cc(SMB2M9,SMB2P15)
       SMB(231) = ten2cc(SMB2M9,SMB2P16)
       SMB(232) = ten2cc(SMB2M9,SMB2P17)
       SMB(237) = ten2cc(SMB2M10,SMB2P14)
       SMB(238) = ten2cc(SMB2M10,SMB2P3)
       SMB(239) = ten2cc(SMB2M10,SMB2P4)
       SMB(240) = ten2cc(SMB2M10,SMB2P5)
       SMB(242) = ten2cc(SMB2M10,SMB2P15)
       SMB(243) = ten2cc(SMB2M10,SMB2P16)
       SMB(244) = ten2cc(SMB2M10,SMB2P17)
       SMB(249) = ten2cc(SMB2M18,SMB2P3)
       SMB(250) = ten2cc(SMB2M18,SMB2P4)
       SMB(254) = ten2cc(SMB2M18,SMB2P5)
       SMB(266) = ten2cc(SMB2M12,SMB2P3)
       SMB(267) = ten2cc(SMB2M12,SMB2P4)
       SMB(268) = ten2cc(SMB2M12,SMB2P5)
       SMB(272) = ten2cc(SMB2M9,SMB2P6)
       SMB(273) = ten2cc(SMB2M9,SMB2P7)
       SMB(274) = ten2cc(SMB2M9,SMB2P8)
       SMB(275) = ten2cc(SMB2M11,SMB2P6)
       SMB(276) = ten2cc(SMB2M11,SMB2P7)
       SMB(277) = ten2cc(SMB2M11,SMB2P8)
       SMB(278) = ten2cc(SMB2M12,SMB2P6)
       SMB(279) = ten2cc(SMB2M12,SMB2P7)
       SMB(280) = ten2cc(SMB2M12,SMB2P8)
       SMB(284) = ten2cc(SMB2M10,SMB2P6)
       SMB(285) = ten2cc(SMB2M10,SMB2P7)
       SMB(286) = ten2cc(SMB2M10,SMB2P8)
       SMB(287) = ten2cc(SMB2M22,SMB2P3)
       SMB(288) = ten2cc(SMB2M22,SMB2P4)
       SMB(290) = ten2cc(SMB2M22,SMB2P5)
       SMB(294) = ten2cc(SMB2M13,SMB2P3)
       SMB(295) = ten2cc(SMB2M13,SMB2P4)
       SMB(296) = ten2cc(SMB2M13,SMB2P5)
       SMB(299) = ten2cc(SMB2M13,SMB2P6)
       SMB(301) = ten2cc(SMB2M13,SMB2P7)
       SMB(302) = ten2cc(SMB2M13,SMB2P8)
       SMB(303) = ten2cc(SMB2M18,SMB2P14)
c************************************************************************************
       SMB(265) = ten3cc(SMB3M2,SMB3P1)
       SMB(269) = ten3cc(SMB3M2,SMB3P4)
       SMB(270) = ten3cc(SMB3M2,SMB3P5)
       SMB(271) = ten3cc(SMB3M2,SMB3P6)
       SMB(293) = ten3cc(SMB3M8,SMB3P1)
       SMB(297) = ten3cc(SMB3M8,SMB3P4)
       SMB(298) = ten3cc(SMB3M8,SMB3P5)
       SMB(300) = ten3cc(SMB3M8,SMB3P6)
c************************************************************************************
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
       resultn=(result(1)+result(2)+result(3)+result(4)+result(5)+res
     &   ult(6)+result(7)+result(8)+result(9))
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
c************************************************************************************
c************************************************************************************
c************************************************************************************
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
       Return
       End
