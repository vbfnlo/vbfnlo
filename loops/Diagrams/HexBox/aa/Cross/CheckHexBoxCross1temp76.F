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
       subroutine CheckHexBoxCross176(p1,p2,p3,p4,p5,p6,psi_p1,barpsi_
     &   p6,barpsi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,resu
     &   lt,resultn,SMB)
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
      Complex*16 SMB(293),Fa(272),F(1461)

     &   

       Complex*16 SMB0P(31),SMB1P1(0:3),SMB1P2(0:3),SMB1P3(0:3),SMB1P4
     &   (0:3),SMB1P5(0:3),SMB1P6(0:3),SMB1P7(0:3),SMB1P8(0:3),SMB1P9
     &   (0:3),SMB1P10(0:3),SMB1P11(0:3),SMB1P12(0:3),SMB1P13(0:3),SM
     &   B1P14(0:3),SMB1P15(0:3),SMB1P16(0:3),SMB1P17(0:3),SMB1P18(0:
     &   3),SMB1P19(0:3),SMB1P20(0:3),SMB1P21(0:3),SMB1P22(0:3),SMB1P
     &   23(0:3),SMB1P24(0:3),SMB1P25(0:3),SMB1P26(0:3),SMB1P27(0:3),
     &   SMB2P1(0:3,0:3),SMB2P2(0:3,0:3),SMB2P3(0:3,0:3),SMB2P4(0:3,0
     &   :3),SMB2P5(0:3,0:3),SMB2P6(0:3,0:3),SMB2P7(0:3,0:3),SMB2P8(0
     &   :3,0:3),SMB2P9(0:3,0:3),SMB2P10(0:3,0:3),SMB2P11(0:3,0:3),SM
     &   B2P12(0:3,0:3),SMB2P13(0:3,0:3),SMB2P14(0:3,0:3),SMB2P15(0:3
     &   ,0:3),SMB2P16(0:3,0:3),SMB2P17(0:3,0:3),SMB2P18(0:3,0:3),SMB
     &   2P19(0:3,0:3),SMB2P20(0:3,0:3),SMB2P21(0:3,0:3),SMB3P1(0:3,0
     &   :3,0:3),SMB3P2(0:3,0:3,0:3),SMB3P3(0:3,0:3,0:3),SMB3P4(0:3,0
     &   :3,0:3),SMB3P5(0:3,0:3,0:3),SMB3P6(0:3,0:3,0:3),SMB3P7(0:3,0
     &   :3,0:3),SMB3P8(0:3,0:3,0:3),SMB3P9(0:3,0:3,0:3),SMB3P10(0:3,
     &   0:3,0:3),SMB3P11(0:3,0:3,0:3),SMB3P12(0:3,0:3,0:3),SMB3P13(0
     &   :3,0:3,0:3),SMB3P14(0:3,0:3,0:3),SMB3P15(0:3,0:3,0:3)
       Complex*16 SMB0M(31),SMB1M1(0:3),SMB1M2(0:3),SMB1M3(0:3),SMB1M4
     &   (0:3),SMB1M5(0:3),SMB1M6(0:3),SMB1M7(0:3),SMB1M8(0:3),SMB1M9
     &   (0:3),SMB1M10(0:3),SMB1M11(0:3),SMB1M12(0:3),SMB1M13(0:3),SM
     &   B1M14(0:3),SMB1M15(0:3),SMB1M16(0:3),SMB1M17(0:3),SMB1M18(0:
     &   3),SMB1M19(0:3),SMB1M20(0:3),SMB1M21(0:3),SMB1M22(0:3),SMB1M
     &   23(0:3),SMB1M24(0:3),SMB1M25(0:3),SMB1M26(0:3),SMB1M27(0:3),
     &   SMB2M1(0:3,0:3),SMB2M2(0:3,0:3),SMB2M3(0:3,0:3),SMB2M4(0:3,0
     &   :3),SMB2M5(0:3,0:3),SMB2M6(0:3,0:3),SMB2M7(0:3,0:3),SMB2M8(0
     &   :3,0:3),SMB2M9(0:3,0:3),SMB2M10(0:3,0:3),SMB2M11(0:3,0:3),SM
     &   B2M12(0:3,0:3),SMB2M13(0:3,0:3),SMB2M14(0:3,0:3),SMB2M15(0:3
     &   ,0:3),SMB2M16(0:3,0:3),SMB2M17(0:3,0:3),SMB2M18(0:3,0:3),SMB
     &   2M19(0:3,0:3),SMB2M20(0:3,0:3),SMB2M21(0:3,0:3),SMB3M1(0:3,0
     &   :3,0:3),SMB3M2(0:3,0:3,0:3),SMB3M3(0:3,0:3,0:3),SMB3M4(0:3,0
     &   :3,0:3),SMB3M5(0:3,0:3,0:3),SMB3M6(0:3,0:3,0:3),SMB3M7(0:3,0
     &   :3,0:3),SMB3M8(0:3,0:3,0:3),SMB3M9(0:3,0:3,0:3),SMB3M10(0:3,
     &   0:3,0:3),SMB3M11(0:3,0:3,0:3),SMB3M12(0:3,0:3,0:3),SMB3M13(0
     &   :3,0:3,0:3),SMB3M14(0:3,0:3,0:3),SMB3M15(0:3,0:3,0:3)

      Complex*16 psi_p1_P(2),barpsi_p6_P(2),barpsi_p2_P(2),psi_p5_P(2) 
       Complex*16 psi_p1_M(2),barpsi_p6_M(2),barpsi_p2_M(2),psi_p5_M(2)
       Complex*16 psi_p1(4),barpsi_p6(4),barpsi_p2(4),psi_p5(4)
       Real*8 P(1825)
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
       Complex*16   mup3(0:3),mup4(0:3)
       Complex*16   p1mup3, p1mup4, p2mup3, p2mup4, p3mup3, p3mup4, 
     -          p4mup3, p4mup4, p5mup3, p5mup4, p6mup3, p6mup4
       Complex*16   mup3mup4
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
      Complex*16 dotrc,dotcc,result(6),resultn
       Real*8 musq
      EXTERNAL   dotrc,dotcc

       COMMON/HexBoxCross1FaFunctions/Fa
       COMMON/HexBoxCross1FhlFunctions/F
      Save/HexBoxCross1FhlFunctions/
       COMMON/HexBoxCross1PFunctions/P
      Save/HexBoxCross1PFunctions/
       COMMON/HexBoxCross1Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12
     &   ,s23,s34,s45,s56,s16,s123,s234,s345
       COMMON/HexBoxCross1FVALFunctions/F0123456R,F123456R,F0123456I,F
     &   123456I
       COMMON/HexBoxCross1EVALFunctions/ E012345R,E12345R,E012345I,E12
     &   345I, E012346R,E12346R,E012346I,E12346I, E012356R,E12356R,E0
     &   12356I,E12356I, E012456R,E12456R,E012456I,E12456I, E013456R,
     &   E13456R,E013456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/HexBoxCross1DVALFunctions/ D01234R,D1234R,D01234I,D1234I
     &   , D01235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D123
     &   6I, D01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1
     &   246I, D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,
     &   D1345I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356
     &   I,D1356I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D023
     &   45I,D2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D0
     &   2356I,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,
     &   D03456I,D3456I
       COMMON/HexBoxCross1CVALFunctions/ C0123R,C123R,C0123I,C123I, C0
     &   124R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C
     &   126R,C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C
     &   0135I,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,
     &   C145I, C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I,
     &    C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236
     &   R,C236R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246
     &   R,C0246I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C034
     &   5I,C345I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C35
     &   6I, C0456R,C456R,C0456I,C456I
       COMMON/HexBoxCross1BVALFunctions/ B012R,B012I, B013R,B013I, B01
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
       call HexBoxCross1FFhl1(F(1))
       call HexBoxCross1FFhl2(F(74))
       call HexBoxCross1FFhl3(F(147))
       call HexBoxCross1FFhl4(F(220))
       call HexBoxCross1FFhl5(F(293))
       call HexBoxCross1FFhl6(F(366))
       call HexBoxCross1FFhl7(F(439))
       call HexBoxCross1FFhl8(F(512))
       call HexBoxCross1FFhl9(F(585))
       call HexBoxCross1FFhl10(F(658))
       call HexBoxCross1FFhl11(F(731))
       call HexBoxCross1FFhl12(F(804))
       call HexBoxCross1FFhl13(F(877))
       call HexBoxCross1FFhl14(F(950))
       call HexBoxCross1FFhl15(F(1023))
       call HexBoxCross1FFhl16(F(1096))
       call HexBoxCross1FFhl17(F(1169))
       call HexBoxCross1FFhl18(F(1242))
       call HexBoxCross1FFhl19(F(1315))
       call HexBoxCross1FFhl20(F(1388))
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
       call  HexBoxCross1Fa1(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(1))
       call  HexBoxCross1Fa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(14))
       call  HexBoxCross1Fa3(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(27))
       call  HexBoxCross1Fa4(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(40))
       call  HexBoxCross1Fa5(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(53))
       call  HexBoxCross1Fa6(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(66))
       call  HexBoxCross1Fa7(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(79))
       call  HexBoxCross1Fa8(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(92))
       call  HexBoxCross1Fa9(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4
     &   ,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(105))
       call  HexBoxCross1Fa10(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(118)
     &   )
       call  HexBoxCross1Fa11(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(131)
     &   )
       call  HexBoxCross1Fa12(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(144)
     &   )
       call  HexBoxCross1Fa13(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(157)
     &   )
       call  HexBoxCross1Fa14(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(170)
     &   )
       call  HexBoxCross1Fa15(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(183)
     &   )
       call  HexBoxCross1Fa16(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(196)
     &   )
       call  HexBoxCross1Fa17(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(209)
     &   )
       call  HexBoxCross1Fa18(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(222)
     &   )
       call  HexBoxCross1Fa19(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(235)
     &   )
       call  HexBoxCross1Fa20(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup
     &   4,p4mup3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(248)
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
       barpsi_p6_P(1)=barpsi_p6(1)
       barpsi_p6_P(2)=barpsi_p6(2)
       barpsi_p6_M(1)=barpsi_p6(3)
       barpsi_p6_M(2)=barpsi_p6(4)
       psi_p5_M(1)=psi_p5(1)
       psi_p5_M(2)=psi_p5(2)
       psi_p5_P(1)=psi_p5(3)
       psi_p5_P(2)=psi_p5(4)
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
       SMB0P(1) = SC1c(barpsi_p2_P,mup3,psi_p5_P,+1)
       SMB0P(2) = SC1c(barpsi_p6_P,mup4,psi_p1_P,+1)
       SMB0P(3) = SC1c(barpsi_p2_P,mup4,psi_p5_P,+1)
       SMB0P(4) = SC1c(barpsi_p6_P,mup3,psi_p1_P,+1)
       SMB0P(5) = SC1r(barpsi_p2_P,p1,psi_p5_P,+1)
       SMB0P(6) = SC1r(barpsi_p6_P,p2,psi_p1_P,+1)
       SMB0P(7) = SC3ccr(barpsi_p2_P,mup4,mup3,p1,psi_p5_P,+1)
       SMB0P(8) = SC1r(barpsi_p6_P,p4,psi_p1_P,+1)
       SMB0P(9) = SC1r(barpsi_p6_P,p5,psi_p1_P,+1)
       SMB0P(10) = SC1r(barpsi_p2_P,p4,psi_p5_P,+1)
       SMB0P(11) = SC3ccr(barpsi_p6_P,mup4,mup3,p4,psi_p1_P,+1)
       SMB0P(12) = SC3ccr(barpsi_p6_P,mup4,mup3,p5,psi_p1_P,+1)
       SMB0P(13) = SC3crr(barpsi_p2_P,mup3,p4,p1,psi_p5_P,+1)
       SMB0P(14) = SC3ccr(barpsi_p2_P,mup4,mup3,p4,psi_p5_P,+1)
       SMB0P(15) = SC3ccr(barpsi_p6_P,mup4,mup3,p2,psi_p1_P,+1)
       SMB0P(16) = SC3crr(barpsi_p6_P,mup4,p5,p4,psi_p1_P,+1)
       SMB0P(17) = SC3crr(barpsi_p6_P,mup3,p5,p4,psi_p1_P,+1)
       SMB0P(18) = SC3crr(barpsi_p6_P,mup4,p5,p2,psi_p1_P,+1)
       SMB0P(19) = SC3crr(barpsi_p6_P,mup3,p5,p2,psi_p1_P,+1)
       SMB0P(20) = SC1r(barpsi_p2_P,p6,psi_p5_P,+1)
       SMB0P(21) = SC3crr(barpsi_p2_P,mup4,p4,p1,psi_p5_P,+1)
       SMB0P(22) = SC3crr(barpsi_p2_P,mup4,p6,p4,psi_p5_P,+1)
       SMB0P(23) = SC3rrr(barpsi_p2_P,p6,p4,p1,psi_p5_P,+1)
       SMB0P(24) = SC3crr(barpsi_p2_P,mup3,p6,p4,psi_p5_P,+1)
       SMB0P(25) = SC3crr(barpsi_p2_P,mup4,p6,p1,psi_p5_P,+1)
       SMB0P(26) = SC3ccr(barpsi_p2_P,mup4,mup3,p6,psi_p5_P,+1)
       SMB0P(27) = SC3crr(barpsi_p2_P,mup3,p6,p1,psi_p5_P,+1)
       SMB0P(28) = SC3crr(barpsi_p6_P,mup4,p4,p2,psi_p1_P,+1)
       SMB0P(29) = SC3rrr(barpsi_p6_P,p5,p4,p2,psi_p1_P,+1)
       SMB0P(30) = SC3crr(barpsi_p6_P,mup3,p4,p2,psi_p1_P,+1)
       SMB0P(31) = SC5ccrrr(barpsi_p2_P,mup4,mup3,p6,p4,p1,psi_p5_P,+1)
       SMB0M(1) = SC1c(barpsi_p2_M,mup3,psi_p5_M,-1)
       SMB0M(2) = SC1c(barpsi_p6_M,mup4,psi_p1_M,-1)
       SMB0M(3) = SC1c(barpsi_p2_M,mup4,psi_p5_M,-1)
       SMB0M(4) = SC1c(barpsi_p6_M,mup3,psi_p1_M,-1)
       SMB0M(5) = SC1r(barpsi_p2_M,p1,psi_p5_M,-1)
       SMB0M(6) = SC1r(barpsi_p6_M,p2,psi_p1_M,-1)
       SMB0M(7) = SC3ccr(barpsi_p2_M,mup4,mup3,p1,psi_p5_M,-1)
       SMB0M(8) = SC1r(barpsi_p6_M,p4,psi_p1_M,-1)
       SMB0M(9) = SC1r(barpsi_p6_M,p5,psi_p1_M,-1)
       SMB0M(10) = SC1r(barpsi_p2_M,p4,psi_p5_M,-1)
       SMB0M(11) = SC3ccr(barpsi_p6_M,mup4,mup3,p4,psi_p1_M,-1)
       SMB0M(12) = SC3ccr(barpsi_p6_M,mup4,mup3,p5,psi_p1_M,-1)
       SMB0M(13) = SC3crr(barpsi_p2_M,mup3,p4,p1,psi_p5_M,-1)
       SMB0M(14) = SC3ccr(barpsi_p2_M,mup4,mup3,p4,psi_p5_M,-1)
       SMB0M(15) = SC3ccr(barpsi_p6_M,mup4,mup3,p2,psi_p1_M,-1)
       SMB0M(16) = SC3crr(barpsi_p6_M,mup4,p5,p4,psi_p1_M,-1)
       SMB0M(17) = SC3crr(barpsi_p6_M,mup3,p5,p4,psi_p1_M,-1)
       SMB0M(18) = SC3crr(barpsi_p6_M,mup4,p5,p2,psi_p1_M,-1)
       SMB0M(19) = SC3crr(barpsi_p6_M,mup3,p5,p2,psi_p1_M,-1)
       SMB0M(20) = SC1r(barpsi_p2_M,p6,psi_p5_M,-1)
       SMB0M(21) = SC3crr(barpsi_p2_M,mup4,p4,p1,psi_p5_M,-1)
       SMB0M(22) = SC3crr(barpsi_p2_M,mup4,p6,p4,psi_p5_M,-1)
       SMB0M(23) = SC3rrr(barpsi_p2_M,p6,p4,p1,psi_p5_M,-1)
       SMB0M(24) = SC3crr(barpsi_p2_M,mup3,p6,p4,psi_p5_M,-1)
       SMB0M(25) = SC3crr(barpsi_p2_M,mup4,p6,p1,psi_p5_M,-1)
       SMB0M(26) = SC3ccr(barpsi_p2_M,mup4,mup3,p6,psi_p5_M,-1)
       SMB0M(27) = SC3crr(barpsi_p2_M,mup3,p6,p1,psi_p5_M,-1)
       SMB0M(28) = SC3crr(barpsi_p6_M,mup4,p4,p2,psi_p1_M,-1)
       SMB0M(29) = SC3rrr(barpsi_p6_M,p5,p4,p2,psi_p1_M,-1)
       SMB0M(30) = SC3crr(barpsi_p6_M,mup3,p4,p2,psi_p1_M,-1)
       SMB0M(31) = SC5ccrrr(barpsi_p2_M,mup4,mup3,p6,p4,p1,psi_p5_M,-1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB1P1(i) = SC1c(barpsi_p2_P,v1,psi_p5_P,+1)
       SMB1P2(i) = SC1c(barpsi_p6_P,v1,psi_p1_P,+1)
       SMB1P3(i) = SC3ccc(barpsi_p2_P,mup4,mup3,v1,psi_p5_P,+1)
       SMB1P4(i) = SC3crc(barpsi_p6_P,mup4,p4,v1,psi_p1_P,+1)
       SMB1P5(i) = SC3crc(barpsi_p6_P,mup4,p5,v1,psi_p1_P,+1)
       SMB1P6(i) = SC3crc(barpsi_p2_P,mup3,p4,v1,psi_p5_P,+1)
       SMB1P7(i) = SC3crc(barpsi_p6_P,mup4,p2,v1,psi_p1_P,+1)
       SMB1P8(i) = SC3rrc(barpsi_p6_P,p5,p4,v1,psi_p1_P,+1)
       SMB1P9(i) = SC3rrc(barpsi_p6_P,p5,p2,v1,psi_p1_P,+1)
       SMB1P10(i) = SC3ccc(barpsi_p6_P,mup4,mup3,v1,psi_p1_P,+1)
       SMB1P11(i) = SC3rrc(barpsi_p2_P,p4,p1,v1,psi_p5_P,+1)
       SMB1P12(i) = SC3rrc(barpsi_p2_P,p6,p4,v1,psi_p5_P,+1)
       SMB1P13(i) = SC3crc(barpsi_p6_P,mup3,p5,v1,psi_p1_P,+1)
       SMB1P14(i) = SC3crc(barpsi_p2_P,mup4,p1,v1,psi_p5_P,+1)
       SMB1P15(i) = SC3crc(barpsi_p2_P,mup4,p6,v1,psi_p5_P,+1)
       SMB1P16(i) = SC3crc(barpsi_p6_P,mup3,p4,v1,psi_p1_P,+1)
       SMB1P17(i) = SC3crc(barpsi_p2_P,mup4,p4,v1,psi_p5_P,+1)
       SMB1P18(i) = SC3crc(barpsi_p2_P,mup3,p1,v1,psi_p5_P,+1)
       SMB1P19(i) = SC3crc(barpsi_p2_P,mup3,p6,v1,psi_p5_P,+1)
       SMB1P20(i) = SC3crc(barpsi_p6_P,mup3,p2,v1,psi_p1_P,+1)
       SMB1P21(i) = SC3rrc(barpsi_p6_P,p4,p2,v1,psi_p1_P,+1)
       SMB1P22(i) = SC3rrc(barpsi_p2_P,p6,p1,v1,psi_p5_P,+1)
       SMB1P23(i) = SC5ccrrc(barpsi_p2_P,mup4,mup3,p4,p1,v1,psi_p5_P,+1)
       SMB1P24(i) = SC5ccrrc(barpsi_p2_P,mup4,mup3,p6,p1,v1,psi_p5_P,+1)
       SMB1P25(i) = SC5ccrrc(barpsi_p2_P,mup4,mup3,p6,p4,v1,psi_p5_P,+1)
       SMB1P26(i) = SC5crrrc(barpsi_p2_P,mup3,p6,p4,p1,v1,psi_p5_P,+1)
       SMB1P27(i) = SC5crrrc(barpsi_p2_P,mup4,p6,p4,p1,v1,psi_p5_P,+1)
       SMB1M1(i) = SC1c(barpsi_p2_M,v1,psi_p5_M,-1)
       SMB1M2(i) = SC1c(barpsi_p6_M,v1,psi_p1_M,-1)
       SMB1M3(i) = SC3ccc(barpsi_p2_M,mup4,mup3,v1,psi_p5_M,-1)
       SMB1M4(i) = SC3crc(barpsi_p6_M,mup4,p4,v1,psi_p1_M,-1)
       SMB1M5(i) = SC3crc(barpsi_p6_M,mup4,p5,v1,psi_p1_M,-1)
       SMB1M6(i) = SC3crc(barpsi_p2_M,mup3,p4,v1,psi_p5_M,-1)
       SMB1M7(i) = SC3crc(barpsi_p6_M,mup4,p2,v1,psi_p1_M,-1)
       SMB1M8(i) = SC3rrc(barpsi_p6_M,p5,p4,v1,psi_p1_M,-1)
       SMB1M9(i) = SC3rrc(barpsi_p6_M,p5,p2,v1,psi_p1_M,-1)
       SMB1M10(i) = SC3ccc(barpsi_p6_M,mup4,mup3,v1,psi_p1_M,-1)
       SMB1M11(i) = SC3rrc(barpsi_p2_M,p4,p1,v1,psi_p5_M,-1)
       SMB1M12(i) = SC3rrc(barpsi_p2_M,p6,p4,v1,psi_p5_M,-1)
       SMB1M13(i) = SC3crc(barpsi_p6_M,mup3,p5,v1,psi_p1_M,-1)
       SMB1M14(i) = SC3crc(barpsi_p2_M,mup4,p1,v1,psi_p5_M,-1)
       SMB1M15(i) = SC3crc(barpsi_p2_M,mup4,p6,v1,psi_p5_M,-1)
       SMB1M16(i) = SC3crc(barpsi_p6_M,mup3,p4,v1,psi_p1_M,-1)
       SMB1M17(i) = SC3crc(barpsi_p2_M,mup4,p4,v1,psi_p5_M,-1)
       SMB1M18(i) = SC3crc(barpsi_p2_M,mup3,p1,v1,psi_p5_M,-1)
       SMB1M19(i) = SC3crc(barpsi_p2_M,mup3,p6,v1,psi_p5_M,-1)
       SMB1M20(i) = SC3crc(barpsi_p6_M,mup3,p2,v1,psi_p1_M,-1)
       SMB1M21(i) = SC3rrc(barpsi_p6_M,p4,p2,v1,psi_p1_M,-1)
       SMB1M22(i) = SC3rrc(barpsi_p2_M,p6,p1,v1,psi_p5_M,-1)
       SMB1M23(i) = SC5ccrrc(barpsi_p2_M,mup4,mup3,p4,p1,v1,psi_p5_M,-1)
       SMB1M24(i) = SC5ccrrc(barpsi_p2_M,mup4,mup3,p6,p1,v1,psi_p5_M,-1)
       SMB1M25(i) = SC5ccrrc(barpsi_p2_M,mup4,mup3,p6,p4,v1,psi_p5_M,-1)
       SMB1M26(i) = SC5crrrc(barpsi_p2_M,mup3,p6,p4,p1,v1,psi_p5_M,-1)
       SMB1M27(i) = SC5crrrc(barpsi_p2_M,mup4,p6,p4,p1,v1,psi_p5_M,-1)
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
       SMB2P1(i,j) = SC3ccc(barpsi_p2_P,v1,mup4,v2,psi_p5_P,+1)
       SMB2P2(i,j) = SC3ccc(barpsi_p6_P,v1,mup3,v2,psi_p1_P,+1)
       SMB2P3(i,j) = SC3ccc(barpsi_p2_P,v1,mup3,v2,psi_p5_P,+1)
       SMB2P4(i,j) = SC3ccc(barpsi_p6_P,v1,mup4,v2,psi_p1_P,+1)
       SMB2P5(i,j) = SC3crc(barpsi_p6_P,v1,p2,v2,psi_p1_P,+1)
       SMB2P6(i,j) = SC3crc(barpsi_p6_P,v1,p4,v2,psi_p1_P,+1)
       SMB2P7(i,j) = SC3crc(barpsi_p6_P,v1,p5,v2,psi_p1_P,+1)
       SMB2P8(i,j) = SC3crc(barpsi_p2_P,v1,p1,v2,psi_p5_P,+1)
       SMB2P9(i,j) = SC3crc(barpsi_p2_P,v1,p4,v2,psi_p5_P,+1)
       SMB2P10(i,j) = SC3crc(barpsi_p2_P,v1,p6,v2,psi_p5_P,+1)
       SMB2P11(i,j) = SC5cccrc(barpsi_p2_P,mup4,mup3,v1,p1,v2,psi_p5_P,+1)
       SMB2P12(i,j) = SC5cccrc(barpsi_p2_P,mup4,mup3,v1,p4,v2,psi_p5_P,+1)
       SMB2P13(i,j) = SC5cccrc(barpsi_p2_P,mup4,mup3,v1,p6,v2,psi_p5_P,+1)
       SMB2P14(i,j) = SC5crcrc(barpsi_p2_P,mup3,p4,v1,p1,v2,psi_p5_P,+1)
       SMB2P15(i,j) = SC5crcrc(barpsi_p2_P,mup3,p6,v1,p4,v2,psi_p5_P,+1)
       SMB2P16(i,j) = SC5crcrc(barpsi_p2_P,mup4,p4,v1,p1,v2,psi_p5_P,+1)
       SMB2P17(i,j) = SC5crcrc(barpsi_p2_P,mup4,p6,v1,p4,v2,psi_p5_P,+1)
       SMB2P18(i,j) = SC5crcrc(barpsi_p2_P,mup4,p6,v1,p1,v2,psi_p5_P,+1)
       SMB2P19(i,j) = SC5crcrc(barpsi_p2_P,mup3,p6,v1,p1,v2,psi_p5_P,+1)
       SMB2P20(i,j) = SC5rrcrc(barpsi_p2_P,p6,p4,v1,p1,v2,psi_p5_P,+1)
       SMB2P21(i,j) = SC7ccrrcrc(barpsi_p2_P,mup4,mup3,p6,p4,v1,p1,v2,psi_p5_P,+1)
       SMB2M1(i,j) = SC3ccc(barpsi_p2_M,v1,mup4,v2,psi_p5_M,-1)
       SMB2M2(i,j) = SC3ccc(barpsi_p6_M,v1,mup3,v2,psi_p1_M,-1)
       SMB2M3(i,j) = SC3ccc(barpsi_p2_M,v1,mup3,v2,psi_p5_M,-1)
       SMB2M4(i,j) = SC3ccc(barpsi_p6_M,v1,mup4,v2,psi_p1_M,-1)
       SMB2M5(i,j) = SC3crc(barpsi_p6_M,v1,p2,v2,psi_p1_M,-1)
       SMB2M6(i,j) = SC3crc(barpsi_p6_M,v1,p4,v2,psi_p1_M,-1)
       SMB2M7(i,j) = SC3crc(barpsi_p6_M,v1,p5,v2,psi_p1_M,-1)
       SMB2M8(i,j) = SC3crc(barpsi_p2_M,v1,p1,v2,psi_p5_M,-1)
       SMB2M9(i,j) = SC3crc(barpsi_p2_M,v1,p4,v2,psi_p5_M,-1)
       SMB2M10(i,j) = SC3crc(barpsi_p2_M,v1,p6,v2,psi_p5_M,-1)
       SMB2M11(i,j) = SC5cccrc(barpsi_p2_M,mup4,mup3,v1,p1,v2,psi_p5_M,-1)
       SMB2M12(i,j) = SC5cccrc(barpsi_p2_M,mup4,mup3,v1,p4,v2,psi_p5_M,-1)
       SMB2M13(i,j) = SC5cccrc(barpsi_p2_M,mup4,mup3,v1,p6,v2,psi_p5_M,-1)
       SMB2M14(i,j) = SC5crcrc(barpsi_p2_M,mup3,p4,v1,p1,v2,psi_p5_M,-1)
       SMB2M15(i,j) = SC5crcrc(barpsi_p2_M,mup3,p6,v1,p4,v2,psi_p5_M,-1)
       SMB2M16(i,j) = SC5crcrc(barpsi_p2_M,mup4,p4,v1,p1,v2,psi_p5_M,-1)
       SMB2M17(i,j) = SC5crcrc(barpsi_p2_M,mup4,p6,v1,p4,v2,psi_p5_M,-1)
       SMB2M18(i,j) = SC5crcrc(barpsi_p2_M,mup4,p6,v1,p1,v2,psi_p5_M,-1)
       SMB2M19(i,j) = SC5crcrc(barpsi_p2_M,mup3,p6,v1,p1,v2,psi_p5_M,-1)
       SMB2M20(i,j) = SC5rrcrc(barpsi_p2_M,p6,p4,v1,p1,v2,psi_p5_M,-1)
       SMB2M21(i,j) = SC7ccrrcrc(barpsi_p2_M,mup4,mup3,p6,p4,v1,p1,v2,psi_p5_M,-1)
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
       SMB3P1(i,j,k) = SC3ccc(barpsi_p2_P,v1,v2,v3,psi_p5_P,+1)
       SMB3P2(i,j,k) = SC3ccc(barpsi_p6_P,v1,v2,v3,psi_p1_P,+1)
       SMB3P3(i,j,k) = SC5ccccc(barpsi_p2_P,mup4,mup3,v1,v2,v3,psi_p5_P,+1)
       SMB3P4(i,j,k) = SC5crccc(barpsi_p2_P,mup3,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P5(i,j,k) = SC5crccc(barpsi_p2_P,mup3,p4,v1,v2,v3,psi_p5_P,+1)
       SMB3P6(i,j,k) = SC5crccc(barpsi_p2_P,mup4,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P7(i,j,k) = SC5crccc(barpsi_p2_P,mup4,p4,v1,v2,v3,psi_p5_P,+1)
       SMB3P8(i,j,k) = SC5rrccc(barpsi_p2_P,p4,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P9(i,j,k) = SC7ccrrccc(barpsi_p2_P,mup4,mup3,p4,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P10(i,j,k) = SC5crccc(barpsi_p2_P,mup3,p6,v1,v2,v3,psi_p5_P,+1)
       SMB3P11(i,j,k) = SC5crccc(barpsi_p2_P,mup4,p6,v1,v2,v3,psi_p5_P,+1)
       SMB3P12(i,j,k) = SC5rrccc(barpsi_p2_P,p6,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P13(i,j,k) = SC7ccrrccc(barpsi_p2_P,mup4,mup3,p6,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P14(i,j,k) = SC5rrccc(barpsi_p2_P,p6,p4,v1,v2,v3,psi_p5_P,+1)
       SMB3P15(i,j,k) = SC7ccrrccc(barpsi_p2_P,mup4,mup3,p6,p4,v1,v2,v3,psi_p5_P,+1)
       SMB3M1(i,j,k) = SC3ccc(barpsi_p2_M,v1,v2,v3,psi_p5_M,-1)
       SMB3M2(i,j,k) = SC3ccc(barpsi_p6_M,v1,v2,v3,psi_p1_M,-1)
       SMB3M3(i,j,k) = SC5ccccc(barpsi_p2_M,mup4,mup3,v1,v2,v3,psi_p5_M,-1)
       SMB3M4(i,j,k) = SC5crccc(barpsi_p2_M,mup3,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M5(i,j,k) = SC5crccc(barpsi_p2_M,mup3,p4,v1,v2,v3,psi_p5_M,-1)
       SMB3M6(i,j,k) = SC5crccc(barpsi_p2_M,mup4,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M7(i,j,k) = SC5crccc(barpsi_p2_M,mup4,p4,v1,v2,v3,psi_p5_M,-1)
       SMB3M8(i,j,k) = SC5rrccc(barpsi_p2_M,p4,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M9(i,j,k) = SC7ccrrccc(barpsi_p2_M,mup4,mup3,p4,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M10(i,j,k) = SC5crccc(barpsi_p2_M,mup3,p6,v1,v2,v3,psi_p5_M,-1)
       SMB3M11(i,j,k) = SC5crccc(barpsi_p2_M,mup4,p6,v1,v2,v3,psi_p5_M,-1)
       SMB3M12(i,j,k) = SC5rrccc(barpsi_p2_M,p6,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M13(i,j,k) = SC7ccrrccc(barpsi_p2_M,mup4,mup3,p6,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M14(i,j,k) = SC5rrccc(barpsi_p2_M,p6,p4,v1,v2,v3,psi_p5_M,-1)
       SMB3M15(i,j,k) = SC7ccrrccc(barpsi_p2_M,mup4,mup3,p6,p4,v1,v2,v3,psi_p5_M,-1)
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
       SMB(2) = SMB0M(4)*SMB0P(3)
       SMB(5) = SMB0M(6)*SMB0P(5)
       SMB(6) = SMB0M(6)*SMB0P(7)
       SMB(7) = SMB0M(8)*SMB0P(1)
       SMB(8) = SMB0M(9)*SMB0P(1)
       SMB(9) = SMB0M(8)*SMB0P(3)
       SMB(10) = SMB0M(9)*SMB0P(3)
       SMB(13) = SMB0M(2)*SMB0P(5)
       SMB(14) = SMB0M(4)*SMB0P(10)
       SMB(15) = SMB0M(2)*SMB0P(10)
       SMB(16) = SMB0M(11)*SMB0P(10)
       SMB(17) = SMB0M(12)*SMB0P(10)
       SMB(18) = SMB0M(2)*SMB0P(13)
       SMB(20) = SMB0M(8)*SMB0P(14)
       SMB(21) = SMB0M(9)*SMB0P(14)
       SMB(24) = SMB0M(6)*SMB0P(1)
       SMB(25) = SMB0M(6)*SMB0P(3)
       SMB(27) = SMB0M(15)*SMB0P(10)
       SMB(28) = SMB0M(6)*SMB0P(14)
       SMB(30) = SMB0M(9)*SMB0P(5)
       SMB(31) = SMB0M(9)*SMB0P(7)
       SMB(32) = SMB0M(16)*SMB0P(1)
       SMB(33) = SMB0M(17)*SMB0P(3)
       SMB(36) = SMB0M(18)*SMB0P(1)
       SMB(37) = SMB0M(19)*SMB0P(3)
       SMB(40) = SMB0M(8)*SMB0P(10)
       SMB(41) = SMB0M(9)*SMB0P(10)
       SMB(43) = SMB0M(4)*SMB0P(5)
       SMB(44) = SMB0M(8)*SMB0P(5)
       SMB(45) = SMB0M(11)*SMB0P(5)
       SMB(46) = SMB0M(12)*SMB0P(5)
       SMB(47) = SMB0M(4)*SMB0P(20)
       SMB(48) = SMB0M(2)*SMB0P(20)
       SMB(49) = SMB0M(8)*SMB0P(20)
       SMB(50) = SMB0M(9)*SMB0P(20)
       SMB(51) = SMB0M(11)*SMB0P(20)
       SMB(52) = SMB0M(12)*SMB0P(20)
       SMB(53) = SMB0M(4)*SMB0P(21)
       SMB(54) = SMB0M(4)*SMB0P(22)
       SMB(55) = SMB0M(4)*SMB0P(23)
       SMB(56) = SMB0M(2)*SMB0P(24)
       SMB(57) = SMB0M(2)*SMB0P(23)
       SMB(60) = SMB0M(8)*SMB0P(23)
       SMB(61) = SMB0M(9)*SMB0P(23)
       SMB(64) = SMB0M(11)*SMB0P(23)
       SMB(65) = SMB0M(12)*SMB0P(23)
       SMB(67) = SMB0M(4)*SMB0P(25)
       SMB(79) = SMB0M(6)*SMB0P(20)
       SMB(80) = SMB0M(6)*SMB0P(26)
       SMB(81) = SMB0M(17)*SMB0P(5)
       SMB(82) = SMB0M(17)*SMB0P(10)
       SMB(83) = SMB0M(17)*SMB0P(20)
       SMB(84) = SMB0M(9)*SMB0P(26)
       SMB(85) = SMB0M(9)*SMB0P(21)
       SMB(86) = SMB0M(9)*SMB0P(25)
       SMB(87) = SMB0M(9)*SMB0P(22)
       SMB(88) = SMB0M(17)*SMB0P(21)
       SMB(89) = SMB0M(17)*SMB0P(25)
       SMB(90) = SMB0M(17)*SMB0P(22)
       SMB(91) = SMB0M(2)*SMB0P(27)
       SMB(97) = SMB0M(9)*SMB0P(13)
       SMB(98) = SMB0M(9)*SMB0P(24)
       SMB(99) = SMB0M(6)*SMB0P(10)
       SMB(100) = SMB0M(6)*SMB0P(21)
       SMB(101) = SMB0M(6)*SMB0P(22)
       SMB(102) = SMB0M(28)*SMB0P(1)
       SMB(103) = SMB0M(28)*SMB0P(5)
       SMB(104) = SMB0M(28)*SMB0P(10)
       SMB(105) = SMB0M(18)*SMB0P(10)
       SMB(106) = SMB0M(28)*SMB0P(20)
       SMB(107) = SMB0M(6)*SMB0P(13)
       SMB(108) = SMB0M(6)*SMB0P(24)
       SMB(109) = SMB0M(8)*SMB0P(13)
       SMB(110) = SMB0M(8)*SMB0P(24)
       SMB(111) = SMB0M(28)*SMB0P(13)
       SMB(112) = SMB0M(18)*SMB0P(13)
       SMB(113) = SMB0M(28)*SMB0P(24)
       SMB(114) = SMB0M(18)*SMB0P(24)
       SMB(115) = SMB0M(29)*SMB0P(1)
       SMB(116) = SMB0M(29)*SMB0P(3)
       SMB(117) = SMB0M(29)*SMB0P(5)
       SMB(118) = SMB0M(29)*SMB0P(10)
       SMB(119) = SMB0M(29)*SMB0P(20)
       SMB(120) = SMB0M(8)*SMB0P(7)
       SMB(121) = SMB0M(8)*SMB0P(26)
       SMB(122) = SMB0M(29)*SMB0P(7)
       SMB(123) = SMB0M(29)*SMB0P(14)
       SMB(124) = SMB0M(29)*SMB0P(26)
       SMB(125) = SMB0M(30)*SMB0P(3)
       SMB(127) = SMB0M(30)*SMB0P(5)
       SMB(128) = SMB0M(19)*SMB0P(5)
       SMB(129) = SMB0M(30)*SMB0P(10)
       SMB(130) = SMB0M(19)*SMB0P(10)
       SMB(131) = SMB0M(30)*SMB0P(20)
       SMB(132) = SMB0M(19)*SMB0P(20)
       SMB(133) = SMB0M(8)*SMB0P(21)
       SMB(134) = SMB0M(8)*SMB0P(22)
       SMB(136) = SMB0M(30)*SMB0P(21)
       SMB(137) = SMB0M(19)*SMB0P(21)
       SMB(138) = SMB0M(30)*SMB0P(22)
       SMB(139) = SMB0M(19)*SMB0P(22)
       SMB(142) = SMB0M(6)*SMB0P(25)
       SMB(143) = SMB0M(8)*SMB0P(25)
       SMB(144) = SMB0M(30)*SMB0P(25)
       SMB(145) = SMB0M(16)*SMB0P(5)
       SMB(146) = SMB0M(16)*SMB0P(10)
       SMB(147) = SMB0M(16)*SMB0P(20)
       SMB(148) = SMB0M(6)*SMB0P(27)
       SMB(149) = SMB0M(8)*SMB0P(27)
       SMB(150) = SMB0M(9)*SMB0P(27)
       SMB(152) = SMB0M(16)*SMB0P(13)
       SMB(153) = SMB0M(28)*SMB0P(27)
       SMB(154) = SMB0M(16)*SMB0P(27)
       SMB(190) = SMB0M(16)*SMB0P(24)
       SMB(217) = SMB0M(6)*SMB0P(23)
       SMB(218) = SMB0M(6)*SMB0P(31)
       SMB(219) = SMB0M(8)*SMB0P(31)
       SMB(220) = SMB0M(9)*SMB0P(31)
       SMB(221) = SMB0M(15)*SMB0P(5)
       SMB(222) = SMB0M(15)*SMB0P(20)
       SMB(232) = SMB0M(15)*SMB0P(23)
c************************************************************************************
       SMB(3) = dotcc(SMB1M2,SMB1P1)
       SMB(4) = dotcc(SMB1M2,SMB1P3)
       SMB(11) = dotcc(SMB1M4,SMB1P1)
       SMB(12) = dotcc(SMB1M5,SMB1P1)
       SMB(19) = dotcc(SMB1M2,SMB1P6)
       SMB(22) = dotcc(SMB1M4,SMB1P6)
       SMB(23) = dotcc(SMB1M5,SMB1P6)
       SMB(26) = dotcc(SMB1M7,SMB1P1)
       SMB(29) = dotcc(SMB1M7,SMB1P6)
       SMB(34) = dotcc(SMB1M8,SMB1P1)
       SMB(35) = dotcc(SMB1M8,SMB1P3)
       SMB(38) = dotcc(SMB1M9,SMB1P1)
       SMB(39) = dotcc(SMB1M9,SMB1P3)
       SMB(42) = dotcc(SMB1M10,SMB1P1)
       SMB(58) = dotcc(SMB1M2,SMB1P11)
       SMB(59) = dotcc(SMB1M2,SMB1P12)
       SMB(62) = dotcc(SMB1M10,SMB1P11)
       SMB(63) = dotcc(SMB1M10,SMB1P12)
       SMB(66) = dotcc(SMB1M13,SMB1P1)
       SMB(68) = dotcc(SMB1M2,SMB1P14)
       SMB(69) = dotcc(SMB1M2,SMB1P15)
       SMB(71) = dotcc(SMB1M13,SMB1P14)
       SMB(72) = dotcc(SMB1M13,SMB1P15)
       SMB(73) = dotcc(SMB1M16,SMB1P1)
       SMB(74) = dotcc(SMB1M2,SMB1P17)
       SMB(75) = dotcc(SMB1M16,SMB1P14)
       SMB(76) = dotcc(SMB1M16,SMB1P17)
       SMB(77) = dotcc(SMB1M13,SMB1P17)
       SMB(78) = dotcc(SMB1M16,SMB1P15)
       SMB(92) = dotcc(SMB1M2,SMB1P18)
       SMB(93) = dotcc(SMB1M2,SMB1P19)
       SMB(95) = dotcc(SMB1M5,SMB1P18)
       SMB(96) = dotcc(SMB1M5,SMB1P19)
       SMB(126) = dotcc(SMB1M20,SMB1P1)
       SMB(135) = dotcc(SMB1M20,SMB1P17)
       SMB(140) = dotcc(SMB1M21,SMB1P1)
       SMB(141) = dotcc(SMB1M21,SMB1P3)
       SMB(151) = dotcc(SMB1M4,SMB1P18)
       SMB(155) = dotcc(SMB1M2,SMB1P22)
       SMB(156) = dotcc(SMB1M2,SMB1P23)
       SMB(157) = dotcc(SMB1M2,SMB1P24)
       SMB(161) = dotcc(SMB1M7,SMB1P18)
       SMB(162) = dotcc(SMB1M7,SMB1P19)
       SMB(163) = dotcc(SMB1M4,SMB1P19)
       SMB(167) = dotcc(SMB1M20,SMB1P14)
       SMB(168) = dotcc(SMB1M20,SMB1P15)
       SMB(185) = dotcc(SMB1M2,SMB1P25)
       SMB(191) = dotcc(SMB1M9,SMB1P18)
       SMB(192) = dotcc(SMB1M8,SMB1P18)
       SMB(193) = dotcc(SMB1M9,SMB1P6)
       SMB(194) = dotcc(SMB1M8,SMB1P6)
       SMB(195) = dotcc(SMB1M9,SMB1P19)
       SMB(196) = dotcc(SMB1M8,SMB1P19)
       SMB(197) = dotcc(SMB1M8,SMB1P14)
       SMB(198) = dotcc(SMB1M9,SMB1P17)
       SMB(199) = dotcc(SMB1M8,SMB1P17)
       SMB(200) = dotcc(SMB1M8,SMB1P15)
       SMB(201) = dotcc(SMB1M9,SMB1P11)
       SMB(202) = dotcc(SMB1M8,SMB1P11)
       SMB(203) = dotcc(SMB1M8,SMB1P22)
       SMB(204) = dotcc(SMB1M9,SMB1P23)
       SMB(205) = dotcc(SMB1M8,SMB1P23)
       SMB(206) = dotcc(SMB1M8,SMB1P24)
       SMB(207) = dotcc(SMB1M21,SMB1P18)
       SMB(208) = dotcc(SMB1M21,SMB1P6)
       SMB(209) = dotcc(SMB1M21,SMB1P19)
       SMB(210) = dotcc(SMB1M21,SMB1P14)
       SMB(211) = dotcc(SMB1M21,SMB1P17)
       SMB(212) = dotcc(SMB1M21,SMB1P15)
       SMB(213) = dotcc(SMB1M21,SMB1P11)
       SMB(214) = dotcc(SMB1M21,SMB1P22)
       SMB(215) = dotcc(SMB1M21,SMB1P23)
       SMB(216) = dotcc(SMB1M21,SMB1P24)
       SMB(223) = dotcc(SMB1M2,SMB1P26)
       SMB(226) = dotcc(SMB1M7,SMB1P11)
       SMB(227) = dotcc(SMB1M4,SMB1P11)
       SMB(228) = dotcc(SMB1M5,SMB1P11)
       SMB(229) = dotcc(SMB1M7,SMB1P12)
       SMB(230) = dotcc(SMB1M4,SMB1P12)
       SMB(231) = dotcc(SMB1M5,SMB1P12)
       SMB(234) = dotcc(SMB1M7,SMB1P26)
       SMB(235) = dotcc(SMB1M4,SMB1P26)
       SMB(236) = dotcc(SMB1M5,SMB1P26)
       SMB(237) = dotcc(SMB1M9,SMB1P12)
       SMB(238) = dotcc(SMB1M9,SMB1P25)
       SMB(239) = dotcc(SMB1M8,SMB1P12)
       SMB(240) = dotcc(SMB1M8,SMB1P25)
       SMB(241) = dotcc(SMB1M21,SMB1P12)
       SMB(242) = dotcc(SMB1M21,SMB1P25)
       SMB(245) = dotcc(SMB1M2,SMB1P27)
       SMB(248) = dotcc(SMB1M20,SMB1P11)
       SMB(250) = dotcc(SMB1M20,SMB1P22)
       SMB(251) = dotcc(SMB1M20,SMB1P12)
       SMB(254) = dotcc(SMB1M20,SMB1P27)
       SMB(255) = dotcc(SMB1M16,SMB1P11)
       SMB(256) = dotcc(SMB1M13,SMB1P11)
       SMB(257) = dotcc(SMB1M16,SMB1P22)
       SMB(258) = dotcc(SMB1M13,SMB1P22)
       SMB(259) = dotcc(SMB1M16,SMB1P12)
       SMB(260) = dotcc(SMB1M13,SMB1P12)
       SMB(261) = dotcc(SMB1M16,SMB1P27)
       SMB(262) = dotcc(SMB1M13,SMB1P27)
       SMB(263) = dotcc(SMB1M10,SMB1P22)
c************************************************************************************
       SMB(70) = ten2cc(SMB2M2,SMB2P1)
       SMB(94) = ten2cc(SMB2M4,SMB2P3)
       SMB(158) = ten2cc(SMB2M5,SMB2P3)
       SMB(159) = ten2cc(SMB2M6,SMB2P3)
       SMB(160) = ten2cc(SMB2M7,SMB2P3)
       SMB(164) = ten2cc(SMB2M5,SMB2P1)
       SMB(165) = ten2cc(SMB2M6,SMB2P1)
       SMB(166) = ten2cc(SMB2M7,SMB2P1)
       SMB(170) = ten2cc(SMB2M5,SMB2P8)
       SMB(171) = ten2cc(SMB2M6,SMB2P8)
       SMB(172) = ten2cc(SMB2M7,SMB2P8)
       SMB(173) = ten2cc(SMB2M5,SMB2P9)
       SMB(174) = ten2cc(SMB2M6,SMB2P9)
       SMB(175) = ten2cc(SMB2M7,SMB2P9)
       SMB(176) = ten2cc(SMB2M5,SMB2P10)
       SMB(177) = ten2cc(SMB2M6,SMB2P10)
       SMB(178) = ten2cc(SMB2M7,SMB2P10)
       SMB(180) = ten2cc(SMB2M5,SMB2P11)
       SMB(181) = ten2cc(SMB2M5,SMB2P12)
       SMB(182) = ten2cc(SMB2M5,SMB2P13)
       SMB(183) = ten2cc(SMB2M6,SMB2P11)
       SMB(184) = ten2cc(SMB2M6,SMB2P12)
       SMB(186) = ten2cc(SMB2M6,SMB2P13)
       SMB(187) = ten2cc(SMB2M7,SMB2P11)
       SMB(188) = ten2cc(SMB2M7,SMB2P12)
       SMB(189) = ten2cc(SMB2M7,SMB2P13)
       SMB(224) = ten2cc(SMB2M4,SMB2P8)
       SMB(225) = ten2cc(SMB2M4,SMB2P9)
       SMB(233) = ten2cc(SMB2M4,SMB2P14)
       SMB(243) = ten2cc(SMB2M4,SMB2P10)
       SMB(244) = ten2cc(SMB2M4,SMB2P15)
       SMB(246) = ten2cc(SMB2M2,SMB2P8)
       SMB(247) = ten2cc(SMB2M2,SMB2P9)
       SMB(249) = ten2cc(SMB2M2,SMB2P10)
       SMB(252) = ten2cc(SMB2M2,SMB2P16)
       SMB(253) = ten2cc(SMB2M2,SMB2P17)
       SMB(269) = ten2cc(SMB2M5,SMB2P14)
       SMB(270) = ten2cc(SMB2M5,SMB2P19)
       SMB(271) = ten2cc(SMB2M5,SMB2P15)
       SMB(272) = ten2cc(SMB2M5,SMB2P16)
       SMB(273) = ten2cc(SMB2M5,SMB2P17)
       SMB(274) = ten2cc(SMB2M5,SMB2P20)
       SMB(275) = ten2cc(SMB2M5,SMB2P21)
       SMB(276) = ten2cc(SMB2M6,SMB2P14)
       SMB(277) = ten2cc(SMB2M6,SMB2P19)
       SMB(278) = ten2cc(SMB2M6,SMB2P15)
       SMB(279) = ten2cc(SMB2M6,SMB2P16)
       SMB(280) = ten2cc(SMB2M6,SMB2P17)
       SMB(281) = ten2cc(SMB2M6,SMB2P20)
       SMB(282) = ten2cc(SMB2M6,SMB2P21)
       SMB(283) = ten2cc(SMB2M7,SMB2P14)
       SMB(284) = ten2cc(SMB2M7,SMB2P19)
       SMB(285) = ten2cc(SMB2M7,SMB2P15)
       SMB(286) = ten2cc(SMB2M7,SMB2P16)
       SMB(287) = ten2cc(SMB2M7,SMB2P17)
       SMB(288) = ten2cc(SMB2M7,SMB2P20)
       SMB(289) = ten2cc(SMB2M7,SMB2P21)
       SMB(290) = ten2cc(SMB2M4,SMB2P19)
c************************************************************************************
       SMB(169) = ten3cc(SMB3M2,SMB3P1)
       SMB(179) = ten3cc(SMB3M2,SMB3P3)
       SMB(264) = ten3cc(SMB3M2,SMB3P4)
       SMB(265) = ten3cc(SMB3M2,SMB3P5)
       SMB(266) = ten3cc(SMB3M2,SMB3P7)
       SMB(267) = ten3cc(SMB3M2,SMB3P8)
       SMB(268) = ten3cc(SMB3M2,SMB3P9)
       SMB(291) = ten3cc(SMB3M2,SMB3P10)
       SMB(292) = ten3cc(SMB3M2,SMB3P14)
       SMB(293) = ten3cc(SMB3M2,SMB3P15)
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


c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       result(1) = (F(2)*SMB(4)+F(4)*SMB(21)+F(6)*SMB(31)+mup3mup4*(F(
     &   1)*SMB(3)+8*(F(5)*SMB(30)+F(3)*SMB(41)+F(10)*SMB(61))+4*(F(7
     &   )*SMB(34)-F(8)*SMB(58)+F(9)*SMB(59)-s234*F(11)*SMB(155)))+Fa
     &   (1)*(4*(F(7)*SMB(66)+F(8)*SMB(92)+F(9)*SMB(93))-8*(F(3)*SMB(
     &   8)-F(10)*SMB(150))-2*F(7)*SMB(160))+F(14)*SMB(187)+F(15)*SMB
     &   (188)+s234*F(16)*SMB(189)-4*F(10)*SMB(220)-2*(F(7)*SMB(35)-F
     &   (8)*SMB(156)-s234*F(11)*SMB(157)+F(9)*SMB(185)-F(17)*Fa(1)*S
     &   MB(284)-mup3mup4*(F(12)*SMB(172)+F(7)*SMB(175)+s234*F(13)*SM
     &   B(178)+F(17)*SMB(288)))+F(18)*SMB(289)+F(1454)*(P(33)*SMB(17
     &   9)-SMB(268)-2*(Fa(1)*(SMB(264)+SMB(291))-mup3mup4*(P(22)*SMB
     &   (169)+SMB(267)-SMB(292)))+SMB(293)))/(s23*s234)
       result(2) = (F(21)*SMB(4)+Fa(2)*SMB(7)+Fa(4)*SMB(19)+F(46)*SMB(
     &   84)+F(47)*SMB(91)+8*(Fa(3)*SMB(8)-p1mup4*(2*F(41)*SMB(98)+F(
     &   54)*SMB(110))-F(55)*SMB(154))+F(57)*SMB(156)+F(60)*SMB(183)-
     &   F(62)*SMB(185)+F(58)*SMB(186)+s45*F(63)*SMB(187)-F(61)*(SMB(
     &   184)+SMB(188))+F(53)*(4*SMB(96)+SMB(189))+4*(F(19)*SMB(1)+F(
     &   22)*SMB(20)+F(31)*SMB(21)+F(38)*SMB(56)+Fa(5)*SMB(66)+Fa(6)*
     &   SMB(73)+Fa(7)*SMB(92)+Fa(8)*SMB(93)-s45*F(32)*(SMB(31)+SMB(9
     &   5))+F(27)*(SMB(18)+SMB(120))+F(24)*SMB(121)-F(44)*SMB(151)+F
     &   (59)*SMB(163)+F(64)*SMB(220))+2*(F(30)*(-2*(SMB(22)+SMB(23)-
     &   2*SMB(32))+SMB(35))+F(50)*SMB(157)-Fa(10)*SMB(159)-Fa(11)*SM
     &   B(160)-4*F(41)*SMB(219)-p6mup4*(4*F(42)*SMB(97)+8*F(41)*SMB(
     &   109)+F(55)*(4*SMB(192)+SMB(276)+SMB(283)))+Fa(9)*(8*F(41)*SM
     &   B(149)+4*F(54)*SMB(150)+F(55)*(2*SMB(223)+SMB(277)+SMB(284))
     &   )-F(55)*(2*(SMB(206)+SMB(235)+SMB(236))+p1mup4*(-2*SMB(196)+
     &   SMB(278)+SMB(285))))+mup3mup4*(-4*F(20)*SMB(3)+8*(F(33)*SMB(
     &   40)+F(34)*SMB(41)+F(35)*SMB(44)-F(36)*SMB(49)+F(37)*SMB(50))
     &   +F(39)*SMB(58)+2*(8*F(41)*SMB(60)+4*F(42)*SMB(61)+2*(F(40)*S
     &   MB(59)+F(56)*SMB(155))+F(44)*SMB(171)+s45*F(32)*(4*SMB(30)+S
     &   MB(172))+F(30)*(-2*SMB(34)+SMB(174)+SMB(175))+F(45)*SMB(177)
     &   +F(43)*SMB(178)+F(55)*(4*SMB(203)+SMB(281)+SMB(288))))+F(65)
     &   *(SMB(282)+SMB(289))+F(1453)*(P(52)*(12*SMB(42)+6*SMB(94)-2*
     &   mup3mup4*SMB(169)+SMB(179))-6*(SMB(233)-SMB(244))-SMB(268)+2
     &   *(Fa(12)*SMB(265)+mup3mup4*SMB(267)-Fa(9)*(SMB(264)+SMB(291)
     &   )-mup3mup4*SMB(292))+SMB(293)))/(s23*s45)
       result(3) = (Fa(13)*SMB(3)+F(74)*SMB(4)+F(81)*SMB(28)+F(84)*SMB
     &   (53)+Fa(24)*SMB(58)+Fa(25)*SMB(59)+Fa(26)*SMB(68)+F(78)*(SMB
     &   (18)+SMB(120))+Fa(33)*SMB(140)+Fa(34)*SMB(155)+4*(F(66)*SMB(
     &   1)+F(67)*SMB(2)-Fa(17)*SMB(19)+Fa(27)*SMB(74)+Fa(29)*SMB(92)
     &   +Fa(30)*SMB(93)+F(80)*(P(45)*SMB(80)-2*(p1mup3*SMB(101)-p1mu
     &   p4*SMB(108)))+F(77)*(SMB(54)-SMB(56)+SMB(121))+Fa(14)*(8*F(7
     &   5)*SMB(5)+F(72)*SMB(170)))+Fa(37)*SMB(173)+F(73)*(-8*(SMB(27
     &   )-SMB(102)+SMB(125))+4*(SMB(29)-SMB(135))+2*SMB(141)+SMB(181
     &   ))+F(90)*SMB(182)+8*(Fa(15)*(SMB(7)-SMB(14))+Fa(16)*(-SMB(9)
     &   +SMB(15))-Fa(18)*SMB(24)+Fa(19)*SMB(43)-Fa(20)*SMB(44)+Fa(21
     &   )*SMB(47)-Fa(23)*SMB(49)+F(71)*P(115)*SMB(67)-F(87)*Fa(28)*S
     &   MB(79)+Fa(31)*SMB(99)+F(89)*(2*p6mup4*SMB(107)+SMB(218)))+2*
     &   (8*(p2mup3*F(72)*SMB(57)-p1mup3*F(72)*SMB(106)-p6mup4*F(72)*
     &   SMB(127)+p1mup4*F(72)*SMB(131)+p2mup3*F(72)*SMB(143)-p2mup4*
     &   F(72)*(SMB(55)+SMB(149)))+F(88)*SMB(156)+F(72)*P(115)*SMB(15
     &   7)+Fa(32)*(2*SMB(126)+SMB(158))-p6mup3*(8*(F(89)*SMB(100)+F(
     &   72)*(P(115)*SMB(13)-SMB(103)))+P(45)*(4*F(80)*SMB(25)-F(71)*
     &   (2*SMB(26)+SMB(164))))+Fa(38)*SMB(176)-F(70)*SMB(185)+4*(F(7
     &   2)*P(115)*SMB(91)+F(72)*SMB(144)-F(72)*SMB(153)+F(71)*P(45)*
     &   SMB(222))+2*(F(72)*Fa(22)*(4*SMB(48)+SMB(69))+F(72)*P(45)*SM
     &   B(162)+F(71)*P(45)*SMB(168)-F(72)*(-(p1mup4*SMB(209))+p1mup3
     &   *SMB(212)+SMB(216)+p2mup4*SMB(223)-SMB(234)-p2mup3*SMB(245)-
     &   2*(p6mup3*SMB(210)+mup3mup4*SMB(214)-p6mup3*SMB(226)-p1mup3*
     &   SMB(229)-SMB(232)-p6mup4*(SMB(207)-SMB(248)))))-Fa(9)*(8*F(7
     &   5)*SMB(148)+F(72)*(4*SMB(250)+SMB(270)))-F(72)*(2*SMB(254)-p
     &   6mup4*SMB(269)-p1mup4*(4*SMB(251)+SMB(271))+p6mup3*SMB(272)+
     &   p1mup3*SMB(273))-mup3mup4*(8*F(75)*SMB(217)+F(72)*SMB(274)))
     &   +F(72)*SMB(275)+F(1452)*(48*(SMB(45)+SMB(51))+24*(SMB(62)-SM
     &   B(63)+SMB(75)-SMB(151))+12*(SMB(78)-SMB(163))+P(49)*SMB(179)
     &   +6*(P(45)*(SMB(70)-SMB(94))-SMB(183)-SMB(186)-SMB(233)+SMB(2
     &   44)+SMB(252)-SMB(253))-SMB(268)+2*(Fa(36)*SMB(169)+Fa(12)*(-
     &   6*(SMB(159)-SMB(247))+SMB(265))+Fa(35)*(6*(SMB(165)-SMB(225)
     &   )-SMB(266))-Fa(9)*(6*(SMB(246)+SMB(249))+SMB(264)+SMB(291))+
     &   mup3mup4*(6*(SMB(171)+SMB(177))+SMB(267)-SMB(292)))+SMB(293)
     &   ))/(s345*s45)
       result(4) = (F(91)*SMB(1)+F(93)*SMB(4)+Fa(39)*SMB(7)+Fa(40)*SMB
     &   (8)+Fa(41)*SMB(19)+F(107)*SMB(20)+F(108)*SMB(21)+F(109)*SMB(
     &   22)+F(112)*SMB(31)+Fa(43)*SMB(73)+F(134)*SMB(84)+Fa(44)*SMB(
     &   92)+Fa(45)*SMB(93)+F(160)*SMB(120)+F(161)*SMB(121)+8*(Fa(46)
     &   *SMB(97)-Fa(47)*SMB(98)-Fa(48)*SMB(109)+Fa(49)*SMB(110)+Fa(5
     &   0)*SMB(149)+Fa(51)*SMB(150))+F(172)*SMB(156)+F(173)*SMB(157)
     &   +Fa(52)*SMB(159)+F(189)*SMB(183)+F(190)*SMB(184)+F(191)*SMB(
     &   185)+F(192)*SMB(186)+F(193)*SMB(187)+F(194)*SMB(188)+F(195)*
     &   SMB(189)+4*(F(102)*SMB(18)-F(113)*SMB(32)-F(121)*SMB(56)+Fa(
     &   42)*SMB(66)+F(135)*SMB(91)+F(154)*SMB(152)-F(196)*SMB(190)+F
     &   a(54)*SMB(192)+Fa(55)*SMB(194)+Fa(56)*SMB(196)+F(211)*SMB(21
     &   9)+F(212)*SMB(220)+Fa(57)*SMB(223))+F(238)*SMB(282)+mup3mup4
     &   *(F(92)*SMB(3)+F(111)*SMB(30)+4*F(114)*SMB(34)+F(115)*SMB(40
     &   )+F(116)*SMB(41)+F(119)*SMB(49)+F(122)*SMB(58)+F(123)*SMB(59
     &   )-8*(F(118)*SMB(44)+F(120)*SMB(50)-F(124)*SMB(60)-F(125)*SMB
     &   (61))+F(171)*SMB(155)+F(185)*SMB(174)-2*(F(183)*SMB(171)+F(1
     &   84)*SMB(172)+F(186)*SMB(175)+F(187)*SMB(177)+F(188)*SMB(178)
     &   +2*(F(209)*SMB(202)-F(218)*SMB(239))-F(155)*SMB(281)-F(210)*
     &   (2*SMB(203)+SMB(288))))+F(251)*SMB(289)+2*(F(110)*SMB(23)-F(
     &   114)*SMB(35)-F(117)*(2*SMB(42)+SMB(94))+F(144)*SMB(95)+F(145
     &   )*SMB(96)+F(170)*SMB(151)+Fa(53)*SMB(160)+F(182)*SMB(163)+F(
     &   209)*SMB(205)-F(217)*SMB(233)-F(155)*SMB(235)-F(210)*(SMB(20
     &   6)+SMB(236))-F(218)*SMB(240)+F(219)*SMB(244)-Fa(59)*SMB(265)
     &   +Fa(60)*SMB(276)+Fa(61)*SMB(277)+Fa(62)*SMB(278)+Fa(63)*SMB(
     &   283)+Fa(64)*SMB(284)-Fa(65)*SMB(285)+F(252)*SMB(290)-Fa(58)*
     &   (SMB(264)+SMB(291)))+F(1456)*(P(52)*SMB(179)-SMB(268)-2*mup3
     &   mup4*(P(52)*SMB(169)-SMB(267)+SMB(292))+SMB(293)))/s23
       result(5) = (F(253)*SMB(1)+F(254)*SMB(2)+Fa(66)*SMB(3)+F(272)*S
     &   MB(4)+Fa(67)*SMB(5)+Fa(68)*SMB(7)+Fa(69)*SMB(8)+Fa(70)*SMB(9
     &   )+Fa(71)*SMB(10)+Fa(72)*(SMB(11)+SMB(12))+Fa(73)*SMB(13)+Fa(
     &   75)*SMB(15)+F(309)*SMB(18)+Fa(76)*SMB(19)+F(314)*SMB(20)+F(3
     &   15)*SMB(21)+Fa(77)*SMB(24)+Fa(78)*SMB(25)+F(324)*SMB(27)+F(3
     &   25)*SMB(28)+F(326)*SMB(29)+Fa(80)*SMB(30)+F(332)*SMB(31)+Fa(
     &   81)*SMB(34)+Fa(82)*SMB(38)+Fa(83)*SMB(40)+Fa(84)*SMB(41)+Fa(
     &   85)*SMB(43)-Fa(86)*SMB(44)+Fa(87)*SMB(47)+Fa(88)*SMB(48)+F(3
     &   86)*SMB(54)+F(388)*SMB(56)+Fa(93)*SMB(58)+F(399)*SMB(67)+Fa(
     &   95)*SMB(68)+Fa(96)*SMB(69)+Fa(98)*SMB(74)+Fa(99)*SMB(79)+F(4
     &   20)*SMB(80)+F(425)*SMB(91)+Fa(102)*SMB(92)+Fa(103)*SMB(93)+F
     &   a(106)*SMB(99)+16*Fa(111)*SMB(106)+F(445)*SMB(120)+F(446)*SM
     &   B(125)+F(449)*SMB(135)+Fa(122)*SMB(140)-F(451)*SMB(141)+Fa(1
     &   26)*SMB(155)+F(455)*SMB(156)+F(456)*SMB(157)+Fa(127)*SMB(158
     &   )+Fa(130)*SMB(165)+Fa(131)*SMB(166)+Fa(133)*SMB(170)-Fa(136)
     &   *SMB(173)+Fa(138)*SMB(176)+Fa(139)*SMB(177)+Fa(140)*SMB(178)
     &   +F(472)*SMB(181)+F(473)*SMB(182)+F(475)*SMB(186)+F(474)*(SMB
     &   (184)+SMB(188))+F(477)*SMB(189)+Fa(145)*SMB(212)+s45*(F(276)
     &   *SMB(6)+2*F(367)*(SMB(71)+2*(SMB(46)-SMB(95)))+F(453)*(4*SMB
     &   (151)+SMB(183))+F(476)*SMB(187)+F(274)*(-2*SMB(167)+SMB(180)
     &   +4*(SMB(161)-SMB(221))))+Fa(147)*SMB(229)-8*(Fa(74)*SMB(14)-
     &   Fa(89)*SMB(49)-Fa(90)*SMB(50)+F(397)*SMB(63)+Fa(97)*SMB(73)-
     &   Fa(100)*SMB(85)-Fa(101)*SMB(87)+Fa(104)*SMB(97)-Fa(105)*SMB(
     &   98)+Fa(107)*SMB(100)+Fa(108)*SMB(101)-Fa(112)*SMB(107)+Fa(11
     &   3)*SMB(108)+Fa(114)*SMB(109)-Fa(115)*SMB(110)+Fa(110)*(SMB(1
     &   04)+SMB(105)+SMB(116))+Fa(117)*SMB(127)-Fa(118)*(SMB(129)+SM
     &   B(130))+Fa(119)*SMB(131)-Fa(120)*SMB(133)-Fa(121)*SMB(134)+F
     &   a(92)*(SMB(57)+SMB(143))-Fa(125)*SMB(149)-Fa(124)*(-2*SMB(14
     &   7)+SMB(230)+SMB(231)))+4*(F(354)*SMB(42)-F(383)*SMB(51)+F(38
     &   5)*SMB(53)+Fa(94)*SMB(59)-F(284)*SMB(84)-F(280)*SMB(102)+F(2
     &   79)*SMB(121)-F(431)*(SMB(111)+SMB(112)-SMB(122))+F(382)*(SMB
     &   (113)+SMB(114)+SMB(124))+F(274)*(P(45)*SMB(37)-2*Fa(9)*SMB(1
     &   28)+SMB(136)+SMB(137))+F(459)*SMB(163)-Fa(142)*SMB(200)-Fa(1
     &   43)*SMB(207)+Fa(116)*(2*SMB(115)+SMB(193)+SMB(208))-Fa(144)*
     &   SMB(209)+Fa(141)*(SMB(198)+SMB(211))+F(442)*SMB(218)+F(479)*
     &   SMB(219)-F(432)*SMB(220)+Fa(91)*(2*SMB(55)+SMB(223))+F(275)*
     &   (SMB(144)-2*SMB(153)-SMB(216)-SMB(232)+SMB(234))+mup3mup4*(F
     &   (431)*(-2*SMB(117)+SMB(201)+SMB(213))+2*(F(394)*SMB(60)+F(39
     &   5)*SMB(61)+F(387)*SMB(203)+F(275)*SMB(214)-F(443)*SMB(217))-
     &   F(382)*(2*SMB(119)+SMB(237)+SMB(241)))-Fa(148)*SMB(245)+Fa(1
     &   49)*SMB(247)+p6mup4*F(275)*SMB(248)+F(387)*(SMB(64)+SMB(65)+
     &   SMB(89)-2*SMB(154)+p1mup4*(2*SMB(83)+SMB(196))-SMB(206)-SMB(
     &   235)-SMB(236)-p6mup4*(2*(SMB(81)+SMB(192))+SMB(255)+SMB(256)
     &   )))+Fa(150)*SMB(273)+F(275)*SMB(275)+Fa(151)*(SMB(280)+SMB(2
     &   87))-2*(-(F(313)*(2*(SMB(22)+SMB(23)-2*SMB(32))-SMB(35)))+F(
     &   283)*(2*SMB(36)+SMB(39))-F(407)*SMB(70)-F(366)*(2*SMB(45)+SM
     &   B(75))+F(304)*(2*(SMB(16)+SMB(17)-SMB(33))+SMB(76)+SMB(77))-
     &   F(414)*SMB(78)+F(384)*(SMB(72)+2*(SMB(52)-SMB(96)))-Fa(128)*
     &   SMB(159)-Fa(129)*SMB(160)-Fa(79)*(2*SMB(26)+SMB(164))+Fa(132
     &   )*SMB(169)-Fa(134)*SMB(171)-Fa(135)*SMB(172)+Fa(137)*(SMB(17
     &   4)+SMB(175))-F(311)*SMB(185)+F(431)*(SMB(204)+SMB(215))-F(37
     &   0)*(-SMB(168)+2*(SMB(162)-SMB(222)))-F(382)*(SMB(238)+SMB(24
     &   2))-F(396)*(4*SMB(62)-2*Fa(9)*SMB(246)+SMB(252))+F(397)*SMB(
     &   253)+F(275)*SMB(254)+2*(F(448)*(SMB(138)+SMB(139))+Fa(9)*(2*
     &   (F(398)*SMB(66)+F(447)*SMB(126)+F(448)*SMB(132)-F(452)*SMB(1
     &   48)-F(394)*SMB(150))+F(431)*SMB(191)+F(382)*SMB(195))-p1mup4
     &   *F(275)*SMB(251)+Fa(9)*(F(397)*SMB(249)+F(275)*SMB(250)-F(38
     &   7)*(SMB(257)+SMB(258)))+p1mup4*F(387)*SMB(259)+p1mup4*F(387)
     &   *SMB(260))-F(387)*SMB(261)-F(387)*SMB(262)-Fa(146)*(6*SMB(22
     &   5)+SMB(266))-p6mup4*F(275)*SMB(269)-p1mup4*F(275)*SMB(271)-F
     &   a(109)*(8*SMB(103)+4*(SMB(210)-SMB(226))-SMB(272))+mup3mup4*
     &   F(275)*SMB(274)+p6mup4*F(387)*SMB(276)+p1mup4*F(387)*SMB(278
     &   )-mup3mup4*F(387)*SMB(281)+p6mup4*F(387)*SMB(283)+Fa(9)*(F(2
     &   75)*SMB(270)-F(387)*(SMB(277)+SMB(284)))+p1mup4*F(387)*SMB(2
     &   85)+Fa(123)*(8*SMB(145)+4*(SMB(197)+SMB(227)+SMB(228))+SMB(2
     &   79)+SMB(286))-mup3mup4*F(387)*SMB(288))+F(483)*(SMB(282)+SMB
     &   (289))+F(1455)*(P(52)*(6*SMB(94)+SMB(179))-6*(SMB(233)-SMB(2
     &   44))-SMB(268)+2*(Fa(12)*SMB(265)+mup3mup4*SMB(267)-Fa(9)*(SM
     &   B(264)+SMB(291))-mup3mup4*SMB(292))+SMB(293)))/s45
       result(6) = Fa(152)*SMB(3)+Fa(160)*SMB(13)+Fa(163)*SMB(19)+Fa(1
     &   65)*SMB(25)+F(657)*SMB(36)+Fa(169)*SMB(38)+F(685)*SMB(39)+Fa
     &   (172)*SMB(43)+Fa(174)*SMB(47)+F(817)*SMB(54)+Fa(180)*SMB(58)
     &   +Fa(181)*SMB(59)+Fa(183)*SMB(68)+Fa(184)*SMB(69)+F(896)*SMB(
     &   71)+16*(Fa(158)*SMB(11)+mup3mup4*F(875)*SMB(60)+Fa(185)*SMB(
     &   73))+Fa(186)*SMB(74)+F(906)*SMB(75)+F(934)*SMB(80)+Fa(193)*S
     &   MB(92)+Fa(194)*SMB(93)+F(969)*SMB(94)+F(970)*SMB(95)+F(1045)
     &   *SMB(120)+F(1046)*SMB(121)+F(1070)*SMB(135)+Fa(212)*SMB(140)
     &   +F(1083)*SMB(141)+F(1104)*SMB(151)+Fa(218)*SMB(155)+F(1131)*
     &   SMB(156)+Fa(219)*SMB(158)+Fa(220)*SMB(159)+Fa(221)*SMB(160)+
     &   F(1148)*SMB(161)+Fa(222)*SMB(164)+Fa(223)*SMB(165)+Fa(224)*S
     &   MB(166)+F(1166)*SMB(167)+F(1167)*SMB(168)+Fa(225)*SMB(169)+F
     &   a(226)*SMB(170)+Fa(227)*SMB(171)+Fa(228)*SMB(172)+Fa(229)*SM
     &   B(173)+Fa(230)*SMB(174)+Fa(231)*SMB(175)+Fa(232)*SMB(176)+Fa
     &   (233)*SMB(177)+Fa(234)*SMB(178)+F(1311)*SMB(179)+F(1312)*SMB
     &   (180)+F(1313)*SMB(181)+F(1314)*SMB(182)+F(1315)*SMB(183)+F(1
     &   316)*SMB(184)+F(1317)*SMB(185)+F(1318)*SMB(186)+F(1319)*SMB(
     &   187)+F(1320)*SMB(188)+F(1321)*SMB(189)+8*(Fa(153)*SMB(5)-Fa(
     &   154)*SMB(7)+Fa(155)*SMB(8)+Fa(156)*SMB(9)+Fa(157)*SMB(10)-Fa
     &   (159)*SMB(12)+Fa(161)*SMB(14)-Fa(162)*SMB(15)+F(576)*SMB(17)
     &   -Fa(164)*SMB(24)+Fa(166)*SMB(26)+Fa(167)*SMB(30)+Fa(170)*SMB
     &   (40)+Fa(171)*SMB(41)+F(735)*SMB(42)+Fa(173)*SMB(44)+F(761)*S
     &   MB(46)+Fa(175)*SMB(48)+Fa(176)*SMB(49)+Fa(177)*SMB(50)-F(146
     &   1)*SMB(52)+Fa(178)*SMB(55)+Fa(179)*SMB(57)-Fa(182)*SMB(66)+F
     &   a(187)*SMB(79)+Fa(188)*SMB(81)+Fa(189)*SMB(82)-Fa(190)*SMB(8
     &   3)+Fa(191)*SMB(85)-p4mup3*F(952)*SMB(86)+Fa(192)*SMB(87)-Fa(
     &   195)*SMB(97)+Fa(196)*SMB(98)+Fa(197)*SMB(99)+Fa(198)*SMB(100
     &   )+Fa(199)*SMB(101)-Fa(201)*SMB(107)-Fa(202)*SMB(108)-Fa(203)
     &   *SMB(109)+Fa(204)*SMB(110)-Fa(205)*SMB(115)+Fa(200)*(SMB(105
     &   )+SMB(116))-Fa(206)*SMB(126)-Fa(207)*SMB(128)-Fa(208)*SMB(13
     &   0)+Fa(209)*SMB(132)+Fa(210)*SMB(133)+Fa(211)*SMB(134)+p4mup3
     &   *F(1084)*SMB(142)+Fa(213)*SMB(143)+Fa(214)*SMB(146)-Fa(215)*
     &   SMB(148)+Fa(216)*SMB(149)+Fa(217)*SMB(150)+mup3mup4*(F(876)*
     &   SMB(61)-F(1037)*SMB(117)+F(1044)*SMB(118)+F(1038)*SMB(119)+F
     &   (1366)*SMB(217)))-12*F(875)*SMB(219)+Fa(246)*SMB(223)+Fa(252
     &   )*SMB(230)+Fa(253)*SMB(231)+Fa(255)*SMB(245)+4*(F(484)*SMB(1
     &   )-F(485)*SMB(2)+F(531)*SMB(6)+F(577)*SMB(18)+F(583)*SMB(20)+
     &   F(584)*SMB(21)-F(601)*SMB(28)+F(628)*SMB(31)+F(629)*SMB(32)-
     &   F(630)*SMB(33)-Fa(168)*SMB(34)-F(658)*SMB(37)-F(816)*SMB(53)
     &   +F(823)*SMB(56)+F(877)*SMB(62)-F(878)*SMB(63)+F(884)*SMB(67)
     &   +F(946)*SMB(84)+F(958)*SMB(91)-F(1044)*SMB(123)+F(1037)*(-SM
     &   B(112)+SMB(122)+SMB(137))-F(1038)*(SMB(114)+SMB(124)-SMB(139
     &   ))-Fa(235)*SMB(191)+Fa(236)*SMB(192)+Fa(237)*SMB(193)+Fa(189
     &   )*SMB(194)-Fa(238)*SMB(195)-Fa(190)*SMB(196)+Fa(240)*SMB(198
     &   )+Fa(214)*SMB(199)-F(956)*(SMB(88)-SMB(152)+mup3mup4*SMB(202
     &   ))+Fa(241)*SMB(207)+Fa(242)*SMB(208)-Fa(243)*SMB(209)-Fa(244
     &   )*SMB(210)+Fa(245)*SMB(211)+F(1367)*SMB(218)+F(1368)*SMB(220
     &   )-Fa(247)*SMB(224)+Fa(248)*SMB(225)+Fa(249)*SMB(226)-Fa(250)
     &   *SMB(227)-Fa(239)*(SMB(197)+SMB(228))+Fa(251)*SMB(229)-F(957
     &   )*(SMB(90)-SMB(190)+mup3mup4*SMB(239))+mup3mup4*(F(1345)*SMB
     &   (201)+F(952)*SMB(203)+F(1365)*SMB(213)+F(1084)*SMB(214)-F(13
     &   90)*SMB(237)-F(1391)*SMB(241))-Fa(254)*SMB(243)-Fa(256)*SMB(
     &   246)+Fa(257)*SMB(247)-Fa(258)*SMB(248)-Fa(259)*SMB(249)-Fa(2
     &   60)*SMB(250)-Fa(261)*SMB(251)-Fa(262)*SMB(255)-Fa(263)*SMB(2
     &   56)+Fa(264)*SMB(257)+Fa(265)*SMB(258)+Fa(266)*SMB(259)+Fa(26
     &   7)*SMB(260)-F(1443)*SMB(263))-F(1459)*SMB(268)+F(1084)*SMB(2
     &   75)+Fa(271)*SMB(280)+F(1449)*SMB(282)+Fa(272)*SMB(287)+F(145
     &   1)*SMB(289)+2*(F(512)*SMB(4)+F(585)*SMB(22)+F(586)*SMB(23)+F
     &   (602)*SMB(29)+F(656)*SMB(35)+F(895)*SMB(70)+F(897)*SMB(72)+F
     &   (907)*SMB(76)-F(908)*SMB(77)+F(909)*SMB(78)+F(971)*SMB(96)+F
     &   (1132)*SMB(157)-F(1149)*SMB(162)-F(1150)*SMB(163)-F(1345)*SM
     &   B(204)+F(956)*SMB(205)-F(1365)*SMB(215)-F(1389)*SMB(233)+F(1
     &   390)*SMB(238)+F(957)*SMB(240)+F(1391)*SMB(242)+F(1392)*SMB(2
     &   44)+F(1422)*SMB(252)-F(1423)*SMB(253)-Fa(269)*SMB(265)+Fa(27
     &   0)*SMB(266)-Fa(258)*SMB(269)-Fa(260)*SMB(270)-Fa(261)*SMB(27
     &   1)+Fa(249)*SMB(272)+Fa(251)*SMB(273)-F(1084)*(SMB(216)-SMB(2
     &   34)+SMB(254)+mup3mup4*SMB(274))-Fa(262)*SMB(276)+Fa(264)*SMB
     &   (277)+Fa(266)*SMB(278)-Fa(250)*SMB(279)+F(875)*(-SMB(235)+SM
     &   B(261)+mup3mup4*SMB(281))-Fa(263)*SMB(283)+Fa(265)*SMB(284)+
     &   Fa(267)*SMB(285)-Fa(239)*SMB(286)+F(952)*(-SMB(206)-SMB(236)
     &   +SMB(262)+mup3mup4*SMB(288))+F(1443)*SMB(290)-Fa(268)*(SMB(2
     &   64)+SMB(291))+mup3mup4*F(1459)*(SMB(267)-SMB(292)))+F(1459)*
     &   SMB(293)
       If (ngluon.eq.0) then
       resultn=-(result(1)+result(2)+result(3)+result(4)+result(5)+res
     &   ult(6))
       elseif (ngluon.eq.1) then
       If (posgluon.eq.3) then
       result(1)=-(result(1)+result(2)+result(3)+result(5)+result(6))
       result(2)=-(result(4))
       elseIf (posgluon.eq.4) then
       result(1)=-(result(1)+result(2)+result(3)+result(4)+result(6))
       result(2)=-(result(5))
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
       Return
       End
