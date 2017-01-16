c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:5/2/2010
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
       subroutine CheckHexBox77(p1,p2,p3,p4,p5,p6,barpsi_p1,psi_p6,bar
     &   psi_p2,psi_p5,mup3,mup4,musq,comp,ngluon,posgluon,result,res
     &   ultn,SMB)
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

      Complex*16 barpsi_p1_P(2),psi_p6_P(2),barpsi_p2_P(2),psi_p5_P(2) 
       Complex*16 barpsi_p1_M(2),psi_p6_M(2),barpsi_p2_M(2),psi_p5_M(2)
       Complex*16 barpsi_p1(4),psi_p6(4),barpsi_p2(4),psi_p5(4)
       Real*8 P(1362)
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
       call  HexBoxFa1(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(1))
       call  HexBoxFa2(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(13))
       call  HexBoxFa3(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(25))
       call  HexBoxFa4(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(37))
       call  HexBoxFa5(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(49))
       call  HexBoxFa6(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(61))
       call  HexBoxFa7(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(73))
       call  HexBoxFa8(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(85))
       call  HexBoxFa9(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mup
     &   3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(97))
       call  HexBoxFa10(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(109))
       call  HexBoxFa11(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(121))
       call  HexBoxFa12(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(133))
       call  HexBoxFa13(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(145))
       call  HexBoxFa14(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(157))
       call  HexBoxFa15(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(169))
       call  HexBoxFa16(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(181))
       call  HexBoxFa17(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(193))
       call  HexBoxFa18(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(205))
       call  HexBoxFa19(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(217))
       call  HexBoxFa20(p1mup3,p1mup4,p2mup3,p2mup4,p3mup3,p3mup4,p4mu
     &   p3,p4mup4,p5mup3,p5mup4,p6mup3,p6mup4,mup3mup4,Fa(229))
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
       SMB0P(1) = SC1c(barpsi_p1_P,mup3,psi_p6_P,+1)
       SMB0P(2) = SC1c(barpsi_p2_P,mup4,psi_p5_P,+1)
       SMB0P(3) = SC1c(barpsi_p1_P,mup4,psi_p6_P,+1)
       SMB0P(4) = SC1c(barpsi_p2_P,mup3,psi_p5_P,+1)
       SMB0P(5) = SC1r(barpsi_p1_P,p5,psi_p6_P,+1)
       SMB0P(6) = SC1r(barpsi_p2_P,p1,psi_p5_P,+1)
       SMB0P(7) = SC1r(barpsi_p2_P,p4,psi_p5_P,+1)
       SMB0P(8) = SC3ccr(barpsi_p2_P,mup4,mup3,p1,psi_p5_P,+1)
       SMB0P(9) = SC3ccr(barpsi_p2_P,mup4,mup3,p4,psi_p5_P,+1)
       SMB0P(10) = SC3crr(barpsi_p1_P,mup4,p5,p2,psi_p6_P,+1)
       SMB0P(11) = SC3crr(barpsi_p1_P,mup4,p5,p4,psi_p6_P,+1)
       SMB0P(12) = SC3crr(barpsi_p1_P,mup3,p5,p2,psi_p6_P,+1)
       SMB0P(13) = SC3crr(barpsi_p1_P,mup3,p5,p4,psi_p6_P,+1)
       SMB0P(14) = SC1r(barpsi_p1_P,p2,psi_p6_P,+1)
       SMB0P(15) = SC1r(barpsi_p2_P,p6,psi_p5_P,+1)
       SMB0P(16) = SC3ccr(barpsi_p2_P,mup4,mup3,p6,psi_p5_P,+1)
       SMB0P(17) = SC1r(barpsi_p1_P,p4,psi_p6_P,+1)
       SMB0P(18) = SC3rrr(barpsi_p1_P,p5,p4,p2,psi_p6_P,+1)
       SMB0P(19) = SC3crr(barpsi_p2_P,mup4,p6,p1,psi_p5_P,+1)
       SMB0P(20) = SC3crr(barpsi_p2_P,mup4,p6,p4,psi_p5_P,+1)
       SMB0P(21) = SC3crr(barpsi_p2_P,mup3,p6,p1,psi_p5_P,+1)
       SMB0P(22) = SC3crr(barpsi_p2_P,mup3,p6,p4,psi_p5_P,+1)
       SMB0P(23) = SC3crr(barpsi_p2_P,mup3,p4,p1,psi_p5_P,+1)
       SMB0P(24) = SC3crr(barpsi_p1_P,mup4,p4,p2,psi_p6_P,+1)
       SMB0P(25) = SC3crr(barpsi_p1_P,mup3,p4,p2,psi_p6_P,+1)
       SMB0P(26) = SC3ccr(barpsi_p1_P,mup4,mup3,p2,psi_p6_P,+1)
       SMB0P(27) = SC3ccr(barpsi_p1_P,mup4,mup3,p4,psi_p6_P,+1)
       SMB0P(28) = SC3ccr(barpsi_p1_P,mup4,mup3,p5,psi_p6_P,+1)
       SMB0P(29) = SC3crr(barpsi_p2_P,mup4,p4,p1,psi_p5_P,+1)
       SMB0P(30) = SC3rrr(barpsi_p2_P,p6,p4,p1,psi_p5_P,+1)
       SMB0P(31) = SC5ccrrr(barpsi_p2_P,mup4,mup3,p6,p4,p1,psi_p5_P,+1)
       SMB0M(1) = SC1c(barpsi_p1_M,mup3,psi_p6_M,-1)
       SMB0M(2) = SC1c(barpsi_p2_M,mup4,psi_p5_M,-1)
       SMB0M(3) = SC1c(barpsi_p1_M,mup4,psi_p6_M,-1)
       SMB0M(4) = SC1c(barpsi_p2_M,mup3,psi_p5_M,-1)
       SMB0M(5) = SC1r(barpsi_p1_M,p5,psi_p6_M,-1)
       SMB0M(6) = SC1r(barpsi_p2_M,p1,psi_p5_M,-1)
       SMB0M(7) = SC1r(barpsi_p2_M,p4,psi_p5_M,-1)
       SMB0M(8) = SC3ccr(barpsi_p2_M,mup4,mup3,p1,psi_p5_M,-1)
       SMB0M(9) = SC3ccr(barpsi_p2_M,mup4,mup3,p4,psi_p5_M,-1)
       SMB0M(10) = SC3crr(barpsi_p1_M,mup4,p5,p2,psi_p6_M,-1)
       SMB0M(11) = SC3crr(barpsi_p1_M,mup4,p5,p4,psi_p6_M,-1)
       SMB0M(12) = SC3crr(barpsi_p1_M,mup3,p5,p2,psi_p6_M,-1)
       SMB0M(13) = SC3crr(barpsi_p1_M,mup3,p5,p4,psi_p6_M,-1)
       SMB0M(14) = SC1r(barpsi_p1_M,p2,psi_p6_M,-1)
       SMB0M(15) = SC1r(barpsi_p2_M,p6,psi_p5_M,-1)
       SMB0M(16) = SC3ccr(barpsi_p2_M,mup4,mup3,p6,psi_p5_M,-1)
       SMB0M(17) = SC1r(barpsi_p1_M,p4,psi_p6_M,-1)
       SMB0M(18) = SC3rrr(barpsi_p1_M,p5,p4,p2,psi_p6_M,-1)
       SMB0M(19) = SC3crr(barpsi_p2_M,mup4,p6,p1,psi_p5_M,-1)
       SMB0M(20) = SC3crr(barpsi_p2_M,mup4,p6,p4,psi_p5_M,-1)
       SMB0M(21) = SC3crr(barpsi_p2_M,mup3,p6,p1,psi_p5_M,-1)
       SMB0M(22) = SC3crr(barpsi_p2_M,mup3,p6,p4,psi_p5_M,-1)
       SMB0M(23) = SC3crr(barpsi_p2_M,mup3,p4,p1,psi_p5_M,-1)
       SMB0M(24) = SC3crr(barpsi_p1_M,mup4,p4,p2,psi_p6_M,-1)
       SMB0M(25) = SC3crr(barpsi_p1_M,mup3,p4,p2,psi_p6_M,-1)
       SMB0M(26) = SC3ccr(barpsi_p1_M,mup4,mup3,p2,psi_p6_M,-1)
       SMB0M(27) = SC3ccr(barpsi_p1_M,mup4,mup3,p4,psi_p6_M,-1)
       SMB0M(28) = SC3ccr(barpsi_p1_M,mup4,mup3,p5,psi_p6_M,-1)
       SMB0M(29) = SC3crr(barpsi_p2_M,mup4,p4,p1,psi_p5_M,-1)
       SMB0M(30) = SC3rrr(barpsi_p2_M,p6,p4,p1,psi_p5_M,-1)
       SMB0M(31) = SC5ccrrr(barpsi_p2_M,mup4,mup3,p6,p4,p1,psi_p5_M,-1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB1P1(i) = SC1c(barpsi_p1_P,v1,psi_p6_P,+1)
       SMB1P2(i) = SC1c(barpsi_p2_P,v1,psi_p5_P,+1)
       SMB1P3(i) = SC3ccc(barpsi_p2_P,mup4,mup3,v1,psi_p5_P,+1)
       SMB1P4(i) = SC3rrc(barpsi_p1_P,p5,p2,v1,psi_p6_P,+1)
       SMB1P5(i) = SC3rrc(barpsi_p1_P,p5,p4,v1,psi_p6_P,+1)
       SMB1P6(i) = SC3crc(barpsi_p2_P,mup3,p1,v1,psi_p5_P,+1)
       SMB1P7(i) = SC3crc(barpsi_p2_P,mup3,p4,v1,psi_p5_P,+1)
       SMB1P8(i) = SC3crc(barpsi_p2_P,mup3,p6,v1,psi_p5_P,+1)
       SMB1P9(i) = SC3crc(barpsi_p2_P,mup4,p1,v1,psi_p5_P,+1)
       SMB1P10(i) = SC3crc(barpsi_p2_P,mup4,p4,v1,psi_p5_P,+1)
       SMB1P11(i) = SC3crc(barpsi_p2_P,mup4,p6,v1,psi_p5_P,+1)
       SMB1P12(i) = SC3rrc(barpsi_p2_P,p6,p1,v1,psi_p5_P,+1)
       SMB1P13(i) = SC3rrc(barpsi_p2_P,p6,p4,v1,psi_p5_P,+1)
       SMB1P14(i) = SC5ccrrc(barpsi_p2_P,mup4,mup3,p6,p1,v1,psi_p5_P,+1)
       SMB1P15(i) = SC5ccrrc(barpsi_p2_P,mup4,mup3,p6,p4,v1,psi_p5_P,+1)
       SMB1P16(i) = SC3crc(barpsi_p1_P,mup3,p2,v1,psi_p6_P,+1)
       SMB1P17(i) = SC3crc(barpsi_p1_P,mup3,p4,v1,psi_p6_P,+1)
       SMB1P18(i) = SC3crc(barpsi_p1_P,mup3,p5,v1,psi_p6_P,+1)
       SMB1P19(i) = SC3crc(barpsi_p1_P,mup4,p2,v1,psi_p6_P,+1)
       SMB1P20(i) = SC3crc(barpsi_p1_P,mup4,p4,v1,psi_p6_P,+1)
       SMB1P21(i) = SC3crc(barpsi_p1_P,mup4,p5,v1,psi_p6_P,+1)
       SMB1P22(i) = SC3ccc(barpsi_p1_P,mup4,mup3,v1,psi_p6_P,+1)
       SMB1P23(i) = SC3rrc(barpsi_p1_P,p4,p2,v1,psi_p6_P,+1)
       SMB1P24(i) = SC3rrc(barpsi_p2_P,p4,p1,v1,psi_p5_P,+1)
       SMB1P25(i) = SC5ccrrc(barpsi_p2_P,mup4,mup3,p4,p1,v1,psi_p5_P,+1)
       SMB1P26(i) = SC5crrrc(barpsi_p2_P,mup3,p6,p4,p1,v1,psi_p5_P,+1)
       SMB1P27(i) = SC5crrrc(barpsi_p2_P,mup4,p6,p4,p1,v1,psi_p5_P,+1)
       SMB1M1(i) = SC1c(barpsi_p1_M,v1,psi_p6_M,-1)
       SMB1M2(i) = SC1c(barpsi_p2_M,v1,psi_p5_M,-1)
       SMB1M3(i) = SC3ccc(barpsi_p2_M,mup4,mup3,v1,psi_p5_M,-1)
       SMB1M4(i) = SC3rrc(barpsi_p1_M,p5,p2,v1,psi_p6_M,-1)
       SMB1M5(i) = SC3rrc(barpsi_p1_M,p5,p4,v1,psi_p6_M,-1)
       SMB1M6(i) = SC3crc(barpsi_p2_M,mup3,p1,v1,psi_p5_M,-1)
       SMB1M7(i) = SC3crc(barpsi_p2_M,mup3,p4,v1,psi_p5_M,-1)
       SMB1M8(i) = SC3crc(barpsi_p2_M,mup3,p6,v1,psi_p5_M,-1)
       SMB1M9(i) = SC3crc(barpsi_p2_M,mup4,p1,v1,psi_p5_M,-1)
       SMB1M10(i) = SC3crc(barpsi_p2_M,mup4,p4,v1,psi_p5_M,-1)
       SMB1M11(i) = SC3crc(barpsi_p2_M,mup4,p6,v1,psi_p5_M,-1)
       SMB1M12(i) = SC3rrc(barpsi_p2_M,p6,p1,v1,psi_p5_M,-1)
       SMB1M13(i) = SC3rrc(barpsi_p2_M,p6,p4,v1,psi_p5_M,-1)
       SMB1M14(i) = SC5ccrrc(barpsi_p2_M,mup4,mup3,p6,p1,v1,psi_p5_M,-1)
       SMB1M15(i) = SC5ccrrc(barpsi_p2_M,mup4,mup3,p6,p4,v1,psi_p5_M,-1)
       SMB1M16(i) = SC3crc(barpsi_p1_M,mup3,p2,v1,psi_p6_M,-1)
       SMB1M17(i) = SC3crc(barpsi_p1_M,mup3,p4,v1,psi_p6_M,-1)
       SMB1M18(i) = SC3crc(barpsi_p1_M,mup3,p5,v1,psi_p6_M,-1)
       SMB1M19(i) = SC3crc(barpsi_p1_M,mup4,p2,v1,psi_p6_M,-1)
       SMB1M20(i) = SC3crc(barpsi_p1_M,mup4,p4,v1,psi_p6_M,-1)
       SMB1M21(i) = SC3crc(barpsi_p1_M,mup4,p5,v1,psi_p6_M,-1)
       SMB1M22(i) = SC3ccc(barpsi_p1_M,mup4,mup3,v1,psi_p6_M,-1)
       SMB1M23(i) = SC3rrc(barpsi_p1_M,p4,p2,v1,psi_p6_M,-1)
       SMB1M24(i) = SC3rrc(barpsi_p2_M,p4,p1,v1,psi_p5_M,-1)
       SMB1M25(i) = SC5ccrrc(barpsi_p2_M,mup4,mup3,p4,p1,v1,psi_p5_M,-1)
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
       SMB2P1(i,j) = SC3ccc(barpsi_p1_P,v1,mup3,v2,psi_p6_P,+1)
       SMB2P2(i,j) = SC3ccc(barpsi_p2_P,v1,mup4,v2,psi_p5_P,+1)
       SMB2P3(i,j) = SC3ccc(barpsi_p1_P,v1,mup4,v2,psi_p6_P,+1)
       SMB2P4(i,j) = SC3ccc(barpsi_p2_P,v1,mup3,v2,psi_p5_P,+1)
       SMB2P5(i,j) = SC3crc(barpsi_p1_P,v1,p2,v2,psi_p6_P,+1)
       SMB2P6(i,j) = SC3crc(barpsi_p2_P,v1,p1,v2,psi_p5_P,+1)
       SMB2P7(i,j) = SC3crc(barpsi_p2_P,v1,p4,v2,psi_p5_P,+1)
       SMB2P8(i,j) = SC3crc(barpsi_p2_P,v1,p6,v2,psi_p5_P,+1)
       SMB2P9(i,j) = SC5cccrc(barpsi_p2_P,mup4,mup3,v1,p1,v2,psi_p5_P,+1)
       SMB2P10(i,j) = SC5cccrc(barpsi_p2_P,mup4,mup3,v1,p4,v2,psi_p5_P,+1)
       SMB2P11(i,j) = SC5cccrc(barpsi_p2_P,mup4,mup3,v1,p6,v2,psi_p5_P,+1)
       SMB2P12(i,j) = SC3crc(barpsi_p1_P,v1,p4,v2,psi_p6_P,+1)
       SMB2P13(i,j) = SC3crc(barpsi_p1_P,v1,p5,v2,psi_p6_P,+1)
       SMB2P14(i,j) = SC5crcrc(barpsi_p2_P,mup4,p6,v1,p1,v2,psi_p5_P,+1)
       SMB2P15(i,j) = SC5crcrc(barpsi_p2_P,mup3,p6,v1,p1,v2,psi_p5_P,+1)
       SMB2P16(i,j) = SC5crcrc(barpsi_p2_P,mup3,p6,v1,p4,v2,psi_p5_P,+1)
       SMB2P17(i,j) = SC5rrcrc(barpsi_p2_P,p6,p4,v1,p1,v2,psi_p5_P,+1)
       SMB2P18(i,j) = SC7ccrrcrc(barpsi_p2_P,mup4,mup3,p6,p4,v1,p1,v2,psi_p5_P,+1)
       SMB2P19(i,j) = SC5crcrc(barpsi_p2_P,mup4,p4,v1,p1,v2,psi_p5_P,+1)
       SMB2P20(i,j) = SC5crcrc(barpsi_p2_P,mup3,p4,v1,p1,v2,psi_p5_P,+1)
       SMB2P21(i,j) = SC5crcrc(barpsi_p2_P,mup4,p6,v1,p4,v2,psi_p5_P,+1)
       SMB2M1(i,j) = SC3ccc(barpsi_p1_M,v1,mup3,v2,psi_p6_M,-1)
       SMB2M2(i,j) = SC3ccc(barpsi_p2_M,v1,mup4,v2,psi_p5_M,-1)
       SMB2M3(i,j) = SC3ccc(barpsi_p1_M,v1,mup4,v2,psi_p6_M,-1)
       SMB2M4(i,j) = SC3ccc(barpsi_p2_M,v1,mup3,v2,psi_p5_M,-1)
       SMB2M5(i,j) = SC3crc(barpsi_p1_M,v1,p2,v2,psi_p6_M,-1)
       SMB2M6(i,j) = SC3crc(barpsi_p2_M,v1,p1,v2,psi_p5_M,-1)
       SMB2M7(i,j) = SC3crc(barpsi_p2_M,v1,p4,v2,psi_p5_M,-1)
       SMB2M8(i,j) = SC3crc(barpsi_p2_M,v1,p6,v2,psi_p5_M,-1)
       SMB2M9(i,j) = SC5cccrc(barpsi_p2_M,mup4,mup3,v1,p1,v2,psi_p5_M,-1)
       SMB2M10(i,j) = SC5cccrc(barpsi_p2_M,mup4,mup3,v1,p4,v2,psi_p5_M,-1)
       SMB2M11(i,j) = SC5cccrc(barpsi_p2_M,mup4,mup3,v1,p6,v2,psi_p5_M,-1)
       SMB2M12(i,j) = SC3crc(barpsi_p1_M,v1,p4,v2,psi_p6_M,-1)
       SMB2M13(i,j) = SC3crc(barpsi_p1_M,v1,p5,v2,psi_p6_M,-1)
       SMB2M14(i,j) = SC5crcrc(barpsi_p2_M,mup4,p6,v1,p1,v2,psi_p5_M,-1)
       SMB2M15(i,j) = SC5crcrc(barpsi_p2_M,mup3,p6,v1,p1,v2,psi_p5_M,-1)
       SMB2M16(i,j) = SC5crcrc(barpsi_p2_M,mup3,p6,v1,p4,v2,psi_p5_M,-1)
       SMB2M17(i,j) = SC5rrcrc(barpsi_p2_M,p6,p4,v1,p1,v2,psi_p5_M,-1)
       SMB2M18(i,j) = SC7ccrrcrc(barpsi_p2_M,mup4,mup3,p6,p4,v1,p1,v2,psi_p5_M,-1)
       SMB2M19(i,j) = SC5crcrc(barpsi_p2_M,mup4,p4,v1,p1,v2,psi_p5_M,-1)
       SMB2M20(i,j) = SC5crcrc(barpsi_p2_M,mup3,p4,v1,p1,v2,psi_p5_M,-1)
       SMB2M21(i,j) = SC5crcrc(barpsi_p2_M,mup4,p6,v1,p4,v2,psi_p5_M,-1)
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
       SMB3P1(i,j,k) = SC3ccc(barpsi_p1_P,v1,v2,v3,psi_p6_P,+1)
       SMB3P2(i,j,k) = SC3ccc(barpsi_p2_P,v1,v2,v3,psi_p5_P,+1)
       SMB3P3(i,j,k) = SC5ccccc(barpsi_p2_P,mup4,mup3,v1,v2,v3,psi_p5_P,+1)
       SMB3P4(i,j,k) = SC5crccc(barpsi_p2_P,mup3,p6,v1,v2,v3,psi_p5_P,+1)
       SMB3P5(i,j,k) = SC5rrccc(barpsi_p2_P,p6,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P6(i,j,k) = SC7ccrrccc(barpsi_p2_P,mup4,mup3,p6,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P7(i,j,k) = SC5crccc(barpsi_p2_P,mup3,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P8(i,j,k) = SC5crccc(barpsi_p2_P,mup3,p4,v1,v2,v3,psi_p5_P,+1)
       SMB3P9(i,j,k) = SC5rrccc(barpsi_p2_P,p4,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P10(i,j,k) = SC7ccrrccc(barpsi_p2_P,mup4,mup3,p4,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3P11(i,j,k) = SC5crccc(barpsi_p2_P,mup4,p6,v1,v2,v3,psi_p5_P,+1)
       SMB3P12(i,j,k) = SC5rrccc(barpsi_p2_P,p6,p4,v1,v2,v3,psi_p5_P,+1)
       SMB3P13(i,j,k) = SC7ccrrccc(barpsi_p2_P,mup4,mup3,p6,p4,v1,v2,v3,psi_p5_P,+1)
       SMB3P14(i,j,k) = SC5crccc(barpsi_p2_P,mup4,p4,v1,v2,v3,psi_p5_P,+1)
       SMB3P15(i,j,k) = SC5crccc(barpsi_p2_P,mup4,p1,v1,v2,v3,psi_p5_P,+1)
       SMB3M1(i,j,k) = SC3ccc(barpsi_p1_M,v1,v2,v3,psi_p6_M,-1)
       SMB3M2(i,j,k) = SC3ccc(barpsi_p2_M,v1,v2,v3,psi_p5_M,-1)
       SMB3M3(i,j,k) = SC5ccccc(barpsi_p2_M,mup4,mup3,v1,v2,v3,psi_p5_M,-1)
       SMB3M4(i,j,k) = SC5crccc(barpsi_p2_M,mup3,p6,v1,v2,v3,psi_p5_M,-1)
       SMB3M5(i,j,k) = SC5rrccc(barpsi_p2_M,p6,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M6(i,j,k) = SC7ccrrccc(barpsi_p2_M,mup4,mup3,p6,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M7(i,j,k) = SC5crccc(barpsi_p2_M,mup3,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M8(i,j,k) = SC5crccc(barpsi_p2_M,mup3,p4,v1,v2,v3,psi_p5_M,-1)
       SMB3M9(i,j,k) = SC5rrccc(barpsi_p2_M,p4,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M10(i,j,k) = SC7ccrrccc(barpsi_p2_M,mup4,mup3,p4,p1,v1,v2,v3,psi_p5_M,-1)
       SMB3M11(i,j,k) = SC5crccc(barpsi_p2_M,mup4,p6,v1,v2,v3,psi_p5_M,-1)
       SMB3M12(i,j,k) = SC5rrccc(barpsi_p2_M,p6,p4,v1,v2,v3,psi_p5_M,-1)
       SMB3M13(i,j,k) = SC7ccrrccc(barpsi_p2_M,mup4,mup3,p6,p4,v1,v2,v3,psi_p5_M,-1)
       SMB3M14(i,j,k) = SC5crccc(barpsi_p2_M,mup4,p4,v1,v2,v3,psi_p5_M,-1)
       SMB3M15(i,j,k) = SC5crccc(barpsi_p2_M,mup4,p1,v1,v2,v3,psi_p5_M,-1)
      enddo
       enddo
      enddo
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c   heP-1 heM-1 hMP-1 hMM-1
c************************************************************************************
       SMB(1) = SMB0M(1)*SMB0M(2)
       SMB(2) = SMB0M(3)*SMB0M(4)
       SMB(5) = SMB0M(5)*SMB0M(4)
       SMB(6) = SMB0M(5)*SMB0M(2)
       SMB(7) = SMB0M(5)*SMB0M(6)
       SMB(8) = SMB0M(5)*SMB0M(7)
       SMB(9) = SMB0M(5)*SMB0M(8)
       SMB(10) = SMB0M(5)*SMB0M(9)
       SMB(11) = SMB0M(10)*SMB0M(4)
       SMB(12) = SMB0M(11)*SMB0M(4)
       SMB(13) = SMB0M(12)*SMB0M(2)
       SMB(14) = SMB0M(13)*SMB0M(2)
       SMB(19) = SMB0M(14)*SMB0M(4)
       SMB(20) = SMB0M(14)*SMB0M(2)
       SMB(21) = SMB0M(14)*SMB0M(15)
       SMB(22) = SMB0M(14)*SMB0M(16)
       SMB(23) = SMB0M(14)*SMB0M(6)
       SMB(24) = SMB0M(14)*SMB0M(7)
       SMB(25) = SMB0M(14)*SMB0M(8)
       SMB(26) = SMB0M(14)*SMB0M(9)
       SMB(27) = SMB0M(17)*SMB0M(4)
       SMB(28) = SMB0M(17)*SMB0M(2)
       SMB(29) = SMB0M(17)*SMB0M(6)
       SMB(30) = SMB0M(17)*SMB0M(7)
       SMB(31) = SMB0M(17)*SMB0M(15)
       SMB(32) = SMB0M(17)*SMB0M(8)
       SMB(33) = SMB0M(17)*SMB0M(9)
       SMB(34) = SMB0M(17)*SMB0M(16)
       SMB(35) = SMB0M(5)*SMB0M(15)
       SMB(36) = SMB0M(5)*SMB0M(16)
       SMB(37) = SMB0M(18)*SMB0M(4)
       SMB(38) = SMB0M(18)*SMB0M(2)
       SMB(39) = SMB0M(18)*SMB0M(6)
       SMB(40) = SMB0M(18)*SMB0M(7)
       SMB(41) = SMB0M(18)*SMB0M(15)
       SMB(42) = SMB0M(18)*SMB0M(8)
       SMB(43) = SMB0M(18)*SMB0M(9)
       SMB(44) = SMB0M(18)*SMB0M(16)
       SMB(45) = SMB0M(1)*SMB0M(15)
       SMB(46) = SMB0M(1)*SMB0M(19)
       SMB(47) = SMB0M(1)*SMB0M(20)
       SMB(48) = SMB0M(3)*SMB0M(15)
       SMB(49) = SMB0M(3)*SMB0M(21)
       SMB(50) = SMB0M(3)*SMB0M(22)
       SMB(118) = SMB0M(1)*SMB0M(7)
       SMB(119) = SMB0M(3)*SMB0M(7)
       SMB(120) = SMB0M(3)*SMB0M(23)
       SMB(121) = SMB0M(24)*SMB0M(4)
       SMB(122) = SMB0M(25)*SMB0M(2)
       SMB(123) = SMB0M(26)*SMB0M(7)
       SMB(124) = SMB0M(27)*SMB0M(7)
       SMB(125) = SMB0M(28)*SMB0M(7)
       SMB(126) = SMB0M(1)*SMB0M(6)
       SMB(127) = SMB0M(1)*SMB0M(29)
       SMB(128) = SMB0M(3)*SMB0M(6)
       SMB(129) = SMB0M(5)*SMB0M(23)
       SMB(130) = SMB0M(5)*SMB0M(21)
       SMB(131) = SMB0M(5)*SMB0M(22)
       SMB(132) = SMB0M(5)*SMB0M(29)
       SMB(133) = SMB0M(5)*SMB0M(30)
       SMB(134) = SMB0M(5)*SMB0M(31)
       SMB(135) = SMB0M(12)*SMB0M(6)
       SMB(136) = SMB0M(13)*SMB0M(6)
       SMB(137) = SMB0M(28)*SMB0M(6)
       SMB(138) = SMB0M(11)*SMB0M(6)
       SMB(139) = SMB0M(12)*SMB0M(7)
       SMB(140) = SMB0M(13)*SMB0M(7)
       SMB(141) = SMB0M(12)*SMB0M(15)
       SMB(142) = SMB0M(13)*SMB0M(15)
       SMB(143) = SMB0M(12)*SMB0M(29)
       SMB(144) = SMB0M(13)*SMB0M(29)
       SMB(145) = SMB0M(13)*SMB0M(19)
       SMB(146) = SMB0M(10)*SMB0M(23)
       SMB(147) = SMB0M(11)*SMB0M(23)
       SMB(148) = SMB0M(11)*SMB0M(21)
       SMB(162) = SMB0M(5)*SMB0M(20)
       SMB(163) = SMB0M(10)*SMB0M(7)
       SMB(164) = SMB0M(11)*SMB0M(7)
       SMB(165) = SMB0M(28)*SMB0M(15)
       SMB(166) = SMB0M(11)*SMB0M(15)
       SMB(167) = SMB0M(12)*SMB0M(20)
       SMB(168) = SMB0M(13)*SMB0M(20)
       SMB(169) = SMB0M(10)*SMB0M(22)
       SMB(170) = SMB0M(11)*SMB0M(22)
       SMB(178) = SMB0M(26)*SMB0M(6)
       SMB(179) = SMB0M(26)*SMB0M(15)
       SMB(186) = SMB0M(14)*SMB0M(23)
       SMB(187) = SMB0M(14)*SMB0M(21)
       SMB(188) = SMB0M(14)*SMB0M(22)
       SMB(189) = SMB0M(14)*SMB0M(20)
       SMB(190) = SMB0M(14)*SMB0M(30)
       SMB(191) = SMB0M(14)*SMB0M(31)
       SMB(192) = SMB0M(17)*SMB0M(23)
       SMB(193) = SMB0M(17)*SMB0M(22)
       SMB(194) = SMB0M(24)*SMB0M(7)
       SMB(195) = SMB0M(24)*SMB0M(23)
       SMB(196) = SMB0M(24)*SMB0M(22)
       SMB(203) = SMB0M(17)*SMB0M(21)
       SMB(204) = SMB0M(17)*SMB0M(19)
       SMB(205) = SMB0M(17)*SMB0M(20)
       SMB(206) = SMB0M(17)*SMB0M(30)
       SMB(207) = SMB0M(17)*SMB0M(31)
       SMB(208) = SMB0M(14)*SMB0M(29)
       SMB(239) = SMB0M(1)*SMB0M(30)
       SMB(240) = SMB0M(3)*SMB0M(30)
       SMB(241) = SMB0M(17)*SMB0M(29)
       SMB(242) = SMB0M(25)*SMB0M(6)
       SMB(243) = SMB0M(27)*SMB0M(6)
       SMB(244) = SMB0M(24)*SMB0M(6)
       SMB(245) = SMB0M(25)*SMB0M(7)
       SMB(246) = SMB0M(25)*SMB0M(15)
       SMB(247) = SMB0M(27)*SMB0M(15)
       SMB(253) = SMB0M(25)*SMB0M(29)
       SMB(254) = SMB0M(25)*SMB0M(19)
       SMB(257) = SMB0M(26)*SMB0M(30)
       SMB(258) = SMB0M(27)*SMB0M(30)
       SMB(259) = SMB0M(28)*SMB0M(30)
       SMB(263) = SMB0M(24)*SMB0M(21)
       SMB(287) = SMB0M(24)*SMB0M(15)
       SMB(288) = SMB0M(25)*SMB0M(20)
c************************************************************************************
       SMB(3) = dotcc(SMB1M1,SMB1M2)
       SMB(4) = dotcc(SMB1M1,SMB1M3)
       SMB(15) = dotcc(SMB1M4,SMB1M2)
       SMB(16) = dotcc(SMB1M5,SMB1M2)
       SMB(17) = dotcc(SMB1M4,SMB1M3)
       SMB(18) = dotcc(SMB1M5,SMB1M3)
       SMB(51) = dotcc(SMB1M1,SMB1M6)
       SMB(52) = dotcc(SMB1M1,SMB1M7)
       SMB(53) = dotcc(SMB1M1,SMB1M8)
       SMB(54) = dotcc(SMB1M1,SMB1M9)
       SMB(55) = dotcc(SMB1M1,SMB1M10)
       SMB(56) = dotcc(SMB1M1,SMB1M11)
       SMB(57) = dotcc(SMB1M1,SMB1M12)
       SMB(58) = dotcc(SMB1M1,SMB1M13)
       SMB(59) = dotcc(SMB1M1,SMB1M14)
       SMB(60) = dotcc(SMB1M1,SMB1M15)
       SMB(61) = dotcc(SMB1M16,SMB1M2)
       SMB(62) = dotcc(SMB1M17,SMB1M2)
       SMB(63) = dotcc(SMB1M18,SMB1M2)
       SMB(64) = dotcc(SMB1M19,SMB1M2)
       SMB(65) = dotcc(SMB1M20,SMB1M2)
       SMB(66) = dotcc(SMB1M21,SMB1M2)
       SMB(67) = dotcc(SMB1M22,SMB1M2)
       SMB(69) = dotcc(SMB1M16,SMB1M9)
       SMB(70) = dotcc(SMB1M16,SMB1M10)
       SMB(71) = dotcc(SMB1M16,SMB1M11)
       SMB(72) = dotcc(SMB1M17,SMB1M9)
       SMB(73) = dotcc(SMB1M17,SMB1M10)
       SMB(74) = dotcc(SMB1M17,SMB1M11)
       SMB(75) = dotcc(SMB1M18,SMB1M9)
       SMB(76) = dotcc(SMB1M18,SMB1M10)
       SMB(77) = dotcc(SMB1M18,SMB1M11)
       SMB(79) = dotcc(SMB1M19,SMB1M6)
       SMB(80) = dotcc(SMB1M19,SMB1M7)
       SMB(81) = dotcc(SMB1M19,SMB1M8)
       SMB(82) = dotcc(SMB1M20,SMB1M6)
       SMB(83) = dotcc(SMB1M20,SMB1M7)
       SMB(84) = dotcc(SMB1M20,SMB1M8)
       SMB(85) = dotcc(SMB1M21,SMB1M6)
       SMB(86) = dotcc(SMB1M21,SMB1M7)
       SMB(87) = dotcc(SMB1M21,SMB1M8)
       SMB(93) = dotcc(SMB1M23,SMB1M2)
       SMB(97) = dotcc(SMB1M23,SMB1M3)
       SMB(102) = dotcc(SMB1M1,SMB1M24)
       SMB(106) = dotcc(SMB1M1,SMB1M25)
       SMB(149) = dotcc(SMB1M4,SMB1M6)
       SMB(150) = dotcc(SMB1M4,SMB1M7)
       SMB(151) = dotcc(SMB1M4,SMB1M8)
       SMB(152) = dotcc(SMB1M4,SMB1M24)
       SMB(153) = dotcc(SMB1M4,SMB1M25)
       SMB(154) = dotcc(SMB1M5,SMB1M6)
       SMB(155) = dotcc(SMB1M5,SMB1M7)
       SMB(156) = dotcc(SMB1M5,SMB1M8)
       SMB(157) = dotcc(SMB1M5,SMB1M9)
       SMB(158) = dotcc(SMB1M5,SMB1M24)
       SMB(159) = dotcc(SMB1M5,SMB1M12)
       SMB(160) = dotcc(SMB1M5,SMB1M25)
       SMB(161) = dotcc(SMB1M5,SMB1M14)
       SMB(171) = dotcc(SMB1M4,SMB1M10)
       SMB(172) = dotcc(SMB1M4,SMB1M13)
       SMB(173) = dotcc(SMB1M4,SMB1M15)
       SMB(174) = dotcc(SMB1M5,SMB1M10)
       SMB(175) = dotcc(SMB1M5,SMB1M11)
       SMB(176) = dotcc(SMB1M5,SMB1M13)
       SMB(177) = dotcc(SMB1M5,SMB1M15)
       SMB(180) = dotcc(SMB1M23,SMB1M8)
       SMB(181) = dotcc(SMB1M23,SMB1M11)
       SMB(182) = dotcc(SMB1M23,SMB1M12)
       SMB(183) = dotcc(SMB1M23,SMB1M13)
       SMB(184) = dotcc(SMB1M23,SMB1M14)
       SMB(185) = dotcc(SMB1M23,SMB1M15)
       SMB(197) = dotcc(SMB1M23,SMB1M6)
       SMB(198) = dotcc(SMB1M23,SMB1M7)
       SMB(199) = dotcc(SMB1M23,SMB1M10)
       SMB(200) = dotcc(SMB1M23,SMB1M24)
       SMB(201) = dotcc(SMB1M23,SMB1M25)
       SMB(202) = dotcc(SMB1M23,SMB1M9)
       SMB(209) = dotcc(SMB1M1,SMB1M26)
       SMB(210) = dotcc(SMB1M1,SMB1M27)
       SMB(212) = dotcc(SMB1M16,SMB1M12)
       SMB(213) = dotcc(SMB1M16,SMB1M13)
       SMB(214) = dotcc(SMB1M16,SMB1M27)
       SMB(215) = dotcc(SMB1M17,SMB1M12)
       SMB(216) = dotcc(SMB1M17,SMB1M13)
       SMB(217) = dotcc(SMB1M17,SMB1M27)
       SMB(218) = dotcc(SMB1M18,SMB1M12)
       SMB(219) = dotcc(SMB1M18,SMB1M13)
       SMB(220) = dotcc(SMB1M18,SMB1M27)
       SMB(221) = dotcc(SMB1M22,SMB1M12)
       SMB(223) = dotcc(SMB1M19,SMB1M26)
       SMB(224) = dotcc(SMB1M20,SMB1M26)
       SMB(225) = dotcc(SMB1M21,SMB1M26)
       SMB(251) = dotcc(SMB1M16,SMB1M24)
       SMB(252) = dotcc(SMB1M17,SMB1M24)
       SMB(255) = dotcc(SMB1M18,SMB1M24)
       SMB(256) = dotcc(SMB1M22,SMB1M24)
       SMB(261) = dotcc(SMB1M19,SMB1M24)
       SMB(262) = dotcc(SMB1M20,SMB1M24)
       SMB(264) = dotcc(SMB1M21,SMB1M24)
       SMB(276) = dotcc(SMB1M22,SMB1M13)
       SMB(279) = dotcc(SMB1M19,SMB1M13)
       SMB(280) = dotcc(SMB1M20,SMB1M13)
       SMB(281) = dotcc(SMB1M21,SMB1M13)
c************************************************************************************
       SMB(68) = ten2cc(SMB2M1,SMB2M2)
       SMB(78) = ten2cc(SMB2M3,SMB2M4)
       SMB(90) = ten2cc(SMB2M5,SMB2M4)
       SMB(91) = ten2cc(SMB2M5,SMB2M2)
       SMB(92) = ten2cc(SMB2M5,SMB2M6)
       SMB(94) = ten2cc(SMB2M5,SMB2M7)
       SMB(95) = ten2cc(SMB2M5,SMB2M8)
       SMB(96) = ten2cc(SMB2M5,SMB2M9)
       SMB(98) = ten2cc(SMB2M5,SMB2M10)
       SMB(99) = ten2cc(SMB2M5,SMB2M11)
       SMB(100) = ten2cc(SMB2M12,SMB2M4)
       SMB(101) = ten2cc(SMB2M12,SMB2M2)
       SMB(103) = ten2cc(SMB2M12,SMB2M6)
       SMB(104) = ten2cc(SMB2M12,SMB2M7)
       SMB(105) = ten2cc(SMB2M12,SMB2M8)
       SMB(107) = ten2cc(SMB2M12,SMB2M9)
       SMB(108) = ten2cc(SMB2M12,SMB2M10)
       SMB(109) = ten2cc(SMB2M12,SMB2M11)
       SMB(110) = ten2cc(SMB2M13,SMB2M4)
       SMB(111) = ten2cc(SMB2M13,SMB2M2)
       SMB(112) = ten2cc(SMB2M13,SMB2M6)
       SMB(113) = ten2cc(SMB2M13,SMB2M7)
       SMB(114) = ten2cc(SMB2M13,SMB2M8)
       SMB(115) = ten2cc(SMB2M13,SMB2M9)
       SMB(116) = ten2cc(SMB2M13,SMB2M10)
       SMB(117) = ten2cc(SMB2M13,SMB2M11)
       SMB(211) = ten2cc(SMB2M1,SMB2M8)
       SMB(222) = ten2cc(SMB2M3,SMB2M15)
       SMB(227) = ten2cc(SMB2M5,SMB2M15)
       SMB(228) = ten2cc(SMB2M5,SMB2M16)
       SMB(229) = ten2cc(SMB2M5,SMB2M17)
       SMB(230) = ten2cc(SMB2M5,SMB2M18)
       SMB(231) = ten2cc(SMB2M12,SMB2M15)
       SMB(232) = ten2cc(SMB2M12,SMB2M16)
       SMB(233) = ten2cc(SMB2M12,SMB2M17)
       SMB(234) = ten2cc(SMB2M12,SMB2M18)
       SMB(235) = ten2cc(SMB2M13,SMB2M15)
       SMB(236) = ten2cc(SMB2M13,SMB2M16)
       SMB(237) = ten2cc(SMB2M13,SMB2M17)
       SMB(238) = ten2cc(SMB2M13,SMB2M18)
       SMB(248) = ten2cc(SMB2M1,SMB2M6)
       SMB(249) = ten2cc(SMB2M1,SMB2M7)
       SMB(250) = ten2cc(SMB2M1,SMB2M19)
       SMB(260) = ten2cc(SMB2M3,SMB2M20)
       SMB(269) = ten2cc(SMB2M5,SMB2M20)
       SMB(270) = ten2cc(SMB2M5,SMB2M19)
       SMB(271) = ten2cc(SMB2M12,SMB2M20)
       SMB(272) = ten2cc(SMB2M12,SMB2M19)
       SMB(273) = ten2cc(SMB2M13,SMB2M20)
       SMB(274) = ten2cc(SMB2M13,SMB2M19)
       SMB(275) = ten2cc(SMB2M1,SMB2M21)
       SMB(277) = ten2cc(SMB2M3,SMB2M8)
       SMB(278) = ten2cc(SMB2M3,SMB2M16)
       SMB(284) = ten2cc(SMB2M5,SMB2M21)
       SMB(285) = ten2cc(SMB2M12,SMB2M21)
       SMB(286) = ten2cc(SMB2M13,SMB2M21)
       SMB(289) = ten2cc(SMB2M3,SMB2M7)
       SMB(291) = ten2cc(SMB2M3,SMB2M6)
c************************************************************************************
       SMB(88) = ten3cc(SMB3M1,SMB3M2)
       SMB(89) = ten3cc(SMB3M1,SMB3M3)
       SMB(226) = ten3cc(SMB3M1,SMB3M4)
       SMB(265) = ten3cc(SMB3M1,SMB3M7)
       SMB(266) = ten3cc(SMB3M1,SMB3M8)
       SMB(267) = ten3cc(SMB3M1,SMB3M9)
       SMB(268) = ten3cc(SMB3M1,SMB3M10)
       SMB(282) = ten3cc(SMB3M1,SMB3M12)
       SMB(283) = ten3cc(SMB3M1,SMB3M13)
       SMB(290) = ten3cc(SMB3M1,SMB3M14)
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
