       subroutine PrintIntHexMassLess4g_QUAD(p1,p2,p3,p4,p5,p6,musq)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 27/10/2008
c Modified:17/11/2008
c ************************************************************************************
c determine the Hexagon of
c p1(gluon)    ---->------->--------->--------->---   p6(B) p6^2=M_h
c                $         $          $         $     
c                $         $          $         $     
c                $         $          $         $     
c                V         V          V         V     
c                $         $          $         $     
c                $         $          $         $     
c               p2         p3         p4        p5     
c              mu_p2      mu_p3      mu_p4     mu_p5      
c Note: No couplings are assigned
c Note:m the mass of the fermion in the loop  
c Note: To make it shorter in the promgram: mu_p2,...->mup2,... 
c Notation of External momenta: p1+p2+p3+p4+p5+p6=0 
c mu_p1,mu_p2,mu_p3,mu_p4,mu_p5 should be think as external current 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Real*8 psq(6)
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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c   QUAD Variables
c
c
c
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       Real*16  p1_QUAD(0:3),p2_QUAD(0:3),p3_QUAD(0:3),p4_QUAD(0:3),
     -          p5_QUAD(0:3),p6_QUAD(0:3)
       Real*16 ptemp(0:3,6)

       Real*16   p1sq_QUAD, p1p2_QUAD, p1p3_QUAD, p1p4_QUAD, p1p5_QUAD,p1p6_QUAD
       Real*16   p2sq_QUAD, p2p3_QUAD, p2p4_QUAD, p2p5_QUAD 
       Real*16   p3sq_QUAD, p3p4_QUAD, p3p5_QUAD 
       Real*16   p4sq_QUAD, p4p5_QUAD 
       Real*16   p5sq_QUAD, p5p6_QUAD 
       Real*16   p6sq_QUAD 
       Real*16   s12_QUAD, s13_QUAD, s14_QUAD, s15_QUAD, s16_QUAD 
       Real*16   s23_QUAD, s24_QUAD, s25_QUAD,s26_QUAD 
       Real*16   s34_QUAD, s35_QUAD,s36_QUAD 
       Real*16   s45_QUAD,s46_QUAD 
       Real*16   s56_QUAD 
       Real*16   s123_QUAD,s234_QUAD,s345_QUAD  
       Real*16   p2p6_QUAD,p3p6_QUAD,p4p6_QUAD  
       Complex*32  B012_QUAD,B013_QUAD,B014_QUAD,B015_QUAD,B016_QUAD 
       Complex*32  B023_QUAD,B024_QUAD,B025_QUAD,B026_QUAD 
       Complex*32  B034_QUAD,B035_QUAD,B036_QUAD 
       Complex*32  B045_QUAD,B046_QUAD 
       Complex*32  B056_QUAD   
       Real*16  B012R_QUAD,B013R_QUAD,B014R_QUAD,B015R_QUAD,B016R_QUAD 
       Real*16  B023R_QUAD,B024R_QUAD,B025R_QUAD,B026R_QUAD 
       Real*16  B034R_QUAD,B035R_QUAD,B036R_QUAD 
       Real*16  B045R_QUAD,B046R_QUAD 
       Real*16  B056R_QUAD  
       Real*16  B012I_QUAD,B013I_QUAD,B014I_QUAD,B015I_QUAD,B016I_QUAD 
       Real*16  B023I_QUAD,B024I_QUAD,B025I_QUAD,B026I_QUAD 
       Real*16  B034I_QUAD,B035I_QUAD,B036I_QUAD 
       Real*16  B045I_QUAD,B046I_QUAD 
       Real*16  B056I_QUAD   
       Complex*32 C0123_QUAD,C0124_QUAD,C0125_QUAD,C0126_QUAD 
       Complex*32 C0134_QUAD,C0135_QUAD,C0136_QUAD 
       Complex*32 C0145_QUAD,C0146_QUAD 
       Complex*32 C0156_QUAD 
       Complex*32 C0234_QUAD,C0235_QUAD,C0236_QUAD 
       Complex*32 C0245_QUAD,C0246_QUAD 
       Complex*32 C0256_QUAD 
       Complex*32 C0345_QUAD,C0346_QUAD 
       Complex*32 C0356_QUAD 
       Complex*32 C0456_QUAD   
       Real*16 C0123R_QUAD,C0124R_QUAD,C0125R_QUAD,C0126R_QUAD 
       Real*16 C0134R_QUAD,C0135R_QUAD,C0136R_QUAD 
       Real*16 C0145R_QUAD,C0146R_QUAD 
       Real*16 C0156R_QUAD 
       Real*16 C0234R_QUAD,C0235R_QUAD,C0236R_QUAD 
       Real*16 C0245R_QUAD,C0246R_QUAD 
       Real*16 C0256R_QUAD 
       Real*16 C0345R_QUAD,C0346R_QUAD 
       Real*16 C0356R_QUAD 
       Real*16 C0456R_QUAD    
       Real*16 C0123I_QUAD,C0124I_QUAD,C0125I_QUAD,C0126I_QUAD 
       Real*16 C0134I_QUAD,C0135I_QUAD,C0136I_QUAD 
       Real*16 C0145I_QUAD,C0146I_QUAD 
       Real*16 C0156I_QUAD 
       Real*16 C0234I_QUAD,C0235I_QUAD,C0236I_QUAD 
       Real*16 C0245I_QUAD,C0246I_QUAD 
       Real*16 C0256I_QUAD 
       Real*16 C0345I_QUAD,C0346I_QUAD 
       Real*16 C0356I_QUAD 
       Real*16 C0456I_QUAD      
       Real*16 C123R_QUAD(4,2),C124R_QUAD(4,2),C125R_QUAD(4,2),C126R_QUAD(4,2) 
       Real*16 C134R_QUAD(4,2),C135R_QUAD(4,2),C136R_QUAD(4,2) 
       Real*16 C145R_QUAD(4,2),C146R_QUAD(4,2) 
       Real*16 C156R_QUAD(4,2) 
       Real*16 C234R_QUAD(4,2),C235R_QUAD(4,2),C236R_QUAD(4,2) 
       Real*16 C245R_QUAD(4,2),C246R_QUAD(4,2) 
       Real*16 C256R_QUAD(4,2) 
       Real*16 C345R_QUAD(4,2),C346R_QUAD(4,2) 
       Real*16 C356R_QUAD(4,2) 
       Real*16 C456R_QUAD(4,2)  
       Real*16 C123I_QUAD(4,2),C124I_QUAD(4,2),C125I_QUAD(4,2),C126I_QUAD(4,2) 
       Real*16 C134I_QUAD(4,2),C135I_QUAD(4,2),C136I_QUAD(4,2) 
       Real*16 C145I_QUAD(4,2),C146I_QUAD(4,2) 
       Real*16 C156I_QUAD(4,2) 
       Real*16 C234I_QUAD(4,2),C235I_QUAD(4,2),C236I_QUAD(4,2) 
       Real*16 C245I_QUAD(4,2),C246I_QUAD(4,2) 
       Real*16 C256I_QUAD(4,2) 
       Real*16 C345I_QUAD(4,2),C346I_QUAD(4,2) 
       Real*16 C356I_QUAD(4,2) 
       Real*16 C456I_QUAD(4,2)  
       Complex*32  D01234_QUAD,D01235_QUAD,D01236_QUAD 
       Complex*32 D01245_QUAD,D01246_QUAD 
       Complex*32 D01256_QUAD 
       Complex*32 D01345_QUAD,D01346_QUAD 
       Complex*32 D01356_QUAD 
       Complex*32 D01456_QUAD 
       Complex*32 D02345_QUAD,D02346_QUAD 
       Complex*32 D02356_QUAD 
       Complex*32 D02456_QUAD 
       Complex*32 D03456_QUAD    
       Real*16 D01234R_QUAD,D01235R_QUAD,D01236R_QUAD 
       Real*16 D01245R_QUAD,D01246R_QUAD 
       Real*16 D01256R_QUAD 
       Real*16 D01345R_QUAD,D01346R_QUAD 
       Real*16 D01356R_QUAD 
       Real*16 D01456R_QUAD 
       Real*16 D02345R_QUAD,D02346R_QUAD 
       Real*16 D02356R_QUAD 
       Real*16 D02456R_QUAD 
       Real*16 D03456R_QUAD  
       Real*16 D01234I_QUAD,D01235I_QUAD,D01236I_QUAD 
       Real*16 D01245I_QUAD,D01246I_QUAD 
       Real*16 D01256I_QUAD 
       Real*16 D01345I_QUAD,D01346I_QUAD 
       Real*16 D01356I_QUAD 
       Real*16 D01456I_QUAD 
       Real*16 D02345I_QUAD,D02346I_QUAD 
       Real*16 D02356I_QUAD 
       Real*16 D02456I_QUAD 
       Real*16 D03456I_QUAD   
       Real*16 D1234R_QUAD(13,3),D1235R_QUAD(13,3),D1236R_QUAD(13,3) 
       Real*16 D1245R_QUAD(13,3),D1246R_QUAD(13,3) 
       Real*16 D1256R_QUAD(13,3) 
       Real*16 D1345R_QUAD(13,3),D1346R_QUAD(13,3) 
       Real*16 D1356R_QUAD(13,3) 
       Real*16 D1456R_QUAD(13,3) 
       Real*16 D2345R_QUAD(13,3),D2346R_QUAD(13,3) 
       Real*16 D2356R_QUAD(13,3) 
       Real*16 D2456R_QUAD(13,3) 
       Real*16 D3456R_QUAD(13,3)   
       Real*16 D1234I_QUAD(13,3),D1235I_QUAD(13,3),D1236I_QUAD(13,3) 
       Real*16 D1245I_QUAD(13,3),D1246I_QUAD(13,3) 
       Real*16 D1256I_QUAD(13,3) 
       Real*16 D1345I_QUAD(13,3),D1346I_QUAD(13,3) 
       Real*16 D1356I_QUAD(13,3) 
       Real*16 D1456I_QUAD(13,3) 
       Real*16 D2345I_QUAD(13,3),D2346I_QUAD(13,3) 
       Real*16 D2356I_QUAD(13,3) 
       Real*16 D2456I_QUAD(13,3) 
       Real*16 D3456I_QUAD(13,3)   
       Complex*32 E012345_QUAD,E012346_QUAD,E012356_QUAD,E012456_QUAD,E013456_QUAD,E023456_QUAD 
       Real*16 E012345R_QUAD,E012346R_QUAD,E012356R_QUAD,E012456R_QUAD,E013456R_QUAD,E023456R_QUAD 
       Real*16 E012345I_QUAD,E012346I_QUAD,E012356I_QUAD,E012456I_QUAD,E013456I_QUAD,E023456I_QUAD  
       Real*16 E12345R_QUAD(46,4),E12346R_QUAD(46,4),E12356R_QUAD(46,4) 
       Real*16 E12456R_QUAD(46,4),E13456R_QUAD(46,4),E23456R_QUAD(46,4)  
       Real*16 E12345I_QUAD(46,4),E12346I_QUAD(46,4),E12356I_QUAD(46,4) 
       Real*16 E12456I_QUAD(46,4),E13456I_QUAD(46,4),E23456I_QUAD(46,4)  
       Complex*32 F0123456_QUAD 
       Real*16 F0123456R_QUAD 
       Real*16 F0123456I_QUAD 
       Real*16 F123456R_QUAD(166,5) 
       Real*16 F123456I_QUAD(166,5) 
       Complex*32 F123456_QUAD(166,5) 
       Real*16 musq_QUAD, dotrr_QUAD 
       External dotrr_QUAD 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      End QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       COMMON/NoAbe4vertHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45
     -   ,s56,s16,s123,s234,s345
       COMMON/NoAbe4vertHexLineFVALFunctions/F0123456R,F123456R,F0123456I,F123456I
       COMMON/NoAbe4vertHexLineEVALFunctions/ E012345R,E12345R,E012345I,E12345I, E01234
     -   6R,E12346R,E012346I,E12346I, E012356R,E12356R,E012356I,E1235
     -   6I, E012456R,E12456R,E012456I,E12456I, E013456R,E13456R,E013
     -   456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/NoAbe4vertHexLineDVALFunctions/ D01234R,D1234R,D01234I,D1234I, D01235R,D1
     -   235R,D01235I,D1235I, D01236R,D1236R,D01236I,D1236I, D01245R,
     -   D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1246I, D01256
     -   R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,D1345I, D013
     -   46R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356I,D1356I, D0
     -   1456R,D1456R,D01456I,D1456I, D02345R,D2345R,D02345I,D2345I, 
     -   D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D02356I,D2356I
     -   , D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,D03456I,D345
     -   6I
       COMMON/NoAbe4vertHexLineCVALFunctions/ C0123R,C123R,C0123I,C123I, C0124R,C124R,C
     -   0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C126R,C0126I,
     -   C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C0135I,C135I,
     -    C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,C145I, C0146
     -   R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I, C0234R,C234
     -   R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236R,C236R,C023
     -   6I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246R,C0246I,C24
     -   6I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C0345I,C345I, C0
     -   346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C356I, C0456R,C
     -   456R,C0456I,C456I
       COMMON/NoAbe4vertHexLineBVALFunctions/ B012R,B012I, B013R,B013I, B014R,B014I, B0
     -   15R,B015I, B016R,B016I, B023R,B023I, B024R,B024I, B025R,B025
     -   I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,B036I, B045R
     -   ,B045I, B046R,B046I, B056R,B056I
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       integer i,jtemp,j
       logical PrintB,PrintC,PrintD,PrintE,PrintF
       real*8 dotrr
       external dotrr
       real*8 musq
       do i=0,3 
       ptemp(i,1)=p1(i)
       ptemp(i,2)=p2(i)
       ptemp(i,3)=p3(i)
       ptemp(i,4)=p4(i)
       ptemp(i,5)=p5(i)
       ptemp(i,6)=p6(i)
       enddo
       psq(1)=abs(dotrr(p1,p1))
       psq(2)=abs(dotrr(p2,p2))
       psq(3)=abs(dotrr(p3,p3))
       psq(4)=abs(dotrr(p4,p4))
       psq(5)=abs(dotrr(p5,p5))
       psq(6)=abs(dotrr(p6,p6))
       if (psq(1).gt.10d-7) then
       j=1
       elseif (psq(2).gt.10d-7)then
       j=2
       elseif (psq(3).gt.10d-7)then
       j=3
       elseif (psq(4).gt.10d-7)then
       j=4
       elseif (psq(5).gt.10d-7)then
       j=5
       elseif (psq(6).gt.10d-7)then
       j=6
       else
       j=7
c$$$       write(*,*) 'No massive invariant'
c$$$       Return
       endif
       if(j.ne.7) then
       do jtemp=1,6
       if (psq(jtemp).le.10d-7) then
       if(abs(ptemp(1,jtemp)).le.10d-7)then
       if(abs(ptemp(2,jtemp)).le.10d-7) then
       ptemp(0,jtemp)=abs(ptemp(3,jtemp))
       else
       if (ptemp(0,jtemp).ge.0d0) then
       ptemp(0,jtemp)=SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
       else
c$$$       print*, 'Minus jtemp', jtemp
      ptemp(0,jtemp)=-SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
      endif
      endif
      else
      if (ptemp(0,jtemp).ge.0d0) then
c$$$       print*, 'jtemp', jtemp
      ptemp(0,jtemp)=SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
      else
c$$$       print*, 'Minus jtemp', jtemp
      ptemp(0,jtemp)=-SQRT(ptemp(1,jtemp)*ptemp(1,jtemp)+
     -                  ptemp(2,jtemp)*ptemp(2,jtemp)+
     -                  ptemp(3,jtemp)*ptemp(3,jtemp))
      endif
      endif
      endif
      enddo
      do i=0,3
      ptemp(i,j)=0q0
      do jtemp=1,6
      if (jtemp.ne.j) then
      ptemp(i,j)=-ptemp(i,jtemp)+ptemp(i,j)
      endif
      enddo
      enddo
      elseif(j.eq.7) then
      do i=0,3 
      ptemp(i,1)=p1(i)
      ptemp(i,2)=p2(i)
      ptemp(i,3)=p3(i)
      ptemp(i,4)=p4(i)
      ptemp(i,5)=p5(i)
      ptemp(i,6)=p6(i)
      enddo 
      endif

      do i=0,3 
      p1_QUAD(i)=ptemp(i,1)
      p2_QUAD(i)=ptemp(i,2)
      p3_QUAD(i)=ptemp(i,3)
      p4_QUAD(i)=ptemp(i,4)
      p5_QUAD(i)=ptemp(i,5)
      p6_QUAD(i)=ptemp(i,6)
      enddo
c$$$      print*, 'p1QUAD',p1_QUAD
c$$$      print*, 'p2QUAD',p2_QUAD
c$$$      print*, 'p3QUAD',p3_QUAD
c$$$      print*, 'p4QUAD',p4_QUAD
c$$$      print*, 'p5QUAD',p5_QUAD
c$$$      print*, 'p6QUAD',p6_QUAD
c$$$      Print*, 'Conser',(p1_QUAD(i)+p2_QUAD(i)+p3_QUAD(i)+
c$$$     -   p4_QUAD(i)+p5_QUAD(i)+p6_QUAD(i),i=0,3)
c$$$      Print*, 'Conser',(p1(i)+p2(i)+p3(i)+
c$$$     -   p4(i)+p5(i)+p6(i),i=0,3)

       p1sq_QUAD = dotrr_QUAD(p1_QUAD,p1_QUAD)
       p1p2_QUAD = dotrr_QUAD(p1_QUAD,p2_QUAD)
       p1p3_QUAD = dotrr_QUAD(p1_QUAD,p3_QUAD)
       p1p4_QUAD = dotrr_QUAD(p1_QUAD,p4_QUAD)
       p1p5_QUAD = dotrr_QUAD(p1_QUAD,p5_QUAD)
       p1p6_QUAD = dotrr_QUAD(p1_QUAD,p6_QUAD)
       p2sq_QUAD = dotrr_QUAD(p2_QUAD,p2_QUAD)
       p2p3_QUAD = dotrr_QUAD(p2_QUAD,p3_QUAD)
       p2p4_QUAD = dotrr_QUAD(p2_QUAD,p4_QUAD)
       p2p5_QUAD = dotrr_QUAD(p2_QUAD,p5_QUAD)
       p2p6_QUAD = dotrr_QUAD(p2_QUAD,p6_QUAD)
       p3sq_QUAD = dotrr_QUAD(p3_QUAD,p3_QUAD)
       p3p4_QUAD = dotrr_QUAD(p3_QUAD,p4_QUAD)
       p3p5_QUAD = dotrr_QUAD(p3_QUAD,p5_QUAD)
       p3p6_QUAD = dotrr_QUAD(p3_QUAD,p6_QUAD)
       p4sq_QUAD = dotrr_QUAD(p4_QUAD,p4_QUAD)
       p4p5_QUAD = dotrr_QUAD(p4_QUAD,p5_QUAD)
       p4p6_QUAD = dotrr_QUAD(p4_QUAD,p6_QUAD)
       p5sq_QUAD = dotrr_QUAD(p5_QUAD,p5_QUAD)
       p5p6_QUAD = dotrr_QUAD(p5_QUAD,p6_QUAD)
       p6sq_QUAD = dotrr_QUAD(p6_QUAD,p6_QUAD)
       s12_QUAD = (p1sq_QUAD +p2sq_QUAD+ 2*p1p2_QUAD) 
       s13_QUAD = (p1sq_QUAD +p3sq_QUAD+ 2*p1p3_QUAD) 
       s14_QUAD = (p1sq_QUAD +p4sq_QUAD+ 2*p1p4_QUAD) 
       s15_QUAD = (p1sq_QUAD +p5sq_QUAD+ 2*p1p5_QUAD) 
       s16_QUAD = (p1sq_QUAD +p6sq_QUAD+ 2*p1p6_QUAD) 
       s23_QUAD = (p2sq_QUAD +p3sq_QUAD+ 2*p2p3_QUAD) 
       s24_QUAD = (p2sq_QUAD +p4sq_QUAD+ 2*p2p4_QUAD) 
       s25_QUAD = (p2sq_QUAD +p5sq_QUAD+ 2*p2p5_QUAD) 
       s26_QUAD = (p2sq_QUAD +p6sq_QUAD+ 2*p2p6_QUAD) 
       s34_QUAD = (p3sq_QUAD +p4sq_QUAD+ 2*p3p4_QUAD) 
       s35_QUAD = (p3sq_QUAD +p5sq_QUAD+ 2*p3p5_QUAD) 
       s36_QUAD = (p3sq_QUAD +p6sq_QUAD+ 2*p3p6_QUAD) 
       s45_QUAD = (p4sq_QUAD +p5sq_QUAD+ 2*p4p5_QUAD) 
       s46_QUAD = (p4sq_QUAD +p6sq_QUAD+ 2*p4p6_QUAD) 
       s56_QUAD = (p5sq_QUAD +p6sq_QUAD+ 2*p5p6_QUAD) 
c       Write(*,'(a5,E44.36)')," p1sq_QUAD ", p1sq_QUAD 
c       Write(*,'(a5,E44.36)')," p1p2_QUAD ", p1p2_QUAD
c       Write(*,'(a5,E44.36)')," p1p3_QUAD ", p1p3_QUAD
c       Write(*,'(a5,E44.36)')," p1p4_QUAD ", p1p4_QUAD
c       Write(*,'(a5,E44.36)')," p1p5_QUAD ", p1p5_QUAD
c       Write(*,'(a5,E44.36)')," p1p6_QUAD ", p1p6_QUAD
c       Write(*,'(a5,E44.36)')," p2sq_QUAD ", p2sq_QUAD 
c       Write(*,'(a5,E44.36)')," p2p3_QUAD ", p2p3_QUAD
c       Write(*,'(a5,E44.36)')," p2p4_QUAD ", p2p4_QUAD
c       Write(*,'(a5,E44.36)')," p2p5_QUAD ", p2p5_QUAD
c       Write(*,'(a5,E44.36)')," p2p6_QUAD ", p2p6_QUAD
c       Write(*,'(a5,E44.36)')," p3sq_QUAD ", p3sq_QUAD 
c       Write(*,'(a5,E44.36)')," p3p4_QUAD ", p3p4_QUAD
c       Write(*,'(a5,E44.36)')," p3p5_QUAD ", p3p5_QUAD
c       Write(*,'(a5,E44.36)')," p3p6_QUAD ", p3p6_QUAD
c       Write(*,'(a5,E44.36)')," p4sq_QUAD ", p4sq_QUAD 
c       Write(*,'(a5,E44.36)')," p4p5_QUAD ", p4p5_QUAD
c       Write(*,'(a5,E44.36)')," p4p6_QUAD ", p4p6_QUAD
c       Write(*,'(a5,E44.36)')," p5sq_QUAD ", p5sq_QUAD 
c       Write(*,'(a5,E44.36)')," p5p6_QUAD ", p5p6_QUAD
c       Write(*,'(a5,E44.36)')," p6sq_QUAD ", p6sq_QUAD 
      s123_QUAD=p1sq_QUAD+p2sq_QUAD+p3sq_QUAD
     -   +2*(p1p2_QUAD+p1p3_QUAD+p2p3_QUAD) 
      s234_QUAD=p2sq_QUAD+p3sq_QUAD+p4sq_QUAD
     -   +2*(p2p3_QUAD+p2p4_QUAD+p3p4_QUAD) 
      s345_QUAD=p3sq_QUAD+p4sq_QUAD+p5sq_QUAD
     -   +2*(p3p4_QUAD+p3p5_QUAD+p4p5_QUAD)
       p1sq=dble(p1sq_QUAD) 
       p1p2=dble(p1p2_QUAD)
       p1p3=dble(p1p3_QUAD)
       p1p4=dble(p1p4_QUAD)
       p1p5=dble(p1p5_QUAD)
       p1p6=dble(p1p6_QUAD)
       p2sq=dble(p2sq_QUAD) 
       p2p3=dble(p2p3_QUAD)
       p2p4=dble(p2p4_QUAD)
       p2p5=dble(p2p5_QUAD)
       p2p6=dble(p2p6_QUAD)
       p3sq=dble(p3sq_QUAD) 
       p3p4=dble(p3p4_QUAD)
       p3p5=dble(p3p5_QUAD)
       p3p6=dble(p3p6_QUAD)
       p4sq=dble(p4sq_QUAD) 
       p4p5=dble(p4p5_QUAD)
       p4p6=dble(p4p6_QUAD)
       p5sq=dble(p5sq_QUAD) 
       p5p6=dble(p5p6_QUAD)
       p6sq=dble(p6sq_QUAD) 
       s12=dble(s12_QUAD) 
       s13=dble(s13_QUAD) 
       s14=dble(s14_QUAD) 
       s15=dble(s15_QUAD) 
       s16=dble(s16_QUAD) 
       s23=dble(s23_QUAD) 
       s24=dble(s24_QUAD) 
       s25=dble(s25_QUAD) 
       s26=dble(s26_QUAD) 
       s34=dble(s34_QUAD) 
       s35=dble(s35_QUAD) 
       s36=dble(s36_QUAD) 
       s45=dble(s45_QUAD) 
       s46=dble(s46_QUAD) 
       s56=dble(s56_QUAD) 
      s123=dble(s123_QUAD) 
      s234=dble(s234_QUAD) 
      s345=dble(s345_QUAD)
       do i=0,3
         p1(i)=dble(p1_QUAD(i))
         p2(i)=dble(p2_QUAD(i))
         p3(i)=dble(p3_QUAD(i))
         p4(i)=dble(p4_QUAD(i))
         p5(i)=dble(p5_QUAD(i))
         p6(i)=dble(p6_QUAD(i))
         enddo

      If(abs(musq-s12).lt.10d-7) then
          musq_QUAD=s12_QUAD
          elseIf(abs(musq-s23).lt.10d-7) then
          musq_QUAD=s23_QUAD
          elseIf(abs(musq-s34).lt.10d-7) then
          musq_QUAD=s34_QUAD
          elseIf(abs(musq-s45).lt.10d-7) then
          musq_QUAD=s45_QUAD
          elseIf(abs(musq-s56).lt.10d-7) then
          musq_QUAD=s56_QUAD
          elseIf(abs(musq-s16).lt.10d-7) then
          musq_QUAD=s16_QUAD
          elseIf(abs(musq-s123).lt.10d-7) then
          musq_QUAD=s123_QUAD
          elseIf(abs(musq-s234).lt.10d-7) then
          musq_QUAD=s234_QUAD
          elseIf(abs(musq-s345).lt.10d-7) then
          musq_QUAD=s345_QUAD
          else
          musq_QUAD=musq
          endif
      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False. 
      PrintF=.False.
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************

c    Calling C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************
       call B0fin_QUAD(p1sq_QUAD,musq_QUAD,B012_QUAD)
       B012R_QUAD=real(B012_QUAD)
       B012I_QUAD=imag(B012_QUAD)
       call From_BijR_QUAD_to_Dbl(B012R,B012I,B012R_QUAD,B012I_QUAD)
       B012=DCMPLX(B012R,B012I)
       call B0fin_QUAD(p2sq_QUAD,musq_QUAD,B023_QUAD)
       B023R_QUAD=real(B023_QUAD)
       B023I_QUAD=imag(B023_QUAD)
       call From_BijR_QUAD_to_Dbl(B023R,B023I,B023R_QUAD,B023I_QUAD)
       B023=DCMPLX(B023R,B023I)
       call B0fin_QUAD(p3sq_QUAD,musq_QUAD,B034_QUAD)
       B034R_QUAD=real(B034_QUAD)
       B034I_QUAD=imag(B034_QUAD)
       call From_BijR_QUAD_to_Dbl(B034R,B034I,B034R_QUAD,B034I_QUAD)
       B034=DCMPLX(B034R,B034I)
       call B0fin_QUAD(p4sq_QUAD,musq_QUAD,B045_QUAD)
       B045R_QUAD=real(B045_QUAD)
       B045I_QUAD=imag(B045_QUAD)
       call From_BijR_QUAD_to_Dbl(B045R,B045I,B045R_QUAD,B045I_QUAD)
       B045=DCMPLX(B045R,B045I)
       call B0fin_QUAD(p5sq_QUAD,musq_QUAD,B056_QUAD)
       B056R_QUAD=real(B056_QUAD)
       B056I_QUAD=imag(B056_QUAD)
       call From_BijR_QUAD_to_Dbl(B056R,B056I,B056R_QUAD,B056I_QUAD)
       B056=DCMPLX(B056R,B056I)
       call B0fin_QUAD(p6sq_QUAD,musq_QUAD,B016_QUAD)
       B016R_QUAD=real(B016_QUAD)
       B016I_QUAD=imag(B016_QUAD)
       call From_BijR_QUAD_to_Dbl(B016R,B016I,B016R_QUAD,B016I_QUAD)
       B016=DCMPLX(B016R,B016I)
       call B0fin_QUAD(s12_QUAD,musq_QUAD,B013_QUAD)
       B013R_QUAD=real(B013_QUAD)
       B013I_QUAD=imag(B013_QUAD)
       call From_BijR_QUAD_to_Dbl(B013R,B013I,B013R_QUAD,B013I_QUAD)
       B013=DCMPLX(B013R,B013I)
       call B0fin_QUAD(s123_QUAD,musq_QUAD,B014_QUAD)
       B014R_QUAD=real(B014_QUAD)
       B014I_QUAD=imag(B014_QUAD)
       call From_BijR_QUAD_to_Dbl(B014R,B014I,B014R_QUAD,B014I_QUAD)
       B014=DCMPLX(B014R,B014I)
       call B0fin_QUAD(s16_QUAD,musq_QUAD,B026_QUAD)
       B026R_QUAD=real(B026_QUAD)
       B026I_QUAD=imag(B026_QUAD)
       call From_BijR_QUAD_to_Dbl(B026R,B026I,B026R_QUAD,B026I_QUAD)
       B026=DCMPLX(B026R,B026I)
       call B0fin_QUAD(s23_QUAD,musq_QUAD,B024_QUAD)
       B024R_QUAD=real(B024_QUAD)
       B024I_QUAD=imag(B024_QUAD)
       call From_BijR_QUAD_to_Dbl(B024R,B024I,B024R_QUAD,B024I_QUAD)
       B024=DCMPLX(B024R,B024I)
       call B0fin_QUAD(s234_QUAD,musq_QUAD,B025_QUAD)
       B025R_QUAD=real(B025_QUAD)
       B025I_QUAD=imag(B025_QUAD)
       call From_BijR_QUAD_to_Dbl(B025R,B025I,B025R_QUAD,B025I_QUAD)
       B025=DCMPLX(B025R,B025I)
       call B0fin_QUAD(s34_QUAD,musq_QUAD,B035_QUAD)
       B035R_QUAD=real(B035_QUAD)
       B035I_QUAD=imag(B035_QUAD)
       call From_BijR_QUAD_to_Dbl(B035R,B035I,B035R_QUAD,B035I_QUAD)
       B035=DCMPLX(B035R,B035I)
       call B0fin_QUAD(s345_QUAD,musq_QUAD,B036_QUAD)
       B036R_QUAD=real(B036_QUAD)
       B036I_QUAD=imag(B036_QUAD)
       call From_BijR_QUAD_to_Dbl(B036R,B036I,B036R_QUAD,B036I_QUAD)
       B036=DCMPLX(B036R,B036I)
       call B0fin_QUAD(s45_QUAD,musq_QUAD,B046_QUAD)
       B046R_QUAD=real(B046_QUAD)
       B046I_QUAD=imag(B046_QUAD)
       call From_BijR_QUAD_to_Dbl(B046R,B046I,B046R_QUAD,B046I_QUAD)
       B046=DCMPLX(B046R,B046I)
       call B0fin_QUAD(s56_QUAD,musq_QUAD,B015_QUAD)
       B015R_QUAD=real(B015_QUAD)
       B015I_QUAD=imag(B015_QUAD)
       call From_BijR_QUAD_to_Dbl(B015R,B015I,B015R_QUAD,B015I_QUAD)
       B015=DCMPLX(B015R,B015I)
       call C0fin_QUAD(p1sq_QUAD,p2sq_QUAD,s12_QUAD,musq_QUAD,C0123_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,s12_QUAD,B02
     -   3_QUAD,B013_QUAD,B012_QUAD,C0123_QUAD,C0123R_QUAD,C0123I_QUA
     -   D,C123R_QUAD,C123I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0123R,C0123I,C123R,C123I,C0123R_QU
     -   AD,C0123I_QUAD,C123R_QUAD,C123I_QUAD)
       C0123=DCMPLX(C0123R,C0123I)
       call C0fin_QUAD(p1sq_QUAD,s16_QUAD,p6sq_QUAD,musq_QUAD,C0126_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p1sq_QUAD,s16_QUAD,p6sq_QUAD,B02
     -   6_QUAD,B016_QUAD,B012_QUAD,C0126_QUAD,C0126R_QUAD,C0126I_QUA
     -   D,C126R_QUAD,C126I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0126R,C0126I,C126R,C126I,C0126R_QU
     -   AD,C0126I_QUAD,C126R_QUAD,C126I_QUAD)
       C0126=DCMPLX(C0126R,C0126I)
       call C0fin_QUAD(p1sq_QUAD,s23_QUAD,s123_QUAD,musq_QUAD,C0124_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p1sq_QUAD,s23_QUAD,s123_QUAD,B02
     -   4_QUAD,B014_QUAD,B012_QUAD,C0124_QUAD,C0124R_QUAD,C0124I_QUA
     -   D,C124R_QUAD,C124I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0124R,C0124I,C124R,C124I,C0124R_QU
     -   AD,C0124I_QUAD,C124R_QUAD,C124I_QUAD)
       C0124=DCMPLX(C0124R,C0124I)
       call C0fin_QUAD(p1sq_QUAD,s234_QUAD,s56_QUAD,musq_QUAD,C0125_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p1sq_QUAD,s234_QUAD,s56_QUAD,B02
     -   5_QUAD,B015_QUAD,B012_QUAD,C0125_QUAD,C0125R_QUAD,C0125I_QUA
     -   D,C125R_QUAD,C125I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0125R,C0125I,C125R,C125I,C0125R_QU
     -   AD,C0125I_QUAD,C125R_QUAD,C125I_QUAD)
       C0125=DCMPLX(C0125R,C0125I)
       call C0fin_QUAD(p2sq_QUAD,p3sq_QUAD,s23_QUAD,musq_QUAD,C0234_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p2sq_QUAD,p3sq_QUAD,s23_QUAD,B03
     -   4_QUAD,B024_QUAD,B023_QUAD,C0234_QUAD,C0234R_QUAD,C0234I_QUA
     -   D,C234R_QUAD,C234I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0234R,C0234I,C234R,C234I,C0234R_QU
     -   AD,C0234I_QUAD,C234R_QUAD,C234I_QUAD)
       C0234=DCMPLX(C0234R,C0234I)
       call C0fin_QUAD(p2sq_QUAD,s34_QUAD,s234_QUAD,musq_QUAD,C0235_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p2sq_QUAD,s34_QUAD,s234_QUAD,B03
     -   5_QUAD,B025_QUAD,B023_QUAD,C0235_QUAD,C0235R_QUAD,C0235I_QUA
     -   D,C235R_QUAD,C235I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0235R,C0235I,C235R,C235I,C0235R_QU
     -   AD,C0235I_QUAD,C235R_QUAD,C235I_QUAD)
       C0235=DCMPLX(C0235R,C0235I)
       call C0fin_QUAD(p2sq_QUAD,s345_QUAD,s16_QUAD,musq_QUAD,C0236_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p2sq_QUAD,s345_QUAD,s16_QUAD,B03
     -   6_QUAD,B026_QUAD,B023_QUAD,C0236_QUAD,C0236R_QUAD,C0236I_QUA
     -   D,C236R_QUAD,C236I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0236R,C0236I,C236R,C236I,C0236R_QU
     -   AD,C0236I_QUAD,C236R_QUAD,C236I_QUAD)
       C0236=DCMPLX(C0236R,C0236I)
       call C0fin_QUAD(p3sq_QUAD,p4sq_QUAD,s34_QUAD,musq_QUAD,C0345_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p3sq_QUAD,p4sq_QUAD,s34_QUAD,B04
     -   5_QUAD,B035_QUAD,B034_QUAD,C0345_QUAD,C0345R_QUAD,C0345I_QUA
     -   D,C345R_QUAD,C345I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0345R,C0345I,C345R,C345I,C0345R_QU
     -   AD,C0345I_QUAD,C345R_QUAD,C345I_QUAD)
       C0345=DCMPLX(C0345R,C0345I)
       call C0fin_QUAD(p3sq_QUAD,s45_QUAD,s345_QUAD,musq_QUAD,C0346_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p3sq_QUAD,s45_QUAD,s345_QUAD,B04
     -   6_QUAD,B036_QUAD,B034_QUAD,C0346_QUAD,C0346R_QUAD,C0346I_QUA
     -   D,C346R_QUAD,C346I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0346R,C0346I,C346R,C346I,C0346R_QU
     -   AD,C0346I_QUAD,C346R_QUAD,C346I_QUAD)
       C0346=DCMPLX(C0346R,C0346I)
       call C0fin_QUAD(p4sq_QUAD,p5sq_QUAD,s45_QUAD,musq_QUAD,C0456_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(p4sq_QUAD,p5sq_QUAD,s45_QUAD,B05
     -   6_QUAD,B046_QUAD,B045_QUAD,C0456_QUAD,C0456R_QUAD,C0456I_QUA
     -   D,C456R_QUAD,C456I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0456R,C0456I,C456R,C456I,C0456R_QU
     -   AD,C0456I_QUAD,C456R_QUAD,C456I_QUAD)
       C0456=DCMPLX(C0456R,C0456I)
       call C0fin_QUAD(s12_QUAD,p3sq_QUAD,s123_QUAD,musq_QUAD,C0134_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s12_QUAD,p3sq_QUAD,s123_QUAD,B03
     -   4_QUAD,B014_QUAD,B013_QUAD,C0134_QUAD,C0134R_QUAD,C0134I_QUA
     -   D,C134R_QUAD,C134I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0134R,C0134I,C134R,C134I,C0134R_QU
     -   AD,C0134I_QUAD,C134R_QUAD,C134I_QUAD)
       C0134=DCMPLX(C0134R,C0134I)
       call C0fin_QUAD(s12_QUAD,s34_QUAD,s56_QUAD,musq_QUAD,C0135_QUAD
     -   )
       call tens_red3_new_Re_Com_QUAD(s12_QUAD,s34_QUAD,s56_QUAD,B035_
     -   QUAD,B015_QUAD,B013_QUAD,C0135_QUAD,C0135R_QUAD,C0135I_QUAD,
     -   C135R_QUAD,C135I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0135R,C0135I,C135R,C135I,C0135R_QU
     -   AD,C0135I_QUAD,C135R_QUAD,C135I_QUAD)
       C0135=DCMPLX(C0135R,C0135I)
       call C0fin_QUAD(s12_QUAD,s345_QUAD,p6sq_QUAD,musq_QUAD,C0136_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s12_QUAD,s345_QUAD,p6sq_QUAD,B03
     -   6_QUAD,B016_QUAD,B013_QUAD,C0136_QUAD,C0136R_QUAD,C0136I_QUA
     -   D,C136R_QUAD,C136I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0136R,C0136I,C136R,C136I,C0136R_QU
     -   AD,C0136I_QUAD,C136R_QUAD,C136I_QUAD)
       C0136=DCMPLX(C0136R,C0136I)
       call C0fin_QUAD(s123_QUAD,p4sq_QUAD,s56_QUAD,musq_QUAD,C0145_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s123_QUAD,p4sq_QUAD,s56_QUAD,B04
     -   5_QUAD,B015_QUAD,B014_QUAD,C0145_QUAD,C0145R_QUAD,C0145I_QUA
     -   D,C145R_QUAD,C145I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0145R,C0145I,C145R,C145I,C0145R_QU
     -   AD,C0145I_QUAD,C145R_QUAD,C145I_QUAD)
       C0145=DCMPLX(C0145R,C0145I)
       call C0fin_QUAD(s123_QUAD,s45_QUAD,p6sq_QUAD,musq_QUAD,C0146_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s123_QUAD,s45_QUAD,p6sq_QUAD,B04
     -   6_QUAD,B016_QUAD,B014_QUAD,C0146_QUAD,C0146R_QUAD,C0146I_QUA
     -   D,C146R_QUAD,C146I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0146R,C0146I,C146R,C146I,C0146R_QU
     -   AD,C0146I_QUAD,C146R_QUAD,C146I_QUAD)
       C0146=DCMPLX(C0146R,C0146I)
       call C0fin_QUAD(s23_QUAD,p4sq_QUAD,s234_QUAD,musq_QUAD,C0245_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s23_QUAD,p4sq_QUAD,s234_QUAD,B04
     -   5_QUAD,B025_QUAD,B024_QUAD,C0245_QUAD,C0245R_QUAD,C0245I_QUA
     -   D,C245R_QUAD,C245I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0245R,C0245I,C245R,C245I,C0245R_QU
     -   AD,C0245I_QUAD,C245R_QUAD,C245I_QUAD)
       C0245=DCMPLX(C0245R,C0245I)
       call C0fin_QUAD(s23_QUAD,s45_QUAD,s16_QUAD,musq_QUAD,C0246_QUAD
     -   )
       call tens_red3_new_Re_Com_QUAD(s23_QUAD,s45_QUAD,s16_QUAD,B046_
     -   QUAD,B026_QUAD,B024_QUAD,C0246_QUAD,C0246R_QUAD,C0246I_QUAD,
     -   C246R_QUAD,C246I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0246R,C0246I,C246R,C246I,C0246R_QU
     -   AD,C0246I_QUAD,C246R_QUAD,C246I_QUAD)
       C0246=DCMPLX(C0246R,C0246I)
       call C0fin_QUAD(s234_QUAD,p5sq_QUAD,s16_QUAD,musq_QUAD,C0256_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s234_QUAD,p5sq_QUAD,s16_QUAD,B05
     -   6_QUAD,B026_QUAD,B025_QUAD,C0256_QUAD,C0256R_QUAD,C0256I_QUA
     -   D,C256R_QUAD,C256I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0256R,C0256I,C256R,C256I,C0256R_QU
     -   AD,C0256I_QUAD,C256R_QUAD,C256I_QUAD)
       C0256=DCMPLX(C0256R,C0256I)
       call C0fin_QUAD(s34_QUAD,p5sq_QUAD,s345_QUAD,musq_QUAD,C0356_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s34_QUAD,p5sq_QUAD,s345_QUAD,B05
     -   6_QUAD,B036_QUAD,B035_QUAD,C0356_QUAD,C0356R_QUAD,C0356I_QUA
     -   D,C356R_QUAD,C356I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0356R,C0356I,C356R,C356I,C0356R_QU
     -   AD,C0356I_QUAD,C356R_QUAD,C356I_QUAD)
       C0356=DCMPLX(C0356R,C0356I)
       call C0fin_QUAD(s56_QUAD,p5sq_QUAD,p6sq_QUAD,musq_QUAD,C0156_QU
     -   AD)
       call tens_red3_new_Re_Com_QUAD(s56_QUAD,p5sq_QUAD,p6sq_QUAD,B05
     -   6_QUAD,B016_QUAD,B015_QUAD,C0156_QUAD,C0156R_QUAD,C0156I_QUA
     -   D,C156R_QUAD,C156I_QUAD)
       call  From_CijR_QUAD_to_Dbl(C0156R,C0156I,C156R,C156I,C0156R_QU
     -   AD,C0156I_QUAD,C156R_QUAD,C156I_QUAD)
       C0156=DCMPLX(C0156R,C0156I)
       call D0fin_QUAD(s12_QUAD,s23_QUAD,p1sq_QUAD,p2sq_QUAD,p3sq_QUAD
     -   ,s123_QUAD,musq_QUAD,D01234_QUAD)
       call tens_red4_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p1
     -   p2_QUAD,p1p3_QUAD,p2p3_QUAD,C0234R_QUAD,C0134R_QUAD,C0124R_Q
     -   UAD,C0123R_QUAD,C234R_QUAD,C134R_QUAD,C124R_QUAD,C123R_QUAD,
     -   C0234I_QUAD,C0134I_QUAD,C0124I_QUAD,C0123I_QUAD,C234I_QUAD,C
     -   134I_QUAD,C124I_QUAD,C123I_QUAD,D01234_QUAD,D01234R_QUAD,D01
     -   234I_QUAD,D1234R_QUAD,D1234I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01234R,D01234I,D1234R,D1234I,D0123
     -   4R_QUAD,D01234I_QUAD,D1234R_QUAD,D1234I_QUAD)
       D01234=DCMPLX(D01234R,D01234I)
       call D0fin_QUAD(s12_QUAD,s234_QUAD,p1sq_QUAD,p2sq_QUAD,s34_QUAD
     -   ,s56_QUAD,musq_QUAD,D01235_QUAD)
       call tens_red4_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,s34_QUAD,p1p
     -   2_QUAD,p1p3_QUAD+ p1p4_QUAD,p2p3_QUAD+ p2p4_QUAD,C0235R_QUAD
     -   ,C0135R_QUAD,C0125R_QUAD,C0123R_QUAD,C235R_QUAD,C135R_QUAD,C
     -   125R_QUAD,C123R_QUAD,C0235I_QUAD,C0135I_QUAD,C0125I_QUAD,C01
     -   23I_QUAD,C235I_QUAD,C135I_QUAD,C125I_QUAD,C123I_QUAD,D01235_
     -   QUAD,D01235R_QUAD,D01235I_QUAD,D1235R_QUAD,D1235I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01235R,D01235I,D1235R,D1235I,D0123
     -   5R_QUAD,D01235I_QUAD,D1235R_QUAD,D1235I_QUAD)
       D01235=DCMPLX(D01235R,D01235I)
       call D0fin_QUAD(s12_QUAD,s16_QUAD,p1sq_QUAD,p2sq_QUAD,s345_QUAD
     -   ,p6sq_QUAD,musq_QUAD,D01236_QUAD)
       call tens_red4_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,s345_QUAD,p1
     -   p2_QUAD,p1p3_QUAD+ p1p4_QUAD+ p1p5_QUAD,p2p3_QUAD+ p2p4_QUAD
     -   + p2p5_QUAD,C0236R_QUAD,C0136R_QUAD,C0126R_QUAD,C0123R_QUAD,
     -   C236R_QUAD,C136R_QUAD,C126R_QUAD,C123R_QUAD,C0236I_QUAD,C013
     -   6I_QUAD,C0126I_QUAD,C0123I_QUAD,C236I_QUAD,C136I_QUAD,C126I_
     -   QUAD,C123I_QUAD,D01236_QUAD,D01236R_QUAD,D01236I_QUAD,D1236R
     -   _QUAD,D1236I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01236R,D01236I,D1236R,D1236I,D0123
     -   6R_QUAD,D01236I_QUAD,D1236R_QUAD,D1236I_QUAD)
       D01236=DCMPLX(D01236R,D01236I)
       call D0fin_QUAD(s123_QUAD,s234_QUAD,p1sq_QUAD,s23_QUAD,p4sq_QUA
     -   D,s56_QUAD,musq_QUAD,D01245_QUAD)
       call tens_red4_new_Re_Com_QUAD(p1sq_QUAD,s23_QUAD,p4sq_QUAD,p1p
     -   2_QUAD+ p1p3_QUAD,p1p4_QUAD,p2p4_QUAD+ p3p4_QUAD,C0245R_QUAD
     -   ,C0145R_QUAD,C0125R_QUAD,C0124R_QUAD,C245R_QUAD,C145R_QUAD,C
     -   125R_QUAD,C124R_QUAD,C0245I_QUAD,C0145I_QUAD,C0125I_QUAD,C01
     -   24I_QUAD,C245I_QUAD,C145I_QUAD,C125I_QUAD,C124I_QUAD,D01245_
     -   QUAD,D01245R_QUAD,D01245I_QUAD,D1245R_QUAD,D1245I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01245R,D01245I,D1245R,D1245I,D0124
     -   5R_QUAD,D01245I_QUAD,D1245R_QUAD,D1245I_QUAD)
       D01245=DCMPLX(D01245R,D01245I)
       call D0fin_QUAD(s123_QUAD,s16_QUAD,p1sq_QUAD,s23_QUAD,s45_QUAD,
     -   p6sq_QUAD,musq_QUAD,D01246_QUAD)
       call tens_red4_new_Re_Com_QUAD(p1sq_QUAD,s23_QUAD,s45_QUAD,p1p2
     -   _QUAD+ p1p3_QUAD,p1p4_QUAD+ p1p5_QUAD,p2p4_QUAD+ p2p5_QUAD+ 
     -   p3p4_QUAD+ p3p5_QUAD,C0246R_QUAD,C0146R_QUAD,C0126R_QUAD,C01
     -   24R_QUAD,C246R_QUAD,C146R_QUAD,C126R_QUAD,C124R_QUAD,C0246I_
     -   QUAD,C0146I_QUAD,C0126I_QUAD,C0124I_QUAD,C246I_QUAD,C146I_QU
     -   AD,C126I_QUAD,C124I_QUAD,D01246_QUAD,D01246R_QUAD,D01246I_QU
     -   AD,D1246R_QUAD,D1246I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01246R,D01246I,D1246R,D1246I,D0124
     -   6R_QUAD,D01246I_QUAD,D1246R_QUAD,D1246I_QUAD)
       D01246=DCMPLX(D01246R,D01246I)
       call D0fin_QUAD(s56_QUAD,s16_QUAD,p1sq_QUAD,s234_QUAD,p5sq_QUAD
     -   ,p6sq_QUAD,musq_QUAD,D01256_QUAD)
       call tens_red4_new_Re_Com_QUAD(p1sq_QUAD,s234_QUAD,p5sq_QUAD,p1
     -   p2_QUAD+ p1p3_QUAD+ p1p4_QUAD,p1p5_QUAD,p2p5_QUAD+ p3p5_QUAD
     -   + p4p5_QUAD,C0256R_QUAD,C0156R_QUAD,C0126R_QUAD,C0125R_QUAD,
     -   C256R_QUAD,C156R_QUAD,C126R_QUAD,C125R_QUAD,C0256I_QUAD,C015
     -   6I_QUAD,C0126I_QUAD,C0125I_QUAD,C256I_QUAD,C156I_QUAD,C126I_
     -   QUAD,C125I_QUAD,D01256_QUAD,D01256R_QUAD,D01256I_QUAD,D1256R
     -   _QUAD,D1256I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01256R,D01256I,D1256R,D1256I,D0125
     -   6R_QUAD,D01256I_QUAD,D1256R_QUAD,D1256I_QUAD)
       D01256=DCMPLX(D01256R,D01256I)
       call D0fin_QUAD(s23_QUAD,s34_QUAD,p2sq_QUAD,p3sq_QUAD,p4sq_QUAD
     -   ,s234_QUAD,musq_QUAD,D02345_QUAD)
       call tens_red4_new_Re_Com_QUAD(p2sq_QUAD,p3sq_QUAD,p4sq_QUAD,p2
     -   p3_QUAD,p2p4_QUAD,p3p4_QUAD,C0345R_QUAD,C0245R_QUAD,C0235R_Q
     -   UAD,C0234R_QUAD,C345R_QUAD,C245R_QUAD,C235R_QUAD,C234R_QUAD,
     -   C0345I_QUAD,C0245I_QUAD,C0235I_QUAD,C0234I_QUAD,C345I_QUAD,C
     -   245I_QUAD,C235I_QUAD,C234I_QUAD,D02345_QUAD,D02345R_QUAD,D02
     -   345I_QUAD,D2345R_QUAD,D2345I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D02345R,D02345I,D2345R,D2345I,D0234
     -   5R_QUAD,D02345I_QUAD,D2345R_QUAD,D2345I_QUAD)
       D02345=DCMPLX(D02345R,D02345I)
       call D0fin_QUAD(s23_QUAD,s345_QUAD,p2sq_QUAD,p3sq_QUAD,s45_QUAD
     -   ,s16_QUAD,musq_QUAD,D02346_QUAD)
       call tens_red4_new_Re_Com_QUAD(p2sq_QUAD,p3sq_QUAD,s45_QUAD,p2p
     -   3_QUAD,p2p4_QUAD+ p2p5_QUAD,p3p4_QUAD+ p3p5_QUAD,C0346R_QUAD
     -   ,C0246R_QUAD,C0236R_QUAD,C0234R_QUAD,C346R_QUAD,C246R_QUAD,C
     -   236R_QUAD,C234R_QUAD,C0346I_QUAD,C0246I_QUAD,C0236I_QUAD,C02
     -   34I_QUAD,C346I_QUAD,C246I_QUAD,C236I_QUAD,C234I_QUAD,D02346_
     -   QUAD,D02346R_QUAD,D02346I_QUAD,D2346R_QUAD,D2346I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D02346R,D02346I,D2346R,D2346I,D0234
     -   6R_QUAD,D02346I_QUAD,D2346R_QUAD,D2346I_QUAD)
       D02346=DCMPLX(D02346R,D02346I)
       call D0fin_QUAD(s234_QUAD,s345_QUAD,p2sq_QUAD,s34_QUAD,p5sq_QUA
     -   D,s16_QUAD,musq_QUAD,D02356_QUAD)
       call tens_red4_new_Re_Com_QUAD(p2sq_QUAD,s34_QUAD,p5sq_QUAD,p2p
     -   3_QUAD+ p2p4_QUAD,p2p5_QUAD,p3p5_QUAD+ p4p5_QUAD,C0356R_QUAD
     -   ,C0256R_QUAD,C0236R_QUAD,C0235R_QUAD,C356R_QUAD,C256R_QUAD,C
     -   236R_QUAD,C235R_QUAD,C0356I_QUAD,C0256I_QUAD,C0236I_QUAD,C02
     -   35I_QUAD,C356I_QUAD,C256I_QUAD,C236I_QUAD,C235I_QUAD,D02356_
     -   QUAD,D02356R_QUAD,D02356I_QUAD,D2356R_QUAD,D2356I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D02356R,D02356I,D2356R,D2356I,D0235
     -   6R_QUAD,D02356I_QUAD,D2356R_QUAD,D2356I_QUAD)
       D02356=DCMPLX(D02356R,D02356I)
       call D0fin_QUAD(s34_QUAD,s45_QUAD,p3sq_QUAD,p4sq_QUAD,p5sq_QUAD
     -   ,s345_QUAD,musq_QUAD,D03456_QUAD)
       call tens_red4_new_Re_Com_QUAD(p3sq_QUAD,p4sq_QUAD,p5sq_QUAD,p3
     -   p4_QUAD,p3p5_QUAD,p4p5_QUAD,C0456R_QUAD,C0356R_QUAD,C0346R_Q
     -   UAD,C0345R_QUAD,C456R_QUAD,C356R_QUAD,C346R_QUAD,C345R_QUAD,
     -   C0456I_QUAD,C0356I_QUAD,C0346I_QUAD,C0345I_QUAD,C456I_QUAD,C
     -   356I_QUAD,C346I_QUAD,C345I_QUAD,D03456_QUAD,D03456R_QUAD,D03
     -   456I_QUAD,D3456R_QUAD,D3456I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D03456R,D03456I,D3456R,D3456I,D0345
     -   6R_QUAD,D03456I_QUAD,D3456R_QUAD,D3456I_QUAD)
       D03456=DCMPLX(D03456R,D03456I)
       call D0fin_QUAD(s123_QUAD,s34_QUAD,s12_QUAD,p3sq_QUAD,p4sq_QUAD
     -   ,s56_QUAD,musq_QUAD,D01345_QUAD)
       call tens_red4_new_Re_Com_QUAD(s12_QUAD,p3sq_QUAD,p4sq_QUAD,p1p
     -   3_QUAD+ p2p3_QUAD,p1p4_QUAD+ p2p4_QUAD,p3p4_QUAD,C0345R_QUAD
     -   ,C0145R_QUAD,C0135R_QUAD,C0134R_QUAD,C345R_QUAD,C145R_QUAD,C
     -   135R_QUAD,C134R_QUAD,C0345I_QUAD,C0145I_QUAD,C0135I_QUAD,C01
     -   34I_QUAD,C345I_QUAD,C145I_QUAD,C135I_QUAD,C134I_QUAD,D01345_
     -   QUAD,D01345R_QUAD,D01345I_QUAD,D1345R_QUAD,D1345I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01345R,D01345I,D1345R,D1345I,D0134
     -   5R_QUAD,D01345I_QUAD,D1345R_QUAD,D1345I_QUAD)
       D01345=DCMPLX(D01345R,D01345I)
       call D0fin_QUAD(s123_QUAD,s345_QUAD,s12_QUAD,p3sq_QUAD,s45_QUAD
     -   ,p6sq_QUAD,musq_QUAD,D01346_QUAD)
       call tens_red4_new_Re_Com_QUAD(s12_QUAD,p3sq_QUAD,s45_QUAD,p1p3
     -   _QUAD+ p2p3_QUAD,p1p4_QUAD+ p1p5_QUAD+ p2p4_QUAD+ p2p5_QUAD,
     -   p3p4_QUAD+ p3p5_QUAD,C0346R_QUAD,C0146R_QUAD,C0136R_QUAD,C01
     -   34R_QUAD,C346R_QUAD,C146R_QUAD,C136R_QUAD,C134R_QUAD,C0346I_
     -   QUAD,C0146I_QUAD,C0136I_QUAD,C0134I_QUAD,C346I_QUAD,C146I_QU
     -   AD,C136I_QUAD,C134I_QUAD,D01346_QUAD,D01346R_QUAD,D01346I_QU
     -   AD,D1346R_QUAD,D1346I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01346R,D01346I,D1346R,D1346I,D0134
     -   6R_QUAD,D01346I_QUAD,D1346R_QUAD,D1346I_QUAD)
       D01346=DCMPLX(D01346R,D01346I)
       call D0fin_QUAD(s56_QUAD,s345_QUAD,s12_QUAD,s34_QUAD,p5sq_QUAD,
     -   p6sq_QUAD,musq_QUAD,D01356_QUAD)
       call tens_red4_new_Re_Com_QUAD(s12_QUAD,s34_QUAD,p5sq_QUAD,p1p3
     -   _QUAD+ p1p4_QUAD+ p2p3_QUAD+ p2p4_QUAD,p1p5_QUAD+ p2p5_QUAD,
     -   p3p5_QUAD+ p4p5_QUAD,C0356R_QUAD,C0156R_QUAD,C0136R_QUAD,C01
     -   35R_QUAD,C356R_QUAD,C156R_QUAD,C136R_QUAD,C135R_QUAD,C0356I_
     -   QUAD,C0156I_QUAD,C0136I_QUAD,C0135I_QUAD,C356I_QUAD,C156I_QU
     -   AD,C136I_QUAD,C135I_QUAD,D01356_QUAD,D01356R_QUAD,D01356I_QU
     -   AD,D1356R_QUAD,D1356I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01356R,D01356I,D1356R,D1356I,D0135
     -   6R_QUAD,D01356I_QUAD,D1356R_QUAD,D1356I_QUAD)
       D01356=DCMPLX(D01356R,D01356I)
       call D0fin_QUAD(s56_QUAD,s45_QUAD,s123_QUAD,p4sq_QUAD,p5sq_QUAD
     -   ,p6sq_QUAD,musq_QUAD,D01456_QUAD)
       call tens_red4_new_Re_Com_QUAD(s123_QUAD,p4sq_QUAD,p5sq_QUAD,p1
     -   p4_QUAD+ p2p4_QUAD+ p3p4_QUAD,p1p5_QUAD+ p2p5_QUAD+ p3p5_QUA
     -   D,p4p5_QUAD,C0456R_QUAD,C0156R_QUAD,C0146R_QUAD,C0145R_QUAD,
     -   C456R_QUAD,C156R_QUAD,C146R_QUAD,C145R_QUAD,C0456I_QUAD,C015
     -   6I_QUAD,C0146I_QUAD,C0145I_QUAD,C456I_QUAD,C156I_QUAD,C146I_
     -   QUAD,C145I_QUAD,D01456_QUAD,D01456R_QUAD,D01456I_QUAD,D1456R
     -   _QUAD,D1456I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D01456R,D01456I,D1456R,D1456I,D0145
     -   6R_QUAD,D01456I_QUAD,D1456R_QUAD,D1456I_QUAD)
       D01456=DCMPLX(D01456R,D01456I)
       call D0fin_QUAD(s234_QUAD,s45_QUAD,s23_QUAD,p4sq_QUAD,p5sq_QUAD
     -   ,s16_QUAD,musq_QUAD,D02456_QUAD)
       call tens_red4_new_Re_Com_QUAD(s23_QUAD,p4sq_QUAD,p5sq_QUAD,p2p
     -   4_QUAD+ p3p4_QUAD,p2p5_QUAD+ p3p5_QUAD,p4p5_QUAD,C0456R_QUAD
     -   ,C0256R_QUAD,C0246R_QUAD,C0245R_QUAD,C456R_QUAD,C256R_QUAD,C
     -   246R_QUAD,C245R_QUAD,C0456I_QUAD,C0256I_QUAD,C0246I_QUAD,C02
     -   45I_QUAD,C456I_QUAD,C256I_QUAD,C246I_QUAD,C245I_QUAD,D02456_
     -   QUAD,D02456R_QUAD,D02456I_QUAD,D2456R_QUAD,D2456I_QUAD)
       call  From_DijR_QUAD_to_Dbl(D02456R,D02456I,D2456R,D2456I,D0245
     -   6R_QUAD,D02456I_QUAD,D2456R_QUAD,D2456I_QUAD)
       D02456=DCMPLX(D02456R,D02456I)
       call E0fin_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p4sq_QUAD,s56_QUA
     -   D,s12_QUAD,s23_QUAD,s34_QUAD,s123_QUAD,s234_QUAD,D02345_QUAD
     -   ,D01345_QUAD,D01245_QUAD,D01235_QUAD,D01234_QUAD,E012345_QUA
     -   D)
       E012345R_QUAD=real(E012345_QUAD)
       E012345I_QUAD=imag(E012345_QUAD)
       call tens_red5_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p4
     -   sq_QUAD,p1p2_QUAD,p1p3_QUAD,p1p4_QUAD,p2p3_QUAD,p2p4_QUAD,p3
     -   p4_QUAD,D02345R_QUAD,D01345R_QUAD,D01245R_QUAD,D01235R_QUAD,
     -   D01234R_QUAD,D2345R_QUAD,D1345R_QUAD,D1245R_QUAD,D1235R_QUAD
     -   ,D1234R_QUAD,D02345I_QUAD,D01345I_QUAD,D01245I_QUAD,D01235I_
     -   QUAD,D01234I_QUAD,D2345I_QUAD,D1345I_QUAD,D1245I_QUAD,D1235I
     -   _QUAD,D1234I_QUAD,E12345R_QUAD,E12345I_QUAD)
       call From_EijR_QUAD_to_Dbl(E012345R,E012345I,E12345R,E12345I,E0
     -   12345R_QUAD,E012345I_QUAD,E12345R_QUAD,E12345I_QUAD)
       E012345=DCMPLX(E012345R,E012345I)
       call E0fin_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,s45_QUAD,p6sq_QUA
     -   D,s12_QUAD,s23_QUAD,s345_QUAD,s123_QUAD,s16_QUAD,D02346_QUAD
     -   ,D01346_QUAD,D01246_QUAD,D01236_QUAD,D01234_QUAD,E012346_QUA
     -   D)
       E012346R_QUAD=real(E012346_QUAD)
       E012346I_QUAD=imag(E012346_QUAD)
       call tens_red5_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,s4
     -   5_QUAD,p1p2_QUAD,p1p3_QUAD,p1p4_QUAD+ p1p5_QUAD,p2p3_QUAD,p2
     -   p4_QUAD+ p2p5_QUAD,p3p4_QUAD+ p3p5_QUAD,D02346R_QUAD,D01346R
     -   _QUAD,D01246R_QUAD,D01236R_QUAD,D01234R_QUAD,D2346R_QUAD,D13
     -   46R_QUAD,D1246R_QUAD,D1236R_QUAD,D1234R_QUAD,D02346I_QUAD,D0
     -   1346I_QUAD,D01246I_QUAD,D01236I_QUAD,D01234I_QUAD,D2346I_QUA
     -   D,D1346I_QUAD,D1246I_QUAD,D1236I_QUAD,D1234I_QUAD,E12346R_QU
     -   AD,E12346I_QUAD)
       call From_EijR_QUAD_to_Dbl(E012346R,E012346I,E12346R,E12346I,E0
     -   12346R_QUAD,E012346I_QUAD,E12346R_QUAD,E12346I_QUAD)
       E012346=DCMPLX(E012346R,E012346I)
       call E0fin_QUAD(p1sq_QUAD,p2sq_QUAD,s34_QUAD,p5sq_QUAD,p6sq_QUA
     -   D,s12_QUAD,s234_QUAD,s345_QUAD,s56_QUAD,s16_QUAD,D02356_QUAD
     -   ,D01356_QUAD,D01256_QUAD,D01236_QUAD,D01235_QUAD,E012356_QUA
     -   D)
       E012356R_QUAD=real(E012356_QUAD)
       E012356I_QUAD=imag(E012356_QUAD)
       call tens_red5_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,s34_QUAD,p5s
     -   q_QUAD,p1p2_QUAD,p1p3_QUAD+ p1p4_QUAD,p1p5_QUAD,p2p3_QUAD+ p
     -   2p4_QUAD,p2p5_QUAD,p3p5_QUAD+ p4p5_QUAD,D02356R_QUAD,D01356R
     -   _QUAD,D01256R_QUAD,D01236R_QUAD,D01235R_QUAD,D2356R_QUAD,D13
     -   56R_QUAD,D1256R_QUAD,D1236R_QUAD,D1235R_QUAD,D02356I_QUAD,D0
     -   1356I_QUAD,D01256I_QUAD,D01236I_QUAD,D01235I_QUAD,D2356I_QUA
     -   D,D1356I_QUAD,D1256I_QUAD,D1236I_QUAD,D1235I_QUAD,E12356R_QU
     -   AD,E12356I_QUAD)
       call From_EijR_QUAD_to_Dbl(E012356R,E012356I,E12356R,E12356I,E0
     -   12356R_QUAD,E012356I_QUAD,E12356R_QUAD,E12356I_QUAD)
       E012356=DCMPLX(E012356R,E012356I)
       call E0fin_QUAD(p1sq_QUAD,s23_QUAD,p4sq_QUAD,p5sq_QUAD,p6sq_QUA
     -   D,s123_QUAD,s234_QUAD,s45_QUAD,s56_QUAD,s16_QUAD,D02456_QUAD
     -   ,D01456_QUAD,D01256_QUAD,D01246_QUAD,D01245_QUAD,E012456_QUA
     -   D)
       E012456R_QUAD=real(E012456_QUAD)
       E012456I_QUAD=imag(E012456_QUAD)
       call tens_red5_new_Re_Com_QUAD(p1sq_QUAD,s23_QUAD,p4sq_QUAD,p5s
     -   q_QUAD,p1p2_QUAD+ p1p3_QUAD,p1p4_QUAD,p1p5_QUAD,p2p4_QUAD+ p
     -   3p4_QUAD,p2p5_QUAD+ p3p5_QUAD,p4p5_QUAD,D02456R_QUAD,D01456R
     -   _QUAD,D01256R_QUAD,D01246R_QUAD,D01245R_QUAD,D2456R_QUAD,D14
     -   56R_QUAD,D1256R_QUAD,D1246R_QUAD,D1245R_QUAD,D02456I_QUAD,D0
     -   1456I_QUAD,D01256I_QUAD,D01246I_QUAD,D01245I_QUAD,D2456I_QUA
     -   D,D1456I_QUAD,D1256I_QUAD,D1246I_QUAD,D1245I_QUAD,E12456R_QU
     -   AD,E12456I_QUAD)
       call From_EijR_QUAD_to_Dbl(E012456R,E012456I,E12456R,E12456I,E0
     -   12456R_QUAD,E012456I_QUAD,E12456R_QUAD,E12456I_QUAD)
       E012456=DCMPLX(E012456R,E012456I)
       call E0fin_QUAD(p2sq_QUAD,p3sq_QUAD,p4sq_QUAD,p5sq_QUAD,s16_QUA
     -   D,s23_QUAD,s34_QUAD,s45_QUAD,s234_QUAD,s345_QUAD,D03456_QUAD
     -   ,D02456_QUAD,D02356_QUAD,D02346_QUAD,D02345_QUAD,E023456_QUA
     -   D)
       E023456R_QUAD=real(E023456_QUAD)
       E023456I_QUAD=imag(E023456_QUAD)
       call tens_red5_new_Re_Com_QUAD(p2sq_QUAD,p3sq_QUAD,p4sq_QUAD,p5
     -   sq_QUAD,p2p3_QUAD,p2p4_QUAD,p2p5_QUAD,p3p4_QUAD,p3p5_QUAD,p4
     -   p5_QUAD,D03456R_QUAD,D02456R_QUAD,D02356R_QUAD,D02346R_QUAD,
     -   D02345R_QUAD,D3456R_QUAD,D2456R_QUAD,D2356R_QUAD,D2346R_QUAD
     -   ,D2345R_QUAD,D03456I_QUAD,D02456I_QUAD,D02356I_QUAD,D02346I_
     -   QUAD,D02345I_QUAD,D3456I_QUAD,D2456I_QUAD,D2356I_QUAD,D2346I
     -   _QUAD,D2345I_QUAD,E23456R_QUAD,E23456I_QUAD)
       call From_EijR_QUAD_to_Dbl(E023456R,E023456I,E23456R,E23456I,E0
     -   23456R_QUAD,E023456I_QUAD,E23456R_QUAD,E23456I_QUAD)
       E023456=DCMPLX(E023456R,E023456I)
       call E0fin_QUAD(s12_QUAD,p3sq_QUAD,p4sq_QUAD,p5sq_QUAD,p6sq_QUA
     -   D,s123_QUAD,s34_QUAD,s45_QUAD,s56_QUAD,s345_QUAD,D03456_QUAD
     -   ,D01456_QUAD,D01356_QUAD,D01346_QUAD,D01345_QUAD,E013456_QUA
     -   D)
       E013456R_QUAD=real(E013456_QUAD)
       E013456I_QUAD=imag(E013456_QUAD)
       call tens_red5_new_Re_Com_QUAD(s12_QUAD,p3sq_QUAD,p4sq_QUAD,p5s
     -   q_QUAD,p1p3_QUAD+ p2p3_QUAD,p1p4_QUAD+ p2p4_QUAD,p1p5_QUAD+ 
     -   p2p5_QUAD,p3p4_QUAD,p3p5_QUAD,p4p5_QUAD,D03456R_QUAD,D01456R
     -   _QUAD,D01356R_QUAD,D01346R_QUAD,D01345R_QUAD,D3456R_QUAD,D14
     -   56R_QUAD,D1356R_QUAD,D1346R_QUAD,D1345R_QUAD,D03456I_QUAD,D0
     -   1456I_QUAD,D01356I_QUAD,D01346I_QUAD,D01345I_QUAD,D3456I_QUA
     -   D,D1456I_QUAD,D1356I_QUAD,D1346I_QUAD,D1345I_QUAD,E13456R_QU
     -   AD,E13456I_QUAD)
       call From_EijR_QUAD_to_Dbl(E013456R,E013456I,E13456R,E13456I,E0
     -   13456R_QUAD,E013456I_QUAD,E13456R_QUAD,E13456I_QUAD)
       E013456=DCMPLX(E013456R,E013456I)
        call F0fin_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p4sq_QUAD,p5sq_Q
     -   UAD,p6sq_QUAD,s12_QUAD,s23_QUAD,s34_QUAD,s45_QUAD,s56_QUAD,s
     -   16_QUAD,s123_QUAD,s234_QUAD,s345_QUAD,E023456_QUAD,E013456_Q
     -   UAD,E012456_QUAD,E012356_QUAD,E012346_QUAD,E012345_QUAD,F012
     -   3456_QUAD)
       F0123456R_QUAD=real(F0123456_QUAD)
       F0123456I_QUAD=imag(F0123456_QUAD)
       call tens_red6_new_Re_Com_QUAD(p1sq_QUAD,p2sq_QUAD,p3sq_QUAD,p4
     -   sq_QUAD,p5sq_QUAD,p1p2_QUAD,p1p3_QUAD,p1p4_QUAD,p1p5_QUAD,p2
     -   p3_QUAD,p2p4_QUAD,p2p5_QUAD,p3p4_QUAD,p3p5_QUAD,p4p5_QUAD,E0
     -   23456R_QUAD,E013456R_QUAD,E012456R_QUAD,E012356R_QUAD,E01234
     -   6R_QUAD,E012345R_QUAD,E23456R_QUAD,E13456R_QUAD,E12456R_QUAD
     -   ,E12356R_QUAD,E12346R_QUAD,E12345R_QUAD,E023456I_QUAD,E01345
     -   6I_QUAD,E012456I_QUAD,E012356I_QUAD,E012346I_QUAD,E012345I_Q
     -   UAD,E23456I_QUAD,E13456I_QUAD,E12456I_QUAD,E12356I_QUAD,E123
     -   46I_QUAD,E12345I_QUAD,F123456R_QUAD,F123456I_QUAD,F123456_QU
     -   AD)
       call  From_FijR_QUAD_to_Dbl(F0123456R,F0123456I,F123456R,F12345
     -   6I,F0123456R_QUAD,F0123456I_QUAD,F123456R_QUAD,F123456I_QUAD
     -   )
       F0123456=DCMPLX(F0123456R,F0123456I)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c   Printing C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************
       If(PrintB) then
       write(*,*)"B12"
       write(*,*) "B12=", DCMPLX(B012R,B012I)
       write(*,*)"B23"
       write(*,*) "B23=", DCMPLX(B023R,B023I)
       write(*,*)"B34"
       write(*,*) "B34=", DCMPLX(B034R,B034I)
       write(*,*)"B45"
       write(*,*) "B45=", DCMPLX(B045R,B045I)
       write(*,*)"B56"
       write(*,*) "B56=", DCMPLX(B056R,B056I)
       write(*,*)"B16"
       write(*,*) "B16=", DCMPLX(B016R,B016I)
       write(*,*)"B13"
       write(*,*) "B13=", DCMPLX(B013R,B013I)
       write(*,*)"B14"
       write(*,*) "B14=", DCMPLX(B014R,B014I)
       write(*,*)"B26"
       write(*,*) "B26=", DCMPLX(B026R,B026I)
       write(*,*)"B24"
       write(*,*) "B24=", DCMPLX(B024R,B024I)
       write(*,*)"B25"
       write(*,*) "B25=", DCMPLX(B025R,B025I)
       write(*,*)"B35"
       write(*,*) "B35=", DCMPLX(B035R,B035I)
       write(*,*)"B36"
       write(*,*) "B36=", DCMPLX(B036R,B036I)
       write(*,*)"B46"
       write(*,*) "B46=", DCMPLX(B046R,B046I)
       write(*,*)"B15"
       write(*,*) "B15=", DCMPLX(B015R,B015I)
      Endif    
       If(PrintC) then
       write(*,*)"C123"
       call print_CijR_QUAD(C0123_QUAD,C123R_QUAD,C123I_QUAD)
       call print_CijR_QUAD_D(C0123,C123R,C123I,C0123_QUAD,C123R_QUAD,
     -   C123I_QUAD)
       call print_CijR(C0123,C123R,C123I)
c       call tens_Int3_Re_Com(0d0,p1sq,p2sq,s12,B023,B013,B012,C0123R,C
c     -   0123I,C123R,C123I)
       write(*,*)"C126"
       call print_CijR_QUAD(C0126_QUAD,C126R_QUAD,C126I_QUAD)
       call print_CijR_QUAD_D(C0126,C126R,C126I,C0126_QUAD,C126R_QUAD,
     -   C126I_QUAD)
       call print_CijR(C0126,C126R,C126I)
c       call tens_Int3_Re_Com(0d0,p1sq,s16,p6sq,B026,B016,B012,C0126R,C
c     -   0126I,C126R,C126I)
       write(*,*)"C124"
       call print_CijR_QUAD(C0124_QUAD,C124R_QUAD,C124I_QUAD)
       call print_CijR_QUAD_D(C0124,C124R,C124I,C0124_QUAD,C124R_QUAD,
     -   C124I_QUAD)
       call print_CijR(C0124,C124R,C124I)
c       call tens_Int3_Re_Com(0d0,p1sq,s23,s123,B024,B014,B012,C0124R,C
c     -   0124I,C124R,C124I)
       write(*,*)"C125"
       call print_CijR_QUAD(C0125_QUAD,C125R_QUAD,C125I_QUAD)
       call print_CijR_QUAD_D(C0125,C125R,C125I,C0125_QUAD,C125R_QUAD,
     -   C125I_QUAD)
       call print_CijR(C0125,C125R,C125I)
c       call tens_Int3_Re_Com(0d0,p1sq,s234,s56,B025,B015,B012,C0125R,C
c     -   0125I,C125R,C125I)
       write(*,*)"C234"
       call print_CijR_QUAD(C0234_QUAD,C234R_QUAD,C234I_QUAD)
       call print_CijR_QUAD_D(C0234,C234R,C234I,C0234_QUAD,C234R_QUAD,
     -   C234I_QUAD)
       call print_CijR(C0234,C234R,C234I)
c       call tens_Int3_Re_Com(0d0,p2sq,p3sq,s23,B034,B024,B023,C0234R,C
c     -   0234I,C234R,C234I)
       write(*,*)"C235"
       call print_CijR_QUAD(C0235_QUAD,C235R_QUAD,C235I_QUAD)
       call print_CijR_QUAD_D(C0235,C235R,C235I,C0235_QUAD,C235R_QUAD,
     -   C235I_QUAD)
       call print_CijR(C0235,C235R,C235I)
c       call tens_Int3_Re_Com(0d0,p2sq,s34,s234,B035,B025,B023,C0235R,C
c     -   0235I,C235R,C235I)
       write(*,*)"C236"
       call print_CijR_QUAD(C0236_QUAD,C236R_QUAD,C236I_QUAD)
       call print_CijR_QUAD_D(C0236,C236R,C236I,C0236_QUAD,C236R_QUAD,
     -   C236I_QUAD)
       call print_CijR(C0236,C236R,C236I)
c       call tens_Int3_Re_Com(0d0,p2sq,s345,s16,B036,B026,B023,C0236R,C
c     -   0236I,C236R,C236I)
       write(*,*)"C345"
       call print_CijR_QUAD(C0345_QUAD,C345R_QUAD,C345I_QUAD)
       call print_CijR_QUAD_D(C0345,C345R,C345I,C0345_QUAD,C345R_QUAD,
     -   C345I_QUAD)
       call print_CijR(C0345,C345R,C345I)
c       call tens_Int3_Re_Com(0d0,p3sq,p4sq,s34,B045,B035,B034,C0345R,C
c     -   0345I,C345R,C345I)
       write(*,*)"C346"
       call print_CijR_QUAD(C0346_QUAD,C346R_QUAD,C346I_QUAD)
       call print_CijR_QUAD_D(C0346,C346R,C346I,C0346_QUAD,C346R_QUAD,
     -   C346I_QUAD)
       call print_CijR(C0346,C346R,C346I)
c       call tens_Int3_Re_Com(0d0,p3sq,s45,s345,B046,B036,B034,C0346R,C
c     -   0346I,C346R,C346I)
       write(*,*)"C456"
       call print_CijR_QUAD(C0456_QUAD,C456R_QUAD,C456I_QUAD)
       call print_CijR_QUAD_D(C0456,C456R,C456I,C0456_QUAD,C456R_QUAD,
     -   C456I_QUAD)
       call print_CijR(C0456,C456R,C456I)
c       call tens_Int3_Re_Com(0d0,p4sq,p5sq,s45,B056,B046,B045,C0456R,C
c     -   0456I,C456R,C456I)
       write(*,*)"C134"
       call print_CijR_QUAD(C0134_QUAD,C134R_QUAD,C134I_QUAD)
       call print_CijR_QUAD_D(C0134,C134R,C134I,C0134_QUAD,C134R_QUAD,
     -   C134I_QUAD)
       call print_CijR(C0134,C134R,C134I)
c       call tens_Int3_Re_Com(0d0,s12,p3sq,s123,B034,B014,B013,C0134R,C
c     -   0134I,C134R,C134I)
       write(*,*)"C135"
       call print_CijR_QUAD(C0135_QUAD,C135R_QUAD,C135I_QUAD)
       call print_CijR_QUAD_D(C0135,C135R,C135I,C0135_QUAD,C135R_QUAD,
     -   C135I_QUAD)
       call print_CijR(C0135,C135R,C135I)
c       call tens_Int3_Re_Com(0d0,s12,s34,s56,B035,B015,B013,C0135R,C01
c     -   35I,C135R,C135I)
       write(*,*)"C136"
       call print_CijR_QUAD(C0136_QUAD,C136R_QUAD,C136I_QUAD)
       call print_CijR_QUAD_D(C0136,C136R,C136I,C0136_QUAD,C136R_QUAD,
     -   C136I_QUAD)
       call print_CijR(C0136,C136R,C136I)
c       call tens_Int3_Re_Com(0d0,s12,s345,p6sq,B036,B016,B013,C0136R,C
c     -   0136I,C136R,C136I)
       write(*,*)"C145"
       call print_CijR_QUAD(C0145_QUAD,C145R_QUAD,C145I_QUAD)
       call print_CijR_QUAD_D(C0145,C145R,C145I,C0145_QUAD,C145R_QUAD,
     -   C145I_QUAD)
       call print_CijR(C0145,C145R,C145I)
c       call tens_Int3_Re_Com(0d0,s123,p4sq,s56,B045,B015,B014,C0145R,C
c     -   0145I,C145R,C145I)
       write(*,*)"C146"
       call print_CijR_QUAD(C0146_QUAD,C146R_QUAD,C146I_QUAD)
       call print_CijR_QUAD_D(C0146,C146R,C146I,C0146_QUAD,C146R_QUAD,
     -   C146I_QUAD)
       call print_CijR(C0146,C146R,C146I)
c       call tens_Int3_Re_Com(0d0,s123,s45,p6sq,B046,B016,B014,C0146R,C
c     -   0146I,C146R,C146I)
       write(*,*)"C245"
       call print_CijR_QUAD(C0245_QUAD,C245R_QUAD,C245I_QUAD)
       call print_CijR_QUAD_D(C0245,C245R,C245I,C0245_QUAD,C245R_QUAD,
     -   C245I_QUAD)
       call print_CijR(C0245,C245R,C245I)
c       call tens_Int3_Re_Com(0d0,s23,p4sq,s234,B045,B025,B024,C0245R,C
c     -   0245I,C245R,C245I)
       write(*,*)"C246"
       call print_CijR_QUAD(C0246_QUAD,C246R_QUAD,C246I_QUAD)
       call print_CijR_QUAD_D(C0246,C246R,C246I,C0246_QUAD,C246R_QUAD,
     -   C246I_QUAD)
       call print_CijR(C0246,C246R,C246I)
c       call tens_Int3_Re_Com(0d0,s23,s45,s16,B046,B026,B024,C0246R,C02
c     -   46I,C246R,C246I)
       write(*,*)"C256"
       call print_CijR_QUAD(C0256_QUAD,C256R_QUAD,C256I_QUAD)
       call print_CijR_QUAD_D(C0256,C256R,C256I,C0256_QUAD,C256R_QUAD,
     -   C256I_QUAD)
       call print_CijR(C0256,C256R,C256I)
c       call tens_Int3_Re_Com(0d0,s234,p5sq,s16,B056,B026,B025,C0256R,C
c     -   0256I,C256R,C256I)
       write(*,*)"C356"
       call print_CijR_QUAD(C0356_QUAD,C356R_QUAD,C356I_QUAD)
       call print_CijR_QUAD_D(C0356,C356R,C356I,C0356_QUAD,C356R_QUAD,
     -   C356I_QUAD)
       call print_CijR(C0356,C356R,C356I)
c       call tens_Int3_Re_Com(0d0,s34,p5sq,s345,B056,B036,B035,C0356R,C
c     -   0356I,C356R,C356I)
       write(*,*)"C156"
       call print_CijR_QUAD(C0156_QUAD,C156R_QUAD,C156I_QUAD)
       call print_CijR_QUAD_D(C0156,C156R,C156I,C0156_QUAD,C156R_QUAD,
     -   C156I_QUAD)
       call print_CijR(C0156,C156R,C156I)
c       call tens_Int3_Re_Com(0d0,s56,p5sq,p6sq,B056,B016,B015,C0156R,C
c     -   0156I,C156R,C156I)
      Endif    
       If(PrintD) then
       write(*,*)"D1234"
       call print_DijR_QUAD(D01234_QUAD,D1234R_QUAD,D1234I_QUAD)
       call print_DijR_QUAD_D(D01234,D1234R,D1234I,D01234_QUAD,D1234R_
     -   QUAD,D1234I_QUAD)
       call print_DijR(D01234,D1234R,D1234I)
cfc       call tens_Int4_Re_Com(0d0,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
cfc     -   C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C0134I,C
cfc     -   0124I,C0123I,C234I,C134I,C124I,C123I,D01234R,D01234I,D1234R,
cfc     -   D1234I)
       write(*,*)"D1235"
       call print_DijR_QUAD(D01235_QUAD,D1235R_QUAD,D1235I_QUAD)
       call print_DijR_QUAD_D(D01235,D1235R,D1235I,D01235_QUAD,D1235R_
     -   QUAD,D1235I_QUAD)
       call print_DijR(D01235,D1235R,D1235I)
cfc       call tens_Int4_Re_Com(0d0,p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
cfc     -   4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0235I
cfc     -   ,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235R,D01235
cfc     -   I,D1235R,D1235I)
       write(*,*)"D1236"
       call print_DijR_QUAD(D01236_QUAD,D1236R_QUAD,D1236I_QUAD)
       call print_DijR_QUAD_D(D01236,D1236R,D1236I,D01236_QUAD,D1236R_
     -   QUAD,D1236I_QUAD)
       call print_DijR(D01236,D1236R,D1236I)
cfc       call tens_Int4_Re_Com(0d0,p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p5,p2
cfc     -   p3+p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C126R,C
cfc     -   123R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123I,D01
cfc     -   236R,D01236I,D1236R,D1236I)
       write(*,*)"D1245"
       call print_DijR_QUAD(D01245_QUAD,D1245R_QUAD,D1245I_QUAD)
       call print_DijR_QUAD_D(D01245,D1245R,D1245I,D01245_QUAD,D1245R_
     -   QUAD,D1245I_QUAD)
       call print_DijR(D01245,D1245R,D1245I)
cfc       call tens_Int4_Re_Com(0d0,p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
cfc     -   4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0245I
cfc     -   ,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245R,D01245
cfc     -   I,D1245R,D1245I)
       write(*,*)"D1246"
       call print_DijR_QUAD(D01246_QUAD,D1246R_QUAD,D1246I_QUAD)
       call print_DijR_QUAD_D(D01246,D1246R,D1246I,D01246_QUAD,D1246R_
     -   QUAD,D1246I_QUAD)
       call print_DijR(D01246,D1246R,D1246I)
cfc       call tens_Int4_Re_Com(0d0,p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,p2p4
cfc     -   +p2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R,C126
cfc     -   R,C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C124I,
cfc     -   D01246R,D01246I,D1246R,D1246I)
       write(*,*)"D1256"
       call print_DijR_QUAD(D01256_QUAD,D1256R_QUAD,D1256I_QUAD)
       call print_DijR_QUAD_D(D01256,D1256R,D1256I,D01256_QUAD,D1256R_
     -   QUAD,D1256I_QUAD)
       call print_DijR(D01256,D1256R,D1256I)
cfc       call tens_Int4_Re_Com(0d0,p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p5,p2
cfc     -   p5+p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C126R,C
cfc     -   125R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125I,D01
cfc     -   256R,D01256I,D1256R,D1256I)
       write(*,*)"D2345"
       call print_DijR_QUAD(D02345_QUAD,D2345R_QUAD,D2345I_QUAD)
       call print_DijR_QUAD_D(D02345,D2345R,D2345I,D02345_QUAD,D2345R_
     -   QUAD,D2345I_QUAD)
       call print_DijR(D02345,D2345R,D2345I)
cfc       call tens_Int4_Re_Com(0d0,p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
cfc     -   C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C0245I,C
cfc     -   0235I,C0234I,C345I,C245I,C235I,C234I,D02345R,D02345I,D2345R,
cfc     -   D2345I)
       write(*,*)"D2346"
       call print_DijR_QUAD(D02346_QUAD,D2346R_QUAD,D2346I_QUAD)
       call print_DijR_QUAD_D(D02346,D2346R,D2346I,D02346_QUAD,D2346R_
     -   QUAD,D2346I_QUAD)
       call print_DijR(D02346,D2346R,D2346I)
cfc       call tens_Int4_Re_Com(0d0,p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4+p3p
cfc     -   5,C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C0346I
cfc     -   ,C0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346R,D02346
cfc     -   I,D2346R,D2346I)
       write(*,*)"D2356"
       call print_DijR_QUAD(D02356_QUAD,D2356R_QUAD,D2356I_QUAD)
       call print_DijR_QUAD_D(D02356,D2356R,D2356I,D02356_QUAD,D2356R_
     -   QUAD,D2356I_QUAD)
       call print_DijR(D02356,D2356R,D2356I)
cfc       call tens_Int4_Re_Com(0d0,p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5+p4p
cfc     -   5,C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C0356I
cfc     -   ,C0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356R,D02356
cfc     -   I,D2356R,D2356I)
       write(*,*)"D3456"
       call print_DijR_QUAD(D03456_QUAD,D3456R_QUAD,D3456I_QUAD)
       call print_DijR_QUAD_D(D03456,D3456R,D3456I,D03456_QUAD,D3456R_
     -   QUAD,D3456I_QUAD)
       call print_DijR(D03456,D3456R,D3456I)
cfc       call tens_Int4_Re_Com(0d0,p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C0456R,
cfc     -   C0356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C0356I,C
cfc     -   0346I,C0345I,C456I,C356I,C346I,C345I,D03456R,D03456I,D3456R,
cfc     -   D3456I)
       write(*,*)"D1345"
       call print_DijR_QUAD(D01345_QUAD,D1345R_QUAD,D1345I_QUAD)
       call print_DijR_QUAD_D(D01345,D1345R,D1345I,D01345_QUAD,D1345R_
     -   QUAD,D1345I_QUAD)
       call print_DijR(D01345,D1345R,D1345I)
cfc       call tens_Int4_Re_Com(0d0,s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
cfc     -   4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0345I
cfc     -   ,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345R,D01345
cfc     -   I,D1345R,D1345I)
       write(*,*)"D1346"
       call print_DijR_QUAD(D01346_QUAD,D1346R_QUAD,D1346I_QUAD)
       call print_DijR_QUAD_D(D01346,D1346R,D1346I,D01346_QUAD,D1346R_
     -   QUAD,D1346I_QUAD)
       call print_DijR(D01346,D1346R,D1346I)
cfc       call tens_Int4_Re_Com(0d0,s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+p2p4
cfc     -   +p2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R,C136
cfc     -   R,C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C134I,
cfc     -   D01346R,D01346I,D1346R,D1346I)
       write(*,*)"D1356"
       call print_DijR_QUAD(D01356_QUAD,D1356R_QUAD,D1356I_QUAD)
       call print_DijR_QUAD_D(D01356,D1356R,D1356I,D01356_QUAD,D1356R_
     -   QUAD,D1356I_QUAD)
       call print_DijR(D01356,D1356R,D1356I)
cfc       call tens_Int4_Re_Com(0d0,s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,p1p5
cfc     -   +p2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R,C136
cfc     -   R,C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C135I,
cfc     -   D01356R,D01356I,D1356R,D1356I)
       write(*,*)"D1456"
       call print_DijR_QUAD(D01456_QUAD,D1456R_QUAD,D1456I_QUAD)
       call print_DijR_QUAD_D(D01456,D1456R,D1456I,D01456_QUAD,D1456R_
     -   QUAD,D1456I_QUAD)
       call print_DijR(D01456,D1456R,D1456I)
cfc       call tens_Int4_Re_Com(0d0,s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p5+p2
cfc     -   p5+p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C146R,C
cfc     -   145R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145I,D01
cfc     -   456R,D01456I,D1456R,D1456I)
       write(*,*)"D2456"
       call print_DijR_QUAD(D02456_QUAD,D2456R_QUAD,D2456I_QUAD)
       call print_DijR_QUAD_D(D02456,D2456R,D2456I,D02456_QUAD,D2456R_
     -   QUAD,D2456I_QUAD)
       call print_DijR(D02456,D2456R,D2456I)
cfc       call tens_Int4_Re_Com(0d0,s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5,p4p
cfc     -   5,C0456R,C0256R,C0246R,C0245R,C456R,C256R,C246R,C245R,C0456I
cfc     -   ,C0256I,C0246I,C0245I,C456I,C256I,C246I,C245I,D02456R,D02456
cfc     -   I,D2456R,D2456I)
      Endif    
       If(PrintE) then
       write(*,*)"E12345"
       call print_EijR_QUAD(E012345_QUAD,E12345R_QUAD,E12345I_QUAD)
       call print_EijR_QUAD_D(E012345,E12345R,E12345I,E012345_QUAD,E12345R_
     -   QUAD,E12345I_QUAD)
       call print_EijR(E012345,E12345R,E12345I)
cfc       call tens_Int5_Re_Com(0d0,p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
cfc     -   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2345R,
cfc     -   D1345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01235I,
cfc     -   D01234I,D2345I,D1345I,D1245I,D1235I,D1234I,E012345R,E012345I
cfc     -   ,E12345R,E12345I)
       write(*,*)"E12346"
       call print_EijR_QUAD(E012346_QUAD,E12346R_QUAD,E12346I_QUAD)
       call print_EijR_QUAD_D(E012346,E12346R,E12346I,E012346_QUAD,E12346R_
     -   QUAD,E12346I_QUAD)
       call print_EijR(E012346,E12346R,E12346I)
cfc       call tens_Int5_Re_Com(0d0,p1sq,p2sq,p3sq,s45,p1p2,p1p3,p1p4+p1p
cfc     -   5,p2p3,p2p4+p2p5,p3p4+p3p5,D02346R,D01346R,D01246R,D01236R,D
cfc     -   01234R,D2346R,D1346R,D1246R,D1236R,D1234R,D02346I,D01346I,D0
cfc     -   1246I,D01236I,D01234I,D2346I,D1346I,D1246I,D1236I,D1234I,E01
cfc     -   2346R,E012346I,E12346R,E12346I)
       write(*,*)"E12356"
       call print_EijR_QUAD(E012356_QUAD,E12356R_QUAD,E12356I_QUAD)
       call print_EijR_QUAD_D(E012356,E12356R,E12356I,E012356_QUAD,E12356R_
     -   QUAD,E12356I_QUAD)
       call print_EijR(E012356,E12356R,E12356I)
cfc       call tens_Int5_Re_Com(0d0,p1sq,p2sq,s34,p5sq,p1p2,p1p3+p1p4,p1p
cfc     -   5,p2p3+p2p4,p2p5,p3p5+p4p5,D02356R,D01356R,D01256R,D01236R,D
cfc     -   01235R,D2356R,D1356R,D1256R,D1236R,D1235R,D02356I,D01356I,D0
cfc     -   1256I,D01236I,D01235I,D2356I,D1356I,D1256I,D1236I,D1235I,E01
cfc     -   2356R,E012356I,E12356R,E12356I)
       write(*,*)"E12456"
       call print_EijR_QUAD(E012456_QUAD,E12456R_QUAD,E12456I_QUAD)
       call print_EijR_QUAD_D(E012456,E12456R,E12456I,E012456_QUAD,E12456R_
     -   QUAD,E12456I_QUAD)
       call print_EijR(E012456,E12456R,E12456I)
cfc       call tens_Int5_Re_Com(0d0,p1sq,s23,p4sq,p5sq,p1p2+p1p3,p1p4,p1p
cfc     -   5,p2p4+p3p4,p2p5+p3p5,p4p5,D02456R,D01456R,D01256R,D01246R,D
cfc     -   01245R,D2456R,D1456R,D1256R,D1246R,D1245R,D02456I,D01456I,D0
cfc     -   1256I,D01246I,D01245I,D2456I,D1456I,D1256I,D1246I,D1245I,E01
cfc     -   2456R,E012456I,E12456R,E12456I)
       write(*,*)"E23456"
       call print_EijR_QUAD(E023456_QUAD,E23456R_QUAD,E23456I_QUAD)
       call print_EijR_QUAD_D(E023456,E23456R,E23456I,E023456_QUAD,E23456R_
     -   QUAD,E23456I_QUAD)
       call print_EijR(E023456,E23456R,E23456I)
cfc       call tens_Int5_Re_Com(0d0,p2sq,p3sq,p4sq,p5sq,p2p3,p2p4,p2p5,p3
cfc     -   p4,p3p5,p4p5,D03456R,D02456R,D02356R,D02346R,D02345R,D3456R,
cfc     -   D2456R,D2356R,D2346R,D2345R,D03456I,D02456I,D02356I,D02346I,
cfc     -   D02345I,D3456I,D2456I,D2356I,D2346I,D2345I,E023456R,E023456I
cfc     -   ,E23456R,E23456I)
       write(*,*)"E13456"
       call print_EijR_QUAD(E013456_QUAD,E13456R_QUAD,E13456I_QUAD)
       call print_EijR_QUAD_D(E013456,E13456R,E13456I,E013456_QUAD,E13456R_
     -   QUAD,E13456I_QUAD)
       call print_EijR(E013456,E13456R,E13456I)
cfc       call tens_Int5_Re_Com(0d0,s12,p3sq,p4sq,p5sq,p1p3+p2p3,p1p4+p2p
cfc     -   4,p1p5+p2p5,p3p4,p3p5,p4p5,D03456R,D01456R,D01356R,D01346R,D
cfc     -   01345R,D3456R,D1456R,D1356R,D1346R,D1345R,D03456I,D01456I,D0
cfc     -   1356I,D01346I,D01345I,D3456I,D1456I,D1356I,D1346I,D1345I,E01
cfc     -   3456R,E013456I,E13456R,E13456I)
      Endif    
       If(PrintF) then
       write(*,*)"F123456"
       call print_FijR_QUAD(F0123456_QUAD,F123456R_QUAD,F123456I_QUAD)
       call print_FijR_QUAD_D(F0123456,F123456R,F123456I,F0123456_QUAD
     -   ,F123456R_QUAD,F123456I_QUAD)
       call print_FijR(F0123456,F123456R,F123456I)
      Endif  


      end
