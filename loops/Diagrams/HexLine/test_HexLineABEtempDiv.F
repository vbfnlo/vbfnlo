       subroutine test_hexlineABEDiv(p1,p2,p3,p4,p5,p6,barpsi_p6,psi_p1,m
     -   up2,mup3,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result,re
     -   sultborn,Div)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 27/10/2008
c Modified:9/9/2009
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
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16   p1mup2, p1mup3, p1mup4, p1mup5, p2mup2, p2mup3, 
     -          p2mup4, p2mup5, p3mup2, p3mup3, p3mup4, p3mup5, 
     -          p4mup2, p4mup3, p4mup4, p4mup5, p5mup2, p5mup3, 
     -          p5mup4, p5mup5, p6mup2, p6mup3, p6mup4, p6mup5
       Complex*16   mup2mup3, mup2mup4, mup2mup5, mup3mup4, mup3mup5, 
     -          mup4mup5
       Real*8 dotrr
       Complex*16 B0finDiv,C0finDiv,D0finDiv,E0fin,F0fin
       EXTERNAL dotrr,B0finDiv,C0finDiv,D0finDiv,E0fin,F0fin
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
       Integer comp,Div 
      Complex*16 SMB(64),F(5930),Fa(445)
      Real*8 P(11674)
      Complex*16 dotrc,dotcc,result(13),resultborn
       Real*8 musq
      EXTERNAL   dotrc,dotcc
      Integer alpha
       COMMON/test_FaFunctionsAbe/Fa
       COMMON/test_FhlFunctionsABE/F
      Save/test_FhlFunctionsABE/
       COMMON/test_PFunctionsABE/P
       COMMON/Invariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45
     -   ,s56,s16,s123,s234,s345
       COMMON/FVALFunctions/F0123456R,F123456R,F0123456I,F123456I
       COMMON/EVALFunctions/ E012345R,E12345R,E012345I,E12345I, E01234
     -   6R,E12346R,E012346I,E12346I, E012356R,E12356R,E012356I,E1235
     -   6I, E012456R,E12456R,E012456I,E12456I, E013456R,E13456R,E013
     -   456I,E13456I, E023456R,E23456R,E023456I,E23456I
       COMMON/DVALFunctions/ D01234R,D1234R,D01234I,D1234I, D01235R,D1
     -   235R,D01235I,D1235I, D01236R,D1236R,D01236I,D1236I, D01245R,
     -   D1245R,D01245I,D1245I, D01246R,D1246R,D01246I,D1246I, D01256
     -   R,D1256R,D01256I,D1256I, D01345R,D1345R,D01345I,D1345I, D013
     -   46R,D1346R,D01346I,D1346I, D01356R,D1356R,D01356I,D1356I, D0
     -   1456R,D1456R,D01456I,D1456I, D02345R,D2345R,D02345I,D2345I, 
     -   D02346R,D2346R,D02346I,D2346I, D02356R,D2356R,D02356I,D2356I
     -   , D02456R,D2456R,D02456I,D2456I, D03456R,D3456R,D03456I,D345
     -   6I
       COMMON/CVALFunctions/ C0123R,C123R,C0123I,C123I, C0124R,C124R,C
     -   0124I,C124I, C0125R,C125R,C0125I,C125I, C0126R,C126R,C0126I,
     -   C126I, C0134R,C134R,C0134I,C134I, C0135R,C135R,C0135I,C135I,
     -    C0136R,C136R,C0136I,C136I, C0145R,C145R,C0145I,C145I, C0146
     -   R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C156I, C0234R,C234
     -   R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0236R,C236R,C023
     -   6I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C246R,C0246I,C24
     -   6I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C0345I,C345I, C0
     -   346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,C356I, C0456R,C
     -   456R,C0456I,C456I
       COMMON/BVALFunctions/ B012R,B012I, B013R,B013I, B014R,B014I, B0
     -   15R,B015I, B016R,B016I, B023R,B023I, B024R,B024I, B025R,B025
     -   I, B026R,B026I, B034R,B034I, B035R,B035I, B036R,B036I, B045R
     -   ,B045I, B046R,B046I, B056R,B056I
       Complex*16 SC1c,SC1r,SC3ccc,SC3ccr,SC3crr,SC3rrr,SC5ccccr,SC5cc
     -   crr,SC5ccrrr,SC7ccccrrr
       External SC1c,SC1r,SC3ccc,SC3ccr,SC3crr,SC3rrr,SC5ccccr,SC5cccr
     -   r,SC5ccrrr,SC7ccccrrr
       Complex*16 barpsi_p6(2),psi_p1(2)
       Integer ngluon, posgluon
       Real*8 musqIn
       Common/musqInv/musq
       musq=musqIn
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
       If(PrintB) then
       write(*,*)"B12"
       write(*,*) "B12=", B012
       endif

       B023=B0finDiv(p2sq,musq,Div)
       B023R=Dble(B023)
       B023I=DImag(B023)
       If(PrintB) then
       write(*,*)"B23"
       write(*,*) "B23=", B023
       endif

       B034=B0finDiv(p3sq,musq,Div)
       B034R=Dble(B034)
       B034I=DImag(B034)
       If(PrintB) then
       write(*,*)"B34"
       write(*,*) "B34=", B034
       endif

       B045=B0finDiv(p4sq,musq,Div)
       B045R=Dble(B045)
       B045I=DImag(B045)
       If(PrintB) then
       write(*,*)"B45"
       write(*,*) "B45=", B045
       endif

       B056=B0finDiv(p5sq,musq,Div)
       B056R=Dble(B056)
       B056I=DImag(B056)
       If(PrintB) then
       write(*,*)"B56"
       write(*,*) "B56=", B056
       endif

       B016=B0finDiv(p6sq,musq,Div)
       B016R=Dble(B016)
       B016I=DImag(B016)
       If(PrintB) then
       write(*,*)"B16"
       write(*,*) "B16=", B016
       endif

       B013=B0finDiv(s12,musq,Div)
       B013R=Dble(B013)
       B013I=DImag(B013)
       If(PrintB) then
       write(*,*)"B13"
       write(*,*) "B13=", B013
       endif

       B014=B0finDiv(s123,musq,Div)
       B014R=Dble(B014)
       B014I=DImag(B014)
       If(PrintB) then
       write(*,*)"B14"
       write(*,*) "B14=", B014
       endif

       B026=B0finDiv(s16,musq,Div)
       B026R=Dble(B026)
       B026I=DImag(B026)
       If(PrintB) then
       write(*,*)"B26"
       write(*,*) "B26=", B026
       endif

       B024=B0finDiv(s23,musq,Div)
       B024R=Dble(B024)
       B024I=DImag(B024)
       If(PrintB) then
       write(*,*)"B24"
       write(*,*) "B24=", B024
       endif

       B025=B0finDiv(s234,musq,Div)
       B025R=Dble(B025)
       B025I=DImag(B025)
       If(PrintB) then
       write(*,*)"B25"
       write(*,*) "B25=", B025
       endif

       B035=B0finDiv(s34,musq,Div)
       B035R=Dble(B035)
       B035I=DImag(B035)
       If(PrintB) then
       write(*,*)"B35"
       write(*,*) "B35=", B035
       endif

       B036=B0finDiv(s345,musq,Div)
       B036R=Dble(B036)
       B036I=DImag(B036)
       If(PrintB) then
       write(*,*)"B36"
       write(*,*) "B36=", B036
       endif

       B046=B0finDiv(s45,musq,Div)
       B046R=Dble(B046)
       B046I=DImag(B046)
       If(PrintB) then
       write(*,*)"B46"
       write(*,*) "B46=", B046
       endif

       B015=B0finDiv(s56,musq,Div)
       B015R=Dble(B015)
       B015I=DImag(B015)
       If(PrintB) then
       write(*,*)"B15"
       write(*,*) "B15=", B015
       endif

       C0123=C0finDiv(p1sq,p2sq,s12,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,p2sq,s12,B023,B013,B012,C0123,C0
     -   123R,C0123I,C123R,C123I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C123"
c$$$       call print_CijR(C0123,C123R,C123I)
c$$$       call tens_Int3_Re_Com(0d0,p1sq,p2sq,s12,B023,B013,B012,C0123R,C
c$$$     -   0123I,C123R,C123I)
c$$$       Endif

       C0126=C0finDiv(p1sq,s16,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,s16,p6sq,B026,B016,B012,C0126,C0
     -   126R,C0126I,C126R,C126I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C126"
c$$$       call print_CijR(C0126,C126R,C126I)
c$$$       call tens_Int3_Re_Com(0d0,p1sq,s16,p6sq,B026,B016,B012,C0126R,C
c$$$     -   0126I,C126R,C126I)
c$$$       Endif

       C0124=C0finDiv(p1sq,s23,s123,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,s23,s123,B024,B014,B012,C0124,C0
     -   124R,C0124I,C124R,C124I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C124"
c$$$       call print_CijR(C0124,C124R,C124I)
c$$$       call tens_Int3_Re_Com(0d0,p1sq,s23,s123,B024,B014,B012,C0124R,C
c$$$     -   0124I,C124R,C124I)
c$$$       Endif

       C0125=C0finDiv(p1sq,s234,s56,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,s234,s56,B025,B015,B012,C0125,C0
     -   125R,C0125I,C125R,C125I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C125"
c$$$       call print_CijR(C0125,C125R,C125I)
c$$$       call tens_Int3_Re_Com(0d0,p1sq,s234,s56,B025,B015,B012,C0125R,C
c$$$     -   0125I,C125R,C125I)
c$$$       Endif

       C0234=C0finDiv(p2sq,p3sq,s23,musq,Div)
       call tens_red3_new_Re_Com_Div(p2sq,p3sq,s23,B034,B024,B023,C0234,C0
     -   234R,C0234I,C234R,C234I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C234"
c$$$       call print_CijR(C0234,C234R,C234I)
c$$$       call tens_Int3_Re_Com(0d0,p2sq,p3sq,s23,B034,B024,B023,C0234R,C
c$$$     -   0234I,C234R,C234I)
c$$$       Endif

       C0235=C0finDiv(p2sq,s34,s234,musq,Div)
       call tens_red3_new_Re_Com_Div(p2sq,s34,s234,B035,B025,B023,C0235,C0
     -   235R,C0235I,C235R,C235I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C235"
c$$$       call print_CijR(C0235,C235R,C235I)
c$$$       call tens_Int3_Re_Com(0d0,p2sq,s34,s234,B035,B025,B023,C0235R,C
c$$$     -   0235I,C235R,C235I)
c$$$       Endif

       C0236=C0finDiv(p2sq,s345,s16,musq,Div)
       call tens_red3_new_Re_Com_Div(p2sq,s345,s16,B036,B026,B023,C0236,C0
     -   236R,C0236I,C236R,C236I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C236"
c$$$       call print_CijR(C0236,C236R,C236I)
c$$$       call tens_Int3_Re_Com(0d0,p2sq,s345,s16,B036,B026,B023,C0236R,C
c$$$     -   0236I,C236R,C236I)
c$$$       Endif

       C0345=C0finDiv(p3sq,p4sq,s34,musq,Div)
       call tens_red3_new_Re_Com_Div(p3sq,p4sq,s34,B045,B035,B034,C0345,C0
     -   345R,C0345I,C345R,C345I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C345"
c$$$       call print_CijR(C0345,C345R,C345I)
c$$$       call tens_Int3_Re_Com(0d0,p3sq,p4sq,s34,B045,B035,B034,C0345R,C
c$$$     -   0345I,C345R,C345I)
c$$$       Endif

       C0346=C0finDiv(p3sq,s45,s345,musq,Div)
       call tens_red3_new_Re_Com_Div(p3sq,s45,s345,B046,B036,B034,C0346,C0
     -   346R,C0346I,C346R,C346I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C346"
c$$$       call print_CijR(C0346,C346R,C346I)
c$$$       call tens_Int3_Re_Com(0d0,p3sq,s45,s345,B046,B036,B034,C0346R,C
c$$$     -   0346I,C346R,C346I)
c$$$       Endif

       C0456=C0finDiv(p4sq,p5sq,s45,musq,Div)
       call tens_red3_new_Re_Com_Div(p4sq,p5sq,s45,B056,B046,B045,C0456,C0
     -   456R,C0456I,C456R,C456I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C456"
c$$$       call print_CijR(C0456,C456R,C456I)
c$$$       call tens_Int3_Re_Com(0d0,p4sq,p5sq,s45,B056,B046,B045,C0456R,C
c$$$     -   0456I,C456R,C456I)
c$$$       Endif

       C0134=C0finDiv(s12,p3sq,s123,musq,Div)
       call tens_red3_new_Re_Com_Div(s12,p3sq,s123,B034,B014,B013,C0134,C0
     -   134R,C0134I,C134R,C134I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C134"
c$$$       call print_CijR(C0134,C134R,C134I)
c$$$       call tens_Int3_Re_Com(0d0,s12,p3sq,s123,B034,B014,B013,C0134R,C
c$$$     -   0134I,C134R,C134I)
c$$$       Endif

       C0135=C0finDiv(s12,s34,s56,musq,Div)
       call tens_red3_new_Re_Com_Div(s12,s34,s56,B035,B015,B013,C0135,C013
     -   5R,C0135I,C135R,C135I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C135"
c$$$       call print_CijR(C0135,C135R,C135I)
c$$$       call tens_Int3_Re_Com(0d0,s12,s34,s56,B035,B015,B013,C0135R,C01
c$$$     -   35I,C135R,C135I)
c$$$       Endif

       C0136=C0finDiv(s12,s345,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(s12,s345,p6sq,B036,B016,B013,C0136,C0
     -   136R,C0136I,C136R,C136I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C136"
c$$$       call print_CijR(C0136,C136R,C136I)
c$$$       call tens_Int3_Re_Com(0d0,s12,s345,p6sq,B036,B016,B013,C0136R,C
c$$$     -   0136I,C136R,C136I)
c$$$       Endif

       C0145=C0finDiv(s123,p4sq,s56,musq,Div)
       call tens_red3_new_Re_Com_Div(s123,p4sq,s56,B045,B015,B014,C0145,C0
     -   145R,C0145I,C145R,C145I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C145"
c$$$       call print_CijR(C0145,C145R,C145I)
c$$$       call tens_Int3_Re_Com(0d0,s123,p4sq,s56,B045,B015,B014,C0145R,C
c$$$     -   0145I,C145R,C145I)
c$$$       Endif

       C0146=C0finDiv(s123,s45,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(s123,s45,p6sq,B046,B016,B014,C0146,C0
     -   146R,C0146I,C146R,C146I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C146"
c$$$       call print_CijR(C0146,C146R,C146I)
c$$$       call tens_Int3_Re_Com(0d0,s123,s45,p6sq,B046,B016,B014,C0146R,C
c$$$     -   0146I,C146R,C146I)
c$$$       Endif

       C0245=C0finDiv(s23,p4sq,s234,musq,Div)
       call tens_red3_new_Re_Com_Div(s23,p4sq,s234,B045,B025,B024,C0245,C0
     -   245R,C0245I,C245R,C245I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C245"
c$$$       call print_CijR(C0245,C245R,C245I)
c$$$       call tens_Int3_Re_Com(0d0,s23,p4sq,s234,B045,B025,B024,C0245R,C
c$$$     -   0245I,C245R,C245I)
c$$$       Endif

       C0246=C0finDiv(s23,s45,s16,musq,Div)
       call tens_red3_new_Re_Com_Div(s23,s45,s16,B046,B026,B024,C0246,C024
     -   6R,C0246I,C246R,C246I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C246"
c$$$       call print_CijR(C0246,C246R,C246I)
c$$$       call tens_Int3_Re_Com(0d0,s23,s45,s16,B046,B026,B024,C0246R,C02
c$$$     -   46I,C246R,C246I)
c$$$       Endif

       C0256=C0finDiv(s234,p5sq,s16,musq,Div)
       call tens_red3_new_Re_Com_Div(s234,p5sq,s16,B056,B026,B025,C0256,C0
     -   256R,C0256I,C256R,C256I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C256"
c$$$       call print_CijR(C0256,C256R,C256I)
c$$$       call tens_Int3_Re_Com(0d0,s234,p5sq,s16,B056,B026,B025,C0256R,C
c$$$     -   0256I,C256R,C256I)
c$$$       Endif

       C0356=C0finDiv(s34,p5sq,s345,musq,Div)
       call tens_red3_new_Re_Com_Div(s34,p5sq,s345,B056,B036,B035,C0356,C0
     -   356R,C0356I,C356R,C356I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C356"
c$$$       call print_CijR(C0356,C356R,C356I)
c$$$       call tens_Int3_Re_Com(0d0,s34,p5sq,s345,B056,B036,B035,C0356R,C
c$$$     -   0356I,C356R,C356I)
c$$$       Endif

       C0156=C0finDiv(s56,p5sq,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(s56,p5sq,p6sq,B056,B016,B015,C0156,C0
     -   156R,C0156I,C156R,C156I)
c$$$       If(PrintC) then
c$$$       write(*,*)"C156"
c$$$       call print_CijR(C0156,C156R,C156I)
c$$$       call tens_Int3_Re_Com(0d0,s56,p5sq,p6sq,B056,B016,B015,C0156R,C
c$$$     -   0156I,C156R,C156I)
c$$$       Endif

       D01234=D0finDiv(s12,s23,p1sq,p2sq,p3sq,s123,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
     -   C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C0134I,C
     -   0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D01234I,
     -   D1234R,D1234I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D1234"
c$$$       call print_DijR(D01234,D1234R,D1234I)
c$$$       call tens_Int4_Re_Com(0d0,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
c$$$     -   C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C0134I,C
c$$$     -   0124I,C0123I,C234I,C134I,C124I,C123I,D01234R,D01234I,D1234R,
c$$$     -   D1234I)
c$$$       endif

       D01235=D0finDiv(s12,s234,p1sq,p2sq,s34,s56,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
     -   4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0235I
     -   ,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235,D01235R
     -   ,D01235I,D1235R,D1235I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D1235"
c$$$       call print_DijR(D01235,D1235R,D1235I)
c$$$       call tens_Int4_Re_Com(0d0,p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
c$$$     -   4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0235I
c$$$     -   ,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235R,D01235
c$$$     -   I,D1235R,D1235I)
c$$$       endif

       D01236=D0finDiv(s12,s16,p1sq,p2sq,s345,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p5,p2
     -   p3+p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C126R,C
     -   123R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123I,D01
     -   236,D01236R,D01236I,D1236R,D1236I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D1236"
c$$$       call print_DijR(D01236,D1236R,D1236I)
c$$$       call tens_Int4_Re_Com(0d0,p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p5,p2
c$$$     -   p3+p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C126R,C
c$$$     -   123R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123I,D01
c$$$     -   236R,D01236I,D1236R,D1236I)
c$$$       endif

       D01245=D0finDiv(s123,s234,p1sq,s23,p4sq,s56,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
     -   4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0245I
     -   ,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245,D01245R
     -   ,D01245I,D1245R,D1245I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D1245"
c$$$       call print_DijR(D01245,D1245R,D1245I)
c$$$       call tens_Int4_Re_Com(0d0,p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
c$$$     -   4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0245I
c$$$     -   ,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245R,D01245
c$$$     -   I,D1245R,D1245I)
c$$$       endif

       D01246=D0finDiv(s123,s16,p1sq,s23,s45,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,p2p4
     -   +p2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R,C126
     -   R,C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C124I,
     -   D01246,D01246R,D01246I,D1246R,D1246I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D1246"
c$$$       call print_DijR(D01246,D1246R,D1246I)
c$$$       call tens_Int4_Re_Com(0d0,p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,p2p4
c$$$     -   +p2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R,C126
c$$$     -   R,C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C124I,
c$$$     -   D01246R,D01246I,D1246R,D1246I)
c$$$       endif

       D01256=D0finDiv(s56,s16,p1sq,s234,p5sq,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p5,p2
     -   p5+p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C126R,C
     -   125R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125I,D01
     -   256,D01256R,D01256I,D1256R,D1256I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D1256"
c$$$       call print_DijR(D01256,D1256R,D1256I)
c$$$       call tens_Int4_Re_Com(0d0,p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p5,p2
c$$$     -   p5+p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C126R,C
c$$$     -   125R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125I,D01
c$$$     -   256R,D01256I,D1256R,D1256I)
c$$$       endif

       D02345=D0finDiv(s23,s34,p2sq,p3sq,p4sq,s234,musq,Div)
       call tens_red4_new_Re_Com_Div(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
     -   C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C0245I,C
     -   0235I,C0234I,C345I,C245I,C235I,C234I,D02345,D02345R,D02345I,
     -   D2345R,D2345I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D2345"
c$$$       call print_DijR(D02345,D2345R,D2345I)
c$$$       call tens_Int4_Re_Com(0d0,p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
c$$$     -   C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C0245I,C
c$$$     -   0235I,C0234I,C345I,C245I,C235I,C234I,D02345R,D02345I,D2345R,
c$$$     -   D2345I)
c$$$       endif

       D02346=D0finDiv(s23,s345,p2sq,p3sq,s45,s16,musq,Div)
       call tens_red4_new_Re_Com_Div(p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4+p3p
     -   5,C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C0346I
     -   ,C0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346,D02346R
     -   ,D02346I,D2346R,D2346I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D2346"
c$$$       call print_DijR(D02346,D2346R,D2346I)
c$$$       call tens_Int4_Re_Com(0d0,p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4+p3p
c$$$     -   5,C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C0346I
c$$$     -   ,C0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346R,D02346
c$$$     -   I,D2346R,D2346I)
c$$$       endif

       D02356=D0finDiv(s234,s345,p2sq,s34,p5sq,s16,musq,Div)
       call tens_red4_new_Re_Com_Div(p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5+p4p
     -   5,C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C0356I
     -   ,C0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356,D02356R
     -   ,D02356I,D2356R,D2356I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D2356"
c$$$       call print_DijR(D02356,D2356R,D2356I)
c$$$       call tens_Int4_Re_Com(0d0,p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5+p4p
c$$$     -   5,C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C0356I
c$$$     -   ,C0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356R,D02356
c$$$     -   I,D2356R,D2356I)
c$$$       endif

       D03456=D0finDiv(s34,s45,p3sq,p4sq,p5sq,s345,musq,Div)
       call tens_red4_new_Re_Com_Div(p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C0456R,
     -   C0356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C0356I,C
     -   0346I,C0345I,C456I,C356I,C346I,C345I,D03456,D03456R,D03456I,
     -   D3456R,D3456I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D3456"
c$$$       call print_DijR(D03456,D3456R,D3456I)
c$$$       call tens_Int4_Re_Com(0d0,p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C0456R,
c$$$     -   C0356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C0356I,C
c$$$     -   0346I,C0345I,C456I,C356I,C346I,C345I,D03456R,D03456I,D3456R,
c$$$     -   D3456I)
c$$$       endif

       D01345=D0finDiv(s123,s34,s12,p3sq,p4sq,s56,musq,Div)
       call tens_red4_new_Re_Com_Div(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
     -   4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0345I
     -   ,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345,D01345R
     -   ,D01345I,D1345R,D1345I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D1345"
c$$$       call print_DijR(D01345,D1345R,D1345I)
c$$$       call tens_Int4_Re_Com(0d0,s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
c$$$     -   4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0345I
c$$$     -   ,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345R,D01345
c$$$     -   I,D1345R,D1345I)
c$$$       endif

       D01346=D0finDiv(s123,s345,s12,p3sq,s45,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+p2p4
     -   +p2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R,C136
     -   R,C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C134I,
     -   D01346,D01346R,D01346I,D1346R,D1346I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D1346"
c$$$       call print_DijR(D01346,D1346R,D1346I)
c$$$       call tens_Int4_Re_Com(0d0,s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+p2p4
c$$$     -   +p2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R,C136
c$$$     -   R,C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C134I,
c$$$     -   D01346R,D01346I,D1346R,D1346I)
c$$$       endif

       D01356=D0finDiv(s56,s345,s12,s34,p5sq,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,p1p5
     -   +p2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R,C136
     -   R,C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C135I,
     -   D01356,D01356R,D01356I,D1356R,D1356I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D1356"
c$$$       call print_DijR(D01356,D1356R,D1356I)
c$$$       call tens_Int4_Re_Com(0d0,s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,p1p5
c$$$     -   +p2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R,C136
c$$$     -   R,C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C135I,
c$$$     -   D01356R,D01356I,D1356R,D1356I)
c$$$       endif

       D01456=D0finDiv(s56,s45,s123,p4sq,p5sq,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p5+p2
     -   p5+p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C146R,C
     -   145R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145I,D01
     -   456,D01456R,D01456I,D1456R,D1456I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D1456"
c$$$       call print_DijR(D01456,D1456R,D1456I)
c$$$       call tens_Int4_Re_Com(0d0,s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p5+p2
c$$$     -   p5+p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C146R,C
c$$$     -   145R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145I,D01
c$$$     -   456R,D01456I,D1456R,D1456I)
c$$$       endif

       D02456=D0finDiv(s234,s45,s23,p4sq,p5sq,s16,musq,Div)
       call tens_red4_new_Re_Com_Div(s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5,p4p
     -   5,C0456R,C0256R,C0246R,C0245R,C456R,C256R,C246R,C245R,C0456I
     -   ,C0256I,C0246I,C0245I,C456I,C256I,C246I,C245I,D02456,D02456R
     -   ,D02456I,D2456R,D2456I)
c$$$       If(PrintD) then
c$$$       write(*,*)"D2456"
c$$$       call print_DijR(D02456,D2456R,D2456I)
c$$$       call tens_Int4_Re_Com(0d0,s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5,p4p
c$$$     -   5,C0456R,C0256R,C0246R,C0245R,C456R,C256R,C246R,C245R,C0456I
c$$$     -   ,C0256I,C0246I,C0245I,C456I,C256I,C246I,C245I,D02456R,D02456
c$$$     -   I,D2456R,D2456I)
c$$$       endif

       E012345=E0fin(p1sq,p2sq,p3sq,p4sq,s56,s12,s23,s34,s123,s234,D02
     -   345,D01345,D01245,D01235,D01234)
       E012345R=Dble(E012345)
       E012345I=DImag(E012345)
       call tens_red5_new_Re_Com(p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
     -   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2345R,
     -   D1345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01235I,
     -   D01234I,D2345I,D1345I,D1245I,D1235I,D1234I,E12345R,E12345I)
c$$$       If(PrintE) then
c$$$       write(*,*)"E12345"
c$$$       call print_EijR(E012345,E12345R,E12345I)
c$$$       call tens_Int5_Re_Com(0d0,p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
c$$$     -   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2345R,
c$$$     -   D1345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01235I,
c$$$     -   D01234I,D2345I,D1345I,D1245I,D1235I,D1234I,E012345R,E012345I
c$$$     -   ,E12345R,E12345I)
c$$$       endif

       E012346=E0fin(p1sq,p2sq,p3sq,s45,p6sq,s12,s23,s345,s123,s16,D02
     -   346,D01346,D01246,D01236,D01234)
       E012346R=Dble(E012346)
       E012346I=DImag(E012346)
       call tens_red5_new_Re_Com(p1sq,p2sq,p3sq,s45,p1p2,p1p3,p1p4+p1p
     -   5,p2p3,p2p4+p2p5,p3p4+p3p5,D02346R,D01346R,D01246R,D01236R,D
     -   01234R,D2346R,D1346R,D1246R,D1236R,D1234R,D02346I,D01346I,D0
     -   1246I,D01236I,D01234I,D2346I,D1346I,D1246I,D1236I,D1234I,E12
     -   346R,E12346I)
c$$$       If(PrintE) then
c$$$       write(*,*)"E12346"
c$$$       call print_EijR(E012346,E12346R,E12346I)
c$$$       call tens_Int5_Re_Com(0d0,p1sq,p2sq,p3sq,s45,p1p2,p1p3,p1p4+p1p
c$$$     -   5,p2p3,p2p4+p2p5,p3p4+p3p5,D02346R,D01346R,D01246R,D01236R,D
c$$$     -   01234R,D2346R,D1346R,D1246R,D1236R,D1234R,D02346I,D01346I,D0
c$$$     -   1246I,D01236I,D01234I,D2346I,D1346I,D1246I,D1236I,D1234I,E01
c$$$     -   2346R,E012346I,E12346R,E12346I)
c$$$       endif

       E012356=E0fin(p1sq,p2sq,s34,p5sq,p6sq,s12,s234,s345,s56,s16,D02
     -   356,D01356,D01256,D01236,D01235)
       E012356R=Dble(E012356)
       E012356I=DImag(E012356)
       call tens_red5_new_Re_Com(p1sq,p2sq,s34,p5sq,p1p2,p1p3+p1p4,p1p
     -   5,p2p3+p2p4,p2p5,p3p5+p4p5,D02356R,D01356R,D01256R,D01236R,D
     -   01235R,D2356R,D1356R,D1256R,D1236R,D1235R,D02356I,D01356I,D0
     -   1256I,D01236I,D01235I,D2356I,D1356I,D1256I,D1236I,D1235I,E12
     -   356R,E12356I)
c$$$       If(PrintE) then
c$$$       write(*,*)"E12356"
c$$$       call print_EijR(E012356,E12356R,E12356I)
c$$$       call tens_Int5_Re_Com(0d0,p1sq,p2sq,s34,p5sq,p1p2,p1p3+p1p4,p1p
c$$$     -   5,p2p3+p2p4,p2p5,p3p5+p4p5,D02356R,D01356R,D01256R,D01236R,D
c$$$     -   01235R,D2356R,D1356R,D1256R,D1236R,D1235R,D02356I,D01356I,D0
c$$$     -   1256I,D01236I,D01235I,D2356I,D1356I,D1256I,D1236I,D1235I,E01
c$$$     -   2356R,E012356I,E12356R,E12356I)
c$$$       endif

       E012456=E0fin(p1sq,s23,p4sq,p5sq,p6sq,s123,s234,s45,s56,s16,D02
     -   456,D01456,D01256,D01246,D01245)
       E012456R=Dble(E012456)
       E012456I=DImag(E012456)
       call tens_red5_new_Re_Com(p1sq,s23,p4sq,p5sq,p1p2+p1p3,p1p4,p1p
     -   5,p2p4+p3p4,p2p5+p3p5,p4p5,D02456R,D01456R,D01256R,D01246R,D
     -   01245R,D2456R,D1456R,D1256R,D1246R,D1245R,D02456I,D01456I,D0
     -   1256I,D01246I,D01245I,D2456I,D1456I,D1256I,D1246I,D1245I,E12
     -   456R,E12456I)
c$$$       If(PrintE) then
c$$$       write(*,*)"E12456"
c$$$       call print_EijR(E012456,E12456R,E12456I)
c$$$       call tens_Int5_Re_Com(0d0,p1sq,s23,p4sq,p5sq,p1p2+p1p3,p1p4,p1p
c$$$     -   5,p2p4+p3p4,p2p5+p3p5,p4p5,D02456R,D01456R,D01256R,D01246R,D
c$$$     -   01245R,D2456R,D1456R,D1256R,D1246R,D1245R,D02456I,D01456I,D0
c$$$     -   1256I,D01246I,D01245I,D2456I,D1456I,D1256I,D1246I,D1245I,E01
c$$$     -   2456R,E012456I,E12456R,E12456I)
c$$$       endif

       E023456=E0fin(p2sq,p3sq,p4sq,p5sq,s16,s23,s34,s45,s234,s345,D03
     -   456,D02456,D02356,D02346,D02345)
       E023456R=Dble(E023456)
       E023456I=DImag(E023456)
       call tens_red5_new_Re_Com(p2sq,p3sq,p4sq,p5sq,p2p3,p2p4,p2p5,p3
     -   p4,p3p5,p4p5,D03456R,D02456R,D02356R,D02346R,D02345R,D3456R,
     -   D2456R,D2356R,D2346R,D2345R,D03456I,D02456I,D02356I,D02346I,
     -   D02345I,D3456I,D2456I,D2356I,D2346I,D2345I,E23456R,E23456I)
c$$$       If(PrintE) then
c$$$       write(*,*)"E23456"
c$$$       call print_EijR(E023456,E23456R,E23456I)
c$$$       call tens_Int5_Re_Com(0d0,p2sq,p3sq,p4sq,p5sq,p2p3,p2p4,p2p5,p3
c$$$     -   p4,p3p5,p4p5,D03456R,D02456R,D02356R,D02346R,D02345R,D3456R,
c$$$     -   D2456R,D2356R,D2346R,D2345R,D03456I,D02456I,D02356I,D02346I,
c$$$     -   D02345I,D3456I,D2456I,D2356I,D2346I,D2345I,E023456R,E023456I
c$$$     -   ,E23456R,E23456I)
c$$$       endif

       E013456=E0fin(s12,p3sq,p4sq,p5sq,p6sq,s123,s34,s45,s56,s345,D03
     -   456,D01456,D01356,D01346,D01345)
       E013456R=Dble(E013456)
       E013456I=DImag(E013456)
       call tens_red5_new_Re_Com(s12,p3sq,p4sq,p5sq,p1p3+p2p3,p1p4+p2p
     -   4,p1p5+p2p5,p3p4,p3p5,p4p5,D03456R,D01456R,D01356R,D01346R,D
     -   01345R,D3456R,D1456R,D1356R,D1346R,D1345R,D03456I,D01456I,D0
     -   1356I,D01346I,D01345I,D3456I,D1456I,D1356I,D1346I,D1345I,E13
     -   456R,E13456I)
c$$$       If(PrintE) then
c$$$       write(*,*)"E13456"
c$$$       call print_EijR(E013456,E13456R,E13456I)
c$$$       call tens_Int5_Re_Com(0d0,s12,p3sq,p4sq,p5sq,p1p3+p2p3,p1p4+p2p
c$$$     -   4,p1p5+p2p5,p3p4,p3p5,p4p5,D03456R,D01456R,D01356R,D01346R,D
c$$$     -   01345R,D3456R,D1456R,D1356R,D1346R,D1345R,D03456I,D01456I,D0
c$$$     -   1356I,D01346I,D01345I,D3456I,D1456I,D1356I,D1346I,D1345I,E01
c$$$     -   3456R,E013456I,E13456R,E13456I)
c$$$       endif

       F0123456=F0fin(p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s45,s5
     -   6,s16,s123,s234,s345,E023456,E013456,E012456,E012356,E012346
     -   ,E012345)
       F0123456R=Dble(F0123456)
       F0123456I=DImag(F0123456)
       call tens_red6_new_Re_Com(p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1
     -   p4,p1p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,E023456R,E013456R,E012
     -   456R,E012356R,E012346R,E012345R,E23456R,E13456R,E12456R,E123
     -   56R,E12346R,E12345R,E023456I,E013456I,E012456I,E012356I,E012
     -   346I,E012345I,E23456I,E13456I,E12456I,E12356I,E12346I,E12345
     -   I,F123456R,F123456I,F123456)
c$$$       If(PrintF) then
c$$$       write(*,*)"F123456"
c$$$       call print_FijR(F0123456,F123456R,F123456I)
c$$$       Endif

c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call test_FFhlABE1(F(1))
       call test_FFhlABE2(F(198))
       call test_FFhlABE3(F(395))
       call test_FFhlABE4(F(592))
       call test_FFhlABE5(F(789))
       call test_FFhlABE6(F(986))
       call test_FFhlABE7(F(1183))
       call test_FFhlABE8(F(1380))
       call test_FFhlABE9(F(1577))
       call test_FFhlABE10(F(1774))
       call test_FFhlABE11(F(1971))
       call test_FFhlABE12(F(2168))
       call test_FFhlABE13(F(2365))
       call test_FFhlABE14(F(2562))
       call test_FFhlABE15(F(2759))
       call test_FFhlABE16(F(2956))
       call test_FFhlABE17(F(3153))
       call test_FFhlABE18(F(3350))
       call test_FFhlABE19(F(3547))
       call test_FFhlABE20(F(3744))
       call test_FFhlABE21(F(3941))
       call test_FFhlABE22(F(4138))
       call test_FFhlABE23(F(4335))
       call test_FFhlABE24(F(4532))
       call test_FFhlABE25(F(4729))
       call test_FFhlABE26(F(4926))
       call test_FFhlABE27(F(5123))
       call test_FFhlABE28(F(5320))
       call test_FFhlABE29(F(5517))
       call test_FFhlABE30(F(5714))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrc(p1,mup3)
       p1mup4 = dotrc(p1,mup4)
       p1mup5 = dotrc(p1,mup5)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup4 = dotrc(p2,mup4)
       p2mup5 = dotrc(p2,mup5)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup4 = dotrc(p3,mup4)
       p3mup5 = dotrc(p3,mup5)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup4 = dotrc(p4,mup4)
       p4mup5 = dotrc(p4,mup5)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrc(p5,mup3)
       p5mup4 = dotrc(p5,mup4)
       p5mup5 = dotrc(p5,mup5)
       p6mup2 = dotrc(p6,mup2)
       p6mup3 = dotrc(p6,mup3)
       p6mup4 = dotrc(p6,mup4)
       p6mup5 = dotrc(p6,mup5)
       mup2mup3 = dotcc(mup2,mup3)
       mup2mup4 = dotcc(mup2,mup4)
       mup2mup5 = dotcc(mup2,mup5)
       mup3mup4 = dotcc(mup3,mup4)
       mup3mup5 = dotcc(mup3,mup5)
       mup4mup5 = dotcc(mup4,mup5)
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call test_FaFunctionhexlineAbe1(p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2,
     -   p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6mu
     -   p3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mup
     -   5,mup4mup5,Fa(1))
       call test_FaFunctionhexlineAbe2(p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2,
     -   p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6mu
     -   p3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mup
     -   5,mup4mup5,Fa(12))
       call test_FaFunctionhexlineAbe3(p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2,
     -   p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6mu
     -   p3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mup
     -   5,mup4mup5,Fa(23))
       call test_FaFunctionhexlineAbe4(p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2,
     -   p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6mu
     -   p3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mup
     -   5,mup4mup5,Fa(34))
       call test_FaFunctionhexlineAbe5(p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2,
     -   p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6mu
     -   p3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mup
     -   5,mup4mup5,Fa(45))
       call test_FaFunctionhexlineAbe6(p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2,
     -   p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6mu
     -   p3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mup
     -   5,mup4mup5,Fa(56))
       call test_FaFunctionhexlineAbe7(p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2,
     -   p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6mu
     -   p3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mup
     -   5,mup4mup5,Fa(67))
       call test_FaFunctionhexlineAbe8(p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2,
     -   p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6mu
     -   p3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mup
     -   5,mup4mup5,Fa(78))
       call test_FaFunctionhexlineAbe9(p1mup2,p1mup3,p1mup4,p1mup5,p2m
     -   up2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2,
     -   p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6mu
     -   p3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mup
     -   5,mup4mup5,Fa(89))
       call test_FaFunctionhexlineAbe10(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(100))
       call test_FaFunctionhexlineAbe11(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(111))
       call test_FaFunctionhexlineAbe12(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(122))
       call test_FaFunctionhexlineAbe13(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(133))
       call test_FaFunctionhexlineAbe14(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(144))
       call test_FaFunctionhexlineAbe15(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(155))
       call test_FaFunctionhexlineAbe16(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(166))
       call test_FaFunctionhexlineAbe17(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(177))
       call test_FaFunctionhexlineAbe18(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(188))
       call test_FaFunctionhexlineAbe19(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(199))
       call test_FaFunctionhexlineAbe20(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(210))
       call test_FaFunctionhexlineAbe21(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(221))
       call test_FaFunctionhexlineAbe22(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(232))
       call test_FaFunctionhexlineAbe23(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(243))
       call test_FaFunctionhexlineAbe24(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(254))
       call test_FaFunctionhexlineAbe25(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(265))
       call test_FaFunctionhexlineAbe26(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(276))
       call test_FaFunctionhexlineAbe27(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(287))
       call test_FaFunctionhexlineAbe28(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(298))
       call test_FaFunctionhexlineAbe29(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(309))
       call test_FaFunctionhexlineAbe30(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(320))
       call test_FaFunctionhexlineAbe31(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(331))
       call test_FaFunctionhexlineAbe32(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(342))
       call test_FaFunctionhexlineAbe33(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(353))
       call test_FaFunctionhexlineAbe34(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(364))
       call test_FaFunctionhexlineAbe35(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(375))
       call test_FaFunctionhexlineAbe36(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(386))
       call test_FaFunctionhexlineAbe37(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(397))
       call test_FaFunctionhexlineAbe38(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(408))
       call test_FaFunctionhexlineAbe39(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(419))
       call test_FaFunctionhexlineAbe40(p1mup2,p1mup3,p1mup4,p1mup5,p2
     -   mup2,p2mup3,p2mup4,p2mup5,p3mup2,p3mup3,p3mup4,p3mup5,p4mup2
     -   ,p4mup3,p4mup4,p4mup5,p5mup2,p5mup3,p5mup4,p5mup5,p6mup2,p6m
     -   up3,p6mup4,p6mup5,mup2mup3,mup2mup4,mup2mup5,mup3mup4,mup3mu
     -   p5,mup4mup5,Fa(430))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the Matrix Element  
c************************************************************************************
c************************************************************************************
c************************************************************************************
       SMB(1) = SC3ccc(barpsi_p6,mup5,mup4,mup3,psi_p1,alpha)
       SMB(2) = SC3ccc(barpsi_p6,mup5,mup4,mup2,psi_p1,alpha)
       SMB(3) = SC3ccr(barpsi_p6,mup5,mup4,p3,psi_p1,alpha)
       SMB(4) = SC3ccr(barpsi_p6,mup5,mup4,p2,psi_p1,alpha)
       SMB(5) = SC5cccrr(barpsi_p6,mup5,mup4,mup2,p3,p2,psi_p1,alpha)
       SMB(6) = SC1c(barpsi_p6,mup5,psi_p1,alpha)
       SMB(7) = SC3ccr(barpsi_p6,mup5,mup4,p4,psi_p1,alpha)
       SMB(8) = SC3ccc(barpsi_p6,mup5,mup3,mup2,psi_p1,alpha)
       SMB(9) = SC3ccr(barpsi_p6,mup5,mup3,p2,psi_p1,alpha)
       SMB(10) = SC3ccr(barpsi_p6,mup5,mup2,p2,psi_p1,alpha)
       SMB(11) = SC3ccr(barpsi_p6,mup5,mup3,p3,psi_p1,alpha)
       SMB(12) = SC5ccccr(barpsi_p6,mup5,mup4,mup3,mup2,p2,psi_p1,alpha)
       SMB(13) = SC5cccrr(barpsi_p6,mup5,mup4,mup3,p3,p2,psi_p1,alpha)
       SMB(14) = SC5ccccr(barpsi_p6,mup5,mup4,mup3,mup2,p3,psi_p1,alpha)
       SMB(15) = SC5ccccr(barpsi_p6,mup5,mup4,mup3,mup2,p4,psi_p1,alpha)
       SMB(16) = SC5cccrr(barpsi_p6,mup5,mup4,mup3,p4,p2,psi_p1,alpha)
       SMB(17) = SC5cccrr(barpsi_p6,mup5,mup4,mup2,p4,p2,psi_p1,alpha)
       SMB(18) = SC5cccrr(barpsi_p6,mup5,mup4,mup3,p4,p3,psi_p1,alpha)
       SMB(19) = SC5cccrr(barpsi_p6,mup5,mup3,mup2,p3,p2,psi_p1,alpha)
       SMB(20) = SC7ccccrrr(barpsi_p6,mup5,mup4,mup3,mup2,p4,p3,p2,psi_p1,alpha)
       SMB(21) = SC1c(barpsi_p6,mup4,psi_p1,alpha)
       SMB(22) = SC3ccc(barpsi_p6,mup4,mup3,mup2,psi_p1,alpha)
       SMB(23) = SC3ccr(barpsi_p6,mup4,mup3,p2,psi_p1,alpha)
       SMB(24) = SC3ccr(barpsi_p6,mup4,mup2,p2,psi_p1,alpha)
       SMB(25) = SC3ccr(barpsi_p6,mup4,mup3,p3,psi_p1,alpha)
       SMB(26) = SC5cccrr(barpsi_p6,mup4,mup3,mup2,p3,p2,psi_p1,alpha)
       SMB(27) = SC1r(barpsi_p6,p2,psi_p1,alpha)
       SMB(28) = SC1r(barpsi_p6,p3,psi_p1,alpha)
       SMB(29) = SC1r(barpsi_p6,p4,psi_p1,alpha)
       SMB(30) = SC1c(barpsi_p6,mup2,psi_p1,alpha)
       SMB(31) = SC1c(barpsi_p6,mup3,psi_p1,alpha)
       SMB(32) = SC3ccr(barpsi_p6,mup3,mup2,p2,psi_p1,alpha)
       SMB(33) = SC3crr(barpsi_p6,mup3,p3,p2,psi_p1,alpha)
       SMB(34) = SC3ccr(barpsi_p6,mup3,mup2,p3,psi_p1,alpha)
       SMB(35) = SC3crr(barpsi_p6,mup2,p3,p2,psi_p1,alpha)
       SMB(36) = SC3ccr(barpsi_p6,mup3,mup2,p4,psi_p1,alpha)
       SMB(37) = SC3crr(barpsi_p6,mup3,p4,p2,psi_p1,alpha)
       SMB(38) = SC3crr(barpsi_p6,mup2,p4,p2,psi_p1,alpha)
       SMB(39) = SC3crr(barpsi_p6,mup3,p4,p3,psi_p1,alpha)
       SMB(40) = SC5ccrrr(barpsi_p6,mup3,mup2,p4,p3,p2,psi_p1,alpha)
       SMB(41) = SC3ccr(barpsi_p6,mup4,mup2,p3,psi_p1,alpha)
       SMB(42) = SC3ccr(barpsi_p6,mup4,mup2,p4,psi_p1,alpha)
       SMB(43) = SC3ccr(barpsi_p6,mup4,mup3,p4,psi_p1,alpha)
       SMB(44) = SC3ccr(barpsi_p6,mup5,mup2,p3,psi_p1,alpha)
       SMB(45) = SC3ccr(barpsi_p6,mup5,mup2,p4,psi_p1,alpha)
       SMB(46) = SC3ccr(barpsi_p6,mup5,mup3,p4,psi_p1,alpha)
       SMB(47) = SC5cccrr(barpsi_p6,mup5,mup3,mup2,p4,p2,psi_p1,alpha)
       SMB(48) = SC5cccrr(barpsi_p6,mup4,mup3,mup2,p4,p2,psi_p1,alpha)
       SMB(49) = SC3crr(barpsi_p6,mup4,p3,p2,psi_p1,alpha)
       SMB(50) = SC3crr(barpsi_p6,mup4,p4,p2,psi_p1,alpha)
       SMB(51) = SC3crr(barpsi_p6,mup4,p4,p3,psi_p1,alpha)
       SMB(52) = SC5ccrrr(barpsi_p6,mup4,mup3,p4,p3,p2,psi_p1,alpha)
       SMB(53) = SC5cccrr(barpsi_p6,mup4,mup3,mup2,p4,p3,psi_p1,alpha)
       SMB(54) = SC5ccrrr(barpsi_p6,mup4,mup2,p4,p3,p2,psi_p1,alpha)
       SMB(55) = SC3crr(barpsi_p6,mup2,p4,p3,psi_p1,alpha)
       SMB(56) = SC3crr(barpsi_p6,mup5,p3,p2,psi_p1,alpha)
       SMB(57) = SC3rrr(barpsi_p6,p4,p3,p2,psi_p1,alpha)
       SMB(58) = SC3crr(barpsi_p6,mup5,p4,p2,psi_p1,alpha)
       SMB(59) = SC3crr(barpsi_p6,mup5,p4,p3,psi_p1,alpha)
       SMB(60) = SC5ccrrr(barpsi_p6,mup5,mup2,p4,p3,p2,psi_p1,alpha)
       SMB(61) = SC5cccrr(barpsi_p6,mup5,mup4,mup2,p4,p3,psi_p1,alpha)
       SMB(62) = SC5ccrrr(barpsi_p6,mup5,mup4,p4,p3,p2,psi_p1,alpha)
       SMB(63) = SC5cccrr(barpsi_p6,mup5,mup3,mup2,p4,p3,psi_p1,alpha)
       SMB(64) = SC5ccrrr(barpsi_p6,mup5,mup3,p4,p3,p2,psi_p1,alpha)
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
       result(1) = -((F(5928)*(8*Fa(4)*SMB(6)+4*(Fa(2)*SMB(3)-Fa(3)*SM
     -   B(4)-Fa(5)*SMB(7)+Fa(7)*SMB(9)-Fa(8)*SMB(10)-Fa(9)*SMB(11))+
     -   s123*SMB(12)+s12*(-2*Fa(6)*SMB(8)+SMB(15))-2*(Fa(1)*SMB(1)-p
     -   4mup3*(s12*SMB(2)+SMB(5))+p4mup2*SMB(13)-p3mup2*SMB(16)+Fa(1
     -   0)*SMB(17)+Fa(11)*SMB(18)+Fa(6)*SMB(19))+SMB(20)))/(s12*s123
     -   *s56))
       result(2) = -((F(5929)*(8*Fa(4)*SMB(6)+4*(Fa(2)*SMB(3)-Fa(3)*SM
     -   B(4)-Fa(5)*SMB(7)+Fa(7)*SMB(9)-Fa(8)*SMB(10)-Fa(9)*SMB(11))+
     -   s123*SMB(12)+s12*(-2*Fa(6)*SMB(8)+SMB(15))-2*(Fa(1)*SMB(1)-p
     -   4mup3*(s12*SMB(2)+SMB(5))+p4mup2*SMB(13)-p3mup2*SMB(16)+Fa(1
     -   0)*SMB(17)+Fa(11)*SMB(18)+Fa(6)*SMB(19))+SMB(20)))/(s12*s123
     -   *s56))
       result(3) = -((F(5930)*(8*Fa(4)*SMB(6)+4*(Fa(2)*SMB(3)-Fa(3)*SM
     -   B(4)-Fa(5)*SMB(7)+Fa(7)*SMB(9)-Fa(8)*SMB(10)-Fa(9)*SMB(11))+
     -   s123*SMB(12)+s12*(-2*Fa(6)*SMB(8)+SMB(15))-2*(Fa(1)*SMB(1)-p
     -   4mup3*(s12*SMB(2)+SMB(5))+p4mup2*SMB(13)-p3mup2*SMB(16)+Fa(1
     -   0)*SMB(17)+Fa(11)*SMB(18)+Fa(6)*SMB(19))+SMB(20)))/(s12*s123
     -   *s56))
       result(4) = (2*(-8*Fa(15)*SMB(6)-4*(Fa(13)*SMB(3)+Fa(14)*SMB(4)
     -   -Fa(16)*SMB(7)-Fa(17)*SMB(9)-Fa(18)*SMB(11))+2*(Fa(12)*SMB(1
     -   )+Fa(19)*SMB(16)+Fa(20)*SMB(18))+F(1)*(-4*Fa(8)*SMB(10)+s123
     -   *SMB(12)+s12*(-2*Fa(6)*SMB(8)+SMB(15))+2*(p4mup3*(s12*SMB(2)
     -   +SMB(5))-p4mup2*SMB(13)-Fa(10)*SMB(17)-Fa(6)*SMB(19))+SMB(20
     -   ))))/(s12*s123*s56)
       result(5) = (2*(-8*Fa(26)*SMB(6)+4*(Fa(23)*SMB(3)-Fa(24)*SMB(4)
     -   +Fa(27)*SMB(7)+Fa(28)*SMB(10))+s123*F(7)*SMB(12)+2*(Fa(21)*S
     -   MB(1)+Fa(25)*SMB(5)+Fa(29)*SMB(17)-F(17)*(2*Fa(9)*SMB(11)+p4
     -   mup2*SMB(13)+Fa(11)*SMB(18)+Fa(6)*SMB(19)))+F(12)*(4*Fa(7)*S
     -   MB(9)+2*p3mup2*SMB(16)+SMB(20))+s12*(-(s123*F(8)*SMB(14))+F(
     -   9)*SMB(15)+4*Fa(30)*SMB(44)+2*(Fa(22)*SMB(2)+F(6)*Fa(6)*SMB(
     -   8)+Fa(31)*SMB(61)))))/(s12*s123*s56)
       result(6) = (2*(-8*Fa(34)*SMB(6)+s123*F(33)*SMB(12)+s12*F(23)*S
     -   MB(15)+F(24)*(4*Fa(2)*SMB(3)+2*p3mup2*SMB(16)+SMB(20))+4*(Fa
     -   (33)*SMB(4)-Fa(36)*SMB(9)+Fa(37)*SMB(10)+Fa(38)*SMB(11)+s123
     -   *Fa(40)*SMB(46))-2*(Fa(32)*SMB(1)-p4mup3*(s12*F(23)*SMB(2)+F
     -   (24)*SMB(5))-s12*Fa(35)*SMB(8)+F(25)*(2*Fa(5)*SMB(7)+p4mup2*
     -   SMB(13)+Fa(10)*SMB(17)+Fa(11)*SMB(18))-Fa(39)*SMB(19)-s123*F
     -   a(41)*SMB(47))))/(s12*s123*s56)
       result(7) = (2*(-(s123*F(34)*SMB(12))-s12*F(34)*SMB(15)-F(34)*S
     -   MB(20)-8*(F(34)*Fa(4)*SMB(6)+s56*Fa(42)*SMB(21))-4*(F(34)*Fa
     -   (2)*SMB(3)-F(34)*Fa(3)*SMB(4)-F(34)*Fa(5)*SMB(7)+F(34)*Fa(7)
     -   *SMB(9)-F(34)*Fa(8)*SMB(10)-F(34)*Fa(9)*SMB(11)+s56*Fa(44)*S
     -   MB(23)-s56*Fa(45)*SMB(24)-s56*Fa(46)*SMB(25))+2*(F(34)*Fa(1)
     -   *SMB(1)-p4mup3*F(34)*(s12*SMB(2)+SMB(5))+s12*F(34)*Fa(6)*SMB
     -   (8)+p4mup2*F(34)*SMB(13)-p3mup2*F(34)*SMB(16)+F(34)*Fa(10)*S
     -   MB(17)+F(34)*Fa(11)*SMB(18)+F(34)*Fa(6)*SMB(19)+s12*s56*Fa(4
     -   3)*SMB(22)+s56*Fa(43)*SMB(26))))/(s12*s123*s56)
       result(8) = (Fa(48)*SMB(2)+Fa(49)*SMB(3)+Fa(50)*SMB(4)+Fa(52)*S
     -   MB(6)+Fa(55)*SMB(10)+8*(Fa(54)*SMB(9)+Fa(56)*SMB(11))+F(101)
     -   *SMB(15)+Fa(58)*SMB(16)+Fa(59)*SMB(17)-4*(Fa(47)*SMB(1)+Fa(5
     -   1)*SMB(5)+Fa(53)*SMB(7)-Fa(57)*SMB(13)-Fa(60)*SMB(18)+Fa(6)*
     -   (F(37)*SMB(8)-F(76)*SMB(19)))+2*(s123*(F(97)*SMB(12)+F(100)*
     -   SMB(14))-F(76)*SMB(20)-4*(Fa(61)*SMB(44)+Fa(62)*SMB(56))-2*F
     -   a(63)*SMB(61))+Fa(64)*SMB(62))/(s123*s56)
       result(9) = (Fa(65)*SMB(1)-Fa(70)*SMB(6)+Fa(73)*SMB(9)+Fa(74)*S
     -   MB(10)-Fa(77)*SMB(19)+2*(F(177)*SMB(12)+F(181)*(2*p3mup2*SMB
     -   (16)+SMB(20)))-4*(Fa(69)*SMB(5)-Fa(76)*SMB(17)-F(178)*(p4mup
     -   2*SMB(13)+Fa(11)*SMB(18))+Fa(81)*SMB(47))-8*(Fa(67)*SMB(3)-F
     -   a(68)*SMB(4)-Fa(71)*SMB(7)+Fa(75)*SMB(11)+Fa(80)*SMB(46)+Fa(
     -   82)*SMB(56)-Fa(83)*SMB(58)+Fa(84)*SMB(59))+Fa(85)*SMB(60)+s1
     -   2*(Fa(66)*SMB(2)+Fa(72)*SMB(8)+F(179)*SMB(14)+F(180)*SMB(15)
     -   +Fa(78)*SMB(44)+Fa(79)*SMB(45)+4*(Fa(86)*SMB(61)-Fa(87)*SMB(
     -   63))))/(s12*s56)
       result(10) = (Fa(88)*SMB(1)+p4mup3*(s12*F(197)*SMB(2)+4*F(198)*
     -   SMB(5))+Fa(90)*SMB(6)+Fa(92)*SMB(9)+Fa(93)*SMB(10)+Fa(94)*SM
     -   B(11)+Fa(95)*SMB(19)+4*Fa(101)*SMB(26)-16*(Fa(96)*SMB(21)-Fa
     -   (420)*SMB(27)-Fa(108)*SMB(28))+s12*(Fa(91)*SMB(8)+F(214)*SMB
     -   (15)+Fa(97)*SMB(22)+8*Fa(117)*SMB(30)+Fa(121)*SMB(36))-8*(Fa
     -   (98)*SMB(23)-Fa(99)*SMB(24)-Fa(100)*SMB(25)+2*Fa(421)*SMB(29
     -   )+Fa(118)*SMB(31)+Fa(120)*SMB(33)-Fa(117)*SMB(35)-Fa(122)*SM
     -   B(37)+Fa(123)*SMB(38)+Fa(124)*SMB(39))+Fa(121)*SMB(40)+2*(F(
     -   198)*(4*(Fa(2)*SMB(3)+Fa(89)*SMB(4))-2*(2*Fa(5)*SMB(7)+p4mup
     -   2*SMB(13)-p3mup2*SMB(16)+Fa(10)*SMB(17)+Fa(11)*SMB(18))+SMB(
     -   20))+s123*(F(213)*SMB(12)+4*(Fa(125)*SMB(43)-Fa(126)*SMB(46)
     -   )+2*(Fa(119)*SMB(32)-Fa(127)*SMB(47)+Fa(128)*SMB(48)))))/(s1
     -   2*s123)
       result(11) = (Fa(129)*SMB(1)+Fa(130)*SMB(2)+Fa(131)*SMB(3)+Fa(1
     -   32)*SMB(4)+Fa(133)*SMB(5)+Fa(134)*SMB(6)+Fa(135)*SMB(7)+Fa(1
     -   36)*SMB(8)+Fa(137)*SMB(9)+Fa(138)*SMB(10)+Fa(139)*SMB(11)+F(
     -   440)*SMB(12)+Fa(140)*SMB(13)+F(445)*SMB(14)+F(446)*SMB(15)+F
     -   a(141)*SMB(16)+Fa(142)*SMB(17)+Fa(143)*SMB(18)+Fa(144)*SMB(1
     -   9)+F(463)*SMB(20)+Fa(145)*SMB(44)+Fa(146)*SMB(45)+Fa(147)*SM
     -   B(46)+Fa(148)*SMB(47)+Fa(149)*SMB(56)+Fa(150)*SMB(58)-Fa(151
     -   )*SMB(59)+Fa(152)*SMB(60)+Fa(153)*SMB(61)+Fa(154)*SMB(62)+Fa
     -   (155)*SMB(63)+Fa(156)*SMB(64))/s56
       result(12) = (Fa(159)*SMB(3)+Fa(160)*SMB(4)+Fa(161)*SMB(5)+Fa(4
     -   22)*SMB(6)+Fa(166)*SMB(7)+Fa(168)*SMB(9)+Fa(169)*SMB(10)+Fa(
     -   170)*SMB(11)+2*F(916)*SMB(12)+p3mup2*F(920)*SMB(16)+Fa(171)*
     -   SMB(17)+F(917)*(p4mup2*SMB(13)+Fa(11)*SMB(18))+Fa(172)*SMB(1
     -   9)+F(928)*SMB(20)+Fa(173)*SMB(21)+Fa(175)*SMB(23)+Fa(176)*SM
     -   B(24)+Fa(423)*SMB(27)-Fa(424)*SMB(28)+Fa(425)*SMB(29)+Fa(218
     -   )*SMB(31)+Fa(219)*SMB(32)+Fa(427)*SMB(35)+Fa(428)*SMB(38)+Fa
     -   (237)*SMB(40)+Fa(240)*SMB(43)+Fa(243)*SMB(46)+Fa(244)*SMB(47
     -   )+Fa(245)*SMB(48)-Fa(250)*SMB(54)+8*(Fa(177)*SMB(25)+Fa(220)
     -   *SMB(33)+Fa(229)*SMB(37)+Fa(236)*SMB(39)+Fa(246)*SMB(49)-Fa(
     -   247)*SMB(50)+Fa(248)*SMB(51)+Fa(252)*SMB(56)-Fa(253)*SMB(58)
     -   +Fa(254)*SMB(59))+4*(Fa(157)*SMB(1)+Fa(178)*SMB(26)+Fa(255)*
     -   SMB(60))+s12*(Fa(158)*SMB(2)+Fa(167)*SMB(8)+F(918)*SMB(14)+2
     -   *F(919)*SMB(15)+Fa(174)*SMB(22)+Fa(221)*SMB(34)+Fa(228)*SMB(
     -   36)+Fa(239)*SMB(42)+Fa(241)*SMB(44)+8*(Fa(426)*SMB(30)+Fa(25
     -   1)*SMB(55))+4*(Fa(238)*SMB(41)-Fa(242)*SMB(45)-Fa(249)*SMB(5
     -   3)-Fa(256)*SMB(61)+Fa(257)*SMB(63))))/s12
       result(13) = Fa(258)*SMB(1)+Fa(259)*SMB(2)+Fa(260)*SMB(3)+Fa(26
     -   1)*SMB(4)+Fa(429)*SMB(6)+Fa(267)*SMB(7)+Fa(268)*SMB(8)+Fa(26
     -   9)*SMB(9)+Fa(270)*SMB(10)+Fa(271)*SMB(11)+Fa(272)*SMB(13)+F(
     -   1668)*SMB(14)+F(1669)*SMB(15)+2*(F(1662)*SMB(12)-F(1690)*SMB
     -   (20))+Fa(430)*SMB(21)+Fa(283)*SMB(22)+Fa(284)*SMB(23)+Fa(285
     -   )*SMB(24)+Fa(286)*SMB(25)+Fa(431)*SMB(27)+Fa(432)*SMB(28)+Fa
     -   (434)*SMB(30)+Fa(435)*SMB(31)+Fa(340)*SMB(32)+Fa(437)*SMB(35
     -   )+Fa(355)*SMB(36)+Fa(371)*SMB(41)+Fa(373)*SMB(43)+Fa(374)*SM
     -   B(44)+Fa(375)*SMB(45)+Fa(442)*SMB(56)+16*(Fa(433)*SMB(29)+Fa
     -   (443)*SMB(57))+8*(Fa(436)*SMB(33)+Fa(356)*SMB(37)-Fa(438)*SM
     -   B(38)+Fa(439)*SMB(39)+Fa(379)*SMB(49)+Fa(440)*SMB(50)+Fa(441
     -   )*SMB(51)+Fa(395)*SMB(55)+Fa(444)*SMB(58)+Fa(445)*SMB(59))-4
     -   *(Fa(262)*SMB(5)+Fa(273)*SMB(16)-Fa(274)*SMB(17)-Fa(275)*SMB
     -   (18)-Fa(276)*SMB(19)+Fa(287)*SMB(26)+Fa(348)*SMB(34)+Fa(370)
     -   *SMB(40)+Fa(372)*SMB(42)+Fa(376)*SMB(46)-Fa(377)*SMB(47)+Fa(
     -   378)*SMB(48)-Fa(392)*SMB(52)-Fa(393)*SMB(53)-Fa(394)*SMB(54)
     -   -Fa(415)*SMB(60)-Fa(416)*SMB(61)-Fa(417)*SMB(62)+Fa(418)*SMB
     -   (63)+Fa(419)*SMB(64))
       If (ngluon.eq.0) then
       result(1)=-(result(1)+result(2)+result(3)+result(4)+result(5)+r
     -   esult(6)+result(7)+result(8)+result(9)+result(10)+result(11)
     -   +result(12)+result(13))
       elseif (ngluon.eq.1) then
       If (posgluon.eq.2) then
       result(1)=-(result(1)+result(2)+result(3)+result(5)+result(6)+r
     -   esult(7)+result(9)+result(10)+result(12))
       result(2)=-(result(4)+result(8)+result(11)+result(13))
       elseIf (posgluon.eq.3) then
       result(1)=-(result(1)+result(2)+result(3)+result(4)+result(6)+r
     -   esult(7)+result(10))
       result(2)=-(result(5)+result(8)+result(9)+result(11)+result(12)
     -   +result(13))
       elseIf (posgluon.eq.4) then
       result(1)=-(result(1)+result(2)+result(3)+result(4)+result(5)+r
     -   esult(7)+result(8))
       result(2)=-(result(6)+result(9)+result(10)+result(11)+result(12
     -   )+result(13))
       elseIf (posgluon.eq.5) then
       result(1)=-(result(1)+result(2)+result(3)+result(4)+result(5)+r
     -   esult(6)+result(8)+result(9)+result(11))
       result(2)=-(result(7)+result(10)+result(12)+result(13))
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     -   . Look to the heading for explanation" 
       endif
       elseif (ngluon.eq.2) then
c posgluon==24 -> gluons on position 2 and 4
c result(1): coefficient of t^b t^a color structure (i.e. the born color structure)
c result(2): coefficient of delta^ab color structure
       if(posgluon.eq.23) then
       result(1)=+(result(4)+result(5)+result(8)+result(9)+result(11)+
     -             result(12)+result(13)-8d0*(result(1)+result(2)+
     -             result(3)+result(6)+result(7)+result(10)))/6d0
       result(2)=-(result(8)+result(11)+result(13))/4d0
       elseif(posgluon.eq.24) then
       result(1)=+(result(4)+result(6)+result(8)+result(9)+result(10)+
     -             result(11)+result(12)+result(13)-8d0*(result(1)+
     -             result(2)+result(3)+result(5)+result(7)))/6d0
       result(2)=-(result(11)+result(13))/4d0
       elseif(posgluon.eq.25) then
       result(1)=+(result(4)+result(7)+result(8)+result(10)+result(11)+
     -             result(12)+result(13)-8d0*(result(1)+result(2)+
     -             result(3)+result(5)+result(6)+result(9)))/6d0
       result(2)=-(result(13))/4d0
       elseif(posgluon.eq.34) then
       result(1)=+(result(5)+result(6)+result(8)+result(9)+result(10)+
     -             result(11)+result(12)+result(13)-8d0*(result(1)+
     -             result(2)+result(3)+result(4)+result(7)))/6d0
       result(2)=-(result(9)+result(11)+result(12)+result(13))/4d0
       elseif(posgluon.eq.35) then
       result(1)=+(result(5)+result(7)+result(8)+result(9)+result(10)+
     -             result(11)+result(12)+result(13)-8d0*(result(1)+
     -             result(2)+result(3)+result(4)+result(6)))/6d0
       result(2)=-(result(12)+result(13))/4d0
       elseif(posgluon.eq.45) then
       result(1)=+(result(6)+result(7)+result(9)+result(10)+result(11)+
     -             result(12)+result(13)-8d0*(result(1)+result(2)+
     -             result(3)+result(4)+result(5)+result(8)))/6d0
       result(2)=-(result(10)+result(12)+result(13))/4d0
       else
       write(*,*) "Error: Wrong position of gluon in HexlineAbe"
       endif
       endif
       resultborn = (2*(p4mup2*s12+(p1mup2+p2mup2)*s123)*SMB(1)-2*p4mu
     -   p3*s12*SMB(2)-4*(p1mup2+p2mup2)*p4mup3*SMB(3)-4*(p1mup3*p4mu
     -   p2+p2mup3*p4mup2+p3mup3*p4mup2-p3mup2*p4mup3)*SMB(4)-2*p4mup
     -   3*SMB(5)-8*(p1mup2+p2mup2)*(p1mup3+p2mup3+p3mup3)*(p1mup4+p2
     -   mup4+p3mup4+p4mup4)*SMB(6)+4*(p1mup2+p2mup2)*(p1mup3+p2mup3+
     -   p3mup3)*SMB(7)+2*(p1mup4+p2mup4+p3mup4+p4mup4)*s12*SMB(8)-4*
     -   p3mup2*(p1mup4+p2mup4+p3mup4+p4mup4)*SMB(9)+4*(p1mup3+p2mup3
     -   +p3mup3)*(p1mup4+p2mup4+p3mup4+p4mup4)*SMB(10)+4*(p1mup2+p2m
     -   up2)*(p1mup4+p2mup4+p3mup4+p4mup4)*SMB(11)-s123*SMB(12)+2*p4
     -   mup2*SMB(13)-s12*SMB(15)-2*p3mup2*SMB(16)+2*(p1mup3+p2mup3+p
     -   3mup3)*SMB(17)+2*(p1mup2+p2mup2)*SMB(18)+2*(p1mup4+p2mup4+p3
     -   mup4+p4mup4)*SMB(19)-SMB(20))/(s12*s123*s56)
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
