       subroutine gluonHexHiggsEven(m,p1,p2,p3,p4,p5,p6,mup1,mup2,mup3
     -   ,mup4,mup5,comp,result)
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
       Complex*16   mup1(0:3),mup2(0:3),mup3(0:3),mup4(0:3),mup5(0:3)
       Complex*16   p1mup1, p1mup2, p1mup3, p1mup4, p1mup5, p2mup1, 
     -          p2mup2, p2mup3, p2mup4, p2mup5, p3mup1, p3mup2, 
     -          p3mup3, p3mup4, p3mup5, p4mup1, p4mup2, p4mup3, 
     -          p4mup4, p4mup5, p5mup1, p5mup2, p5mup3, p5mup4, p5mup5
       Complex*16   mup1mup2, mup1mup3, mup1mup4, mup1mup5, mup2mup3, 
     -          mup2mup4, mup2mup5, mup3mup4, mup3mup5, mup4mup5
       Real*8 dotrr
       Complex*16 B0tM,C0fin1M,D04,E01M,F0fin1M
       EXTERNAL dotrr,B0tM,C0fin1M,D04,E01M,F0fin1M
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
       Integer comp 
      Complex*16 F(3823),Fa(50)
      Real*8 FI(3823),FR(3823),P(581)
      Complex*16 dotrc,dotcc,result
       Real*8 p5t,m,msq
      EXTERNAL   dotrc,dotcc
      Integer alpha
       COMMON/FFunctions/F
      Save/FFunctions/
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
       msq=m*m
      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False. 
      PrintF=.False.
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
        if (comp.gt.0d0) then    
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
cfc
cfc       do alpha=0,3
cfc       print*, "Conse", p1(alpha)+ p2(alpha)+ p3(alpha)+
cfc     &       p4(alpha)+ p5(alpha)+ p6(alpha)
cfc       enddo

      s123=p1sq+p2sq+p3sq+2*(p1p2+p1p3+p2p3) 
      s234=p2sq+p3sq+p4sq+2*(p2p3+p2p4+p3p4) 
      s345=p3sq+p4sq+p5sq+2*(p3p4+p3p5+p4p5)
c    Calling C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************
       B012=B0tM(m,p1sq)
c FC       If(PrintB) then
c FC       write(*,*)"B12"
c FC       write(*,*) "B12=", B012
c FC       endif

       B023=B0tM(m,p2sq)
c FC       If(PrintB) then
c FC       write(*,*)"B23"
c FC       write(*,*) "B23=", B023
c FC       endif

       B034=B0tM(m,p3sq)
c FC       If(PrintB) then
c FC       write(*,*)"B34"
c FC       write(*,*) "B34=", B034
c FC       endif

       B045=B0tM(m,p4sq)
c FC       If(PrintB) then
c FC       write(*,*)"B45"
c FC       write(*,*) "B45=", B045
c FC       endif

       B056=B0tM(m,p5sq)
c FC       If(PrintB) then
c FC       write(*,*)"B56"
c FC       write(*,*) "B56=", B056
c FC       endif

       B016=B0tM(m,p6sq)
c FC       If(PrintB) then
c FC       write(*,*)"B16"
c FC       write(*,*) "B16=", B016
c FC       endif

       B013=B0tM(m,s12)
c FC       If(PrintB) then
c FC       write(*,*)"B13"
c FC       write(*,*) "B13=", B013
c FC       endif

       B014=B0tM(m,s123)
c FC       If(PrintB) then
c FC       write(*,*)"B14"
c FC       write(*,*) "B14=", B014
c FC       endif

       B026=B0tM(m,s16)
c FC       If(PrintB) then
c FC       write(*,*)"B26"
c FC       write(*,*) "B26=", B026
c FC       endif

       B024=B0tM(m,s23)
c FC       If(PrintB) then
c FC       write(*,*)"B24"
c FC       write(*,*) "B24=", B024
c FC       endif

       B025=B0tM(m,s234)
c FC       If(PrintB) then
c FC       write(*,*)"B25"
c FC       write(*,*) "B25=", B025
c FC       endif

       B035=B0tM(m,s34)
c FC       If(PrintB) then
c FC       write(*,*)"B35"
c FC       write(*,*) "B35=", B035
c FC       endif

       B036=B0tM(m,s345)
c FC       If(PrintB) then
c FC       write(*,*)"B36"
c FC       write(*,*) "B36=", B036
c FC       endif

       B046=B0tM(m,s45)
c FC       If(PrintB) then
c FC       write(*,*)"B46"
c FC       write(*,*) "B46=", B046
c FC       endif

       B015=B0tM(m,s56)
c FC       If(PrintB) then
c FC       write(*,*)"B15"
c FC       write(*,*) "B15=", B015
c FC       endif

       C0123=C0fin1M(m,p1sq,p2sq,s12)
       call tens_red3_new_Re_Com_1M(m,p1sq,p2sq,s12,B023,B013,B012,C01
     -   23,C0123R,C0123I,C123R,C123I)
c FC       If(PrintC) then
c FC       write(*,*)"C123"
c FC       call print_CijR(C0123,C123R,C123I)
c FC       call tens_Int3_Re_Com(m,p1sq,p2sq,s12,B023,B013,B012,C0123R,C01
c FC     -   23I,C123R,C123I)
c FC       Endif

       C0126=C0fin1M(m,p1sq,s16,p6sq)
       call tens_red3_new_Re_Com_1M(m,p1sq,s16,p6sq,B026,B016,B012,C01
     -   26,C0126R,C0126I,C126R,C126I)
c FC       If(PrintC) then
c FC       write(*,*)"C126"
c FC       call print_CijR(C0126,C126R,C126I)
c FC       call tens_Int3_Re_Com(m,p1sq,s16,p6sq,B026,B016,B012,C0126R,C01
c FC     -   26I,C126R,C126I)
c FC       Endif

       C0124=C0fin1M(m,p1sq,s23,s123)
       call tens_red3_new_Re_Com_1M(m,p1sq,s23,s123,B024,B014,B012,C01
     -   24,C0124R,C0124I,C124R,C124I)
c FC       If(PrintC) then
c FC       write(*,*)"C124"
c FC       call print_CijR(C0124,C124R,C124I)
c FC       call tens_Int3_Re_Com(m,p1sq,s23,s123,B024,B014,B012,C0124R,C01
c FC     -   24I,C124R,C124I)
c FC       Endif

       C0125=C0fin1M(m,p1sq,s234,s56)
       call tens_red3_new_Re_Com_1M(m,p1sq,s234,s56,B025,B015,B012,C01
     -   25,C0125R,C0125I,C125R,C125I)
c FC       If(PrintC) then
c FC       write(*,*)"C125"
c FC       call print_CijR(C0125,C125R,C125I)
c FC       call tens_Int3_Re_Com(m,p1sq,s234,s56,B025,B015,B012,C0125R,C01
c FC     -   25I,C125R,C125I)
c FC       Endif

       C0234=C0fin1M(m,p2sq,p3sq,s23)
       call tens_red3_new_Re_Com_1M(m,p2sq,p3sq,s23,B034,B024,B023,C02
     -   34,C0234R,C0234I,C234R,C234I)
c FC       If(PrintC) then
c FC       write(*,*)"C234"
c FC       call print_CijR(C0234,C234R,C234I)
c FC       call tens_Int3_Re_Com(m,p2sq,p3sq,s23,B034,B024,B023,C0234R,C02
c FC     -   34I,C234R,C234I)
c FC       Endif

       C0235=C0fin1M(m,p2sq,s34,s234)
       call tens_red3_new_Re_Com_1M(m,p2sq,s34,s234,B035,B025,B023,C02
     -   35,C0235R,C0235I,C235R,C235I)
c FC       If(PrintC) then
c FC       write(*,*)"C235"
c FC       call print_CijR(C0235,C235R,C235I)
c FC       call tens_Int3_Re_Com(m,p2sq,s34,s234,B035,B025,B023,C0235R,C02
c FC     -   35I,C235R,C235I)
c FC       Endif

       C0236=C0fin1M(m,p2sq,s345,s16)
       call tens_red3_new_Re_Com_1M(m,p2sq,s345,s16,B036,B026,B023,C02
     -   36,C0236R,C0236I,C236R,C236I)
c FC       If(PrintC) then
c FC       write(*,*)"C236"
c FC       call print_CijR(C0236,C236R,C236I)
c FC       call tens_Int3_Re_Com(m,p2sq,s345,s16,B036,B026,B023,C0236R,C02
c FC     -   36I,C236R,C236I)
c FC       Endif

       C0345=C0fin1M(m,p3sq,p4sq,s34)
       call tens_red3_new_Re_Com_1M(m,p3sq,p4sq,s34,B045,B035,B034,C03
     -   45,C0345R,C0345I,C345R,C345I)
c FC       If(PrintC) then
c FC       write(*,*)"C345"
c FC       call print_CijR(C0345,C345R,C345I)
c FC       call tens_Int3_Re_Com(m,p3sq,p4sq,s34,B045,B035,B034,C0345R,C03
c FC     -   45I,C345R,C345I)
c FC       Endif

       C0346=C0fin1M(m,p3sq,s45,s345)
       call tens_red3_new_Re_Com_1M(m,p3sq,s45,s345,B046,B036,B034,C03
     -   46,C0346R,C0346I,C346R,C346I)
c FC       If(PrintC) then
c FC       write(*,*)"C346"
c FC       call print_CijR(C0346,C346R,C346I)
c FC       call tens_Int3_Re_Com(m,p3sq,s45,s345,B046,B036,B034,C0346R,C03
c FC     -   46I,C346R,C346I)
c FC       Endif

       C0456=C0fin1M(m,p4sq,p5sq,s45)
       call tens_red3_new_Re_Com_1M(m,p4sq,p5sq,s45,B056,B046,B045,C04
     -   56,C0456R,C0456I,C456R,C456I)
c FC       If(PrintC) then
c FC       write(*,*)"C456"
c FC       call print_CijR(C0456,C456R,C456I)
c FC       call tens_Int3_Re_Com(m,p4sq,p5sq,s45,B056,B046,B045,C0456R,C04
c FC     -   56I,C456R,C456I)
c FC       Endif

       C0134=C0fin1M(m,s12,p3sq,s123)
       call tens_red3_new_Re_Com_1M(m,s12,p3sq,s123,B034,B014,B013,C01
     -   34,C0134R,C0134I,C134R,C134I)
c FC       If(PrintC) then
c FC       write(*,*)"C134"
c FC       call print_CijR(C0134,C134R,C134I)
c FC       call tens_Int3_Re_Com(m,s12,p3sq,s123,B034,B014,B013,C0134R,C01
c FC     -   34I,C134R,C134I)
c FC       Endif

       C0135=C0fin1M(m,s12,s34,s56)
       call tens_red3_new_Re_Com_1M(m,s12,s34,s56,B035,B015,B013,C0135
     -   ,C0135R,C0135I,C135R,C135I)
c FC       If(PrintC) then
c FC       write(*,*)"C135"
c FC       call print_CijR(C0135,C135R,C135I)
c FC       call tens_Int3_Re_Com(m,s12,s34,s56,B035,B015,B013,C0135R,C0135
c FC     -   I,C135R,C135I)
c FC       Endif

       C0136=C0fin1M(m,s12,s345,p6sq)
       call tens_red3_new_Re_Com_1M(m,s12,s345,p6sq,B036,B016,B013,C01
     -   36,C0136R,C0136I,C136R,C136I)
c FC       If(PrintC) then
c FC       write(*,*)"C136"
c FC       call print_CijR(C0136,C136R,C136I)
c FC       call tens_Int3_Re_Com(m,s12,s345,p6sq,B036,B016,B013,C0136R,C01
c FC     -   36I,C136R,C136I)
c FC       Endif

       C0145=C0fin1M(m,s123,p4sq,s56)
       call tens_red3_new_Re_Com_1M(m,s123,p4sq,s56,B045,B015,B014,C01
     -   45,C0145R,C0145I,C145R,C145I)
c FC       If(PrintC) then
c FC       write(*,*)"C145"
c FC       call print_CijR(C0145,C145R,C145I)
c FC       call tens_Int3_Re_Com(m,s123,p4sq,s56,B045,B015,B014,C0145R,C01
c FC     -   45I,C145R,C145I)
c FC       Endif

       C0146=C0fin1M(m,s123,s45,p6sq)
       call tens_red3_new_Re_Com_1M(m,s123,s45,p6sq,B046,B016,B014,C01
     -   46,C0146R,C0146I,C146R,C146I)
c FC       If(PrintC) then
c FC       write(*,*)"C146"
c FC       call print_CijR(C0146,C146R,C146I)
c FC       call tens_Int3_Re_Com(m,s123,s45,p6sq,B046,B016,B014,C0146R,C01
c FC     -   46I,C146R,C146I)
c FC       Endif

       C0245=C0fin1M(m,s23,p4sq,s234)
       call tens_red3_new_Re_Com_1M(m,s23,p4sq,s234,B045,B025,B024,C02
     -   45,C0245R,C0245I,C245R,C245I)
c FC       If(PrintC) then
c FC       write(*,*)"C245"
c FC       call print_CijR(C0245,C245R,C245I)
c FC       call tens_Int3_Re_Com(m,s23,p4sq,s234,B045,B025,B024,C0245R,C02
c FC     -   45I,C245R,C245I)
c FC       Endif

       C0246=C0fin1M(m,s23,s45,s16)
       call tens_red3_new_Re_Com_1M(m,s23,s45,s16,B046,B026,B024,C0246
     -   ,C0246R,C0246I,C246R,C246I)
c FC       If(PrintC) then
c FC       write(*,*)"C246"
c FC       call print_CijR(C0246,C246R,C246I)
c FC       call tens_Int3_Re_Com(m,s23,s45,s16,B046,B026,B024,C0246R,C0246
c FC     -   I,C246R,C246I)
c FC       Endif

       C0256=C0fin1M(m,s234,p5sq,s16)
       call tens_red3_new_Re_Com_1M(m,s234,p5sq,s16,B056,B026,B025,C02
     -   56,C0256R,C0256I,C256R,C256I)
c FC       If(PrintC) then
c FC       write(*,*)"C256"
c FC       call print_CijR(C0256,C256R,C256I)
c FC       call tens_Int3_Re_Com(m,s234,p5sq,s16,B056,B026,B025,C0256R,C02
c FC     -   56I,C256R,C256I)
c FC       Endif

       C0356=C0fin1M(m,s34,p5sq,s345)
       call tens_red3_new_Re_Com_1M(m,s34,p5sq,s345,B056,B036,B035,C03
     -   56,C0356R,C0356I,C356R,C356I)
c FC       If(PrintC) then
c FC       write(*,*)"C356"
c FC       call print_CijR(C0356,C356R,C356I)
c FC       call tens_Int3_Re_Com(m,s34,p5sq,s345,B056,B036,B035,C0356R,C03
c FC     -   56I,C356R,C356I)
c FC       Endif

       C0156=C0fin1M(m,s56,p5sq,p6sq)
       call tens_red3_new_Re_Com_1M(m,s56,p5sq,p6sq,B056,B016,B015,C01
     -   56,C0156R,C0156I,C156R,C156I)
c FC       If(PrintC) then
c FC       write(*,*)"C156"
c FC       call print_CijR(C0156,C156R,C156I)
c FC       call tens_Int3_Re_Com(m,s56,p5sq,p6sq,B056,B016,B015,C0156R,C01
c FC     -   56I,C156R,C156I)
c FC       Endif

       D01234=D04(p1sq,p2sq,p3sq,s123,s12,s23,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0
     -   234R,C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C01
     -   34I,C0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D01
     -   234I,D1234R,D1234I)
c FC       If(PrintD) then
c FC       write(*,*)"D1234"
c FC       call print_DijR(D01234,D1234R,D1234I)
c FC       call tens_Int4_Re_Com(m,p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,C0
c FC     -   134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C0134I,C01
c FC     -   24I,C0123I,C234I,C134I,C124I,C123I,D01234R,D01234I,D1234R,D1
c FC     -   234I)
c FC       endif

       D01235=D04(p1sq,p2sq,s34,s56,s12,s234,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p
     -   3+p2p4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C
     -   0235I,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235,D0
     -   1235R,D01235I,D1235R,D1235I)
c FC       If(PrintD) then
c FC       write(*,*)"D1235"
c FC       call print_DijR(D01235,D1235R,D1235I)
c FC       call tens_Int4_Re_Com(m,p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p4,
c FC     -   C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0235I,C
c FC     -   0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235R,D01235I,
c FC     -   D1235R,D1235I)
c FC       endif

       D01236=D04(p1sq,p2sq,s345,p6sq,s12,s16,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1
     -   p5,p2p3+p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C1
     -   26R,C123R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123
     -   I,D01236,D01236R,D01236I,D1236R,D1236I)
c FC       If(PrintD) then
c FC       write(*,*)"D1236"
c FC       call print_DijR(D01236,D1236R,D1236I)
c FC       call tens_Int4_Re_Com(m,p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p5,p2p3
c FC     -   +p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C126R,C12
c FC     -   3R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123I,D0123
c FC     -   6R,D01236I,D1236R,D1236I)
c FC       endif

       D01245=D04(p1sq,s23,p4sq,s56,s123,s234,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p
     -   4+p3p4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C
     -   0245I,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245,D0
     -   1245R,D01245I,D1245R,D1245I)
c FC       If(PrintD) then
c FC       write(*,*)"D1245"
c FC       call print_DijR(D01245,D1245R,D1245I)
c FC       call tens_Int4_Re_Com(m,p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p4,
c FC     -   C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0245I,C
c FC     -   0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245R,D01245I,
c FC     -   D1245R,D1245I)
c FC       endif

       D01246=D04(p1sq,s23,s45,p6sq,s123,s16,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5
     -   ,p2p4+p2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R
     -   ,C126R,C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C
     -   124I,D01246,D01246R,D01246I,D1246R,D1246I)
c FC       If(PrintD) then
c FC       write(*,*)"D1246"
c FC       call print_DijR(D01246,D1246R,D1246I)
c FC       call tens_Int4_Re_Com(m,p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,p2p4+p
c FC     -   2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R,C126R,
c FC     -   C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C124I,D0
c FC     -   1246R,D01246I,D1246R,D1246I)
c FC       endif

       D01256=D04(p1sq,s234,p5sq,p6sq,s56,s16,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1
     -   p5,p2p5+p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C1
     -   26R,C125R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125
     -   I,D01256,D01256R,D01256I,D1256R,D1256I)
c FC       If(PrintD) then
c FC       write(*,*)"D1256"
c FC       call print_DijR(D01256,D1256R,D1256I)
c FC       call tens_Int4_Re_Com(m,p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p5,p2p5
c FC     -   +p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C126R,C12
c FC     -   5R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125I,D0125
c FC     -   6R,D01256I,D1256R,D1256I)
c FC       endif

       D02345=D04(p2sq,p3sq,p4sq,s234,s23,s34,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0
     -   345R,C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C02
     -   45I,C0235I,C0234I,C345I,C245I,C235I,C234I,D02345,D02345R,D02
     -   345I,D2345R,D2345I)
c FC       If(PrintD) then
c FC       write(*,*)"D2345"
c FC       call print_DijR(D02345,D2345R,D2345I)
c FC       call tens_Int4_Re_Com(m,p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,C0
c FC     -   245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C0245I,C02
c FC     -   35I,C0234I,C345I,C245I,C235I,C234I,D02345R,D02345I,D2345R,D2
c FC     -   345I)
c FC       endif
c FC
       D02346=D04(p2sq,p3sq,s45,s16,s23,s345,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p
     -   4+p3p5,C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C
     -   0346I,C0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346,D0
     -   2346R,D02346I,D2346R,D2346I)
c FC       If(PrintD) then
c FC       write(*,*)"D2346"
c FC       call print_DijR(D02346,D2346R,D2346I)
c FC       call tens_Int4_Re_Com(m,p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4+p3p5,
c FC     -   C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C0346I,C
c FC     -   0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346R,D02346I,
c FC     -   D2346R,D2346I)
c FC       endif

       D02356=D04(p2sq,s34,p5sq,s16,s234,s345,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p
     -   5+p4p5,C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C
     -   0356I,C0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356,D0
     -   2356R,D02356I,D2356R,D2356I)
c FC       If(PrintD) then
c FC       write(*,*)"D2356"
c FC       call print_DijR(D02356,D2356R,D2356I)
c FC       call tens_Int4_Re_Com(m,p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5+p4p5,
c FC     -   C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C0356I,C
c FC     -   0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356R,D02356I,
c FC     -   D2356R,D2356I)
c FC       endif
c FC
       D03456=D04(p3sq,p4sq,p5sq,s345,s34,s45,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C0
     -   456R,C0356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C03
     -   56I,C0346I,C0345I,C456I,C356I,C346I,C345I,D03456,D03456R,D03
     -   456I,D3456R,D3456I)
c FC       If(PrintD) then
c FC       write(*,*)"D3456"
c FC       call print_DijR(D03456,D3456R,D3456I)
c FC       call tens_Int4_Re_Com(m,p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C0456R,C0
c FC     -   356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C0356I,C03
c FC     -   46I,C0345I,C456I,C356I,C346I,C345I,D03456R,D03456I,D3456R,D3
c FC     -   456I)
c FC       endif

       D01345=D04(s12,p3sq,p4sq,s56,s123,s34,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p
     -   4,p3p4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C
     -   0345I,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345,D0
     -   1345R,D01345I,D1345R,D1345I)
c FC       If(PrintD) then
c FC       write(*,*)"D1345"
c FC       call print_DijR(D01345,D1345R,D1345I)
c FC       call tens_Int4_Re_Com(m,s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p4,
c FC     -   C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0345I,C
c FC     -   0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345R,D01345I,
c FC     -   D1345R,D1345I)
c FC       endif

       D01346=D04(s12,p3sq,s45,p6sq,s123,s345,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5
     -   +p2p4+p2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R
     -   ,C136R,C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C
     -   134I,D01346,D01346R,D01346I,D1346R,D1346I)
c FC       If(PrintD) then
c FC       write(*,*)"D1346"
c FC       call print_DijR(D01346,D1346R,D1346I)
c FC       call tens_Int4_Re_Com(m,s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+p2p4+p
c FC     -   2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R,C136R,
c FC     -   C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C134I,D0
c FC     -   1346R,D01346I,D1346R,D1346I)
c FC       endif

       D01356=D04(s12,s34,p5sq,p6sq,s56,s345,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4
     -   ,p1p5+p2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R
     -   ,C136R,C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C
     -   135I,D01356,D01356R,D01356I,D1356R,D1356I)
c FC       If(PrintD) then
c FC       write(*,*)"D1356"
c FC       call print_DijR(D01356,D1356R,D1356I)
c FC       call tens_Int4_Re_Com(m,s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,p1p5+p
c FC     -   2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R,C136R,
c FC     -   C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C135I,D0
c FC     -   1356R,D01356I,D1356R,D1356I)
c FC       endif

       D01456=D04(s123,p4sq,p5sq,p6sq,s56,s45,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1
     -   p5+p2p5+p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C1
     -   46R,C145R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145
     -   I,D01456,D01456R,D01456I,D1456R,D1456I)
c FC       If(PrintD) then
c FC       write(*,*)"D1456"
c FC       call print_DijR(D01456,D1456R,D1456I)
c FC       call tens_Int4_Re_Com(m,s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p5+p2p5
c FC     -   +p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C146R,C14
c FC     -   5R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145I,D0145
c FC     -   6R,D01456I,D1456R,D1456I)
c FC       endif

       D02456=D04(s23,p4sq,p5sq,s16,s234,s45,m,m,m,m)
       call tens_red4_new_Re_Com_1M(m,s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p
     -   5,p4p5,C0456R,C0256R,C0246R,C0245R,C456R,C256R,C246R,C245R,C
     -   0456I,C0256I,C0246I,C0245I,C456I,C256I,C246I,C245I,D02456,D0
     -   2456R,D02456I,D2456R,D2456I)
c FC       If(PrintD) then
c FC       write(*,*)"D2456"
c FC       call print_DijR(D02456,D2456R,D2456I)
c FC       call tens_Int4_Re_Com(m,s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5,p4p5,
c FC     -   C0456R,C0256R,C0246R,C0245R,C456R,C256R,C246R,C245R,C0456I,C
c FC     -   0256I,C0246I,C0245I,C456I,C256I,C246I,C245I,D02456R,D02456I,
c FC     -   D2456R,D2456I)
c FC       endif

       E012345=E01M(m,p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2p3,p2p4,p3p
     -   4,D02345,D01345,D01245,D01235,D01234)
       E012345R=Dble(E012345)
       E012345I=DImag(E012345)
       call tens_red5_new_Re_Com_1M(m,p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1
     -   p4,p2p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2
     -   345R,D1345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01
     -   235I,D01234I,D2345I,D1345I,D1245I,D1235I,D1234I,E12345R,E123
     -   45I)
c FC       If(PrintE) then
c FC       write(*,*)"E12345"
c FC       call print_EijR(E012345,E12345R,E12345I)
c FC       call tens_Int5_Re_Com(m,p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2p3
c FC     -   ,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2345R,D1
c FC     -   345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01235I,D0
c FC     -   1234I,D2345I,D1345I,D1245I,D1235I,D1234I,E012345R,E012345I,E
c FC     -   12345R,E12345I)
c FC       endif

       E012346=E01M(m,p1sq,p2sq,p3sq,s45,p1p2,p1p3,p1p4+p1p5,p2p3,p2p4
     -   +p2p5,p3p4+p3p5,D02346,D01346,D01246,D01236,D01234)
       E012346R=Dble(E012346)
       E012346I=DImag(E012346)
       call tens_red5_new_Re_Com_1M(m,p1sq,p2sq,p3sq,s45,p1p2,p1p3,p1p
     -   4+p1p5,p2p3,p2p4+p2p5,p3p4+p3p5,D02346R,D01346R,D01246R,D012
     -   36R,D01234R,D2346R,D1346R,D1246R,D1236R,D1234R,D02346I,D0134
     -   6I,D01246I,D01236I,D01234I,D2346I,D1346I,D1246I,D1236I,D1234
     -   I,E12346R,E12346I)
c FC       If(PrintE) then
c FC       write(*,*)"E12346"
c FC       call print_EijR(E012346,E12346R,E12346I)
c FC       call tens_Int5_Re_Com(m,p1sq,p2sq,p3sq,s45,p1p2,p1p3,p1p4+p1p5,
c FC     -   p2p3,p2p4+p2p5,p3p4+p3p5,D02346R,D01346R,D01246R,D01236R,D01
c FC     -   234R,D2346R,D1346R,D1246R,D1236R,D1234R,D02346I,D01346I,D012
c FC     -   46I,D01236I,D01234I,D2346I,D1346I,D1246I,D1236I,D1234I,E0123
c FC     -   46R,E012346I,E12346R,E12346I)
c FC       endif

       E012356=E01M(m,p1sq,p2sq,s34,p5sq,p1p2,p1p3+p1p4,p1p5,p2p3+p2p4
     -   ,p2p5,p3p5+p4p5,D02356,D01356,D01256,D01236,D01235)
       E012356R=Dble(E012356)
       E012356I=DImag(E012356)
       call tens_red5_new_Re_Com_1M(m,p1sq,p2sq,s34,p5sq,p1p2,p1p3+p1p
     -   4,p1p5,p2p3+p2p4,p2p5,p3p5+p4p5,D02356R,D01356R,D01256R,D012
     -   36R,D01235R,D2356R,D1356R,D1256R,D1236R,D1235R,D02356I,D0135
     -   6I,D01256I,D01236I,D01235I,D2356I,D1356I,D1256I,D1236I,D1235
     -   I,E12356R,E12356I)
c FC       If(PrintE) then
c FC       write(*,*)"E12356"
c FC       call print_EijR(E012356,E12356R,E12356I)
c FC       call tens_Int5_Re_Com(m,p1sq,p2sq,s34,p5sq,p1p2,p1p3+p1p4,p1p5,
c FC     -   p2p3+p2p4,p2p5,p3p5+p4p5,D02356R,D01356R,D01256R,D01236R,D01
c FC     -   235R,D2356R,D1356R,D1256R,D1236R,D1235R,D02356I,D01356I,D012
c FC     -   56I,D01236I,D01235I,D2356I,D1356I,D1256I,D1236I,D1235I,E0123
c FC     -   56R,E012356I,E12356R,E12356I)
c FC       endif

       E012456=E01M(m,p1sq,s23,p4sq,p5sq,p1p2+p1p3,p1p4,p1p5,p2p4+p3p4
     -   ,p2p5+p3p5,p4p5,D02456,D01456,D01256,D01246,D01245)
       E012456R=Dble(E012456)
       E012456I=DImag(E012456)
       call tens_red5_new_Re_Com_1M(m,p1sq,s23,p4sq,p5sq,p1p2+p1p3,p1p
     -   4,p1p5,p2p4+p3p4,p2p5+p3p5,p4p5,D02456R,D01456R,D01256R,D012
     -   46R,D01245R,D2456R,D1456R,D1256R,D1246R,D1245R,D02456I,D0145
     -   6I,D01256I,D01246I,D01245I,D2456I,D1456I,D1256I,D1246I,D1245
     -   I,E12456R,E12456I)
c FC       If(PrintE) then
c FC       write(*,*)"E12456"
c FC       call print_EijR(E012456,E12456R,E12456I)
c FC       call tens_Int5_Re_Com(m,p1sq,s23,p4sq,p5sq,p1p2+p1p3,p1p4,p1p5,
c FC     -   p2p4+p3p4,p2p5+p3p5,p4p5,D02456R,D01456R,D01256R,D01246R,D01
c FC     -   245R,D2456R,D1456R,D1256R,D1246R,D1245R,D02456I,D01456I,D012
c FC     -   56I,D01246I,D01245I,D2456I,D1456I,D1256I,D1246I,D1245I,E0124
c FC     -   56R,E012456I,E12456R,E12456I)
c FC       endif

       E023456=E01M(m,p2sq,p3sq,p4sq,p5sq,p2p3,p2p4,p2p5,p3p4,p3p5,p4p
     -   5,D03456,D02456,D02356,D02346,D02345)
       E023456R=Dble(E023456)
       E023456I=DImag(E023456)
       call tens_red5_new_Re_Com_1M(m,p2sq,p3sq,p4sq,p5sq,p2p3,p2p4,p2
     -   p5,p3p4,p3p5,p4p5,D03456R,D02456R,D02356R,D02346R,D02345R,D3
     -   456R,D2456R,D2356R,D2346R,D2345R,D03456I,D02456I,D02356I,D02
     -   346I,D02345I,D3456I,D2456I,D2356I,D2346I,D2345I,E23456R,E234
     -   56I)
c FC       If(PrintE) then
c FC       write(*,*)"E23456"
c FC       call print_EijR(E023456,E23456R,E23456I)
c FC       call tens_Int5_Re_Com(m,p2sq,p3sq,p4sq,p5sq,p2p3,p2p4,p2p5,p3p4
c FC     -   ,p3p5,p4p5,D03456R,D02456R,D02356R,D02346R,D02345R,D3456R,D2
c FC     -   456R,D2356R,D2346R,D2345R,D03456I,D02456I,D02356I,D02346I,D0
c FC     -   2345I,D3456I,D2456I,D2356I,D2346I,D2345I,E023456R,E023456I,E
c FC     -   23456R,E23456I)
c FC       endif

       E013456=E01M(m,s12,p3sq,p4sq,p5sq,p1p3+p2p3,p1p4+p2p4,p1p5+p2p5
     -   ,p3p4,p3p5,p4p5,D03456,D01456,D01356,D01346,D01345)
       E013456R=Dble(E013456)
       E013456I=DImag(E013456)
       call tens_red5_new_Re_Com_1M(m,s12,p3sq,p4sq,p5sq,p1p3+p2p3,p1p
     -   4+p2p4,p1p5+p2p5,p3p4,p3p5,p4p5,D03456R,D01456R,D01356R,D013
     -   46R,D01345R,D3456R,D1456R,D1356R,D1346R,D1345R,D03456I,D0145
     -   6I,D01356I,D01346I,D01345I,D3456I,D1456I,D1356I,D1346I,D1345
     -   I,E13456R,E13456I)
c FC       If(PrintE) then
c FC       write(*,*)"E13456"
c FC       call print_EijR(E013456,E13456R,E13456I)
c FC       call tens_Int5_Re_Com(m,s12,p3sq,p4sq,p5sq,p1p3+p2p3,p1p4+p2p4,
c FC     -   p1p5+p2p5,p3p4,p3p5,p4p5,D03456R,D01456R,D01356R,D01346R,D01
c FC     -   345R,D3456R,D1456R,D1356R,D1346R,D1345R,D03456I,D01456I,D013
c FC     -   56I,D01346I,D01345I,D3456I,D1456I,D1356I,D1346I,D1345I,E0134
c FC     -   56R,E013456I,E13456R,E13456I)
c FC       endif

       F0123456=F0fin1M(m,p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,s12,s23,s34,s4
     -   5,s56,s16,s123,s234,s345,E023456,E013456,E012456,E012356,E01
     -   2346,E012345)
       F0123456R=Dble(F0123456)
       F0123456I=DImag(F0123456)
       call tens_red6_new_Re_Com(p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1
     -   p4,p1p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,E023456R,E013456R,E012
     -   456R,E012356R,E012346R,E012345R,E23456R,E13456R,E12456R,E123
     -   56R,E12346R,E12345R,E023456I,E013456I,E012456I,E012356I,E012
     -   346I,E012345I,E23456I,E13456I,E12456I,E12356I,E12346I,E12345
     -   I,F123456R,F123456I,F123456)
c FC       If(PrintF) then
c FC       write(*,*)"F123456"
c FC       call print_FijR(F0123456,F123456R,F123456I)
c FC       call tens_I6_R3(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1p
c FC     -   5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D3456
c FC     -   I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E234
c FC     -   56R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E013456
c FC     -   I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E12456
c FC     -   I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F1234
c FC     -   56I)
c FC       call tens_I6_R41(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R42(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R51(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R52(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R53(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R54(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R55(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R56(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R57(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R58(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R59(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p1
c FC     -   p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D345
c FC     -   6I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E23
c FC     -   456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E01345
c FC     -   6I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E1245
c FC     -   6I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F123
c FC     -   456I)
c FC       call tens_I6_R510(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p
c FC     -   1p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D34
c FC     -   56I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E2
c FC     -   3456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E0134
c FC     -   56I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E124
c FC     -   56I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F12
c FC     -   3456I)
c FC       call tens_I6_R511(msq,p1sq,p2sq,p3sq,p4sq,p5sq,p1p2,p1p3,p1p4,p
c FC     -   1p5,p2p3,p2p4,p2p5,p3p4,p3p5,p4p5,D03456R,D3456R,D03456I,D34
c FC     -   56I,E023456R,E013456R,E012456R,E012356R,E012346R,E012345R,E2
c FC     -   3456R,E13456R,E12456R,E12356R,E12346R,E12345R,E023456I,E0134
c FC     -   56I,E012456I,E012356I,E012346I,E012345I,E23456I,E13456I,E124
c FC     -   56I,E12356I,E12346I,E12345I,F0123456R,F0123456I,F123456R,F12
c FC     -   3456I)
c FC       Endif

c************************************************************************************
c************************************************************************************
c************************************************************************************
 
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call FFHexEven1(F(1))
       call FFHexEven2(F(765))
       call FFHexEven3(F(1529))
       call FFHexEven4(F(2293))
       call FFHexEven5(F(3057))


      endif  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup1 = dotrc(p1,mup1)
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrc(p1,mup3)
       p1mup4 = dotrc(p1,mup4)
       p1mup5 = dotrc(p1,mup5)
       p2mup1 = dotrc(p2,mup1)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup4 = dotrc(p2,mup4)
       p2mup5 = dotrc(p2,mup5)
       p3mup1 = dotrc(p3,mup1)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup4 = dotrc(p3,mup4)
       p3mup5 = dotrc(p3,mup5)
       p4mup1 = dotrc(p4,mup1)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup4 = dotrc(p4,mup4)
       p4mup5 = dotrc(p4,mup5)
       p5mup1 = dotrc(p5,mup1)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrc(p5,mup3)
       p5mup4 = dotrc(p5,mup4)
       p5mup5 = dotrc(p5,mup5)
       mup1mup2 = dotcc(mup1,mup2)
       mup1mup3 = dotcc(mup1,mup3)
       mup1mup4 = dotcc(mup1,mup4)
       mup1mup5 = dotcc(mup1,mup5)
       mup2mup3 = dotcc(mup2,mup3)
       mup2mup4 = dotcc(mup2,mup4)
       mup2mup5 = dotcc(mup2,mup5)
       mup3mup4 = dotcc(mup3,mup4)
       mup3mup5 = dotcc(mup3,mup5)
       mup4mup5 = dotcc(mup4,mup5)

cfc         print*,"p1mup1 ",       p1mup1 
cfc         print*,"p1mup2 ",       p1mup2 
cfc         print*,"p1mup3 ",       p1mup3 
cfc         print*,"p1mup4 ",       p1mup4 
cfc         print*,"p1mup5 ",       p1mup5 
cfc         print*,"p2mup1 ",       p2mup1 
cfc         print*,"p2mup2 ",       p2mup2 
cfc         print*,"p2mup3 ",       p2mup3 
cfc         print*,"p2mup4 ",       p2mup4 
cfc         print*,"p2mup5 ",       p2mup5 
cfc         print*,"p3mup1 ",       p3mup1 
cfc         print*,"p3mup2 ",       p3mup2 
cfc         print*,"p3mup3 ",       p3mup3 
cfc         print*,"p3mup4 ",       p3mup4 
cfc         print*,"p3mup5 ",       p3mup5 
cfc         print*,"p4mup1 ",       p4mup1 
cfc         print*,"p4mup2 ",       p4mup2 
cfc         print*,"p4mup3 ",       p4mup3 
cfc         print*,"p4mup4 ",       p4mup4 
cfc         print*,"p4mup5 ",       p4mup5 
cfc         print*,"p5mup1 ",       p5mup1 
cfc         print*,"p5mup2 ",       p5mup2 
cfc         print*,"p5mup3 ",       p5mup3 
cfc         print*,"p5mup4 ",       p5mup4 
cfc         print*,"p5mup5 ",       p5mup5 
cfc         print*,"mup1mup2 ",       mup1mup2
cfc         print*,"mup1mup3 ",       mup1mup3
cfc         print*,"mup1mup4 ",       mup1mup4
cfc         print*,"mup1mup5 ",       mup1mup5
cfc         print*,"mup2mup3 ",       mup2mup3
cfc         print*,"mup2mup4 ",       mup2mup4
cfc         print*,"mup2mup5 ",       mup2mup5
cfc         print*,"mup3mup4 ",       mup3mup4
cfc         print*,"mup3mup5 ",       mup3mup5
cfc         print*,"mup4mup5 ",       mup4mup5



c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call FaFunctiongggg(p1mup1,p1mup2,p1mup3,p1mup4,p1mup5,p2mup1,p
     -   2mup2,p2mup3,p2mup4,p2mup5,p3mup1,p3mup2,p3mup3,p3mup4,p3mup
     -   5,p4mup1,p4mup2,p4mup3,p4mup4,p4mup5,p5mup1,p5mup2,p5mup3,p5
     -   mup4,p5mup5,mup1mup2,mup1mup3,mup1mup4,mup1mup5,mup2mup3,mup
     -   2mup4,mup2mup5,mup3mup4,mup3mup5,mup4mup5,Fa)
c************************************************************************************
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


       result = msq*(Fa(2)+Fa(3)+Fa(4)+Fa(5)+Fa(6)+Fa(15)+Fa(16)+Fa(17
     -   )+Fa(19)+Fa(20)+2*(Fa(1)+Fa(7)+Fa(14)+Fa(21))+Fa(22)-4*(Fa(1
     -   3)-Fa(18)+Fa(36)+Fa(37))-16*(Fa(10)+Fa(11)+Fa(12)+Fa(23)+Fa(
     -   27)-Fa(29)+Fa(30)+Fa(40)+Fa(44)+Fa(45)+Fa(46)+Fa(48)+Fa(49))
     -   -8*(Fa(8)+Fa(9)+Fa(24)+Fa(25)+Fa(26)+Fa(28)+Fa(31)+Fa(32)+Fa
     -   (33)+Fa(34)+Fa(35)+Fa(38)-Fa(39)+Fa(41)+Fa(42)+Fa(43)+Fa(47)
     -   +Fa(50)))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
