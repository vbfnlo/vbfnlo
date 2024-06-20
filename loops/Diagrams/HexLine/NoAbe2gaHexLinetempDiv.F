       subroutine NoAbe2gaHexLineDiv(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p
     -   1,mup2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result
     -   _graph,result_color,Div)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 27/10/2008
c Modified:5/6/2012
c ************************************************************************************
c determine the non-abelian finite QCD virtual corrections along the quark line i.e 
c The sum of all virtual corrections with color factor CA of the sort showed
c below.Attention!!!!This class of virtual corrections correspond to three 
c different type of born diagram.The p5 leg can be attacht in whaterver
c position among p2 and p3:(p5,p2,p3),(p2,p5,p3),(p2,p3,p5). So,
c to get the full result you should call this subroutine plus
c the three abelian penline subroutines for this combination of momenta(pay
c attention that the conventions for the names are different for this one
c and penlinemmm or penline). To finish, you have the cross terms
c  p2<->p3, so you have to call again the non-abelian part and three abelian penline
c  subroutines that goes with this configuration.(Pay attention to the position of the gluon)  
c                          mu_p6                       
c                            p6                      
c                            $                     
c                            $                     
c                            $                     
c                            V                      
c                   $$$$$$$$$$$$$$$$$$$$$              
c                   $                   $              
c                   $                   $               
c                   $                   $               
c psi(p1)   ------>------->------->------->-----   bar_psi(p5)
c                      $      $      $                
c                      $      $      $                
c                      $      $      $               
c                      V      V      V               
c                      $      $      $              
c                      $      $      $               
c                     p2      p3     p4              
c                   mu_p2   mu_p3   mu_p4           
c Note: The vertices are just Gamma^(mu_p2)..terms. So the correct
c factor should be added by hand for each boson. i.e, Without adding
c anything it represents pp->gluon+2 photons up to the electromagnetic charge 
c (for uu->g+2 photons is missing (2/3)^2) from the coupling. For W and Z,
c we have to account for the Diracgamma_5 that are missing in the vertex
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
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       IMPLICIT NONE
      Complex*16 SMB(63),Fa(872),F(9706)

      Complex*16 barpsi_p4(2),psi_p1(2)
       Real*8 P(9687)
       Complex*16  SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5cccrr,SC5ccccr,SC7ccccrrr
       EXTERNAL    SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr, SC5ccrrr,SC5cccrr,SC5ccccr,SC7ccccrrr 
       Integer alpha,comp
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
       Complex*16   mup2(0:3),mup3(0:3),mup5(0:3),mup6(0:3)
       Complex*16   p1mup2, p1mup3, p1mup5, p1mup6, p2mup2, p2mup3, 
     -          p2mup5, p2mup6, p3mup2, p3mup3, p3mup5, p3mup6, 
     -          p4mup2, p4mup3, p4mup5, p4mup6, p5mup2, p5mup3, 
     -          p5mup5, p5mup6, p6mup2, p6mup3, p6mup5, p6mup6
       Complex*16   mup2mup3, mup2mup5, mup2mup6, mup3mup5, mup3mup6, 
     -          mup5mup6
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
      Complex*16 dotrc,dotcc,result_graph(6),result_color(2)
       Real*8 musq
      EXTERNAL   dotrc,dotcc
       COMMON/NoAbe2gaHexLineFaFunctions/Fa
       COMMON/NoAbe2gaHexLineFhlFunctions/F
      Save/NoAbe2gaHexLineFhlFunctions/
       COMMON/NoAbe2gaHexLinePFunctions/P
      Save/NoAbe2gaHexLinePFunctions/
       COMMON/NoAbe2gaHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6sq,
     -   s12,s23,s34,s45,s56,s16,s123,s234,s345
       COMMON/NoAbe2gaHexLineFVALFunctions/F0123456R,F123456R,F0123456
     -   I,F123456I
       COMMON/NoAbe2gaHexLineEVALFunctions/ E012345R,E12345R,E012345I,
     -   E12345I, E012346R,E12346R,E012346I,E12346I, E012356R,E12356R
     -   ,E012356I,E12356I, E012456R,E12456R,E012456I,E12456I, E01345
     -   6R,E13456R,E013456I,E13456I, E023456R,E23456R,E023456I,E2345
     -   6I
       COMMON/NoAbe2gaHexLineDVALFunctions/ D01234R,D1234R,D01234I,D12
     -   34I, D01235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I,D
     -   1236I, D01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D01246I
     -   ,D1246I, D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D0134
     -   5I,D1345I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D01
     -   356I,D1356I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R,D
     -   02345I,D2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D2356R
     -   ,D02356I,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D345
     -   6R,D03456I,D3456I
       COMMON/NoAbe2gaHexLineCVALFunctions/ C0123R,C123R,C0123I,C123I,
     -    C0124R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C0126
     -   R,C126R,C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C135
     -   R,C0135I,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C014
     -   5I,C145I, C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C15
     -   6I, C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, C0
     -   236R,C236R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R,C
     -   246R,C0246I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R,C
     -   0345I,C345I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356I,
     -   C356I, C0456R,C456R,C0456I,C456I
       COMMON/NoAbe2gaHexLineBVALFunctions/ B012R,B012I, B013R,B013I, 
     -   B014R,B014I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,B0
     -   24I, B025R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B03
     -   6R,B036I, B045R,B045I, B046R,B046I, B056R,B056I
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
     -   3,C0123R,C0123I,C123R,C123I)
       C0126=C0finDiv(p1sq,s16,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,s16,p6sq,B026,B016,B012,C012
     -   6,C0126R,C0126I,C126R,C126I)
       C0124=C0finDiv(p1sq,s23,s123,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,s23,s123,B024,B014,B012,C012
     -   4,C0124R,C0124I,C124R,C124I)
       C0125=C0finDiv(p1sq,s234,s56,musq,Div)
       call tens_red3_new_Re_Com_Div(p1sq,s234,s56,B025,B015,B012,C012
     -   5,C0125R,C0125I,C125R,C125I)
       C0234=C0finDiv(p2sq,p3sq,s23,musq,Div)
       call tens_red3_new_Re_Com_Div(p2sq,p3sq,s23,B034,B024,B023,C023
     -   4,C0234R,C0234I,C234R,C234I)
       C0235=C0finDiv(p2sq,s34,s234,musq,Div)
       call tens_red3_new_Re_Com_Div(p2sq,s34,s234,B035,B025,B023,C023
     -   5,C0235R,C0235I,C235R,C235I)
       C0236=C0finDiv(p2sq,s345,s16,musq,Div)
       call tens_red3_new_Re_Com_Div(p2sq,s345,s16,B036,B026,B023,C023
     -   6,C0236R,C0236I,C236R,C236I)
       C0345=C0finDiv(p3sq,p4sq,s34,musq,Div)
       call tens_red3_new_Re_Com_Div(p3sq,p4sq,s34,B045,B035,B034,C034
     -   5,C0345R,C0345I,C345R,C345I)
       C0346=C0finDiv(p3sq,s45,s345,musq,Div)
       call tens_red3_new_Re_Com_Div(p3sq,s45,s345,B046,B036,B034,C034
     -   6,C0346R,C0346I,C346R,C346I)
       C0456=C0finDiv(p4sq,p5sq,s45,musq,Div)
       call tens_red3_new_Re_Com_Div(p4sq,p5sq,s45,B056,B046,B045,C045
     -   6,C0456R,C0456I,C456R,C456I)
       C0134=C0finDiv(s12,p3sq,s123,musq,Div)
       call tens_red3_new_Re_Com_Div(s12,p3sq,s123,B034,B014,B013,C013
     -   4,C0134R,C0134I,C134R,C134I)
       C0135=C0finDiv(s12,s34,s56,musq,Div)
       call tens_red3_new_Re_Com_Div(s12,s34,s56,B035,B015,B013,C0135,
     -   C0135R,C0135I,C135R,C135I)
       C0136=C0finDiv(s12,s345,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(s12,s345,p6sq,B036,B016,B013,C013
     -   6,C0136R,C0136I,C136R,C136I)
       C0145=C0finDiv(s123,p4sq,s56,musq,Div)
       call tens_red3_new_Re_Com_Div(s123,p4sq,s56,B045,B015,B014,C014
     -   5,C0145R,C0145I,C145R,C145I)
       C0146=C0finDiv(s123,s45,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(s123,s45,p6sq,B046,B016,B014,C014
     -   6,C0146R,C0146I,C146R,C146I)
       C0245=C0finDiv(s23,p4sq,s234,musq,Div)
       call tens_red3_new_Re_Com_Div(s23,p4sq,s234,B045,B025,B024,C024
     -   5,C0245R,C0245I,C245R,C245I)
       C0246=C0finDiv(s23,s45,s16,musq,Div)
       call tens_red3_new_Re_Com_Div(s23,s45,s16,B046,B026,B024,C0246,
     -   C0246R,C0246I,C246R,C246I)
       C0256=C0finDiv(s234,p5sq,s16,musq,Div)
       call tens_red3_new_Re_Com_Div(s234,p5sq,s16,B056,B026,B025,C025
     -   6,C0256R,C0256I,C256R,C256I)
       C0356=C0finDiv(s34,p5sq,s345,musq,Div)
       call tens_red3_new_Re_Com_Div(s34,p5sq,s345,B056,B036,B035,C035
     -   6,C0356R,C0356I,C356R,C356I)
       C0156=C0finDiv(s56,p5sq,p6sq,musq,Div)
       call tens_red3_new_Re_Com_Div(s56,p5sq,p6sq,B056,B016,B015,C015
     -   6,C0156R,C0156I,C156R,C156I)
       D01234=D0finDiv(s12,s23,p1sq,p2sq,p3sq,s123,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C02
     -   34R,C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C013
     -   4I,C0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D012
     -   34I,D1234R,D1234I)
       D01235=D0finDiv(s12,s234,p1sq,p2sq,s34,s56,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3
     -   +p2p4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0
     -   235I,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235,D01
     -   235R,D01235I,D1235R,D1235I)
       D01236=D0finDiv(s12,s16,p1sq,p2sq,s345,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p
     -   5,p2p3+p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C12
     -   6R,C123R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123I
     -   ,D01236,D01236R,D01236I,D1236R,D1236I)
       D01245=D0finDiv(s123,s234,p1sq,s23,p4sq,s56,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4
     -   +p3p4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0
     -   245I,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245,D01
     -   245R,D01245I,D1245R,D1245I)
       D01246=D0finDiv(s123,s16,p1sq,s23,s45,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,
     -   p2p4+p2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R,
     -   C126R,C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C1
     -   24I,D01246,D01246R,D01246I,D1246R,D1246I)
       D01256=D0finDiv(s56,s16,p1sq,s234,p5sq,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p
     -   5,p2p5+p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C12
     -   6R,C125R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125I
     -   ,D01256,D01256R,D01256I,D1256R,D1256I)
       D02345=D0finDiv(s23,s34,p2sq,p3sq,p4sq,s234,musq,Div)
       call tens_red4_new_Re_Com_Div(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C03
     -   45R,C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C024
     -   5I,C0235I,C0234I,C345I,C245I,C235I,C234I,D02345,D02345R,D023
     -   45I,D2345R,D2345I)
       D02346=D0finDiv(s23,s345,p2sq,p3sq,s45,s16,musq,Div)
       call tens_red4_new_Re_Com_Div(p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4
     -   +p3p5,C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C0
     -   346I,C0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346,D02
     -   346R,D02346I,D2346R,D2346I)
       D02356=D0finDiv(s234,s345,p2sq,s34,p5sq,s16,musq,Div)
       call tens_red4_new_Re_Com_Div(p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5
     -   +p4p5,C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C0
     -   356I,C0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356,D02
     -   356R,D02356I,D2356R,D2356I)
       D03456=D0finDiv(s34,s45,p3sq,p4sq,p5sq,s345,musq,Div)
       call tens_red4_new_Re_Com_Div(p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C04
     -   56R,C0356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C035
     -   6I,C0346I,C0345I,C456I,C356I,C346I,C345I,D03456,D03456R,D034
     -   56I,D3456R,D3456I)
       D01345=D0finDiv(s123,s34,s12,p3sq,p4sq,s56,musq,Div)
       call tens_red4_new_Re_Com_Div(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4
     -   ,p3p4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0
     -   345I,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345,D01
     -   345R,D01345I,D1345R,D1345I)
       D01346=D0finDiv(s123,s345,s12,p3sq,s45,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+
     -   p2p4+p2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R,
     -   C136R,C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C1
     -   34I,D01346,D01346R,D01346I,D1346R,D1346I)
       D01356=D0finDiv(s56,s345,s12,s34,p5sq,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,
     -   p1p5+p2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R,
     -   C136R,C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C1
     -   35I,D01356,D01356R,D01356I,D1356R,D1356I)
       D01456=D0finDiv(s56,s45,s123,p4sq,p5sq,p6sq,musq,Div)
       call tens_red4_new_Re_Com_Div(s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p
     -   5+p2p5+p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C14
     -   6R,C145R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145I
     -   ,D01456,D01456R,D01456I,D1456R,D1456I)
       D02456=D0finDiv(s234,s45,s23,p4sq,p5sq,s16,musq,Div)
       call tens_red4_new_Re_Com_Div(s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5
     -   ,p4p5,C0456R,C0256R,C0246R,C0245R,C456R,C256R,C246R,C245R,C0
     -   456I,C0256I,C0246I,C0245I,C456I,C256I,C246I,C245I,D02456,D02
     -   456R,D02456I,D2456R,D2456I)
       E012345=E0fin(p1sq,p2sq,p3sq,p4sq,s56,s12,s23,s34,s123,s234,D02
     -   345,D01345,D01245,D01235,D01234)
       E012345R=Dble(E012345)
       E012345I=DImag(E012345)
       call tens_red5_new_Re_Com(p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
     -   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,D2345R,
     -   D1345R,D1245R,D1235R,D1234R,D02345I,D01345I,D01245I,D01235I,
     -   D01234I,D2345I,D1345I,D1245I,D1235I,D1234I,E12345R,E12345I)
       E012346=E0fin(p1sq,p2sq,p3sq,s45,p6sq,s12,s23,s345,s123,s16,D02
     -   346,D01346,D01246,D01236,D01234)
       E012346R=Dble(E012346)
       E012346I=DImag(E012346)
       call tens_red5_new_Re_Com(p1sq,p2sq,p3sq,s45,p1p2,p1p3,p1p4+p1p
     -   5,p2p3,p2p4+p2p5,p3p4+p3p5,D02346R,D01346R,D01246R,D01236R,D
     -   01234R,D2346R,D1346R,D1246R,D1236R,D1234R,D02346I,D01346I,D0
     -   1246I,D01236I,D01234I,D2346I,D1346I,D1246I,D1236I,D1234I,E12
     -   346R,E12346I)
       E012356=E0fin(p1sq,p2sq,s34,p5sq,p6sq,s12,s234,s345,s56,s16,D02
     -   356,D01356,D01256,D01236,D01235)
       E012356R=Dble(E012356)
       E012356I=DImag(E012356)
       call tens_red5_new_Re_Com(p1sq,p2sq,s34,p5sq,p1p2,p1p3+p1p4,p1p
     -   5,p2p3+p2p4,p2p5,p3p5+p4p5,D02356R,D01356R,D01256R,D01236R,D
     -   01235R,D2356R,D1356R,D1256R,D1236R,D1235R,D02356I,D01356I,D0
     -   1256I,D01236I,D01235I,D2356I,D1356I,D1256I,D1236I,D1235I,E12
     -   356R,E12356I)
       E012456=E0fin(p1sq,s23,p4sq,p5sq,p6sq,s123,s234,s45,s56,s16,D02
     -   456,D01456,D01256,D01246,D01245)
       E012456R=Dble(E012456)
       E012456I=DImag(E012456)
       call tens_red5_new_Re_Com(p1sq,s23,p4sq,p5sq,p1p2+p1p3,p1p4,p1p
     -   5,p2p4+p3p4,p2p5+p3p5,p4p5,D02456R,D01456R,D01256R,D01246R,D
     -   01245R,D2456R,D1456R,D1256R,D1246R,D1245R,D02456I,D01456I,D0
     -   1256I,D01246I,D01245I,D2456I,D1456I,D1256I,D1246I,D1245I,E12
     -   456R,E12456I)
       E023456=E0fin(p2sq,p3sq,p4sq,p5sq,s16,s23,s34,s45,s234,s345,D03
     -   456,D02456,D02356,D02346,D02345)
       E023456R=Dble(E023456)
       E023456I=DImag(E023456)
       call tens_red5_new_Re_Com(p2sq,p3sq,p4sq,p5sq,p2p3,p2p4,p2p5,p3
     -   p4,p3p5,p4p5,D03456R,D02456R,D02356R,D02346R,D02345R,D3456R,
     -   D2456R,D2356R,D2346R,D2345R,D03456I,D02456I,D02356I,D02346I,
     -   D02345I,D3456I,D2456I,D2356I,D2346I,D2345I,E23456R,E23456I)
       E013456=E0fin(s12,p3sq,p4sq,p5sq,p6sq,s123,s34,s45,s56,s345,D03
     -   456,D01456,D01356,D01346,D01345)
       E013456R=Dble(E013456)
       E013456I=DImag(E013456)
       call tens_red5_new_Re_Com(s12,p3sq,p4sq,p5sq,p1p3+p2p3,p1p4+p2p
     -   4,p1p5+p2p5,p3p4,p3p5,p4p5,D03456R,D01456R,D01356R,D01346R,D
     -   01345R,D3456R,D1456R,D1356R,D1346R,D1345R,D03456I,D01456I,D0
     -   1356I,D01346I,D01345I,D3456I,D1456I,D1356I,D1346I,D1345I,E13
     -   456R,E13456I)
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
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call NoAbe2gaHexLineFFhl1(F(1))
       call NoAbe2gaHexLineFFhl2(F(195))
       call NoAbe2gaHexLineFFhl3(F(389))
       call NoAbe2gaHexLineFFhl4(F(583))
       call NoAbe2gaHexLineFFhl5(F(777))
       call NoAbe2gaHexLineFFhl6(F(971))
       call NoAbe2gaHexLineFFhl7(F(1165))
       call NoAbe2gaHexLineFFhl8(F(1359))
       call NoAbe2gaHexLineFFhl9(F(1553))
       call NoAbe2gaHexLineFFhl10(F(1747))
       call NoAbe2gaHexLineFFhl11(F(1941))
       call NoAbe2gaHexLineFFhl12(F(2135))
       call NoAbe2gaHexLineFFhl13(F(2329))
       call NoAbe2gaHexLineFFhl14(F(2523))
       call NoAbe2gaHexLineFFhl15(F(2717))
       call NoAbe2gaHexLineFFhl16(F(2911))
       call NoAbe2gaHexLineFFhl17(F(3105))
       call NoAbe2gaHexLineFFhl18(F(3299))
       call NoAbe2gaHexLineFFhl19(F(3493))
       call NoAbe2gaHexLineFFhl20(F(3687))
       call NoAbe2gaHexLineFFhl21(F(3881))
       call NoAbe2gaHexLineFFhl22(F(4075))
       call NoAbe2gaHexLineFFhl23(F(4269))
       call NoAbe2gaHexLineFFhl24(F(4463))
       call NoAbe2gaHexLineFFhl25(F(4657))
       call NoAbe2gaHexLineFFhl26(F(4851))
       call NoAbe2gaHexLineFFhl27(F(5045))
       call NoAbe2gaHexLineFFhl28(F(5239))
       call NoAbe2gaHexLineFFhl29(F(5433))
       call NoAbe2gaHexLineFFhl30(F(5627))
       call NoAbe2gaHexLineFFhl31(F(5821))
       call NoAbe2gaHexLineFFhl32(F(6015))
       call NoAbe2gaHexLineFFhl33(F(6209))
       call NoAbe2gaHexLineFFhl34(F(6403))
       call NoAbe2gaHexLineFFhl35(F(6597))
       call NoAbe2gaHexLineFFhl36(F(6791))
       call NoAbe2gaHexLineFFhl37(F(6985))
       call NoAbe2gaHexLineFFhl38(F(7179))
       call NoAbe2gaHexLineFFhl39(F(7373))
       call NoAbe2gaHexLineFFhl40(F(7567))
       call NoAbe2gaHexLineFFhl41(F(7761))
       call NoAbe2gaHexLineFFhl42(F(7955))
       call NoAbe2gaHexLineFFhl43(F(8149))
       call NoAbe2gaHexLineFFhl44(F(8343))
       call NoAbe2gaHexLineFFhl45(F(8537))
       call NoAbe2gaHexLineFFhl46(F(8731))
       call NoAbe2gaHexLineFFhl47(F(8925))
       call NoAbe2gaHexLineFFhl48(F(9119))
       call NoAbe2gaHexLineFFhl49(F(9313))
       call NoAbe2gaHexLineFFhl50(F(9507))
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
       p1mup5 = dotrc(p1,mup5)
       p1mup6 = dotrc(p1,mup6)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup5 = dotrc(p2,mup5)
       p2mup6 = dotrc(p2,mup6)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup5 = dotrc(p3,mup5)
       p3mup6 = dotrc(p3,mup6)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup5 = dotrc(p4,mup5)
       p4mup6 = dotrc(p4,mup6)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrc(p5,mup3)
       p5mup5 = dotrc(p5,mup5)
       p5mup6 = dotrc(p5,mup6)
       p6mup2 = dotrc(p6,mup2)
       p6mup3 = dotrc(p6,mup3)
       p6mup5 = dotrc(p6,mup5)
       p6mup6 = dotrc(p6,mup6)
       mup2mup3 = dotcc(mup2,mup3)
       mup2mup5 = dotcc(mup2,mup5)
       mup2mup6 = dotcc(mup2,mup6)
       mup3mup5 = dotcc(mup3,mup5)
       mup3mup6 = dotcc(mup3,mup6)
       mup5mup6 = dotcc(mup5,mup6)
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call NoAbe2gaHexLineFa1(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2mu
     -   p3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,p
     -   4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mup
     -   5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5mu
     -   p6,Fa(1))
       call NoAbe2gaHexLineFa2(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2mu
     -   p3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,p
     -   4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mup
     -   5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5mu
     -   p6,Fa(25))
       call NoAbe2gaHexLineFa3(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2mu
     -   p3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,p
     -   4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mup
     -   5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5mu
     -   p6,Fa(49))
       call NoAbe2gaHexLineFa4(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2mu
     -   p3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,p
     -   4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mup
     -   5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5mu
     -   p6,Fa(73))
       call NoAbe2gaHexLineFa5(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2mu
     -   p3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,p
     -   4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mup
     -   5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5mu
     -   p6,Fa(97))
       call NoAbe2gaHexLineFa6(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2mu
     -   p3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,p
     -   4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mup
     -   5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5mu
     -   p6,Fa(121))
       call NoAbe2gaHexLineFa7(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2mu
     -   p3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,p
     -   4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mup
     -   5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5mu
     -   p6,Fa(145))
       call NoAbe2gaHexLineFa8(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2mu
     -   p3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,p
     -   4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mup
     -   5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5mu
     -   p6,Fa(169))
       call NoAbe2gaHexLineFa9(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2mu
     -   p3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,p
     -   4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mup
     -   5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5mu
     -   p6,Fa(193))
       call NoAbe2gaHexLineFa10(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(217))
       call NoAbe2gaHexLineFa11(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(241))
       call NoAbe2gaHexLineFa12(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(265))
       call NoAbe2gaHexLineFa13(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(289))
       call NoAbe2gaHexLineFa14(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(313))
       call NoAbe2gaHexLineFa15(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(337))
       call NoAbe2gaHexLineFa16(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(361))
       call NoAbe2gaHexLineFa17(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(385))
       call NoAbe2gaHexLineFa18(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(409))
       call NoAbe2gaHexLineFa19(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(433))
       call NoAbe2gaHexLineFa20(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(457))
       call NoAbe2gaHexLineFa21(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(481))
       call NoAbe2gaHexLineFa22(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(505))
       call NoAbe2gaHexLineFa23(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(529))
       call NoAbe2gaHexLineFa24(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(553))
       call NoAbe2gaHexLineFa25(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(577))
       call NoAbe2gaHexLineFa26(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(601))
       call NoAbe2gaHexLineFa27(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(625))
       call NoAbe2gaHexLineFa28(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(649))
       call NoAbe2gaHexLineFa29(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(673))
       call NoAbe2gaHexLineFa30(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(697))
       call NoAbe2gaHexLineFa31(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(721))
       call NoAbe2gaHexLineFa32(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(745))
       call NoAbe2gaHexLineFa33(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(769))
       call NoAbe2gaHexLineFa34(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(793))
       call NoAbe2gaHexLineFa35(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2m
     -   up3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3,
     -   p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6mu
     -   p5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5m
     -   up6,Fa(817))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the Matrix Element  
c************************************************************************************
c************************************************************************************
c************************************************************************************
       SMB(1) = SC1c(barpsi_p4,mup2,psi_p1,alpha)
       SMB(2) = SC1c(barpsi_p4,mup5,psi_p1,alpha)
       SMB(3) = SC3ccc(barpsi_p4,mup6,mup5,mup2,psi_p1,alpha)
       SMB(4) = SC1c(barpsi_p4,mup3,psi_p1,alpha)
       SMB(5) = SC3ccc(barpsi_p4,mup6,mup5,mup3,psi_p1,alpha)
       SMB(6) = SC1r(barpsi_p4,p2,psi_p1,alpha)
       SMB(7) = SC3ccr(barpsi_p4,mup6,mup5,p2,psi_p1,alpha)
       SMB(8) = SC1c(barpsi_p4,mup6,psi_p1,alpha)
       SMB(9) = SC1r(barpsi_p4,p5,psi_p1,alpha)
       SMB(10) = SC1r(barpsi_p4,p6,psi_p1,alpha)
       SMB(11) = SC3ccr(barpsi_p4,mup6,mup5,p5,psi_p1,alpha)
       SMB(12) = SC3ccr(barpsi_p4,mup6,mup5,p6,psi_p1,alpha)
       SMB(13) = SC3ccr(barpsi_p4,mup3,mup2,p2,psi_p1,alpha)
       SMB(14) = SC3ccr(barpsi_p4,mup5,mup2,p2,psi_p1,alpha)
       SMB(15) = SC5ccccr(barpsi_p4,mup6,mup5,mup3,mup2,p2,psi_p1,alpha)
       SMB(16) = SC3crr(barpsi_p4,mup2,p5,p2,psi_p1,alpha)
       SMB(17) = SC3crr(barpsi_p4,mup2,p6,p2,psi_p1,alpha)
       SMB(18) = SC3ccr(barpsi_p4,mup6,mup2,p2,psi_p1,alpha)
       SMB(19) = SC5cccrr(barpsi_p4,mup6,mup5,mup2,p5,p2,psi_p1,alpha)
       SMB(20) = SC5cccrr(barpsi_p4,mup6,mup5,mup2,p6,p2,psi_p1,alpha)
       SMB(21) = SC3ccr(barpsi_p4,mup3,mup2,p5,psi_p1,alpha)
       SMB(22) = SC3ccr(barpsi_p4,mup3,mup2,p6,psi_p1,alpha)
       SMB(23) = SC3ccr(barpsi_p4,mup5,mup2,p5,psi_p1,alpha)
       SMB(24) = SC3ccr(barpsi_p4,mup5,mup2,p6,psi_p1,alpha)
       SMB(25) = SC3ccc(barpsi_p4,mup5,mup3,mup2,psi_p1,alpha)
       SMB(26) = SC3ccr(barpsi_p4,mup5,mup3,p2,psi_p1,alpha)
       SMB(27) = SC3ccr(barpsi_p4,mup5,mup3,p5,psi_p1,alpha)
       SMB(28) = SC3ccr(barpsi_p4,mup5,mup3,p6,psi_p1,alpha)
       SMB(29) = SC3ccr(barpsi_p4,mup6,mup2,p5,psi_p1,alpha)
       SMB(30) = SC3ccr(barpsi_p4,mup6,mup2,p6,psi_p1,alpha)
       SMB(31) = SC3ccc(barpsi_p4,mup6,mup3,mup2,psi_p1,alpha)
       SMB(32) = SC5ccccr(barpsi_p4,mup6,mup5,mup3,mup2,p5,psi_p1,alpha)
       SMB(33) = SC5ccccr(barpsi_p4,mup6,mup5,mup3,mup2,p6,psi_p1,alpha)
       SMB(34) = SC3crr(barpsi_p4,mup3,p5,p2,psi_p1,alpha)
       SMB(35) = SC3crr(barpsi_p4,mup3,p6,p2,psi_p1,alpha)
       SMB(36) = SC3crr(barpsi_p4,mup5,p5,p2,psi_p1,alpha)
       SMB(37) = SC3crr(barpsi_p4,mup5,p6,p2,psi_p1,alpha)
       SMB(38) = SC3ccr(barpsi_p4,mup6,mup3,p2,psi_p1,alpha)
       SMB(39) = SC3crr(barpsi_p4,mup6,p5,p2,psi_p1,alpha)
       SMB(40) = SC3crr(barpsi_p4,mup6,p6,p2,psi_p1,alpha)
       SMB(41) = SC5cccrr(barpsi_p4,mup6,mup5,mup3,p5,p2,psi_p1,alpha)
       SMB(42) = SC5cccrr(barpsi_p4,mup6,mup5,mup3,p6,p2,psi_p1,alpha)
       SMB(43) = SC3ccr(barpsi_p4,mup6,mup3,p5,psi_p1,alpha)
       SMB(44) = SC3ccr(barpsi_p4,mup6,mup3,p6,psi_p1,alpha)
       SMB(45) = SC3crr(barpsi_p4,mup3,p6,p5,psi_p1,alpha)
       SMB(46) = SC3crr(barpsi_p4,mup5,p6,p5,psi_p1,alpha)
       SMB(47) = SC3crr(barpsi_p4,mup6,p6,p5,psi_p1,alpha)
       SMB(48) = SC5cccrr(barpsi_p4,mup6,mup5,mup3,p6,p5,psi_p1,alpha)
       SMB(49) = SC5ccrrr(barpsi_p4,mup3,mup2,p6,p5,p2,psi_p1,alpha)
       SMB(50) = SC5cccrr(barpsi_p4,mup6,mup3,mup2,p5,p2,psi_p1,alpha)
       SMB(51) = SC5cccrr(barpsi_p4,mup6,mup3,mup2,p6,p2,psi_p1,alpha)
       SMB(52) = SC5cccrr(barpsi_p4,mup5,mup3,mup2,p5,p2,psi_p1,alpha)
       SMB(53) = SC5cccrr(barpsi_p4,mup5,mup3,mup2,p6,p2,psi_p1,alpha)
       SMB(54) = SC7ccccrrr(barpsi_p4,mup6,mup5,mup3,mup2,p6,p5,p2,psi_p1,alpha)
       SMB(55) = SC5ccrrr(barpsi_p4,mup5,mup2,p6,p5,p2,psi_p1,alpha)
       SMB(56) = SC5ccrrr(barpsi_p4,mup6,mup2,p6,p5,p2,psi_p1,alpha)
       SMB(57) = SC3crr(barpsi_p4,mup2,p6,p5,psi_p1,alpha)
       SMB(58) = SC5cccrr(barpsi_p4,mup6,mup5,mup2,p6,p5,psi_p1,alpha)
       SMB(59) = SC3rrr(barpsi_p4,p6,p5,p2,psi_p1,alpha)
       SMB(60) = SC5cccrr(barpsi_p4,mup5,mup3,mup2,p6,p5,psi_p1,alpha)
       SMB(61) = SC5ccrrr(barpsi_p4,mup5,mup3,p6,p5,p2,psi_p1,alpha)
       SMB(62) = SC5cccrr(barpsi_p4,mup6,mup3,mup2,p6,p5,psi_p1,alpha)
       SMB(63) = SC5ccrrr(barpsi_p4,mup6,mup3,p6,p5,p2,psi_p1,alpha)
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
       result_graph(1) = (Fa(1)*SMB(1)+Fa(827)*SMB(4)+Fa(15)*SMB(8)+Fa
     -   (828)*SMB(9)+Fa(829)*SMB(10)+Fa(56)*SMB(13)+F(177)*SMB(15)+F
     -   a(57)*SMB(16)+Fa(58)*SMB(17)+Fa(59)*SMB(21)+Fa(60)*SMB(22)+F
     -   a(61)*SMB(25)-Fa(62)*SMB(26)+Fa(63)*SMB(27)+Fa(64)*SMB(28)+F
     -   a(65)*SMB(31)+F(269)*SMB(32)+F(270)*SMB(33)+Fa(830)*SMB(34)+
     -   Fa(831)*SMB(35)+Fa(88)*SMB(43)+Fa(89)*SMB(44)+Fa(832)*SMB(45
     -   )+Fa(100)*SMB(46)-4*(Fa(2)*SMB(2)-Fa(14)*SMB(6)+Fa(55)*(F(15
     -   8)*SMB(11)+F(159)*SMB(12))-Fa(84)*(F(158)*SMB(36)+F(159)*SMB
     -   (37))+F(158)*Fa(86)*SMB(39)+F(159)*Fa(86)*SMB(40)-Fa(101)*SM
     -   B(47))+Fa(103)*SMB(49)+Fa(104)*SMB(50)+Fa(105)*SMB(51)+Fa(10
     -   6)*SMB(52)+Fa(107)*SMB(53)+F(351)*SMB(54)+2*(Fa(13)*SMB(5)+F
     -   a(85)*SMB(38)+F(158)*Fa(87)*SMB(41)+F(159)*Fa(87)*SMB(42)+F(
     -   335)*Fa(102)*SMB(48)+mup3mup6*(-(F(84)*SMB(14))+F(96)*SMB(23
     -   )+F(132)*SMB(24)+F(335)*SMB(55))+mup3mup5*(F(84)*SMB(18)+F(1
     -   05)*SMB(29)+F(141)*SMB(30)-F(335)*SMB(56)))+Fa(108)*SMB(57)+
     -   Fa(109)*SMB(59)+Fa(110)*SMB(60)+Fa(111)*SMB(61)+Fa(112)*SMB(
     -   62)+Fa(113)*SMB(63))/(s234*s34)
       result_graph(2) = (Fa(115)*SMB(2)+Fa(833)*SMB(4)+Fa(126)*SMB(5)
     -   +Fa(127)*SMB(6)+4*Fa(128)*SMB(8)+Fa(129)*SMB(9)+Fa(130)*SMB(
     -   10)+Fa(131)*SMB(13)+F(475)*SMB(15)+Fa(132)*SMB(16)+Fa(133)*S
     -   MB(17)+Fa(137)*SMB(26)+Fa(138)*SMB(27)+Fa(139)*SMB(28)+Fa(14
     -   2)*SMB(35)+Fa(143)*SMB(38)+Fa(144)*SMB(43)+Fa(145)*SMB(44)+F
     -   a(150)*SMB(49)+Fa(151)*SMB(50)+Fa(152)*SMB(51)+Fa(153)*SMB(5
     -   2)+Fa(154)*SMB(53)+F(562)*SMB(54)+2*(mup3mup6*F(443)*SMB(14)
     -   +mup3mup5*F(372)*SMB(18)+Fa(141)*SMB(34)+Fa(146)*SMB(45)-F(5
     -   35)*(p6mup2*SMB(41)-p5mup2*SMB(42)+2*(Fa(147)*SMB(46)-Fa(148
     -   )*SMB(47))+Fa(149)*SMB(48)+mup3mup6*(2*(p6mup2*SMB(36)-p5mup
     -   2*SMB(37))-SMB(55))+mup3mup5*(-2*(p6mup2*SMB(39)-p5mup2*SMB(
     -   40))+SMB(56))))+s12*(Fa(114)*SMB(1)-Fa(134)*SMB(21)-Fa(135)*
     -   SMB(22)+mup3mup6*F(498)*SMB(23)+Fa(136)*SMB(25)+mup3mup5*F(4
     -   26)*SMB(29)+Fa(140)*SMB(31)+F(514)*SMB(32)+F(515)*SMB(33)+2*
     -   (mup3mup6*F(445)*SMB(24)+mup3mup5*F(374)*SMB(30)-Fa(155)*SMB
     -   (57))+Fa(156)*SMB(60)+Fa(157)*SMB(62)))/(s12*s34)
       result_graph(3) = (Fa(834)*SMB(4)+Fa(835)*SMB(6)+4*(Fa(159)*SMB
     -   (2)+Fa(183)*SMB(7)-Fa(184)*SMB(8))+Fa(836)*SMB(9)+Fa(837)*SM
     -   B(10)+Fa(212)*SMB(13)+Fa(213)*SMB(14)+F(683)*SMB(15)+Fa(838)
     -   *SMB(16)+Fa(839)*SMB(17)+Fa(235)*SMB(18)+s12*(Fa(158)*SMB(1)
     -   +Fa(237)*SMB(22)+2*(Fa(236)*SMB(21)-F(676)*Fa(238)*SMB(25)))
     -   +Fa(239)*SMB(26)+Fa(240)*SMB(27)+Fa(241)*SMB(28)+8*s12*Fa(24
     -   2)*F(9573)*SMB(31)+Fa(243)*SMB(34)-Fa(244)*SMB(35)+Fa(245)*S
     -   MB(38)-2*F(652)*(2*Fa(211)*SMB(12)+p5mup3*SMB(20)+p6mup2*SMB
     -   (41))+Fa(246)*SMB(43)+Fa(247)*SMB(44)+Fa(248)*SMB(45)-Fa(149
     -   )*(F(629)*SMB(5)+2*F(652)*SMB(48))-Fa(249)*SMB(49)+Fa(250)*S
     -   MB(50)+Fa(251)*SMB(51)+Fa(252)*SMB(52)+Fa(253)*SMB(53)+F(651
     -   )*(4*Fa(210)*SMB(11)+2*(p6mup3*SMB(19)+p5mup2*SMB(42))+SMB(5
     -   4)))/(s12*s123)
       result_graph(4) = (Fa(254)*SMB(1)+Fa(255)*SMB(2)+Fa(840)*SMB(4)
     -   +Fa(266)*SMB(5)+Fa(267)*SMB(6)+Fa(268)*SMB(8)+Fa(269)*SMB(9)
     -   +Fa(270)*SMB(10)+Fa(271)*SMB(13)+F(1079)*SMB(15)+Fa(272)*SMB
     -   (16)+Fa(273)*SMB(17)+Fa(274)*SMB(21)+Fa(275)*SMB(22)+Fa(276)
     -   *SMB(25)+Fa(277)*SMB(26)+Fa(278)*SMB(27)+Fa(279)*SMB(28)+Fa(
     -   280)*SMB(31)+F(1230)*SMB(32)+F(1231)*SMB(33)+Fa(841)*SMB(34)
     -   +Fa(842)*SMB(35)+Fa(301)*SMB(36)+Fa(303)*SMB(38)+Fa(306)*SMB
     -   (41)+Fa(308)*SMB(43)+Fa(309)*SMB(44)+Fa(843)*SMB(45)+4*(Fa(3
     -   02)*SMB(37)+Fa(304)*SMB(39)+Fa(305)*SMB(40)+Fa(320)*SMB(46)+
     -   Fa(321)*SMB(47))+Fa(323)*SMB(49)+Fa(324)*SMB(50)+Fa(325)*SMB
     -   (51)+Fa(326)*SMB(52)+Fa(327)*SMB(53)+F(1723)*SMB(54)+2*(Fa(3
     -   07)*SMB(42)+Fa(322)*SMB(48)+mup3mup6*(-(F(1078)*SMB(14))+F(1
     -   152)*SMB(23)+F(1153)*SMB(24)+F(1724)*SMB(55))+mup3mup5*(F(11
     -   00)*SMB(18)+F(1223)*SMB(29)+F(1224)*SMB(30)+F(1725)*SMB(56))
     -   )+Fa(328)*SMB(57)+Fa(329)*SMB(59)+Fa(330)*SMB(60)+Fa(331)*SM
     -   B(61)+Fa(332)*SMB(62)+Fa(333)*SMB(63))/s34
       result_graph(5) = (Fa(845)*SMB(2)+Fa(846)*SMB(4)+Fa(369)*SMB(5)
     -   +Fa(847)*SMB(6)+Fa(394)*SMB(7)+Fa(848)*SMB(8)+Fa(849)*SMB(9)
     -   +Fa(850)*SMB(10)-Fa(450)*SMB(12)+Fa(451)*SMB(13)+Fa(452)*SMB
     -   (14)+F(2458)*SMB(15)+Fa(851)*SMB(16)+Fa(852)*SMB(17)+Fa(485)
     -   *SMB(18)+Fa(486)*SMB(19)+Fa(487)*SMB(20)+Fa(493)*SMB(26)+Fa(
     -   494)*SMB(27)+Fa(495)*SMB(28)+Fa(499)*SMB(34)+Fa(500)*SMB(35)
     -   +Fa(501)*SMB(36)+Fa(502)*SMB(37)+Fa(503)*SMB(38)+Fa(504)*SMB
     -   (39)+4*(Fa(449)*SMB(11)-Fa(505)*SMB(40))+Fa(506)*SMB(43)+Fa(
     -   507)*SMB(44)+Fa(508)*SMB(45)+Fa(509)*SMB(46)+Fa(510)*SMB(47)
     -   -2*F(2647)*(p6mup2*SMB(41)-p5mup2*SMB(42)+Fa(149)*SMB(48))+F
     -   a(511)*SMB(49)+Fa(512)*SMB(50)+Fa(513)*SMB(51)+Fa(514)*SMB(5
     -   2)+Fa(515)*SMB(53)+F(2688)*SMB(54)+Fa(516)*SMB(55)+Fa(517)*S
     -   MB(56)+s12*(Fa(844)*SMB(1)+Fa(358)*SMB(3)+Fa(488)*SMB(21)+Fa
     -   (489)*SMB(22)+Fa(490)*SMB(23)+Fa(491)*SMB(24)+Fa(492)*SMB(25
     -   )+Fa(496)*SMB(29)+Fa(497)*SMB(30)+Fa(498)*SMB(31)+F(2579)*SM
     -   B(32)+F(2580)*SMB(33)+Fa(518)*SMB(57)+Fa(519)*SMB(60)+Fa(520
     -   )*SMB(62)))/s12
       result_graph(6) = Fa(853)*SMB(1)+Fa(854)*SMB(2)+Fa(547)*SMB(3)+
     -   Fa(855)*SMB(4)+Fa(558)*SMB(5)+Fa(856)*SMB(6)+Fa(584)*SMB(7)+
     -   Fa(857)*SMB(8)+Fa(858)*SMB(9)+Fa(859)*SMB(10)+Fa(648)*SMB(12
     -   )+Fa(649)*SMB(13)+Fa(650)*SMB(14)+F(7370)*SMB(15)+Fa(860)*SM
     -   B(16)+Fa(861)*SMB(17)+Fa(683)*SMB(18)+Fa(684)*SMB(19)+Fa(685
     -   )*SMB(20)+Fa(686)*SMB(21)+Fa(687)*SMB(22)+Fa(688)*SMB(23)+Fa
     -   (689)*SMB(24)+Fa(690)*SMB(25)+Fa(691)*SMB(26)+Fa(692)*SMB(27
     -   )+Fa(693)*SMB(28)+Fa(694)*SMB(29)+Fa(695)*SMB(30)+Fa(696)*SM
     -   B(31)+F(8137)*SMB(32)+F(8138)*SMB(33)+Fa(862)*SMB(34)+Fa(863
     -   )*SMB(35)+Fa(864)*SMB(36)+Fa(865)*SMB(37)+Fa(733)*SMB(38)+4*
     -   (Fa(647)*SMB(11)+Fa(866)*SMB(39)+Fa(867)*SMB(40))+Fa(756)*SM
     -   B(41)+Fa(758)*SMB(43)+Fa(759)*SMB(44)+Fa(868)*SMB(45)+Fa(869
     -   )*SMB(46)+Fa(870)*SMB(47)+2*(Fa(757)*SMB(42)-Fa(788)*SMB(48)
     -   )+Fa(789)*SMB(49)+Fa(790)*SMB(50)+Fa(791)*SMB(51)+Fa(792)*SM
     -   B(52)+Fa(793)*SMB(53)+F(9432)*SMB(54)+Fa(794)*SMB(55)+Fa(795
     -   )*SMB(56)+Fa(871)*SMB(57)+Fa(872)*SMB(59)+Fa(823)*SMB(60)+Fa
     -   (824)*SMB(61)+Fa(825)*SMB(62)+Fa(826)*SMB(63)
       If (ngluon.eq.2) then
       result_color(1)=-(result_graph(1)+result_graph(2)+result_graph(
     -   3)+result_graph(4)+result_graph(5)+result_graph(6))
       elseif (ngluon.eq.3) then
       If (posgluon.eq.2) then
       result_color(1)=-(result_graph(4)+result_graph(6))
       result_color(2)=-(result_graph(1)+result_graph(2)+result_graph(
     -   3)+result_graph(5))
       elseIf (posgluon.eq.3) then
       result_color(1)=-(result_graph(5)+result_graph(6))
       result_color(2)=-(result_graph(1)+result_graph(2)+result_graph(
     -   3)+result_graph(4))
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     -   . Look to the heading for explanation" 
       endif
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     -   . Look to the heading for explanation" 
       endif
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
