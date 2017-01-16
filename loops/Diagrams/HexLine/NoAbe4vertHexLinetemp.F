       subroutine NoAbe4vertHexLine(p1,p2,p3,p4,p5,p6,barpsi_p4,psi_p1
     -   ,mup2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_
     -   graph,result_color)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 27/10/2008
c Modified:28/6/2012
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
      Complex*16 SMB(32),Fa(199),F(140)

      Complex*16 barpsi_p4(2),psi_p1(2)
       Real*8 P(39)
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
      Complex*16 dotrc,dotcc,result_graph(6,3),result_color(3)
       Real*8 musq
      EXTERNAL   dotrc,dotcc
       COMMON/NoAbe4vertHexLineFaFunctions/Fa
       COMMON/NoAbe4vertHexLineFhlFunctions/F
      Save/NoAbe4vertHexLineFhlFunctions/
       COMMON/NoAbe4vertHexLinePFunctions/P
      Save/NoAbe4vertHexLinePFunctions/
       COMMON/NoAbe4vertHexLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,p6s
     -   q,s12,s23,s34,s45,s56,s16,s123,s234,s345
       COMMON/NoAbe4vertHexLineFVALFunctions/F0123456R,F123456R,F01234
     -   56I,F123456I
       COMMON/NoAbe4vertHexLineEVALFunctions/ E012345R,E12345R,E012345
     -   I,E12345I, E012346R,E12346R,E012346I,E12346I, E012356R,E1235
     -   6R,E012356I,E12356I, E012456R,E12456R,E012456I,E12456I, E013
     -   456R,E13456R,E013456I,E13456I, E023456R,E23456R,E023456I,E23
     -   456I
       COMMON/NoAbe4vertHexLineDVALFunctions/ D01234R,D1234R,D01234I,D
     -   1234I, D01235R,D1235R,D01235I,D1235I, D01236R,D1236R,D01236I
     -   ,D1236I, D01245R,D1245R,D01245I,D1245I, D01246R,D1246R,D0124
     -   6I,D1246I, D01256R,D1256R,D01256I,D1256I, D01345R,D1345R,D01
     -   345I,D1345I, D01346R,D1346R,D01346I,D1346I, D01356R,D1356R,D
     -   01356I,D1356I, D01456R,D1456R,D01456I,D1456I, D02345R,D2345R
     -   ,D02345I,D2345I, D02346R,D2346R,D02346I,D2346I, D02356R,D235
     -   6R,D02356I,D2356I, D02456R,D2456R,D02456I,D2456I, D03456R,D3
     -   456R,D03456I,D3456I
       COMMON/NoAbe4vertHexLineCVALFunctions/ C0123R,C123R,C0123I,C123
     -   I, C0124R,C124R,C0124I,C124I, C0125R,C125R,C0125I,C125I, C01
     -   26R,C126R,C0126I,C126I, C0134R,C134R,C0134I,C134I, C0135R,C1
     -   35R,C0135I,C135I, C0136R,C136R,C0136I,C136I, C0145R,C145R,C0
     -   145I,C145I, C0146R,C146R,C0146I,C146I, C0156R,C156R,C0156I,C
     -   156I, C0234R,C234R,C0234I,C234I, C0235R,C235R,C0235I,C235I, 
     -   C0236R,C236R,C0236I,C236I, C0245R,C245R,C0245I,C245I, C0246R
     -   ,C246R,C0246I,C246I, C0256R,C256R,C0256I,C256I, C0345R,C345R
     -   ,C0345I,C345I, C0346R,C346R,C0346I,C346I, C0356R,C356R,C0356
     -   I,C356I, C0456R,C456R,C0456I,C456I
       COMMON/NoAbe4vertHexLineBVALFunctions/ B012R,B012I, B013R,B013I
     -   , B014R,B014I, B015R,B015I, B016R,B016I, B023R,B023I, B024R,
     -   B024I, B025R,B025I, B026R,B026I, B034R,B034I, B035R,B035I, B
     -   036R,B036I, B045R,B045I, B046R,B046I, B056R,B056I
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
     -   123R,C0123I,C123R,C123I)
       C0126=C0fin(p1sq,s16,p6sq,musq)
       call tens_red3_new_Re_Com(p1sq,s16,p6sq,B026,B016,B012,C0126,C0
     -   126R,C0126I,C126R,C126I)
       C0124=C0fin(p1sq,s23,s123,musq)
       call tens_red3_new_Re_Com(p1sq,s23,s123,B024,B014,B012,C0124,C0
     -   124R,C0124I,C124R,C124I)
       C0125=C0fin(p1sq,s234,s56,musq)
       call tens_red3_new_Re_Com(p1sq,s234,s56,B025,B015,B012,C0125,C0
     -   125R,C0125I,C125R,C125I)
       C0234=C0fin(p2sq,p3sq,s23,musq)
       call tens_red3_new_Re_Com(p2sq,p3sq,s23,B034,B024,B023,C0234,C0
     -   234R,C0234I,C234R,C234I)
       C0235=C0fin(p2sq,s34,s234,musq)
       call tens_red3_new_Re_Com(p2sq,s34,s234,B035,B025,B023,C0235,C0
     -   235R,C0235I,C235R,C235I)
       C0236=C0fin(p2sq,s345,s16,musq)
       call tens_red3_new_Re_Com(p2sq,s345,s16,B036,B026,B023,C0236,C0
     -   236R,C0236I,C236R,C236I)
       C0345=C0fin(p3sq,p4sq,s34,musq)
       call tens_red3_new_Re_Com(p3sq,p4sq,s34,B045,B035,B034,C0345,C0
     -   345R,C0345I,C345R,C345I)
       C0346=C0fin(p3sq,s45,s345,musq)
       call tens_red3_new_Re_Com(p3sq,s45,s345,B046,B036,B034,C0346,C0
     -   346R,C0346I,C346R,C346I)
       C0456=C0fin(p4sq,p5sq,s45,musq)
       call tens_red3_new_Re_Com(p4sq,p5sq,s45,B056,B046,B045,C0456,C0
     -   456R,C0456I,C456R,C456I)
       C0134=C0fin(s12,p3sq,s123,musq)
       call tens_red3_new_Re_Com(s12,p3sq,s123,B034,B014,B013,C0134,C0
     -   134R,C0134I,C134R,C134I)
       C0135=C0fin(s12,s34,s56,musq)
       call tens_red3_new_Re_Com(s12,s34,s56,B035,B015,B013,C0135,C013
     -   5R,C0135I,C135R,C135I)
       C0136=C0fin(s12,s345,p6sq,musq)
       call tens_red3_new_Re_Com(s12,s345,p6sq,B036,B016,B013,C0136,C0
     -   136R,C0136I,C136R,C136I)
       C0145=C0fin(s123,p4sq,s56,musq)
       call tens_red3_new_Re_Com(s123,p4sq,s56,B045,B015,B014,C0145,C0
     -   145R,C0145I,C145R,C145I)
       C0146=C0fin(s123,s45,p6sq,musq)
       call tens_red3_new_Re_Com(s123,s45,p6sq,B046,B016,B014,C0146,C0
     -   146R,C0146I,C146R,C146I)
       C0245=C0fin(s23,p4sq,s234,musq)
       call tens_red3_new_Re_Com(s23,p4sq,s234,B045,B025,B024,C0245,C0
     -   245R,C0245I,C245R,C245I)
       C0246=C0fin(s23,s45,s16,musq)
       call tens_red3_new_Re_Com(s23,s45,s16,B046,B026,B024,C0246,C024
     -   6R,C0246I,C246R,C246I)
       C0256=C0fin(s234,p5sq,s16,musq)
       call tens_red3_new_Re_Com(s234,p5sq,s16,B056,B026,B025,C0256,C0
     -   256R,C0256I,C256R,C256I)
       C0356=C0fin(s34,p5sq,s345,musq)
       call tens_red3_new_Re_Com(s34,p5sq,s345,B056,B036,B035,C0356,C0
     -   356R,C0356I,C356R,C356I)
       C0156=C0fin(s56,p5sq,p6sq,musq)
       call tens_red3_new_Re_Com(s56,p5sq,p6sq,B056,B016,B015,C0156,C0
     -   156R,C0156I,C156R,C156I)
       D01234=D0fin(s12,s23,p1sq,p2sq,p3sq,s123,musq)
       call tens_red4_new_Re_Com(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
     -   C0134R,C0124R,C0123R,C234R,C134R,C124R,C123R,C0234I,C0134I,C
     -   0124I,C0123I,C234I,C134I,C124I,C123I,D01234,D01234R,D01234I,
     -   D1234R,D1234I)
       D01235=D0fin(s12,s234,p1sq,p2sq,s34,s56,musq)
       call tens_red4_new_Re_Com(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
     -   4,C0235R,C0135R,C0125R,C0123R,C235R,C135R,C125R,C123R,C0235I
     -   ,C0135I,C0125I,C0123I,C235I,C135I,C125I,C123I,D01235,D01235R
     -   ,D01235I,D1235R,D1235I)
       D01236=D0fin(s12,s16,p1sq,p2sq,s345,p6sq,musq)
       call tens_red4_new_Re_Com(p1sq,p2sq,s345,p1p2,p1p3+p1p4+p1p5,p2
     -   p3+p2p4+p2p5,C0236R,C0136R,C0126R,C0123R,C236R,C136R,C126R,C
     -   123R,C0236I,C0136I,C0126I,C0123I,C236I,C136I,C126I,C123I,D01
     -   236,D01236R,D01236I,D1236R,D1236I)
       D01245=D0fin(s123,s234,p1sq,s23,p4sq,s56,musq)
       call tens_red4_new_Re_Com(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
     -   4,C0245R,C0145R,C0125R,C0124R,C245R,C145R,C125R,C124R,C0245I
     -   ,C0145I,C0125I,C0124I,C245I,C145I,C125I,C124I,D01245,D01245R
     -   ,D01245I,D1245R,D1245I)
       D01246=D0fin(s123,s16,p1sq,s23,s45,p6sq,musq)
       call tens_red4_new_Re_Com(p1sq,s23,s45,p1p2+p1p3,p1p4+p1p5,p2p4
     -   +p2p5+p3p4+p3p5,C0246R,C0146R,C0126R,C0124R,C246R,C146R,C126
     -   R,C124R,C0246I,C0146I,C0126I,C0124I,C246I,C146I,C126I,C124I,
     -   D01246,D01246R,D01246I,D1246R,D1246I)
       D01256=D0fin(s56,s16,p1sq,s234,p5sq,p6sq,musq)
       call tens_red4_new_Re_Com(p1sq,s234,p5sq,p1p2+p1p3+p1p4,p1p5,p2
     -   p5+p3p5+p4p5,C0256R,C0156R,C0126R,C0125R,C256R,C156R,C126R,C
     -   125R,C0256I,C0156I,C0126I,C0125I,C256I,C156I,C126I,C125I,D01
     -   256,D01256R,D01256I,D1256R,D1256I)
       D02345=D0fin(s23,s34,p2sq,p3sq,p4sq,s234,musq)
       call tens_red4_new_Re_Com(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
     -   C0245R,C0235R,C0234R,C345R,C245R,C235R,C234R,C0345I,C0245I,C
     -   0235I,C0234I,C345I,C245I,C235I,C234I,D02345,D02345R,D02345I,
     -   D2345R,D2345I)
       D02346=D0fin(s23,s345,p2sq,p3sq,s45,s16,musq)
       call tens_red4_new_Re_Com(p2sq,p3sq,s45,p2p3,p2p4+p2p5,p3p4+p3p
     -   5,C0346R,C0246R,C0236R,C0234R,C346R,C246R,C236R,C234R,C0346I
     -   ,C0246I,C0236I,C0234I,C346I,C246I,C236I,C234I,D02346,D02346R
     -   ,D02346I,D2346R,D2346I)
       D02356=D0fin(s234,s345,p2sq,s34,p5sq,s16,musq)
       call tens_red4_new_Re_Com(p2sq,s34,p5sq,p2p3+p2p4,p2p5,p3p5+p4p
     -   5,C0356R,C0256R,C0236R,C0235R,C356R,C256R,C236R,C235R,C0356I
     -   ,C0256I,C0236I,C0235I,C356I,C256I,C236I,C235I,D02356,D02356R
     -   ,D02356I,D2356R,D2356I)
       D03456=D0fin(s34,s45,p3sq,p4sq,p5sq,s345,musq)
       call tens_red4_new_Re_Com(p3sq,p4sq,p5sq,p3p4,p3p5,p4p5,C0456R,
     -   C0356R,C0346R,C0345R,C456R,C356R,C346R,C345R,C0456I,C0356I,C
     -   0346I,C0345I,C456I,C356I,C346I,C345I,D03456,D03456R,D03456I,
     -   D3456R,D3456I)
       D01345=D0fin(s123,s34,s12,p3sq,p4sq,s56,musq)
       call tens_red4_new_Re_Com(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
     -   4,C0345R,C0145R,C0135R,C0134R,C345R,C145R,C135R,C134R,C0345I
     -   ,C0145I,C0135I,C0134I,C345I,C145I,C135I,C134I,D01345,D01345R
     -   ,D01345I,D1345R,D1345I)
       D01346=D0fin(s123,s345,s12,p3sq,s45,p6sq,musq)
       call tens_red4_new_Re_Com(s12,p3sq,s45,p1p3+p2p3,p1p4+p1p5+p2p4
     -   +p2p5,p3p4+p3p5,C0346R,C0146R,C0136R,C0134R,C346R,C146R,C136
     -   R,C134R,C0346I,C0146I,C0136I,C0134I,C346I,C146I,C136I,C134I,
     -   D01346,D01346R,D01346I,D1346R,D1346I)
       D01356=D0fin(s56,s345,s12,s34,p5sq,p6sq,musq)
       call tens_red4_new_Re_Com(s12,s34,p5sq,p1p3+p1p4+p2p3+p2p4,p1p5
     -   +p2p5,p3p5+p4p5,C0356R,C0156R,C0136R,C0135R,C356R,C156R,C136
     -   R,C135R,C0356I,C0156I,C0136I,C0135I,C356I,C156I,C136I,C135I,
     -   D01356,D01356R,D01356I,D1356R,D1356I)
       D01456=D0fin(s56,s45,s123,p4sq,p5sq,p6sq,musq)
       call tens_red4_new_Re_Com(s123,p4sq,p5sq,p1p4+p2p4+p3p4,p1p5+p2
     -   p5+p3p5,p4p5,C0456R,C0156R,C0146R,C0145R,C456R,C156R,C146R,C
     -   145R,C0456I,C0156I,C0146I,C0145I,C456I,C156I,C146I,C145I,D01
     -   456,D01456R,D01456I,D1456R,D1456I)
       D02456=D0fin(s234,s45,s23,p4sq,p5sq,s16,musq)
       call tens_red4_new_Re_Com(s23,p4sq,p5sq,p2p4+p3p4,p2p5+p3p5,p4p
     -   5,C0456R,C0256R,C0246R,C0245R,C456R,C256R,C246R,C245R,C0456I
     -   ,C0256I,C0246I,C0245I,C456I,C256I,C246I,C245I,D02456,D02456R
     -   ,D02456I,D2456R,D2456I)
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
       call NoAbe4vertHexLineFFhl1(F(1))
       call NoAbe4vertHexLineFFhl2(F(15))
       call NoAbe4vertHexLineFFhl3(F(29))
       call NoAbe4vertHexLineFFhl4(F(43))
       call NoAbe4vertHexLineFFhl5(F(57))
       call NoAbe4vertHexLineFFhl6(F(71))
       call NoAbe4vertHexLineFFhl7(F(85))
       call NoAbe4vertHexLineFFhl8(F(99))
       call NoAbe4vertHexLineFFhl9(F(113))
       call NoAbe4vertHexLineFFhl10(F(127))
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
       call NoAbe4vertHexLineFa1(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2
     -   mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3
     -   ,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6m
     -   up5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5
     -   mup6,Fa(1))
       call NoAbe4vertHexLineFa2(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2
     -   mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3
     -   ,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6m
     -   up5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5
     -   mup6,Fa(20))
       call NoAbe4vertHexLineFa3(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2
     -   mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3
     -   ,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6m
     -   up5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5
     -   mup6,Fa(39))
       call NoAbe4vertHexLineFa4(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2
     -   mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3
     -   ,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6m
     -   up5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5
     -   mup6,Fa(58))
       call NoAbe4vertHexLineFa5(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2
     -   mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3
     -   ,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6m
     -   up5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5
     -   mup6,Fa(77))
       call NoAbe4vertHexLineFa6(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2
     -   mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3
     -   ,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6m
     -   up5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5
     -   mup6,Fa(96))
       call NoAbe4vertHexLineFa7(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2
     -   mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3
     -   ,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6m
     -   up5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5
     -   mup6,Fa(115))
       call NoAbe4vertHexLineFa8(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2
     -   mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3
     -   ,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6m
     -   up5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5
     -   mup6,Fa(134))
       call NoAbe4vertHexLineFa9(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p2
     -   mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup3
     -   ,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6m
     -   up5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup5
     -   mup6,Fa(153))
       call NoAbe4vertHexLineFa10(p1mup2,p1mup3,p1mup5,p1mup6,p2mup2,p
     -   2mup3,p2mup5,p2mup6,p3mup2,p3mup3,p3mup5,p3mup6,p4mup2,p4mup
     -   3,p4mup5,p4mup6,p5mup2,p5mup3,p5mup5,p5mup6,p6mup2,p6mup3,p6
     -   mup5,p6mup6,mup2mup3,mup2mup5,mup2mup6,mup3mup5,mup3mup6,mup
     -   5mup6,Fa(172))
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
       SMB(4) = SC1c(barpsi_p4,mup6,psi_p1,alpha)
       SMB(5) = SC1c(barpsi_p4,mup3,psi_p1,alpha)
       SMB(6) = SC1r(barpsi_p4,p2,psi_p1,alpha)
       SMB(7) = SC1r(barpsi_p4,p3,psi_p1,alpha)
       SMB(8) = SC3ccc(barpsi_p4,mup6,mup5,mup3,psi_p1,alpha)
       SMB(9) = SC3ccr(barpsi_p4,mup6,mup5,p2,psi_p1,alpha)
       SMB(10) = SC3ccr(barpsi_p4,mup6,mup5,p3,psi_p1,alpha)
       SMB(11) = SC3ccr(barpsi_p4,mup3,mup2,p2,psi_p1,alpha)
       SMB(12) = SC3ccr(barpsi_p4,mup3,mup2,p3,psi_p1,alpha)
       SMB(13) = SC3ccr(barpsi_p4,mup5,mup2,p2,psi_p1,alpha)
       SMB(14) = SC3ccc(barpsi_p4,mup5,mup3,mup2,psi_p1,alpha)
       SMB(15) = SC3ccr(barpsi_p4,mup5,mup3,p2,psi_p1,alpha)
       SMB(16) = SC5ccccr(barpsi_p4,mup6,mup5,mup3,mup2,p2,psi_p1,alpha)
       SMB(17) = SC3ccr(barpsi_p4,mup5,mup2,p3,psi_p1,alpha)
       SMB(18) = SC3ccr(barpsi_p4,mup5,mup3,p3,psi_p1,alpha)
       SMB(19) = SC5ccccr(barpsi_p4,mup6,mup5,mup3,mup2,p3,psi_p1,alpha)
       SMB(20) = SC3ccc(barpsi_p4,mup6,mup3,mup2,psi_p1,alpha)
       SMB(21) = SC3ccr(barpsi_p4,mup6,mup2,p2,psi_p1,alpha)
       SMB(22) = SC3ccr(barpsi_p4,mup6,mup3,p2,psi_p1,alpha)
       SMB(23) = SC3ccr(barpsi_p4,mup6,mup2,p3,psi_p1,alpha)
       SMB(24) = SC3ccr(barpsi_p4,mup6,mup3,p3,psi_p1,alpha)
       SMB(25) = SC3crr(barpsi_p4,mup2,p3,p2,psi_p1,alpha)
       SMB(26) = SC5cccrr(barpsi_p4,mup6,mup5,mup2,p3,p2,psi_p1,alpha)
       SMB(27) = SC3crr(barpsi_p4,mup3,p3,p2,psi_p1,alpha)
       SMB(28) = SC3crr(barpsi_p4,mup5,p3,p2,psi_p1,alpha)
       SMB(29) = SC5cccrr(barpsi_p4,mup6,mup5,mup3,p3,p2,psi_p1,alpha)
       SMB(30) = SC3crr(barpsi_p4,mup6,p3,p2,psi_p1,alpha)
       SMB(31) = SC5cccrr(barpsi_p4,mup5,mup3,mup2,p3,p2,psi_p1,alpha)
       SMB(32) = SC5cccrr(barpsi_p4,mup6,mup3,mup2,p3,p2,psi_p1,alpha)
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
       result_graph(1,1) = Fa(1)*SMB(1)+Fa(9)*SMB(5)-8*(Fa(4)*SMB(2)-F
     -   a(12)*SMB(6)+Fa(15)*SMB(7))+Fa(18)*SMB(8)+F(84)*SMB(16)+F(86
     -   )*SMB(19)+Fa(31)*SMB(20)+Fa(32)*SMB(21)+Fa(36)*SMB(24)-4*(Fa
     -   (7)*SMB(4)-Fa(20)*SMB(9)-Fa(21)*SMB(10)-Fa(24)*SMB(13)-Fa(29
     -   )*SMB(18)+Fa(45)*SMB(30))-2*(Fa(6)*SMB(3)-p4mup6*F(1)*SMB(14
     -   )+Fa(34)*SMB(22)+Fa(35)*SMB(23)-Fa(40)*SMB(26)+Fa(44)*SMB(29
     -   )-p4mup6*F(132)*(-2*(p3mup5*SMB(11)-p2mup5*SMB(12)+p3mup2*SM
     -   B(15)+mup3mup5*SMB(25)-mup2mup5*SMB(27))+SMB(31))-p1mup5*F(8
     -   1)*SMB(32))
       result_graph(1,2) = Fa(10)*SMB(5)+Fa(13)*SMB(6)+Fa(16)*SMB(7)+2
     -   *(mup5mup6*F(70)*SMB(11)+mup5mup6*F(72)*SMB(12))+4*(Fa(2)*SM
     -   B(1)-Fa(38)*SMB(25)+Fa(41)*SMB(27))
       result_graph(1,3) = Fa(11)*SMB(5)-8*(Fa(8)*SMB(4)-Fa(14)*SMB(6)
     -   -Fa(17)*SMB(7))+Fa(23)*SMB(12)+Fa(25)*SMB(13)+Fa(26)*SMB(14)
     -   +F(85)*SMB(16)+Fa(30)*SMB(18)+F(87)*SMB(19)-4*(Fa(3)*SMB(1)+
     -   Fa(5)*SMB(2)+Fa(20)*SMB(9)+Fa(21)*SMB(10)-Fa(33)*SMB(21)+p3m
     -   up2*p4mup5*F(132)*SMB(22)-Fa(37)*SMB(24)-Fa(39)*SMB(25)+Fa(4
     -   2)*SMB(27)+Fa(43)*SMB(28))+2*(Fa(6)*SMB(3)+Fa(19)*SMB(8)-Fa(
     -   22)*SMB(11)-Fa(27)*SMB(15)-Fa(28)*SMB(17)+p4mup5*F(1)*SMB(20
     -   )-Fa(40)*SMB(26)+Fa(44)*SMB(29)+p1mup6*F(81)*SMB(31)+p4mup5*
     -   F(132)*SMB(32))
       result_graph(2,1) = (Fa(48)*SMB(2)+Fa(52)*SMB(5)+Fa(55)*SMB(6)+
     -   s34*Fa(75)*SMB(21)+Fa(77)*SMB(22)+F(92)*(-SMB(19)+2*(p4mup3*
     -   SMB(3)-p3mup6*SMB(14)-mup3mup6*SMB(17)+mup3mup5*SMB(23)))+Fa
     -   (79)*SMB(24)+4*(Fa(46)*SMB(1)+Fa(50)*SMB(4)-Fa(58)*SMB(7)-Fa
     -   (62)*SMB(9)-Fa(63)*SMB(11)+Fa(70)*SMB(15)+Fa(72)*SMB(18)-Fa(
     -   81)*SMB(25)-Fa(86)*SMB(28)+Fa(88)*SMB(30))+F(126)*(s34*SMB(1
     -   6)+4*(Fa(65)*SMB(12)+Fa(83)*SMB(27))+2*(Fa(67)*SMB(13)+Fa(89
     -   )*SMB(31)))-2*(s34*Fa(60)*SMB(8)+s34*Fa(74)*SMB(20)+Fa(87)*S
     -   MB(29)-p1mup5*F(103)*SMB(32)))/s34
       result_graph(2,2) = (-2*(mup5mup6*s34*F(126)*SMB(11)+mup5mup6*F
     -   (93)*(-2*p4mup3*SMB(1)+SMB(12))+2*(s34*Fa(53)*SMB(5)+2*Fa(56
     -   )*SMB(6)+Fa(84)*SMB(27))))/s34
       result_graph(2,3) = (Fa(51)*SMB(4)+Fa(54)*SMB(5)+Fa(57)*SMB(6)+
     -   Fa(64)*SMB(11)+Fa(66)*SMB(12)+s34*Fa(68)*SMB(13)+Fa(71)*SMB(
     -   15)+Fa(73)*SMB(18)+F(92)*(SMB(19)-2*(p4mup3*SMB(3)-mup3mup6*
     -   SMB(17)+p3mup5*SMB(20)+mup3mup5*SMB(23)))+4*(Fa(47)*SMB(1)+F
     -   a(49)*SMB(2)-Fa(59)*SMB(7)+Fa(62)*SMB(9)+Fa(78)*SMB(22)+Fa(8
     -   0)*SMB(24)-Fa(82)*SMB(25)-Fa(85)*SMB(27)+Fa(86)*SMB(28)-Fa(8
     -   8)*SMB(30))+2*(s34*Fa(61)*SMB(8)-s34*Fa(69)*SMB(14)+Fa(87)*S
     -   MB(29)+p1mup6*F(103)*SMB(31))+F(126)*(-(s34*SMB(16))+2*(Fa(7
     -   6)*SMB(21)+Fa(90)*SMB(32))))/s34
       result_graph(3,1) = (-8*Fa(93)*SMB(2)-Fa(96)*SMB(4)+Fa(106)*SMB
     -   (8)-4*(Fa(108)*SMB(9)-Fa(109)*SMB(10)-Fa(111)*SMB(13))+F(117
     -   )*SMB(16)+Fa(116)*SMB(21)+s12*(2*Fa(95)*SMB(3)+F(110)*(SMB(1
     -   9)-2*(p3mup5*SMB(20)+mup3mup5*SMB(23))))+2*(Fa(123)*SMB(26)+
     -   F(127)*(4*(Fa(98)*SMB(5)-Fa(103)*SMB(7))+2*Fa(114)*SMB(18)+p
     -   4mup6*(4*mup3mup5*p3mup2*SMB(6)-2*(p3mup5*SMB(11)+p3mup2*SMB
     -   (15)+mup3mup5*SMB(25))+SMB(31)))+F(110)*(-2*(Fa(118)*SMB(22)
     -   -Fa(119)*SMB(24))+Fa(125)*SMB(32))))/s12
       result_graph(3,2) = (Fa(99)*SMB(5)-8*(Fa(101)*SMB(6)-Fa(104)*SM
     -   B(7))+mup5mup6*F(116)*SMB(11)-2*mup5mup6*s12*F(110)*SMB(12)+
     -   4*(s12*Fa(91)*SMB(1)+Fa(121)*SMB(25)))/s12
       result_graph(3,3) = (-(Fa(94)*SMB(2))+Fa(100)*SMB(5)-8*(Fa(97)*
     -   SMB(4)+Fa(102)*SMB(6)-Fa(105)*SMB(7))+Fa(107)*SMB(8)+Fa(110)
     -   *SMB(11)+Fa(112)*SMB(13)+F(118)*SMB(16)+s12*(4*Fa(92)*SMB(1)
     -   -2*(Fa(95)*SMB(3)+F(110)*(-(mup5mup6*SMB(12))+p3mup6*SMB(14)
     -   +mup3mup6*SMB(17)))+F(119)*SMB(19))+4*(Fa(108)*SMB(9)-Fa(109
     -   )*SMB(10)+Fa(117)*SMB(21)+Fa(122)*SMB(25))-2*(Fa(123)*SMB(26
     -   )-F(110)*(-2*(Fa(113)*SMB(15)-Fa(115)*SMB(18))+Fa(124)*SMB(3
     -   1)))+F(127)*(-4*(p3mup2*p4mup5*SMB(22)-Fa(120)*SMB(24))+2*p4
     -   mup5*SMB(32)))/s12
       result_graph(4,1) = (-4*(Fa(126)*SMB(1)+Fa(130)*SMB(4)-Fa(132)*
     -   SMB(5))+8*Fa(134)*SMB(6)+Fa(136)*SMB(7)+F(123)*(4*Fa(128)*SM
     -   B(2)+s234*(-SMB(19)+2*(p4mup3*SMB(3)+p3mup2*SMB(8)-mup3mup6*
     -   SMB(17)+mup3mup5*SMB(23)))+2*(s34*Fa(143)*SMB(14)-2*(Fa(142)
     -   *SMB(13)+Fa(145)*SMB(15))+Fa(146)*SMB(18)+Fa(158)*SMB(31)))+
     -   2*(s234*s34*Fa(148)*SMB(20)-Fa(152)*SMB(24)-2*(Fa(138)*SMB(1
     -   1)-Fa(140)*SMB(12)+Fa(154)*SMB(25)-Fa(156)*SMB(27))+p1mup5*F
     -   (120)*(-2*(p4mup3*SMB(21)+p3mup2*SMB(22))+SMB(32))))/(s234*s
     -   34)
       result_graph(4,2) = (2*mup5mup6*F(123)*(2*(p4mup3*SMB(1)+p3mup2
     -   *SMB(5))-SMB(12)))/s34
       result_graph(4,3) = (8*Fa(135)*SMB(6)+Fa(137)*SMB(7)+Fa(141)*SM
     -   B(12)-4*(Fa(127)*SMB(1)+Fa(129)*SMB(2)-Fa(133)*SMB(5)+Fa(139
     -   )*SMB(11)+Fa(155)*SMB(25)-Fa(157)*SMB(27))+2*(s234*s34*Fa(14
     -   4)*SMB(14)-Fa(147)*SMB(18)+p1mup6*F(120)*(-2*(p4mup3*SMB(13)
     -   +p3mup2*SMB(15))+SMB(31)))+F(123)*(s234*(SMB(19)-2*(p4mup3*S
     -   MB(3)+p3mup2*SMB(8)-mup3mup6*SMB(17)+mup3mup5*SMB(23)))+2*(s
     -   34*Fa(149)*SMB(20)+2*(Fa(131)*SMB(4)-Fa(150)*SMB(21)-Fa(151)
     -   *SMB(22))+Fa(153)*SMB(24)+Fa(159)*SMB(32))))/(s234*s34)
       result_graph(5,1) = (Fa(162)*SMB(2)+4*Fa(164)*SMB(4)+8*(Fa(166)
     -   *SMB(5)+Fa(169)*SMB(6))-2*(4*Fa(171)*SMB(7)+Fa(173)*SMB(8)+F
     -   a(179)*SMB(21)+2*(Fa(174)*SMB(11)+Fa(182)*SMB(25)))+F(124)*(
     -   s34*SMB(16)-4*(Fa(177)*SMB(15)-Fa(178)*SMB(18))+2*(Fa(67)*SM
     -   B(13)+Fa(89)*SMB(31)))+F(121)*(s12*(4*Fa(160)*SMB(1)+SMB(19)
     -   -2*(p4mup3*SMB(3)-p3mup6*SMB(14)-mup3mup6*SMB(17)+p3mup5*SMB
     -   (20)+mup3mup5*SMB(23)))-4*(Fa(118)*SMB(22)-Fa(119)*SMB(24))+
     -   2*Fa(125)*SMB(32)))/(s12*s34)
       result_graph(5,2) = (-4*Fa(167)*SMB(5)+2*(mup5mup6*s34*F(124)*S
     -   MB(11)+mup5mup6*s12*F(121)*(-2*p4mup3*SMB(1)+SMB(12))))/(s12
     -   *s34)
       result_graph(5,3) = (Fa(165)*SMB(4)+Fa(168)*SMB(5)+8*(Fa(170)*S
     -   MB(6)-Fa(172)*SMB(7))+Fa(175)*SMB(11)+2*(Fa(173)*SMB(8)-Fa(1
     -   76)*SMB(13))+4*(Fa(163)*SMB(2)-Fa(183)*SMB(25))+F(121)*(-4*(
     -   Fa(113)*SMB(15)-Fa(115)*SMB(18))-s12*(4*Fa(161)*SMB(1)+SMB(1
     -   9)-2*(p4mup3*SMB(3)+mup5mup6*SMB(12)-p3mup6*SMB(14)-mup3mup6
     -   *SMB(17)+p3mup5*SMB(20)+mup3mup5*SMB(23)))+2*Fa(124)*SMB(31)
     -   )+F(124)*(-(s34*SMB(16))-4*(Fa(180)*SMB(22)-Fa(181)*SMB(24))
     -   +2*(Fa(76)*SMB(21)+Fa(90)*SMB(32))))/(s12*s34)
       result_graph(6,1) = (2*F(125)*(-4*Fa(184)*SMB(2)+2*(Fa(190)*SMB
     -   (13)+Fa(114)*SMB(18))+p4mup6*(s12*SMB(14)-2*p3mup2*SMB(15)+S
     -   MB(31)))+F(122)*(-8*Fa(186)*SMB(4)-s123*SMB(16)+4*(Fa(196)*S
     -   MB(21)-Fa(198)*SMB(22)+Fa(199)*SMB(24))+2*(s123*Fa(189)*SMB(
     -   8)+s12*Fa(195)*SMB(20)+Fa(195)*SMB(32))))/(s12*s123)
       result_graph(6,2) = (2*F(122)*(2*Fa(188)*SMB(5)-mup5mup6*SMB(11
     -   )))/s12
       result_graph(6,3) = (F(122)*(-8*Fa(185)*SMB(2)+s123*(4*Fa(188)*
     -   SMB(5)-2*(Fa(189)*SMB(8)+mup5mup6*SMB(11))+SMB(16))+2*(2*(Fa
     -   (191)*SMB(13)-Fa(193)*SMB(15)+Fa(194)*SMB(18))+Fa(192)*(s12*
     -   SMB(14)+SMB(31))))+2*F(125)*(-4*Fa(187)*SMB(4)+2*(Fa(197)*SM
     -   B(21)+Fa(120)*SMB(24))+p4mup5*(s12*SMB(20)-2*p3mup2*SMB(22)+
     -   SMB(32))))/(s12*s123)
       result_graph(1,1) =-result_graph(1,1)
       result_graph(1,2) =-result_graph(1,2)
       result_graph(1,3) =-result_graph(1,3)
       result_graph(2,1) =-result_graph(2,1)
       result_graph(2,2) =-result_graph(2,2)
       result_graph(2,3) =-result_graph(2,3)
       result_graph(3,1) =-result_graph(3,1)
       result_graph(3,2) =-result_graph(3,2)
       result_graph(3,3) =-result_graph(3,3)
       result_graph(4,1) =-result_graph(4,1)
       result_graph(4,2) =-result_graph(4,2)
       result_graph(4,3) =-result_graph(4,3)
       result_graph(5,1) =-result_graph(5,1)
       result_graph(5,2) =-result_graph(5,2)
       result_graph(5,3) =-result_graph(5,3)
       result_graph(6,1) =-result_graph(6,1)
       result_graph(6,2) =-result_graph(6,2)
       result_graph(6,3) =-result_graph(6,3)
       If (ngluon.eq.2) then
       result_color(1)=(result_graph(1,1)+result_graph(2,1)+result_gra
     -   ph(3,1)+result_graph(4,1)+result_graph(5,1)+result_graph(6,1
     -   ))
       result_color(2)=(result_graph(1,2)+result_graph(2,2)+result_gra
     -   ph(3,2)+result_graph(4,2)+result_graph(5,2)+result_graph(6,2
     -   ))
       result_color(3)=(result_graph(1,3)+result_graph(2,3)+result_gra
     -   ph(3,3)+result_graph(4,3)+result_graph(5,3)+result_graph(6,3
     -   ))
       elseif (ngluon.eq.3) then
       If (posgluon.eq.2) then
       result_color(1)=(result_graph(1,1))
       result_color(2)=(result_graph(1,1))
       elseIf (posgluon.eq.3) then
       result_color(1)=(result_graph(1,1))
       result_color(2)=(result_graph(1,1))
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
