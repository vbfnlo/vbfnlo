       subroutine NoAbe4vertHexLine_QUAD(p1t,p2t,p3t,p4t,p5t,p6t,barpsi_p4,psi_p1
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
      Real*8   p1t(0:3),p2t(0:3),p3t(0:3),p4t(0:3),p5t(0:3),p6t(0:3)
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
        Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   p6sq 
       Real*8   s12, s16 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45 
       Real*8   s56 
       Real*8   s123,s234,s345


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

       Real*8 E012345R,E012346R,E012356R,E012456R,E013456R,E023456R 
       Real*8 E012345I,E012346I,E012356I,E012456I,E013456I,E023456I  
       Real*8 E12345R(46,4),E12346R(46,4),E12356R(46,4) 
       Real*8 E12456R(46,4),E13456R(46,4),E23456R(46,4)  
       Real*8 E12345I(46,4),E12346I(46,4),E12356I(46,4) 
       Real*8 E12456I(46,4),E13456I(46,4),E23456I(46,4)  

       Real*8 F0123456R 
       Real*8 F0123456I 
       Real*8 F123456R(166,5) 
       Real*8 F123456I(166,5) 

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
       Integer ngluon, posgluon,jtemp
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

      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False. 
      PrintF=.False.
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      do jtemp=0,3
         p1(jtemp)=p1t(jtemp)
         p2(jtemp)=p2t(jtemp)
         p3(jtemp)=p3t(jtemp)
         p4(jtemp)=p4t(jtemp)
         p5(jtemp)=p5t(jtemp)
         p6(jtemp)=p6t(jtemp)
      enddo

        if (comp.gt.0) then    
c    Calling C_ij,D_ij,E_ij,F_ij Functions    
c************************************************************************************
c************************************************************************************

       call PrintIntHexMassLess4g_QUAD(p1,p2,p3,p4,p5,p6,musq)

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
