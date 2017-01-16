       subroutine NoAbe2gaHexLine_QUAD(p1t,p2t,p3t,p4t,p5t,p6t,barpsi_p4,psi_p1,m
     -   up2,mup3,mup5,mup6,alpha,musq,ngluon,posgluon,comp,result_gr
     -   aph,result_color)
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
       Integer ngluon, posgluon
       integer jtemp
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

c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
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


       call PrintIntHexMassLess2g_QUAD(p1,p2,p3,p4,p5,p6,musq)

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
