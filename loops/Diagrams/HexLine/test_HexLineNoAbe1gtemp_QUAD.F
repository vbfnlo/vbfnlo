       subroutine test_hexlineNoAbe1g_QUAD(p1t,p2t,p3t,p4t,p5t,p6t,barpsi_p5,psi_
     -   p1,mup2,mup3,mup4,mup6,alpha,musqIn,ngluon,posgluon,comp,resul
     -   t,resultborn)
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 08/09/2009
c Modified:8/9/2009
c ************************************************************************************
c determine the rational terms for the non-abelian finite QCD virtual corrections along the quark line i.e 
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
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
       Real*8   p1t(0:3),p2t(0:3),p3t(0:3),p4t(0:3),p5t(0:3),p6t(0:3)
       Complex*16   mup2(0:3),mup3(0:3),mup4(0:3),mup6(0:3)
       Complex*16   p1mup2, p1mup3, p1mup4, p1mup6, p2mup2, p2mup3, 
     -          p2mup4, p2mup6, p3mup2, p3mup3, p3mup4, p3mup6, 
     -          p4mup2, p4mup3, p4mup4, p4mup6, p5mup2, p5mup3, 
     -          p5mup4, p5mup6, p6mup2, p6mup3, p6mup4, p6mup6
       Complex*16   mup2mup3, mup2mup4, mup2mup6, mup3mup4, mup3mup6, 
     -          mup4mup6

c       Complex*16 B0fin,C0fin,D0fin,E0fin,F0fin
c       EXTERNAL dotrr,B0fin,C0fin,D0fin,E0fin,F0fin
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


       Integer comp 
      Complex*16 SMB(64),F(6395),Fa(587)
      Real*8 P(7129)
      Complex*16 dotrc,dotcc,result(10),resultborn,resulttemp(3)
       Real*8 musq
      EXTERNAL   dotrc,dotcc
      Integer alpha
       COMMON/test_FaFunctionsNoAbe1g/Fa
       COMMON/test_FhlFunctionsNoABE1g/F
      Save/test_FhlFunctionsNoABE1g/
       COMMON/test_PFunctionsNoABE1g/P
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
       Complex*16 barpsi_p5(2),psi_p1(2)
       Integer ngluon, posgluon,jtemp
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
       do jtemp=0,3
         p1(jtemp)=p1t(jtemp)
         p2(jtemp)=p2t(jtemp)
         p3(jtemp)=p3t(jtemp)
         p4(jtemp)=p4t(jtemp)
         p5(jtemp)=p5t(jtemp)
         p6(jtemp)=p6t(jtemp)
      enddo

       if (comp.gt.0) then


         call PrintIntHexMassLess_QUAD(p1,p2,p3,p4,p5,p6,musq)

c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call test_FFhlNoABE1g1(F(1))
       call test_FFhlNoABE1g2(F(214))
       call test_FFhlNoABE1g3(F(427))
       call test_FFhlNoABE1g4(F(640))
       call test_FFhlNoABE1g5(F(853))
       call test_FFhlNoABE1g6(F(1066))
       call test_FFhlNoABE1g7(F(1279))
       call test_FFhlNoABE1g8(F(1492))
       call test_FFhlNoABE1g9(F(1705))
       call test_FFhlNoABE1g10(F(1918))
       call test_FFhlNoABE1g11(F(2131))
       call test_FFhlNoABE1g12(F(2344))
       call test_FFhlNoABE1g13(F(2557))
       call test_FFhlNoABE1g14(F(2770))
       call test_FFhlNoABE1g15(F(2983))
       call test_FFhlNoABE1g16(F(3196))
       call test_FFhlNoABE1g17(F(3409))
       call test_FFhlNoABE1g18(F(3622))
       call test_FFhlNoABE1g19(F(3835))
       call test_FFhlNoABE1g20(F(4048))
       call test_FFhlNoABE1g21(F(4261))
       call test_FFhlNoABE1g22(F(4474))
       call test_FFhlNoABE1g23(F(4687))
       call test_FFhlNoABE1g24(F(4900))
       call test_FFhlNoABE1g25(F(5113))
       call test_FFhlNoABE1g26(F(5326))
       call test_FFhlNoABE1g27(F(5539))
       call test_FFhlNoABE1g28(F(5752))
       call test_FFhlNoABE1g29(F(5965))
       call test_FFhlNoABE1g30(F(6178))
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
       p1mup6 = dotrc(p1,mup6)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup4 = dotrc(p2,mup4)
       p2mup6 = dotrc(p2,mup6)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup4 = dotrc(p3,mup4)
       p3mup6 = dotrc(p3,mup6)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup4 = dotrc(p4,mup4)
       p4mup6 = dotrc(p4,mup6)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrc(p5,mup3)
       p5mup4 = dotrc(p5,mup4)
       p5mup6 = dotrc(p5,mup6)
       p6mup2 = dotrc(p6,mup2)
       p6mup3 = dotrc(p6,mup3)
       p6mup4 = dotrc(p6,mup4)
       p6mup6 = dotrc(p6,mup6)
       mup2mup3 = dotcc(mup2,mup3)
       mup2mup4 = dotcc(mup2,mup4)
       mup2mup6 = dotcc(mup2,mup6)
       mup3mup4 = dotcc(mup3,mup4)
       mup3mup6 = dotcc(mup3,mup6)
       mup4mup6 = dotcc(mup4,mup6)
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call test_FaFunctionhexlineNoAbe1g1(p1mup2,p1mup3,p1mup4,p1mup6
     -   ,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4m
     -   up2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2,
     -   p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mup
     -   3mup6,mup4mup6,Fa(1))
       call test_FaFunctionhexlineNoAbe1g2(p1mup2,p1mup3,p1mup4,p1mup6
     -   ,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4m
     -   up2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2,
     -   p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mup
     -   3mup6,mup4mup6,Fa(15))
       call test_FaFunctionhexlineNoAbe1g3(p1mup2,p1mup3,p1mup4,p1mup6
     -   ,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4m
     -   up2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2,
     -   p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mup
     -   3mup6,mup4mup6,Fa(29))
       call test_FaFunctionhexlineNoAbe1g4(p1mup2,p1mup3,p1mup4,p1mup6
     -   ,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4m
     -   up2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2,
     -   p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mup
     -   3mup6,mup4mup6,Fa(43))
       call test_FaFunctionhexlineNoAbe1g5(p1mup2,p1mup3,p1mup4,p1mup6
     -   ,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4m
     -   up2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2,
     -   p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mup
     -   3mup6,mup4mup6,Fa(57))
       call test_FaFunctionhexlineNoAbe1g6(p1mup2,p1mup3,p1mup4,p1mup6
     -   ,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4m
     -   up2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2,
     -   p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mup
     -   3mup6,mup4mup6,Fa(71))
       call test_FaFunctionhexlineNoAbe1g7(p1mup2,p1mup3,p1mup4,p1mup6
     -   ,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4m
     -   up2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2,
     -   p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mup
     -   3mup6,mup4mup6,Fa(85))
       call test_FaFunctionhexlineNoAbe1g8(p1mup2,p1mup3,p1mup4,p1mup6
     -   ,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4m
     -   up2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2,
     -   p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mup
     -   3mup6,mup4mup6,Fa(99))
       call test_FaFunctionhexlineNoAbe1g9(p1mup2,p1mup3,p1mup4,p1mup6
     -   ,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4m
     -   up2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2,
     -   p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mup
     -   3mup6,mup4mup6,Fa(113))
       call test_FaFunctionhexlineNoAbe1g10(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(127))
       call test_FaFunctionhexlineNoAbe1g11(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(141))
       call test_FaFunctionhexlineNoAbe1g12(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(155))
       call test_FaFunctionhexlineNoAbe1g13(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(169))
       call test_FaFunctionhexlineNoAbe1g14(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(183))
       call test_FaFunctionhexlineNoAbe1g15(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(197))
       call test_FaFunctionhexlineNoAbe1g16(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(211))
       call test_FaFunctionhexlineNoAbe1g17(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(225))
       call test_FaFunctionhexlineNoAbe1g18(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(239))
       call test_FaFunctionhexlineNoAbe1g19(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(253))
       call test_FaFunctionhexlineNoAbe1g20(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(267))
       call test_FaFunctionhexlineNoAbe1g21(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(281))
       call test_FaFunctionhexlineNoAbe1g22(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(295))
       call test_FaFunctionhexlineNoAbe1g23(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(309))
       call test_FaFunctionhexlineNoAbe1g24(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(323))
       call test_FaFunctionhexlineNoAbe1g25(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(337))
       call test_FaFunctionhexlineNoAbe1g26(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(351))
       call test_FaFunctionhexlineNoAbe1g27(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(365))
       call test_FaFunctionhexlineNoAbe1g28(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(379))
       call test_FaFunctionhexlineNoAbe1g29(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(393))
       call test_FaFunctionhexlineNoAbe1g30(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(407))
       call test_FaFunctionhexlineNoAbe1g31(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(421))
       call test_FaFunctionhexlineNoAbe1g32(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(435))
       call test_FaFunctionhexlineNoAbe1g33(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(449))
       call test_FaFunctionhexlineNoAbe1g34(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(463))
       call test_FaFunctionhexlineNoAbe1g35(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(477))
       call test_FaFunctionhexlineNoAbe1g36(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(491))
       call test_FaFunctionhexlineNoAbe1g37(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(505))
       call test_FaFunctionhexlineNoAbe1g38(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(519))
       call test_FaFunctionhexlineNoAbe1g39(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(533))
       call test_FaFunctionhexlineNoAbe1g40(p1mup2,p1mup3,p1mup4,p1mup
     -   6,p2mup2,p2mup3,p2mup4,p2mup6,p3mup2,p3mup3,p3mup4,p3mup6,p4
     -   mup2,p4mup3,p4mup4,p4mup6,p5mup2,p5mup3,p5mup4,p5mup6,p6mup2
     -   ,p6mup3,p6mup4,p6mup6,mup2mup3,mup2mup4,mup2mup6,mup3mup4,mu
     -   p3mup6,mup4mup6,Fa(547))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the Matrix Element  
c************************************************************************************
c************************************************************************************
c************************************************************************************
       SMB(1) = SC3ccr(barpsi_p5,mup3,mup2,p2,psi_p1,alpha)
       SMB(2) = SC3ccr(barpsi_p5,mup4,mup2,p2,psi_p1,alpha)
       SMB(3) = SC3ccc(barpsi_p5,mup4,mup3,mup2,psi_p1,alpha)
       SMB(4) = SC3ccr(barpsi_p5,mup4,mup3,p2,psi_p1,alpha)
       SMB(5) = SC5ccccr(barpsi_p5,mup6,mup4,mup3,mup2,p2,psi_p1,alpha)
       SMB(6) = SC3ccr(barpsi_p5,mup3,mup2,p3,psi_p1,alpha)
       SMB(7) = SC3ccr(barpsi_p5,mup4,mup2,p3,psi_p1,alpha)
       SMB(8) = SC3ccr(barpsi_p5,mup4,mup3,p3,psi_p1,alpha)
       SMB(9) = SC5ccccr(barpsi_p5,mup6,mup4,mup3,mup2,p3,psi_p1,alpha)
       SMB(10) = SC3ccr(barpsi_p5,mup3,mup2,p4,psi_p1,alpha)
       SMB(11) = SC3ccr(barpsi_p5,mup4,mup2,p4,psi_p1,alpha)
       SMB(12) = SC3ccr(barpsi_p5,mup4,mup3,p4,psi_p1,alpha)
       SMB(13) = SC5ccccr(barpsi_p5,mup6,mup4,mup3,mup2,p4,psi_p1,alpha)
       SMB(14) = SC1c(barpsi_p5,mup2,psi_p1,alpha)
       SMB(15) = SC1c(barpsi_p5,mup3,psi_p1,alpha)
       SMB(16) = SC1c(barpsi_p5,mup4,psi_p1,alpha)
       SMB(17) = SC3ccc(barpsi_p5,mup6,mup3,mup2,psi_p1,alpha)
       SMB(18) = SC3ccc(barpsi_p5,mup6,mup4,mup2,psi_p1,alpha)
       SMB(19) = SC3ccc(barpsi_p5,mup6,mup4,mup3,psi_p1,alpha)
       SMB(20) = SC1r(barpsi_p5,p2,psi_p1,alpha)
       SMB(21) = SC3ccr(barpsi_p5,mup6,mup4,p2,psi_p1,alpha)
       SMB(22) = SC1r(barpsi_p5,p3,psi_p1,alpha)
       SMB(23) = SC3ccr(barpsi_p5,mup6,mup4,p3,psi_p1,alpha)
       SMB(24) = SC1r(barpsi_p5,p4,psi_p1,alpha)
       SMB(25) = SC3ccr(barpsi_p5,mup6,mup4,p4,psi_p1,alpha)
       SMB(26) = SC1c(barpsi_p5,mup6,psi_p1,alpha)
       SMB(27) = SC3crr(barpsi_p5,mup2,p3,p2,psi_p1,alpha)
       SMB(28) = SC5cccrr(barpsi_p5,mup6,mup4,mup2,p3,p2,psi_p1,alpha)
       SMB(29) = SC3crr(barpsi_p5,mup2,p4,p2,psi_p1,alpha)
       SMB(30) = SC5cccrr(barpsi_p5,mup6,mup4,mup2,p4,p2,psi_p1,alpha)
       SMB(31) = SC3ccr(barpsi_p5,mup6,mup2,p2,psi_p1,alpha)
       SMB(32) = SC3crr(barpsi_p5,mup3,p4,p2,psi_p1,alpha)
       SMB(33) = SC5cccrr(barpsi_p5,mup6,mup4,mup3,p4,p2,psi_p1,alpha)
       SMB(34) = SC3crr(barpsi_p5,mup3,p4,p3,psi_p1,alpha)
       SMB(35) = SC5cccrr(barpsi_p5,mup6,mup4,mup3,p4,p3,psi_p1,alpha)
       SMB(36) = SC3ccr(barpsi_p5,mup6,mup3,p2,psi_p1,alpha)
       SMB(37) = SC3ccr(barpsi_p5,mup6,mup3,p3,psi_p1,alpha)
       SMB(38) = SC3crr(barpsi_p5,mup3,p3,p2,psi_p1,alpha)
       SMB(39) = SC5ccrrr(barpsi_p5,mup3,mup2,p4,p3,p2,psi_p1,alpha)
       SMB(40) = SC5cccrr(barpsi_p5,mup4,mup3,mup2,p3,p2,psi_p1,alpha)
       SMB(41) = SC5cccrr(barpsi_p5,mup6,mup4,mup3,p3,p2,psi_p1,alpha)
       SMB(42) = SC7ccccrrr(barpsi_p5,mup6,mup4,mup3,mup2,p4,p3,p2,psi_p1,alpha)
       SMB(43) = SC5cccrr(barpsi_p5,mup6,mup3,mup2,p3,p2,psi_p1,alpha)
       SMB(44) = SC3ccr(barpsi_p5,mup6,mup3,p4,psi_p1,alpha)
       SMB(45) = SC5cccrr(barpsi_p5,mup6,mup3,mup2,p4,p2,psi_p1,alpha)
       SMB(46) = SC3ccr(barpsi_p5,mup6,mup2,p3,psi_p1,alpha)
       SMB(47) = SC3ccr(barpsi_p5,mup6,mup2,p4,psi_p1,alpha)
       SMB(48) = SC3crr(barpsi_p5,mup4,p4,p2,psi_p1,alpha)
       SMB(49) = SC3crr(barpsi_p5,mup4,p4,p3,psi_p1,alpha)
       SMB(50) = SC3crr(barpsi_p5,mup2,p4,p3,psi_p1,alpha)
       SMB(51) = SC5cccrr(barpsi_p5,mup6,mup4,mup2,p4,p3,psi_p1,alpha)
       SMB(52) = SC3crr(barpsi_p5,mup4,p3,p2,psi_p1,alpha)
       SMB(53) = SC5ccrrr(barpsi_p5,mup4,mup2,p4,p3,p2,psi_p1,alpha)
       SMB(54) = SC5cccrr(barpsi_p5,mup4,mup3,mup2,p4,p2,psi_p1,alpha)
       SMB(55) = SC5cccrr(barpsi_p5,mup6,mup3,mup2,p4,p3,psi_p1,alpha)
       SMB(56) = SC5cccrr(barpsi_p5,mup4,mup3,mup2,p4,p3,psi_p1,alpha)
       SMB(57) = SC5ccrrr(barpsi_p5,mup4,mup3,p4,p3,p2,psi_p1,alpha)
       SMB(58) = SC3crr(barpsi_p5,mup6,p3,p2,psi_p1,alpha)
       SMB(59) = SC3rrr(barpsi_p5,p4,p3,p2,psi_p1,alpha)
       SMB(60) = SC5ccrrr(barpsi_p5,mup6,mup4,p4,p3,p2,psi_p1,alpha)
       SMB(61) = SC3crr(barpsi_p5,mup6,p4,p2,psi_p1,alpha)
       SMB(62) = SC3crr(barpsi_p5,mup6,p4,p3,psi_p1,alpha)
       SMB(63) = SC5ccrrr(barpsi_p5,mup6,mup2,p4,p3,p2,psi_p1,alpha)
       SMB(64) = SC5ccrrr(barpsi_p5,mup6,mup3,p4,p3,p2,psi_p1,alpha)
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
       result(1) = (Fa(2)*SMB(2)+Fa(3)*SMB(3)+s45*F(1)*SMB(5)+Fa(8)*SM
     -   B(11)+Fa(9)*SMB(12)+s345*F(1)*SMB(13)-Fa(11)*SMB(15)+p5mup4*
     -   s345*F(2)*SMB(17)+p4mup3*s345*F(2)*SMB(18)+Fa(13)*SMB(19)+8*
     -   (F(1)*Fa(14)*SMB(20)+F(1)*Fa(16)*SMB(22)+F(1)*Fa(18)*SMB(24)
     -   -F(1)*Fa(20)*SMB(26))+p4mup3*F(2)*SMB(28)+F(2)*Fa(23)*SMB(30
     -   )+p3mup2*F(2)*SMB(33)+F(2)*Fa(26)*SMB(35)+p4mup6*F(2)*SMB(40
     -   )+F(1)*SMB(42)+p5mup4*F(2)*SMB(43)+4*(Fa(5)*SMB(6)+Fa(6)*SMB
     -   (7)-Fa(7)*SMB(8)-Fa(10)*SMB(14)-Fa(12)*SMB(16)-F(1)*Fa(15)*S
     -   MB(21)-F(1)*Fa(17)*SMB(23)-F(1)*Fa(19)*SMB(25)-F(1)*Fa(21)*S
     -   MB(27)-F(1)*Fa(22)*SMB(29)-F(1)*Fa(24)*SMB(31)-mup4mup6*p3mu
     -   p2*F(1)*SMB(32)-F(1)*Fa(25)*SMB(34)-p3mup2*p5mup4*F(1)*SMB(3
     -   6)-F(1)*Fa(27)*SMB(37)+F(1)*Fa(28)*SMB(38)+F(1)*Fa(29)*SMB(4
     -   8)+F(1)*Fa(30)*SMB(49)+F(1)*Fa(31)*SMB(52))+mup3mup6*F(2)*SM
     -   B(53)-2*(F(1)*Fa(1)*SMB(1)+F(1)*Fa(4)*SMB(4)+mup4mup6*s345*F
     -   (1)*SMB(10)+mup4mup6*F(1)*SMB(39)+p4mup2*F(1)*SMB(41)+p3mup6
     -   *F(1)*SMB(54)-Fa(32)*SMB(56)+mup2mup6*F(1)*SMB(57)))/(s16*s3
     -   45*s45)
       result(2) = (Fa(34)*SMB(2)+4*(Fa(36)*SMB(4)+Fa(39)*SMB(8)-Fa(41
     -   )*SMB(12)+Fa(43)*SMB(15)-Fa(44)*SMB(16))+F(10)*(s45*SMB(5)+8
     -   *(Fa(46)*SMB(20)-Fa(47)*SMB(22)-Fa(49)*SMB(24)+Fa(51)*SMB(26
     -   ))+SMB(42)-4*(Fa(15)*SMB(21)-Fa(48)*SMB(23)-Fa(50)*SMB(25)+F
     -   a(21)*SMB(27)+Fa(22)*SMB(29)+Fa(24)*SMB(31)-Fa(52)*SMB(34)-F
     -   a(54)*SMB(37)+Fa(55)*SMB(49))+2*(p4mup3*SMB(28)+Fa(23)*SMB(3
     -   0)-Fa(53)*SMB(35)-mup4mup6*(2*(p3mup2*SMB(32)-p4mup2*SMB(38)
     -   )+SMB(39))+p4mup6*SMB(40)+p5mup4*SMB(43)+p3mup2*(SMB(33)-2*(
     -   p5mup4*SMB(36)-mup3mup6*SMB(48)))-p4mup2*(SMB(41)+2*mup3mup6
     -   *SMB(52))+mup3mup6*SMB(53)))-2*(Fa(33)*SMB(1)+Fa(45)*SMB(19)
     -   -Fa(56)*SMB(54))+s12*(s45*Fa(35)*SMB(3)+Fa(40)*SMB(11)+4*(Fa
     -   (37)*SMB(6)+Fa(38)*SMB(7)+Fa(42)*SMB(14))+3*s345*F(16)*(-SMB
     -   (13)+2*(mup4mup6*SMB(10)-p5mup4*SMB(17)-p4mup3*SMB(18)))+2*F
     -   a(57)*SMB(56)))/(s12*s345*s45)
       result(3) = (-4*(Fa(59)*SMB(2)-Fa(61)*SMB(4)-Fa(62)*SMB(8)-Fa(6
     -   4)*SMB(15))+8*Fa(65)*SMB(16)+s12*(Fa(60)*SMB(3)+F(23)*(SMB(1
     -   3)-2*(mup4mup6*(SMB(10)+2*p4mup3*SMB(14))-p5mup4*SMB(17)-p4m
     -   up3*SMB(18))))+Fa(66)*SMB(19)+Fa(60)*SMB(40)+F(23)*(8*(Fa(67
     -   )*SMB(20)-Fa(69)*SMB(22)+Fa(70)*SMB(24)-Fa(72)*SMB(26))+4*(F
     -   a(68)*SMB(21)+Fa(48)*SMB(23)-Fa(71)*SMB(25)+Fa(73)*SMB(29)+F
     -   a(75)*SMB(31)+Fa(52)*SMB(34)-p3mup2*p5mup4*SMB(36)+Fa(54)*SM
     -   B(37))+SMB(42)+2*(p4mup3*SMB(28)-Fa(74)*SMB(30)+p3mup2*SMB(3
     -   3)-Fa(53)*SMB(35)-mup4mup6*(2*(p4mup3*SMB(27)+p3mup2*SMB(32)
     -   -p4mup2*SMB(38))+SMB(39))-p4mup2*SMB(41)+p5mup4*SMB(43)))+s1
     -   23*(Fa(58)*SMB(1)-3*s45*F(17)*SMB(5)+4*Fa(63)*SMB(12)+2*Fa(7
     -   6)*SMB(54)))/(s12*s123*s45)
       result(4) = (s123*F(26)*SMB(5)+s12*(F(26)*SMB(13)-2*F(27)*Fa(82
     -   )*SMB(17)+p4mup3*F(28)*SMB(18))+Fa(84)*SMB(21)+F(30)*Fa(85)*
     -   SMB(26)+F(28)*(p4mup3*SMB(28)+p3mup2*SMB(33))+F(29)*(Fa(48)*
     -   SMB(23)+Fa(87)*SMB(36))+F(26)*SMB(42)-2*(s56*(2*(Fa(77)*SMB(
     -   2)-Fa(79)*SMB(4)+Fa(80)*SMB(8))-4*Fa(81)*SMB(16)+Fa(78)*(s12
     -   *SMB(3)+SMB(40)))+F(27)*(Fa(83)*SMB(19)+Fa(74)*SMB(30)+Fa(53
     -   )*SMB(35)+2*(Fa(71)*SMB(25)+Fa(86)*SMB(31)+Fa(88)*SMB(37))+p
     -   4mup2*SMB(41)+Fa(82)*SMB(43))))/(s12*s123*s56)
       result(5) = (Fa(89)*SMB(1)+Fa(90)*SMB(2)+Fa(92)*SMB(4)+Fa(93)*S
     -   MB(6)+Fa(94)*SMB(7)+Fa(95)*SMB(8)+Fa(96)*SMB(11)-Fa(97)*SMB(
     -   12)+Fa(98)*SMB(14)+Fa(99)*SMB(15)+Fa(100)*SMB(16)+s345*(F(88
     -   )*SMB(13)+F(70)*(p5mup4*SMB(17)+p4mup3*SMB(18)))+s45*(Fa(91)
     -   *SMB(3)+F(61)*SMB(5)+Fa(101)*SMB(19))+Fa(102)*SMB(20)+8*(Fa(
     -   104)*SMB(22)+Fa(106)*SMB(24)+Fa(108)*SMB(26))+Fa(110)*SMB(33
     -   )+mup4mup6*(s345*F(69)*SMB(10)+2*F(31)*SMB(39))+F(61)*SMB(42
     -   )-4*(Fa(103)*SMB(21)-Fa(105)*SMB(23)-Fa(107)*SMB(25)+Fa(109)
     -   *SMB(32)-Fa(111)*SMB(34)+Fa(113)*SMB(36)-Fa(114)*SMB(37)-Fa(
     -   115)*SMB(38)-Fa(116)*SMB(48)-Fa(117)*SMB(49)-Fa(118)*SMB(52)
     -   )+2*(Fa(112)*SMB(35)+p4mup2*F(31)*SMB(41)+F(37)*(p4mup3*SMB(
     -   28)+Fa(23)*SMB(30)-2*(Fa(21)*SMB(27)+Fa(22)*SMB(29)+Fa(24)*S
     -   MB(31))+p4mup6*SMB(40)+p5mup4*SMB(43)+mup3mup6*SMB(53))+Fa(1
     -   19)*SMB(54))+Fa(120)*SMB(56)+Fa(121)*SMB(57))/(s345*s45)
       result(6) = (Fa(122)*SMB(1)+Fa(123)*SMB(2)+Fa(125)*SMB(4)+s45*F
     -   (147)*SMB(5)+Fa(128)*SMB(8)+Fa(130)*SMB(12)+Fa(132)*SMB(15)+
     -   Fa(133)*SMB(16)+Fa(136)*SMB(20)-Fa(137)*SMB(21)-Fa(138)*SMB(
     -   22)+Fa(139)*SMB(23)-8*(Fa(140)*SMB(24)-Fa(142)*SMB(26))-Fa(1
     -   43)*SMB(27)+Fa(144)*SMB(28)+Fa(145)*SMB(29)+Fa(146)*SMB(30)+
     -   Fa(147)*SMB(31)+p3mup2*F(210)*SMB(33)+4*(Fa(141)*SMB(25)-p3m
     -   up2*p5mup4*F(209)*SMB(36))+F(211)*(Fa(52)*SMB(34)+Fa(54)*SMB
     -   (37))+mup4mup6*(-4*p3mup2*F(209)*SMB(32)+p4mup2*F(211)*SMB(3
     -   8)-2*F(212)*SMB(39))+Fa(148)*SMB(40)+F(212)*(-2*(Fa(53)*SMB(
     -   35)+p4mup2*SMB(41))+SMB(42))+p5mup4*F(210)*SMB(43)+Fa(150)*S
     -   MB(48)-Fa(151)*SMB(49)-Fa(154)*SMB(52)+Fa(155)*SMB(53)+Fa(15
     -   6)*SMB(54)+s12*(Fa(126)*SMB(6)+Fa(127)*SMB(7)+mup4mup6*F(157
     -   )*SMB(10)+Fa(129)*SMB(11)+F(166)*SMB(13)-Fa(131)*SMB(14)+p5m
     -   up4*F(139)*SMB(17)+s45*(Fa(124)*SMB(3)+3*F(148)*SMB(9)+Fa(13
     -   4)*SMB(18)+Fa(135)*SMB(19))-4*(Fa(149)*SMB(46)+Fa(152)*SMB(5
     -   0))+Fa(153)*SMB(51)+2*Fa(157)*SMB(56)))/(s12*s45)
       result(7) = (Fa(159)*SMB(2)+Fa(161)*SMB(4)+Fa(162)*SMB(8)+Fa(16
     -   6)*SMB(15)+s12*(Fa(160)*SMB(3)-Fa(163)*SMB(10)+F(254)*SMB(13
     -   )-Fa(165)*SMB(14)+Fa(168)*SMB(17)+p4mup3*F(265)*SMB(18))-Fa(
     -   169)*SMB(19)+Fa(557)*SMB(20)+Fa(179)*SMB(21)-Fa(180)*SMB(22)
     -   +F(270)*Fa(48)*SMB(23)+Fa(558)*SMB(24)+8*(Fa(167)*SMB(16)+Fa
     -   (187)*SMB(26))-Fa(188)*SMB(27)+Fa(189)*SMB(29)-Fa(191)*SMB(3
     -   2)+F(265)*(p4mup3*SMB(28)+p3mup2*SMB(33))+Fa(192)*SMB(34)+Fa
     -   (193)*SMB(36)-4*(Fa(190)*SMB(31)+Fa(194)*SMB(37))+Fa(195)*SM
     -   B(38)-Fa(196)*SMB(39)+Fa(197)*SMB(40)+F(254)*(-2*(2*Fa(71)*S
     -   MB(25)+Fa(74)*SMB(30)+Fa(53)*SMB(35)+p4mup2*SMB(41))+SMB(42)
     -   )+Fa(198)*SMB(43)+s123*(F(243)*SMB(5)-4*(Fa(164)*SMB(12)-Fa(
     -   199)*SMB(44))+Fa(200)*SMB(45)-2*(Fa(158)*SMB(1)+Fa(201)*SMB(
     -   54))))/(s12*s123)
       result(8) = (-(Fa(202)*SMB(1))-Fa(203)*SMB(2)+Fa(204)*SMB(3)+Fa
     -   (205)*SMB(4)-Fa(206)*SMB(6)+Fa(207)*SMB(7)+Fa(208)*SMB(8)+mu
     -   p4mup6*F(375)*SMB(10)+Fa(209)*SMB(11)+Fa(210)*SMB(12)+F(410)
     -   *SMB(13)-Fa(211)*SMB(14)+Fa(212)*SMB(15)+Fa(559)*SMB(16)+p5m
     -   up4*F(584)*SMB(17)+Fa(560)*SMB(20)+Fa(229)*SMB(21)+Fa(561)*S
     -   MB(22)+Fa(239)*SMB(23)+Fa(240)*SMB(24)+Fa(241)*SMB(25)+Fa(24
     -   2)*SMB(26)+Fa(245)*SMB(29)+Fa(246)*SMB(30)+Fa(247)*SMB(31)+F
     -   a(248)*SMB(32)+Fa(249)*SMB(33)+Fa(250)*SMB(34)+Fa(251)*SMB(3
     -   5)+Fa(252)*SMB(36)+Fa(253)*SMB(37)-Fa(254)*SMB(38)+mup4mup6*
     -   F(937)*SMB(39)+s45*(F(333)*SMB(5)+F(374)*SMB(9)+Fa(218)*SMB(
     -   18)+Fa(219)*SMB(19)+Fa(244)*SMB(28)+Fa(255)*SMB(40)+Fa(256)*
     -   SMB(41))+F(944)*SMB(42)+p5mup4*F(945)*SMB(43)+Fa(562)*SMB(48
     -   )+Fa(563)*SMB(49)+4*(Fa(243)*SMB(27)+Fa(257)*SMB(46)+Fa(268)
     -   *SMB(50))-Fa(269)*SMB(51)+Fa(270)*SMB(52)-Fa(271)*SMB(53)-Fa
     -   (272)*SMB(54)-Fa(273)*SMB(56)+Fa(274)*SMB(57)+16*(Fa(275)*SM
     -   B(58)-Fa(276)*SMB(59))+8*Fa(277)*SMB(60))/s45
       result(9) = (-(Fa(278)*SMB(1))+Fa(279)*SMB(2)-Fa(281)*SMB(4)+F(
     -   1053)*SMB(5)+Fa(284)*SMB(8)+Fa(287)*SMB(12)+Fa(295)*SMB(15)+
     -   Fa(296)*SMB(16)+Fa(299)*SMB(19)+Fa(565)*SMB(20)+Fa(314)*SMB(
     -   21)+Fa(566)*SMB(22)+Fa(326)*SMB(23)+Fa(567)*SMB(24)+Fa(338)*
     -   SMB(25)+Fa(568)*SMB(26)+Fa(569)*SMB(27)+Fa(351)*SMB(28)+Fa(5
     -   70)*SMB(29)+Fa(359)*SMB(30)-Fa(360)*SMB(31)+Fa(361)*SMB(32)+
     -   p3mup2*F(1515)*SMB(33)+4*Fa(362)*SMB(34)+Fa(363)*SMB(36)-Fa(
     -   364)*SMB(37)+Fa(365)*SMB(38)+Fa(366)*SMB(39)+Fa(367)*SMB(40)
     -   +F(1516)*(Fa(53)*SMB(35)+p4mup2*SMB(41))+F(1541)*SMB(42)-Fa(
     -   368)*SMB(43)-Fa(369)*SMB(44)-Fa(370)*SMB(45)-Fa(373)*SMB(48)
     -   +Fa(374)*SMB(49)+Fa(377)*SMB(52)-Fa(378)*SMB(53)+Fa(379)*SMB
     -   (54)+s12*(Fa(280)*SMB(3)+Fa(282)*SMB(6)+Fa(283)*SMB(7)+F(108
     -   5)*SMB(9)+F(1109)*SMB(13)+Fa(564)*SMB(14)+Fa(297)*SMB(17)+Fa
     -   (298)*SMB(18)+4*(Fa(285)*SMB(10)-Fa(286)*SMB(11)+Fa(371)*SMB
     -   (46)+Fa(372)*SMB(47)-Fa(375)*SMB(50))+Fa(376)*SMB(51)+Fa(380
     -   )*SMB(55)+Fa(381)*SMB(56))-16*(Fa(382)*SMB(58)-Fa(383)*SMB(6
     -   1)+Fa(384)*SMB(62))+8*Fa(385)*SMB(63))/s12
       result(10) = Fa(386)*SMB(1)+Fa(387)*SMB(2)+Fa(388)*SMB(3)+Fa(38
     -   9)*SMB(4)+F(1675)*SMB(5)+Fa(390)*SMB(6)+Fa(391)*SMB(7)+Fa(39
     -   2)*SMB(8)+F(1754)*SMB(9)+Fa(393)*SMB(10)+Fa(394)*SMB(11)+Fa(
     -   395)*SMB(12)+F(1833)*SMB(13)+Fa(571)*SMB(14)+Fa(572)*SMB(15)
     -   +Fa(573)*SMB(16)+Fa(414)*SMB(17)+Fa(415)*SMB(18)+Fa(416)*SMB
     -   (19)+Fa(574)*SMB(20)+Fa(432)*SMB(21)+Fa(575)*SMB(22)+Fa(448)
     -   *SMB(23)+Fa(576)*SMB(24)+Fa(464)*SMB(25)+Fa(577)*SMB(26)+Fa(
     -   578)*SMB(27)+Fa(477)*SMB(28)+Fa(579)*SMB(29)+Fa(485)*SMB(30)
     -   +Fa(486)*SMB(31)+Fa(487)*SMB(32)+Fa(488)*SMB(33)+Fa(495)*SMB
     -   (35)+Fa(496)*SMB(36)+Fa(497)*SMB(37)+Fa(581)*SMB(38)+Fa(504)
     -   *SMB(39)+Fa(505)*SMB(40)+Fa(506)*SMB(41)+F(5375)*SMB(42)+Fa(
     -   507)*SMB(43)+4*(Fa(580)*SMB(34)+Fa(508)*SMB(44)-Fa(509)*SMB(
     -   45))+Fa(510)*SMB(46)+Fa(511)*SMB(47)+Fa(582)*SMB(48)+Fa(583)
     -   *SMB(49)+Fa(522)*SMB(50)+Fa(523)*SMB(51)+Fa(524)*SMB(52)+Fa(
     -   525)*SMB(53)+Fa(526)*SMB(54)+Fa(527)*SMB(55)+Fa(528)*SMB(56)
     -   +Fa(529)*SMB(57)+8*(Fa(585)*SMB(59)+Fa(544)*SMB(60)-2*(Fa(58
     -   4)*SMB(58)+Fa(586)*SMB(61)+Fa(587)*SMB(62))+Fa(555)*SMB(63)-
     -   Fa(556)*SMB(64))
c The Result must be multiplied by 1/2 to 
c recover the proper CA factor
       result(1) = (1/2.d0)*result(1)
       result(2) = (1/2.d0)*result(2)
       result(3) = (1/2.d0)*result(3)
       result(4) = (1/2.d0)*result(4)
       result(5) = (1/2.d0)*result(5)
       result(6) = (1/2.d0)*result(6)
       result(7) = (1/2.d0)*result(7)
       result(8) = (1/2.d0)*result(8)
       result(9) = (1/2.d0)*result(9)
       result(10) = (1/2.d0)*result(10)
c$$$       print*, ""
c$$$       print*, "test dbe", result
c$$$       print*, ""
       If (ngluon.eq.1) then
       result(1)=-(result(1)+result(2)+result(3)+result(4)+result(5)+r
     -   esult(6)+result(7)+result(8)+result(9)+result(10))
       elseif (ngluon.eq.2) then
c  One gluon case is:  result(i)*CA. True color factor is f^abc T^c T^b= (-I)/2 CA.
c  To restore I should multiply for 2*I
c The color factor should be provided from outside
c result 1: posglu inside loop
c result 2: posglu left of loop
c result 3: posglu right of loop
       If (posgluon.eq.2) then
	resulttemp(1)=-(0d0,2d0)*(result(10)+result(8)+result(5))
	resulttemp(2)=-(0d0,2d0)*(result(9)+result(6)+result(7)
     -                           +result(2)+result(3)+result(4))
	resulttemp(3)=-(0d0,2d0)*(result(1))
       elseif (posgluon.eq.3) then
        resulttemp(1)=-(0d0,2d0)*(result(10)+result(8)+result(9)+result(6))
        resulttemp(2)=-(0d0,2d0)*(result(7)+result(3)+result(4))
        resulttemp(3)=-(0d0,2d0)*(result(5)+result(1)+result(2))
       elseif (posgluon.eq.4) then
        resulttemp(1)=-(0d0,2d0)*(result(10)+result(9)+result(7))
        resulttemp(2)=-(0d0,2d0)*(result(4))
        resulttemp(3)=-(0d0,2d0)*(result(8)+result(5)+result(6)
     -                           +result(1)+result(2)+result(3))
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     -   . Look to the heading for explanation" 
       endif
       result(1) = resulttemp(1)
       result(2) = resulttemp(2)
       result(3) = resulttemp(3)
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     -   . Look to the heading for explanation" 
       endif
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
