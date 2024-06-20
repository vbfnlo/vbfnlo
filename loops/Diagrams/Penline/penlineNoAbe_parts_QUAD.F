       subroutine penlineNoAbe_pmmm_QUAD(p1t,p2t,p3t,p4t,p5t,barpsi_p4,psi_p1,mu
     -   p2,mup3,mup5,alpha,musqIn,ngluon,posgluon,gaugetest,comp,resul
     -   tgauge,result)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 3/7/2008
c Modified:31/5/2012
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
c                        mu_p5                       
c                          p5                      
c                          $                     
c                          $                     
c                          $                     
c                          V                      
c                   $$$$$$$$$$$$$$$$$              
c                   $               $              
c                   $               $               
c                   $               $               
c psi(p1)    ---->------->------->------->---   bar_psi(p4)
c                      $       $                     
c                      $       $                     
c                      $       $                    
c                      V       V                    
c                      $       $                    
c                      $       $                    
c                     p2       p3                  
c                   mu_p2    mu_p3              
c Note: The vertices are just Gamma^(mu_p2)..terms. So the correct
c factor should be added by hand for each boson. i.e, Without adding
c anything it represents pp->gluon+2 photons up to the electromagnetic charge 
c (for uu->g+2 photons is missing (2/3)^2) from the coupling. For W and Z,
c we have to account for the Diracgamma_5 that are missing in the vertex
c Note: To make it shorter in the promgram: mu_p2,...->mup2,... 
c Notation of External momenta: p1+p2+p3+p4+p5=0 
c mu_p2,mu_p3,mu_p5, should be think as external current 
c alpha is the helicity of the initial spinor 
c musq is the renormalization scale energy  
c gaugetest,integer value.Different gauge test:
c gaugetest=1 give the result after replacing mu_p2 with p2 
c gaugetest=2 give the result after replacing mu_p3 with p3 
c gaugetest=3 give the result after replacing mu_p5 with p5 
c gaugetest=4 calculates the three different gaugue tests 
c gaugetest>5 no calculation of gauge test in case you
c have already done it before(it safes times)
c comp: integer value.The first time called with p1...p5, comp=1
c ATTENTION: ONLY!!!If you have to call the subroutine consecutively with the same arguments
c(p1,p2,p3,p4,p5). Then, comp=-1 (it safes 4000 lines of code) 
c This applies when you have for examples the same diagram for an off-shell photon
c and a Z boson. The differences are in the coupling and  the part that depends on the
c polarization vector that are calculated at the end of this program.
c resultgauge is an array of dimension three. 
c In case you use gaugetest=(1,3). The result is given in the first argument:
c resultgauge(1).The argument 2,and 3 is set to zero 
c In case you use gaugetest=4. resultgauge(1) is the result of gaugetest=1,
c resultgauge(2) of gaugetest=2 and resultgauge(3) of gaugetest=3
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
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
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
      Real*8   p1t(0:3),p2t(0:3),p3t(0:3),p4t(0:3),p5t(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3),mup5(0:3) 
       Real*8   p1sq, p2sq, 
     -          p3sq, p4sq, p5sq
       Real*8   s12, s15, s23, s34, s45
       Complex*16   p1mup2, p1mup3, p1mup5, p2mup2, p2mup3, p2mup5, 
     -          p3mup2, p3mup3, p3mup5, p4mup2, p4mup3, p4mup5, 
     -          p5mup2, p5mup3, p5mup5
       Complex*16 mup2mup3,mup2mup5,mup3mup5
       Real*8 dotrr
       Complex*16 B0fin,C0fin,D0fin,E0fin
       EXTERNAL dotrr,B0fin,C0fin,D0fin,E0fin

      Real*8 B012R,B013R,B014R,B015R,B023R,B024R,B025R,B034R,B035R,B045R
      Real*8 B012I,B013I,B014I,B015I,B023I,B024I,B025I,B034I,B035I,B045I
       Real*8 C0123R,C0124R,C0125R,C0134R,C0135R,C0145R,C0234R,C0235R,
     -   C0245R,C0345R
       Real*8 C0123I,C0124I,C0125I,C0134I,C0135I,C0145I,C0234I,C0235I,
     -   C0245I,C0345I
       Real*8 Cij123R(4,2),Cij124R(4,2),Cij125R(4,2),Cij134R(4,2),Cij1
     -   35R(4,2),Cij145R(4,2),Cij234R(4,2),Cij235R(4,2),Cij245R(4,2)
     -   ,Cij345R(4,2)
       Real*8 Cij123I(4,2),Cij124I(4,2),Cij125I(4,2),Cij134I(4,2),Cij1
     -   35I(4,2),Cij145I(4,2),Cij234I(4,2),Cij235I(4,2),Cij245I(4,2)
     -   ,Cij345I(4,2)

        Real*8 D01234R,D01235R,D01245R,D01345R,D02345R,EE0R
        Real*8 D01234I,D01235I,D01245I,D01345I,D02345I,EE0I
        Real*8 Dij1234R(13,3),Dij1235R(13,3),Dij1245R(13,3),Dij1345R(13
     -   ,3),Dij2345R(13,3),EijR(46,4)
        Real*8 Dij1234I(13,3),Dij1235I(13,3),Dij1245I(13,3),Dij1345I(13
     -   ,3),Dij2345I(13,3),EijI(46,4)
       Complex*16 SMB(16),SMG(16) ,Fa(50),F(438),K(7)
       Real*8 FI(438),FR(438),KI(7),  KR(7)
       Complex*16 SC1c,SC3ccc,SC3rcc,SC3rrc,SC5rrccc,SC1r,dotrc,dotcc,result(6)
       Complex*16 resultgauge(3) 
       Real*8 musq,Is12,Is34,Is15s34,Is45,P(759) 
       EXTERNAL   dotrc,dotcc,SC1c,SC1r,SC3ccc,SC3rcc,SC3rrc,SC5rrccc
       Integer alpha,comp,gaugetest,ngluon,posgluon
       common/PartsFfunctionsNOABEmmm/F
       common/PartsPfunctionsNOABEmmm/P
       common/PartsKfunctionsNOABEmmm/K
       SAVE/PartsFfunctionsNOABEmmm/
       SAVE/PartsPfunctionsNOABEmmm/
       SAVE/PartsKfunctionsNOABEmmm/
       common/PartsKinematicsNOABEmmm/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
       real*8 musqIn
       Common/musqInv/musq

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      End QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       COMMON/EVALFunctionsP/ EE0R,EijR,EE0I,EijI
       COMMON/DVALFunctionsP/ D01234R,DIJ1234R,D01234I,DIJ1234I, 
     &    D01235R,DIJ1235R,D01235I,DIJ1235I,
     &    D01245R,DIJ1245R,D01245I,DIJ1245I, 
     &    D01345R,DIJ1345R,D01345I,DIJ1345I, 
     &    D02345R,Dij2345R,D02345I,Dij2345I 

       COMMON/CVALFunctionsP/ C0123R,CIJ123R,C0123I,CIJ123I, 
     &    C0124R,CIJ124R,C0124I,CIJ124I, C0125R,CIJ125R,C0125I,CIJ125I,
     &    C0134R,CIJ134R,C0134I,CIJ134I, C0135R,CIJ135R,C0135I,CIJ135I,
     &    C0145R,CIJ145R,C0145I,CIJ145I, C0234R,CIJ234R,C0234I,CIJ234I, 
     &    C0235R,CIJ235R,C0235I,CIJ235I, C0245R,CIJ245R,C0245I,CIJ245I, 
     &    C0345R,CIJ345R,C0345I,CIJ345I
       COMMON/BVALFunctionsP/ B012R,B012I, B013R,B013I, B014R,B014I, 
     &    B015R,B015I,  B023R,B023I,  B024R,B024I, B025R,B025I,
     &    B034R,B034I,  B035R,B035I,   B045R,B045I
       real*8 p1sqT,p2sqT,p3sqT,p4sqT,p5sqT,s12T,s23T,s34T,s45T,
     &     s15T,Is12T,Is45T,Is12s45T
       COMMON/InvariantsP/p1sqT,p2sqT,p3sqT,p4sqT,p5sqT,s12T,s23T,s34T,s45T,
     &     s15T,Is12T,Is45T,Is12s45T
       integer jtemp
       real*8 Is12s45

       musq=musqIn


       do jtemp=0,3
         p1(jtemp)=p1t(jtemp)
         p2(jtemp)=p2t(jtemp)
         p3(jtemp)=p3t(jtemp)
         p4(jtemp)=p4t(jtemp)
         p5(jtemp)=p5t(jtemp)
      enddo
       If (comp.eq.1) then    
c    Calling C_ij,D_ij Functions    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
          
          call PrintIntPenMassless_QUAD(p1,p2,p3,p4,p5,musq)
          p1sq=p1sqT
          p2sq=p2sqT
          p3sq=p3sqT
          p4sq=p4sqT
          p5sq=p5sqT
          s12=s12T
          s23=s23T
          s34=s34T
          s45=s45T
          s15=s15T
          Is12=Is12T
          Is45=Is45T
          Is12s45=Is12s45T
          is34=1d0/s34
          Is15s34=1d0/(s15*s34)
          Is45=1d0/s45
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
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       P(1) = 2*p5sq+s12-s34
       P(2) = 2*p5sq+5*s12-s34
       P(3) = s12+s34
       P(4) = 2*p5sq-s34
       P(5) = 2*p5sq+4*s12-s34
       P(6) = p5sq+2*s12-s34
       P(7) = 3*s12-s34
       P(8) = p5sq+s12-s34
       P(9) = -3*s12+s34
       P(10) = p5sq-s34
       P(11) = -s12+2*P(10)
       P(12) = s12+2*P(10)
       P(13) = s12-2*P(10)
       P(14) = 2*p5sq+s45
       P(15) = 2*p5sq+5*s45
       P(16) = s12-s34
       P(17) = 2*p5sq-s45
       P(18) = p5sq+2*s45
       P(19) = -2*p5sq+s45
       P(20) = -2*p5sq+s15
       P(21) = p5sq-s15
       P(22) = 2*p5sq-s15
       P(23) = p3sq+s12+2*P(10)
       P(24) = s12+2*s34
       P(25) = -s45+2*P(8)
       P(26) = p3sq+6*s12-2*s45+4*P(10)
       P(27) = p3sq+3*s12+2*P(10)
       P(28) = p3sq+3*s12-2*s34
       P(29) = p3sq-2*s45+4*P(8)
       P(30) = -1+3*Is12*s34
       P(31) = -(p5sq*s34)+P(16)**2
       P(32) = -p5sq+s34
       P(33) = 5*s12**2-7*s12*s34+2*s34*P(32)
       P(34) = p3sq+2*s12
       P(35) = -2*s12*s45+3*s34*s45+P(16)*P(34)
       P(36) = -s12+s34
       P(37) = p5sq*s12-s12*s45
       P(38) = p5sq*s34-s34*s45
       P(39) = s12*s34-s34**2+p3sq*P(36)-2*P(37)+3*P(38)
       P(40) = p3sq+s12-s45
       P(41) = p3sq+s12+2*s34-s45
       P(42) = -p3sq-p5sq+s34+s45
       P(43) = p3sq+p5sq-s34-s45
       P(44) = 1-2*Is12*s34
       P(45) = p3sq+s12-2*s34
       P(46) = -s45+2*P(10)
       P(47) = s12-s45+2*P(10)
       P(48) = p5sq+s12-s45
       P(49) = p3sq-3*s34+2*P(48)
       P(50) = p5sq+s12
       P(51) = 2*p3sq-5*s34-3*s45+4*P(50)
       P(52) = p5sq-s45
       P(53) = p3sq+3*P(16)+2*P(52)
       P(54) = 2*p3sq+4*p5sq+7*s12-5*s34-3*s45
       P(55) = p3sq+s12-3*s34+2*P(52)
       P(56) = s12-s45
       P(57) = 2*p3sq+4*p5sq-5*s34+3*P(56)
       P(58) = p3sq+2*p5sq
       P(59) = 3-5*Is12*s34
       P(60) = s12+Is12*s34*P(32)
       P(61) = -4*s34+s45+2*P(60)
       P(62) = p3sq*s12+6*s12**2-9*s12*s34+2*s34*P(32)
       P(63) = -2*s12+s34
       P(64) = -4*s34*s45+p3sq*P(63)
       P(65) = 4*s34-3*P(56)+Is12*P(64)
       P(66) = p3sq-s34+3*P(52)
       P(67) = 2*s12+s34
       P(68) = s34*P(56)+p3sq*P(67)
       P(69) = p3sq-s34
       P(70) = -p3sq+s34
       P(71) = 3-7*Is12*s34
       P(72) = -3+5*Is12*s34
       P(73) = s12-Is12*s34*P(10)
       P(74) = 4*s34-s45-2*P(73)
       P(75) = 2*s12-s34
       P(76) = 3*s12-4*s34
       P(77) = p3sq*P(75)+P(56)*P(76)
       P(78) = -p3sq+s34-3*P(52)
       P(79) = 3*s12+2*P(10)
       P(80) = p3sq+s45+2*P(16)
       P(81) = p3sq+2*p5sq+5*s12-4*s34+s45
       P(82) = p3sq-p5sq-s34+s45
       P(83) = -p3sq+s12
       P(84) = p3sq-s12
       P(85) = p5sq+2*P(69)
       P(86) = s45-2*P(8)
       P(87) = 4*p3sq+p5sq-3*s45-2*P(3)
       P(88) = -2*p3sq+s12+s34-s45
       P(89) = p5sq+s12-s34-s45
       P(90) = p3sq+2*P(89)
       P(91) = -s12+s45
       P(92) = s45+2*P(16)
       P(93) = 3*s12+2*P(69)
       P(94) = -3*p5sq+s34+2*P(84)
       P(95) = s34-s45+2*P(84)
       P(96) = 2*s34-3*s45
       P(97) = p5sq-s45+2*P(69)
       P(98) = s34+s45
       P(99) = p5sq-P(98)
       P(100) = p3sq+s12+2*P(99)
       P(101) = p3sq+p5sq+s12-s34
       P(102) = -3*s45+2*P(101)
       P(103) = s12+4*P(10)
       P(104) = p3sq+s12-2*s34-s45
       P(105) = p3sq-s12-2*s34+s45
       P(106) = s45+2*P(8)
       P(107) = p5sq+3*s12-s34+s45
       P(108) = -p5sq+s34+s45
       P(109) = 2*s12+s45
       P(110) = p3sq-5*s12-2*s45
       P(111) = p3sq-s12+s45
       P(112) = p3sq-2*s45
       P(113) = 3*s45-2*P(84)
       P(114) = p3sq+2*P(56)
       P(115) = p3sq-s34+2*P(52)
       P(116) = -p3sq+s12+s45
       P(117) = 5*s12+4*P(10)
       P(118) = p5sq+s12-s34-2*s45
       P(119) = p3sq+2*P(118)
       P(120) = p5sq+s45
       P(121) = -3*s12-4*s34+2*P(120)
       P(122) = p3sq-2*s34+3*P(56)
       P(123) = p3sq+p5sq+s12-s34-2*s45
       P(124) = p3sq+s12
       P(125) = -3*p5sq+s34+2*P(124)
       P(126) = 3*s12-5*s34
       P(127) = p5sq*s34-P(16)**2
       P(128) = s12*s45-2*P(127)
       P(129) = p3sq*P(63)-P(56)*P(76)
       P(130) = -3+7*Is12*s34
       P(131) = s12+s34-s45
       P(132) = -p3sq+p5sq+s34-s45
       P(133) = -s34+s45-2*P(84)
       P(134) = 5*s12+2*P(10)
       P(135) = 7*s12+2*P(10)
       P(136) = p3sq+4*P(56)
       P(137) = -3*s45+2*P(84)
       P(138) = -p3sq+s34-2*P(52)
       P(139) = p3sq+3*P(56)
       P(140) = 3*p5sq+s12-5*s45+2*P(69)
       P(141) = 4*p3sq+p5sq-2*s34-3*s45
       P(142) = 3*s12+4*s34-2*P(120)
       P(143) = p2sq-s12
       P(144) = -p2sq+p5sq+s12+s15
       P(145) = p2sq+s12
       P(146) = p2sq+p5sq-s12-s15-2*s34
       P(147) = -p5sq+s15+2*P(3)
       P(148) = p2sq-2*s12+s34
       P(149) = p2sq-2*s12
       P(150) = p5sq-s12+2*s34
       P(151) = -p2sq+s12
       P(152) = p2sq-s15
       P(153) = -p5sq+s12+s34-2*P(152)
       P(154) = p2sq-s12+2*P(21)
       P(155) = p5sq-s12
       P(156) = s34+2*P(155)
       P(157) = p2sq+p5sq-s12+s15
       P(158) = p2sq+p5sq-s12-s15
       P(159) = p5sq-s15+2*P(143)
       P(160) = 2*p5sq-s12+s34
       P(161) = p2sq-s12+2*s15+s34
       P(162) = p2sq-s12-2*P(21)
       P(163) = s12-s15
       P(164) = p2sq+p5sq+3*s34-2*P(163)
       P(165) = p2sq+s12-2*s15-s34
       P(166) = p2sq+2*p5sq+s12-4*s15-s34
       P(167) = p2sq-p5sq+s12-s15-s34
       P(168) = p2sq-s12+s15
       P(169) = p2sq-s15+s34
       P(170) = 3*p2sq-s15+s34
       P(171) = 2*p5sq-3*s12+s34
       P(172) = p2sq-p5sq-s12-s15
       P(173) = p5sq-s15-2*P(3)
       P(174) = p5sq-s12-s34+2*P(152)
       P(175) = p2sq-3*s12
       P(176) = p2sq-s15+3*P(155)
       P(177) = -p5sq+s15
       P(178) = -p2sq-p5sq+s12+s15+2*s34
       P(179) = p2sq-4*s12
       P(180) = -p2sq+s15+s34-3*P(155)
       P(181) = -p5sq+s12-2*s34
       P(182) = -p2sq+s12+2*P(21)
       P(183) = -p2sq+p5sq-s12+s15+s34
       P(184) = -p2sq+s15-3*P(155)
       P(185) = -p2sq-p5sq+s12+s15
       P(186) = p5sq-s15+s23-s45
       P(187) = p2sq-s12-2*s15+s23-s45
       P(188) = p2sq-s12-s23-2*s34+s45
       P(189) = p2sq-s12-s15
       P(190) = -s23+s45-2*P(189)
       P(191) = p5sq+s15-s23-s45
       P(192) = p2sq+3*s12+s23+2*s34-s45
       P(193) = s23-s45
       P(194) = p2sq-3*s12-4*s34-2*P(193)
       FR(1) = 2*(Dij2345R(1,1)-Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(2
     -   ,2)+2*Dij2345R(4,2)+4*(EijR(37,4)+EijR(38,4))-8*(EijR(22,3)-
     -   EijR(23,3)+EijR(43,4))+(EijR(2,2)+2*(EijR(8,3)+EijR(11,3)-2*
     -   EijR(17,3)))*P(143)-EijR(2,3)*P(145)+(-EijR(10,3)-EijR(13,3)
     -   +2*(EijR(2,1)-EijR(3,1)+EijR(5,2)-EijR(6,2)+EijR(20,3)))*P(1
     -   86)+EijR(3,2)*P(187)-EijR(3,3)*P(188)+EijR(8,2)*P(190)+(EijR
     -   (9,2)-EijR(10,2))*P(191)+EijR(9,3)*P(192)+EijR(12,3)*P(194))
       FI(1) = 2*(Dij2345I(1,1)-Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(2
     -   ,2)+2*Dij2345I(4,2)+4*(EijI(37,4)+EijI(38,4))-8*(EijI(22,3)-
     -   EijI(23,3)+EijI(43,4))+(EijI(2,2)+2*(EijI(8,3)+EijI(11,3)-2*
     -   EijI(17,3)))*P(143)-EijI(2,3)*P(145)+(-EijI(10,3)-EijI(13,3)
     -   +2*(EijI(2,1)-EijI(3,1)+EijI(5,2)-EijI(6,2)+EijI(20,3)))*P(1
     -   86)+EijI(3,2)*P(187)-EijI(3,3)*P(188)+EijI(8,2)*P(190)+(EijI
     -   (9,2)-EijI(10,2))*P(191)+EijI(9,3)*P(192)+EijI(12,3)*P(194))
       F(1)=DCMPLX(FR(1),FI(1))
       P(195) = s15-s23+s45
       P(196) = p2sq-s15-s34
       P(197) = 4*s12+s23-s45-2*P(196)
       P(198) = -p5sq+s15-3*s23+s45+2*P(16)
       P(199) = s23-s45+2*P(3)
       P(200) = s15-s23
       P(201) = p5sq-s45-2*P(16)-3*P(200)
       P(202) = p2sq-p5sq+s12+s15-s23+s45
       FR(2) = 2*(-Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(2,1)+Dij2345R(
     -   4,2)+Dij2345R(5,2)-Dij2345R(6,2)-4*(EijR(22,3)-EijR(24,3)-Ei
     -   jR(37,4)+EijR(43,4)+EijR(44,4)-EijR(45,4))-EijR(2,3)*P(145)+
     -   EijR(20,3)*P(173)+EijR(2,2)*P(175)+(2*(EijR(2,1)+EijR(5,2))+
     -   EijR(15,3)-EijR(16,3))*P(186)+(EijR(3,2)+EijR(12,3)-EijR(13,
     -   3))*P(188)-2*(p5sq*EijR(4,1)+EijR(7,2)*P(21)-(EijR(8,3)-EijR
     -   (17,3)-EijR(18,3)+EijR(19,3))*P(143)+EijR(6,2)*P(193)-EijR(3
     -   ,1)*P(195))+EijR(8,2)*P(197)+EijR(9,2)*P(198)+EijR(9,3)*P(19
     -   9)+EijR(10,2)*P(201)+EijR(10,3)*P(202))
       FI(2) = 2*(-Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(2,1)+Dij2345I(
     -   4,2)+Dij2345I(5,2)-Dij2345I(6,2)-4*(EijI(22,3)-EijI(24,3)-Ei
     -   jI(37,4)+EijI(43,4)+EijI(44,4)-EijI(45,4))-EijI(2,3)*P(145)+
     -   EijI(20,3)*P(173)+EijI(2,2)*P(175)+(2*(EijI(2,1)+EijI(5,2))+
     -   EijI(15,3)-EijI(16,3))*P(186)+(EijI(3,2)+EijI(12,3)-EijI(13,
     -   3))*P(188)-2*(p5sq*EijI(4,1)+EijI(7,2)*P(21)-(EijI(8,3)-EijI
     -   (17,3)-EijI(18,3)+EijI(19,3))*P(143)+EijI(6,2)*P(193)-EijI(3
     -   ,1)*P(195))+EijI(8,2)*P(197)+EijI(9,2)*P(198)+EijI(9,3)*P(19
     -   9)+EijI(10,2)*P(201)+EijI(10,3)*P(202))
       F(2)=DCMPLX(FR(2),FI(2))
       P(203) = p2sq+p3sq-s23-s34
       P(204) = p3sq+s15-s23-s34
       FR(3) = -4*(D02345R+Dij2345R(1,1)+s12*(EijR(2,2)-EijR(5,2))+s34
     -   *(EijR(6,2)-EijR(8,2))-p2sq*(EijR(2,3)+EijR(5,3)-2*EijR(8,3)
     -   )+4*EijR(21,3)-2*(EijR(36,4)+EijR(37,4)+2*(EijR(22,3)-EijR(4
     -   0,4)))+(EijR(6,3)+EijR(9,3)-2*EijR(17,3))*P(203)-(EijR(7,3)+
     -   EijR(10,3)-2*EijR(18,3))*P(204))
       FI(3) = -4*(D02345I+Dij2345I(1,1)+s12*(EijI(2,2)-EijI(5,2))+s34
     -   *(EijI(6,2)-EijI(8,2))-p2sq*(EijI(2,3)+EijI(5,3)-2*EijI(8,3)
     -   )+4*EijI(21,3)-2*(EijI(36,4)+EijI(37,4)+2*(EijI(22,3)-EijI(4
     -   0,4)))+(EijI(6,3)+EijI(9,3)-2*EijI(17,3))*P(203)-(EijI(7,3)+
     -   EijI(10,3)-2*EijI(18,3))*P(204))
       F(3)=DCMPLX(FR(3),FI(3))
       P(205) = p5sq-s12-s15+s34
       P(206) = p5sq-s15+s34
       P(207) = p2sq-3*s12+2*P(206)
       P(208) = s15+s34
       P(209) = p2sq-3*s12-s23+s45+2*P(208)
       P(210) = p2sq+3*s12+2*s15-s23-s45
       P(211) = p2sq+p3sq+s12-s34
       P(212) = -s23-s45+2*P(211)
       P(213) = p5sq+s15-s23-s45+2*P(69)
       P(214) = p2sq+s12-s23-s45+2*P(69)
       FR(4) = 2*(D02345R+Dij2345R(1,2)+Dij2345R(2,1)-Dij2345R(4,2)-6*
     -   EijR(11,2)+4*(EijR(37,4)-EijR(40,4)+EijR(41,4)-EijR(43,4))-(
     -   EijR(2,2)-EijR(2,3)+EijR(8,3))*P(145)-EijR(7,2)*P(186)-EijR(
     -   9,2)*P(191)+2*(EijR(22,3)-EijR(23,3)+EijR(3,1)*P(16)-(EijR(1
     -   ,1)+EijR(1,2))*P(21)+EijR(10,2)*P(52)-EijR(3,2)*P(56)+EijR(2
     -   ,1)*P(205))+EijR(5,2)*P(207)-EijR(6,2)*P(209)+EijR(8,2)*P(21
     -   0)+(-EijR(9,3)+EijR(17,3))*P(212)+(EijR(10,3)-EijR(18,3)+Eij
     -   R(19,3)-EijR(20,3))*P(213)+(-EijR(11,3)+EijR(12,3))*P(214))
       FI(4) = 2*(D02345I+Dij2345I(1,2)+Dij2345I(2,1)-Dij2345I(4,2)-6*
     -   EijI(11,2)+4*(EijI(37,4)-EijI(40,4)+EijI(41,4)-EijI(43,4))-(
     -   EijI(2,2)-EijI(2,3)+EijI(8,3))*P(145)-EijI(7,2)*P(186)-EijI(
     -   9,2)*P(191)+2*(EijI(22,3)-EijI(23,3)+EijI(3,1)*P(16)-(EijI(1
     -   ,1)+EijI(1,2))*P(21)+EijI(10,2)*P(52)-EijI(3,2)*P(56)+EijI(2
     -   ,1)*P(205))+EijI(5,2)*P(207)-EijI(6,2)*P(209)+EijI(8,2)*P(21
     -   0)+(-EijI(9,3)+EijI(17,3))*P(212)+(EijI(10,3)-EijI(18,3)+Eij
     -   I(19,3)-EijI(20,3))*P(213)+(-EijI(11,3)+EijI(12,3))*P(214))
       F(4)=DCMPLX(FR(4),FI(4))
       P(215) = p2sq-s12-s23+s45
       P(216) = p5sq+s15+s23-s45-2*P(16)
       P(217) = p2sq+s12-s23-s45
       P(218) = -p5sq+2*s12+s15+s23+s45
       P(219) = p2sq-p5sq+s12-s15+s23+s45-2*P(69)
       FR(5) = 2*(D02345R+Dij2345R(1,2)+Dij2345R(3,1)-Dij2345R(5,2)-6*
     -   EijR(11,2)+4*(EijR(37,4)-EijR(40,4)+EijR(42,4)-EijR(44,4))-(
     -   EijR(2,2)-EijR(2,3)+EijR(8,3))*P(145)+2*(EijR(22,3)-EijR(24,
     -   3)+EijR(4,1)*P(16)-(EijR(1,1)+EijR(1,2))*P(21)+EijR(4,2)*P(5
     -   2)-EijR(10,2)*P(56)+EijR(2,1)*P(205))+EijR(5,2)*P(207)+(EijR
     -   (14,3)-EijR(15,3))*P(213)+(-EijR(9,3)+EijR(17,3)-EijR(19,3)+
     -   EijR(20,3))*P(214)-EijR(6,2)*P(215)-EijR(7,2)*P(216)+EijR(8,
     -   2)*P(217)+EijR(9,2)*P(218)+(-EijR(10,3)+EijR(18,3))*P(219))
       FI(5) = 2*(D02345I+Dij2345I(1,2)+Dij2345I(3,1)-Dij2345I(5,2)-6*
     -   EijI(11,2)+4*(EijI(37,4)-EijI(40,4)+EijI(42,4)-EijI(44,4))-(
     -   EijI(2,2)-EijI(2,3)+EijI(8,3))*P(145)+2*(EijI(22,3)-EijI(24,
     -   3)+EijI(4,1)*P(16)-(EijI(1,1)+EijI(1,2))*P(21)+EijI(4,2)*P(5
     -   2)-EijI(10,2)*P(56)+EijI(2,1)*P(205))+EijI(5,2)*P(207)+(EijI
     -   (14,3)-EijI(15,3))*P(213)+(-EijI(9,3)+EijI(17,3)-EijI(19,3)+
     -   EijI(20,3))*P(214)-EijI(6,2)*P(215)-EijI(7,2)*P(216)+EijI(8,
     -   2)*P(217)+EijI(9,2)*P(218)+(-EijI(10,3)+EijI(18,3))*P(219))
       F(5)=DCMPLX(FR(5),FI(5))
       P(220) = p5sq-P(195)
       P(221) = -p2sq+s12+2*P(220)
       P(222) = p2sq-s12-s23+s45-2*P(21)
       P(223) = p5sq+s15-s23+s45
       P(224) = p2sq-s12-2*P(193)
       P(225) = p5sq+s15-s45
       P(226) = p2sq-s12-2*P(225)
       FR(6) = 2*(D02345R-Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(3,1)+Di
     -   j2345R(4,2)+Dij2345R(5,2)-Dij2345R(6,2)-8*EijR(22,3)-4*(EijR
     -   (11,2)-EijR(23,3)-EijR(24,3)-EijR(37,4)+EijR(43,4)+EijR(44,4
     -   )-EijR(45,4))-EijR(2,3)*P(145)+EijR(10,2)*P(172)+EijR(20,3)*
     -   P(173)+(2*EijR(5,2)+EijR(15,3)-EijR(16,3))*P(186)+(EijR(12,3
     -   )-EijR(13,3))*P(188)-EijR(4,2)*P(191)+2*(Dij2345R(1,1)+(EE0R
     -   +EijR(1,1)-EijR(6,2))*P(21)+(EijR(8,3)-EijR(17,3)-EijR(18,3)
     -   +EijR(19,3))*P(143)-EijR(7,2)*P(193))+EijR(9,3)*P(199)+EijR(
     -   10,3)*P(202)-EijR(3,2)*P(215)+EijR(2,1)*P(221)+EijR(3,1)*P(2
     -   22)+EijR(4,1)*P(223)+EijR(8,2)*P(224)-EijR(9,2)*P(226))
       FI(6) = 2*(D02345I-Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(3,1)+Di
     -   j2345I(4,2)+Dij2345I(5,2)-Dij2345I(6,2)-8*EijI(22,3)-4*(EijI
     -   (11,2)-EijI(23,3)-EijI(24,3)-EijI(37,4)+EijI(43,4)+EijI(44,4
     -   )-EijI(45,4))-EijI(2,3)*P(145)+EijI(10,2)*P(172)+EijI(20,3)*
     -   P(173)+(2*EijI(5,2)+EijI(15,3)-EijI(16,3))*P(186)+(EijI(12,3
     -   )-EijI(13,3))*P(188)-EijI(4,2)*P(191)+2*(Dij2345I(1,1)+(EE0I
     -   +EijI(1,1)-EijI(6,2))*P(21)+(EijI(8,3)-EijI(17,3)-EijI(18,3)
     -   +EijI(19,3))*P(143)-EijI(7,2)*P(193))+EijI(9,3)*P(199)+EijI(
     -   10,3)*P(202)-EijI(3,2)*P(215)+EijI(2,1)*P(221)+EijI(3,1)*P(2
     -   22)+EijI(4,1)*P(223)+EijI(8,2)*P(224)-EijI(9,2)*P(226))
       F(6)=DCMPLX(FR(6),FI(6))
       P(227) = p5sq-3*s15+s23-s45
       P(228) = s12+s15-s23-s34
       P(229) = 2*s12+s15-s23-s34
       P(230) = -p5sq+s15-s23+s45+2*P(145)
       P(231) = p2sq+s12-2*P(220)
       FR(7) = 2*(D02345R-Dij2345R(1,2)-Dij2345R(3,2)-4*(EijR(11,2)+Ei
     -   jR(22,3)-EijR(24,3)-EijR(37,4)-EijR(39,4))-8*EijR(44,4)-EijR
     -   (2,3)*P(145)-(EijR(4,3)-2*(EijR(5,2)-EijR(7,2)))*P(186)-(Eij
     -   R(9,3)+EijR(16,3)-2*EijR(20,3))*P(188)+EijR(3,1)*P(215)+EijR
     -   (2,1)*P(221)-EijR(4,1)*P(227)+2*(Dij2345R(5,2)-s12*EijR(2,2)
     -   +s34*EijR(8,2)-s34*EijR(10,2)+(EE0R+EijR(1,1))*P(21)+(EijR(8
     -   ,3)+EijR(14,3)-2*EijR(18,3))*P(143)-EijR(4,2)*P(228)+EijR(9,
     -   2)*P(229))+EijR(10,3)*P(230)-EijR(15,3)*P(231))
       FI(7) = 2*(D02345I-Dij2345I(1,2)-Dij2345I(3,2)-4*(EijI(11,2)+Ei
     -   jI(22,3)-EijI(24,3)-EijI(37,4)-EijI(39,4))-8*EijI(44,4)-EijI
     -   (2,3)*P(145)-(EijI(4,3)-2*(EijI(5,2)-EijI(7,2)))*P(186)-(Eij
     -   I(9,3)+EijI(16,3)-2*EijI(20,3))*P(188)+EijI(3,1)*P(215)+EijI
     -   (2,1)*P(221)-EijI(4,1)*P(227)+2*(Dij2345I(5,2)-s12*EijI(2,2)
     -   +s34*EijI(8,2)-s34*EijI(10,2)+(EE0I+EijI(1,1))*P(21)+(EijI(8
     -   ,3)+EijI(14,3)-2*EijI(18,3))*P(143)-EijI(4,2)*P(228)+EijI(9,
     -   2)*P(229))+EijI(10,3)*P(230)-EijI(15,3)*P(231))
       F(7)=DCMPLX(FR(7),FI(7))
       P(232) = s23-2*s45
       P(233) = p3sq+s12+s23-2*P(98)
       P(234) = p5sq-2*s45
       P(235) = p2sq+P(98)
       P(236) = p3sq+5*s12+s23-2*P(235)
       P(237) = s23+s34
       P(238) = p2sq-P(237)
       P(239) = -p3sq-s12+s45+2*P(238)
       P(240) = p5sq-2*s15+3*s23-4*P(16)
       P(241) = -p2sq+s23+2*s34
       P(242) = s12-s23
       P(243) = s15-s34
       P(244) = p3sq-3*P(242)-2*P(243)
       P(245) = p2sq-s12-s23-s34+s45
       P(246) = -p5sq+s12+s15-s23+s45
       P(247) = p2sq-p5sq+s15-s23+s45
       FR(8) = 2*(EE0R*s45-4*Dij2345R(1,1)+3*Dij2345R(3,1)-EijR(4,2)*P
     -   (52)+(-EijR(3,1)+EijR(10,2))*P(104)-EijR(1,2)*P(193)-EijR(1,
     -   1)*P(232)+EijR(2,1)*P(233)+EijR(4,1)*P(234)+EijR(5,2)*P(236)
     -   +EijR(6,2)*P(239)+EijR(7,2)*P(240)-EijR(9,2)*P(244)-2*(D0234
     -   5R+Dij2345R(1,2)-Dij2345R(5,2)+s12*EijR(2,3)-p2sq*EijR(8,3)-
     -   EijR(11,2)-2*(EijR(22,3)-EijR(24,3)+EijR(37,4)-EijR(40,4)+Ei
     -   jR(42,4)-EijR(44,4))-(-EijR(5,3)+EijR(7,3))*P(143)-EijR(2,2)
     -   *P(149)-(-EijR(14,3)+EijR(15,3))*P(186)-EijR(8,2)*P(241)-(-E
     -   ijR(9,3)+EijR(17,3)-EijR(19,3)+EijR(20,3))*P(245)-EijR(10,3)
     -   *P(246)+EijR(18,3)*P(247)))
       FI(8) = 2*(EE0I*s45-4*Dij2345I(1,1)+3*Dij2345I(3,1)-EijI(4,2)*P
     -   (52)+(-EijI(3,1)+EijI(10,2))*P(104)-EijI(1,2)*P(193)-EijI(1,
     -   1)*P(232)+EijI(2,1)*P(233)+EijI(4,1)*P(234)+EijI(5,2)*P(236)
     -   +EijI(6,2)*P(239)+EijI(7,2)*P(240)-EijI(9,2)*P(244)-2*(D0234
     -   5I+Dij2345I(1,2)-Dij2345I(5,2)+s12*EijI(2,3)-p2sq*EijI(8,3)-
     -   EijI(11,2)-2*(EijI(22,3)-EijI(24,3)+EijI(37,4)-EijI(40,4)+Ei
     -   jI(42,4)-EijI(44,4))-(-EijI(5,3)+EijI(7,3))*P(143)-EijI(2,2)
     -   *P(149)-(-EijI(14,3)+EijI(15,3))*P(186)-EijI(8,2)*P(241)-(-E
     -   ijI(9,3)+EijI(17,3)-EijI(19,3)+EijI(20,3))*P(245)-EijI(10,3)
     -   *P(246)+EijI(18,3)*P(247)))
       F(8)=DCMPLX(FR(8),FI(8))
       P(248) = p2sq+s34+s45
       P(249) = -p3sq-3*s12-s23+2*P(248)
       P(250) = p3sq+3*s12+s23-2*s34
       P(251) = 2*p2sq-P(250)
       P(252) = p5sq-s45+2*P(200)
       P(253) = p2sq-s34-s45
       P(254) = p3sq-s12-s23+2*P(253)
       P(255) = p2sq-2*s12-s23-s34+s45
       P(256) = s12+s23+s34-s45
       FR(9) = 2*(EE0R*s45+3*Dij2345R(2,1)-EijR(3,1)*P(45)+(EijR(4,1)-
     -   EijR(10,2))*P(52)+EijR(3,2)*P(104)-EijR(1,2)*P(193)-EijR(1,1
     -   )*P(232)+EijR(2,1)*P(233)-EijR(5,2)*P(249)+EijR(6,2)*P(251)-
     -   EijR(7,2)*P(252)-EijR(8,2)*P(254)-2*(Dij2345R(1,1)+Dij2345R(
     -   1,2)-Dij2345R(4,2)+s12*EijR(2,3)-p2sq*EijR(8,3)-2*(EijR(11,2
     -   )+EijR(21,3)-EijR(23,3)+EijR(37,4)-EijR(40,4)+EijR(41,4)-Eij
     -   R(43,4))-(EijR(2,2)-EijR(5,3)+EijR(6,3))*P(143)-(-EijR(10,3)
     -   +EijR(18,3)-EijR(19,3)+EijR(20,3))*P(186)-EijR(9,2)*P(191)-(
     -   -EijR(11,3)+EijR(12,3))*P(245)+EijR(9,3)*P(255)+EijR(17,3)*P
     -   (256)))
       FI(9) = 2*(EE0I*s45+3*Dij2345I(2,1)-EijI(3,1)*P(45)+(EijI(4,1)-
     -   EijI(10,2))*P(52)+EijI(3,2)*P(104)-EijI(1,2)*P(193)-EijI(1,1
     -   )*P(232)+EijI(2,1)*P(233)-EijI(5,2)*P(249)+EijI(6,2)*P(251)-
     -   EijI(7,2)*P(252)-EijI(8,2)*P(254)-2*(Dij2345I(1,1)+Dij2345I(
     -   1,2)-Dij2345I(4,2)+s12*EijI(2,3)-p2sq*EijI(8,3)-2*(EijI(11,2
     -   )+EijI(21,3)-EijI(23,3)+EijI(37,4)-EijI(40,4)+EijI(41,4)-Eij
     -   I(43,4))-(EijI(2,2)-EijI(5,3)+EijI(6,3))*P(143)-(-EijI(10,3)
     -   +EijI(18,3)-EijI(19,3)+EijI(20,3))*P(186)-EijI(9,2)*P(191)-(
     -   -EijI(11,3)+EijI(12,3))*P(245)+EijI(9,3)*P(255)+EijI(17,3)*P
     -   (256)))
       F(9)=DCMPLX(FR(9),FI(9))
       P(257) = 3*p2sq-s12
       P(258) = -p2sq+p5sq+s12-s15+s23-s45
       P(259) = p3sq-s15-s34
       P(260) = 4*p2sq-s23-s45+2*P(259)
       FR(10) = -2*(-Dij2345R(1,1)+Dij2345R(1,2)+Dij2345R(2,1)-Dij2345
     -   R(4,2)+4*(-3*EijR(11,2)+EijR(22,3)-EijR(37,4)+EijR(43,4))+Ei
     -   jR(2,3)*P(145)+(-2*EijR(2,1)+EijR(10,3)-EijR(20,3))*P(186)-E
     -   ijR(12,3)*P(188)+(-EijR(9,2)+EijR(10,2))*P(191)-EijR(9,3)*P(
     -   199)-EijR(3,2)*P(214)-EijR(2,2)*P(257)+2*(D02345R+EijR(4,1)*
     -   P(52)+(-EijR(8,3)+EijR(17,3))*P(143)-EijR(3,1)*P(200)-EijR(6
     -   ,2)*P(215)-EijR(5,2)*P(258))+EijR(8,2)*P(260))
       FI(10) = -2*(-Dij2345I(1,1)+Dij2345I(1,2)+Dij2345I(2,1)-Dij2345
     -   I(4,2)+4*(-3*EijI(11,2)+EijI(22,3)-EijI(37,4)+EijI(43,4))+Ei
     -   jI(2,3)*P(145)+(-2*EijI(2,1)+EijI(10,3)-EijI(20,3))*P(186)-E
     -   ijI(12,3)*P(188)+(-EijI(9,2)+EijI(10,2))*P(191)-EijI(9,3)*P(
     -   199)-EijI(3,2)*P(214)-EijI(2,2)*P(257)+2*(D02345I+EijI(4,1)*
     -   P(52)+(-EijI(8,3)+EijI(17,3))*P(143)-EijI(3,1)*P(200)-EijI(6
     -   ,2)*P(215)-EijI(5,2)*P(258))+EijI(8,2)*P(260))
       F(10)=DCMPLX(FR(10),FI(10))
       P(261) = p2sq-3*s12-2*P(220)
       P(262) = p2sq-2*p3sq-s12-s23+s45
       P(263) = p5sq+s15-s23+s45+2*P(84)
       P(264) = p3sq+s12-s34-s45
       P(265) = -3*p2sq+s12+2*P(225)
       FR(11) = 2*(D02345R-Dij2345R(1,2)-Dij2345R(3,1)+Dij2345R(5,2)+4
     -   *(EijR(11,2)-EijR(22,3)+EijR(37,4)-EijR(44,4))-EijR(2,3)*P(1
     -   45)+EijR(15,3)*P(186)+(-EijR(9,3)+EijR(20,3))*P(188)-EijR(4,
     -   2)*P(191)+EijR(10,3)*P(202)+EijR(10,2)*P(214)-EijR(2,1)*P(26
     -   1)+EijR(3,1)*P(262)+EijR(4,1)*P(263)+2*(Dij2345R(1,1)+p2sq*E
     -   ijR(2,2)+(EE0R+EijR(1,1))*P(21)+(EijR(8,3)-EijR(18,3))*P(143
     -   )+EijR(7,2)*P(215)+EijR(5,2)*P(258)-EijR(8,2)*P(264))+EijR(9
     -   ,2)*P(265))
       FI(11) = 2*(D02345I-Dij2345I(1,2)-Dij2345I(3,1)+Dij2345I(5,2)+4
     -   *(EijI(11,2)-EijI(22,3)+EijI(37,4)-EijI(44,4))-EijI(2,3)*P(1
     -   45)+EijI(15,3)*P(186)+(-EijI(9,3)+EijI(20,3))*P(188)-EijI(4,
     -   2)*P(191)+EijI(10,3)*P(202)+EijI(10,2)*P(214)-EijI(2,1)*P(26
     -   1)+EijI(3,1)*P(262)+EijI(4,1)*P(263)+2*(Dij2345I(1,1)+p2sq*E
     -   ijI(2,2)+(EE0I+EijI(1,1))*P(21)+(EijI(8,3)-EijI(18,3))*P(143
     -   )+EijI(7,2)*P(215)+EijI(5,2)*P(258)-EijI(8,2)*P(264))+EijI(9
     -   ,2)*P(265))
       F(11)=DCMPLX(FR(11),FI(11))
       P(266) = s12-s34+2*P(21)
       P(267) = -3*s12+2*P(206)
       P(268) = -5*s12+2*P(206)
       P(269) = s12-s45+3*P(69)
       P(270) = p5sq-s45+3*P(69)
       FR(12) = 2*(Dij2345R(1,2)+s12*EijR(2,2)-4*(EijR(11,2)-EijR(37,4
     -   )+EijR(40,4))+2*(D02345R+Dij2345R(1,1)+EijR(21,3)-EijR(1,2)*
     -   P(21))+EE0R*P(36)-EijR(4,1)*P(52)+EijR(3,1)*P(56)+(EijR(2,3)
     -   -EijR(8,3))*P(145)+(EijR(10,3)-EijR(18,3))*P(213)+(-EijR(9,3
     -   )+EijR(17,3))*P(214)-3*(EijR(7,2)*P(43)-EijR(6,2)*P(264))-Ei
     -   jR(1,1)*P(266)+EijR(2,1)*P(267)+EijR(5,2)*P(268)-EijR(8,2)*P
     -   (269)+EijR(9,2)*P(270))
       FI(12) = 2*(Dij2345I(1,2)+s12*EijI(2,2)-4*(EijI(11,2)-EijI(37,4
     -   )+EijI(40,4))+2*(D02345I+Dij2345I(1,1)+EijI(21,3)-EijI(1,2)*
     -   P(21))+EE0I*P(36)-EijI(4,1)*P(52)+EijI(3,1)*P(56)+(EijI(2,3)
     -   -EijI(8,3))*P(145)+(EijI(10,3)-EijI(18,3))*P(213)+(-EijI(9,3
     -   )+EijI(17,3))*P(214)-3*(EijI(7,2)*P(43)-EijI(6,2)*P(264))-Ei
     -   jI(1,1)*P(266)+EijI(2,1)*P(267)+EijI(5,2)*P(268)-EijI(8,2)*P
     -   (269)+EijI(9,2)*P(270))
       F(12)=DCMPLX(FR(12),FI(12))
       P(271) = -p2sq+s34
       P(272) = p2sq+s12+2*s15-s34
       P(273) = -p5sq+s45
       P(274) = -s12+s34+2*P(21)
       P(275) = s12*P(1)+p2sq*P(274)
       P(276) = -s12-2*s15+s34
       P(277) = 2*s15+s34
       P(278) = p2sq**2+s12*P(276)-p2sq*P(277)
       P(279) = p3sq*p5sq+p5sq*s12-s12*s15-p5sq*s45+s15*s45
       P(280) = s12**2+s12*s23-s12*s34-s23*s34-s12*s45+s34*s45+p2sq*P(
     -   36)+2*P(279)
       P(281) = s12-2*s15-s34
       P(282) = s15+s45
       P(283) = -s12+3*s34+2*P(282)
       P(284) = p3sq*p5sq+p5sq**2
       P(285) = P(16)*P(200)+s45*P(281)+p5sq*P(283)-2*P(284)
       P(286) = 2*p5sq+s12
       P(287) = p5sq-s12-2*s15
       P(288) = p2sq**2+s12*P(286)+2*p2sq*P(287)
       P(289) = s12+s23+2*P(21)
       P(290) = s12+2*P(21)
       P(291) = -2*s12-s23+s45
       P(292) = p2sq**2+2*p3sq*P(21)+s12*P(289)-s45*P(290)+p2sq*P(291)
       P(293) = p3sq*p5sq+p5sq**2-p3sq*s15-p5sq*s15-p5sq*s34+s15*s34
       P(294) = -(p5sq*s12)+s12*s15-s12*s23+p2sq*P(186)+s45*P(290)-2*P
     -   (293)
       P(295) = s12-s23-s34
       P(296) = s15-s45
       P(297) = -s23+s34-2*P(296)
       P(298) = s12**2+p3sq*s15-s15*s45
       P(299) = p2sq**2+s23*s34-s34*s45+p2sq*P(295)+s12*P(297)-2*P(298
     -   )
       P(300) = s12+s15
       P(301) = s34-2*P(300)
       P(302) = s15-s23+2*s45
       P(303) = p3sq*s15-s15*s45
       P(304) = s15*s34+s23*s34-s34*s45+p2sq*P(200)+p5sq*P(301)+s12*P(
     -   302)-2*P(303)
       P(305) = p5sq+s12+s23
       P(306) = p2sq+s15-2*P(305)
       P(307) = p5sq*s12+s45**2
       P(308) = -(p2sq*p5sq)-s12*s15+p5sq*s23+s12*s23+s45*P(306)+2*P(3
     -   07)
       P(309) = p5sq+P(208)
       P(310) = s12+2*P(309)
       P(311) = 3*p2sq-P(310)
       FR(13) = -C0345R+8*EijR(46,4)+P(56)*(Dij2345R(2,1)+EijR(3,2)*P(
     -   215))+D02345R*P(271)-Dij2345R(1,1)*P(272)+(Dij2345R(3,1)+Eij
     -   R(4,2)*P(186))*P(273)-EijR(2,1)*P(275)+EijR(2,2)*P(278)+EijR
     -   (3,1)*P(280)+EijR(4,1)*P(285)-EijR(5,2)*P(288)+EijR(6,2)*P(2
     -   92)+EijR(7,2)*P(294)-EijR(8,2)*P(299)+EijR(9,2)*P(304)+EijR(
     -   10,2)*P(308)+2*(Dij2345R(7,2)+EijR(22,3)*P(145)+P(21)*(-(EE0
     -   R*P(16))+EijR(1,2)*P(143)+EijR(1,1)*P(148))+EijR(24,3)*P(213
     -   )-EijR(23,3)*P(214)+EijR(11,2)*P(311))
       FI(13) = -C0345I+8*EijI(46,4)+P(56)*(Dij2345I(2,1)+EijI(3,2)*P(
     -   215))+D02345I*P(271)-Dij2345I(1,1)*P(272)+(Dij2345I(3,1)+Eij
     -   I(4,2)*P(186))*P(273)-EijI(2,1)*P(275)+EijI(2,2)*P(278)+EijI
     -   (3,1)*P(280)+EijI(4,1)*P(285)-EijI(5,2)*P(288)+EijI(6,2)*P(2
     -   92)+EijI(7,2)*P(294)-EijI(8,2)*P(299)+EijI(9,2)*P(304)+EijI(
     -   10,2)*P(308)+2*(Dij2345I(7,2)+EijI(22,3)*P(145)+P(21)*(-(EE0
     -   I*P(16))+EijI(1,2)*P(143)+EijI(1,1)*P(148))+EijI(24,3)*P(213
     -   )-EijI(23,3)*P(214)+EijI(11,2)*P(311))
       F(13)=DCMPLX(FR(13),FI(13))
       P(312) = s23-s45+2*P(21)
       P(313) = 2*p2sq+s12
       P(314) = s23-3*s45+2*P(21)
       P(315) = 3*P(143)+2*P(220)
       P(316) = 2*P(21)+3*P(143)
       P(317) = 2*p2sq+s12-s23+3*s34+s45
       FR(14) = 2*(Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(4,2)-14*EijR(1
     -   1,2)-10*EijR(22,3)+6*EijR(23,3)+4*(EijR(37,4)-EijR(43,4))+3*
     -   (D02345R-s34*EijR(3,2)+EijR(4,1)*P(52))+2*(EijR(8,3)-EijR(17
     -   ,3))*P(143)-EijR(2,3)*P(145)+(2*EijR(2,1)-EijR(10,3)+EijR(20
     -   ,3))*P(186)+EijR(12,3)*P(188)+EijR(9,2)*P(191)+EijR(9,3)*P(1
     -   99)+(EE0R+EijR(1,1))*P(312)-EijR(2,2)*P(313)-EijR(3,1)*P(314
     -   )+EijR(5,2)*P(315)-EijR(6,2)*P(316)+EijR(8,2)*P(317))
       FI(14) = 2*(Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(4,2)-14*EijI(1
     -   1,2)-10*EijI(22,3)+6*EijI(23,3)+4*(EijI(37,4)-EijI(43,4))+3*
     -   (D02345I-s34*EijI(3,2)+EijI(4,1)*P(52))+2*(EijI(8,3)-EijI(17
     -   ,3))*P(143)-EijI(2,3)*P(145)+(2*EijI(2,1)-EijI(10,3)+EijI(20
     -   ,3))*P(186)+EijI(12,3)*P(188)+EijI(9,2)*P(191)+EijI(9,3)*P(1
     -   99)+(EE0I+EijI(1,1))*P(312)-EijI(2,2)*P(313)-EijI(3,1)*P(314
     -   )+EijI(5,2)*P(315)-EijI(6,2)*P(316)+EijI(8,2)*P(317))
       F(14)=DCMPLX(FR(14),FI(14))
       P(318) = s15-s23+s34
       P(319) = p3sq-P(318)
       P(320) = p5sq-s45+2*P(319)
       P(321) = p3sq-s15+s23-s34
       P(322) = -p5sq+s15-s23+s45
       FR(15) = 2*(EE0R*s45+EijR(4,1)*P(52)-EijR(3,1)*P(104)-EijR(1,2)
     -   *P(193)-EijR(1,1)*P(232)+EijR(2,1)*P(233)+EijR(5,2)*P(236)+E
     -   ijR(6,2)*P(239)+EijR(7,2)*P(320)-2*(D02345R+Dij2345R(1,2)+s1
     -   2*EijR(2,3)-p2sq*EijR(8,3)-EijR(11,2)+2*(Dij2345R(1,1)-EijR(
     -   22,3)-EijR(37,4)+EijR(40,4))+EijR(5,3)*P(143)-EijR(2,2)*P(14
     -   9)-EijR(18,3)*P(186)-EijR(8,2)*P(241)-(-EijR(9,3)+EijR(17,3)
     -   )*P(245)+EijR(9,2)*P(321)-EijR(10,3)*P(322)))
       FI(15) = 2*(EE0I*s45+EijI(4,1)*P(52)-EijI(3,1)*P(104)-EijI(1,2)
     -   *P(193)-EijI(1,1)*P(232)+EijI(2,1)*P(233)+EijI(5,2)*P(236)+E
     -   ijI(6,2)*P(239)+EijI(7,2)*P(320)-2*(D02345I+Dij2345I(1,2)+s1
     -   2*EijI(2,3)-p2sq*EijI(8,3)-EijI(11,2)+2*(Dij2345I(1,1)-EijI(
     -   22,3)-EijI(37,4)+EijI(40,4))+EijI(5,3)*P(143)-EijI(2,2)*P(14
     -   9)-EijI(18,3)*P(186)-EijI(8,2)*P(241)-(-EijI(9,3)+EijI(17,3)
     -   )*P(245)+EijI(9,2)*P(321)-EijI(10,3)*P(322)))
       F(15)=DCMPLX(FR(15),FI(15))
       P(323) = -s12+2*P(220)
       P(324) = 2*p2sq+3*s12
       P(325) = -3*s34+2*P(56)
       P(326) = s23+3*s34+2*s45
       P(327) = 3*s12+2*s15-P(326)
       P(328) = -p2sq+s12+s23+5*s34-s45
       P(329) = p2sq-s23
       P(330) = -p5sq+s15+s45+2*P(16)+3*P(329)
       FR(16) = 2*(-Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(5,2)-3*s34*Ei
     -   jR(10,2)-12*EijR(11,2)-6*(EijR(22,3)-EijR(24,3))+4*(EijR(37,
     -   4)-EijR(44,4))+2*(D02345R+(EijR(8,3)-EijR(18,3))*P(143))-Eij
     -   R(2,3)*P(145)+EijR(15,3)*P(186)+(-EijR(9,3)+EijR(20,3))*P(18
     -   8)+EijR(10,3)*P(202)+(EE0R+EijR(1,1))*P(312)+EijR(5,2)*P(315
     -   )-EijR(7,2)*P(316)+EijR(2,1)*P(323)-EijR(2,2)*P(324)-EijR(3,
     -   1)*P(325)+EijR(4,1)*P(327)+EijR(8,2)*P(328)+EijR(9,2)*P(330)
     -   )
       FI(16) = 2*(-Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(5,2)-3*s34*Ei
     -   jI(10,2)-12*EijI(11,2)-6*(EijI(22,3)-EijI(24,3))+4*(EijI(37,
     -   4)-EijI(44,4))+2*(D02345I+(EijI(8,3)-EijI(18,3))*P(143))-Eij
     -   I(2,3)*P(145)+EijI(15,3)*P(186)+(-EijI(9,3)+EijI(20,3))*P(18
     -   8)+EijI(10,3)*P(202)+(EE0I+EijI(1,1))*P(312)+EijI(5,2)*P(315
     -   )-EijI(7,2)*P(316)+EijI(2,1)*P(323)-EijI(2,2)*P(324)-EijI(3,
     -   1)*P(325)+EijI(4,1)*P(327)+EijI(8,2)*P(328)+EijI(9,2)*P(330)
     -   )
       F(16)=DCMPLX(FR(16),FI(16))
       P(331) = s12+2*P(220)
       P(332) = -3*s34+2*P(40)
       P(333) = p2sq-s12+2*P(220)
       P(334) = p2sq+2*p3sq+s12-s23-5*s34-s45
       FR(17) = 2*(Dij2345R(1,1)-Dij2345R(1,2)-s12*EijR(2,2)-6*EijR(22
     -   ,3)-4*(EijR(11,2)-EijR(37,4))+EijR(4,1)*P(115)+2*(D02345R+Ei
     -   jR(8,3)*P(143))-EijR(2,3)*P(145)-EijR(10,3)*P(186)-EijR(9,3)
     -   *P(188)+EijR(9,2)*P(191)+(EE0R+EijR(1,1))*P(312)+EijR(2,1)*P
     -   (331)-EijR(3,1)*P(332)+EijR(5,2)*P(333)-EijR(8,2)*P(334))
       FI(17) = 2*(Dij2345I(1,1)-Dij2345I(1,2)-s12*EijI(2,2)-6*EijI(22
     -   ,3)-4*(EijI(11,2)-EijI(37,4))+EijI(4,1)*P(115)+2*(D02345I+Ei
     -   jI(8,3)*P(143))-EijI(2,3)*P(145)-EijI(10,3)*P(186)-EijI(9,3)
     -   *P(188)+EijI(9,2)*P(191)+(EE0I+EijI(1,1))*P(312)+EijI(2,1)*P
     -   (331)-EijI(3,1)*P(332)+EijI(5,2)*P(333)-EijI(8,2)*P(334))
       F(17)=DCMPLX(FR(17),FI(17))
       P(335) = 2*p3sq-3*s45+4*P(16)
       P(336) = p2sq-s34
       P(337) = p3sq+s12-s23+2*P(336)
       P(338) = -s23+s45+2*P(143)
       P(339) = -3*s45+4*P(16)
       P(340) = 4*p2sq+p3sq-s12+s45-2*P(237)
       P(341) = p3sq-p5sq+s12-s23+s45+2*P(336)
       P(342) = s23+s34-s45
       P(343) = -4*p2sq-p3sq+p5sq+s12+2*P(342)
       FR(18) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(3,1)-Dij2345R(
     -   5,2)+2*(p3sq*EijR(6,2)-p3sq*EijR(8,2))-6*(EijR(21,3)-EijR(22
     -   ,3))+4*(EijR(37,4)-EijR(40,4)+EijR(42,4)-EijR(44,4))+(-EijR(
     -   9,3)+EijR(17,3)-EijR(19,3)+EijR(20,3))*P(40)+EijR(14,3)*P(52
     -   )+EijR(15,3)*P(273)+(EijR(2,2)-EijR(5,2))*P(335)+EijR(2,3)*P
     -   (337)+(EijR(5,3)-EijR(7,3))*P(338)+(EijR(7,2)-EijR(9,2))*P(3
     -   39)-EijR(8,3)*P(340)-EijR(10,3)*P(341)-EijR(18,3)*P(343))
       FI(18) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(3,1)-Dij2345I(
     -   5,2)+2*(p3sq*EijI(6,2)-p3sq*EijI(8,2))-6*(EijI(21,3)-EijI(22
     -   ,3))+4*(EijI(37,4)-EijI(40,4)+EijI(42,4)-EijI(44,4))+(-EijI(
     -   9,3)+EijI(17,3)-EijI(19,3)+EijI(20,3))*P(40)+EijI(14,3)*P(52
     -   )+EijI(15,3)*P(273)+(EijI(2,2)-EijI(5,2))*P(335)+EijI(2,3)*P
     -   (337)+(EijI(5,3)-EijI(7,3))*P(338)+(EijI(7,2)-EijI(9,2))*P(3
     -   39)-EijI(8,3)*P(340)-EijI(10,3)*P(341)-EijI(18,3)*P(343))
       F(18)=DCMPLX(FR(18),FI(18))
       P(344) = -p3sq+s15+s34+s45
       P(345) = s15+s45-2*P(69)
       P(346) = 3*p3sq+2*P(238)
       P(347) = p2sq+3*p3sq-s23-s34
       FR(19) = -4*(Dij2345R(1,1)-Dij2345R(2,1)+p5sq*(EijR(2,1)-EijR(3
     -   ,1))+p3sq*EijR(3,3)-2*(3*(EijR(22,3)-EijR(23,3))+EijR(37,4)+
     -   EijR(38,4)-2*EijR(43,4))+(EijR(5,2)-EijR(6,2))*P(21)+(EijR(9
     -   ,2)-EijR(10,2))*P(52)-EijR(3,2)*P(69)-EijR(2,3)*P(203)+(EijR
     -   (8,3)+EijR(11,3)-2*EijR(17,3))*P(215)+EijR(2,2)*P(344)-EijR(
     -   8,2)*P(345)+EijR(9,3)*P(346)-EijR(12,3)*P(347))
       FI(19) = -4*(Dij2345I(1,1)-Dij2345I(2,1)+p5sq*(EijI(2,1)-EijI(3
     -   ,1))+p3sq*EijI(3,3)-2*(3*(EijI(22,3)-EijI(23,3))+EijI(37,4)+
     -   EijI(38,4)-2*EijI(43,4))+(EijI(5,2)-EijI(6,2))*P(21)+(EijI(9
     -   ,2)-EijI(10,2))*P(52)-EijI(3,2)*P(69)-EijI(2,3)*P(203)+(EijI
     -   (8,3)+EijI(11,3)-2*EijI(17,3))*P(215)+EijI(2,2)*P(344)-EijI(
     -   8,2)*P(345)+EijI(9,3)*P(346)-EijI(12,3)*P(347))
       F(19)=DCMPLX(FR(19),FI(19))
       P(348) = p3sq+s12-s15-s34-s45
       P(349) = 2*p3sq+s12-s34-s45
       P(350) = -s12+s15+s34+s45
       P(351) = p2sq+2*p3sq-s23-s34
       P(352) = s12-s34-s45
       FR(20) = 4*(p5sq*(-EijR(2,1)+EijR(4,1))+p3sq*(EijR(3,2)+EijR(12
     -   ,3)-EijR(13,3))+2*(2*EijR(22,3)-EijR(23,3)-EijR(24,3)+EijR(3
     -   7,4)-EijR(43,4)-EijR(44,4)+EijR(45,4))+(-EijR(5,2)+EijR(7,2)
     -   )*P(21)+(EijR(2,3)-EijR(10,3))*P(203)+(-EijR(8,3)+EijR(17,3)
     -   +EijR(18,3)-EijR(19,3))*P(215)+EijR(2,2)*P(348)-EijR(8,2)*P(
     -   349)+EijR(9,2)*P(350)+(-EijR(9,3)+EijR(20,3))*P(351)+EijR(10
     -   ,2)*P(352))
       FI(20) = 4*(p5sq*(-EijI(2,1)+EijI(4,1))+p3sq*(EijI(3,2)+EijI(12
     -   ,3)-EijI(13,3))+2*(2*EijI(22,3)-EijI(23,3)-EijI(24,3)+EijI(3
     -   7,4)-EijI(43,4)-EijI(44,4)+EijI(45,4))+(-EijI(5,2)+EijI(7,2)
     -   )*P(21)+(EijI(2,3)-EijI(10,3))*P(203)+(-EijI(8,3)+EijI(17,3)
     -   +EijI(18,3)-EijI(19,3))*P(215)+EijI(2,2)*P(348)-EijI(8,2)*P(
     -   349)+EijI(9,2)*P(350)+(-EijI(9,3)+EijI(20,3))*P(351)+EijI(10
     -   ,2)*P(352))
       F(20)=DCMPLX(FR(20),FI(20))
       P(353) = s12+s23-s45
       P(354) = -p3sq+s12+s23+s34
       P(355) = s12-s23-s34+s45
       P(356) = -p3sq+s12+2*s23+s34-s45
       P(357) = p5sq-s12-s23
       P(358) = -p5sq+s12+s23+s45-2*P(69)
       P(359) = p3sq-P(237)
       P(360) = p5sq-s12+2*P(359)
       FR(21) = -4*(Dij2345R(1,2)-Dij2345R(5,2)-2*(EijR(11,2)+EijR(21,
     -   3)-2*EijR(22,3)+EijR(24,3)+EijR(37,4)-EijR(40,4)+EijR(42,4)-
     -   EijR(44,4))+(-EijR(2,1)+EijR(4,1))*P(16)+(EijR(14,3)-EijR(15
     -   ,3))*P(43)+EijR(4,2)*P(52)+(EijR(8,2)-EijR(9,3)-EijR(10,2)+E
     -   ijR(17,3)-EijR(19,3)+EijR(20,3))*P(56)+(EijR(5,3)-EijR(7,3))
     -   *P(193)-EijR(2,2)*P(353)+EijR(2,3)*P(354)+(-EijR(5,2)+EijR(7
     -   ,2))*P(355)-EijR(8,3)*P(356)-EijR(9,2)*P(357)-EijR(10,3)*P(3
     -   58)-EijR(18,3)*P(360))
       FI(21) = -4*(Dij2345I(1,2)-Dij2345I(5,2)-2*(EijI(11,2)+EijI(21,
     -   3)-2*EijI(22,3)+EijI(24,3)+EijI(37,4)-EijI(40,4)+EijI(42,4)-
     -   EijI(44,4))+(-EijI(2,1)+EijI(4,1))*P(16)+(EijI(14,3)-EijI(15
     -   ,3))*P(43)+EijI(4,2)*P(52)+(EijI(8,2)-EijI(9,3)-EijI(10,2)+E
     -   ijI(17,3)-EijI(19,3)+EijI(20,3))*P(56)+(EijI(5,3)-EijI(7,3))
     -   *P(193)-EijI(2,2)*P(353)+EijI(2,3)*P(354)+(-EijI(5,2)+EijI(7
     -   ,2))*P(355)-EijI(8,3)*P(356)-EijI(9,2)*P(357)-EijI(10,3)*P(3
     -   58)-EijI(18,3)*P(360))
       F(21)=DCMPLX(FR(21),FI(21))
       P(361) = p3sq-s34-s45
       P(362) = p3sq+p5sq-s34-2*s45
       FR(22) = 4*(-Dij2345R(1,1)+Dij2345R(3,1)+p3sq*(EijR(12,3)-EijR(
     -   13,3))+2*(2*(EijR(22,3)-EijR(24,3))+EijR(37,4)-EijR(43,4)-Ei
     -   jR(44,4)+EijR(45,4))+EijR(4,2)*P(52)+(-EijR(8,2)+EijR(10,2))
     -   *P(69)+(EijR(2,3)-EijR(10,3))*P(203)+(-EijR(8,3)+EijR(17,3)+
     -   EijR(18,3)-EijR(19,3))*P(215)+(-EijR(9,3)+EijR(20,3))*P(351)
     -   +EijR(2,2)*P(361)-EijR(9,2)*P(362))
       FI(22) = 4*(-Dij2345I(1,1)+Dij2345I(3,1)+p3sq*(EijI(12,3)-EijI(
     -   13,3))+2*(2*(EijI(22,3)-EijI(24,3))+EijI(37,4)-EijI(43,4)-Ei
     -   jI(44,4)+EijI(45,4))+EijI(4,2)*P(52)+(-EijI(8,2)+EijI(10,2))
     -   *P(69)+(EijI(2,3)-EijI(10,3))*P(203)+(-EijI(8,3)+EijI(17,3)+
     -   EijI(18,3)-EijI(19,3))*P(215)+(-EijI(9,3)+EijI(20,3))*P(351)
     -   +EijI(2,2)*P(361)-EijI(9,2)*P(362))
       F(22)=DCMPLX(FR(22),FI(22))
       P(363) = p3sq+2*P(352)
       FR(23) = -4*(p3sq*(EijR(8,2)+EijR(9,3)-EijR(10,2)+EijR(16,3)-2*
     -   EijR(20,3))-2*(EijR(22,3)-EijR(24,3)+EijR(37,4)+EijR(39,4)-2
     -   *EijR(44,4))-(EijR(2,3)-2*EijR(10,3)+EijR(15,3))*P(203)+(Eij
     -   R(8,3)+EijR(14,3)-2*EijR(18,3))*P(215)-EijR(2,2)*P(264)-EijR
     -   (4,2)*P(352)+EijR(9,2)*P(363))
       FI(23) = -4*(p3sq*(EijI(8,2)+EijI(9,3)-EijI(10,2)+EijI(16,3)-2*
     -   EijI(20,3))-2*(EijI(22,3)-EijI(24,3)+EijI(37,4)+EijI(39,4)-2
     -   *EijI(44,4))-(EijI(2,3)-2*EijI(10,3)+EijI(15,3))*P(203)+(Eij
     -   I(8,3)+EijI(14,3)-2*EijI(18,3))*P(215)-EijI(2,2)*P(264)-EijI
     -   (4,2)*P(352)+EijI(9,2)*P(363))
       F(23)=DCMPLX(FR(23),FI(23))
       P(364) = p2sq-p3sq+s34+s45
       P(365) = p2sq-p3sq-p5sq+s34+2*s45
       FR(24) = 4*(-Dij2345R(1,1)+Dij2345R(3,1)-s12*EijR(2,1)+p3sq*(Ei
     -   jR(3,1)-EijR(9,3)+EijR(20,3))-2*(2*EijR(11,2)-EijR(22,3)+Eij
     -   R(24,3)-EijR(37,4)+EijR(44,4))+EijR(4,2)*P(52)+(EijR(8,2)-Ei
     -   jR(10,2))*P(56)-EijR(4,1)*P(84)+(EijR(5,2)-EijR(7,2))*P(143)
     -   +(EijR(2,3)-EijR(10,3))*P(203)+(-EijR(8,3)+EijR(18,3))*P(215
     -   )-EijR(2,2)*P(364)+EijR(9,2)*P(365))
       FI(24) = 4*(-Dij2345I(1,1)+Dij2345I(3,1)-s12*EijI(2,1)+p3sq*(Ei
     -   jI(3,1)-EijI(9,3)+EijI(20,3))-2*(2*EijI(11,2)-EijI(22,3)+Eij
     -   I(24,3)-EijI(37,4)+EijI(44,4))+EijI(4,2)*P(52)+(EijI(8,2)-Ei
     -   jI(10,2))*P(56)-EijI(4,1)*P(84)+(EijI(5,2)-EijI(7,2))*P(143)
     -   +(EijI(2,3)-EijI(10,3))*P(203)+(-EijI(8,3)+EijI(18,3))*P(215
     -   )-EijI(2,2)*P(364)+EijI(9,2)*P(365))
       F(24)=DCMPLX(FR(24),FI(24))
       P(366) = p3sq-s45+3*P(16)
       P(367) = p2sq+p3sq
       P(368) = 4*s12-s23-5*s34-2*s45+3*P(367)
       P(369) = -s45+3*P(16)
       P(370) = -s23+s45+3*P(143)
       P(371) = 3*p3sq+2*P(56)
       P(372) = -3*p2sq-p3sq+2*p5sq-4*s12+s23+5*s34
       FR(25) = 2*(-(p3sq*(EijR(3,1)+2*(EijR(9,3)-EijR(20,3))))+EijR(1
     -   0,2)*P(114)+2*(Dij2345R(1,1)-Dij2345R(3,1)+3*EijR(11,2)+EijR
     -   (22,3)+EijR(24,3)+2*(EijR(37,4)-EijR(44,4))-EijR(4,2)*P(52)+
     -   (EijR(2,3)-EijR(10,3))*P(203)+(-EijR(8,3)+EijR(18,3))*P(215)
     -   )+EijR(2,1)*P(366)+EijR(2,2)*P(368)-EijR(4,1)*P(369)+(-EijR(
     -   5,2)+EijR(7,2))*P(370)-EijR(8,2)*P(371)+EijR(9,2)*P(372))
       FI(25) = 2*(-(p3sq*(EijI(3,1)+2*(EijI(9,3)-EijI(20,3))))+EijI(1
     -   0,2)*P(114)+2*(Dij2345I(1,1)-Dij2345I(3,1)+3*EijI(11,2)+EijI
     -   (22,3)+EijI(24,3)+2*(EijI(37,4)-EijI(44,4))-EijI(4,2)*P(52)+
     -   (EijI(2,3)-EijI(10,3))*P(203)+(-EijI(8,3)+EijI(18,3))*P(215)
     -   )+EijI(2,1)*P(366)+EijI(2,2)*P(368)-EijI(4,1)*P(369)+(-EijI(
     -   5,2)+EijI(7,2))*P(370)-EijI(8,2)*P(371)+EijI(9,2)*P(372))
       F(25)=DCMPLX(FR(25),FI(25))
       P(373) = -3*s45+2*P(45)
       P(374) = s45-2*P(45)
       P(375) = -2*p2sq-p3sq+s23+s34
       FR(26) = -2*(Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(2,1)+Dij2345R
     -   (4,2)+4*(2*(EijR(21,3)-EijR(22,3))-EijR(37,4)+EijR(40,4)-Eij
     -   R(41,4)+EijR(43,4))+(EijR(11,3)-EijR(12,3))*P(40)-(2*(EijR(7
     -   ,2)-EijR(9,2))+EijR(10,3)-EijR(18,3)+EijR(19,3))*P(52)+EijR(
     -   9,3)*P(212)-EijR(20,3)*P(273)-EijR(2,3)*P(337)+(-EijR(5,3)+E
     -   ijR(6,3))*P(338)+EijR(8,3)*P(340)+(-EijR(2,2)+EijR(5,2))*P(3
     -   73)+(EijR(6,2)-EijR(8,2))*P(374)+2*(D02345R+EijR(17,3)*P(375
     -   )))
       FI(26) = -2*(Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(2,1)+Dij2345I
     -   (4,2)+4*(2*(EijI(21,3)-EijI(22,3))-EijI(37,4)+EijI(40,4)-Eij
     -   I(41,4)+EijI(43,4))+(EijI(11,3)-EijI(12,3))*P(40)-(2*(EijI(7
     -   ,2)-EijI(9,2))+EijI(10,3)-EijI(18,3)+EijI(19,3))*P(52)+EijI(
     -   9,3)*P(212)-EijI(20,3)*P(273)-EijI(2,3)*P(337)+(-EijI(5,3)+E
     -   ijI(6,3))*P(338)+EijI(8,3)*P(340)+(-EijI(2,2)+EijI(5,2))*P(3
     -   73)+(EijI(6,2)-EijI(8,2))*P(374)+2*(D02345I+EijI(17,3)*P(375
     -   )))
       F(26)=DCMPLX(FR(26),FI(26))
       P(376) = -p5sq+s12+s15-s34
       P(377) = p2sq+s12-2*P(195)
       P(378) = -2*s15+s23+s34-s45
       P(379) = p5sq+P(378)
       P(380) = p2sq-3*s12+2*P(379)
       P(381) = p2sq-3*s12+s23+2*s34-s45
       P(382) = p2sq+s23+3*P(56)
       P(383) = -p3sq+2*s12+s23+s34-s45
       P(384) = p3sq-s34-2*P(353)
       FR(27) = 2*(-D02345R-Dij2345R(1,1)+6*EijR(11,2)-4*(EijR(22,3)-E
     -   ijR(23,3)-EijR(37,4)+EijR(40,4)-EijR(41,4)+EijR(43,4))+EijR(
     -   7,2)*P(186)+EijR(9,2)*P(191)+EijR(2,2)*P(377)-EijR(5,2)*P(38
     -   0)+EijR(6,2)*P(381)-EijR(8,2)*P(382)-2*(Dij2345R(1,2)-Dij234
     -   5R(4,2)+EijR(3,1)*P(16)-(EijR(1,1)+EijR(1,2))*P(21)-(-EijR(1
     -   0,3)+EijR(18,3)-EijR(19,3)+EijR(20,3))*P(43)+EijR(10,2)*P(52
     -   )-(EijR(3,2)+EijR(11,3)-EijR(12,3))*P(56)-(-EijR(5,3)+EijR(6
     -   ,3))*P(193)+EijR(2,3)*P(354)-EijR(8,3)*P(356)-EijR(2,1)*P(37
     -   6)-EijR(9,3)*P(383)-EijR(17,3)*P(384)))
       FI(27) = 2*(-D02345I-Dij2345I(1,1)+6*EijI(11,2)-4*(EijI(22,3)-E
     -   ijI(23,3)-EijI(37,4)+EijI(40,4)-EijI(41,4)+EijI(43,4))+EijI(
     -   7,2)*P(186)+EijI(9,2)*P(191)+EijI(2,2)*P(377)-EijI(5,2)*P(38
     -   0)+EijI(6,2)*P(381)-EijI(8,2)*P(382)-2*(Dij2345I(1,2)-Dij234
     -   5I(4,2)+EijI(3,1)*P(16)-(EijI(1,1)+EijI(1,2))*P(21)-(-EijI(1
     -   0,3)+EijI(18,3)-EijI(19,3)+EijI(20,3))*P(43)+EijI(10,2)*P(52
     -   )-(EijI(3,2)+EijI(11,3)-EijI(12,3))*P(56)-(-EijI(5,3)+EijI(6
     -   ,3))*P(193)+EijI(2,3)*P(354)-EijI(8,3)*P(356)-EijI(2,1)*P(37
     -   6)-EijI(9,3)*P(383)-EijI(17,3)*P(384)))
       F(27)=DCMPLX(FR(27),FI(27))
       P(385) = p3sq+s23-2*s34+3*P(56)
       P(386) = p3sq-s12-s23-2*s34
       P(387) = s23-s34
       P(388) = p3sq+3*s12-5*s45+2*P(387)
       P(389) = s34-s45
       P(390) = p3sq-s12-3*s23-2*P(389)
       P(391) = -3*s23-4*s34+s45+2*P(84)
       P(392) = -p5sq+s45-2*P(69)
       FR(28) = -2*(Dij2345R(1,1)+Dij2345R(1,2)+s45*(EijR(1,1)-EijR(2,
     -   1))-4*(EijR(21,3)-EijR(22,3)+EijR(36,4)+EijR(37,4)-2*EijR(40
     -   ,4))-(EijR(6,2)+EijR(6,3)-EijR(8,2)+EijR(9,3)-2*EijR(17,3))*
     -   P(40)+(EijR(7,2)-EijR(9,2))*P(52)+(EijR(7,3)-2*EijR(18,3))*P
     -   (97)-EijR(1,3)*P(193)-EijR(1,2)*P(232)-EijR(2,2)*P(385)-EijR
     -   (2,3)*P(386)+EijR(5,2)*P(388)-EijR(5,3)*P(390)+EijR(8,3)*P(3
     -   91)-EijR(10,3)*P(392))
       FI(28) = -2*(Dij2345I(1,1)+Dij2345I(1,2)+s45*(EijI(1,1)-EijI(2,
     -   1))-4*(EijI(21,3)-EijI(22,3)+EijI(36,4)+EijI(37,4)-2*EijI(40
     -   ,4))-(EijI(6,2)+EijI(6,3)-EijI(8,2)+EijI(9,3)-2*EijI(17,3))*
     -   P(40)+(EijI(7,2)-EijI(9,2))*P(52)+(EijI(7,3)-2*EijI(18,3))*P
     -   (97)-EijI(1,3)*P(193)-EijI(1,2)*P(232)-EijI(2,2)*P(385)-EijI
     -   (2,3)*P(386)+EijI(5,2)*P(388)-EijI(5,3)*P(390)+EijI(8,3)*P(3
     -   91)-EijI(10,3)*P(392))
       F(28)=DCMPLX(FR(28),FI(28))
       FR(29) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-6*(EijR(21,3)-EijR(22,3
     -   ))+4*(EijR(37,4)-EijR(40,4))+(-EijR(9,3)+EijR(17,3))*P(40)+(
     -   EijR(10,3)-EijR(18,3))*P(52)+2*(p3sq*EijR(6,2)-p3sq*EijR(8,2
     -   )+(-EijR(7,2)+EijR(9,2))*P(69))+(EijR(2,2)-EijR(5,2))*P(335)
     -   +EijR(2,3)*P(337)+EijR(5,3)*P(338)-EijR(8,3)*P(340))
       FI(29) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-6*(EijI(21,3)-EijI(22,3
     -   ))+4*(EijI(37,4)-EijI(40,4))+(-EijI(9,3)+EijI(17,3))*P(40)+(
     -   EijI(10,3)-EijI(18,3))*P(52)+2*(p3sq*EijI(6,2)-p3sq*EijI(8,2
     -   )+(-EijI(7,2)+EijI(9,2))*P(69))+(EijI(2,2)-EijI(5,2))*P(335)
     -   +EijI(2,3)*P(337)+EijI(5,3)*P(338)-EijI(8,3)*P(340))
       F(29)=DCMPLX(FR(29),FI(29))
       P(393) = p2sq-p3sq+s15+s34+s45
       P(394) = p2sq-p5sq-s12+s15
       P(395) = p2sq-p3sq+s12+s34
       FR(30) = 4*(D02345R-Dij2345R(1,1)+Dij2345R(2,1)-p5sq*EijR(2,1)+
     -   s45*EijR(3,1)-6*EijR(11,2)+p3sq*EijR(12,3)+4*EijR(22,3)-2*(E
     -   ijR(23,3)-EijR(37,4)+EijR(43,4))+(EijR(4,1)-EijR(9,2)+EijR(1
     -   0,2))*P(52)-EijR(3,2)*P(56)-EijR(6,2)*P(143)+EijR(2,3)*P(203
     -   )+(-EijR(8,3)+EijR(17,3))*P(215)-EijR(9,3)*P(351)-EijR(2,2)*
     -   P(393)+EijR(5,2)*P(394)+EijR(8,2)*P(395))
       FI(30) = 4*(D02345I-Dij2345I(1,1)+Dij2345I(2,1)-p5sq*EijI(2,1)+
     -   s45*EijI(3,1)-6*EijI(11,2)+p3sq*EijI(12,3)+4*EijI(22,3)-2*(E
     -   ijI(23,3)-EijI(37,4)+EijI(43,4))+(EijI(4,1)-EijI(9,2)+EijI(1
     -   0,2))*P(52)-EijI(3,2)*P(56)-EijI(6,2)*P(143)+EijI(2,3)*P(203
     -   )+(-EijI(8,3)+EijI(17,3))*P(215)-EijI(9,3)*P(351)-EijI(2,2)*
     -   P(393)+EijI(5,2)*P(394)+EijI(8,2)*P(395))
       F(30)=DCMPLX(FR(30),FI(30))
       P(396) = p3sq+2*s12-3*s34-s45
       P(397) = -s23-5*s34+2*P(56)+3*P(367)
       P(398) = p3sq+2*s12-3*s34
       P(399) = 3*p2sq-s23-5*s34-2*s45+4*P(124)
       FR(31) = -2*(D02345R+(EijR(4,1)+2*EijR(10,2))*P(52)-EijR(3,2)*P
     -   (114)-2*(-Dij2345R(2,1)+4*EijR(11,2)+p3sq*EijR(12,3)+3*EijR(
     -   22,3)+EijR(23,3)+2*(EijR(37,4)-EijR(43,4))+EijR(2,3)*P(203)-
     -   (EijR(8,3)-EijR(17,3))*P(215)-EijR(9,3)*P(351))+(EijR(5,2)-E
     -   ijR(6,2))*P(370)-EijR(2,1)*P(396)-EijR(2,2)*P(397)+EijR(3,1)
     -   *P(398)+EijR(8,2)*P(399))
       FI(31) = -2*(D02345I+(EijI(4,1)+2*EijI(10,2))*P(52)-EijI(3,2)*P
     -   (114)-2*(-Dij2345I(2,1)+4*EijI(11,2)+p3sq*EijI(12,3)+3*EijI(
     -   22,3)+EijI(23,3)+2*(EijI(37,4)-EijI(43,4))+EijI(2,3)*P(203)-
     -   (EijI(8,3)-EijI(17,3))*P(215)-EijI(9,3)*P(351))+(EijI(5,2)-E
     -   ijI(6,2))*P(370)-EijI(2,1)*P(396)-EijI(2,2)*P(397)+EijI(3,1)
     -   *P(398)+EijI(8,2)*P(399))
       F(31)=DCMPLX(FR(31),FI(31))
       P(400) = 3*s12+s23-2*s34-s45
       P(401) = p3sq+s12+s23-3*s34-2*s45
       P(402) = p3sq-3*s12-s34-s45
       P(403) = 2*p3sq+s12-s45
       P(404) = 2*p3sq-s12+s45
       P(405) = p5sq-s45-3*P(69)
       FR(32) = 2*(-Dij2345R(1,1)+6*EijR(21,3)-10*EijR(22,3)+4*(EijR(3
     -   7,4)-EijR(40,4))+EE0R*P(16)+EijR(4,1)*P(52)-(EijR(3,1)-2*(Ei
     -   jR(9,3)-EijR(17,3)))*P(56)-(EijR(1,2)+2*EijR(5,3))*P(193)+Ei
     -   jR(7,2)*P(270)+EijR(1,1)*P(355)-2*(Dij2345R(1,2)-EijR(11,2)-
     -   EijR(10,3)*P(42)-EijR(18,3)*P(43)+EijR(2,3)*P(354)-EijR(8,3)
     -   *P(356))+EijR(2,1)*P(400)+EijR(2,2)*P(401)-EijR(5,2)*P(402)-
     -   EijR(6,2)*P(403)+EijR(8,2)*P(404)+EijR(9,2)*P(405))
       FI(32) = 2*(-Dij2345I(1,1)+6*EijI(21,3)-10*EijI(22,3)+4*(EijI(3
     -   7,4)-EijI(40,4))+EE0I*P(16)+EijI(4,1)*P(52)-(EijI(3,1)-2*(Ei
     -   jI(9,3)-EijI(17,3)))*P(56)-(EijI(1,2)+2*EijI(5,3))*P(193)+Ei
     -   jI(7,2)*P(270)+EijI(1,1)*P(355)-2*(Dij2345I(1,2)-EijI(11,2)-
     -   EijI(10,3)*P(42)-EijI(18,3)*P(43)+EijI(2,3)*P(354)-EijI(8,3)
     -   *P(356))+EijI(2,1)*P(400)+EijI(2,2)*P(401)-EijI(5,2)*P(402)-
     -   EijI(6,2)*P(403)+EijI(8,2)*P(404)+EijI(9,2)*P(405))
       F(32)=DCMPLX(FR(32),FI(32))
       P(406) = -p2sq+s12+s23-s45
       P(407) = s12-3*s34-2*P(282)
       P(408) = p3sq*p5sq+p5sq**2+s15*s45
       P(409) = -(s12*s15)+s12*s23+s15*s34-s23*s34-s12*s45+s34*s45+p5s
     -   q*P(407)+2*P(408)
       P(410) = p5sq**2-p5sq*s15-p5sq*s34+s15*s34+p3sq*P(21)
       P(411) = p5sq*s12-s12*s15+s12*s23-s45*P(290)+p2sq*P(322)+2*P(41
     -   0)
       P(412) = p2sq*p5sq+s12*s15-p5sq*s23-s12*s23-s45*P(306)-2*P(307)
       P(413) = s15+s23+s34-s45
       P(414) = p5sq+P(413)
       P(415) = -3*p2sq+s12+2*P(414)
       P(416) = p2sq-3*s12+2*P(359)
       P(417) = p2sq-s23-3*P(56)
       P(418) = -s15+s23+3*P(52)+2*P(69)
       FR(33) = C0345R-6*Dij2345R(7,2)+8*EijR(46,4)+Dij2345R(2,1)*P(91
     -   )+P(52)*(Dij2345R(3,1)+EijR(4,2)*P(186))+4*EijR(21,3)*P(193)
     -   +Dij2345R(1,1)*P(272)+EijR(2,1)*P(275)-EijR(2,2)*P(278)-EijR
     -   (3,1)*P(280)+EijR(5,2)*P(288)-EijR(6,2)*P(292)+EijR(8,2)*P(2
     -   99)-EijR(9,2)*P(304)+D02345R*P(336)+EijR(3,2)*P(56)*P(406)+E
     -   ijR(4,1)*P(409)+EijR(7,2)*P(411)+EijR(10,2)*P(412)+2*(P(21)*
     -   (EE0R*P(16)-EijR(1,2)*P(143)-EijR(1,1)*P(148))+EijR(11,2)*P(
     -   415)+EijR(22,3)*P(416)-EijR(23,3)*P(417)-EijR(24,3)*P(418))
       FI(33) = C0345I-6*Dij2345I(7,2)+8*EijI(46,4)+Dij2345I(2,1)*P(91
     -   )+P(52)*(Dij2345I(3,1)+EijI(4,2)*P(186))+4*EijI(21,3)*P(193)
     -   +Dij2345I(1,1)*P(272)+EijI(2,1)*P(275)-EijI(2,2)*P(278)-EijI
     -   (3,1)*P(280)+EijI(5,2)*P(288)-EijI(6,2)*P(292)+EijI(8,2)*P(2
     -   99)-EijI(9,2)*P(304)+D02345I*P(336)+EijI(3,2)*P(56)*P(406)+E
     -   ijI(4,1)*P(409)+EijI(7,2)*P(411)+EijI(10,2)*P(412)+2*(P(21)*
     -   (EE0I*P(16)-EijI(1,2)*P(143)-EijI(1,1)*P(148))+EijI(11,2)*P(
     -   415)+EijI(22,3)*P(416)-EijI(23,3)*P(417)-EijI(24,3)*P(418))
       F(33)=DCMPLX(FR(33),FI(33))
       P(419) = p3sq+s12-3*s34-s45
       P(420) = p2sq+3*p3sq-s23-5*s34+2*P(56)
       FR(34) = 2*(p3sq*(EijR(3,1)-EijR(8,2)-2*EijR(9,3))+4*EijR(37,4)
     -   -EijR(4,1)*P(69)-2*(EijR(11,2)-EijR(22,3)-EijR(2,3)*P(203))-
     -   (EijR(5,2)+2*EijR(8,3))*P(215)+EijR(2,1)*P(419)+EijR(2,2)*P(
     -   420))
       FI(34) = 2*(p3sq*(EijI(3,1)-EijI(8,2)-2*EijI(9,3))+4*EijI(37,4)
     -   -EijI(4,1)*P(69)-2*(EijI(11,2)-EijI(22,3)-EijI(2,3)*P(203))-
     -   (EijI(5,2)+2*EijI(8,3))*P(215)+EijI(2,1)*P(419)+EijI(2,2)*P(
     -   420))
       F(34)=DCMPLX(FR(34),FI(34))
       P(421) = p5sq-s12+s34
       P(422) = -s12+s15+s34
       P(423) = p5sq+s12-2*s15-s34
       P(424) = s15+s45-2*P(16)
       P(425) = s12-P(208)
       P(426) = p5sq-s45+2*P(425)
       FR(35) = 4*(-2*(EijR(21,3)-EijR(22,3)-EijR(37,4)+EijR(40,4)-Eij
     -   R(41,4)+EijR(43,4))+(-EijR(5,3)+EijR(6,3))*P(21)+(-EijR(10,3
     -   )+EijR(18,3)-EijR(19,3)+EijR(20,3))*P(52)+(EijR(11,3)-EijR(1
     -   2,3))*P(352)+(EijR(2,2)-EijR(5,2)+EijR(6,2)-EijR(8,2))*P(421
     -   )+EijR(2,3)*P(422)+EijR(8,3)*P(423)-EijR(9,3)*P(424)-EijR(17
     -   ,3)*P(426))
       FI(35) = 4*(-2*(EijI(21,3)-EijI(22,3)-EijI(37,4)+EijI(40,4)-Eij
     -   I(41,4)+EijI(43,4))+(-EijI(5,3)+EijI(6,3))*P(21)+(-EijI(10,3
     -   )+EijI(18,3)-EijI(19,3)+EijI(20,3))*P(52)+(EijI(11,3)-EijI(1
     -   2,3))*P(352)+(EijI(2,2)-EijI(5,2)+EijI(6,2)-EijI(8,2))*P(421
     -   )+EijI(2,3)*P(422)+EijI(8,3)*P(423)-EijI(9,3)*P(424)-EijI(17
     -   ,3)*P(426))
       F(35)=DCMPLX(FR(35),FI(35))
       P(427) = p5sq-s12+s15+s34-s45
       P(428) = s12-2*s15-s34+s45
       FR(36) = 4*(-2*(EijR(21,3)-EijR(22,3)-EijR(37,4)+EijR(40,4)-Eij
     -   R(42,4)+EijR(44,4))+(-EijR(5,3)+EijR(7,3))*P(21)+(-EijR(14,3
     -   )+EijR(15,3))*P(52)+(EijR(9,3)-EijR(17,3)+EijR(19,3)-EijR(20
     -   ,3))*P(352)+(EijR(2,2)-EijR(5,2)+EijR(7,2)-EijR(9,2))*P(421)
     -   +EijR(2,3)*P(422)+EijR(8,3)*P(423)-EijR(10,3)*P(427)-EijR(18
     -   ,3)*P(428))
       FI(36) = 4*(-2*(EijI(21,3)-EijI(22,3)-EijI(37,4)+EijI(40,4)-Eij
     -   I(42,4)+EijI(44,4))+(-EijI(5,3)+EijI(7,3))*P(21)+(-EijI(14,3
     -   )+EijI(15,3))*P(52)+(EijI(9,3)-EijI(17,3)+EijI(19,3)-EijI(20
     -   ,3))*P(352)+(EijI(2,2)-EijI(5,2)+EijI(7,2)-EijI(9,2))*P(421)
     -   +EijI(2,3)*P(422)+EijI(8,3)*P(423)-EijI(10,3)*P(427)-EijI(18
     -   ,3)*P(428))
       F(36)=DCMPLX(FR(36),FI(36))
       P(429) = p5sq+s15+s34
       P(430) = p3sq+s12-s23+s45-2*P(429)
       P(431) = p3sq+s12-s23-2*P(208)
       P(432) = p5sq+s34
       P(433) = p3sq+s12-2*P(432)
       P(434) = -s23+s45+2*P(21)
       P(435) = p3sq+s12-s15
       P(436) = -s23+s45-4*P(432)+2*P(435)
       P(437) = -4*s15-2*s23-6*s34-s45+3*P(124)
       P(438) = -s23-6*s34+3*P(124)-2*P(282)
       FR(37) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(2,1)+Dij2345R(
     -   2,2)-2*Dij2345R(4,2)+4*(EijR(22,3)-EijR(23,3)+EijR(37,4)+Eij
     -   R(38,4)-2*EijR(43,4))+(EijR(2,1)-EijR(3,1))*P(17)-(EijR(9,2)
     -   -EijR(10,2)+EijR(10,3)+EijR(13,3)-2*EijR(20,3))*P(52)+EijR(3
     -   ,3)*P(104)-EijR(2,2)*P(430)-EijR(2,3)*P(431)-EijR(3,2)*P(433
     -   )+(EijR(5,2)-EijR(6,2)+EijR(8,3)+EijR(11,3)-2*EijR(17,3))*P(
     -   434)+EijR(8,2)*P(436)+EijR(9,3)*P(437)-EijR(12,3)*P(438))
       FI(37) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(2,1)+Dij2345I(
     -   2,2)-2*Dij2345I(4,2)+4*(EijI(22,3)-EijI(23,3)+EijI(37,4)+Eij
     -   I(38,4)-2*EijI(43,4))+(EijI(2,1)-EijI(3,1))*P(17)-(EijI(9,2)
     -   -EijI(10,2)+EijI(10,3)+EijI(13,3)-2*EijI(20,3))*P(52)+EijI(3
     -   ,3)*P(104)-EijI(2,2)*P(430)-EijI(2,3)*P(431)-EijI(3,2)*P(433
     -   )+(EijI(5,2)-EijI(6,2)+EijI(8,3)+EijI(11,3)-2*EijI(17,3))*P(
     -   434)+EijI(8,2)*P(436)+EijI(9,3)*P(437)-EijI(12,3)*P(438))
       F(37)=DCMPLX(FR(37),FI(37))
       P(439) = s15+s34-s45
       P(440) = p3sq-3*p5sq+s12-s23-2*P(439)
       P(441) = -s23-4*s34-s45+2*P(435)
       P(442) = p3sq-p5sq+s12-s23+s45-2*P(208)
       P(443) = -p5sq-s23-4*s34+2*P(435)
       FR(38) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(3,1)-Dij2345R(
     -   4,2)-Dij2345R(5,2)+Dij2345R(6,2)+6*(EijR(11,2)+EijR(22,3)-Ei
     -   jR(23,3))+4*(EijR(37,4)-EijR(43,4)-EijR(44,4)+EijR(45,4))+(E
     -   ijR(2,1)-EijR(4,1))*P(17)+(EijR(4,2)+EijR(15,3)-EijR(16,3))*
     -   P(52)+(-EijR(12,3)+EijR(13,3))*P(104)-EijR(2,2)*P(430)-EijR(
     -   2,3)*P(431)+(EijR(8,2)-EijR(10,2))*P(433)+(EijR(5,2)-EijR(7,
     -   2)+EijR(8,3)-EijR(17,3)-EijR(18,3)+EijR(19,3))*P(434)+EijR(9
     -   ,2)*P(440)+EijR(9,3)*P(441)+EijR(10,3)*P(442)-EijR(20,3)*P(4
     -   43))
       FI(38) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(3,1)-Dij2345I(
     -   4,2)-Dij2345I(5,2)+Dij2345I(6,2)+6*(EijI(11,2)+EijI(22,3)-Ei
     -   jI(23,3))+4*(EijI(37,4)-EijI(43,4)-EijI(44,4)+EijI(45,4))+(E
     -   ijI(2,1)-EijI(4,1))*P(17)+(EijI(4,2)+EijI(15,3)-EijI(16,3))*
     -   P(52)+(-EijI(12,3)+EijI(13,3))*P(104)-EijI(2,2)*P(430)-EijI(
     -   2,3)*P(431)+(EijI(8,2)-EijI(10,2))*P(433)+(EijI(5,2)-EijI(7,
     -   2)+EijI(8,3)-EijI(17,3)-EijI(18,3)+EijI(19,3))*P(434)+EijI(9
     -   ,2)*P(440)+EijI(9,3)*P(441)+EijI(10,3)*P(442)-EijI(20,3)*P(4
     -   43))
       F(38)=DCMPLX(FR(38),FI(38))
       P(444) = p3sq-p5sq+s12-s15
       P(445) = s23+4*s34-2*P(444)
       P(446) = 3*p5sq-2*s45
       FR(39) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(
     -   4,2)-Dij2345R(5,2)+Dij2345R(6,2)+10*EijR(22,3)-6*EijR(23,3)-
     -   4*(EijR(24,3)-EijR(37,4)+EijR(43,4)+EijR(44,4)-EijR(45,4))+(
     -   EijR(2,1)-EijR(3,1))*P(17)+(EijR(15,3)-EijR(16,3))*P(52)-(Ei
     -   jR(3,2)+EijR(12,3)-EijR(13,3))*P(104)-EijR(2,2)*P(430)-EijR(
     -   2,3)*P(431)+(EijR(5,2)-EijR(6,2)+EijR(8,3)-EijR(17,3)-EijR(1
     -   8,3)+EijR(19,3))*P(434)+EijR(9,3)*P(441)+EijR(10,3)*P(442)-E
     -   ijR(20,3)*P(443)-EijR(8,2)*P(445)+(-EijR(9,2)+EijR(10,2))*P(
     -   446))
       FI(39) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(
     -   4,2)-Dij2345I(5,2)+Dij2345I(6,2)+10*EijI(22,3)-6*EijI(23,3)-
     -   4*(EijI(24,3)-EijI(37,4)+EijI(43,4)+EijI(44,4)-EijI(45,4))+(
     -   EijI(2,1)-EijI(3,1))*P(17)+(EijI(15,3)-EijI(16,3))*P(52)-(Ei
     -   jI(3,2)+EijI(12,3)-EijI(13,3))*P(104)-EijI(2,2)*P(430)-EijI(
     -   2,3)*P(431)+(EijI(5,2)-EijI(6,2)+EijI(8,3)-EijI(17,3)-EijI(1
     -   8,3)+EijI(19,3))*P(434)+EijI(9,3)*P(441)+EijI(10,3)*P(442)-E
     -   ijI(20,3)*P(443)-EijI(8,2)*P(445)+(-EijI(9,2)+EijI(10,2))*P(
     -   446))
       F(39)=DCMPLX(FR(39),FI(39))
       P(447) = p3sq-5*p5sq+s12-s23+3*s45-2*P(208)
       P(448) = s23+2*P(208)
       P(449) = p3sq+s12-P(448)
       P(450) = -p5sq+s45+2*P(449)
       P(451) = p5sq+s15+s34-s45
       P(452) = p3sq+s12-s23-2*P(451)
       FR(40) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(3,1)+Dij2345R(
     -   3,2)-2*Dij2345R(5,2)+6*EijR(11,2)+12*(EijR(22,3)-EijR(24,3))
     -   +4*(EijR(37,4)+EijR(39,4))-8*EijR(44,4)+(EijR(2,1)-EijR(4,1)
     -   )*P(17)-EijR(4,3)*P(52)+(EijR(8,2)+EijR(9,3)-EijR(10,2)+EijR
     -   (16,3)-2*EijR(20,3))*P(104)-EijR(2,2)*P(430)-EijR(2,3)*P(431
     -   )+(EijR(5,2)-EijR(7,2)+EijR(8,3)+EijR(14,3)-2*EijR(18,3))*P(
     -   434)+EijR(4,2)*P(446)+EijR(9,2)*P(447)+EijR(10,3)*P(450)-Eij
     -   R(15,3)*P(452))
       FI(40) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(3,1)+Dij2345I(
     -   3,2)-2*Dij2345I(5,2)+6*EijI(11,2)+12*(EijI(22,3)-EijI(24,3))
     -   +4*(EijI(37,4)+EijI(39,4))-8*EijI(44,4)+(EijI(2,1)-EijI(4,1)
     -   )*P(17)-EijI(4,3)*P(52)+(EijI(8,2)+EijI(9,3)-EijI(10,2)+EijI
     -   (16,3)-2*EijI(20,3))*P(104)-EijI(2,2)*P(430)-EijI(2,3)*P(431
     -   )+(EijI(5,2)-EijI(7,2)+EijI(8,3)+EijI(14,3)-2*EijI(18,3))*P(
     -   434)+EijI(4,2)*P(446)+EijI(9,2)*P(447)+EijI(10,3)*P(450)-Eij
     -   I(15,3)*P(452))
       F(40)=DCMPLX(FR(40),FI(40))
       P(453) = s45-2*P(421)
       P(454) = p2sq-p5sq-s15
       P(455) = p3sq+3*s12-s23-4*s34+s45+2*P(454)
       P(456) = -s23+s45+2*P(158)
       P(457) = p2sq+p3sq-p5sq-s15-s45
       P(458) = -s23+6*P(16)+2*P(457)
       FR(41) = -2*(Dij2345R(1,1)-Dij2345R(1,2)-Dij2345R(2,1)+Dij2345R
     -   (4,2)-14*EijR(22,3)+10*EijR(23,3)+4*(EijR(11,2)-EijR(37,4)+E
     -   ijR(43,4))+(3*(EijR(9,2)-EijR(10,2))+EijR(10,3)-EijR(20,3))*
     -   P(52)+EijR(12,3)*P(104)+EijR(3,2)*P(122)+EijR(2,3)*P(431)+(-
     -   EijR(8,3)+EijR(17,3))*P(434)-EijR(9,3)*P(441)+(EijR(2,1)-Eij
     -   R(3,1))*P(453)+EijR(2,2)*P(455)+(-EijR(5,2)+EijR(6,2))*P(456
     -   )-EijR(8,2)*P(458))
       FI(41) = -2*(Dij2345I(1,1)-Dij2345I(1,2)-Dij2345I(2,1)+Dij2345I
     -   (4,2)-14*EijI(22,3)+10*EijI(23,3)+4*(EijI(11,2)-EijI(37,4)+E
     -   ijI(43,4))+(3*(EijI(9,2)-EijI(10,2))+EijI(10,3)-EijI(20,3))*
     -   P(52)+EijI(12,3)*P(104)+EijI(3,2)*P(122)+EijI(2,3)*P(431)+(-
     -   EijI(8,3)+EijI(17,3))*P(434)-EijI(9,3)*P(441)+(EijI(2,1)-Eij
     -   I(3,1))*P(453)+EijI(2,2)*P(455)+(-EijI(5,2)+EijI(6,2))*P(456
     -   )-EijI(8,2)*P(458))
       F(41)=DCMPLX(FR(41),FI(41))
       P(459) = p3sq-5*p5sq+3*s12-s23+2*P(152)-4*P(389)
       FR(42) = -2*(Dij2345R(1,1)-Dij2345R(1,2)-Dij2345R(3,1)+Dij2345R
     -   (5,2)-2*(EijR(11,2)+8*EijR(22,3)-5*EijR(24,3)+2*(EijR(37,4)-
     -   EijR(44,4)))-(3*EijR(4,2)+EijR(15,3))*P(52)+(-EijR(9,3)+EijR
     -   (20,3))*P(104)+(-EijR(8,2)+EijR(10,2))*P(122)+EijR(2,3)*P(43
     -   1)+(-EijR(8,3)+EijR(18,3))*P(434)-EijR(10,3)*P(442)+(EijR(2,
     -   1)-EijR(4,1))*P(453)+EijR(2,2)*P(455)+(-EijR(5,2)+EijR(7,2))
     -   *P(456)-EijR(9,2)*P(459))
       FI(42) = -2*(Dij2345I(1,1)-Dij2345I(1,2)-Dij2345I(3,1)+Dij2345I
     -   (5,2)-2*(EijI(11,2)+8*EijI(22,3)-5*EijI(24,3)+2*(EijI(37,4)-
     -   EijI(44,4)))-(3*EijI(4,2)+EijI(15,3))*P(52)+(-EijI(9,3)+EijI
     -   (20,3))*P(104)+(-EijI(8,2)+EijI(10,2))*P(122)+EijI(2,3)*P(43
     -   1)+(-EijI(8,3)+EijI(18,3))*P(434)-EijI(10,3)*P(442)+(EijI(2,
     -   1)-EijI(4,1))*P(453)+EijI(2,2)*P(455)+(-EijI(5,2)+EijI(7,2))
     -   *P(456)-EijI(9,2)*P(459))
       F(42)=DCMPLX(FR(42),FI(42))
       P(460) = p5sq-2*s12+s34
       FR(43) = -4*(D02345R+Dij2345R(1,1)+4*(EijR(21,3)-EijR(22,3))-2*
     -   (EijR(37,4)-EijR(40,4))+EijR(5,3)*P(21)+(-EijR(7,2)+EijR(9,2
     -   ))*P(43)+(EijR(10,3)-EijR(18,3))*P(52)+(EijR(6,2)-EijR(8,2))
     -   *P(264)+(-EijR(9,3)+EijR(17,3))*P(352)-EijR(2,3)*P(422)-EijR
     -   (8,3)*P(423)+(-EijR(2,2)+EijR(5,2))*P(460))
       FI(43) = -4*(D02345I+Dij2345I(1,1)+4*(EijI(21,3)-EijI(22,3))-2*
     -   (EijI(37,4)-EijI(40,4))+EijI(5,3)*P(21)+(-EijI(7,2)+EijI(9,2
     -   ))*P(43)+(EijI(10,3)-EijI(18,3))*P(52)+(EijI(6,2)-EijI(8,2))
     -   *P(264)+(-EijI(9,3)+EijI(17,3))*P(352)-EijI(2,3)*P(422)-EijI
     -   (8,3)*P(423)+(-EijI(2,2)+EijI(5,2))*P(460))
       F(43)=DCMPLX(FR(43),FI(43))
       P(461) = p2sq-2*p5sq-s12+s45
       P(462) = p2sq+p3sq-s23-2*P(208)
       P(463) = p2sq-s12-2*s15+s45
       P(464) = -p2sq+s12+s23+2*P(21)
       P(465) = 2*s15+s23+s34
       P(466) = p5sq-P(465)
       P(467) = p2sq+p3sq+s45+2*P(466)
       P(468) = p2sq+p3sq-s15-s23-2*s34
       P(469) = p2sq+p3sq-s23-2*s34
       P(470) = p2sq+p3sq+p5sq
       P(471) = -3*s23+s45-4*P(208)+2*P(470)
       FR(44) = 2*(D02345R+3*Dij2345R(1,1)-6*(EijR(11,2)+EijR(22,3)-Ei
     -   jR(23,3))+4*(EijR(37,4)-EijR(40,4)+EijR(41,4)-EijR(43,4))+(E
     -   ijR(7,2)-EijR(9,2))*P(191)+(-EijR(10,3)+EijR(18,3)-EijR(19,3
     -   )+EijR(20,3))*P(200)+(-EijR(5,3)+EijR(6,3))*P(434)-EijR(2,2)
     -   *P(461)-EijR(2,3)*P(462)+EijR(5,2)*P(463)+(EijR(6,2)-EijR(8,
     -   2))*P(464)+EijR(8,3)*P(467)+2*(Dij2345R(1,2)-Dij2345R(2,1)-D
     -   ij2345R(4,2)-(EijR(1,1)+EijR(1,2)-EijR(2,1))*P(21)+EijR(9,3)
     -   *P(468))+(EijR(11,3)-EijR(12,3))*P(469)-EijR(17,3)*P(471))
       FI(44) = 2*(D02345I+3*Dij2345I(1,1)-6*(EijI(11,2)+EijI(22,3)-Ei
     -   jI(23,3))+4*(EijI(37,4)-EijI(40,4)+EijI(41,4)-EijI(43,4))+(E
     -   ijI(7,2)-EijI(9,2))*P(191)+(-EijI(10,3)+EijI(18,3)-EijI(19,3
     -   )+EijI(20,3))*P(200)+(-EijI(5,3)+EijI(6,3))*P(434)-EijI(2,2)
     -   *P(461)-EijI(2,3)*P(462)+EijI(5,2)*P(463)+(EijI(6,2)-EijI(8,
     -   2))*P(464)+EijI(8,3)*P(467)+2*(Dij2345I(1,2)-Dij2345I(2,1)-D
     -   ij2345I(4,2)-(EijI(1,1)+EijI(1,2)-EijI(2,1))*P(21)+EijI(9,3)
     -   *P(468))+(EijI(11,3)-EijI(12,3))*P(469)-EijI(17,3)*P(471))
       F(44)=DCMPLX(FR(44),FI(44))
       P(472) = p2sq+s12-s23-2*s34-s45
       P(473) = s15+s23+2*P(98)
       P(474) = 3*p5sq+2*s12-P(473)
       P(475) = s12-P(98)
       P(476) = -3*p5sq+s15+s23-2*P(475)
       P(477) = p2sq+p3sq-3*s15-2*s34
       P(478) = p2sq+p3sq-5*s15-s23+s45+2*P(10)
       FR(45) = 2*(D02345R+3*Dij2345R(1,1)-6*(EijR(11,2)+EijR(21,3)-Ei
     -   jR(24,3))+4*(EijR(37,4)-EijR(40,4)+EijR(42,4)-EijR(44,4))+2*
     -   (Dij2345R(1,2)-Dij2345R(3,1)-Dij2345R(5,2)-(EijR(1,1)+EijR(1
     -   ,2)-EijR(2,1))*P(21))+(-EijR(14,3)+EijR(15,3))*P(200)+(-EijR
     -   (5,3)+EijR(7,3))*P(434)-EijR(2,2)*P(461)-EijR(2,3)*P(462)+Ei
     -   jR(5,2)*P(463)+EijR(8,3)*P(467)+(EijR(9,3)-EijR(17,3)+EijR(1
     -   9,3)-EijR(20,3))*P(469)+(-EijR(6,2)+EijR(8,2))*P(472)+EijR(7
     -   ,2)*P(474)+EijR(9,2)*P(476)+EijR(10,3)*P(477)-EijR(18,3)*P(4
     -   78))
       FI(45) = 2*(D02345I+3*Dij2345I(1,1)-6*(EijI(11,2)+EijI(21,3)-Ei
     -   jI(24,3))+4*(EijI(37,4)-EijI(40,4)+EijI(42,4)-EijI(44,4))+2*
     -   (Dij2345I(1,2)-Dij2345I(3,1)-Dij2345I(5,2)-(EijI(1,1)+EijI(1
     -   ,2)-EijI(2,1))*P(21))+(-EijI(14,3)+EijI(15,3))*P(200)+(-EijI
     -   (5,3)+EijI(7,3))*P(434)-EijI(2,2)*P(461)-EijI(2,3)*P(462)+Ei
     -   jI(5,2)*P(463)+EijI(8,3)*P(467)+(EijI(9,3)-EijI(17,3)+EijI(1
     -   9,3)-EijI(20,3))*P(469)+(-EijI(6,2)+EijI(8,2))*P(472)+EijI(7
     -   ,2)*P(474)+EijI(9,2)*P(476)+EijI(10,3)*P(477)-EijI(18,3)*P(4
     -   78))
       F(45)=DCMPLX(FR(45),FI(45))
       P(479) = p2sq+s12-2*P(208)
       P(480) = p2sq+4*p5sq+s12-6*s15-2*s34
       P(481) = p2sq+p5sq+s12-3*s15-2*s34
       FR(46) = 2*(D02345R+Dij2345R(1,2)+4*(EijR(36,4)+EijR(37,4)+2*(E
     -   ijR(21,3)-EijR(22,3)-EijR(40,4)))-(EijR(7,3)+EijR(10,3)-2*Ei
     -   jR(18,3))*P(191)+(EijR(6,3)+EijR(9,3)-2*EijR(17,3))*P(472)-E
     -   ijR(2,3)*P(479)-EijR(5,3)*P(480)+2*(Dij2345R(1,1)+(EijR(1,2)
     -   +EijR(1,3)+EijR(2,2)-2*EijR(5,2))*P(21)+EijR(8,3)*P(481)))
       FI(46) = 2*(D02345I+Dij2345I(1,2)+4*(EijI(36,4)+EijI(37,4)+2*(E
     -   ijI(21,3)-EijI(22,3)-EijI(40,4)))-(EijI(7,3)+EijI(10,3)-2*Ei
     -   jI(18,3))*P(191)+(EijI(6,3)+EijI(9,3)-2*EijI(17,3))*P(472)-E
     -   ijI(2,3)*P(479)-EijI(5,3)*P(480)+2*(Dij2345I(1,1)+(EijI(1,2)
     -   +EijI(1,3)+EijI(2,2)-2*EijI(5,2))*P(21)+EijI(8,3)*P(481)))
       F(46)=DCMPLX(FR(46),FI(46))
       P(482) = p2sq-s12+s45-2*P(432)
       P(483) = p2sq-s12+s45-2*P(208)
       P(484) = p3sq+s34
       P(485) = -p2sq-s12+s23+s45+2*P(484)
       P(486) = -p5sq-s15+s23+s45+2*P(69)
       P(487) = p2sq+s12-s23-s45-2*P(484)
       FR(47) = 2*(D02345R+3*Dij2345R(1,1)+2*(Dij2345R(1,2)-3*EijR(11,
     -   2)-5*EijR(21,3)+2*(EijR(22,3)+EijR(37,4)-EijR(40,4))-(EijR(1
     -   ,1)+EijR(1,2)-EijR(2,1))*P(21))+(-EijR(10,3)+EijR(18,3))*P(2
     -   00)-EijR(5,3)*P(434)-EijR(2,3)*P(462)+EijR(8,3)*P(467)+(EijR
     -   (9,3)-EijR(17,3))*P(469)-EijR(2,2)*P(482)+EijR(5,2)*P(483)+E
     -   ijR(6,2)*P(485)+(-EijR(7,2)+EijR(9,2))*P(486)+EijR(8,2)*P(48
     -   7))
       FI(47) = 2*(D02345I+3*Dij2345I(1,1)+2*(Dij2345I(1,2)-3*EijI(11,
     -   2)-5*EijI(21,3)+2*(EijI(22,3)+EijI(37,4)-EijI(40,4))-(EijI(1
     -   ,1)+EijI(1,2)-EijI(2,1))*P(21))+(-EijI(10,3)+EijI(18,3))*P(2
     -   00)-EijI(5,3)*P(434)-EijI(2,3)*P(462)+EijI(8,3)*P(467)+(EijI
     -   (9,3)-EijI(17,3))*P(469)-EijI(2,2)*P(482)+EijI(5,2)*P(483)+E
     -   ijI(6,2)*P(485)+(-EijI(7,2)+EijI(9,2))*P(486)+EijI(8,2)*P(48
     -   7))
       F(47)=DCMPLX(FR(47),FI(47))
       P(488) = p2sq+s23-3*s34+2*P(56)
       P(489) = p2sq+s15
       P(490) = -p3sq+s23+s34+2*P(489)
       P(491) = p2sq-s23+s34
       P(492) = -p3sq-s15+s23+s34
       P(493) = 2*p5sq+s12-s34-s45
       P(494) = p2sq*P(274)+s12*P(493)
       P(495) = -3*s34+2*P(163)
       P(496) = p3sq-2*s15-s23-s34
       P(497) = p2sq**2+p2sq*P(495)+s12*P(496)
       P(498) = p5sq+s12-s15-s34
       P(499) = p3sq+2*P(498)
       P(500) = p3sq*p5sq+p5sq*s12-s12*s15
       P(501) = -(p2sq*s12)+s12**2+s12*s23+p2sq*s34-s12*s34-s23*s34+s4
     -   5**2-s45*P(499)+2*P(500)
       P(502) = p2sq-s23-s34
       P(503) = p3sq**2+p3sq*P(325)+P(56)*P(502)
       P(504) = p3sq+3*p5sq+s12-2*P(208)
       P(505) = p3sq*p5sq+p5sq**2-p5sq*s15
       P(506) = p5sq*s12-s12*s15+s12*s23-3*p5sq*s34+s15*s34-s23*s34+s4
     -   5**2-s45*P(504)+2*P(505)
       P(507) = p5sq-2*s15-s34
       P(508) = -s12-s23+s45+2*P(432)
       P(509) = p2sq**2+2*p2sq*P(507)+s12*P(508)
       P(510) = -p5sq+s12+s15
       P(511) = -s23+s45
       P(512) = s23+2*s34+s45
       P(513) = p3sq*p5sq+p5sq*s12+s12*s34-s15*P(124)+s45*P(510)
       P(514) = p2sq**2-s12**2+P(98)*P(193)+p3sq*P(511)-p2sq*P(512)+2*
     -   P(513)
       P(515) = p3sq-p5sq+s12+2*s15+s23-s34
       P(516) = -p5sq-s15+s23+s45
       P(517) = p5sq*s12+s12*s15+p3sq*s23+p5sq*s23-s12*s23-s23*s34+s45
     -   **2-2*P(293)-s45*P(515)+p2sq*P(516)
       P(518) = s15+s23+s34
       P(519) = -(s23*s34)+s12*P(518)
       P(520) = -s23-2*s45+3*P(16)
       P(521) = -2*s15-s23-s45+3*P(16)
       P(522) = p2sq**2+p3sq**2+s45*P(465)-2*P(519)+p2sq*P(520)+p3sq*P
     -   (521)
       P(523) = s15-s23+2*P(52)
       P(524) = p5sq+s12-2*s15-s23-3*s34-s45
       P(525) = s23*s34+s34**2
       P(526) = p3sq**2+s12*s15-s12*s23-s12*s34+s15*s34-p5sq*P(465)+s4
     -   5*P(465)+p2sq*P(523)+p3sq*P(524)+2*P(525)
       P(527) = -s15+s23+s34
       P(528) = s15-2*s23
       P(529) = 2*p5sq+s12-3*P(98)
       P(530) = -p3sq**2-2*s34*P(98)+p5sq*P(237)+p2sq*P(273)+s12*P(527
     -   )+s45*P(528)-p3sq*P(529)
       P(531) = p2sq+p3sq+s12
       P(532) = s15+s23+s45
       P(533) = -9*s34+3*P(531)-2*P(532)
       FR(48) = 4*Dij2345R(7,2)+8*EijR(46,4)+EijR(4,2)*P(52)*P(204)-D0
     -   2345R*P(488)-Dij2345R(1,1)*P(490)+Dij2345R(2,1)*P(491)+Dij23
     -   45R(3,1)*P(492)-EijR(2,1)*P(494)+EijR(2,2)*P(497)+EijR(3,1)*
     -   P(501)+EijR(3,2)*P(503)-EijR(4,1)*P(506)-EijR(5,2)*P(509)+Ei
     -   jR(6,2)*P(514)+EijR(7,2)*P(517)-EijR(8,2)*P(522)+EijR(9,2)*P
     -   (526)+EijR(10,2)*P(530)-2*(C0345R-P(21)*(-(EE0R*P(16))+EijR(
     -   1,2)*P(143)+EijR(1,1)*P(148))+EijR(24,3)*P(200)-EijR(21,3)*P
     -   (434)+EijR(22,3)*P(462)-EijR(23,3)*P(469)-EijR(11,2)*P(533))
       FI(48) = 4*Dij2345I(7,2)+8*EijI(46,4)+EijI(4,2)*P(52)*P(204)-D0
     -   2345I*P(488)-Dij2345I(1,1)*P(490)+Dij2345I(2,1)*P(491)+Dij23
     -   45I(3,1)*P(492)-EijI(2,1)*P(494)+EijI(2,2)*P(497)+EijI(3,1)*
     -   P(501)+EijI(3,2)*P(503)-EijI(4,1)*P(506)-EijI(5,2)*P(509)+Ei
     -   jI(6,2)*P(514)+EijI(7,2)*P(517)-EijI(8,2)*P(522)+EijI(9,2)*P
     -   (526)+EijI(10,2)*P(530)-2*(C0345I-P(21)*(-(EE0I*P(16))+EijI(
     -   1,2)*P(143)+EijI(1,1)*P(148))+EijI(24,3)*P(200)-EijI(21,3)*P
     -   (434)+EijI(22,3)*P(462)-EijI(23,3)*P(469)-EijI(11,2)*P(533))
       F(48)=DCMPLX(FR(48),FI(48))
       P(534) = p2sq+p5sq+s15+s34
       P(535) = p3sq+s12-s23+s45-2*P(534)
       P(536) = -s23+s45-2*P(394)
       P(537) = p2sq+p5sq
       P(538) = p3sq+s12-4*s34-2*P(537)
       FR(49) = -2*(-Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(4,2)-EijR(2,
     -   1)*P(17)+(2*EijR(4,1)+EijR(9,2)+EijR(10,3)-EijR(20,3))*P(52)
     -   +EijR(12,3)*P(104)+2*(D02345R+s45*EijR(3,1)-s34*EijR(3,2)-7*
     -   EijR(11,2)-5*(EijR(22,3)-EijR(23,3))-2*(EijR(37,4)-EijR(43,4
     -   ))-EijR(6,2)*P(143))+EijR(2,3)*P(431)+(-EijR(8,3)+EijR(17,3)
     -   )*P(434)-EijR(9,3)*P(441)+EijR(2,2)*P(535)-EijR(5,2)*P(536)-
     -   EijR(8,2)*P(538))
       FI(49) = -2*(-Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(4,2)-EijI(2,
     -   1)*P(17)+(2*EijI(4,1)+EijI(9,2)+EijI(10,3)-EijI(20,3))*P(52)
     -   +EijI(12,3)*P(104)+2*(D02345I+s45*EijI(3,1)-s34*EijI(3,2)-7*
     -   EijI(11,2)-5*(EijI(22,3)-EijI(23,3))-2*(EijI(37,4)-EijI(43,4
     -   ))-EijI(6,2)*P(143))+EijI(2,3)*P(431)+(-EijI(8,3)+EijI(17,3)
     -   )*P(434)-EijI(9,3)*P(441)+EijI(2,2)*P(535)-EijI(5,2)*P(536)-
     -   EijI(8,2)*P(538))
       F(49)=DCMPLX(FR(49),FI(49))
       P(539) = -p3sq-s12+s23-s45+2*P(534)
       P(540) = p3sq+s12-4*s34-s45
       P(541) = p2sq-s45
       P(542) = 3*p5sq+2*P(541)
       FR(50) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(5,2)+EijR(2,1)
     -   *P(17)+EijR(15,3)*P(52)+(EijR(9,3)-EijR(20,3))*P(104)-2*(D02
     -   345R-s34*EijR(10,2)-7*EijR(11,2)-8*EijR(22,3)+5*EijR(24,3)-2
     -   *(EijR(37,4)-EijR(44,4))+EijR(4,1)*P(89)-EijR(7,2)*P(143)-Ei
     -   jR(3,1)*P(352))-EijR(2,3)*P(431)+(EijR(8,3)-EijR(18,3))*P(43
     -   4)+EijR(10,3)*P(442)+EijR(5,2)*P(536)+EijR(2,2)*P(539)+EijR(
     -   8,2)*P(540)-EijR(9,2)*P(542))
       FI(50) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(5,2)+EijI(2,1)
     -   *P(17)+EijI(15,3)*P(52)+(EijI(9,3)-EijI(20,3))*P(104)-2*(D02
     -   345I-s34*EijI(10,2)-7*EijI(11,2)-8*EijI(22,3)+5*EijI(24,3)-2
     -   *(EijI(37,4)-EijI(44,4))+EijI(4,1)*P(89)-EijI(7,2)*P(143)-Ei
     -   jI(3,1)*P(352))-EijI(2,3)*P(431)+(EijI(8,3)-EijI(18,3))*P(43
     -   4)+EijI(10,3)*P(442)+EijI(5,2)*P(536)+EijI(2,2)*P(539)+EijI(
     -   8,2)*P(540)-EijI(9,2)*P(542))
       F(50)=DCMPLX(FR(50),FI(50))
       P(543) = -s45+2*P(421)
       P(544) = p5sq+s15
       P(545) = p3sq+3*s12-s23-4*s34+s45-2*P(544)
       P(546) = p3sq-4*s34+3*P(56)
       FR(51) = 2*(-Dij2345R(1,1)+Dij2345R(1,2)+10*EijR(11,2)+20*EijR(
     -   22,3)+4*EijR(37,4)-(2*EijR(4,1)+3*EijR(9,2)+EijR(10,3))*P(52
     -   )+EijR(9,3)*P(104)-2*(D02345R-EijR(3,1)*P(352))-EijR(2,3)*P(
     -   431)+(EijR(5,2)+EijR(8,3))*P(434)+EijR(2,1)*P(543)-EijR(2,2)
     -   *P(545)+EijR(8,2)*P(546))
       FI(51) = 2*(-Dij2345I(1,1)+Dij2345I(1,2)+10*EijI(11,2)+20*EijI(
     -   22,3)+4*EijI(37,4)-(2*EijI(4,1)+3*EijI(9,2)+EijI(10,3))*P(52
     -   )+EijI(9,3)*P(104)-2*(D02345I-EijI(3,1)*P(352))-EijI(2,3)*P(
     -   431)+(EijI(5,2)+EijI(8,3))*P(434)+EijI(2,1)*P(543)-EijI(2,2)
     -   *P(545)+EijI(8,2)*P(546))
       F(51)=DCMPLX(FR(51),FI(51))
       FR(52) = 4*(p5sq*(EijR(2,1)-EijR(3,1)+EijR(5,2)-EijR(6,2))+s15*
     -   (EijR(2,2)-EijR(5,2)+EijR(6,2)-EijR(8,2))+2*(EijR(37,4)+EijR
     -   (38,4)-2*EijR(43,4)))
       FI(52) = 4*(p5sq*(EijI(2,1)-EijI(3,1)+EijI(5,2)-EijI(6,2))+s15*
     -   (EijI(2,2)-EijI(5,2)+EijI(6,2)-EijI(8,2))+2*(EijI(37,4)+EijI
     -   (38,4)-2*EijI(43,4)))
       F(52)=DCMPLX(FR(52),FI(52))
       FR(53) = 4*(p5sq*(EijR(2,1)-EijR(4,1)+EijR(5,2)-EijR(7,2))+s15*
     -   (EijR(2,2)-EijR(5,2)+EijR(7,2)-EijR(9,2))-2*(EijR(23,3)-EijR
     -   (24,3)-EijR(37,4)+EijR(43,4)+EijR(44,4)-EijR(45,4)))
       FI(53) = 4*(p5sq*(EijI(2,1)-EijI(4,1)+EijI(5,2)-EijI(7,2))+s15*
     -   (EijI(2,2)-EijI(5,2)+EijI(7,2)-EijI(9,2))-2*(EijI(23,3)-EijI
     -   (24,3)-EijI(37,4)+EijI(43,4)+EijI(44,4)-EijI(45,4)))
       F(53)=DCMPLX(FR(53),FI(53))
       FR(54) = -2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(2,1)+Dij2345R
     -   (2,2)-2*Dij2345R(4,2)+p5sq*(EijR(9,2)-EijR(10,2)+EijR(10,3)+
     -   EijR(13,3)+2*(EijR(2,1)-EijR(3,1)+EijR(5,2)-EijR(6,2)-EijR(2
     -   0,3)))-s15*(EijR(9,2)-EijR(10,2)+EijR(10,3)+EijR(13,3)+2*(Ei
     -   jR(2,1)-EijR(3,1)+EijR(3,2)+EijR(5,2)-EijR(6,2)-EijR(8,2)-Ei
     -   jR(20,3)))-4*(EijR(22,3)-EijR(23,3)+EijR(37,4)+EijR(38,4)-2*
     -   EijR(43,4))-(EijR(2,2)+EijR(2,3)+EijR(3,2)-EijR(3,3)-2*EijR(
     -   8,2)-3*(EijR(9,3)-EijR(12,3)))*P(143)+(EijR(3,2)-EijR(3,3)-E
     -   ijR(8,2)+EijR(9,2)-EijR(9,3)-EijR(10,2)+EijR(10,3)+EijR(13,3
     -   )+2*(EijR(12,3)-EijR(20,3)))*P(193))
       FI(54) = -2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(2,1)+Dij2345I
     -   (2,2)-2*Dij2345I(4,2)+p5sq*(EijI(9,2)-EijI(10,2)+EijI(10,3)+
     -   EijI(13,3)+2*(EijI(2,1)-EijI(3,1)+EijI(5,2)-EijI(6,2)-EijI(2
     -   0,3)))-s15*(EijI(9,2)-EijI(10,2)+EijI(10,3)+EijI(13,3)+2*(Ei
     -   jI(2,1)-EijI(3,1)+EijI(3,2)+EijI(5,2)-EijI(6,2)-EijI(8,2)-Ei
     -   jI(20,3)))-4*(EijI(22,3)-EijI(23,3)+EijI(37,4)+EijI(38,4)-2*
     -   EijI(43,4))-(EijI(2,2)+EijI(2,3)+EijI(3,2)-EijI(3,3)-2*EijI(
     -   8,2)-3*(EijI(9,3)-EijI(12,3)))*P(143)+(EijI(3,2)-EijI(3,3)-E
     -   ijI(8,2)+EijI(9,2)-EijI(9,3)-EijI(10,2)+EijI(10,3)+EijI(13,3
     -   )+2*(EijI(12,3)-EijI(20,3)))*P(193))
       F(54)=DCMPLX(FR(54),FI(54))
       P(547) = -s23+s45+2*P(168)
       P(548) = s23-s45-2*P(143)
       P(549) = p5sq+s15+s23-s45
       P(550) = p2sq+p5sq-s12-s15+s23-s45
       FR(55) = 2*(-Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(2,1)+Dij2345R
     -   (4,2)+Dij2345R(5,2)-Dij2345R(6,2)-4*(-2*EijR(22,3)+EijR(23,3
     -   )+EijR(24,3)-EijR(37,4)+EijR(43,4)+EijR(44,4)-EijR(45,4))-2*
     -   (s15*EijR(3,1)-p5sq*EijR(4,1)+(EijR(2,1)+EijR(5,2)-EijR(7,2)
     -   )*P(21))+(EijR(2,2)+EijR(2,3))*P(143)+EijR(20,3)*P(159)-(Eij
     -   R(9,2)-EijR(15,3)+EijR(16,3))*P(186)+(EijR(3,2)+EijR(12,3)-E
     -   ijR(13,3))*P(215)-EijR(8,2)*P(547)+EijR(9,3)*P(548)+EijR(10,
     -   2)*P(549)-EijR(10,3)*P(550))
       FI(55) = 2*(-Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(2,1)+Dij2345I
     -   (4,2)+Dij2345I(5,2)-Dij2345I(6,2)-4*(-2*EijI(22,3)+EijI(23,3
     -   )+EijI(24,3)-EijI(37,4)+EijI(43,4)+EijI(44,4)-EijI(45,4))-2*
     -   (s15*EijI(3,1)-p5sq*EijI(4,1)+(EijI(2,1)+EijI(5,2)-EijI(7,2)
     -   )*P(21))+(EijI(2,2)+EijI(2,3))*P(143)+EijI(20,3)*P(159)-(Eij
     -   I(9,2)-EijI(15,3)+EijI(16,3))*P(186)+(EijI(3,2)+EijI(12,3)-E
     -   ijI(13,3))*P(215)-EijI(8,2)*P(547)+EijI(9,3)*P(548)+EijI(10,
     -   2)*P(549)-EijI(10,3)*P(550))
       F(55)=DCMPLX(FR(55),FI(55))
       FR(56) = -8*(EijR(2,3)-EijR(2,4)+EijR(5,3)+EijR(5,4)-EijR(6,3)-
     -   EijR(6,4)-EijR(9,3)+EijR(9,4)-2*(EijR(8,3)-EijR(17,3))+3*(Ei
     -   jR(8,4)-EijR(17,4)+EijR(23,4)-EijR(26,4)))
       FI(56) = -8*(EijI(2,3)-EijI(2,4)+EijI(5,3)+EijI(5,4)-EijI(6,3)-
     -   EijI(6,4)-EijI(9,3)+EijI(9,4)-2*(EijI(8,3)-EijI(17,3))+3*(Ei
     -   jI(8,4)-EijI(17,4)+EijI(23,4)-EijI(26,4)))
       F(56)=DCMPLX(FR(56),FI(56))
       FR(57) = 8*(EijR(2,4)-EijR(6,3)+EijR(7,3)-EijR(9,3)-EijR(9,4)+E
     -   ijR(10,3)-EijR(10,4)+EijR(17,4)-EijR(23,4)-EijR(24,4)+EijR(2
     -   5,4)+EijR(28,4)-2*(EijR(8,4)-EijR(17,3)+EijR(18,3)-EijR(26,4
     -   )-EijR(27,4)+EijR(35,4)))
       FI(57) = 8*(EijI(2,4)-EijI(6,3)+EijI(7,3)-EijI(9,3)-EijI(9,4)+E
     -   ijI(10,3)-EijI(10,4)+EijI(17,4)-EijI(23,4)-EijI(24,4)+EijI(2
     -   5,4)+EijI(28,4)-2*(EijI(8,4)-EijI(17,3)+EijI(18,3)-EijI(26,4
     -   )-EijI(27,4)+EijI(35,4)))
       F(57)=DCMPLX(FR(57),FI(57))
       FR(58) = 8*(EijR(2,4)-EijR(6,3)-EijR(9,3)-EijR(9,4)+EijR(17,4)-
     -   EijR(23,4)-2*(EijR(8,4)-EijR(17,3)-EijR(26,4)))
       FI(58) = 8*(EijI(2,4)-EijI(6,3)-EijI(9,3)-EijI(9,4)+EijI(17,4)-
     -   EijI(23,4)-2*(EijI(8,4)-EijI(17,3)-EijI(26,4)))
       F(58)=DCMPLX(FR(58),FI(58))
       FR(59) = -8*(EijR(2,2)-EijR(2,4)-EijR(5,2)-EijR(5,3)+EijR(6,2)+
     -   EijR(6,3)-EijR(8,2)+EijR(8,3)-EijR(9,3)-EijR(11,3)+EijR(12,3
     -   )-EijR(17,4)-EijR(18,4)-EijR(20,4)-4*EijR(26,4)+2*(EijR(8,4)
     -   +EijR(9,4)+EijR(23,4)+EijR(29,4)))
       FI(59) = -8*(EijI(2,2)-EijI(2,4)-EijI(5,2)-EijI(5,3)+EijI(6,2)+
     -   EijI(6,3)-EijI(8,2)+EijI(8,3)-EijI(9,3)-EijI(11,3)+EijI(12,3
     -   )-EijI(17,4)-EijI(18,4)-EijI(20,4)-4*EijI(26,4)+2*(EijI(8,4)
     -   +EijI(9,4)+EijI(23,4)+EijI(29,4)))
       F(59)=DCMPLX(FR(59),FI(59))
       FR(60) = 8*(EijR(2,3)+EijR(2,4)-EijR(8,3)-EijR(8,4)-EijR(11,3)+
     -   EijR(11,4)+EijR(12,3)-EijR(12,4)-2*(EijR(9,3)-EijR(17,3))-3*
     -   (EijR(9,4)-EijR(20,4)-EijR(26,4)+EijR(29,4)))
       FI(60) = 8*(EijI(2,3)+EijI(2,4)-EijI(8,3)-EijI(8,4)-EijI(11,3)+
     -   EijI(11,4)+EijI(12,3)-EijI(12,4)-2*(EijI(9,3)-EijI(17,3))-3*
     -   (EijI(9,4)-EijI(20,4)-EijI(26,4)+EijI(29,4)))
       F(60)=DCMPLX(FR(60),FI(60))
       FR(61) = 8*(EijR(2,3)+EijR(2,4)-EijR(8,3)-EijR(8,4)-EijR(10,4)-
     -   EijR(11,3)+EijR(12,3)+EijR(20,4)+EijR(27,4)-EijR(29,4)+EijR(
     -   30,4)-EijR(31,4)-2*(EijR(9,3)+EijR(9,4)-EijR(17,3)-EijR(26,4
     -   )-EijR(28,4)+EijR(35,4)))
       FI(61) = 8*(EijI(2,3)+EijI(2,4)-EijI(8,3)-EijI(8,4)-EijI(10,4)-
     -   EijI(11,3)+EijI(12,3)+EijI(20,4)+EijI(27,4)-EijI(29,4)+EijI(
     -   30,4)-EijI(31,4)-2*(EijI(9,3)+EijI(9,4)-EijI(17,3)-EijI(26,4
     -   )-EijI(28,4)+EijI(35,4)))
       F(61)=DCMPLX(FR(61),FI(61))
       FR(62) = 8*(EijR(2,3)+EijR(2,4)-EijR(8,3)-EijR(8,4)-EijR(11,3)+
     -   EijR(12,3)+EijR(20,4)-2*(EijR(9,3)+EijR(9,4)-EijR(17,3)-EijR
     -   (26,4))-EijR(29,4))
       FI(62) = 8*(EijI(2,3)+EijI(2,4)-EijI(8,3)-EijI(8,4)-EijI(11,3)+
     -   EijI(12,3)+EijI(20,4)-2*(EijI(9,3)+EijI(9,4)-EijI(17,3)-EijI
     -   (26,4))-EijI(29,4))
       F(62)=DCMPLX(FR(62),FI(62))
       FR(63) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)-EijR(5,2)+EijR(6,2)+
     -   EijR(6,3)+EijR(7,3)-EijR(8,2)+3*EijR(8,3)+EijR(9,4)+EijR(10,
     -   4)-EijR(17,3)-EijR(17,4)-EijR(18,3)-EijR(19,3)+EijR(20,3)+Ei
     -   jR(23,4)+EijR(24,4)-EijR(25,4)-EijR(28,4)-2*(EijR(5,3)-EijR(
     -   8,4)+EijR(26,4)+EijR(27,4)-EijR(35,4)))
       FI(63) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)-EijI(5,2)+EijI(6,2)+
     -   EijI(6,3)+EijI(7,3)-EijI(8,2)+3*EijI(8,3)+EijI(9,4)+EijI(10,
     -   4)-EijI(17,3)-EijI(17,4)-EijI(18,3)-EijI(19,3)+EijI(20,3)+Ei
     -   jI(23,4)+EijI(24,4)-EijI(25,4)-EijI(28,4)-2*(EijI(5,3)-EijI(
     -   8,4)+EijI(26,4)+EijI(27,4)-EijI(35,4)))
       F(63)=DCMPLX(FR(63),FI(63))
       FR(64) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-E
     -   ijR(8,4)-EijR(10,3)-EijR(10,4)-EijR(11,3)+EijR(12,3)-3*(EijR
     -   (9,3)-EijR(17,3))+EijR(18,3)-EijR(19,3)+EijR(20,3)+EijR(20,4
     -   )+EijR(27,4)-EijR(29,4)+EijR(30,4)-EijR(31,4)+2*(EijR(2,3)-E
     -   ijR(8,3)-EijR(9,4)+EijR(26,4)+EijR(28,4)-EijR(35,4)))
       FI(64) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-E
     -   ijI(8,4)-EijI(10,3)-EijI(10,4)-EijI(11,3)+EijI(12,3)-3*(EijI
     -   (9,3)-EijI(17,3))+EijI(18,3)-EijI(19,3)+EijI(20,3)+EijI(20,4
     -   )+EijI(27,4)-EijI(29,4)+EijI(30,4)-EijI(31,4)+2*(EijI(2,3)-E
     -   ijI(8,3)-EijI(9,4)+EijI(26,4)+EijI(28,4)-EijI(35,4)))
       F(64)=DCMPLX(FR(64),FI(64))
       FR(65) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-E
     -   ijR(8,4)-EijR(9,4)+EijR(21,4)+EijR(26,4)-EijR(32,4)+EijR(33,
     -   4)-EijR(34,4)+2*(EijR(2,3)-EijR(8,3)-EijR(9,3)-EijR(10,3)-Ei
     -   jR(10,4)+EijR(17,3)+EijR(18,3)-EijR(19,3)+EijR(20,3)+EijR(27
     -   ,4)+EijR(28,4)-EijR(35,4)))
       FI(65) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-E
     -   ijI(8,4)-EijI(9,4)+EijI(21,4)+EijI(26,4)-EijI(32,4)+EijI(33,
     -   4)-EijI(34,4)+2*(EijI(2,3)-EijI(8,3)-EijI(9,3)-EijI(10,3)-Ei
     -   jI(10,4)+EijI(17,3)+EijI(18,3)-EijI(19,3)+EijI(20,3)+EijI(27
     -   ,4)+EijI(28,4)-EijI(35,4)))
       F(65)=DCMPLX(FR(65),FI(65))
       FR(66) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(7,2)-E
     -   ijR(8,2)-EijR(8,4)+EijR(9,2)-EijR(9,4)-EijR(10,3)-EijR(10,4)
     -   +2*(EijR(2,3)-EijR(8,3)-EijR(9,3)+EijR(17,3))+EijR(18,3)-Eij
     -   R(19,3)+EijR(20,3)+EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,
     -   4))
       FI(66) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(7,2)-E
     -   ijI(8,2)-EijI(8,4)+EijI(9,2)-EijI(9,4)-EijI(10,3)-EijI(10,4)
     -   +2*(EijI(2,3)-EijI(8,3)-EijI(9,3)+EijI(17,3))+EijI(18,3)-Eij
     -   I(19,3)+EijI(20,3)+EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,
     -   4))
       F(66)=DCMPLX(FR(66),FI(66))
       FR(67) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)-EijR(5,2)+EijR(6,2)+
     -   EijR(6,3)-EijR(8,2)+3*EijR(8,3)+EijR(9,4)-EijR(17,3)-EijR(17
     -   ,4)+EijR(23,4)-2*(EijR(5,3)-EijR(8,4)+EijR(26,4)))
       FI(67) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)-EijI(5,2)+EijI(6,2)+
     -   EijI(6,3)-EijI(8,2)+3*EijI(8,3)+EijI(9,4)-EijI(17,3)-EijI(17
     -   ,4)+EijI(23,4)-2*(EijI(5,3)-EijI(8,4)+EijI(26,4)))
       F(67)=DCMPLX(FR(67),FI(67))
       FR(68) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-E
     -   ijR(8,4)-EijR(11,3)+EijR(12,3)-3*(EijR(9,3)-EijR(17,3))+EijR
     -   (20,4)+2*(EijR(2,3)-EijR(8,3)-EijR(9,4)+EijR(26,4))-EijR(29,
     -   4))
       FI(68) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-E
     -   ijI(8,4)-EijI(11,3)+EijI(12,3)-3*(EijI(9,3)-EijI(17,3))+EijI
     -   (20,4)+2*(EijI(2,3)-EijI(8,3)-EijI(9,4)+EijI(26,4))-EijI(29,
     -   4))
       F(68)=DCMPLX(FR(68),FI(68))
       FR(69) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)+EijR(7,2)-E
     -   ijR(8,2)-EijR(8,4)-EijR(9,2)-EijR(9,4)-EijR(10,3)-EijR(10,4)
     -   +2*(EijR(2,3)-EijR(8,3)-EijR(9,3)+EijR(17,3))+EijR(18,3)-Eij
     -   R(19,3)+EijR(20,3)+EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,
     -   4))
       FI(69) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)+EijI(7,2)-E
     -   ijI(8,2)-EijI(8,4)-EijI(9,2)-EijI(9,4)-EijI(10,3)-EijI(10,4)
     -   +2*(EijI(2,3)-EijI(8,3)-EijI(9,3)+EijI(17,3))+EijI(18,3)-Eij
     -   I(19,3)+EijI(20,3)+EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,
     -   4))
       F(69)=DCMPLX(FR(69),FI(69))
       FR(70) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-E
     -   ijR(8,4)-EijR(9,4)+2*(EijR(2,3)-EijR(8,3)-EijR(9,3)+EijR(17,
     -   3))+EijR(26,4))
       FI(70) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-E
     -   ijI(8,4)-EijI(9,4)+2*(EijI(2,3)-EijI(8,3)-EijI(9,3)+EijI(17,
     -   3))+EijI(26,4))
       F(70)=DCMPLX(FR(70),FI(70))
       FR(71) = -4*(s45*(-EijR(2,1)+EijR(3,1)-EijR(5,2)+EijR(6,2))+s15
     -   *(EijR(9,2)-EijR(10,2))+s23*(EijR(2,1)-EijR(3,1)+EijR(5,2)-E
     -   ijR(6,2)-EijR(9,2)+EijR(10,2))-2*(EijR(22,3)-EijR(23,3)+EijR
     -   (37,4)+EijR(38,4)-2*EijR(43,4)))
       FI(71) = -4*(s45*(-EijI(2,1)+EijI(3,1)-EijI(5,2)+EijI(6,2))+s15
     -   *(EijI(9,2)-EijI(10,2))+s23*(EijI(2,1)-EijI(3,1)+EijI(5,2)-E
     -   ijI(6,2)-EijI(9,2)+EijI(10,2))-2*(EijI(22,3)-EijI(23,3)+EijI
     -   (37,4)+EijI(38,4)-2*EijI(43,4)))
       F(71)=DCMPLX(FR(71),FI(71))
       FR(72) = -8*(EijR(2,3)-EijR(2,4)-EijR(8,3)-EijR(11,3)+EijR(12,3
     -   )-EijR(17,4)-EijR(18,4)-EijR(20,4)-4*EijR(26,4)+2*(EijR(8,4)
     -   -EijR(9,3)+EijR(9,4)+EijR(17,3)+EijR(23,4)+EijR(29,4)))
       FI(72) = -8*(EijI(2,3)-EijI(2,4)-EijI(8,3)-EijI(11,3)+EijI(12,3
     -   )-EijI(17,4)-EijI(18,4)-EijI(20,4)-4*EijI(26,4)+2*(EijI(8,4)
     -   -EijI(9,3)+EijI(9,4)+EijI(17,3)+EijI(23,4)+EijI(29,4)))
       F(72)=DCMPLX(FR(72),FI(72))
       FR(73) = 8*(EijR(2,4)-EijR(8,4)-EijR(9,3)+EijR(10,3)-EijR(10,4)
     -   -EijR(11,3)+EijR(12,3)+EijR(17,3)-EijR(18,3)+EijR(19,3)-EijR
     -   (20,3)+EijR(20,4)+EijR(27,4)-EijR(29,4)+EijR(30,4)-EijR(31,4
     -   )-2*(EijR(9,4)-EijR(26,4)-EijR(28,4)+EijR(35,4)))
       FI(73) = 8*(EijI(2,4)-EijI(8,4)-EijI(9,3)+EijI(10,3)-EijI(10,4)
     -   -EijI(11,3)+EijI(12,3)+EijI(17,3)-EijI(18,3)+EijI(19,3)-EijI
     -   (20,3)+EijI(20,4)+EijI(27,4)-EijI(29,4)+EijI(30,4)-EijI(31,4
     -   )-2*(EijI(9,4)-EijI(26,4)-EijI(28,4)+EijI(35,4)))
       F(73)=DCMPLX(FR(73),FI(73))
       FR(74) = 8*(EijR(2,4)-EijR(8,4)-EijR(9,3)-EijR(11,3)+EijR(12,3)
     -   +EijR(17,3)+EijR(20,4)-2*(EijR(9,4)-EijR(26,4))-EijR(29,4))
       FI(74) = 8*(EijI(2,4)-EijI(8,4)-EijI(9,3)-EijI(11,3)+EijI(12,3)
     -   +EijI(17,3)+EijI(20,4)-2*(EijI(9,4)-EijI(26,4))-EijI(29,4))
       F(74)=DCMPLX(FR(74),FI(74))
       FR(75) = -8*(EijR(2,2)-EijR(2,4)+EijR(3,2)-EijR(3,3)+EijR(8,3)+
     -   EijR(8,4)-EijR(9,3)+EijR(11,3)-EijR(11,4)+EijR(12,4)-2*(EijR
     -   (8,2)-EijR(12,3)+EijR(17,3))+3*(EijR(9,4)-EijR(20,4)-EijR(26
     -   ,4)+EijR(29,4)))
       FI(75) = -8*(EijI(2,2)-EijI(2,4)+EijI(3,2)-EijI(3,3)+EijI(8,3)+
     -   EijI(8,4)-EijI(9,3)+EijI(11,3)-EijI(11,4)+EijI(12,4)-2*(EijI
     -   (8,2)-EijI(12,3)+EijI(17,3))+3*(EijI(9,4)-EijI(20,4)-EijI(26
     -   ,4)+EijI(29,4)))
       F(75)=DCMPLX(FR(75),FI(75))
       FR(76) = 8*(EijR(2,3)+EijR(2,4)-EijR(3,3)+EijR(3,4)-3*(EijR(9,3
     -   )-EijR(12,3))-4*(EijR(9,4)+EijR(12,4))+6*EijR(20,4))
       FI(76) = 8*(EijI(2,3)+EijI(2,4)-EijI(3,3)+EijI(3,4)-3*(EijI(9,3
     -   )-EijI(12,3))-4*(EijI(9,4)+EijI(12,4))+6*EijI(20,4))
       F(76)=DCMPLX(FR(76),FI(76))
       FR(77) = 8*(EijR(2,3)+EijR(2,4)-EijR(3,3)-EijR(10,4)-EijR(12,4)
     -   +EijR(13,4)-3*(EijR(9,3)+EijR(9,4)-EijR(12,3)-EijR(20,4)-Eij
     -   R(28,4)+EijR(31,4)))
       FI(77) = 8*(EijI(2,3)+EijI(2,4)-EijI(3,3)-EijI(10,4)-EijI(12,4)
     -   +EijI(13,4)-3*(EijI(9,3)+EijI(9,4)-EijI(12,3)-EijI(20,4)-Eij
     -   I(28,4)+EijI(31,4)))
       F(77)=DCMPLX(FR(77),FI(77))
       FR(78) = 8*(EijR(2,3)+EijR(2,4)-EijR(3,3)-EijR(12,4)-3*(EijR(9,
     -   3)+EijR(9,4)-EijR(12,3)-EijR(20,4)))
       FI(78) = 8*(EijI(2,3)+EijI(2,4)-EijI(3,3)-EijI(12,4)-3*(EijI(9,
     -   3)+EijI(9,4)-EijI(12,3)-EijI(20,4)))
       F(78)=DCMPLX(FR(78),FI(78))
       FR(79) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)+EijR(3,2)+EijR(8,4)+
     -   EijR(9,3)+EijR(10,4)+EijR(11,3)-EijR(13,3)-3*EijR(17,3)-EijR
     -   (18,3)+EijR(19,3)+EijR(20,3)-EijR(20,4)-EijR(27,4)+EijR(29,4
     -   )-EijR(30,4)+EijR(31,4)-2*(EijR(8,2)-EijR(8,3)-EijR(9,4)+Eij
     -   R(26,4)+EijR(28,4)-EijR(35,4)))
       FI(79) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)+EijI(3,2)+EijI(8,4)+
     -   EijI(9,3)+EijI(10,4)+EijI(11,3)-EijI(13,3)-3*EijI(17,3)-EijI
     -   (18,3)+EijI(19,3)+EijI(20,3)-EijI(20,4)-EijI(27,4)+EijI(29,4
     -   )-EijI(30,4)+EijI(31,4)-2*(EijI(8,2)-EijI(8,3)-EijI(9,4)+Eij
     -   I(26,4)+EijI(28,4)-EijI(35,4)))
       F(79)=DCMPLX(FR(79),FI(79))
       FR(80) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,2)-EijR(3,3)-5*EijR(9,3)
     -   -EijR(10,3)-EijR(10,4)+4*EijR(12,3)-EijR(12,4)-EijR(13,3)+Ei
     -   jR(13,4)+2*(EijR(2,3)-EijR(8,2)+EijR(20,3))-3*(EijR(9,4)-Eij
     -   R(20,4)-EijR(28,4)+EijR(31,4)))
       FI(80) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,2)-EijI(3,3)-5*EijI(9,3)
     -   -EijI(10,3)-EijI(10,4)+4*EijI(12,3)-EijI(12,4)-EijI(13,3)+Ei
     -   jI(13,4)+2*(EijI(2,3)-EijI(8,2)+EijI(20,3))-3*(EijI(9,4)-Eij
     -   I(20,4)-EijI(28,4)+EijI(31,4)))
       F(80)=DCMPLX(FR(80),FI(80))
       FR(81) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,2)+EijR(20,4)+EijR(21,4)
     -   +EijR(22,4)-4*(EijR(9,3)-EijR(20,3)-EijR(28,4))+2*(EijR(2,3)
     -   -EijR(8,2)-EijR(9,4)-EijR(10,3)-EijR(10,4)+EijR(12,3)-EijR(1
     -   3,3)-EijR(31,4)-EijR(34,4)))
       FI(81) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,2)+EijI(20,4)+EijI(21,4)
     -   +EijI(22,4)-4*(EijI(9,3)-EijI(20,3)-EijI(28,4))+2*(EijI(2,3)
     -   -EijI(8,2)-EijI(9,4)-EijI(10,3)-EijI(10,4)+EijI(12,3)-EijI(1
     -   3,3)-EijI(31,4)-EijI(34,4)))
       F(81)=DCMPLX(FR(81),FI(81))
       FR(82) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,2)+EijR(9,2)-4*EijR(9,3)
     -   -EijR(10,2)-EijR(10,3)-EijR(10,4)-EijR(13,3)+EijR(20,4)+2*(E
     -   ijR(2,3)-EijR(8,2)-EijR(9,4)+EijR(12,3)+EijR(20,3)+EijR(28,4
     -   ))-EijR(31,4))
       FI(82) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,2)+EijI(9,2)-4*EijI(9,3)
     -   -EijI(10,2)-EijI(10,3)-EijI(10,4)-EijI(13,3)+EijI(20,4)+2*(E
     -   ijI(2,3)-EijI(8,2)-EijI(9,4)+EijI(12,3)+EijI(20,3)+EijI(28,4
     -   ))-EijI(31,4))
       F(82)=DCMPLX(FR(82),FI(82))
       FR(83) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)+EijR(3,2)+EijR(8,4)+
     -   EijR(9,3)+EijR(11,3)-3*EijR(17,3)-EijR(20,4)-2*(EijR(8,2)-Ei
     -   jR(8,3)-EijR(9,4)+EijR(26,4))+EijR(29,4))
       FI(83) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)+EijI(3,2)+EijI(8,4)+
     -   EijI(9,3)+EijI(11,3)-3*EijI(17,3)-EijI(20,4)-2*(EijI(8,2)-Ei
     -   jI(8,3)-EijI(9,4)+EijI(26,4))+EijI(29,4))
       F(83)=DCMPLX(FR(83),FI(83))
       FR(84) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,2)-EijR(3,3)+2*(EijR(2,3
     -   )-EijR(8,2))-5*EijR(9,3)+4*EijR(12,3)-EijR(12,4)-3*(EijR(9,4
     -   )-EijR(20,4)))
       FI(84) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,2)-EijI(3,3)+2*(EijI(2,3
     -   )-EijI(8,2))-5*EijI(9,3)+4*EijI(12,3)-EijI(12,4)-3*(EijI(9,4
     -   )-EijI(20,4)))
       F(84)=DCMPLX(FR(84),FI(84))
       FR(85) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,2)-EijR(9,2)-4*EijR(9,3)
     -   +EijR(10,2)-EijR(10,3)-EijR(10,4)-EijR(13,3)+EijR(20,4)+2*(E
     -   ijR(2,3)-EijR(8,2)-EijR(9,4)+EijR(12,3)+EijR(20,3)+EijR(28,4
     -   ))-EijR(31,4))
       FI(85) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,2)-EijI(9,2)-4*EijI(9,3)
     -   +EijI(10,2)-EijI(10,3)-EijI(10,4)-EijI(13,3)+EijI(20,4)+2*(E
     -   ijI(2,3)-EijI(8,2)-EijI(9,4)+EijI(12,3)+EijI(20,3)+EijI(28,4
     -   ))-EijI(31,4))
       F(85)=DCMPLX(FR(85),FI(85))
       FR(86) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,2)+2*(EijR(2,3)-EijR(8,2
     -   )-2*EijR(9,3)-EijR(9,4)+EijR(12,3))+EijR(20,4))
       FI(86) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,2)+2*(EijI(2,3)-EijI(8,2
     -   )-2*EijI(9,3)-EijI(9,4)+EijI(12,3))+EijI(20,4))
       F(86)=DCMPLX(FR(86),FI(86))
       FR(87) = -8*(EijR(2,3)-EijR(2,4)-EijR(8,3)-EijR(9,3)+EijR(9,4)-
     -   EijR(10,3)+EijR(10,4)+EijR(17,3)-EijR(17,4)+EijR(18,3)-EijR(
     -   19,3)+EijR(20,3)+EijR(23,4)+EijR(24,4)-EijR(25,4)-EijR(28,4)
     -   +2*(EijR(8,4)-EijR(26,4)-EijR(27,4)+EijR(35,4)))
       FI(87) = -8*(EijI(2,3)-EijI(2,4)-EijI(8,3)-EijI(9,3)+EijI(9,4)-
     -   EijI(10,3)+EijI(10,4)+EijI(17,3)-EijI(17,4)+EijI(18,3)-EijI(
     -   19,3)+EijI(20,3)+EijI(23,4)+EijI(24,4)-EijI(25,4)-EijI(28,4)
     -   +2*(EijI(8,4)-EijI(26,4)-EijI(27,4)+EijI(35,4)))
       F(87)=DCMPLX(FR(87),FI(87))
       FR(88) = 8*(EijR(2,4)-EijR(8,4)-EijR(9,3)-EijR(9,4)+EijR(10,3)+
     -   EijR(14,3)-EijR(15,3)+EijR(17,3)-EijR(18,3)-EijR(19,3)+EijR(
     -   20,3)+EijR(21,4)+EijR(26,4)-EijR(32,4)+EijR(33,4)-EijR(34,4)
     -   -2*(EijR(10,4)-EijR(27,4)-EijR(28,4)+EijR(35,4)))
       FI(88) = 8*(EijI(2,4)-EijI(8,4)-EijI(9,3)-EijI(9,4)+EijI(10,3)+
     -   EijI(14,3)-EijI(15,3)+EijI(17,3)-EijI(18,3)-EijI(19,3)+EijI(
     -   20,3)+EijI(21,4)+EijI(26,4)-EijI(32,4)+EijI(33,4)-EijI(34,4)
     -   -2*(EijI(10,4)-EijI(27,4)-EijI(28,4)+EijI(35,4)))
       F(88)=DCMPLX(FR(88),FI(88))
       FR(89) = 8*(EijR(2,4)-EijR(6,2)+EijR(7,2)+EijR(8,2)-EijR(8,4)-E
     -   ijR(9,2)-EijR(9,3)-EijR(9,4)-EijR(10,4)+EijR(17,3)-EijR(19,3
     -   )+EijR(20,3)+EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4))
       FI(89) = 8*(EijI(2,4)-EijI(6,2)+EijI(7,2)+EijI(8,2)-EijI(8,4)-E
     -   ijI(9,2)-EijI(9,3)-EijI(9,4)-EijI(10,4)+EijI(17,3)-EijI(19,3
     -   )+EijI(20,3)+EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4))
       F(89)=DCMPLX(FR(89),FI(89))
       FR(90) = -8*(EijR(2,2)-EijR(2,4)-EijR(8,2)+EijR(8,3)+EijR(8,4)-
     -   EijR(9,2)-EijR(9,3)+EijR(10,2)+EijR(10,4)+EijR(12,3)-EijR(13
     -   ,3)-EijR(17,3)-EijR(18,3)+EijR(19,3)+EijR(20,3)-EijR(20,4)-E
     -   ijR(27,4)+EijR(29,4)-EijR(30,4)+EijR(31,4)+2*(EijR(9,4)-EijR
     -   (26,4)-EijR(28,4)+EijR(35,4)))
       FI(90) = -8*(EijI(2,2)-EijI(2,4)-EijI(8,2)+EijI(8,3)+EijI(8,4)-
     -   EijI(9,2)-EijI(9,3)+EijI(10,2)+EijI(10,4)+EijI(12,3)-EijI(13
     -   ,3)-EijI(17,3)-EijI(18,3)+EijI(19,3)+EijI(20,3)-EijI(20,4)-E
     -   ijI(27,4)+EijI(29,4)-EijI(30,4)+EijI(31,4)+2*(EijI(9,4)-EijI
     -   (26,4)-EijI(28,4)+EijI(35,4)))
       F(90)=DCMPLX(FR(90),FI(90))
       FR(91) = 8*(EijR(2,3)+EijR(2,4)-EijR(10,3)-EijR(10,4)+EijR(12,3
     -   )-EijR(12,4)-EijR(13,3)+EijR(13,4)-2*(EijR(9,3)-EijR(20,3))-
     -   3*(EijR(9,4)-EijR(20,4)-EijR(28,4)+EijR(31,4)))
       FI(91) = 8*(EijI(2,3)+EijI(2,4)-EijI(10,3)-EijI(10,4)+EijI(12,3
     -   )-EijI(12,4)-EijI(13,3)+EijI(13,4)-2*(EijI(9,3)-EijI(20,3))-
     -   3*(EijI(9,4)-EijI(20,4)-EijI(28,4)+EijI(31,4)))
       F(91)=DCMPLX(FR(91),FI(91))
       FR(92) = 8*(EijR(2,3)+EijR(2,4)-EijR(10,3)+EijR(12,3)-EijR(13,3
     -   )+EijR(20,4)+EijR(21,4)+EijR(22,4)+4*EijR(28,4)-2*(EijR(9,3)
     -   +EijR(9,4)+EijR(10,4)-EijR(20,3)+EijR(31,4)+EijR(34,4)))
       FI(92) = 8*(EijI(2,3)+EijI(2,4)-EijI(10,3)+EijI(12,3)-EijI(13,3
     -   )+EijI(20,4)+EijI(21,4)+EijI(22,4)+4*EijI(28,4)-2*(EijI(9,3)
     -   +EijI(9,4)+EijI(10,4)-EijI(20,3)+EijI(31,4)+EijI(34,4)))
       F(92)=DCMPLX(FR(92),FI(92))
       FR(93) = 8*(EijR(2,3)+EijR(2,4)-EijR(3,2)+EijR(8,2)-EijR(9,2)+E
     -   ijR(10,2)-EijR(10,3)-EijR(10,4)+EijR(12,3)-EijR(13,3)+EijR(2
     -   0,4)-2*(EijR(9,3)+EijR(9,4)-EijR(20,3)-EijR(28,4))-EijR(31,4
     -   ))
       FI(93) = 8*(EijI(2,3)+EijI(2,4)-EijI(3,2)+EijI(8,2)-EijI(9,2)+E
     -   ijI(10,2)-EijI(10,3)-EijI(10,4)+EijI(12,3)-EijI(13,3)+EijI(2
     -   0,4)-2*(EijI(9,3)+EijI(9,4)-EijI(20,3)-EijI(28,4))-EijI(31,4
     -   ))
       F(93)=DCMPLX(FR(93),FI(93))
       FR(94) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)-EijR(8,2)+EijR(8,4)-
     -   EijR(9,2)+EijR(9,4)+EijR(10,2)+EijR(10,3)+EijR(14,3)-EijR(16
     -   ,3)-EijR(17,3)-3*EijR(18,3)+EijR(19,3)+EijR(20,3)-EijR(21,4)
     -   -EijR(26,4)+EijR(32,4)-EijR(33,4)+EijR(34,4)+2*(EijR(8,3)+Ei
     -   jR(10,4)-EijR(27,4)-EijR(28,4)+EijR(35,4)))
       FI(94) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)-EijI(8,2)+EijI(8,4)-
     -   EijI(9,2)+EijI(9,4)+EijI(10,2)+EijI(10,3)+EijI(14,3)-EijI(16
     -   ,3)-EijI(17,3)-3*EijI(18,3)+EijI(19,3)+EijI(20,3)-EijI(21,4)
     -   -EijI(26,4)+EijI(32,4)-EijI(33,4)+EijI(34,4)+2*(EijI(8,3)+Ei
     -   jI(10,4)-EijI(27,4)-EijI(28,4)+EijI(35,4)))
       F(94)=DCMPLX(FR(94),FI(94))
       FR(95) = 8*(EijR(2,2)+EijR(2,4)-EijR(8,2)-EijR(9,2)+EijR(10,2)-
     -   3*(EijR(9,3)+EijR(10,3))+EijR(12,3)-EijR(13,3)+EijR(15,3)-Ei
     -   jR(16,3)+EijR(20,4)+EijR(21,4)+EijR(22,4)+4*(EijR(20,3)+EijR
     -   (28,4))+2*(EijR(2,3)-EijR(9,4)-EijR(10,4)-EijR(31,4)-EijR(34
     -   ,4)))
       FI(95) = 8*(EijI(2,2)+EijI(2,4)-EijI(8,2)-EijI(9,2)+EijI(10,2)-
     -   3*(EijI(9,3)+EijI(10,3))+EijI(12,3)-EijI(13,3)+EijI(15,3)-Ei
     -   jI(16,3)+EijI(20,4)+EijI(21,4)+EijI(22,4)+4*(EijI(20,3)+EijI
     -   (28,4))+2*(EijI(2,3)-EijI(9,4)-EijI(10,4)-EijI(31,4)-EijI(34
     -   ,4)))
       F(95)=DCMPLX(FR(95),FI(95))
       FR(96) = 8*(EijR(2,2)+EijR(2,4)-EijR(8,2)-EijR(9,2)-EijR(9,4)+E
     -   ijR(10,2)-EijR(15,4)+2*(EijR(2,3)-EijR(9,3)+EijR(15,3)-EijR(
     -   16,3))+EijR(16,4)-4*(EijR(10,3)-EijR(20,3))-3*(EijR(10,4)-Ei
     -   jR(21,4)-EijR(28,4)+EijR(34,4)))
       FI(96) = 8*(EijI(2,2)+EijI(2,4)-EijI(8,2)-EijI(9,2)-EijI(9,4)+E
     -   ijI(10,2)-EijI(15,4)+2*(EijI(2,3)-EijI(9,3)+EijI(15,3)-EijI(
     -   16,3))+EijI(16,4)-4*(EijI(10,3)-EijI(20,3))-3*(EijI(10,4)-Ei
     -   jI(21,4)-EijI(28,4)+EijI(34,4)))
       F(96)=DCMPLX(FR(96),FI(96))
       FR(97) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,1)-EijR(4,1)-EijR(9,2)-E
     -   ijR(9,4)+EijR(15,3)-EijR(16,3)-3*(EijR(10,3)-EijR(20,3))+Eij
     -   R(21,4)+2*(EijR(2,3)-EijR(9,3)-EijR(10,4)+EijR(28,4))-EijR(3
     -   4,4))
       FI(97) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,1)-EijI(4,1)-EijI(9,2)-E
     -   ijI(9,4)+EijI(15,3)-EijI(16,3)-3*(EijI(10,3)-EijI(20,3))+Eij
     -   I(21,4)+2*(EijI(2,3)-EijI(9,3)-EijI(10,4)+EijI(28,4))-EijI(3
     -   4,4))
       F(97)=DCMPLX(FR(97),FI(97))
       FR(98) = -8*(-EijR(2,3)-EijR(2,4)-EijR(5,2)+EijR(7,2)-EijR(8,2)
     -   +EijR(8,4)+EijR(9,4)+EijR(10,2)+EijR(10,3)+EijR(10,4)-EijR(1
     -   7,3)+2*(EijR(2,2)+EijR(8,3)-EijR(9,2)-EijR(18,3))+EijR(19,3)
     -   -EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(35,4))
       FI(98) = -8*(-EijI(2,3)-EijI(2,4)-EijI(5,2)+EijI(7,2)-EijI(8,2)
     -   +EijI(8,4)+EijI(9,4)+EijI(10,2)+EijI(10,3)+EijI(10,4)-EijI(1
     -   7,3)+2*(EijI(2,2)+EijI(8,3)-EijI(9,2)-EijI(18,3))+EijI(19,3)
     -   -EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(35,4))
       F(98)=DCMPLX(FR(98),FI(98))
       FR(99) = -8*(EijR(2,1)-EijR(2,4)-EijR(3,1)+EijR(10,4)-EijR(12,3
     -   )+EijR(13,3)+3*(EijR(9,3)-EijR(20,3))-EijR(20,4)-2*(EijR(2,3
     -   )-EijR(9,4)-EijR(10,3)+EijR(28,4))+EijR(31,4))
       FI(99) = -8*(EijI(2,1)-EijI(2,4)-EijI(3,1)+EijI(10,4)-EijI(12,3
     -   )+EijI(13,3)+3*(EijI(9,3)-EijI(20,3))-EijI(20,4)-2*(EijI(2,3
     -   )-EijI(9,4)-EijI(10,3)+EijI(28,4))+EijI(31,4))
       F(99)=DCMPLX(FR(99),FI(99))
       FR(100) = -8*(EijR(2,1)-EijR(2,4)-EijR(4,1)+EijR(8,2)+EijR(9,4)
     -   -EijR(10,2)-EijR(15,3)+EijR(16,3)+3*(EijR(10,3)-EijR(20,3))-
     -   EijR(21,4)-2*(EijR(2,3)-EijR(9,3)-EijR(10,4)+EijR(28,4))+Eij
     -   R(34,4))
       FI(100) = -8*(EijI(2,1)-EijI(2,4)-EijI(4,1)+EijI(8,2)+EijI(9,4)
     -   -EijI(10,2)-EijI(15,3)+EijI(16,3)+3*(EijI(10,3)-EijI(20,3))-
     -   EijI(21,4)-2*(EijI(2,3)-EijI(9,3)-EijI(10,4)+EijI(28,4))+Eij
     -   I(34,4))
       F(100)=DCMPLX(FR(100),FI(100))
       FR(101) = -8*(EijR(2,1)-EijR(2,4)-EijR(3,1)+EijR(9,2)+EijR(9,4)
     -   -EijR(10,2)+EijR(10,4)-2*(EijR(2,3)-EijR(9,3)-EijR(10,3)+Eij
     -   R(20,3))-EijR(28,4))
       FI(101) = -8*(EijI(2,1)-EijI(2,4)-EijI(3,1)+EijI(9,2)+EijI(9,4)
     -   -EijI(10,2)+EijI(10,4)-2*(EijI(2,3)-EijI(9,3)-EijI(10,3)+Eij
     -   I(20,3))-EijI(28,4))
       F(101)=DCMPLX(FR(101),FI(101))
       FR(102) = -4*(EijR(2,2)-EijR(5,2)+EijR(5,3)+EijR(6,2)-EijR(6,3)
     -   -EijR(8,2)-4*(EijR(8,3)-EijR(8,4))+3*(EijR(2,3)-EijR(9,3))-2
     -   *(EijR(2,4)-EijR(9,4)+EijR(17,4)-EijR(23,4)-2*(EijR(17,3)-Ei
     -   jR(26,4))))
       FI(102) = -4*(EijI(2,2)-EijI(5,2)+EijI(5,3)+EijI(6,2)-EijI(6,3)
     -   -EijI(8,2)-4*(EijI(8,3)-EijI(8,4))+3*(EijI(2,3)-EijI(9,3))-2
     -   *(EijI(2,4)-EijI(9,4)+EijI(17,4)-EijI(23,4)-2*(EijI(17,3)-Ei
     -   jI(26,4))))
       F(102)=DCMPLX(FR(102),FI(102))
       FR(103) = -4*(EijR(2,3)-EijR(8,3)-EijR(11,3)+EijR(12,3)+4*(EijR
     -   (9,4)-EijR(26,4))+2*(EijR(2,2)-EijR(2,4)-EijR(5,2)+EijR(6,2)
     -   -EijR(8,2)+EijR(8,4)-EijR(9,3)+EijR(17,3)-EijR(20,4)+EijR(29
     -   ,4)))
       FI(103) = -4*(EijI(2,3)-EijI(8,3)-EijI(11,3)+EijI(12,3)+4*(EijI
     -   (9,4)-EijI(26,4))+2*(EijI(2,2)-EijI(2,4)-EijI(5,2)+EijI(6,2)
     -   -EijI(8,2)+EijI(8,4)-EijI(9,3)+EijI(17,3)-EijI(20,4)+EijI(29
     -   ,4)))
       F(103)=DCMPLX(FR(103),FI(103))
       FR(104) = -4*(EijR(2,3)-EijR(6,2)+EijR(8,2)-EijR(8,3)+EijR(9,3)
     -   -EijR(17,3)+3*(EijR(7,2)-EijR(9,2)-EijR(10,3)+EijR(18,3))-Ei
     -   jR(19,3)+EijR(20,3)+2*(EijR(2,2)-EijR(2,4)-EijR(5,2)+EijR(8,
     -   4)+EijR(9,4)+EijR(10,4)-EijR(26,4)-EijR(27,4)-EijR(28,4)+Eij
     -   R(35,4)))
       FI(104) = -4*(EijI(2,3)-EijI(6,2)+EijI(8,2)-EijI(8,3)+EijI(9,3)
     -   -EijI(17,3)+3*(EijI(7,2)-EijI(9,2)-EijI(10,3)+EijI(18,3))-Ei
     -   jI(19,3)+EijI(20,3)+2*(EijI(2,2)-EijI(2,4)-EijI(5,2)+EijI(8,
     -   4)+EijI(9,4)+EijI(10,4)-EijI(26,4)-EijI(27,4)-EijI(28,4)+Eij
     -   I(35,4)))
       F(104)=DCMPLX(FR(104),FI(104))
       FR(105) = -4*(EijR(2,3)+EijR(6,2)-EijR(8,2)-EijR(8,3)+EijR(9,3)
     -   -EijR(17,3)+2*(EijR(2,2)-EijR(2,4)-EijR(5,2)+EijR(8,4)+EijR(
     -   9,4)-EijR(26,4)))
       FI(105) = -4*(EijI(2,3)+EijI(6,2)-EijI(8,2)-EijI(8,3)+EijI(9,3)
     -   -EijI(17,3)+2*(EijI(2,2)-EijI(2,4)-EijI(5,2)+EijI(8,4)+EijI(
     -   9,4)-EijI(26,4)))
       F(105)=DCMPLX(FR(105),FI(105))
       FR(106) = -4*(EijR(2,1)+EijR(2,3)-EijR(3,1)+EijR(3,2)+EijR(5,2)
     -   -EijR(6,2)+EijR(8,3)-EijR(11,3)-3*(EijR(8,2)-EijR(12,3))-4*(
     -   EijR(9,3)-EijR(9,4)+EijR(26,4))+2*(EijR(2,2)-EijR(2,4)+EijR(
     -   8,4)-EijR(20,4)+EijR(29,4)))
       FI(106) = -4*(EijI(2,1)+EijI(2,3)-EijI(3,1)+EijI(3,2)+EijI(5,2)
     -   -EijI(6,2)+EijI(8,3)-EijI(11,3)-3*(EijI(8,2)-EijI(12,3))-4*(
     -   EijI(9,3)-EijI(9,4)+EijI(26,4))+2*(EijI(2,2)-EijI(2,4)+EijI(
     -   8,4)-EijI(20,4)+EijI(29,4)))
       F(106)=DCMPLX(FR(106),FI(106))
       FR(107) = -4*(EijR(2,2)-EijR(2,3)+EijR(3,2)-EijR(3,3)+EijR(9,3)
     -   +EijR(12,3)-2*(EijR(2,4)+EijR(8,2)-EijR(12,4))+6*(EijR(9,4)-
     -   EijR(20,4)))
       FI(107) = -4*(EijI(2,2)-EijI(2,3)+EijI(3,2)-EijI(3,3)+EijI(9,3)
     -   +EijI(12,3)-2*(EijI(2,4)+EijI(8,2)-EijI(12,4))+6*(EijI(9,4)-
     -   EijI(20,4)))
       F(107)=DCMPLX(FR(107),FI(107))
       FR(108) = 4*(-EijR(2,2)+EijR(2,3)+EijR(3,2)+EijR(10,3)+EijR(12,
     -   3)+EijR(13,3)-4*(EijR(9,4)-EijR(28,4))+2*(EijR(2,4)+EijR(9,2
     -   )-EijR(9,3)-EijR(10,2)-EijR(10,4)-EijR(20,3)+EijR(20,4)-EijR
     -   (31,4)))
       FI(108) = 4*(-EijI(2,2)+EijI(2,3)+EijI(3,2)+EijI(10,3)+EijI(12,
     -   3)+EijI(13,3)-4*(EijI(9,4)-EijI(28,4))+2*(EijI(2,4)+EijI(9,2
     -   )-EijI(9,3)-EijI(10,2)-EijI(10,4)-EijI(20,3)+EijI(20,4)-EijI
     -   (31,4)))
       F(108)=DCMPLX(FR(108),FI(108))
       FR(109) = -4*(EijR(2,2)-EijR(2,3)+EijR(3,2)+4*EijR(9,4)-EijR(12
     -   ,3)-2*(EijR(2,4)+EijR(8,2)-EijR(9,3)+EijR(20,4)))
       FI(109) = -4*(EijI(2,2)-EijI(2,3)+EijI(3,2)+4*EijI(9,4)-EijI(12
     -   ,3)-2*(EijI(2,4)+EijI(8,2)-EijI(9,3)+EijI(20,4)))
       F(109)=DCMPLX(FR(109),FI(109))
       FR(110) = -4*(EijR(2,1)+EijR(2,2)-EijR(2,3)-EijR(3,1)+EijR(6,2)
     -   -4*EijR(8,2)-EijR(9,3)+EijR(10,2)-EijR(10,3)-EijR(17,3)-EijR
     -   (18,3)-EijR(19,3)-3*(EijR(7,2)-EijR(8,3)-EijR(20,3))-2*(EijR
     -   (2,4)-EijR(5,2)-EijR(8,4)-EijR(9,2)-EijR(9,4)-EijR(10,4)+Eij
     -   R(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4)))
       FI(110) = -4*(EijI(2,1)+EijI(2,2)-EijI(2,3)-EijI(3,1)+EijI(6,2)
     -   -4*EijI(8,2)-EijI(9,3)+EijI(10,2)-EijI(10,3)-EijI(17,3)-EijI
     -   (18,3)-EijI(19,3)-3*(EijI(7,2)-EijI(8,3)-EijI(20,3))-2*(EijI
     -   (2,4)-EijI(5,2)-EijI(8,4)-EijI(9,2)-EijI(9,4)-EijI(10,4)+Eij
     -   I(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4)))
       F(110)=DCMPLX(FR(110),FI(110))
       FR(111) = 4*(EijR(2,1)+3*EijR(2,3)-EijR(3,1)-EijR(3,2)-EijR(8,2
     -   )-EijR(9,2)+EijR(10,2)-EijR(10,3)+EijR(12,3)+EijR(13,3)-4*(E
     -   ijR(9,3)+EijR(9,4)-EijR(28,4))+2*(EijR(2,2)+EijR(2,4)-EijR(1
     -   0,4)+EijR(20,4)-EijR(31,4)))
       FI(111) = 4*(EijI(2,1)+3*EijI(2,3)-EijI(3,1)-EijI(3,2)-EijI(8,2
     -   )-EijI(9,2)+EijI(10,2)-EijI(10,3)+EijI(12,3)+EijI(13,3)-4*(E
     -   ijI(9,3)+EijI(9,4)-EijI(28,4))+2*(EijI(2,2)+EijI(2,4)-EijI(1
     -   0,4)+EijI(20,4)-EijI(31,4)))
       F(111)=DCMPLX(FR(111),FI(111))
       FR(112) = 4*(EijR(2,1)-EijR(3,1)+3*(EijR(2,3)-EijR(9,3))-EijR(1
     -   5,3)+EijR(16,3)-4*(EijR(10,4)-EijR(28,4))+2*(EijR(2,2)+EijR(
     -   2,4)-EijR(8,2)-EijR(9,4)-EijR(10,3)+EijR(20,3)+EijR(21,4)-Ei
     -   jR(34,4)))
       FI(112) = 4*(EijI(2,1)-EijI(3,1)+3*(EijI(2,3)-EijI(9,3))-EijI(1
     -   5,3)+EijI(16,3)-4*(EijI(10,4)-EijI(28,4))+2*(EijI(2,2)+EijI(
     -   2,4)-EijI(8,2)-EijI(9,4)-EijI(10,3)+EijI(20,3)+EijI(21,4)-Ei
     -   jI(34,4)))
       F(112)=DCMPLX(FR(112),FI(112))
       FR(113) = -4*(EijR(2,1)-EijR(2,3)-EijR(3,1)+EijR(6,2)-4*EijR(8,
     -   2)+3*(EijR(2,2)+EijR(8,3))-EijR(9,3)-EijR(17,3)-2*(EijR(2,4)
     -   -EijR(8,4)-EijR(9,4)+EijR(26,4)))
       FI(113) = -4*(EijI(2,1)-EijI(2,3)-EijI(3,1)+EijI(6,2)-4*EijI(8,
     -   2)+3*(EijI(2,2)+EijI(8,3))-EijI(9,3)-EijI(17,3)-2*(EijI(2,4)
     -   -EijI(8,4)-EijI(9,4)+EijI(26,4)))
       F(113)=DCMPLX(FR(113),FI(113))
       FR(114) = -4*(EijR(2,1)-3*EijR(2,3)-EijR(3,1)+EijR(3,2)-EijR(8,
     -   2)+4*(EijR(9,3)+EijR(9,4))-EijR(12,3)-2*(EijR(2,4)+EijR(20,4
     -   )))
       FI(114) = -4*(EijI(2,1)-3*EijI(2,3)-EijI(3,1)+EijI(3,2)-EijI(8,
     -   2)+4*(EijI(9,3)+EijI(9,4))-EijI(12,3)-2*(EijI(2,4)+EijI(20,4
     -   )))
       F(114)=DCMPLX(FR(114),FI(114))
       FR(115) = -4*(EijR(2,1)-EijR(3,1)-3*(EijR(2,3)-EijR(9,3))-2*(Ei
     -   jR(2,4)-EijR(9,4)))
       FI(115) = -4*(EijI(2,1)-EijI(3,1)-3*(EijI(2,3)-EijI(9,3))-2*(Ei
     -   jI(2,4)-EijI(9,4)))
       F(115)=DCMPLX(FR(115),FI(115))
       FR(116) = -2*(D02345R+Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(4,2)
     -   +p5sq*(-2*(EijR(1,1)+EijR(1,2)-EijR(2,1)-EijR(5,2))-EijR(7,2
     -   )+EijR(9,2)+EijR(10,3)-EijR(18,3)+EijR(19,3)-EijR(20,3))+s15
     -   *(EijR(7,2)+2*(EijR(1,1)+EijR(1,2)-EijR(2,1)-EijR(5,2)-EijR(
     -   6,2)+EijR(8,2))-EijR(9,2)-EijR(10,3)+EijR(18,3)-EijR(19,3)+E
     -   ijR(20,3))+2*(Dij2345R(1,1)-3*EijR(11,2)-EijR(22,3)+EijR(23,
     -   3)-2*(EijR(37,4)-EijR(40,4)+EijR(41,4)-EijR(43,4)))+(-EijR(2
     -   ,2)-EijR(2,3)+EijR(5,2)-EijR(6,2)+EijR(8,2)+EijR(8,3)+EijR(1
     -   1,3)-EijR(12,3)+2*(EijR(9,3)-EijR(17,3)))*P(143)+(EijR(6,2)-
     -   EijR(7,2)-EijR(8,2)+EijR(9,2)-EijR(9,3)+EijR(10,3)-EijR(11,3
     -   )+EijR(12,3)+EijR(17,3)-EijR(18,3)+EijR(19,3)-EijR(20,3))*P(
     -   193))
       FI(116) = -2*(D02345I+Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(4,2)
     -   +p5sq*(-2*(EijI(1,1)+EijI(1,2)-EijI(2,1)-EijI(5,2))-EijI(7,2
     -   )+EijI(9,2)+EijI(10,3)-EijI(18,3)+EijI(19,3)-EijI(20,3))+s15
     -   *(EijI(7,2)+2*(EijI(1,1)+EijI(1,2)-EijI(2,1)-EijI(5,2)-EijI(
     -   6,2)+EijI(8,2))-EijI(9,2)-EijI(10,3)+EijI(18,3)-EijI(19,3)+E
     -   ijI(20,3))+2*(Dij2345I(1,1)-3*EijI(11,2)-EijI(22,3)+EijI(23,
     -   3)-2*(EijI(37,4)-EijI(40,4)+EijI(41,4)-EijI(43,4)))+(-EijI(2
     -   ,2)-EijI(2,3)+EijI(5,2)-EijI(6,2)+EijI(8,2)+EijI(8,3)+EijI(1
     -   1,3)-EijI(12,3)+2*(EijI(9,3)-EijI(17,3)))*P(143)+(EijI(6,2)-
     -   EijI(7,2)-EijI(8,2)+EijI(9,2)-EijI(9,3)+EijI(10,3)-EijI(11,3
     -   )+EijI(12,3)+EijI(17,3)-EijI(18,3)+EijI(19,3)-EijI(20,3))*P(
     -   193))
       F(116)=DCMPLX(FR(116),FI(116))
       P(551) = -p2sq+s12+2*s15
       P(552) = p2sq+2*p5sq-s12-4*s15
       FR(117) = 2*(D02345R+Dij2345R(1,1)-6*EijR(11,2)-4*(EijR(22,3)-E
     -   ijR(23,3)-EijR(37,4)+EijR(40,4)-EijR(41,4)+EijR(43,4))-2*(Ei
     -   jR(1,1)+EijR(1,2)-EijR(2,1))*P(21)+(-EijR(7,2)+EijR(9,2))*P(
     -   186)+(-EijR(6,2)+EijR(8,2))*P(215)+EijR(2,2)*P(551)+EijR(5,2
     -   )*P(552))
       FI(117) = 2*(D02345I+Dij2345I(1,1)-6*EijI(11,2)-4*(EijI(22,3)-E
     -   ijI(23,3)-EijI(37,4)+EijI(40,4)-EijI(41,4)+EijI(43,4))-2*(Ei
     -   jI(1,1)+EijI(1,2)-EijI(2,1))*P(21)+(-EijI(7,2)+EijI(9,2))*P(
     -   186)+(-EijI(6,2)+EijI(8,2))*P(215)+EijI(2,2)*P(551)+EijI(5,2
     -   )*P(552))
       F(117)=DCMPLX(FR(117),FI(117))
       P(553) = p2sq+p5sq-s12+s15-s23-s45
       FR(118) = -2*(D02345R+Dij2345R(1,1)-Dij2345R(3,1)-4*(EijR(11,2)
     -   +EijR(22,3)-EijR(24,3)+EijR(37,4)-EijR(43,4)-EijR(44,4)+EijR
     -   (45,4))-EijR(2,2)*P(143)-EijR(4,2)*P(191)+2*((EE0R+EijR(1,1)
     -   )*P(21)+(EijR(5,2)-EijR(7,2))*P(193))+(EijR(3,1)+EijR(8,2)-E
     -   ijR(10,2))*P(215)+EijR(4,1)*P(223)-EijR(2,1)*P(224)+EijR(9,2
     -   )*P(553))
       FI(118) = -2*(D02345I+Dij2345I(1,1)-Dij2345I(3,1)-4*(EijI(11,2)
     -   +EijI(22,3)-EijI(24,3)+EijI(37,4)-EijI(43,4)-EijI(44,4)+EijI
     -   (45,4))-EijI(2,2)*P(143)-EijI(4,2)*P(191)+2*((EE0I+EijI(1,1)
     -   )*P(21)+(EijI(5,2)-EijI(7,2))*P(193))+(EijI(3,1)+EijI(8,2)-E
     -   ijI(10,2))*P(215)+EijI(4,1)*P(223)-EijI(2,1)*P(224)+EijI(9,2
     -   )*P(553))
       F(118)=DCMPLX(FR(118),FI(118))
       FR(119) = -2*(Dij2345R(2,1)+s23*(-EijR(1,1)-EijR(1,2)+EijR(2,1)
     -   +EijR(5,2)-EijR(6,2)+EijR(8,2)+2*(EijR(7,2)-EijR(9,2)))+2*s1
     -   5*(-EijR(7,2)+EijR(9,2))+s45*(EE0R+EijR(1,2)+EijR(3,2)-EijR(
     -   4,1)+2*(EijR(1,1)-EijR(2,1)-EijR(5,2))+EijR(7,2)-EijR(10,2))
     -   +p5sq*(EijR(4,1)-EijR(7,2)+EijR(10,2))+4*(EijR(11,2)+EijR(21
     -   ,3)-EijR(22,3)-EijR(37,4)+EijR(40,4)-EijR(41,4)+EijR(43,4))+
     -   (EijR(2,1)-EijR(3,1)+EijR(3,2)+EijR(5,2)-EijR(6,2)-EijR(8,2)
     -   )*P(84))
       FI(119) = -2*(Dij2345I(2,1)+s23*(-EijI(1,1)-EijI(1,2)+EijI(2,1)
     -   +EijI(5,2)-EijI(6,2)+EijI(8,2)+2*(EijI(7,2)-EijI(9,2)))+2*s1
     -   5*(-EijI(7,2)+EijI(9,2))+s45*(EE0I+EijI(1,2)+EijI(3,2)-EijI(
     -   4,1)+2*(EijI(1,1)-EijI(2,1)-EijI(5,2))+EijI(7,2)-EijI(10,2))
     -   +p5sq*(EijI(4,1)-EijI(7,2)+EijI(10,2))+4*(EijI(11,2)+EijI(21
     -   ,3)-EijI(22,3)-EijI(37,4)+EijI(40,4)-EijI(41,4)+EijI(43,4))+
     -   (EijI(2,1)-EijI(3,1)+EijI(3,2)+EijI(5,2)-EijI(6,2)-EijI(8,2)
     -   )*P(84))
       F(119)=DCMPLX(FR(119),FI(119))
       FR(120) = -4*(EijR(2,1)+EijR(3,1)-EijR(4,1)-EijR(9,2)-3*(EijR(2
     -   ,3)-EijR(9,3))+EijR(10,2)+EijR(10,3)-EijR(20,3)-2*(EijR(2,4)
     -   -EijR(8,2)-EijR(9,4)-EijR(10,4)+EijR(28,4)))
       FI(120) = -4*(EijI(2,1)+EijI(3,1)-EijI(4,1)-EijI(9,2)-3*(EijI(2
     -   ,3)-EijI(9,3))+EijI(10,2)+EijI(10,3)-EijI(20,3)-2*(EijI(2,4)
     -   -EijI(8,2)-EijI(9,4)-EijI(10,4)+EijI(28,4)))
       F(120)=DCMPLX(FR(120),FI(120))
       FR(121) = 4*(EijR(2,1)+EijR(3,1)-EijR(4,1)+EijR(9,2)+3*(EijR(2,
     -   3)-EijR(9,3))-EijR(10,2)-EijR(10,3)+EijR(20,3)+2*(EijR(2,2)+
     -   EijR(2,4)-EijR(9,4)-EijR(10,4)+EijR(28,4)))
       FI(121) = 4*(EijI(2,1)+EijI(3,1)-EijI(4,1)+EijI(9,2)+3*(EijI(2,
     -   3)-EijI(9,3))-EijI(10,2)-EijI(10,3)+EijI(20,3)+2*(EijI(2,2)+
     -   EijI(2,4)-EijI(9,4)-EijI(10,4)+EijI(28,4)))
       F(121)=DCMPLX(FR(121),FI(121))
       P(554) = p2sq+3*s12
       P(555) = p2sq-2*p3sq-s23-3*P(56)
       P(556) = p3sq-s15
       P(557) = 6*s12+s23-3*s45+2*P(556)
       FR(122) = 2*(-Dij2345R(1,2)-3*(Dij2345R(1,1)-Dij2345R(2,1))+Dij
     -   2345R(4,2)-4*(3*EijR(11,2)-2*EijR(22,3)+EijR(23,3)-EijR(37,4
     -   )+EijR(43,4))+(EijR(2,3)-2*EijR(6,2))*P(143)+(-EijR(10,3)+Ei
     -   jR(20,3))*P(186)+EijR(12,3)*P(215)+2*(D02345R+EijR(4,1)*P(52
     -   )+EijR(2,1)*P(177)-EijR(3,1)*P(296)+EijR(5,2)*P(394))+(-EijR
     -   (9,2)+EijR(10,2))*P(418)+EijR(9,3)*P(548)-EijR(2,2)*P(554)+E
     -   ijR(3,2)*P(555)+EijR(8,2)*P(557))
       FI(122) = 2*(-Dij2345I(1,2)-3*(Dij2345I(1,1)-Dij2345I(2,1))+Dij
     -   2345I(4,2)-4*(3*EijI(11,2)-2*EijI(22,3)+EijI(23,3)-EijI(37,4
     -   )+EijI(43,4))+(EijI(2,3)-2*EijI(6,2))*P(143)+(-EijI(10,3)+Ei
     -   jI(20,3))*P(186)+EijI(12,3)*P(215)+2*(D02345I+EijI(4,1)*P(52
     -   )+EijI(2,1)*P(177)-EijI(3,1)*P(296)+EijI(5,2)*P(394))+(-EijI
     -   (9,2)+EijI(10,2))*P(418)+EijI(9,3)*P(548)-EijI(2,2)*P(554)+E
     -   ijI(3,2)*P(555)+EijI(8,2)*P(557))
       F(122)=DCMPLX(FR(122),FI(122))
       P(558) = p2sq+s12+s15
       P(559) = -p2sq+p5sq+s12-s15
       P(560) = p2sq+p3sq+2*s12-s45
       FR(123) = -4*(D02345R-Dij2345R(1,1)+Dij2345R(2,1)-p5sq*EijR(2,1
     -   )+s45*EijR(3,1)+2*(-3*EijR(11,2)+EijR(23,3)-EijR(37,4)+EijR(
     -   43,4))-EijR(3,2)*P(40)+(-EijR(9,2)+EijR(10,2))*P(43)+EijR(4,
     -   1)*P(52)-EijR(6,2)*P(143)-EijR(2,2)*P(558)-EijR(5,2)*P(559)+
     -   EijR(8,2)*P(560))
       FI(123) = -4*(D02345I-Dij2345I(1,1)+Dij2345I(2,1)-p5sq*EijI(2,1
     -   )+s45*EijI(3,1)+2*(-3*EijI(11,2)+EijI(23,3)-EijI(37,4)+EijI(
     -   43,4))-EijI(3,2)*P(40)+(-EijI(9,2)+EijI(10,2))*P(43)+EijI(4,
     -   1)*P(52)-EijI(6,2)*P(143)-EijI(2,2)*P(558)-EijI(5,2)*P(559)+
     -   EijI(8,2)*P(560))
       F(123)=DCMPLX(FR(123),FI(123))
       P(561) = 3*p2sq+s12
       P(562) = -p2sq+s12+s23-3*s45
       P(563) = p2sq-s12+s23-s45
       P(564) = 2*p3sq-s23-s45+3*P(145)
       P(565) = p5sq-s15+s23-s45+2*P(69)
       FR(124) = 2*(Dij2345R(1,1)-3*(D02345R+EijR(4,1)*P(52))+EijR(2,1
     -   )*P(224)-(EE0R+EijR(1,1))*P(312)+EijR(2,2)*P(561)+EijR(3,1)*
     -   P(562)+2*(-Dij2345R(2,1)+7*EijR(11,2)+EijR(22,3)+EijR(23,3)+
     -   2*(EijR(37,4)-EijR(43,4))+EijR(3,2)*P(40)-EijR(10,2)*P(43)+E
     -   ijR(6,2)*P(143)-EijR(5,2)*P(563))-EijR(8,2)*P(564)+EijR(9,2)
     -   *P(565))
       FI(124) = 2*(Dij2345I(1,1)-3*(D02345I+EijI(4,1)*P(52))+EijI(2,1
     -   )*P(224)-(EE0I+EijI(1,1))*P(312)+EijI(2,2)*P(561)+EijI(3,1)*
     -   P(562)+2*(-Dij2345I(2,1)+7*EijI(11,2)+EijI(22,3)+EijI(23,3)+
     -   2*(EijI(37,4)-EijI(43,4))+EijI(3,2)*P(40)-EijI(10,2)*P(43)+E
     -   ijI(6,2)*P(143)-EijI(5,2)*P(563))-EijI(8,2)*P(564)+EijI(9,2)
     -   *P(565))
       F(124)=DCMPLX(FR(124),FI(124))
       FR(125) = 4*(s34*(EijR(2,2)-EijR(5,2)+EijR(7,2)-EijR(9,2))+p3sq
     -   *(EijR(6,2)-EijR(7,2)-EijR(8,2)+EijR(9,2))-2*(EijR(11,2)+2*E
     -   ijR(21,3)-EijR(22,3)-EijR(37,4)+EijR(40,4)))
       FI(125) = 4*(s34*(EijI(2,2)-EijI(5,2)+EijI(7,2)-EijI(9,2))+p3sq
     -   *(EijI(6,2)-EijI(7,2)-EijI(8,2)+EijI(9,2))-2*(EijI(11,2)+2*E
     -   ijI(21,3)-EijI(22,3)-EijI(37,4)+EijI(40,4)))
       F(125)=DCMPLX(FR(125),FI(125))
       P(566) = p2sq+s12-s34
       P(567) = p2sq+2*s12-s34-s45
       FR(126) = -4*(Dij2345R(1,1)-Dij2345R(2,1)+2*(EijR(11,2)-2*EijR(
     -   22,3)+EijR(23,3)-EijR(37,4)+EijR(43,4))+(EijR(2,1)-EijR(3,1)
     -   )*P(16)+(EijR(9,2)-EijR(10,2))*P(52)+EijR(3,2)*P(56)+(-EijR(
     -   5,2)+EijR(6,2))*P(143)+EijR(2,2)*P(566)-EijR(8,2)*P(567))
       FI(126) = -4*(Dij2345I(1,1)-Dij2345I(2,1)+2*(EijI(11,2)-2*EijI(
     -   22,3)+EijI(23,3)-EijI(37,4)+EijI(43,4))+(EijI(2,1)-EijI(3,1)
     -   )*P(16)+(EijI(9,2)-EijI(10,2))*P(52)+EijI(3,2)*P(56)+(-EijI(
     -   5,2)+EijI(6,2))*P(143)+EijI(2,2)*P(566)-EijI(8,2)*P(567))
       F(126)=DCMPLX(FR(126),FI(126))
       P(568) = p2sq-p5sq+s12-s34+s45
       FR(127) = -4*(Dij2345R(1,1)-Dij2345R(3,1)+2*(EijR(11,2)-2*EijR(
     -   22,3)+EijR(24,3)-EijR(37,4)+EijR(44,4))+(EijR(2,1)-EijR(4,1)
     -   )*P(16)-EijR(4,2)*P(52)+(-EijR(8,2)+EijR(10,2))*P(56)+(-EijR
     -   (5,2)+EijR(7,2))*P(143)+EijR(2,2)*P(566)-EijR(9,2)*P(568))
       FI(127) = -4*(Dij2345I(1,1)-Dij2345I(3,1)+2*(EijI(11,2)-2*EijI(
     -   22,3)+EijI(24,3)-EijI(37,4)+EijI(44,4))+(EijI(2,1)-EijI(4,1)
     -   )*P(16)-EijI(4,2)*P(52)+(-EijI(8,2)+EijI(10,2))*P(56)+(-EijI
     -   (5,2)+EijI(7,2))*P(143)+EijI(2,2)*P(566)-EijI(9,2)*P(568))
       F(127)=DCMPLX(FR(127),FI(127))
       FR(128) = 4*(s34*(EijR(2,1)+EijR(2,2)-EijR(4,1)-EijR(9,2))+p3sq
     -   *(-EijR(3,1)+EijR(4,1)-EijR(8,2)+EijR(9,2))+2*(EijR(11,2)+2*
     -   EijR(22,3)+EijR(37,4)))
       FI(128) = 4*(s34*(EijI(2,1)+EijI(2,2)-EijI(4,1)-EijI(9,2))+p3sq
     -   *(-EijI(3,1)+EijI(4,1)-EijI(8,2)+EijI(9,2))+2*(EijI(11,2)+2*
     -   EijI(22,3)+EijI(37,4)))
       F(128)=DCMPLX(FR(128),FI(128))
       FR(129) = -8*(EijR(2,3)-EijR(2,4)+EijR(5,3)+EijR(5,4)-2*EijR(8,
     -   3)+3*(EijR(8,4)-EijR(17,4)))
       FI(129) = -8*(EijI(2,3)-EijI(2,4)+EijI(5,3)+EijI(5,4)-2*EijI(8,
     -   3)+3*(EijI(8,4)-EijI(17,4)))
       F(129)=DCMPLX(FR(129),FI(129))
       FR(130) = 8*(EijR(2,4)+EijR(6,3)+EijR(9,3)-EijR(9,4)+EijR(17,4)
     -   -EijR(23,4)-2*(EijR(8,4)+EijR(17,3)-EijR(26,4)))
       FI(130) = 8*(EijI(2,4)+EijI(6,3)+EijI(9,3)-EijI(9,4)+EijI(17,4)
     -   -EijI(23,4)-2*(EijI(8,4)+EijI(17,3)-EijI(26,4)))
       F(130)=DCMPLX(FR(130),FI(130))
       FR(131) = 8*(EijR(2,4)+EijR(7,3)+EijR(10,3)-EijR(10,4)+EijR(17,
     -   4)-EijR(24,4)-2*(EijR(8,4)+EijR(18,3)-EijR(27,4)))
       FI(131) = 8*(EijI(2,4)+EijI(7,3)+EijI(10,3)-EijI(10,4)+EijI(17,
     -   4)-EijI(24,4)-2*(EijI(8,4)+EijI(18,3)-EijI(27,4)))
       F(131)=DCMPLX(FR(131),FI(131))
       FR(132) = -8*(EijR(2,2)+EijR(2,3)-EijR(2,4)-EijR(5,2)-EijR(8,3)
     -   -EijR(9,3)+EijR(9,4)+EijR(17,3)-EijR(17,4)+EijR(23,4)+2*(Eij
     -   R(8,4)-EijR(26,4)))
       FI(132) = -8*(EijI(2,2)+EijI(2,3)-EijI(2,4)-EijI(5,2)-EijI(8,3)
     -   -EijI(9,3)+EijI(9,4)+EijI(17,3)-EijI(17,4)+EijI(23,4)+2*(Eij
     -   I(8,4)-EijI(26,4)))
       F(132)=DCMPLX(FR(132),FI(132))
       FR(133) = -8*(EijR(2,2)-EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)
     -   +EijR(8,4)-EijR(9,3)-EijR(11,3)+EijR(12,3)+EijR(17,3)-EijR(2
     -   0,4)+2*(EijR(9,4)-EijR(26,4))+EijR(29,4))
       FI(133) = -8*(EijI(2,2)-EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)
     -   +EijI(8,4)-EijI(9,3)-EijI(11,3)+EijI(12,3)+EijI(17,3)-EijI(2
     -   0,4)+2*(EijI(9,4)-EijI(26,4))+EijI(29,4))
       F(133)=DCMPLX(FR(133),FI(133))
       FR(134) = -8*(EijR(2,2)-EijR(2,4)-EijR(5,2)+EijR(8,4)+2*(EijR(7
     -   ,2)-EijR(9,2))+EijR(9,4)-EijR(10,3)+EijR(10,4)+EijR(18,3)-Ei
     -   jR(19,3)+EijR(20,3)-EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(35
     -   ,4))
       FI(134) = -8*(EijI(2,2)-EijI(2,4)-EijI(5,2)+EijI(8,4)+2*(EijI(7
     -   ,2)-EijI(9,2))+EijI(9,4)-EijI(10,3)+EijI(10,4)+EijI(18,3)-Ei
     -   jI(19,3)+EijI(20,3)-EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(35
     -   ,4))
       F(134)=DCMPLX(FR(134),FI(134))
       FR(135) = -8*(EijR(2,2)-EijR(2,4)-EijR(5,2)+EijR(8,4)+EijR(9,4)
     -   -EijR(26,4))
       FI(135) = -8*(EijI(2,2)-EijI(2,4)-EijI(5,2)+EijI(8,4)+EijI(9,4)
     -   -EijI(26,4))
       F(135)=DCMPLX(FR(135),FI(135))
       FR(136) = -8*(EijR(2,2)-EijR(2,4)-EijR(5,2)-EijR(5,3)+EijR(8,3)
     -   -EijR(10,3)+EijR(10,4)-EijR(17,4)+EijR(18,3)+EijR(24,4)+2*(E
     -   ijR(8,4)-EijR(27,4)))
       FI(136) = -8*(EijI(2,2)-EijI(2,4)-EijI(5,2)-EijI(5,3)+EijI(8,3)
     -   -EijI(10,3)+EijI(10,4)-EijI(17,4)+EijI(18,3)+EijI(24,4)+2*(E
     -   ijI(8,4)-EijI(27,4)))
       F(136)=DCMPLX(FR(136),FI(136))
       FR(137) = 8*(EijR(2,3)+EijR(2,4)-EijR(6,2)+EijR(7,2)+EijR(8,2)-
     -   EijR(8,3)-EijR(8,4)-EijR(9,2)-EijR(9,4)-EijR(10,4)+EijR(19,3
     -   )-EijR(20,3)+EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4))
       FI(137) = 8*(EijI(2,3)+EijI(2,4)-EijI(6,2)+EijI(7,2)+EijI(8,2)-
     -   EijI(8,3)-EijI(8,4)-EijI(9,2)-EijI(9,4)-EijI(10,4)+EijI(19,3
     -   )-EijI(20,3)+EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4))
       F(137)=DCMPLX(FR(137),FI(137))
       FR(138) = 8*(EijR(2,3)+EijR(2,4)-EijR(7,2)-EijR(8,3)-EijR(8,4)+
     -   EijR(9,2)+EijR(14,3)-EijR(15,3)+EijR(21,4)-2*(EijR(10,4)-Eij
     -   R(27,4))-EijR(32,4))
       FI(138) = 8*(EijI(2,3)+EijI(2,4)-EijI(7,2)-EijI(8,3)-EijI(8,4)+
     -   EijI(9,2)+EijI(14,3)-EijI(15,3)+EijI(21,4)-2*(EijI(10,4)-Eij
     -   I(27,4))-EijI(32,4))
       F(138)=DCMPLX(FR(138),FI(138))
       FR(139) = 8*(EijR(2,3)+EijR(2,4)-EijR(8,3)-EijR(8,4)-EijR(10,4)
     -   +EijR(27,4))
       FI(139) = 8*(EijI(2,3)+EijI(2,4)-EijI(8,3)-EijI(8,4)-EijI(10,4)
     -   +EijI(27,4))
       F(139)=DCMPLX(FR(139),FI(139))
       FR(140) = -8*(EijR(2,2)-EijR(2,4)-EijR(5,2)-EijR(5,3)+EijR(8,3)
     -   +2*EijR(8,4)-EijR(17,4))
       FI(140) = -8*(EijI(2,2)-EijI(2,4)-EijI(5,2)-EijI(5,3)+EijI(8,3)
     -   +2*EijI(8,4)-EijI(17,4))
       F(140)=DCMPLX(FR(140),FI(140))
       FR(141) = 8*(EijR(2,3)+EijR(2,4)-EijR(6,2)+EijR(8,2)-EijR(8,3)-
     -   EijR(8,4)-EijR(9,4)+EijR(26,4))
       FI(141) = 8*(EijI(2,3)+EijI(2,4)-EijI(6,2)+EijI(8,2)-EijI(8,3)-
     -   EijI(8,4)-EijI(9,4)+EijI(26,4))
       F(141)=DCMPLX(FR(141),FI(141))
       FR(142) = 8*(EijR(2,3)+EijR(2,4)-EijR(7,2)-EijR(8,3)-EijR(8,4)+
     -   EijR(9,2)-EijR(10,4)+EijR(27,4))
       FI(142) = 8*(EijI(2,3)+EijI(2,4)-EijI(7,2)-EijI(8,3)-EijI(8,4)+
     -   EijI(9,2)-EijI(10,4)+EijI(27,4))
       F(142)=DCMPLX(FR(142),FI(142))
       FR(143) = 8*(EijR(2,3)+EijR(2,4)-EijR(8,3)-EijR(8,4))
       FI(143) = 8*(EijI(2,3)+EijI(2,4)-EijI(8,3)-EijI(8,4))
       F(143)=DCMPLX(FR(143),FI(143))
       FR(144) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)-EijR(5,2)+EijR(6,2)
     -   -EijR(8,2)+EijR(9,3)+EijR(9,4)+3*(EijR(8,3)-EijR(17,3))-EijR
     -   (17,4)+EijR(23,4)-2*(EijR(5,3)-EijR(6,3)-EijR(8,4)+EijR(26,4
     -   )))
       FI(144) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)-EijI(5,2)+EijI(6,2)
     -   -EijI(8,2)+EijI(9,3)+EijI(9,4)+3*(EijI(8,3)-EijI(17,3))-EijI
     -   (17,4)+EijI(23,4)-2*(EijI(5,3)-EijI(6,3)-EijI(8,4)+EijI(26,4
     -   )))
       F(144)=DCMPLX(FR(144),FI(144))
       FR(145) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(7,2)-EijR(8,4)-
     -   EijR(9,2)-EijR(9,4)-EijR(10,3)-EijR(10,4)+2*(EijR(2,3)-EijR(
     -   8,3)-EijR(9,3)+EijR(17,3))+EijR(18,3)-EijR(19,3)+EijR(20,3)+
     -   EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4))
       FI(145) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(7,2)-EijI(8,4)-
     -   EijI(9,2)-EijI(9,4)-EijI(10,3)-EijI(10,4)+2*(EijI(2,3)-EijI(
     -   8,3)-EijI(9,3)+EijI(17,3))+EijI(18,3)-EijI(19,3)+EijI(20,3)+
     -   EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4))
       F(145)=DCMPLX(FR(145),FI(145))
       FR(146) = -8*(EijR(2,1)+EijR(2,2)-EijR(2,3)-EijR(2,4)-EijR(3,1)
     -   +EijR(3,2)+EijR(5,2)-EijR(6,2)+EijR(8,4)-EijR(12,3)-4*EijR(1
     -   7,3)-EijR(20,4)-2*(EijR(8,2)-EijR(8,3)-EijR(9,3)-EijR(9,4)-E
     -   ijR(11,3)+EijR(26,4))+EijR(29,4))
       FI(146) = -8*(EijI(2,1)+EijI(2,2)-EijI(2,3)-EijI(2,4)-EijI(3,1)
     -   +EijI(3,2)+EijI(5,2)-EijI(6,2)+EijI(8,4)-EijI(12,3)-4*EijI(1
     -   7,3)-EijI(20,4)-2*(EijI(8,2)-EijI(8,3)-EijI(9,3)-EijI(9,4)-E
     -   ijI(11,3)+EijI(26,4))+EijI(29,4))
       F(146)=DCMPLX(FR(146),FI(146))
       FR(147) = 8*(EijR(2,2)+EijR(2,4)-EijR(8,2)-4*EijR(9,3)-EijR(10,
     -   3)-EijR(10,4)-EijR(13,3)+EijR(20,4)+2*(EijR(2,3)-EijR(9,4)+E
     -   ijR(12,3)+EijR(20,3)+EijR(28,4))-EijR(31,4))
       FI(147) = 8*(EijI(2,2)+EijI(2,4)-EijI(8,2)-4*EijI(9,3)-EijI(10,
     -   3)-EijI(10,4)-EijI(13,3)+EijI(20,4)+2*(EijI(2,3)-EijI(9,4)+E
     -   ijI(12,3)+EijI(20,3)+EijI(28,4))-EijI(31,4))
       F(147)=DCMPLX(FR(147),FI(147))
       FR(148) = -8*(EijR(2,1)+EijR(2,2)-EijR(2,4)-EijR(3,1)+EijR(5,2)
     -   +EijR(7,2)+EijR(8,4)+EijR(9,4)+EijR(10,2)+EijR(10,3)+EijR(10
     -   ,4)+3*(EijR(8,3)-EijR(17,3))-2*(EijR(2,3)+EijR(6,2)+EijR(9,2
     -   )-EijR(9,3)+EijR(18,3)-EijR(19,3))-EijR(20,3)-EijR(26,4)-Eij
     -   R(27,4)-EijR(28,4)+EijR(35,4))
       FI(148) = -8*(EijI(2,1)+EijI(2,2)-EijI(2,4)-EijI(3,1)+EijI(5,2)
     -   +EijI(7,2)+EijI(8,4)+EijI(9,4)+EijI(10,2)+EijI(10,3)+EijI(10
     -   ,4)+3*(EijI(8,3)-EijI(17,3))-2*(EijI(2,3)+EijI(6,2)+EijI(9,2
     -   )-EijI(9,3)+EijI(18,3)-EijI(19,3))-EijI(20,3)-EijI(26,4)-Eij
     -   I(27,4)-EijI(28,4)+EijI(35,4))
       F(148)=DCMPLX(FR(148),FI(148))
       FR(149) = 8*(EijR(2,4)+EijR(3,2)-EijR(9,2)-5*EijR(9,3)+EijR(10,
     -   2)-EijR(10,4)-EijR(13,3)+3*(EijR(2,3)-EijR(8,2)+EijR(20,3))+
     -   EijR(20,4)+2*(EijR(2,2)-EijR(9,4)-EijR(10,3)+EijR(12,3)+EijR
     -   (28,4))-EijR(31,4))
       FI(149) = 8*(EijI(2,4)+EijI(3,2)-EijI(9,2)-5*EijI(9,3)+EijI(10,
     -   2)-EijI(10,4)-EijI(13,3)+3*(EijI(2,3)-EijI(8,2)+EijI(20,3))+
     -   EijI(20,4)+2*(EijI(2,2)-EijI(9,4)-EijI(10,3)+EijI(12,3)+EijI
     -   (28,4))-EijI(31,4))
       F(149)=DCMPLX(FR(149),FI(149))
       FR(150) = 8*(EijR(2,4)-EijR(9,2)+3*(EijR(2,3)-EijR(9,3))-EijR(9
     -   ,4)+EijR(10,2)+EijR(15,3)-EijR(16,3)-4*(EijR(10,3)-EijR(20,3
     -   ))+EijR(21,4)+2*(EijR(2,2)-EijR(8,2)-EijR(10,4)+EijR(28,4))-
     -   EijR(34,4))
       FI(150) = 8*(EijI(2,4)-EijI(9,2)+3*(EijI(2,3)-EijI(9,3))-EijI(9
     -   ,4)+EijI(10,2)+EijI(15,3)-EijI(16,3)-4*(EijI(10,3)-EijI(20,3
     -   ))+EijI(21,4)+2*(EijI(2,2)-EijI(8,2)-EijI(10,4)+EijI(28,4))-
     -   EijI(34,4))
       F(150)=DCMPLX(FR(150),FI(150))
       FR(151) = 8*(EijR(2,4)-EijR(3,1)+EijR(4,1)+3*(EijR(2,3)-EijR(8,
     -   2)-EijR(9,3))-EijR(9,4)+EijR(10,2)-EijR(10,4)+2*(EijR(2,2)-E
     -   ijR(10,3)+EijR(20,3))+EijR(28,4))
       FI(151) = 8*(EijI(2,4)-EijI(3,1)+EijI(4,1)+3*(EijI(2,3)-EijI(8,
     -   2)-EijI(9,3))-EijI(9,4)+EijI(10,2)-EijI(10,4)+2*(EijI(2,2)-E
     -   ijI(10,3)+EijI(20,3))+EijI(28,4))
       F(151)=DCMPLX(FR(151),FI(151))
       FR(152) = -8*(EijR(2,1)-EijR(2,4)-EijR(3,1)+EijR(8,4)-2*(EijR(2
     -   ,3)-EijR(5,2)+EijR(6,2)-EijR(9,3))+EijR(9,4)+3*(EijR(8,3)-Ei
     -   jR(17,3))-EijR(26,4))
       FI(152) = -8*(EijI(2,1)-EijI(2,4)-EijI(3,1)+EijI(8,4)-2*(EijI(2
     -   ,3)-EijI(5,2)+EijI(6,2)-EijI(9,3))+EijI(9,4)+3*(EijI(8,3)-Ei
     -   jI(17,3))-EijI(26,4))
       F(152)=DCMPLX(FR(152),FI(152))
       FR(153) = 8*(EijR(2,1)+3*(EijR(2,2)+EijR(2,3))+EijR(2,4)-EijR(3
     -   ,1)+EijR(3,2)-4*EijR(8,2)-5*EijR(9,3)-2*(EijR(9,4)-EijR(12,3
     -   ))+EijR(20,4))
       FI(153) = 8*(EijI(2,1)+3*(EijI(2,2)+EijI(2,3))+EijI(2,4)-EijI(3
     -   ,1)+EijI(3,2)-4*EijI(8,2)-5*EijI(9,3)-2*(EijI(9,4)-EijI(12,3
     -   ))+EijI(20,4))
       F(153)=DCMPLX(FR(153),FI(153))
       FR(154) = 8*(EijR(2,1)+EijR(2,4)-EijR(4,1)+3*(EijR(2,2)+EijR(2,
     -   3)-EijR(9,3))-EijR(9,4)+EijR(10,2)-EijR(10,4)-2*(EijR(8,2)+E
     -   ijR(9,2)+EijR(10,3)-EijR(20,3))+EijR(28,4))
       FI(154) = 8*(EijI(2,1)+EijI(2,4)-EijI(4,1)+3*(EijI(2,2)+EijI(2,
     -   3)-EijI(9,3))-EijI(9,4)+EijI(10,2)-EijI(10,4)-2*(EijI(8,2)+E
     -   ijI(9,2)+EijI(10,3)-EijI(20,3))+EijI(28,4))
       F(154)=DCMPLX(FR(154),FI(154))
       FR(155) = 8*(EijR(2,1)+EijR(2,4)-EijR(3,1)+3*(EijR(2,2)+EijR(2,
     -   3)-EijR(8,2)-EijR(9,3))-EijR(9,4))
       FI(155) = 8*(EijI(2,1)+EijI(2,4)-EijI(3,1)+3*(EijI(2,2)+EijI(2,
     -   3)-EijI(8,2)-EijI(9,3))-EijI(9,4))
       F(155)=DCMPLX(FR(155),FI(155))
       P(569) = p2sq-2*s12-s23+s45
       P(570) = p2sq+p3sq-2*s12-s23+s45
       P(571) = p2sq-p5sq+s12+s45
       FR(156) = 4*(-Dij2345R(1,1)+Dij2345R(3,1)+p3sq*EijR(3,1)-2*(Eij
     -   R(11,2)-2*(EijR(22,3)-EijR(24,3))-EijR(37,4)+EijR(44,4))+(Ei
     -   jR(8,2)-EijR(10,2))*P(40)+EijR(4,2)*P(52)-EijR(2,2)*P(145)+(
     -   EijR(5,2)-EijR(7,2))*P(338)+EijR(2,1)*P(569)-EijR(4,1)*P(570
     -   )+EijR(9,2)*P(571))
       FI(156) = 4*(-Dij2345I(1,1)+Dij2345I(3,1)+p3sq*EijI(3,1)-2*(Eij
     -   I(11,2)-2*(EijI(22,3)-EijI(24,3))-EijI(37,4)+EijI(44,4))+(Ei
     -   jI(8,2)-EijI(10,2))*P(40)+EijI(4,2)*P(52)-EijI(2,2)*P(145)+(
     -   EijI(5,2)-EijI(7,2))*P(338)+EijI(2,1)*P(569)-EijI(4,1)*P(570
     -   )+EijI(9,2)*P(571))
       F(156)=DCMPLX(FR(156),FI(156))
       FR(157) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)-EijR(5,2)+EijR(7,2)
     -   -EijR(9,2)+EijR(10,3)+EijR(10,4)-EijR(17,4)+3*(EijR(8,3)-Eij
     -   R(18,3))+EijR(24,4)-2*(EijR(5,3)-EijR(7,3)-EijR(8,4)+EijR(27
     -   ,4)))
       FI(157) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)-EijI(5,2)+EijI(7,2)
     -   -EijI(9,2)+EijI(10,3)+EijI(10,4)-EijI(17,4)+3*(EijI(8,3)-Eij
     -   I(18,3))+EijI(24,4)-2*(EijI(5,3)-EijI(7,3)-EijI(8,4)+EijI(27
     -   ,4)))
       F(157)=DCMPLX(FR(157),FI(157))
       FR(158) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-
     -   EijR(8,4)-EijR(9,3)-EijR(9,4)-EijR(10,4)+EijR(17,3)+2*(EijR(
     -   2,3)-EijR(8,3)-EijR(10,3)+EijR(18,3))-EijR(19,3)+EijR(20,3)+
     -   EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4))
       FI(158) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-
     -   EijI(8,4)-EijI(9,3)-EijI(9,4)-EijI(10,4)+EijI(17,3)+2*(EijI(
     -   2,3)-EijI(8,3)-EijI(10,3)+EijI(18,3))-EijI(19,3)+EijI(20,3)+
     -   EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4))
       F(158)=DCMPLX(FR(158),FI(158))
       FR(159) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(7,2)-EijR(8,4)-
     -   EijR(9,2)-EijR(14,3)+EijR(15,3)-3*(EijR(10,3)-EijR(18,3))+Ei
     -   jR(21,4)+2*(EijR(2,3)-EijR(8,3)-EijR(10,4)+EijR(27,4))-EijR(
     -   32,4))
       FI(159) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(7,2)-EijI(8,4)-
     -   EijI(9,2)-EijI(14,3)+EijI(15,3)-3*(EijI(10,3)-EijI(18,3))+Ei
     -   jI(21,4)+2*(EijI(2,3)-EijI(8,3)-EijI(10,4)+EijI(27,4))-EijI(
     -   32,4))
       F(159)=DCMPLX(FR(159),FI(159))
       FR(160) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(7,2)-EijR(8,4)-
     -   EijR(9,2)-EijR(10,4)+2*(EijR(2,3)-EijR(8,3)-EijR(10,3)+EijR(
     -   18,3))+EijR(27,4))
       FI(160) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(7,2)-EijI(8,4)-
     -   EijI(9,2)-EijI(10,4)+2*(EijI(2,3)-EijI(8,3)-EijI(10,3)+EijI(
     -   18,3))+EijI(27,4))
       F(160)=DCMPLX(FR(160),FI(160))
       FR(161) = -8*(EijR(2,1)-EijR(2,3)-EijR(2,4)-EijR(4,1)-EijR(8,2)
     -   +EijR(8,4)+EijR(9,3)+EijR(9,4)+EijR(10,2)+EijR(10,3)+EijR(10
     -   ,4)+2*(EijR(5,2)-EijR(7,2)+EijR(8,3)-EijR(17,3)-EijR(18,3)+E
     -   ijR(19,3))-EijR(20,3)-EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(
     -   35,4))
       FI(161) = -8*(EijI(2,1)-EijI(2,3)-EijI(2,4)-EijI(4,1)-EijI(8,2)
     -   +EijI(8,4)+EijI(9,3)+EijI(9,4)+EijI(10,2)+EijI(10,3)+EijI(10
     -   ,4)+2*(EijI(5,2)-EijI(7,2)+EijI(8,3)-EijI(17,3)-EijI(18,3)+E
     -   ijI(19,3))-EijI(20,3)-EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(
     -   35,4))
       F(161)=DCMPLX(FR(161),FI(161))
       FR(162) = 8*(EijR(2,1)+EijR(2,4)-EijR(3,1)+EijR(3,2)-EijR(9,2)+
     -   EijR(10,2)-EijR(10,4)+EijR(12,3)-EijR(13,3)-3*(EijR(8,2)+Eij
     -   R(9,3)-EijR(20,3))+EijR(20,4)+2*(EijR(2,2)+EijR(2,3)-EijR(9,
     -   4)-EijR(10,3)+EijR(28,4))-EijR(31,4))
       FI(162) = 8*(EijI(2,1)+EijI(2,4)-EijI(3,1)+EijI(3,2)-EijI(9,2)+
     -   EijI(10,2)-EijI(10,4)+EijI(12,3)-EijI(13,3)-3*(EijI(8,2)+Eij
     -   I(9,3)-EijI(20,3))+EijI(20,4)+2*(EijI(2,2)+EijI(2,3)-EijI(9,
     -   4)-EijI(10,3)+EijI(28,4))-EijI(31,4))
       F(162)=DCMPLX(FR(162),FI(162))
       FR(163) = 8*(EijR(2,1)+EijR(2,4)-EijR(4,1)-EijR(8,2)-EijR(9,4)+
     -   EijR(10,2)+EijR(15,3)-EijR(16,3)-3*(EijR(10,3)-EijR(20,3))+E
     -   ijR(21,4)+2*(EijR(2,2)+EijR(2,3)-EijR(9,2)-EijR(9,3)-EijR(10
     -   ,4)+EijR(28,4))-EijR(34,4))
       FI(163) = 8*(EijI(2,1)+EijI(2,4)-EijI(4,1)-EijI(8,2)-EijI(9,4)+
     -   EijI(10,2)+EijI(15,3)-EijI(16,3)-3*(EijI(10,3)-EijI(20,3))+E
     -   ijI(21,4)+2*(EijI(2,2)+EijI(2,3)-EijI(9,2)-EijI(9,3)-EijI(10
     -   ,4)+EijI(28,4))-EijI(34,4))
       F(163)=DCMPLX(FR(163),FI(163))
       FR(164) = 8*(EijR(2,1)+EijR(2,4)-EijR(4,1)-EijR(8,2)-EijR(9,4)+
     -   EijR(10,2)-EijR(10,4)+2*(EijR(2,2)+EijR(2,3)-EijR(9,2)-EijR(
     -   9,3)-EijR(10,3)+EijR(20,3))+EijR(28,4))
       FI(164) = 8*(EijI(2,1)+EijI(2,4)-EijI(4,1)-EijI(8,2)-EijI(9,4)+
     -   EijI(10,2)-EijI(10,4)+2*(EijI(2,2)+EijI(2,3)-EijI(9,2)-EijI(
     -   9,3)-EijI(10,3)+EijI(20,3))+EijI(28,4))
       F(164)=DCMPLX(FR(164),FI(164))
       FR(165) = -8*(EijR(2,1)-EijR(2,4)-EijR(4,1)+EijR(4,2)+EijR(8,4)
     -   -EijR(9,2)+3*(EijR(8,3)+EijR(10,3))-EijR(15,3)-5*EijR(18,3)-
     -   EijR(21,4)-2*(EijR(2,3)-EijR(5,2)+EijR(7,2)-EijR(10,4)-EijR(
     -   14,3)+EijR(27,4))+EijR(32,4))
       FI(165) = -8*(EijI(2,1)-EijI(2,4)-EijI(4,1)+EijI(4,2)+EijI(8,4)
     -   -EijI(9,2)+3*(EijI(8,3)+EijI(10,3))-EijI(15,3)-5*EijI(18,3)-
     -   EijI(21,4)-2*(EijI(2,3)-EijI(5,2)+EijI(7,2)-EijI(10,4)-EijI(
     -   14,3)+EijI(27,4))+EijI(32,4))
       F(165)=DCMPLX(FR(165),FI(165))
       FR(166) = 8*(EijR(2,1)+EijR(2,4)-EijR(3,1)+EijR(4,2)-4*EijR(9,2
     -   )-EijR(9,4)-5*EijR(10,3)-EijR(16,3)+3*(EijR(2,2)+EijR(2,3)+E
     -   ijR(20,3))+EijR(21,4)-2*(EijR(8,2)+EijR(9,3)-EijR(10,2)+EijR
     -   (10,4)-EijR(15,3)-EijR(28,4))-EijR(34,4))
       FI(166) = 8*(EijI(2,1)+EijI(2,4)-EijI(3,1)+EijI(4,2)-4*EijI(9,2
     -   )-EijI(9,4)-5*EijI(10,3)-EijI(16,3)+3*(EijI(2,2)+EijI(2,3)+E
     -   ijI(20,3))+EijI(21,4)-2*(EijI(8,2)+EijI(9,3)-EijI(10,2)+EijI
     -   (10,4)-EijI(15,3)-EijI(28,4))-EijI(34,4))
       F(166)=DCMPLX(FR(166),FI(166))
       FR(167) = 8*(EijR(2,1)+EijR(2,4)-EijR(4,1)+2*EijR(4,2)-EijR(4,3
     -   )-7*EijR(10,3)-5*(EijR(9,2)-EijR(15,3))-EijR(15,4)+3*(EijR(2
     -   ,2)+EijR(2,3)-EijR(10,4)+EijR(21,4)))
       FI(167) = 8*(EijI(2,1)+EijI(2,4)-EijI(4,1)+2*EijI(4,2)-EijI(4,3
     -   )-7*EijI(10,3)-5*(EijI(9,2)-EijI(15,3))-EijI(15,4)+3*(EijI(2
     -   ,2)+EijI(2,3)-EijI(10,4)+EijI(21,4)))
       F(167)=DCMPLX(FR(167),FI(167))
       FR(168) = 8*(EijR(2,1)+3*(EijR(2,2)+EijR(2,3))+EijR(2,4)-EijR(4
     -   ,1)+EijR(4,2)-4*EijR(9,2)-5*EijR(10,3)-2*(EijR(10,4)-EijR(15
     -   ,3))+EijR(21,4))
       FI(168) = 8*(EijI(2,1)+3*(EijI(2,2)+EijI(2,3))+EijI(2,4)-EijI(4
     -   ,1)+EijI(4,2)-4*EijI(9,2)-5*EijI(10,3)-2*(EijI(10,4)-EijI(15
     -   ,3))+EijI(21,4))
       F(168)=DCMPLX(FR(168),FI(168))
       FR(169) = -8*(EijR(2,1)-EijR(2,4)-EijR(4,1)+EijR(8,4)-2*(EijR(2
     -   ,3)-EijR(5,2)+EijR(7,2)-EijR(10,3))+EijR(10,4)+3*(EijR(8,3)-
     -   EijR(18,3))-EijR(27,4))
       FI(169) = -8*(EijI(2,1)-EijI(2,4)-EijI(4,1)+EijI(8,4)-2*(EijI(2
     -   ,3)-EijI(5,2)+EijI(7,2)-EijI(10,3))+EijI(10,4)+3*(EijI(8,3)-
     -   EijI(18,3))-EijI(27,4))
       F(169)=DCMPLX(FR(169),FI(169))
       FR(170) = 8*(EijR(2,1)+EijR(2,4)-EijR(3,1)-EijR(9,4)+EijR(10,2)
     -   +3*(EijR(2,2)+EijR(2,3)-EijR(10,3))-EijR(10,4)-2*(EijR(8,2)+
     -   EijR(9,2)+EijR(9,3)-EijR(20,3))+EijR(28,4))
       FI(170) = 8*(EijI(2,1)+EijI(2,4)-EijI(3,1)-EijI(9,4)+EijI(10,2)
     -   +3*(EijI(2,2)+EijI(2,3)-EijI(10,3))-EijI(10,4)-2*(EijI(8,2)+
     -   EijI(9,2)+EijI(9,3)-EijI(20,3))+EijI(28,4))
       F(170)=DCMPLX(FR(170),FI(170))
       FR(171) = 8*(EijR(2,1)+EijR(2,4)-EijR(4,1)+3*(EijR(2,2)+EijR(2,
     -   3)-EijR(9,2)-EijR(10,3))-EijR(10,4))
       FI(171) = 8*(EijI(2,1)+EijI(2,4)-EijI(4,1)+3*(EijI(2,2)+EijI(2,
     -   3)-EijI(9,2)-EijI(10,3))-EijI(10,4))
       F(171)=DCMPLX(FR(171),FI(171))
       FR(172) = 4*(EijR(1,1)+EijR(1,2)-EijR(2,1)+EijR(2,3)+EijR(5,2)+
     -   3*EijR(5,3)-4*(EijR(8,3)+EijR(8,4))-2*(EijR(2,2)-EijR(2,4)-E
     -   ijR(17,4)))
       FI(172) = 4*(EijI(1,1)+EijI(1,2)-EijI(2,1)+EijI(2,3)+EijI(5,2)+
     -   3*EijI(5,3)-4*(EijI(8,3)+EijI(8,4))-2*(EijI(2,2)-EijI(2,4)-E
     -   ijI(17,4)))
       F(172)=DCMPLX(FR(172),FI(172))
       FR(173) = 4*(EijR(2,2)-EijR(5,2)+3*(EijR(2,3)-EijR(8,3))-EijR(9
     -   ,3)+EijR(17,3)+2*(EijR(2,4)-EijR(8,4)-EijR(9,4)+EijR(26,4)))
       FI(173) = 4*(EijI(2,2)-EijI(5,2)+3*(EijI(2,3)-EijI(8,3))-EijI(9
     -   ,3)+EijI(17,3)+2*(EijI(2,4)-EijI(8,4)-EijI(9,4)+EijI(26,4)))
       F(173)=DCMPLX(FR(173),FI(173))
       FR(174) = 4*(EijR(2,2)-EijR(5,2)+3*(EijR(2,3)-EijR(8,3))-EijR(1
     -   0,3)+EijR(18,3)+2*(EijR(2,4)-EijR(8,4)-EijR(10,4)+EijR(27,4)
     -   ))
       FI(174) = 4*(EijI(2,2)-EijI(5,2)+3*(EijI(2,3)-EijI(8,3))-EijI(1
     -   0,3)+EijI(18,3)+2*(EijI(2,4)-EijI(8,4)-EijI(10,4)+EijI(27,4)
     -   ))
       F(174)=DCMPLX(FR(174),FI(174))
       FR(175) = 4*(EijR(2,2)-EijR(5,2)+3*(EijR(2,3)-EijR(8,3))+2*(Eij
     -   R(2,4)-EijR(8,4)))
       FI(175) = 4*(EijI(2,2)-EijI(5,2)+3*(EijI(2,3)-EijI(8,3))+2*(Eij
     -   I(2,4)-EijI(8,4)))
       F(175)=DCMPLX(FR(175),FI(175))
       FR(176) = 4*(EijR(2,2)+3*EijR(2,3)-EijR(8,2)-4*(EijR(9,3)+EijR(
     -   9,4))+EijR(12,3)+2*(EijR(2,4)+EijR(20,4)))
       FI(176) = 4*(EijI(2,2)+3*EijI(2,3)-EijI(8,2)-4*(EijI(9,3)+EijI(
     -   9,4))+EijI(12,3)+2*(EijI(2,4)+EijI(20,4)))
       F(176)=DCMPLX(FR(176),FI(176))
       FR(177) = 4*(EijR(2,2)-EijR(8,2)+3*(EijR(2,3)-EijR(9,3))+2*(Eij
     -   R(2,4)-EijR(9,4)))
       FI(177) = 4*(EijI(2,2)-EijI(8,2)+3*(EijI(2,3)-EijI(9,3))+2*(Eij
     -   I(2,4)-EijI(9,4)))
       F(177)=DCMPLX(FR(177),FI(177))
       FR(178) = 4*(EijR(2,1)+4*EijR(2,2)+5*EijR(2,3)-EijR(4,1)-EijR(8
     -   ,2)-3*(EijR(9,2)+EijR(9,3)+EijR(10,3))+EijR(20,3)+2*(EijR(2,
     -   4)-EijR(9,4)-EijR(10,4)+EijR(28,4)))
       FI(178) = 4*(EijI(2,1)+4*EijI(2,2)+5*EijI(2,3)-EijI(4,1)-EijI(8
     -   ,2)-3*(EijI(9,2)+EijI(9,3)+EijI(10,3))+EijI(20,3)+2*(EijI(2,
     -   4)-EijI(9,4)-EijI(10,4)+EijI(28,4)))
       F(178)=DCMPLX(FR(178),FI(178))
       P(572) = p2sq-3*s12+2*P(21)
       P(573) = p2sq-s23-3*P(56)-2*P(69)
       FR(179) = -2*(3*D02345R+Dij2345R(1,2)+8*EijR(21,3)-10*EijR(22,3
     -   )+4*(Dij2345R(1,1)-EijR(37,4)+EijR(40,4))-2*(EijR(11,2)+(Eij
     -   R(1,1)+EijR(1,2)-EijR(2,1))*P(21))+(-EijR(2,3)+EijR(8,3))*P(
     -   143)-EijR(2,2)*P(175)+(EijR(10,3)-EijR(18,3))*P(186)+(EijR(9
     -   ,3)-EijR(17,3))*P(215)+(-EijR(7,2)+EijR(9,2))*P(418)+EijR(5,
     -   2)*P(572)+(-EijR(6,2)+EijR(8,2))*P(573))
       FI(179) = -2*(3*D02345I+Dij2345I(1,2)+8*EijI(21,3)-10*EijI(22,3
     -   )+4*(Dij2345I(1,1)-EijI(37,4)+EijI(40,4))-2*(EijI(11,2)+(Eij
     -   I(1,1)+EijI(1,2)-EijI(2,1))*P(21))+(-EijI(2,3)+EijI(8,3))*P(
     -   143)-EijI(2,2)*P(175)+(EijI(10,3)-EijI(18,3))*P(186)+(EijI(9
     -   ,3)-EijI(17,3))*P(215)+(-EijI(7,2)+EijI(9,2))*P(418)+EijI(5,
     -   2)*P(572)+(-EijI(6,2)+EijI(8,2))*P(573))
       F(179)=DCMPLX(FR(179),FI(179))
       FR(180) = -4*(EE0R+EijR(1,1)-EijR(2,3)-EijR(3,1)-EijR(6,2)+EijR
     -   (9,3)+3*(EijR(2,1)+EijR(5,2)+EijR(8,3)-EijR(17,3))+2*(EijR(2
     -   ,2)-EijR(2,4)-EijR(8,2)+EijR(8,4)+EijR(9,4)-EijR(26,4)))
       FI(180) = -4*(EE0I+EijI(1,1)-EijI(2,3)-EijI(3,1)-EijI(6,2)+EijI
     -   (9,3)+3*(EijI(2,1)+EijI(5,2)+EijI(8,3)-EijI(17,3))+2*(EijI(2
     -   ,2)-EijI(2,4)-EijI(8,2)+EijI(8,4)+EijI(9,4)-EijI(26,4)))
       F(180)=DCMPLX(FR(180),FI(180))
       FR(181) = -4*(EE0R+EijR(1,1)+EijR(2,2)-EijR(4,1)+4*EijR(5,2)-Ei
     -   jR(7,2)+5*EijR(8,3)+EijR(10,3)+3*(EijR(2,1)-EijR(2,3)-EijR(1
     -   8,3))-2*(EijR(2,4)-EijR(8,4)+EijR(9,2)-EijR(10,4)+EijR(27,4)
     -   ))
       FI(181) = -4*(EE0I+EijI(1,1)+EijI(2,2)-EijI(4,1)+4*EijI(5,2)-Ei
     -   jI(7,2)+5*EijI(8,3)+EijI(10,3)+3*(EijI(2,1)-EijI(2,3)-EijI(1
     -   8,3))-2*(EijI(2,4)-EijI(8,4)+EijI(9,2)-EijI(10,4)+EijI(27,4)
     -   ))
       F(181)=DCMPLX(FR(181),FI(181))
       FR(182) = 4*(EijR(2,2)+EijR(4,1)-EijR(8,2)+3*(EijR(2,3)-EijR(9,
     -   3))-EijR(10,3)+EijR(20,3)+2*(EijR(2,4)+EijR(9,2)-EijR(9,4)-E
     -   ijR(10,4)+EijR(28,4)))
       FI(182) = 4*(EijI(2,2)+EijI(4,1)-EijI(8,2)+3*(EijI(2,3)-EijI(9,
     -   3))-EijI(10,3)+EijI(20,3)+2*(EijI(2,4)+EijI(9,2)-EijI(9,4)-E
     -   ijI(10,4)+EijI(28,4)))
       F(182)=DCMPLX(FR(182),FI(182))
       FR(183) = 4*(EijR(2,1)+5*EijR(2,3)-6*EijR(10,3)+4*(EijR(2,2)-Ei
     -   jR(10,4))+EijR(15,3)+2*(EijR(2,4)-EijR(9,2)+EijR(21,4)))
       FI(183) = 4*(EijI(2,1)+5*EijI(2,3)-6*EijI(10,3)+4*(EijI(2,2)-Ei
     -   jI(10,4))+EijI(15,3)+2*(EijI(2,4)-EijI(9,2)+EijI(21,4)))
       F(183)=DCMPLX(FR(183),FI(183))
       P(574) = -p2sq+s12+s23+s45
       P(575) = -s15+s23+3*P(52)
       FR(184) = 2*(3*D02345R-Dij2345R(1,1)-12*EijR(11,2)+8*(EijR(22,3
     -   )-EijR(23,3))+4*(EijR(37,4)-EijR(43,4))-EijR(9,2)*P(191)+Eij
     -   R(2,1)*P(224)+2*(Dij2345R(2,1)+(EE0R+EijR(1,1))*P(21)-EijR(3
     -   ,2)*P(40)+EijR(10,2)*P(52)+(EijR(5,2)-EijR(6,2))*P(338))-Eij
     -   R(2,2)*P(561)+EijR(8,2)*P(564)+EijR(3,1)*P(574)+EijR(4,1)*P(
     -   575))
       FI(184) = 2*(3*D02345I-Dij2345I(1,1)-12*EijI(11,2)+8*(EijI(22,3
     -   )-EijI(23,3))+4*(EijI(37,4)-EijI(43,4))-EijI(9,2)*P(191)+Eij
     -   I(2,1)*P(224)+2*(Dij2345I(2,1)+(EE0I+EijI(1,1))*P(21)-EijI(3
     -   ,2)*P(40)+EijI(10,2)*P(52)+(EijI(5,2)-EijI(6,2))*P(338))-Eij
     -   I(2,2)*P(561)+EijI(8,2)*P(564)+EijI(3,1)*P(574)+EijI(4,1)*P(
     -   575))
       F(184)=DCMPLX(FR(184),FI(184))
       P(576) = 2*p2sq-s12
       P(577) = -p2sq+s12+s23-3*s45+2*P(21)
       P(578) = p2sq-s12-s23+2*s34+s45
       P(579) = -s23+s34+s45
       P(580) = 5*p2sq-3*s12+2*P(579)
       FR(185) = -2*(3*D02345R+Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(4,
     -   2)-12*EijR(11,2)-8*(EijR(22,3)-EijR(23,3))-4*(EijR(37,4)-Eij
     -   R(43,4))-EijR(2,3)*P(143)-EijR(2,1)*P(162)+(2*EijR(9,2)-EijR
     -   (10,2)+EijR(10,3)-EijR(20,3))*P(186)-EijR(12,3)*P(215)-EijR(
     -   9,3)*P(548)+EijR(4,1)*P(575)+2*(Dij2345R(1,1)+(EE0R+EijR(1,1
     -   ))*P(21)+(EijR(5,2)-EijR(6,2))*P(158)-EijR(2,2)*P(576))-EijR
     -   (3,1)*P(577)-EijR(3,2)*P(578)+EijR(8,2)*P(580))
       FI(185) = -2*(3*D02345I+Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(4,
     -   2)-12*EijI(11,2)-8*(EijI(22,3)-EijI(23,3))-4*(EijI(37,4)-Eij
     -   I(43,4))-EijI(2,3)*P(143)-EijI(2,1)*P(162)+(2*EijI(9,2)-EijI
     -   (10,2)+EijI(10,3)-EijI(20,3))*P(186)-EijI(12,3)*P(215)-EijI(
     -   9,3)*P(548)+EijI(4,1)*P(575)+2*(Dij2345I(1,1)+(EE0I+EijI(1,1
     -   ))*P(21)+(EijI(5,2)-EijI(6,2))*P(158)-EijI(2,2)*P(576))-EijI
     -   (3,1)*P(577)-EijI(3,2)*P(578)+EijI(8,2)*P(580))
       F(185)=DCMPLX(FR(185),FI(185))
       P(581) = -2*p2sq+s12
       P(582) = p2sq-s23+2*s34-3*P(56)
       P(583) = p5sq+s15+s23-3*s45+2*P(16)
       P(584) = p2sq-s12-s23+s34+s45
       P(585) = 3*p2sq-s12+2*P(220)
       FR(186) = -2*(3*D02345R+Dij2345R(1,2)-Dij2345R(3,1)-Dij2345R(5,
     -   2)-12*(EijR(11,2)+EijR(22,3))+8*EijR(24,3)-4*(EijR(37,4)-Eij
     -   R(44,4))-EijR(2,3)*P(143)-EijR(2,1)*P(162)-(EijR(4,2)+EijR(1
     -   5,3))*P(186)+(EijR(9,3)-EijR(20,3))*P(215)+EijR(10,3)*P(550)
     -   -EijR(10,2)*P(578)+EijR(3,1)*P(582)+EijR(4,1)*P(583)+2*(Dij2
     -   345R(1,1)+(EE0R+EijR(1,1))*P(21)+(EijR(5,2)-EijR(7,2))*P(158
     -   )+EijR(2,2)*P(581)+EijR(8,2)*P(584))+EijR(9,2)*P(585))
       FI(186) = -2*(3*D02345I+Dij2345I(1,2)-Dij2345I(3,1)-Dij2345I(5,
     -   2)-12*(EijI(11,2)+EijI(22,3))+8*EijI(24,3)-4*(EijI(37,4)-Eij
     -   I(44,4))-EijI(2,3)*P(143)-EijI(2,1)*P(162)-(EijI(4,2)+EijI(1
     -   5,3))*P(186)+(EijI(9,3)-EijI(20,3))*P(215)+EijI(10,3)*P(550)
     -   -EijI(10,2)*P(578)+EijI(3,1)*P(582)+EijI(4,1)*P(583)+2*(Dij2
     -   345I(1,1)+(EE0I+EijI(1,1))*P(21)+(EijI(5,2)-EijI(7,2))*P(158
     -   )+EijI(2,2)*P(581)+EijI(8,2)*P(584))+EijI(9,2)*P(585))
       F(186)=DCMPLX(FR(186),FI(186))
       P(586) = p3sq+s12+s23-2*s45
       P(587) = p3sq+s12+s23+2*P(541)
       P(588) = 3*p3sq+s12-s45
       P(589) = p5sq-s45+4*P(69)
       FR(187) = -2*(EE0R*s45-6*EijR(11,2)-4*(EijR(37,4)-EijR(40,4))-(
     -   EijR(3,1)-2*EijR(8,2))*P(40)+EijR(4,1)*P(52)-2*(Dij2345R(1,1
     -   )+EijR(9,2)*P(97)+EijR(2,2)*P(145))-EijR(1,2)*P(193)-EijR(1,
     -   1)*P(232)+EijR(2,1)*P(586)+EijR(5,2)*P(587)-EijR(6,2)*P(588)
     -   +EijR(7,2)*P(589))
       FI(187) = -2*(EE0I*s45-6*EijI(11,2)-4*(EijI(37,4)-EijI(40,4))-(
     -   EijI(3,1)-2*EijI(8,2))*P(40)+EijI(4,1)*P(52)-2*(Dij2345I(1,1
     -   )+EijI(9,2)*P(97)+EijI(2,2)*P(145))-EijI(1,2)*P(193)-EijI(1,
     -   1)*P(232)+EijI(2,1)*P(586)+EijI(5,2)*P(587)-EijI(6,2)*P(588)
     -   +EijI(7,2)*P(589))
       F(187)=DCMPLX(FR(187),FI(187))
       FR(188) = 2*(2*(EijR(11,2)+3*EijR(22,3)+2*EijR(37,4))+(EijR(4,1
     -   )+2*EijR(9,2))*P(69)+(EE0R+EijR(1,1)+2*(EijR(2,1)+EijR(5,2))
     -   )*P(215))
       FI(188) = 2*(2*(EijI(11,2)+3*EijI(22,3)+2*EijI(37,4))+(EijI(4,1
     -   )+2*EijI(9,2))*P(69)+(EE0I+EijI(1,1)+2*(EijI(2,1)+EijI(5,2))
     -   )*P(215))
       F(188)=DCMPLX(FR(188),FI(188))
       FR(189) = -4*(EE0R+EijR(1,1)+EijR(2,2)+3*(EijR(2,1)-EijR(2,3))+
     -   4*EijR(5,2)+5*EijR(8,3)-2*(EijR(2,4)-EijR(8,4)))
       FI(189) = -4*(EE0I+EijI(1,1)+EijI(2,2)+3*(EijI(2,1)-EijI(2,3))+
     -   4*EijI(5,2)+5*EijI(8,3)-2*(EijI(2,4)-EijI(8,4)))
       F(189)=DCMPLX(FR(189),FI(189))
       FR(190) = 4*(EijR(2,1)+4*EijR(2,2)+5*EijR(2,3)-EijR(8,2)-3*EijR
     -   (9,3)+2*(EijR(2,4)-EijR(9,4)))
       FI(190) = 4*(EijI(2,1)+4*EijI(2,2)+5*EijI(2,3)-EijI(8,2)-3*EijI
     -   (9,3)+2*(EijI(2,4)-EijI(9,4)))
       F(190)=DCMPLX(FR(190),FI(190))
       FR(191) = 4*(EijR(2,1)+4*EijR(2,2)+5*EijR(2,3)-EijR(9,2)-3*EijR
     -   (10,3)+2*(EijR(2,4)-EijR(10,4)))
       FI(191) = 4*(EijI(2,1)+4*EijI(2,2)+5*EijI(2,3)-EijI(9,2)-3*EijI
     -   (10,3)+2*(EijI(2,4)-EijI(10,4)))
       F(191)=DCMPLX(FR(191),FI(191))
       P(590) = p2sq-3*s12-2*P(21)
       P(591) = p2sq-p3sq-s23+s34-2*P(56)
       P(592) = p3sq-s15+s23-s34+2*P(52)
       FR(192) = -2*(3*D02345R+Dij2345R(1,2)+4*(Dij2345R(1,1)-EijR(11,
     -   2)-3*EijR(22,3)-EijR(37,4))-EijR(2,3)*P(143)+EijR(10,3)*P(18
     -   6)+EijR(9,3)*P(215)+EijR(4,1)*P(418)+EijR(3,1)*P(573)-EijR(2
     -   ,1)*P(590)+2*((EE0R+EijR(1,1)+EijR(5,2))*P(21)-EijR(2,2)*P(1
     -   49)+EijR(8,2)*P(591)+EijR(9,2)*P(592)))
       FI(192) = -2*(3*D02345I+Dij2345I(1,2)+4*(Dij2345I(1,1)-EijI(11,
     -   2)-3*EijI(22,3)-EijI(37,4))-EijI(2,3)*P(143)+EijI(10,3)*P(18
     -   6)+EijI(9,3)*P(215)+EijI(4,1)*P(418)+EijI(3,1)*P(573)-EijI(2
     -   ,1)*P(590)+2*((EE0I+EijI(1,1)+EijI(5,2))*P(21)-EijI(2,2)*P(1
     -   49)+EijI(8,2)*P(591)+EijI(9,2)*P(592)))
       F(192)=DCMPLX(FR(192),FI(192))
       FR(193) = 4*(EijR(2,1)+4*EijR(2,2)+5*EijR(2,3)+2*EijR(2,4))
       FI(193) = 4*(EijI(2,1)+4*EijI(2,2)+5*EijI(2,3)+2*EijI(2,4))
       F(193)=DCMPLX(FR(193),FI(193))
       P(593) = p5sq-s12-s15
       P(594) = p3sq+s23-s34-2*P(593)
       P(595) = p3sq*p5sq+p5sq*s12+s12*s34-s15*P(124)
       P(596) = -p2sq**2+s12**2+p3sq*s23-s23*s34+s45**2+p2sq*P(512)-s4
     -   5*P(594)-2*P(595)
       P(597) = -(p5sq*s12)-s12*s15-p3sq*s23-p5sq*s23+s12*s23+s23*s34-
     -   s45**2+p2sq*P(191)+2*P(293)+s45*P(515)
       P(598) = s15-s23-s34
       P(599) = p3sq**2+p2sq*P(52)+2*s34*P(98)-p5sq*P(237)-s45*P(528)+
     -   p3sq*P(529)+s12*P(598)
       P(600) = p5sq+s15+s23+s45
       P(601) = -7*s12+13*s34-3*P(367)+2*P(600)
       P(602) = s23-s45+2*P(143)
       P(603) = p2sq-p3sq+4*s12+s23
       P(604) = -p2sq+p3sq+s23+4*P(56)
       P(605) = s15-s23-4*P(52)-2*P(69)
       FR(194) = (Dij2345R(3,1)*P(204)+D02345R*P(488)+Dij2345R(1,1)*P(
     -   490)-Dij2345R(2,1)*P(491)+EijR(4,2)*P(52)*P(492)+EijR(2,1)*P
     -   (494)-EijR(2,2)*P(497)-EijR(3,1)*P(501)-EijR(3,2)*P(503)+Eij
     -   R(4,1)*P(506)+EijR(5,2)*P(509)+EijR(8,2)*P(522)-EijR(9,2)*P(
     -   526)+EijR(6,2)*P(596)+EijR(7,2)*P(597)+EijR(10,2)*P(599)+2*(
     -   C0345R-4*Dij2345R(7,2)+P(21)*(EE0R*P(16)+EijR(1,2)*P(151))+E
     -   ijR(1,1)*P(148)*P(177)+EijR(11,2)*P(601)+EijR(21,3)*P(602)-E
     -   ijR(22,3)*P(603)+EijR(23,3)*P(604)+EijR(24,3)*P(605)))/2.
       FI(194) = (Dij2345I(3,1)*P(204)+D02345I*P(488)+Dij2345I(1,1)*P(
     -   490)-Dij2345I(2,1)*P(491)+EijI(4,2)*P(52)*P(492)+EijI(2,1)*P
     -   (494)-EijI(2,2)*P(497)-EijI(3,1)*P(501)-EijI(3,2)*P(503)+Eij
     -   I(4,1)*P(506)+EijI(5,2)*P(509)+EijI(8,2)*P(522)-EijI(9,2)*P(
     -   526)+EijI(6,2)*P(596)+EijI(7,2)*P(597)+EijI(10,2)*P(599)+2*(
     -   C0345I-4*Dij2345I(7,2)+P(21)*(EE0I*P(16)+EijI(1,2)*P(151))+E
     -   ijI(1,1)*P(148)*P(177)+EijI(11,2)*P(601)+EijI(21,3)*P(602)-E
     -   ijI(22,3)*P(603)+EijI(23,3)*P(604)+EijI(24,3)*P(605)))/2.
       F(194)=DCMPLX(FR(194),FI(194))
       P(606) = 2*p2sq+p3sq+3*s12-s45
       P(607) = p3sq-s45
       P(608) = p2sq+3*s12+2*P(607)
       FR(195) = -2*(-Dij2345R(1,1)+Dij2345R(1,2)+Dij2345R(2,1)+Dij234
     -   5R(2,2)-2*Dij2345R(4,2)+s45*EijR(3,2)-s45*EijR(8,2)+8*(EijR(
     -   22,3)-EijR(23,3))-EijR(3,3)*P(40)-EijR(13,3)*P(42)+(EijR(10,
     -   3)-2*EijR(20,3))*P(43)-(EijR(8,3)+EijR(11,3)-2*EijR(17,3))*P
     -   (143)+EijR(2,3)*P(145)+(-EijR(9,2)+EijR(10,2))*P(191)+(-EijR
     -   (2,1)+EijR(3,1)-EijR(5,2)+EijR(6,2))*P(193)-EijR(9,3)*P(606)
     -   +EijR(12,3)*P(608))
       FI(195) = -2*(-Dij2345I(1,1)+Dij2345I(1,2)+Dij2345I(2,1)+Dij234
     -   5I(2,2)-2*Dij2345I(4,2)+s45*EijI(3,2)-s45*EijI(8,2)+8*(EijI(
     -   22,3)-EijI(23,3))-EijI(3,3)*P(40)-EijI(13,3)*P(42)+(EijI(10,
     -   3)-2*EijI(20,3))*P(43)-(EijI(8,3)+EijI(11,3)-2*EijI(17,3))*P
     -   (143)+EijI(2,3)*P(145)+(-EijI(9,2)+EijI(10,2))*P(191)+(-EijI
     -   (2,1)+EijI(3,1)-EijI(5,2)+EijI(6,2))*P(193)-EijI(9,3)*P(606)
     -   +EijI(12,3)*P(608))
       F(195)=DCMPLX(FR(195),FI(195))
       P(609) = -s15+s23-3*P(52)
       P(610) = p2sq-s12-s23+3*s45
       P(611) = p2sq-s12+s15-s23+3*P(52)
       P(612) = p2sq-p3sq-p5sq+s12+s34+s45
       P(613) = p2sq-p5sq+2*s12+s34
       FR(196) = D02345R+3*(Dij2345R(1,1)-Dij2345R(3,1))+EijR(2,2)*P(1
     -   51)+EijR(3,1)*P(215)+EijR(4,1)*P(223)-EijR(2,1)*P(224)+EijR(
     -   10,2)*P(562)+EijR(4,2)*P(609)+EijR(8,2)*P(610)+EijR(9,2)*P(6
     -   11)-2*(Dij2345R(1,2)-Dij2345R(4,2)-Dij2345R(5,2)+Dij2345R(6,
     -   2)+8*EijR(22,3)+2*(EijR(11,2)-EijR(23,3))-6*EijR(24,3)-(EE0R
     -   +EijR(1,1))*P(21)-(-EijR(12,3)+EijR(13,3))*P(40)-(EijR(15,3)
     -   -EijR(16,3))*P(43)-(EijR(8,3)-EijR(17,3)-EijR(18,3)+EijR(19,
     -   3))*P(143)+EijR(2,3)*P(145)-(EijR(5,2)-EijR(7,2))*P(193)-Eij
     -   R(9,3)*P(560)-EijR(10,3)*P(612)+EijR(20,3)*P(613))
       FI(196) = D02345I+3*(Dij2345I(1,1)-Dij2345I(3,1))+EijI(2,2)*P(1
     -   51)+EijI(3,1)*P(215)+EijI(4,1)*P(223)-EijI(2,1)*P(224)+EijI(
     -   10,2)*P(562)+EijI(4,2)*P(609)+EijI(8,2)*P(610)+EijI(9,2)*P(6
     -   11)-2*(Dij2345I(1,2)-Dij2345I(4,2)-Dij2345I(5,2)+Dij2345I(6,
     -   2)+8*EijI(22,3)+2*(EijI(11,2)-EijI(23,3))-6*EijI(24,3)-(EE0I
     -   +EijI(1,1))*P(21)-(-EijI(12,3)+EijI(13,3))*P(40)-(EijI(15,3)
     -   -EijI(16,3))*P(43)-(EijI(8,3)-EijI(17,3)-EijI(18,3)+EijI(19,
     -   3))*P(143)+EijI(2,3)*P(145)-(EijI(5,2)-EijI(7,2))*P(193)-Eij
     -   I(9,3)*P(560)-EijI(10,3)*P(612)+EijI(20,3)*P(613))
       F(196)=DCMPLX(FR(196),FI(196))
       P(614) = p3sq-s12+s23-2*s45
       P(615) = p3sq-s12+s23+2*s45
       P(616) = -3*P(52)-2*P(200)
       P(617) = p3sq-s12-s23-2*s45
       P(618) = 2*p2sq+p3sq+s12-s45
       FR(197) = EE0R*s45+3*Dij2345R(2,1)+(EijR(4,1)+EijR(10,2))*P(52)
     -   +EijR(3,1)*P(83)+EijR(3,2)*P(111)-EijR(1,1)*P(232)+EijR(1,2)
     -   *P(511)+(EijR(2,1)+EijR(5,2))*P(614)-EijR(6,2)*P(615)+EijR(7
     -   ,2)*P(616)-EijR(8,2)*P(617)+2*(D02345R-Dij2345R(1,2)+Dij2345
     -   R(4,2)+2*(p2sq*EijR(8,3)+EijR(11,2)+3*(EijR(21,3)-EijR(22,3)
     -   ))+(EijR(11,3)-EijR(12,3))*P(40)+(-EijR(10,3)+EijR(18,3)-Eij
     -   R(19,3)+EijR(20,3))*P(43)+(-EijR(5,3)+EijR(6,3))*P(143)-EijR
     -   (2,3)*P(145)+EijR(9,2)*P(191)+EijR(9,3)*P(560)-EijR(17,3)*P(
     -   618))
       FI(197) = EE0I*s45+3*Dij2345I(2,1)+(EijI(4,1)+EijI(10,2))*P(52)
     -   +EijI(3,1)*P(83)+EijI(3,2)*P(111)-EijI(1,1)*P(232)+EijI(1,2)
     -   *P(511)+(EijI(2,1)+EijI(5,2))*P(614)-EijI(6,2)*P(615)+EijI(7
     -   ,2)*P(616)-EijI(8,2)*P(617)+2*(D02345I-Dij2345I(1,2)+Dij2345
     -   I(4,2)+2*(p2sq*EijI(8,3)+EijI(11,2)+3*(EijI(21,3)-EijI(22,3)
     -   ))+(EijI(11,3)-EijI(12,3))*P(40)+(-EijI(10,3)+EijI(18,3)-Eij
     -   I(19,3)+EijI(20,3))*P(43)+(-EijI(5,3)+EijI(6,3))*P(143)-EijI
     -   (2,3)*P(145)+EijI(9,2)*P(191)+EijI(9,3)*P(560)-EijI(17,3)*P(
     -   618))
       F(197)=DCMPLX(FR(197),FI(197))
       P(619) = p2sq-s12-s23+4*s45
       P(620) = 3*P(143)+2*P(193)
       P(621) = 4*p2sq+3*p3sq+5*s12-s23
       P(622) = s15-s23-3*P(69)
       FR(198) = 3*Dij2345R(2,1)-22*EijR(11,2)-16*EijR(22,3)-(3*EijR(3
     -   ,2)+2*EijR(12,3))*P(40)+(3*EijR(10,2)-2*(EijR(10,3)-EijR(20,
     -   3)))*P(43)+4*(D02345R+EijR(23,3)+EijR(4,1)*P(52))+(-3*EijR(6
     -   ,2)+2*(EijR(8,3)-EijR(17,3)))*P(143)-EijR(2,1)*P(224)+(EE0R+
     -   EijR(1,1))*P(312)-2*(Dij2345R(1,2)-Dij2345R(4,2)+EijR(2,3)*P
     -   (145)+EijR(2,2)*P(313)-EijR(9,3)*P(560))+EijR(3,1)*P(619)+Ei
     -   jR(5,2)*P(620)+EijR(8,2)*P(621)+EijR(9,2)*P(622)
       FI(198) = 3*Dij2345I(2,1)-22*EijI(11,2)-16*EijI(22,3)-(3*EijI(3
     -   ,2)+2*EijI(12,3))*P(40)+(3*EijI(10,2)-2*(EijI(10,3)-EijI(20,
     -   3)))*P(43)+4*(D02345I+EijI(23,3)+EijI(4,1)*P(52))+(-3*EijI(6
     -   ,2)+2*(EijI(8,3)-EijI(17,3)))*P(143)-EijI(2,1)*P(224)+(EE0I+
     -   EijI(1,1))*P(312)-2*(Dij2345I(1,2)-Dij2345I(4,2)+EijI(2,3)*P
     -   (145)+EijI(2,2)*P(313)-EijI(9,3)*P(560))+EijI(3,1)*P(619)+Ei
     -   jI(5,2)*P(620)+EijI(8,2)*P(621)+EijI(9,2)*P(622)
       F(198)=DCMPLX(FR(198),FI(198))
       P(623) = -p2sq+s12-2*s15+s23+s45
       P(624) = -s23-s45+2*P(168)
       P(625) = -3*P(143)+2*P(193)
       FR(199) = Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(2,1)+Dij2345R(2,
     -   2)-4*(EijR(22,3)-EijR(23,3))-2*(Dij2345R(4,2)-(EijR(2,1)-Eij
     -   R(3,1)+EijR(5,2)-EijR(6,2))*P(21))+(EijR(2,2)+EijR(2,3))*P(1
     -   51)+(EijR(10,3)+EijR(13,3)-2*EijR(20,3))*P(186)+EijR(10,2)*P
     -   (191)+EijR(3,3)*P(215)+EijR(9,3)*P(370)+EijR(9,2)*P(516)+Eij
     -   R(3,2)*P(623)+EijR(8,2)*P(624)+EijR(12,3)*P(625)
       FI(199) = Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(2,1)+Dij2345I(2,
     -   2)-4*(EijI(22,3)-EijI(23,3))-2*(Dij2345I(4,2)-(EijI(2,1)-Eij
     -   I(3,1)+EijI(5,2)-EijI(6,2))*P(21))+(EijI(2,2)+EijI(2,3))*P(1
     -   51)+(EijI(10,3)+EijI(13,3)-2*EijI(20,3))*P(186)+EijI(10,2)*P
     -   (191)+EijI(3,3)*P(215)+EijI(9,3)*P(370)+EijI(9,2)*P(516)+Eij
     -   I(3,2)*P(623)+EijI(8,2)*P(624)+EijI(12,3)*P(625)
       F(199)=DCMPLX(FR(199),FI(199))
       P(626) = 2*s15-s45
       P(627) = -p2sq+p3sq+s23
       P(628) = s23-s45-2*P(21)
       P(629) = p3sq+s12+s45-2*P(489)
       P(630) = p3sq-p5sq-s12-s15+s45
       P(631) = p3sq-s12+s15+s23-s45
       FR(200) = Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(3,1)-Dij2345R(4,
     -   2)-Dij2345R(5,2)+Dij2345R(6,2)+6*EijR(11,2)-4*(EijR(22,3)-Ei
     -   jR(24,3))+EijR(4,1)*P(19)+2*(EijR(2,1)+EijR(5,2))*P(21)+EijR
     -   (4,2)*P(52)+(EijR(2,2)+EijR(2,3))*P(151)-EijR(20,3)*P(159)+E
     -   ijR(16,3)*P(186)+EijR(13,3)*P(215)+EijR(15,3)*P(322)+EijR(9,
     -   3)*P(338)+EijR(12,3)*P(406)+EijR(6,2)*P(511)+EijR(10,3)*P(55
     -   0)+EijR(3,1)*P(626)+EijR(3,2)*P(627)+EijR(7,2)*P(628)-EijR(8
     -   ,2)*P(629)+EijR(9,2)*P(630)-EijR(10,2)*P(631)
       FI(200) = Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(3,1)-Dij2345I(4,
     -   2)-Dij2345I(5,2)+Dij2345I(6,2)+6*EijI(11,2)-4*(EijI(22,3)-Ei
     -   jI(24,3))+EijI(4,1)*P(19)+2*(EijI(2,1)+EijI(5,2))*P(21)+EijI
     -   (4,2)*P(52)+(EijI(2,2)+EijI(2,3))*P(151)-EijI(20,3)*P(159)+E
     -   ijI(16,3)*P(186)+EijI(13,3)*P(215)+EijI(15,3)*P(322)+EijI(9,
     -   3)*P(338)+EijI(12,3)*P(406)+EijI(6,2)*P(511)+EijI(10,3)*P(55
     -   0)+EijI(3,1)*P(626)+EijI(3,2)*P(627)+EijI(7,2)*P(628)-EijI(8
     -   ,2)*P(629)+EijI(9,2)*P(630)-EijI(10,2)*P(631)
       F(200)=DCMPLX(FR(200),FI(200))
       P(632) = p2sq-s12+2*s15-s23-s45
       FR(201) = D02345R+Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(4,2)+4*E
     -   ijR(21,3)-6*(EijR(11,2)+EijR(22,3))+2*(Dij2345R(1,1)+EijR(23
     -   ,3)-(EijR(1,1)+EijR(1,2)-EijR(2,1))*P(21))+EijR(8,3)*P(143)+
     -   (EijR(2,2)+EijR(2,3))*P(151)+EijR(5,2)*P(154)+(EijR(10,3)+Ei
     -   jR(19,3))*P(186)+EijR(7,2)*P(191)+EijR(11,3)*P(215)+(EijR(18
     -   ,3)+EijR(20,3))*P(322)+EijR(9,3)*P(338)+EijR(12,3)*P(406)+Ei
     -   jR(9,2)*P(516)+EijR(17,3)*P(548)+EijR(6,2)*P(623)+EijR(8,2)*
     -   P(632)
       FI(201) = D02345I+Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(4,2)+4*E
     -   ijI(21,3)-6*(EijI(11,2)+EijI(22,3))+2*(Dij2345I(1,1)+EijI(23
     -   ,3)-(EijI(1,1)+EijI(1,2)-EijI(2,1))*P(21))+EijI(8,3)*P(143)+
     -   (EijI(2,2)+EijI(2,3))*P(151)+EijI(5,2)*P(154)+(EijI(10,3)+Ei
     -   jI(19,3))*P(186)+EijI(7,2)*P(191)+EijI(11,3)*P(215)+(EijI(18
     -   ,3)+EijI(20,3))*P(322)+EijI(9,3)*P(338)+EijI(12,3)*P(406)+Ei
     -   jI(9,2)*P(516)+EijI(17,3)*P(548)+EijI(6,2)*P(623)+EijI(8,2)*
     -   P(632)
       F(201)=DCMPLX(FR(201),FI(201))
       P(633) = 2*s15-3*s45
       P(634) = -p2sq+3*p3sq+s23+2*P(56)
       P(635) = 3*p3sq+5*s12-2*s15-s45
       P(636) = p3sq+p5sq-P(98)
       P(637) = s15-s23-2*P(636)
       FR(202) = 3*Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(4,2)+18*EijR(1
     -   1,2)-4*EijR(22,3)+EijR(2,3)*P(151)+EijR(10,3)*P(186)+EijR(20
     -   ,3)*P(322)+(EijR(6,2)+EijR(9,3))*P(338)-2*(D02345R+Dij2345R(
     -   2,1)-EijR(2,1)*P(21)+EijR(4,1)*P(52)+EijR(5,2)*P(394))+EijR(
     -   12,3)*P(406)+EijR(2,2)*P(554)+EijR(9,2)*P(565)+EijR(3,1)*P(6
     -   33)+EijR(3,2)*P(634)-EijR(8,2)*P(635)+EijR(10,2)*P(637)
       FI(202) = 3*Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(4,2)+18*EijI(1
     -   1,2)-4*EijI(22,3)+EijI(2,3)*P(151)+EijI(10,3)*P(186)+EijI(20
     -   ,3)*P(322)+(EijI(6,2)+EijI(9,3))*P(338)-2*(D02345I+Dij2345I(
     -   2,1)-EijI(2,1)*P(21)+EijI(4,1)*P(52)+EijI(5,2)*P(394))+EijI(
     -   12,3)*P(406)+EijI(2,2)*P(554)+EijI(9,2)*P(565)+EijI(3,1)*P(6
     -   33)+EijI(3,2)*P(634)-EijI(8,2)*P(635)+EijI(10,2)*P(637)
       F(202)=DCMPLX(FR(202),FI(202))
       P(638) = p3sq-s12-2*s15-s23+s45
       P(639) = p3sq-s12-s23
       P(640) = p3sq-s12-s15
       P(641) = s23-s45-2*P(640)
       P(642) = -2*s23+s45+3*P(84)
       P(643) = -s23+2*s45+3*P(84)
       FR(203) = -Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(2,1)-Dij2345R(2
     -   ,2)+2*Dij2345R(4,2)+4*(EijR(22,3)-EijR(23,3))+EijR(3,1)*P(17
     -   )+EijR(2,1)*P(19)+(EijR(9,2)+2*EijR(20,3))*P(52)+EijR(3,2)*P
     -   (84)-EijR(3,3)*P(111)+(EijR(8,3)+EijR(11,3)-2*EijR(17,3))*P(
     -   193)+(EijR(10,2)+EijR(10,3)+EijR(13,3))*P(273)+EijR(6,2)*P(4
     -   34)+EijR(5,2)*P(628)+EijR(2,2)*P(638)+EijR(2,3)*P(639)+EijR(
     -   8,2)*P(641)-EijR(9,3)*P(642)+EijR(12,3)*P(643)
       FI(203) = -Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(2,1)-Dij2345I(2
     -   ,2)+2*Dij2345I(4,2)+4*(EijI(22,3)-EijI(23,3))+EijI(3,1)*P(17
     -   )+EijI(2,1)*P(19)+(EijI(9,2)+2*EijI(20,3))*P(52)+EijI(3,2)*P
     -   (84)-EijI(3,3)*P(111)+(EijI(8,3)+EijI(11,3)-2*EijI(17,3))*P(
     -   193)+(EijI(10,2)+EijI(10,3)+EijI(13,3))*P(273)+EijI(6,2)*P(4
     -   34)+EijI(5,2)*P(628)+EijI(2,2)*P(638)+EijI(2,3)*P(639)+EijI(
     -   8,2)*P(641)-EijI(9,3)*P(642)+EijI(12,3)*P(643)
       F(203)=DCMPLX(FR(203),FI(203))
       P(644) = p3sq-p5sq-s12-s23-2*P(296)
       P(645) = s23-s45-2*P(84)
       P(646) = p3sq+p5sq-s12-s23-s45
       P(647) = p5sq-s23+2*P(84)
       FR(204) = -Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(3,1)+Dij2345R(4
     -   ,2)+Dij2345R(5,2)-Dij2345R(6,2)+2*(-3*EijR(11,2)+EijR(22,3)+
     -   EijR(23,3)-2*EijR(24,3))+EijR(4,1)*P(17)+EijR(2,1)*P(19)+Eij
     -   R(15,3)*P(52)+EijR(8,2)*P(83)+EijR(10,2)*P(84)+(EijR(12,3)-E
     -   ijR(13,3))*P(111)+(EijR(8,3)+EijR(19,3))*P(193)+(EijR(4,2)+E
     -   ijR(16,3))*P(273)+EijR(7,2)*P(434)+(EijR(17,3)+EijR(18,3))*P
     -   (511)+EijR(5,2)*P(628)+EijR(2,2)*P(638)+EijR(2,3)*P(639)-Eij
     -   R(9,2)*P(644)+EijR(9,3)*P(645)-EijR(10,3)*P(646)+EijR(20,3)*
     -   P(647)
       FI(204) = -Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(3,1)+Dij2345I(4
     -   ,2)+Dij2345I(5,2)-Dij2345I(6,2)+2*(-3*EijI(11,2)+EijI(22,3)+
     -   EijI(23,3)-2*EijI(24,3))+EijI(4,1)*P(17)+EijI(2,1)*P(19)+Eij
     -   I(15,3)*P(52)+EijI(8,2)*P(83)+EijI(10,2)*P(84)+(EijI(12,3)-E
     -   ijI(13,3))*P(111)+(EijI(8,3)+EijI(19,3))*P(193)+(EijI(4,2)+E
     -   ijI(16,3))*P(273)+EijI(7,2)*P(434)+(EijI(17,3)+EijI(18,3))*P
     -   (511)+EijI(5,2)*P(628)+EijI(2,2)*P(638)+EijI(2,3)*P(639)-Eij
     -   I(9,2)*P(644)+EijI(9,3)*P(645)-EijI(10,3)*P(646)+EijI(20,3)*
     -   P(647)
       F(204)=DCMPLX(FR(204),FI(204))
       P(648) = p2sq+2*p5sq-s12-4*s15+s45
       P(649) = p2sq-s12-s23
       P(650) = -p2sq+s12+s23
       P(651) = p3sq-s12-2*s23+s45
       P(652) = p3sq-s12-s23+s45
       FR(205) = -D02345R-Dij2345R(1,2)+Dij2345R(2,1)+Dij2345R(4,2)+6*
     -   EijR(11,2)-4*(EijR(21,3)-EijR(22,3))+(EijR(18,3)+EijR(20,3))
     -   *P(52)+(-EijR(11,3)+EijR(12,3))*P(111)+EijR(9,2)*P(191)+EijR
     -   (6,3)*P(193)+(EijR(10,3)+EijR(19,3))*P(273)+EijR(2,2)*P(463)
     -   +EijR(5,3)*P(511)+EijR(7,2)*P(516)+EijR(2,3)*P(639)+EijR(9,3
     -   )*P(645)-EijR(5,2)*P(648)+EijR(6,2)*P(649)+EijR(8,2)*P(650)-
     -   EijR(8,3)*P(651)-2*(Dij2345R(1,1)-(EijR(1,1)+EijR(1,2)-EijR(
     -   2,1))*P(21)-EijR(17,3)*P(652))
       FI(205) = -D02345I-Dij2345I(1,2)+Dij2345I(2,1)+Dij2345I(4,2)+6*
     -   EijI(11,2)-4*(EijI(21,3)-EijI(22,3))+(EijI(18,3)+EijI(20,3))
     -   *P(52)+(-EijI(11,3)+EijI(12,3))*P(111)+EijI(9,2)*P(191)+EijI
     -   (6,3)*P(193)+(EijI(10,3)+EijI(19,3))*P(273)+EijI(2,2)*P(463)
     -   +EijI(5,3)*P(511)+EijI(7,2)*P(516)+EijI(2,3)*P(639)+EijI(9,3
     -   )*P(645)-EijI(5,2)*P(648)+EijI(6,2)*P(649)+EijI(8,2)*P(650)-
     -   EijI(8,3)*P(651)-2*(Dij2345I(1,1)-(EijI(1,1)+EijI(1,2)-EijI(
     -   2,1))*P(21)-EijI(17,3)*P(652))
       F(205)=DCMPLX(FR(205),FI(205))
       P(653) = -p3sq+3*s12+s23-s45+2*P(489)
       P(654) = s23-s45+2*P(394)
       P(655) = p3sq+5*s12+2*P(541)
       FR(206) = -3*Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(4,2)+EijR(2,1
     -   )*P(19)+(2*EijR(4,1)+EijR(20,3))*P(52)-EijR(9,2)*P(97)+EijR(
     -   12,3)*P(111)+2*(D02345R+Dij2345R(2,1)+s45*EijR(3,1)-9*EijR(1
     -   1,2)+EijR(22,3)+EijR(23,3)-EijR(3,2)*P(40)+EijR(10,2)*P(43)-
     -   EijR(6,2)*P(143))+EijR(8,3)*P(193)+EijR(10,3)*P(273)+EijR(17
     -   ,3)*P(511)+EijR(2,3)*P(639)+EijR(9,3)*P(645)-EijR(2,2)*P(653
     -   )+EijR(5,2)*P(654)+EijR(8,2)*P(655)
       FI(206) = -3*Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(4,2)+EijI(2,1
     -   )*P(19)+(2*EijI(4,1)+EijI(20,3))*P(52)-EijI(9,2)*P(97)+EijI(
     -   12,3)*P(111)+2*(D02345I+Dij2345I(2,1)+s45*EijI(3,1)-9*EijI(1
     -   1,2)+EijI(22,3)+EijI(23,3)-EijI(3,2)*P(40)+EijI(10,2)*P(43)-
     -   EijI(6,2)*P(143))+EijI(8,3)*P(193)+EijI(10,3)*P(273)+EijI(17
     -   ,3)*P(511)+EijI(2,3)*P(639)+EijI(9,3)*P(645)-EijI(2,2)*P(653
     -   )+EijI(5,2)*P(654)+EijI(8,2)*P(655)
       F(206)=DCMPLX(FR(206),FI(206))
       FR(207) = 4*(-EijR(6,2)+EijR(6,3)+EijR(8,2)+EijR(9,3)-2*EijR(17
     -   ,3))
       FI(207) = 4*(-EijI(6,2)+EijI(6,3)+EijI(8,2)+EijI(9,3)-2*EijI(17
     -   ,3))
       F(207)=DCMPLX(FR(207),FI(207))
       FR(208) = -4*(EijR(2,1)-EijR(3,1)+EijR(3,2)+EijR(5,2)-EijR(6,2)
     -   -EijR(8,2)-EijR(9,3)-EijR(11,3)+EijR(12,3)+EijR(17,3))
       FI(208) = -4*(EijI(2,1)-EijI(3,1)+EijI(3,2)+EijI(5,2)-EijI(6,2)
     -   -EijI(8,2)-EijI(9,3)-EijI(11,3)+EijI(12,3)+EijI(17,3))
       F(208)=DCMPLX(FR(208),FI(208))
       FR(209) = -4*(EijR(2,1)-EijR(3,1)+EijR(5,2)+EijR(6,2)-3*EijR(7,
     -   2)-2*(EijR(8,2)-EijR(9,2))-EijR(9,3)+EijR(10,2)+EijR(17,3)-E
     -   ijR(19,3)+EijR(20,3))
       FI(209) = -4*(EijI(2,1)-EijI(3,1)+EijI(5,2)+EijI(6,2)-3*EijI(7,
     -   2)-2*(EijI(8,2)-EijI(9,2))-EijI(9,3)+EijI(10,2)+EijI(17,3)-E
     -   ijI(19,3)+EijI(20,3))
       F(209)=DCMPLX(FR(209),FI(209))
       FR(210) = -4*(EijR(2,1)-EijR(3,1)+EijR(5,2)+EijR(6,2)-2*EijR(8,
     -   2)-EijR(9,3)+EijR(17,3))
       FI(210) = -4*(EijI(2,1)-EijI(3,1)+EijI(5,2)+EijI(6,2)-2*EijI(8,
     -   2)-EijI(9,3)+EijI(17,3))
       F(210)=DCMPLX(FR(210),FI(210))
       FR(211) = 4*(EijR(2,1)-EijR(3,1)+EijR(5,2)-EijR(6,2)+EijR(9,3)+
     -   EijR(11,3)-EijR(12,3)-EijR(17,3))
       FI(211) = 4*(EijI(2,1)-EijI(3,1)+EijI(5,2)-EijI(6,2)+EijI(9,3)+
     -   EijI(11,3)-EijI(12,3)-EijI(17,3))
       F(211)=DCMPLX(FR(211),FI(211))
       FR(212) = 4*(EijR(3,3)+EijR(9,3)-2*EijR(12,3))
       FI(212) = 4*(EijI(3,3)+EijI(9,3)-2*EijI(12,3))
       F(212)=DCMPLX(FR(212),FI(212))
       FR(213) = -4*(EijR(3,2)-EijR(8,2)-EijR(9,3)+2*(EijR(9,2)-EijR(1
     -   0,2))+EijR(12,3)-EijR(13,3)+EijR(20,3))
       FI(213) = -4*(EijI(3,2)-EijI(8,2)-EijI(9,3)+2*(EijI(9,2)-EijI(1
     -   0,2))+EijI(12,3)-EijI(13,3)+EijI(20,3))
       F(213)=DCMPLX(FR(213),FI(213))
       FR(214) = -4*(EijR(3,2)-EijR(8,2)-EijR(9,3)+EijR(12,3))
       FI(214) = -4*(EijI(3,2)-EijI(8,2)-EijI(9,3)+EijI(12,3))
       F(214)=DCMPLX(FR(214),FI(214))
       FR(215) = 4*(EijR(2,1)-EijR(4,1)+EijR(5,2)+EijR(6,2)-2*EijR(7,2
     -   )-EijR(8,2)+EijR(9,2)+EijR(9,3)-EijR(17,3)+EijR(19,3)-EijR(2
     -   0,3))
       FI(215) = 4*(EijI(2,1)-EijI(4,1)+EijI(5,2)+EijI(6,2)-2*EijI(7,2
     -   )-EijI(8,2)+EijI(9,2)+EijI(9,3)-EijI(17,3)+EijI(19,3)-EijI(2
     -   0,3))
       F(215)=DCMPLX(FR(215),FI(215))
       FR(216) = 4*(EijR(3,2)-EijR(8,2)+EijR(9,2)+EijR(9,3)-EijR(10,2)
     -   -EijR(12,3)+EijR(13,3)-EijR(20,3))
       FI(216) = 4*(EijI(3,2)-EijI(8,2)+EijI(9,2)+EijI(9,3)-EijI(10,2)
     -   -EijI(12,3)+EijI(13,3)-EijI(20,3))
       F(216)=DCMPLX(FR(216),FI(216))
       FR(217) = 4*(-EijR(3,1)+EijR(4,1)+EijR(4,2)-EijR(9,2)+EijR(9,3)
     -   +EijR(16,3)-2*EijR(20,3))
       FI(217) = -4*(EijI(3,1)-EijI(4,1)-EijI(4,2)+EijI(9,2)-EijI(9,3)
     -   -EijI(16,3)+2*EijI(20,3))
       F(217)=DCMPLX(FR(217),FI(217))
       FR(218) = -4*(EijR(3,1)-EijR(4,1)-EijR(9,2)-EijR(9,3)+EijR(10,2
     -   )+EijR(20,3))
       FI(218) = -4*(EijI(3,1)-EijI(4,1)-EijI(9,2)-EijI(9,3)+EijI(10,2
     -   )+EijI(20,3))
       F(218)=DCMPLX(FR(218),FI(218))
       FR(219) = -2*(EijR(3,2)-EijR(8,2))+4*(EijR(9,3)-EijR(12,3))
       FI(219) = -2*(EijI(3,2)-EijI(8,2))+4*(EijI(9,3)-EijI(12,3))
       F(219)=DCMPLX(FR(219),FI(219))
       FR(220) = 2*(EijR(3,1)+3*EijR(8,2)-4*EijR(9,2)-EijR(10,2)-2*(Ei
     -   jR(4,1)-EijR(9,3)+EijR(20,3)))
       FI(220) = 2*(EijI(3,1)+3*EijI(8,2)-4*EijI(9,2)-EijI(10,2)-2*(Ei
     -   jI(4,1)-EijI(9,3)+EijI(20,3)))
       F(220)=DCMPLX(FR(220),FI(220))
       FR(221) = 2*(EE0R+EijR(1,1)-EijR(6,2)+EijR(8,2)+2*(EijR(2,1)+Ei
     -   jR(5,2)+EijR(9,3)-EijR(17,3)))
       FI(221) = 2*(EE0I+EijI(1,1)-EijI(6,2)+EijI(8,2)+2*(EijI(2,1)+Ei
     -   jI(5,2)+EijI(9,3)-EijI(17,3)))
       F(221)=DCMPLX(FR(221),FI(221))
       FR(222) = 3*D02345R-Dij2345R(1,1)-16*EijR(11,2)+4*EijR(23,3)+2*
     -   (Dij2345R(2,1)+(EE0R+EijR(1,1))*P(21)-EijR(3,2)*P(40)+EijR(1
     -   0,2)*P(43))-(EijR(2,1)-2*(EijR(5,2)-EijR(6,2)))*P(143)-EijR(
     -   9,2)*P(213)-EijR(2,2)*P(561)+EijR(8,2)*P(564)+EijR(4,1)*P(57
     -   5)+EijR(3,1)*P(610)
       FI(222) = 3*D02345I-Dij2345I(1,1)-16*EijI(11,2)+4*EijI(23,3)+2*
     -   (Dij2345I(2,1)+(EE0I+EijI(1,1))*P(21)-EijI(3,2)*P(40)+EijI(1
     -   0,2)*P(43))-(EijI(2,1)-2*(EijI(5,2)-EijI(6,2)))*P(143)-EijI(
     -   9,2)*P(213)-EijI(2,2)*P(561)+EijI(8,2)*P(564)+EijI(4,1)*P(57
     -   5)+EijI(3,1)*P(610)
       F(222)=DCMPLX(FR(222),FI(222))
       FR(223) = 2*(EijR(3,1)+3*EijR(8,2)+2*EijR(9,3))
       FI(223) = 2*(EijI(3,1)+3*EijI(8,2)+2*EijI(9,3))
       F(223)=DCMPLX(FR(223),FI(223))
       FR(224) = -4*(EijR(2,2)+EijR(2,3)-EijR(5,2)+EijR(5,3)-2*EijR(8,
     -   3))
       FI(224) = -4*(EijI(2,2)+EijI(2,3)-EijI(5,2)+EijI(5,3)-2*EijI(8,
     -   3))
       F(224)=DCMPLX(FR(224),FI(224))
       FR(225) = -4*(EijR(2,3)+2*(EijR(2,2)-EijR(5,2)+EijR(6,2)-EijR(8
     -   ,2))-EijR(8,3)-EijR(9,3)+EijR(17,3))
       FI(225) = -4*(EijI(2,3)+2*(EijI(2,2)-EijI(5,2)+EijI(6,2)-EijI(8
     -   ,2))-EijI(8,3)-EijI(9,3)+EijI(17,3))
       F(225)=DCMPLX(FR(225),FI(225))
       FR(226) = -4*(EijR(2,3)+2*(EijR(2,2)-EijR(5,2))-EijR(8,3)+3*(Ei
     -   jR(7,2)-EijR(9,2))-EijR(10,3)+EijR(18,3))
       FI(226) = -4*(EijI(2,3)+2*(EijI(2,2)-EijI(5,2))-EijI(8,3)+3*(Ei
     -   jI(7,2)-EijI(9,2))-EijI(10,3)+EijI(18,3))
       F(226)=DCMPLX(FR(226),FI(226))
       FR(227) = -8*(EijR(2,2)-EijR(5,2))-4*(EijR(2,3)-EijR(8,3))
       FI(227) = -8*(EijI(2,2)-EijI(5,2))-4*(EijI(2,3)-EijI(8,3))
       F(227)=DCMPLX(FR(227),FI(227))
       FR(228) = -4*(EijR(2,1)+EijR(2,3)-EijR(3,1)+EijR(5,2)-EijR(6,2)
     -   -EijR(8,3)-EijR(9,3)+EijR(17,3))
       FI(228) = -4*(EijI(2,1)+EijI(2,3)-EijI(3,1)+EijI(5,2)-EijI(6,2)
     -   -EijI(8,3)-EijI(9,3)+EijI(17,3))
       F(228)=DCMPLX(FR(228),FI(228))
       FR(229) = -4*(EijR(2,2)+EijR(2,3)+EijR(3,2)-2*(EijR(8,2)+EijR(9
     -   ,3))+EijR(12,3))
       FI(229) = -4*(EijI(2,2)+EijI(2,3)+EijI(3,2)-2*(EijI(8,2)+EijI(9
     -   ,3))+EijI(12,3))
       F(229)=DCMPLX(FR(229),FI(229))
       FR(230) = -4*(EijR(2,2)+EijR(2,3)-EijR(8,2)-EijR(9,3)-2*(EijR(9
     -   ,2)-EijR(10,2))-EijR(10,3)+EijR(20,3))
       FI(230) = -4*(EijI(2,2)+EijI(2,3)-EijI(8,2)-EijI(9,3)-2*(EijI(9
     -   ,2)-EijI(10,2))-EijI(10,3)+EijI(20,3))
       F(230)=DCMPLX(FR(230),FI(230))
       FR(231) = -4*(EijR(2,2)+EijR(2,3)-EijR(8,2)-EijR(9,3))
       FI(231) = -4*(EijI(2,2)+EijI(2,3)-EijI(8,2)-EijI(9,3))
       F(231)=DCMPLX(FR(231),FI(231))
       FR(232) = -4*(EijR(2,1)-EijR(2,2)+EijR(2,3)-EijR(4,1)+2*(EijR(5
     -   ,2)-EijR(7,2))-EijR(8,3)+EijR(9,2)-EijR(10,3)+EijR(18,3))
       FI(232) = -4*(EijI(2,1)-EijI(2,2)+EijI(2,3)-EijI(4,1)+2*(EijI(5
     -   ,2)-EijI(7,2))-EijI(8,3)+EijI(9,2)-EijI(10,3)+EijI(18,3))
       F(232)=DCMPLX(FR(232),FI(232))
       FR(233) = 4*(EijR(2,1)-EijR(2,3)-EijR(3,1)+EijR(9,3)+EijR(10,3)
     -   -EijR(20,3))
       FI(233) = 4*(EijI(2,1)-EijI(2,3)-EijI(3,1)+EijI(9,3)+EijI(10,3)
     -   -EijI(20,3))
       F(233)=DCMPLX(FR(233),FI(233))
       FR(234) = 4*(EijR(2,1)-EijR(2,3)-EijR(4,1)-EijR(4,2)+EijR(9,2)+
     -   2*EijR(10,3)-EijR(15,3))
       FI(234) = 4*(EijI(2,1)-EijI(2,3)-EijI(4,1)-EijI(4,2)+EijI(9,2)+
     -   2*EijI(10,3)-EijI(15,3))
       F(234)=DCMPLX(FR(234),FI(234))
       FR(235) = 4*(EijR(2,1)-EijR(2,3)-EijR(4,1)+EijR(10,3))
       FI(235) = 4*(EijI(2,1)-EijI(2,3)-EijI(4,1)+EijI(10,3))
       F(235)=DCMPLX(FR(235),FI(235))
       FR(236) = -2*(EijR(2,1)-EijR(3,1))-6*(EijR(2,2)-EijR(8,2))-4*(E
     -   ijR(2,3)-EijR(9,3))
       FI(236) = -2*(EijI(2,1)-EijI(3,1))-6*(EijI(2,2)-EijI(8,2))-4*(E
     -   ijI(2,3)-EijI(9,3))
       F(236)=DCMPLX(FR(236),FI(236))
       FR(237) = -2*(EE0R+EijR(1,1)+EijR(2,2)+EijR(5,2)+2*(EijR(2,1)+E
     -   ijR(2,3)-EijR(8,3)))
       FI(237) = -2*(EE0I+EijI(1,1)+EijI(2,2)+EijI(5,2)+2*(EijI(2,1)+E
     -   ijI(2,3)-EijI(8,3)))
       F(237)=DCMPLX(FR(237),FI(237))
       FR(238) = -2*(EijR(2,1)+3*EijR(2,2)-5*EijR(9,2)+2*(EijR(2,3)-Ei
     -   jR(4,1)-EijR(10,3)))
       FI(238) = -2*(EijI(2,1)+3*EijI(2,2)-5*EijI(9,2)+2*(EijI(2,3)-Ei
     -   jI(4,1)-EijI(10,3)))
       F(238)=DCMPLX(FR(238),FI(238))
       P(656) = s15-s23-3*P(52)
       FR(239) = -3*D02345R+Dij2345R(1,1)+16*EijR(11,2)-4*EijR(22,3)-2
     -   *(Dij2345R(2,1)+(EE0R+EijR(1,1))*P(21)-EijR(3,2)*P(40)+EijR(
     -   10,2)*P(43))+(EijR(2,1)-2*(EijR(5,2)-EijR(6,2)))*P(143)+EijR
     -   (9,2)*P(213)+EijR(2,2)*P(561)+EijR(3,1)*P(562)-EijR(8,2)*P(5
     -   64)+EijR(4,1)*P(656)
       FI(239) = -3*D02345I+Dij2345I(1,1)+16*EijI(11,2)-4*EijI(22,3)-2
     -   *(Dij2345I(2,1)+(EE0I+EijI(1,1))*P(21)-EijI(3,2)*P(40)+EijI(
     -   10,2)*P(43))+(EijI(2,1)-2*(EijI(5,2)-EijI(6,2)))*P(143)+EijI
     -   (9,2)*P(213)+EijI(2,2)*P(561)+EijI(3,1)*P(562)-EijI(8,2)*P(5
     -   64)+EijI(4,1)*P(656)
       F(239)=DCMPLX(FR(239),FI(239))
       FR(240) = -2*(EijR(2,1)+3*EijR(2,2)+2*EijR(2,3))
       FI(240) = -2*(EijI(2,1)+3*EijI(2,2)+2*EijI(2,3))
       F(240)=DCMPLX(FR(240),FI(240))
       FR(241) = 8*(EijR(2,2)+EijR(2,3)-EijR(5,2)+EijR(6,2)-EijR(8,2)-
     -   EijR(8,3)-EijR(11,3)+EijR(12,3)-2*(EijR(9,3)-EijR(17,3)))
       FI(241) = 8*(EijI(2,2)+EijI(2,3)-EijI(5,2)+EijI(6,2)-EijI(8,2)-
     -   EijI(8,3)-EijI(11,3)+EijI(12,3)-2*(EijI(9,3)-EijI(17,3)))
       F(241)=DCMPLX(FR(241),FI(241))
       FR(242) = 8*(EijR(2,2)+EijR(2,3)-EijR(5,2)+EijR(7,2)-EijR(8,3)-
     -   EijR(9,2)-EijR(9,3)-EijR(10,3)+EijR(17,3)+EijR(18,3)-EijR(19
     -   ,3)+EijR(20,3))
       FI(242) = 8*(EijI(2,2)+EijI(2,3)-EijI(5,2)+EijI(7,2)-EijI(8,3)-
     -   EijI(9,2)-EijI(9,3)-EijI(10,3)+EijI(17,3)+EijI(18,3)-EijI(19
     -   ,3)+EijI(20,3))
       F(242)=DCMPLX(FR(242),FI(242))
       FR(243) = 4*(EijR(6,2)-EijR(8,2)+2*(EijR(2,2)+EijR(2,3)-EijR(5,
     -   2)-EijR(8,3)-EijR(9,3)+EijR(17,3)))
       FI(243) = 4*(EijI(6,2)-EijI(8,2)+2*(EijI(2,2)+EijI(2,3)-EijI(5,
     -   2)-EijI(8,3)-EijI(9,3)+EijI(17,3)))
       F(243)=DCMPLX(FR(243),FI(243))
       FR(244) = 8*(EijR(2,2)+EijR(2,3)+EijR(3,2)-EijR(3,3)-2*EijR(8,2
     -   )-3*(EijR(9,3)-EijR(12,3)))
       FI(244) = 8*(EijI(2,2)+EijI(2,3)+EijI(3,2)-EijI(3,3)-2*EijI(8,2
     -   )-3*(EijI(9,3)-EijI(12,3)))
       F(244)=DCMPLX(FR(244),FI(244))
       FR(245) = 8*(EijR(2,2)+EijR(2,3)-EijR(8,2)-EijR(9,2)+EijR(10,2)
     -   -EijR(10,3)+EijR(12,3)-EijR(13,3)-2*(EijR(9,3)-EijR(20,3)))
       FI(245) = 8*(EijI(2,2)+EijI(2,3)-EijI(8,2)-EijI(9,2)+EijI(10,2)
     -   -EijI(10,3)+EijI(12,3)-EijI(13,3)-2*(EijI(9,3)-EijI(20,3)))
       F(245)=DCMPLX(FR(245),FI(245))
       FR(246) = 4*(EijR(3,2)-3*EijR(8,2)-4*EijR(9,3)+2*(EijR(2,2)+Eij
     -   R(2,3)+EijR(12,3)))
       FI(246) = 4*(EijI(3,2)-3*EijI(8,2)-4*EijI(9,3)+2*(EijI(2,2)+Eij
     -   I(2,3)+EijI(12,3)))
       F(246)=DCMPLX(FR(246),FI(246))
       FR(247) = 8*(EijR(2,3)+EijR(6,2)-EijR(7,2)-EijR(8,2)-EijR(8,3)+
     -   EijR(9,2)-EijR(9,3)-EijR(10,3)+EijR(17,3)+EijR(18,3)-EijR(19
     -   ,3)+EijR(20,3))
       FI(247) = 8*(EijI(2,3)+EijI(6,2)-EijI(7,2)-EijI(8,2)-EijI(8,3)+
     -   EijI(9,2)-EijI(9,3)-EijI(10,3)+EijI(17,3)+EijI(18,3)-EijI(19
     -   ,3)+EijI(20,3))
       F(247)=DCMPLX(FR(247),FI(247))
       FR(248) = 8*(EijR(2,2)+EijR(2,3)+EijR(3,2)-EijR(10,3)+EijR(12,3
     -   )-EijR(13,3)-2*(EijR(8,2)+EijR(9,3)-EijR(20,3)))
       FI(248) = 8*(EijI(2,2)+EijI(2,3)+EijI(3,2)-EijI(10,3)+EijI(12,3
     -   )-EijI(13,3)-2*(EijI(8,2)+EijI(9,3)-EijI(20,3)))
       F(248)=DCMPLX(FR(248),FI(248))
       FR(249) = 8*(EijR(2,2)+EijR(2,3)-EijR(8,2)-EijR(9,2)-EijR(9,3)+
     -   EijR(10,2)+EijR(15,3)-EijR(16,3)-2*(EijR(10,3)-EijR(20,3)))
       FI(249) = 8*(EijI(2,2)+EijI(2,3)-EijI(8,2)-EijI(9,2)-EijI(9,3)+
     -   EijI(10,2)+EijI(15,3)-EijI(16,3)-2*(EijI(10,3)-EijI(20,3)))
       F(249)=DCMPLX(FR(249),FI(249))
       FR(250) = 4*(-EijR(3,1)+EijR(4,1)-3*EijR(8,2)+EijR(10,2)+2*(Eij
     -   R(2,2)+EijR(2,3)-EijR(9,3)-EijR(10,3)+EijR(20,3)))
       FI(250) = 4*(-EijI(3,1)+EijI(4,1)-3*EijI(8,2)+EijI(10,2)+2*(Eij
     -   I(2,2)+EijI(2,3)-EijI(9,3)-EijI(10,3)+EijI(20,3)))
       F(250)=DCMPLX(FR(250),FI(250))
       FR(251) = 4*(EijR(2,2)-EijR(5,2)+2*(EijR(2,3)+EijR(6,2)-EijR(8,
     -   2)-EijR(8,3)-EijR(9,3)+EijR(17,3)))
       FI(251) = 4*(EijI(2,2)-EijI(5,2)+2*(EijI(2,3)+EijI(6,2)-EijI(8,
     -   2)-EijI(8,3)-EijI(9,3)+EijI(17,3)))
       F(251)=DCMPLX(FR(251),FI(251))
       FR(252) = 4*(EijR(2,1)+3*EijR(2,2)-EijR(3,1)-5*EijR(8,2)-4*EijR
     -   (9,3)+2*(EijR(2,3)+EijR(3,2)+EijR(12,3)))
       FI(252) = 4*(EijI(2,1)+3*EijI(2,2)-EijI(3,1)-5*EijI(8,2)-4*EijI
     -   (9,3)+2*(EijI(2,3)+EijI(3,2)+EijI(12,3)))
       F(252)=DCMPLX(FR(252),FI(252))
       FR(253) = 4*(EijR(2,1)-EijR(4,1)+3*(EijR(2,2)-EijR(9,2))+2*(Eij
     -   R(2,3)-EijR(8,2)-EijR(9,3)+EijR(10,2)-EijR(10,3)+EijR(20,3))
     -   )
       FI(253) = 4*(EijI(2,1)-EijI(4,1)+3*(EijI(2,2)-EijI(9,2))+2*(Eij
     -   I(2,3)-EijI(8,2)-EijI(9,3)+EijI(10,2)-EijI(10,3)+EijI(20,3))
     -   )
       F(253)=DCMPLX(FR(253),FI(253))
       FR(254) = 4*(EijR(2,1)-EijR(3,1)+3*(EijR(2,2)-EijR(8,2))+2*(Eij
     -   R(2,3)-EijR(9,3)))
       FI(254) = 4*(EijI(2,1)-EijI(3,1)+3*(EijI(2,2)-EijI(8,2))+2*(Eij
     -   I(2,3)-EijI(9,3)))
       F(254)=DCMPLX(FR(254),FI(254))
       P(657) = 5*p2sq+3*s12
       P(658) = -5*p3sq-3*P(56)
       P(659) = s23-s45-4*P(143)
       P(660) = 6*s12-2*s23-3*s45+5*P(367)
       P(661) = -s15+s23-3*P(52)-4*P(69)
       P(662) = 3*P(52)+4*P(69)
       FR(255) = 5*D02345R-3*(Dij2345R(1,1)-Dij2345R(2,1))-26*EijR(11,
     -   2)+8*(EijR(22,3)-EijR(23,3))+2*(EE0R+EijR(1,1))*P(21)-(EijR(
     -   2,1)-4*EijR(5,2))*P(143)+EijR(4,1)*P(575)+EijR(3,1)*P(619)-E
     -   ijR(2,2)*P(657)+EijR(3,2)*P(658)+EijR(6,2)*P(659)+EijR(8,2)*
     -   P(660)+EijR(9,2)*P(661)+EijR(10,2)*P(662)
       FI(255) = 5*D02345I-3*(Dij2345I(1,1)-Dij2345I(2,1))-26*EijI(11,
     -   2)+8*(EijI(22,3)-EijI(23,3))+2*(EE0I+EijI(1,1))*P(21)-(EijI(
     -   2,1)-4*EijI(5,2))*P(143)+EijI(4,1)*P(575)+EijI(3,1)*P(619)-E
     -   ijI(2,2)*P(657)+EijI(3,2)*P(658)+EijI(6,2)*P(659)+EijI(8,2)*
     -   P(660)+EijI(9,2)*P(661)+EijI(10,2)*P(662)
       F(255)=DCMPLX(FR(255),FI(255))
       P(663) = p2sq+p3sq-p5sq-s34+s45
       FR(256) = 2*(Dij2345R(1,1)-Dij2345R(1,2)-Dij2345R(3,1)+Dij2345R
     -   (5,2)+p2sq*EijR(2,2)-p3sq*EijR(3,1)-6*EijR(22,3)+4*(EijR(11,
     -   2)+EijR(24,3))+(EijR(9,3)-EijR(20,3))*P(40)+EijR(15,3)*P(43)
     -   +(-EijR(8,2)+EijR(10,2))*P(56)+EijR(4,2)*P(82)+(EijR(8,3)-Ei
     -   jR(18,3))*P(143)-EijR(2,3)*P(145)+(-EijR(5,2)+EijR(7,2))*P(3
     -   38)-EijR(2,1)*P(569)+EijR(4,1)*P(570)+EijR(10,3)*P(612)-EijR
     -   (9,2)*P(663))
       FI(256) = 2*(Dij2345I(1,1)-Dij2345I(1,2)-Dij2345I(3,1)+Dij2345I
     -   (5,2)+p2sq*EijI(2,2)-p3sq*EijI(3,1)-6*EijI(22,3)+4*(EijI(11,
     -   2)+EijI(24,3))+(EijI(9,3)-EijI(20,3))*P(40)+EijI(15,3)*P(43)
     -   +(-EijI(8,2)+EijI(10,2))*P(56)+EijI(4,2)*P(82)+(EijI(8,3)-Ei
     -   jI(18,3))*P(143)-EijI(2,3)*P(145)+(-EijI(5,2)+EijI(7,2))*P(3
     -   38)-EijI(2,1)*P(569)+EijI(4,1)*P(570)+EijI(10,3)*P(612)-EijI
     -   (9,2)*P(663))
       F(256)=DCMPLX(FR(256),FI(256))
       P(664) = -p2sq+s12+2*P(193)
       P(665) = p2sq-s12-s23-s45
       P(666) = 3*p2sq+s12-s23-s45
       P(667) = p5sq+s15-s23-s45-2*P(69)
       FR(257) = -3*D02345R+Dij2345R(1,1)+EijR(2,2)*P(257)-2*(Dij2345R
     -   (1,2)+Dij2345R(2,1)-Dij2345R(4,2)-8*EijR(11,2)+6*EijR(22,3)-
     -   4*EijR(23,3)+(EE0R+EijR(1,1))*P(21)+EijR(12,3)*P(40)-(-EijR(
     -   10,3)+EijR(20,3))*P(43)-EijR(3,2)*P(56)-EijR(10,2)*P(82)-(Ei
     -   jR(8,3)-EijR(17,3))*P(143)+EijR(2,3)*P(145)-(-EijR(5,2)+EijR
     -   (6,2))*P(338)-EijR(9,3)*P(560))+EijR(4,1)*P(656)+EijR(2,1)*P
     -   (664)+EijR(3,1)*P(665)-EijR(8,2)*P(666)+EijR(9,2)*P(667)
       FI(257) = -3*D02345I+Dij2345I(1,1)+EijI(2,2)*P(257)-2*(Dij2345I
     -   (1,2)+Dij2345I(2,1)-Dij2345I(4,2)-8*EijI(11,2)+6*EijI(22,3)-
     -   4*EijI(23,3)+(EE0I+EijI(1,1))*P(21)+EijI(12,3)*P(40)-(-EijI(
     -   10,3)+EijI(20,3))*P(43)-EijI(3,2)*P(56)-EijI(10,2)*P(82)-(Ei
     -   jI(8,3)-EijI(17,3))*P(143)+EijI(2,3)*P(145)-(-EijI(5,2)+EijI
     -   (6,2))*P(338)-EijI(9,3)*P(560))+EijI(4,1)*P(656)+EijI(2,1)*P
     -   (664)+EijI(3,1)*P(665)-EijI(8,2)*P(666)+EijI(9,2)*P(667)
       F(257)=DCMPLX(FR(257),FI(257))
       P(668) = p2sq+2*s12
       P(669) = p3sq+3*s12+s23+2*P(541)
       P(670) = -5*p3sq-s12+s45
       P(671) = p5sq-s45+6*P(69)
       FR(258) = EE0R*s45-4*(Dij2345R(1,1)-p2sq*EijR(8,3))-6*EijR(11,2
     -   )+8*(EijR(21,3)-EijR(22,3))-(EijR(3,1)-2*(EijR(9,3)-EijR(17,
     -   3)))*P(40)+EijR(4,1)*P(52)-EijR(1,1)*P(232)+EijR(1,2)*P(511)
     -   +EijR(2,1)*P(586)-2*(Dij2345R(1,2)-EijR(10,3)*P(42)-EijR(18,
     -   3)*P(43)+EijR(5,3)*P(143)+EijR(2,3)*P(145)+EijR(9,2)*P(270)-
     -   EijR(8,2)*P(403)+EijR(2,2)*P(668))+EijR(5,2)*P(669)+EijR(6,2
     -   )*P(670)+EijR(7,2)*P(671)
       FI(258) = EE0I*s45-4*(Dij2345I(1,1)-p2sq*EijI(8,3))-6*EijI(11,2
     -   )+8*(EijI(21,3)-EijI(22,3))-(EijI(3,1)-2*(EijI(9,3)-EijI(17,
     -   3)))*P(40)+EijI(4,1)*P(52)-EijI(1,1)*P(232)+EijI(1,2)*P(511)
     -   +EijI(2,1)*P(586)-2*(Dij2345I(1,2)-EijI(10,3)*P(42)-EijI(18,
     -   3)*P(43)+EijI(5,3)*P(143)+EijI(2,3)*P(145)+EijI(9,2)*P(270)-
     -   EijI(8,2)*P(403)+EijI(2,2)*P(668))+EijI(5,2)*P(669)+EijI(6,2
     -   )*P(670)+EijI(7,2)*P(671)
       F(258)=DCMPLX(FR(258),FI(258))
       P(672) = p2sq-s23+s45
       P(673) = s12-2*P(672)
       P(674) = -p5sq+s45-5*P(69)
       FR(259) = -Dij2345R(1,1)+p3sq*EijR(3,1)-4*EijR(11,2)-12*EijR(22
     -   ,3)-2*(Dij2345R(1,2)-EijR(9,3)*P(40)+EijR(10,3)*P(43)+EijR(4
     -   ,1)*P(69)-EijR(8,3)*P(143)+EijR(2,3)*P(145))+(EE0R+EijR(1,1)
     -   )*P(406)-EijR(2,2)*P(554)+EijR(8,2)*P(588)+EijR(5,2)*P(664)+
     -   EijR(2,1)*P(673)+EijR(9,2)*P(674)
       FI(259) = -Dij2345I(1,1)+p3sq*EijI(3,1)-4*EijI(11,2)-12*EijI(22
     -   ,3)-2*(Dij2345I(1,2)-EijI(9,3)*P(40)+EijI(10,3)*P(43)+EijI(4
     -   ,1)*P(69)-EijI(8,3)*P(143)+EijI(2,3)*P(145))+(EE0I+EijI(1,1)
     -   )*P(406)-EijI(2,2)*P(554)+EijI(8,2)*P(588)+EijI(5,2)*P(664)+
     -   EijI(2,1)*P(673)+EijI(9,2)*P(674)
       F(259)=DCMPLX(FR(259),FI(259))
       P(675) = p2sq-5*s12
       P(676) = p2sq-5*s12+2*P(21)
       P(677) = -p2sq+4*p3sq+s23-2*s34+3*P(56)
       P(678) = s15-s23-3*P(52)-4*P(69)
       P(679) = p2sq-4*p3sq-s23+2*s34-3*P(56)
       P(680) = -s15+s23+3*P(52)+4*P(69)
       FR(260) = 5*D02345R+6*Dij2345R(1,1)+Dij2345R(1,2)-2*(EijR(11,2)
     -   -6*EijR(21,3)+7*EijR(22,3)+(EijR(1,1)+EijR(1,2)-EijR(2,1))*P
     -   (21))+EijR(8,3)*P(143)+EijR(2,3)*P(151)+EijR(10,3)*P(186)+Ei
     -   jR(9,3)*P(215)+EijR(18,3)*P(322)+EijR(17,3)*P(406)-EijR(2,2)
     -   *P(675)+EijR(5,2)*P(676)+EijR(6,2)*P(677)+EijR(7,2)*P(678)+E
     -   ijR(8,2)*P(679)+EijR(9,2)*P(680)
       FI(260) = 5*D02345I+6*Dij2345I(1,1)+Dij2345I(1,2)-2*(EijI(11,2)
     -   -6*EijI(21,3)+7*EijI(22,3)+(EijI(1,1)+EijI(1,2)-EijI(2,1))*P
     -   (21))+EijI(8,3)*P(143)+EijI(2,3)*P(151)+EijI(10,3)*P(186)+Ei
     -   jI(9,3)*P(215)+EijI(18,3)*P(322)+EijI(17,3)*P(406)-EijI(2,2)
     -   *P(675)+EijI(5,2)*P(676)+EijI(6,2)*P(677)+EijI(7,2)*P(678)+E
     -   ijI(8,2)*P(679)+EijI(9,2)*P(680)
       F(260)=DCMPLX(FR(260),FI(260))
       P(681) = p2sq-s12-s23+4*s45-2*P(21)
       P(682) = p2sq+p3sq-s23+2*s34
       P(683) = s23-s45-2*P(21)-4*P(143)
       P(684) = 7*p2sq+p3sq-3*s23-2*P(16)
       FR(261) = 5*D02345R+Dij2345R(1,2)-Dij2345R(4,2)-26*EijR(11,2)-1
     -   2*(EijR(22,3)-EijR(23,3))+EijR(2,3)*P(151)-EijR(2,1)*P(162)+
     -   EijR(10,3)*P(186)+(-2*EijR(9,2)+EijR(10,2))*P(200)+2*(Dij234
     -   5R(1,1)-Dij2345R(2,1)+(EE0R+EijR(1,1))*P(21)+EijR(5,2)*P(159
     -   )-EijR(2,2)*P(257))+EijR(20,3)*P(322)+EijR(9,3)*P(338)+EijR(
     -   12,3)*P(406)+EijR(4,1)*P(575)+EijR(3,1)*P(681)-EijR(3,2)*P(6
     -   82)+EijR(6,2)*P(683)+EijR(8,2)*P(684)
       FI(261) = 5*D02345I+Dij2345I(1,2)-Dij2345I(4,2)-26*EijI(11,2)-1
     -   2*(EijI(22,3)-EijI(23,3))+EijI(2,3)*P(151)-EijI(2,1)*P(162)+
     -   EijI(10,3)*P(186)+(-2*EijI(9,2)+EijI(10,2))*P(200)+2*(Dij234
     -   5I(1,1)-Dij2345I(2,1)+(EE0I+EijI(1,1))*P(21)+EijI(5,2)*P(159
     -   )-EijI(2,2)*P(257))+EijI(20,3)*P(322)+EijI(9,3)*P(338)+EijI(
     -   12,3)*P(406)+EijI(4,1)*P(575)+EijI(3,1)*P(681)-EijI(3,2)*P(6
     -   82)+EijI(6,2)*P(683)+EijI(8,2)*P(684)
       F(261)=DCMPLX(FR(261),FI(261))
       P(685) = p5sq+s15+s23+2*P(475)
       P(686) = 5*p2sq+p3sq+s23-2*P(300)
       FR(262) = 5*D02345R+Dij2345R(1,2)-Dij2345R(5,2)-20*EijR(11,2)-1
     -   2*(EijR(22,3)-EijR(24,3))+EijR(2,3)*P(151)-EijR(2,1)*P(162)+
     -   EijR(4,2)*P(200)+EijR(9,3)*P(215)+EijR(15,3)*P(322)+EijR(20,
     -   3)*P(406)+2*(Dij2345R(1,1)-Dij2345R(3,1)+(EE0R+EijR(1,1))*P(
     -   21)+EijR(5,2)*P(159)-EijR(2,2)*P(257)+EijR(8,2)*P(491))+EijR
     -   (10,3)*P(550)+EijR(3,1)*P(582)-EijR(10,2)*P(682)+EijR(7,2)*P
     -   (683)+EijR(4,1)*P(685)+EijR(9,2)*P(686)
       FI(262) = 5*D02345I+Dij2345I(1,2)-Dij2345I(5,2)-20*EijI(11,2)-1
     -   2*(EijI(22,3)-EijI(24,3))+EijI(2,3)*P(151)-EijI(2,1)*P(162)+
     -   EijI(4,2)*P(200)+EijI(9,3)*P(215)+EijI(15,3)*P(322)+EijI(20,
     -   3)*P(406)+2*(Dij2345I(1,1)-Dij2345I(3,1)+(EE0I+EijI(1,1))*P(
     -   21)+EijI(5,2)*P(159)-EijI(2,2)*P(257)+EijI(8,2)*P(491))+EijI
     -   (10,3)*P(550)+EijI(3,1)*P(582)-EijI(10,2)*P(682)+EijI(7,2)*P
     -   (683)+EijI(4,1)*P(685)+EijI(9,2)*P(686)
       F(262)=DCMPLX(FR(262),FI(262))
       P(687) = p2sq-p3sq-s12-s23+s34+s45
       P(688) = p3sq+p5sq-s15+s23-s34-s45
       FR(263) = 5*D02345R+Dij2345R(1,2)-12*(EijR(11,2)+EijR(22,3))+4*
     -   (Dij2345R(1,1)-EijR(2,2)*P(143))+EijR(2,3)*P(151)+EijR(10,3)
     -   *P(186)+EijR(9,3)*P(215)+EijR(4,1)*P(418)+EijR(3,1)*P(573)-E
     -   ijR(2,1)*P(590)+2*((EE0R+EijR(1,1))*P(21)+EijR(5,2)*P(158)+E
     -   ijR(8,2)*P(687)+EijR(9,2)*P(688))
       FI(263) = 5*D02345I+Dij2345I(1,2)-12*(EijI(11,2)+EijI(22,3))+4*
     -   (Dij2345I(1,1)-EijI(2,2)*P(143))+EijI(2,3)*P(151)+EijI(10,3)
     -   *P(186)+EijI(9,3)*P(215)+EijI(4,1)*P(418)+EijI(3,1)*P(573)-E
     -   ijI(2,1)*P(590)+2*((EE0I+EijI(1,1))*P(21)+EijI(5,2)*P(158)+E
     -   ijI(8,2)*P(687)+EijI(9,2)*P(688))
       F(263)=DCMPLX(FR(263),FI(263))
       P(689) = 4*p2sq+p3sq+s12-s23-2*s34+s45
       P(690) = -s23+s45+4*P(143)
       P(691) = p3sq-P(98)
       P(692) = -s23+4*P(145)+2*P(691)
       FR(264) = -Dij2345R(1,2)+3*(Dij2345R(1,1)-Dij2345R(2,1))+Dij234
     -   5R(4,2)+12*EijR(11,2)-10*(EijR(22,3)-EijR(23,3))+(3*(EijR(9,
     -   2)-EijR(10,2))+EijR(20,3))*P(52)+(EijR(2,1)-EijR(3,1))*P(92)
     -   +EijR(12,3)*P(111)+EijR(3,2)*P(139)+EijR(8,3)*P(193)+EijR(10
     -   ,3)*P(273)+EijR(17,3)*P(511)+EijR(2,3)*P(639)+EijR(9,3)*P(64
     -   5)+EijR(5,2)*P(659)+EijR(2,2)*P(689)+EijR(6,2)*P(690)-EijR(8
     -   ,2)*P(692)
       FI(264) = -Dij2345I(1,2)+3*(Dij2345I(1,1)-Dij2345I(2,1))+Dij234
     -   5I(4,2)+12*EijI(11,2)-10*(EijI(22,3)-EijI(23,3))+(3*(EijI(9,
     -   2)-EijI(10,2))+EijI(20,3))*P(52)+(EijI(2,1)-EijI(3,1))*P(92)
     -   +EijI(12,3)*P(111)+EijI(3,2)*P(139)+EijI(8,3)*P(193)+EijI(10
     -   ,3)*P(273)+EijI(17,3)*P(511)+EijI(2,3)*P(639)+EijI(9,3)*P(64
     -   5)+EijI(5,2)*P(659)+EijI(2,2)*P(689)+EijI(6,2)*P(690)-EijI(8
     -   ,2)*P(692)
       F(264)=DCMPLX(FR(264),FI(264))
       P(693) = p2sq+s45
       P(694) = p3sq-3*p5sq+s12-s23-2*s34+4*P(693)
       FR(265) = -Dij2345R(1,2)+3*(Dij2345R(1,1)-Dij2345R(3,1))+Dij234
     -   5R(5,2)+6*EijR(11,2)-12*EijR(22,3)+10*EijR(24,3)+(-3*EijR(4,
     -   2)+EijR(15,3))*P(52)+(EijR(2,1)-EijR(4,1))*P(92)+(-EijR(9,3)
     -   +EijR(20,3))*P(111)+(-EijR(8,2)+EijR(10,2))*P(139)+EijR(8,3)
     -   *P(193)+EijR(18,3)*P(511)+EijR(2,3)*P(639)-EijR(10,3)*P(646)
     -   +EijR(5,2)*P(659)+EijR(2,2)*P(689)+EijR(7,2)*P(690)-EijR(9,2
     -   )*P(694)
       FI(265) = -Dij2345I(1,2)+3*(Dij2345I(1,1)-Dij2345I(3,1))+Dij234
     -   5I(5,2)+6*EijI(11,2)-12*EijI(22,3)+10*EijI(24,3)+(-3*EijI(4,
     -   2)+EijI(15,3))*P(52)+(EijI(2,1)-EijI(4,1))*P(92)+(-EijI(9,3)
     -   +EijI(20,3))*P(111)+(-EijI(8,2)+EijI(10,2))*P(139)+EijI(8,3)
     -   *P(193)+EijI(18,3)*P(511)+EijI(2,3)*P(639)-EijI(10,3)*P(646)
     -   +EijI(5,2)*P(659)+EijI(2,2)*P(689)+EijI(7,2)*P(690)-EijI(9,2
     -   )*P(694)
       F(265)=DCMPLX(FR(265),FI(265))
       P(695) = -s45+2*P(3)
       FR(266) = -Dij2345R(1,1)-Dij2345R(1,2)+10*(EijR(21,3)-EijR(22,3
     -   ))+EijR(18,3)*P(52)-4*(p3sq*EijR(6,2)-p3sq*EijR(8,2)-(EijR(7
     -   ,2)-EijR(9,2))*P(69))+(-EijR(9,3)+EijR(17,3))*P(111)+EijR(10
     -   ,3)*P(273)+EijR(5,3)*P(511)+EijR(2,3)*P(639)-EijR(8,3)*P(651
     -   )-(EijR(2,2)-EijR(5,2))*P(695)
       FI(266) = -Dij2345I(1,1)-Dij2345I(1,2)+10*(EijI(21,3)-EijI(22,3
     -   ))+EijI(18,3)*P(52)-4*(p3sq*EijI(6,2)-p3sq*EijI(8,2)-(EijI(7
     -   ,2)-EijI(9,2))*P(69))+(-EijI(9,3)+EijI(17,3))*P(111)+EijI(10
     -   ,3)*P(273)+EijI(5,3)*P(511)+EijI(2,3)*P(639)-EijI(8,3)*P(651
     -   )-(EijI(2,2)-EijI(5,2))*P(695)
       F(266)=DCMPLX(FR(266),FI(266))
       P(696) = -2*s34+s45
       P(697) = p3sq-s12-s23+s45+2*P(336)
       P(698) = p5sq-s45-2*P(69)
       FR(267) = Dij2345R(1,1)-Dij2345R(1,2)-12*EijR(22,3)+2*(p3sq*Eij
     -   R(3,1)-EijR(11,2)-EijR(4,1)*P(69))+(EijR(8,2)-EijR(9,3))*P(1
     -   11)+EijR(8,3)*P(193)+EijR(10,3)*P(273)+EijR(5,2)*P(548)+EijR
     -   (2,3)*P(639)+EijR(2,1)*P(696)+EijR(2,2)*P(697)+EijR(9,2)*P(6
     -   98)
       FI(267) = Dij2345I(1,1)-Dij2345I(1,2)-12*EijI(22,3)+2*(p3sq*Eij
     -   I(3,1)-EijI(11,2)-EijI(4,1)*P(69))+(EijI(8,2)-EijI(9,3))*P(1
     -   11)+EijI(8,3)*P(193)+EijI(10,3)*P(273)+EijI(5,2)*P(548)+EijI
     -   (2,3)*P(639)+EijI(2,1)*P(696)+EijI(2,2)*P(697)+EijI(9,2)*P(6
     -   98)
       F(267)=DCMPLX(FR(267),FI(267))
       P(699) = -p2sq+s12+s23-4*s45
       P(700) = 5*p3sq+3*P(56)
       P(701) = s15-s23+3*P(52)+4*P(69)
       P(702) = -3*P(52)-4*P(69)
       FR(268) = (-5*D02345R+3*(Dij2345R(1,1)-Dij2345R(2,1))+38*EijR(1
     -   1,2)+(EijR(2,1)-4*EijR(5,2))*P(143)+2*(EE0R+EijR(1,1))*P(177
     -   )+EijR(4,1)*P(656)+EijR(2,2)*P(657)-EijR(8,2)*P(660)+EijR(6,
     -   2)*P(690)+EijR(3,1)*P(699)+EijR(3,2)*P(700)+EijR(9,2)*P(701)
     -   +EijR(10,2)*P(702))/2.
       FI(268) = (-5*D02345I+3*(Dij2345I(1,1)-Dij2345I(2,1))+38*EijI(1
     -   1,2)+(EijI(2,1)-4*EijI(5,2))*P(143)+2*(EE0I+EijI(1,1))*P(177
     -   )+EijI(4,1)*P(656)+EijI(2,2)*P(657)-EijI(8,2)*P(660)+EijI(6,
     -   2)*P(690)+EijI(3,1)*P(699)+EijI(3,2)*P(700)+EijI(9,2)*P(701)
     -   +EijI(10,2)*P(702))/2.
       F(268)=DCMPLX(FR(268),FI(268))
       FR(269) = -4*(Dij1235R(3,2)-Dij1235R(6,2))
       FI(269) = -4*(Dij1235I(3,2)-Dij1235I(6,2))
       F(269)=DCMPLX(FR(269),FI(269))
       FR(270) = -2*(D01235R+Dij1235R(1,1)-2*(Dij1235R(3,1)-Dij1235R(5
     -   ,2)+Dij1235R(6,2)))
       FI(270) = -2*(D01235I+Dij1235I(1,1)-2*(Dij1235I(3,1)-Dij1235I(5
     -   ,2)+Dij1235I(6,2)))
       F(270)=DCMPLX(FR(270),FI(270))
       FR(271) = 2*Is34*(C0235R-Cij235R(1,1)+Cij235R(2,1)-4*(Dij1235R(
     -   7,2)-Dij1235R(12,3)+Dij1235R(13,3))+2*(D01235R+Dij1235R(1,1)
     -   )*P(21)+(Dij1235R(2,1)+2*(Dij1235R(4,2)-Dij1235R(5,2)))*P(14
     -   3)+Dij1235R(3,1)*P(144)-Dij1235R(2,2)*P(145)+Dij1235R(3,2)*P
     -   (146)+Dij1235R(6,2)*P(147))
       FI(271) = 2*Is34*(C0235I-Cij235I(1,1)+Cij235I(2,1)-4*(Dij1235I(
     -   7,2)-Dij1235I(12,3)+Dij1235I(13,3))+2*(D01235I+Dij1235I(1,1)
     -   )*P(21)+(Dij1235I(2,1)+2*(Dij1235I(4,2)-Dij1235I(5,2)))*P(14
     -   3)+Dij1235I(3,1)*P(144)-Dij1235I(2,2)*P(145)+Dij1235I(3,2)*P
     -   (146)+Dij1235I(6,2)*P(147))
       F(271)=DCMPLX(FR(271),FI(271))
       FR(272) = 2*Is34*(-2*C0235R-Cij235R(1,1)+4*(Dij1235R(7,2)-Dij12
     -   35R(11,3)+Dij1235R(12,3))+D01235R*P(16)+Dij1235R(2,2)*P(143)
     -   -Dij1235R(1,1)*P(148)+Dij1235R(2,1)*P(149)+Dij1235R(3,1)*P(1
     -   50)+Dij1235R(1,2)*P(151)+(-Dij1235R(5,2)+Dij1235R(6,2))*P(15
     -   3))
       FI(272) = 2*Is34*(-2*C0235I-Cij235I(1,1)+4*(Dij1235I(7,2)-Dij12
     -   35I(11,3)+Dij1235I(12,3))+D01235I*P(16)+Dij1235I(2,2)*P(143)
     -   -Dij1235I(1,1)*P(148)+Dij1235I(2,1)*P(149)+Dij1235I(3,1)*P(1
     -   50)+Dij1235I(1,2)*P(151)+(-Dij1235I(5,2)+Dij1235I(6,2))*P(15
     -   3))
       F(272)=DCMPLX(FR(272),FI(272))
       FR(273) = 2*(Dij1235R(3,1)+2*Dij1235R(6,2))
       FI(273) = 2*(Dij1235I(3,1)+2*Dij1235I(6,2))
       F(273)=DCMPLX(FR(273),FI(273))
       FR(274) = 2*Is34*(-Cij235R(1,1)+p2sq*Dij1235R(2,1)-6*Dij1235R(7
     -   ,2)+4*Dij1235R(12,3)+2*(C0235R+Dij1235R(4,2)*P(143))-Dij1235
     -   R(2,2)*P(145)-Dij1235R(6,2)*P(146)+(D01235R+Dij1235R(1,1))*P
     -   (154)+Dij1235R(3,1)*P(156))
       FI(274) = 2*Is34*(-Cij235I(1,1)+p2sq*Dij1235I(2,1)-6*Dij1235I(7
     -   ,2)+4*Dij1235I(12,3)+2*(C0235I+Dij1235I(4,2)*P(143))-Dij1235
     -   I(2,2)*P(145)-Dij1235I(6,2)*P(146)+(D01235I+Dij1235I(1,1))*P
     -   (154)+Dij1235I(3,1)*P(156))
       F(274)=DCMPLX(FR(274),FI(274))
       FR(275) = 2*(Dij1235R(1,1)+Dij1235R(1,2)-Dij1235R(2,1)-Dij1235R
     -   (2,2)+2*(Dij1235R(2,3)+Dij1235R(4,3)-Dij1235R(5,2)-Dij1235R(
     -   5,3)+Dij1235R(6,2)-Dij1235R(8,3)-2*(Dij1235R(6,3)-Dij1235R(1
     -   0,3))))
       FI(275) = 2*(Dij1235I(1,1)+Dij1235I(1,2)-Dij1235I(2,1)-Dij1235I
     -   (2,2)+2*(Dij1235I(2,3)+Dij1235I(4,3)-Dij1235I(5,2)-Dij1235I(
     -   5,3)+Dij1235I(6,2)-Dij1235I(8,3)-2*(Dij1235I(6,3)-Dij1235I(1
     -   0,3))))
       F(275)=DCMPLX(FR(275),FI(275))
       FR(276) = 2*(Dij1235R(2,2)-Dij1235R(4,2)+Dij1235R(5,2)-Dij1235R
     -   (6,2)+2*(Dij1235R(2,3)-Dij1235R(6,3)-Dij1235R(7,3)+Dij1235R(
     -   9,3)-2*(Dij1235R(8,3)-Dij1235R(10,3))))
       FI(276) = 2*(Dij1235I(2,2)-Dij1235I(4,2)+Dij1235I(5,2)-Dij1235I
     -   (6,2)+2*(Dij1235I(2,3)-Dij1235I(6,3)-Dij1235I(7,3)+Dij1235I(
     -   9,3)-2*(Dij1235I(8,3)-Dij1235I(10,3))))
       F(276)=DCMPLX(FR(276),FI(276))
       FR(277) = 2*(Dij1235R(2,2)-Dij1235R(4,2)+2*(Dij1235R(2,3)+Dij12
     -   35R(5,2)-Dij1235R(6,2)-Dij1235R(6,3)-Dij1235R(8,3)+Dij1235R(
     -   10,3)))
       FI(277) = 2*(Dij1235I(2,2)-Dij1235I(4,2)+2*(Dij1235I(2,3)+Dij12
     -   35I(5,2)-Dij1235I(6,2)-Dij1235I(6,3)-Dij1235I(8,3)+Dij1235I(
     -   10,3)))
       F(277)=DCMPLX(FR(277),FI(277))
       FR(278) = -4*(Dij1235R(2,1)-Dij1235R(2,3)-Dij1235R(3,1)+Dij1235
     -   R(3,2)+Dij1235R(4,2)-Dij1235R(5,2)-Dij1235R(6,2)+Dij1235R(6,
     -   3)+Dij1235R(7,3)-Dij1235R(9,3)+2*(Dij1235R(8,3)-Dij1235R(10,
     -   3)))
       FI(278) = -4*(Dij1235I(2,1)-Dij1235I(2,3)-Dij1235I(3,1)+Dij1235
     -   I(3,2)+Dij1235I(4,2)-Dij1235I(5,2)-Dij1235I(6,2)+Dij1235I(6,
     -   3)+Dij1235I(7,3)-Dij1235I(9,3)+2*(Dij1235I(8,3)-Dij1235I(10,
     -   3)))
       F(278)=DCMPLX(FR(278),FI(278))
       FR(279) = 4*(Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(3,2)-Dij1235R
     -   (3,3)-2*Dij1235R(6,2)-3*(Dij1235R(8,3)-Dij1235R(9,3)))
       FI(279) = 4*(Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(3,2)-Dij1235I
     -   (3,3)-2*Dij1235I(6,2)-3*(Dij1235I(8,3)-Dij1235I(9,3)))
       F(279)=DCMPLX(FR(279),FI(279))
       FR(280) = 4*(Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(3,2)-2*(Dij12
     -   35R(6,2)+Dij1235R(8,3))+Dij1235R(9,3))
       FI(280) = 4*(Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(3,2)-2*(Dij12
     -   35I(6,2)+Dij1235I(8,3))+Dij1235I(9,3))
       F(280)=DCMPLX(FR(280),FI(280))
       FR(281) = 2*(-Dij1235R(2,1)+Dij1235R(2,2)+Dij1235R(3,1)-Dij1235
     -   R(3,2)-4*Dij1235R(8,3)+2*(Dij1235R(2,3)+Dij1235R(9,3)))
       FI(281) = -2*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(3,1)+Dij1235
     -   I(3,2)+4*Dij1235I(8,3)-2*(Dij1235I(2,3)+Dij1235I(9,3)))
       F(281)=DCMPLX(FR(281),FI(281))
       FR(282) = -2*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(3,1)+Dij1235
     -   R(6,2)-2*(Dij1235R(2,3)-Dij1235R(8,3)))
       FI(282) = -2*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(3,1)+Dij1235
     -   I(6,2)-2*(Dij1235I(2,3)-Dij1235I(8,3)))
       F(282)=DCMPLX(FR(282),FI(282))
       FR(283) = -2*(D01235R+Dij1235R(1,1)+Dij1235R(2,2)-Dij1235R(3,1)
     -   +Dij1235R(4,2)+Dij1235R(5,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,1
     -   )-Dij1235R(2,3)+Dij1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)))
       FI(283) = -2*(D01235I+Dij1235I(1,1)+Dij1235I(2,2)-Dij1235I(3,1)
     -   +Dij1235I(4,2)+Dij1235I(5,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,1
     -   )-Dij1235I(2,3)+Dij1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)))
       F(283)=DCMPLX(FR(283),FI(283))
       FR(284) = -C0235R-Cij235R(1,1)+Cij235R(2,1)+4*(Dij1235R(7,2)+Di
     -   j1235R(12,3)-Dij1235R(13,3))-2*(D01235R+Dij1235R(1,1))*P(21)
     -   +(Dij1235R(2,1)+Dij1235R(2,2))*P(143)-Dij1235R(3,1)*P(157)+D
     -   ij1235R(3,2)*P(158)-Dij1235R(6,2)*P(159)
       FI(284) = -C0235I-Cij235I(1,1)+Cij235I(2,1)+4*(Dij1235I(7,2)+Di
     -   j1235I(12,3)-Dij1235I(13,3))-2*(D01235I+Dij1235I(1,1))*P(21)
     -   +(Dij1235I(2,1)+Dij1235I(2,2))*P(143)-Dij1235I(3,1)*P(157)+D
     -   ij1235I(3,2)*P(158)-Dij1235I(6,2)*P(159)
       F(284)=DCMPLX(FR(284),FI(284))
       FR(285) = 2*Is34*(Cij235R(1,1)-Cij235R(2,1)+6*Dij1235R(7,2)+4*(
     -   Dij1235R(12,3)-Dij1235R(13,3))+Dij1235R(3,2)*P(150)+(Dij1235
     -   R(2,1)-Dij1235R(3,1))*P(160)+Dij1235R(2,2)*P(161)+(-Dij1235R
     -   (4,2)+Dij1235R(5,2))*P(162)-Dij1235R(6,2)*P(164))
       FI(285) = 2*Is34*(Cij235I(1,1)-Cij235I(2,1)+6*Dij1235I(7,2)+4*(
     -   Dij1235I(12,3)-Dij1235I(13,3))+Dij1235I(3,2)*P(150)+(Dij1235
     -   I(2,1)-Dij1235I(3,1))*P(160)+Dij1235I(2,2)*P(161)+(-Dij1235I
     -   (4,2)+Dij1235I(5,2))*P(162)-Dij1235I(6,2)*P(164))
       F(285)=DCMPLX(FR(285),FI(285))
       FR(286) = 6*(Dij1235R(2,2)-Dij1235R(4,2))+8*(Dij1235R(5,2)-Dij1
     -   235R(6,2))
       FI(286) = 6*(Dij1235I(2,2)-Dij1235I(4,2))+8*(Dij1235I(5,2)-Dij1
     -   235I(6,2))
       F(286)=DCMPLX(FR(286),FI(286))
       FR(287) = 6*(Dij1235R(2,1)+Dij1235R(2,2)-Dij1235R(3,1))+8*Dij12
     -   35R(3,2)-14*Dij1235R(6,2)
       FI(287) = 6*(Dij1235I(2,1)+Dij1235I(2,2)-Dij1235I(3,1))+8*Dij12
     -   35I(3,2)-14*Dij1235I(6,2)
       F(287)=DCMPLX(FR(287),FI(287))
       FR(288) = 2*Is34*(C0235R+Cij235R(1,1)-6*Dij1235R(7,2)-4*(Dij123
     -   5R(11,3)-Dij1235R(12,3))-2*(Dij1235R(1,1)+Dij1235R(1,2)-Dij1
     -   235R(2,1))*P(21)-Dij1235R(2,2)*P(165)+Dij1235R(4,2)*P(166)+(
     -   -Dij1235R(5,2)+Dij1235R(6,2))*P(167))
       FI(288) = 2*Is34*(C0235I+Cij235I(1,1)-6*Dij1235I(7,2)-4*(Dij123
     -   5I(11,3)-Dij1235I(12,3))-2*(Dij1235I(1,1)+Dij1235I(1,2)-Dij1
     -   235I(2,1))*P(21)-Dij1235I(2,2)*P(165)+Dij1235I(4,2)*P(166)+(
     -   -Dij1235I(5,2)+Dij1235I(6,2))*P(167))
       F(288)=DCMPLX(FR(288),FI(288))
       FR(289) = -8*Dij1235R(7,2)+2*(C0235R-p2sq*Dij1235R(2,2)+(D01235
     -   R+Dij1235R(1,1))*P(21))+(-Dij1235R(2,1)+Dij1235R(4,2))*P(143
     -   )+Dij1235R(5,2)*P(151)+Dij1235R(3,1)*P(168)-Dij1235R(3,2)*P(
     -   169)+Dij1235R(6,2)*P(170)
       FI(289) = -8*Dij1235I(7,2)+2*(C0235I-p2sq*Dij1235I(2,2)+(D01235
     -   I+Dij1235I(1,1))*P(21))+(-Dij1235I(2,1)+Dij1235I(4,2))*P(143
     -   )+Dij1235I(5,2)*P(151)+Dij1235I(3,1)*P(168)-Dij1235I(3,2)*P(
     -   169)+Dij1235I(6,2)*P(170)
       F(289)=DCMPLX(FR(289),FI(289))
       FR(290) = 6*(Dij1235R(2,1)+Dij1235R(2,2))-4*(Dij1235R(3,1)+2*Di
     -   j1235R(6,2))
       FI(290) = 6*(Dij1235I(2,1)+Dij1235I(2,2))-4*(Dij1235I(3,1)+2*Di
     -   j1235I(6,2))
       F(290)=DCMPLX(FR(290),FI(290))
       FR(291) = 2*Is34*(Cij235R(1,1)+10*Dij1235R(7,2)+4*Dij1235R(12,3
     -   )-Dij1235R(6,2)*P(150)-2*(C0235R+Dij1235R(3,1)*P(155))+Dij12
     -   35R(2,2)*P(161)-Dij1235R(4,2)*P(162)+Dij1235R(2,1)*P(171))
       FI(291) = 2*Is34*(Cij235I(1,1)+10*Dij1235I(7,2)+4*Dij1235I(12,3
     -   )-Dij1235I(6,2)*P(150)-2*(C0235I+Dij1235I(3,1)*P(155))+Dij12
     -   35I(2,2)*P(161)-Dij1235I(4,2)*P(162)+Dij1235I(2,1)*P(171))
       F(291)=DCMPLX(FR(291),FI(291))
       FR(292) = 8*Is34*(Dij1235R(2,2)-Dij1235R(2,3)-Dij1235R(4,2)-Dij
     -   1235R(4,3)+2*Dij1235R(6,3))
       FI(292) = 8*Is34*(Dij1235I(2,2)-Dij1235I(2,3)-Dij1235I(4,2)-Dij
     -   1235I(4,3)+2*Dij1235I(6,3))
       F(292)=DCMPLX(FR(292),FI(292))
       FR(293) = -8*Is34*(Dij1235R(2,3)-Dij1235R(5,2)+Dij1235R(6,2)-Di
     -   j1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3))
       FI(293) = -8*Is34*(Dij1235I(2,3)-Dij1235I(5,2)+Dij1235I(6,2)-Di
     -   j1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3))
       F(293)=DCMPLX(FR(293),FI(293))
       FR(294) = 8*Is34*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(2,3)-Dij
     -   1235R(3,1)+2*(Dij1235R(4,2)-Dij1235R(5,2))+Dij1235R(6,2)+Dij
     -   1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3))
       FI(294) = 8*Is34*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(2,3)-Dij
     -   1235I(3,1)+2*(Dij1235I(4,2)-Dij1235I(5,2))+Dij1235I(6,2)+Dij
     -   1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3))
       F(294)=DCMPLX(FR(294),FI(294))
       FR(295) = -8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Di
     -   j1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+
     -   Dij1235R(9,3))
       FI(295) = -8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Di
     -   j1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+
     -   Dij1235I(9,3))
       F(295)=DCMPLX(FR(295),FI(295))
       FR(296) = -8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+2*
     -   (Dij1235R(2,2)-Dij1235R(6,2))-Dij1235R(8,3))
       FI(296) = -8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+2*
     -   (Dij1235I(2,2)-Dij1235I(6,2))-Dij1235I(8,3))
       F(296)=DCMPLX(FR(296),FI(296))
       FR(297) = -2*Is34*(C0235R-Cij235R(1,1)+Cij235R(2,1)-4*(Dij1235R
     -   (7,2)-Dij1235R(12,3)+Dij1235R(13,3))+2*(D01235R+Dij1235R(1,1
     -   ))*P(21)+(Dij1235R(2,1)+2*(Dij1235R(4,2)-Dij1235R(5,2)))*P(1
     -   43)-Dij1235R(2,2)*P(145)+Dij1235R(3,2)*P(146)-Dij1235R(3,1)*
     -   P(172)-Dij1235R(6,2)*P(173))
       FI(297) = -2*Is34*(C0235I-Cij235I(1,1)+Cij235I(2,1)-4*(Dij1235I
     -   (7,2)-Dij1235I(12,3)+Dij1235I(13,3))+2*(D01235I+Dij1235I(1,1
     -   ))*P(21)+(Dij1235I(2,1)+2*(Dij1235I(4,2)-Dij1235I(5,2)))*P(1
     -   43)-Dij1235I(2,2)*P(145)+Dij1235I(3,2)*P(146)-Dij1235I(3,1)*
     -   P(172)-Dij1235I(6,2)*P(173))
       F(297)=DCMPLX(FR(297),FI(297))
       FR(298) = 2*Is34*(Cij235R(1,1)+p2sq*Dij1235R(1,2)-s12*Dij1235R(
     -   1,2)-p2sq*Dij1235R(2,1)-p2sq*Dij1235R(2,2)+s12*Dij1235R(2,2)
     -   -p5sq*Dij1235R(3,1)+s12*Dij1235R(3,1)+2*(C0235R+s12*Dij1235R
     -   (2,1)-s34*Dij1235R(3,1))+(-p5sq+s12-2*(p2sq-s15)+s34)*Dij123
     -   5R(5,2)-4*(Dij1235R(7,2)-Dij1235R(11,3)+Dij1235R(12,3))+D012
     -   35R*P(36)+Dij1235R(1,1)*P(148)+Dij1235R(6,2)*P(174))
       FI(298) = 2*Is34*(Cij235I(1,1)+p2sq*Dij1235I(1,2)-s12*Dij1235I(
     -   1,2)-p2sq*Dij1235I(2,1)-p2sq*Dij1235I(2,2)+s12*Dij1235I(2,2)
     -   -p5sq*Dij1235I(3,1)+s12*Dij1235I(3,1)+2*(C0235I+s12*Dij1235I
     -   (2,1)-s34*Dij1235I(3,1))+(-p5sq+s12-2*(p2sq-s15)+s34)*Dij123
     -   5I(5,2)-4*(Dij1235I(7,2)-Dij1235I(11,3)+Dij1235I(12,3))+D012
     -   35I*P(36)+Dij1235I(1,1)*P(148)+Dij1235I(6,2)*P(174))
       F(298)=DCMPLX(FR(298),FI(298))
       FR(299) = -4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+
     -   2*(Dij1235R(2,3)-Dij1235R(8,3)))
       FI(299) = -4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+
     -   2*(Dij1235I(2,3)-Dij1235I(8,3)))
       F(299)=DCMPLX(FR(299),FI(299))
       FR(300) = 4*Is34*(D01235R+Dij1235R(1,1)-Dij1235R(2,2)+3*Dij1235
     -   R(4,2)+2*(Dij1235R(2,1)-Dij1235R(2,3)+Dij1235R(6,3)))
       FI(300) = 4*Is34*(D01235I+Dij1235I(1,1)-Dij1235I(2,2)+3*Dij1235
     -   I(4,2)+2*(Dij1235I(2,1)-Dij1235I(2,3)+Dij1235I(6,3)))
       F(300)=DCMPLX(FR(300),FI(300))
       FR(301) = 2*Is34*(3*C0235R+Cij235R(1,1)-8*Dij1235R(7,2)-4*Dij12
     -   35R(12,3)+2*(D01235R+Dij1235R(1,1))*P(21)-Dij1235R(2,2)*P(14
     -   3)+Dij1235R(6,2)*P(158)-Dij1235R(2,1)*P(175)+Dij1235R(3,1)*P
     -   (176))
       FI(301) = 2*Is34*(3*C0235I+Cij235I(1,1)-8*Dij1235I(7,2)-4*Dij12
     -   35I(12,3)+2*(D01235I+Dij1235I(1,1))*P(21)-Dij1235I(2,2)*P(14
     -   3)+Dij1235I(6,2)*P(158)-Dij1235I(2,1)*P(175)+Dij1235I(3,1)*P
     -   (176))
       F(301)=DCMPLX(FR(301),FI(301))
       FR(302) = -2*Is34*(-Cij235R(1,1)+p2sq*Dij1235R(2,1)-6*Dij1235R(
     -   7,2)+4*Dij1235R(12,3)+2*(C0235R+Dij1235R(4,2)*P(143))-Dij123
     -   5R(2,2)*P(145)-Dij1235R(6,2)*P(146)+(D01235R+Dij1235R(1,1))*
     -   P(154)+Dij1235R(3,1)*P(156))
       FI(302) = -2*Is34*(-Cij235I(1,1)+p2sq*Dij1235I(2,1)-6*Dij1235I(
     -   7,2)+4*Dij1235I(12,3)+2*(C0235I+Dij1235I(4,2)*P(143))-Dij123
     -   5I(2,2)*P(145)-Dij1235I(6,2)*P(146)+(D01235I+Dij1235I(1,1))*
     -   P(154)+Dij1235I(3,1)*P(156))
       F(302)=DCMPLX(FR(302),FI(302))
       FR(303) = -4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)+2*Dij1235R(2,3
     -   ))
       FI(303) = -4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)+2*Dij1235I(2,3
     -   ))
       F(303)=DCMPLX(FR(303),FI(303))
       FR(304) = 8*Is34*(-Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(4,2)+Di
     -   j1235R(4,3)-2*Dij1235R(6,3))
       FI(304) = 8*Is34*(-Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(4,2)+Di
     -   j1235I(4,3)-2*Dij1235I(6,3))
       F(304)=DCMPLX(FR(304),FI(304))
       FR(305) = 8*Is34*(Dij1235R(2,3)-Dij1235R(5,2)+Dij1235R(6,2)-Dij
     -   1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3))
       FI(305) = 8*Is34*(Dij1235I(2,3)-Dij1235I(5,2)+Dij1235I(6,2)-Dij
     -   1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3))
       F(305)=DCMPLX(FR(305),FI(305))
       FR(306) = 4*(Dij1235R(3,2)-Dij1235R(6,2))
       FI(306) = 4*(Dij1235I(3,2)-Dij1235I(6,2))
       F(306)=DCMPLX(FR(306),FI(306))
       FR(307) = -8*Is34*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(2,3)-Di
     -   j1235R(3,1)+2*(Dij1235R(4,2)-Dij1235R(5,2))+Dij1235R(6,2)+Di
     -   j1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3))
       FI(307) = -8*Is34*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(2,3)-Di
     -   j1235I(3,1)+2*(Dij1235I(4,2)-Dij1235I(5,2))+Dij1235I(6,2)+Di
     -   j1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3))
       F(307)=DCMPLX(FR(307),FI(307))
       FR(308) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3))
       FI(308) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3))
       F(308)=DCMPLX(FR(308),FI(308))
       FR(309) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+2*(
     -   Dij1235R(2,2)-Dij1235R(6,2))-Dij1235R(8,3))
       FI(309) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+2*(
     -   Dij1235I(2,2)-Dij1235I(6,2))-Dij1235I(8,3))
       F(309)=DCMPLX(FR(309),FI(309))
       FR(310) = 2*(D01235R+Dij1235R(1,1)-2*(Dij1235R(3,1)-Dij1235R(5,
     -   2)+Dij1235R(6,2)))
       FI(310) = 2*(D01235I+Dij1235I(1,1)-2*(Dij1235I(3,1)-Dij1235I(5,
     -   2)+Dij1235I(6,2)))
       F(310)=DCMPLX(FR(310),FI(310))
       FR(311) = -2*(Dij1235R(3,1)+2*Dij1235R(6,2))
       FI(311) = -2*(Dij1235I(3,1)+2*Dij1235I(6,2))
       F(311)=DCMPLX(FR(311),FI(311))
       FR(312) = -4*Is34*(D01235R+Dij1235R(1,1)-Dij1235R(2,2)+3*Dij123
     -   5R(4,2)+2*(Dij1235R(2,1)-Dij1235R(2,3)+Dij1235R(6,3)))
       FI(312) = -4*Is34*(D01235I+Dij1235I(1,1)-Dij1235I(2,2)+3*Dij123
     -   5I(4,2)+2*(Dij1235I(2,1)-Dij1235I(2,3)+Dij1235I(6,3)))
       F(312)=DCMPLX(FR(312),FI(312))
       FR(313) = 4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+2
     -   *(Dij1235R(2,3)-Dij1235R(8,3)))
       FI(313) = 4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+2
     -   *(Dij1235I(2,3)-Dij1235I(8,3)))
       F(313)=DCMPLX(FR(313),FI(313))
       FR(314) = -2*Is34*(3*C0235R+Cij235R(1,1)-8*Dij1235R(7,2)-4*Dij1
     -   235R(12,3)+2*(D01235R+Dij1235R(1,1))*P(21)-Dij1235R(2,2)*P(1
     -   43)+Dij1235R(6,2)*P(158)-Dij1235R(2,1)*P(175)+Dij1235R(3,1)*
     -   P(176))
       FI(314) = -2*Is34*(3*C0235I+Cij235I(1,1)-8*Dij1235I(7,2)-4*Dij1
     -   235I(12,3)+2*(D01235I+Dij1235I(1,1))*P(21)-Dij1235I(2,2)*P(1
     -   43)+Dij1235I(6,2)*P(158)-Dij1235I(2,1)*P(175)+Dij1235I(3,1)*
     -   P(176))
       F(314)=DCMPLX(FR(314),FI(314))
       FR(315) = 4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)+2*Dij1235R(2,3)
     -   )
       FI(315) = 4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)+2*Dij1235I(2,3)
     -   )
       F(315)=DCMPLX(FR(315),FI(315))
       FR(316) = (8*Dij1235R(7,2)+(Dij1235R(2,1)+Dij1235R(5,2))*P(143)
     -   +Dij1235R(4,2)*P(151)-Dij1235R(3,1)*P(168)+Dij1235R(3,2)*P(1
     -   69)-Dij1235R(6,2)*P(170)-2*(C0235R-p2sq*Dij1235R(2,2)-(D0123
     -   5R+Dij1235R(1,1))*P(177)))/2.
       FI(316) = (8*Dij1235I(7,2)+(Dij1235I(2,1)+Dij1235I(5,2))*P(143)
     -   +Dij1235I(4,2)*P(151)-Dij1235I(3,1)*P(168)+Dij1235I(3,2)*P(1
     -   69)-Dij1235I(6,2)*P(170)-2*(C0235I-p2sq*Dij1235I(2,2)-(D0123
     -   5I+Dij1235I(1,1))*P(177)))/2.
       F(316)=DCMPLX(FR(316),FI(316))
       FR(317) = Is34*(C0235R-Cij235R(1,1)+Cij235R(2,1)-4*(Dij1235R(7,
     -   2)-Dij1235R(12,3)+Dij1235R(13,3))+2*(D01235R+Dij1235R(1,1))*
     -   P(21)+(Dij1235R(2,1)+2*(Dij1235R(4,2)-Dij1235R(5,2)))*P(143)
     -   +Dij1235R(3,1)*P(144)-Dij1235R(2,2)*P(145)+Dij1235R(3,2)*P(1
     -   46)+Dij1235R(6,2)*P(147))
       FI(317) = Is34*(C0235I-Cij235I(1,1)+Cij235I(2,1)-4*(Dij1235I(7,
     -   2)-Dij1235I(12,3)+Dij1235I(13,3))+2*(D01235I+Dij1235I(1,1))*
     -   P(21)+(Dij1235I(2,1)+2*(Dij1235I(4,2)-Dij1235I(5,2)))*P(143)
     -   +Dij1235I(3,1)*P(144)-Dij1235I(2,2)*P(145)+Dij1235I(3,2)*P(1
     -   46)+Dij1235I(6,2)*P(147))
       F(317)=DCMPLX(FR(317),FI(317))
       FR(318) = Is34*(-Cij235R(1,1)+s12*Dij1235R(1,2)-s12*Dij1235R(2,
     -   2)-s12*Dij1235R(3,1)-2*(C0235R+s12*Dij1235R(2,1)-s34*Dij1235
     -   R(3,1))-(s12+2*s15+s34)*Dij1235R(5,2)+p5sq*(Dij1235R(3,1)+Di
     -   j1235R(5,2))+p2sq*(-Dij1235R(1,2)+Dij1235R(2,1)+Dij1235R(2,2
     -   )+2*Dij1235R(5,2))+4*(Dij1235R(7,2)-Dij1235R(11,3)+Dij1235R(
     -   12,3))-D01235R*P(36)-Dij1235R(1,1)*P(148)-Dij1235R(6,2)*P(17
     -   4))
       FI(318) = Is34*(-Cij235I(1,1)+s12*Dij1235I(1,2)-s12*Dij1235I(2,
     -   2)-s12*Dij1235I(3,1)-2*(C0235I+s12*Dij1235I(2,1)-s34*Dij1235
     -   I(3,1))-(s12+2*s15+s34)*Dij1235I(5,2)+p5sq*(Dij1235I(3,1)+Di
     -   j1235I(5,2))+p2sq*(-Dij1235I(1,2)+Dij1235I(2,1)+Dij1235I(2,2
     -   )+2*Dij1235I(5,2))+4*(Dij1235I(7,2)-Dij1235I(11,3)+Dij1235I(
     -   12,3))-D01235I*P(36)-Dij1235I(1,1)*P(148)-Dij1235I(6,2)*P(17
     -   4))
       F(318)=DCMPLX(FR(318),FI(318))
       FR(319) = Is34*(-Cij235R(1,1)+p2sq*Dij1235R(2,1)-6*Dij1235R(7,2
     -   )+4*Dij1235R(12,3)+2*(C0235R+Dij1235R(4,2)*P(143))-Dij1235R(
     -   2,2)*P(145)+(D01235R+Dij1235R(1,1))*P(154)+Dij1235R(3,1)*P(1
     -   56)+Dij1235R(6,2)*P(178))
       FI(319) = Is34*(-Cij235I(1,1)+p2sq*Dij1235I(2,1)-6*Dij1235I(7,2
     -   )+4*Dij1235I(12,3)+2*(C0235I+Dij1235I(4,2)*P(143))-Dij1235I(
     -   2,2)*P(145)+(D01235I+Dij1235I(1,1))*P(154)+Dij1235I(3,1)*P(1
     -   56)+Dij1235I(6,2)*P(178))
       F(319)=DCMPLX(FR(319),FI(319))
       FR(320) = Is34*(-5*C0235R+14*Dij1235R(7,2)+2*(p2sq*Dij1235R(2,2
     -   )-(D01235R+Dij1235R(1,1))*P(21))+Dij1235R(4,2)*P(151)-Dij123
     -   5R(6,2)*P(169)+Dij1235R(2,1)*P(179)+Dij1235R(3,1)*P(180))
       FI(320) = Is34*(-5*C0235I+14*Dij1235I(7,2)+2*(p2sq*Dij1235I(2,2
     -   )-(D01235I+Dij1235I(1,1))*P(21))+Dij1235I(4,2)*P(151)-Dij123
     -   5I(6,2)*P(169)+Dij1235I(2,1)*P(179)+Dij1235I(3,1)*P(180))
       F(320)=DCMPLX(FR(320),FI(320))
       FR(321) = Is34*(-Cij235R(1,1)+Cij235R(2,1)-6*Dij1235R(7,2)-4*(D
     -   ij1235R(12,3)-Dij1235R(13,3))+(-Dij1235R(2,1)+Dij1235R(3,1))
     -   *P(160)-Dij1235R(2,2)*P(161)+Dij1235R(4,2)*P(162)+Dij1235R(6
     -   ,2)*P(164)+Dij1235R(3,2)*P(181)+Dij1235R(5,2)*P(182))
       FI(321) = Is34*(-Cij235I(1,1)+Cij235I(2,1)-6*Dij1235I(7,2)-4*(D
     -   ij1235I(12,3)-Dij1235I(13,3))+(-Dij1235I(2,1)+Dij1235I(3,1))
     -   *P(160)-Dij1235I(2,2)*P(161)+Dij1235I(4,2)*P(162)+Dij1235I(6
     -   ,2)*P(164)+Dij1235I(3,2)*P(181)+Dij1235I(5,2)*P(182))
       F(321)=DCMPLX(FR(321),FI(321))
       FR(322) = Is34*(-C0235R-Cij235R(1,1)+6*Dij1235R(7,2)+4*(Dij1235
     -   R(11,3)-Dij1235R(12,3))+2*(Dij1235R(1,1)+Dij1235R(1,2)-Dij12
     -   35R(2,1))*P(21)+Dij1235R(2,2)*P(165)-Dij1235R(4,2)*P(166)-(D
     -   ij1235R(5,2)-Dij1235R(6,2))*P(183))
       FI(322) = Is34*(-C0235I-Cij235I(1,1)+6*Dij1235I(7,2)+4*(Dij1235
     -   I(11,3)-Dij1235I(12,3))+2*(Dij1235I(1,1)+Dij1235I(1,2)-Dij12
     -   35I(2,1))*P(21)+Dij1235I(2,2)*P(165)-Dij1235I(4,2)*P(166)-(D
     -   ij1235I(5,2)-Dij1235I(6,2))*P(183))
       F(322)=DCMPLX(FR(322),FI(322))
       FR(323) = Is34*(-Cij235R(1,1)-10*Dij1235R(7,2)-4*Dij1235R(12,3)
     -   +Dij1235R(6,2)*P(150)+2*(C0235R+Dij1235R(3,1)*P(155))-Dij123
     -   5R(2,2)*P(161)+Dij1235R(4,2)*P(162)-Dij1235R(2,1)*P(171))
       FI(323) = Is34*(-Cij235I(1,1)-10*Dij1235I(7,2)-4*Dij1235I(12,3)
     -   +Dij1235I(6,2)*P(150)+2*(C0235I+Dij1235I(3,1)*P(155))-Dij123
     -   5I(2,2)*P(161)+Dij1235I(4,2)*P(162)-Dij1235I(2,1)*P(171))
       F(323)=DCMPLX(FR(323),FI(323))
       FR(324) = 4*Is34*(-Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(4,2)+Di
     -   j1235R(4,3)-2*Dij1235R(6,3))
       FI(324) = 4*Is34*(-Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(4,2)+Di
     -   j1235I(4,3)-2*Dij1235I(6,3))
       F(324)=DCMPLX(FR(324),FI(324))
       FR(325) = 4*Is34*(Dij1235R(2,3)-Dij1235R(5,2)+Dij1235R(6,2)-Dij
     -   1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3))
       FI(325) = 4*Is34*(Dij1235I(2,3)-Dij1235I(5,2)+Dij1235I(6,2)-Dij
     -   1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3))
       F(325)=DCMPLX(FR(325),FI(325))
       FR(326) = -4*Is34*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(2,3)-Di
     -   j1235R(3,1)+2*(Dij1235R(4,2)-Dij1235R(5,2))+Dij1235R(6,2)+Di
     -   j1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3))
       FI(326) = -4*Is34*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(2,3)-Di
     -   j1235I(3,1)+2*(Dij1235I(4,2)-Dij1235I(5,2))+Dij1235I(6,2)+Di
     -   j1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3))
       F(326)=DCMPLX(FR(326),FI(326))
       FR(327) = 4*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3))
       FI(327) = 4*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3))
       F(327)=DCMPLX(FR(327),FI(327))
       FR(328) = 4*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+2*(
     -   Dij1235R(2,2)-Dij1235R(6,2))-Dij1235R(8,3))
       FI(328) = 4*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+2*(
     -   Dij1235I(2,2)-Dij1235I(6,2))-Dij1235I(8,3))
       F(328)=DCMPLX(FR(328),FI(328))
       FR(329) = -2*Is34*(D01235R+Dij1235R(1,1)-Dij1235R(2,2)+3*Dij123
     -   5R(4,2)+2*(Dij1235R(2,1)-Dij1235R(2,3)+Dij1235R(6,3)))
       FI(329) = -2*Is34*(D01235I+Dij1235I(1,1)-Dij1235I(2,2)+3*Dij123
     -   5I(4,2)+2*(Dij1235I(2,1)-Dij1235I(2,3)+Dij1235I(6,3)))
       F(329)=DCMPLX(FR(329),FI(329))
       FR(330) = 2*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+2
     -   *(Dij1235R(2,3)-Dij1235R(8,3)))
       FI(330) = 2*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+2
     -   *(Dij1235I(2,3)-Dij1235I(8,3)))
       F(330)=DCMPLX(FR(330),FI(330))
       FR(331) = Is34*(-3*C0235R-Cij235R(1,1)+8*Dij1235R(7,2)+4*Dij123
     -   5R(12,3)-2*(D01235R+Dij1235R(1,1))*P(21)+Dij1235R(2,2)*P(143
     -   )+Dij1235R(2,1)*P(175)+Dij1235R(3,1)*P(184)+Dij1235R(6,2)*P(
     -   185))
       FI(331) = Is34*(-3*C0235I-Cij235I(1,1)+8*Dij1235I(7,2)+4*Dij123
     -   5I(12,3)-2*(D01235I+Dij1235I(1,1))*P(21)+Dij1235I(2,2)*P(143
     -   )+Dij1235I(2,1)*P(175)+Dij1235I(3,1)*P(184)+Dij1235I(6,2)*P(
     -   185))
       F(331)=DCMPLX(FR(331),FI(331))
       FR(332) = 2*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)+2*Dij1235R(2,3)
     -   )
       FI(332) = 2*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)+2*Dij1235I(2,3)
     -   )
       F(332)=DCMPLX(FR(332),FI(332))
       FR(333) = Is34*(5*C0235R-14*Dij1235R(7,2)-2*(p2sq*Dij1235R(2,2)
     -   -(D01235R+Dij1235R(1,1))*P(21))+Dij1235R(4,2)*P(143)+Dij1235
     -   R(6,2)*P(169)-Dij1235R(2,1)*P(179)-Dij1235R(3,1)*P(180))
       FI(333) = Is34*(5*C0235I-14*Dij1235I(7,2)-2*(p2sq*Dij1235I(2,2)
     -   -(D01235I+Dij1235I(1,1))*P(21))+Dij1235I(4,2)*P(143)+Dij1235
     -   I(6,2)*P(169)-Dij1235I(2,1)*P(179)-Dij1235I(3,1)*P(180))
       F(333)=DCMPLX(FR(333),FI(333))
       FR(334) = -2*(Dij1235R(3,2)-Dij1235R(6,2))
       FI(334) = -2*(Dij1235I(3,2)-Dij1235I(6,2))
       F(334)=DCMPLX(FR(334),FI(334))
       FR(335) = -D01235R-Dij1235R(1,1)+2*(Dij1235R(3,1)-Dij1235R(5,2)
     -   +Dij1235R(6,2))
       FI(335) = -D01235I-Dij1235I(1,1)+2*(Dij1235I(3,1)-Dij1235I(5,2)
     -   +Dij1235I(6,2))
       F(335)=DCMPLX(FR(335),FI(335))
       FR(336) = Dij1235R(3,1)+2*Dij1235R(6,2)
       FI(336) = Dij1235I(3,1)+2*Dij1235I(6,2)
       F(336)=DCMPLX(FR(336),FI(336))
       FR(337) = -3*(Dij1235R(2,2)-Dij1235R(4,2))-4*(Dij1235R(5,2)-Dij
     -   1235R(6,2))
       FI(337) = -3*(Dij1235I(2,2)-Dij1235I(4,2))-4*(Dij1235I(5,2)-Dij
     -   1235I(6,2))
       F(337)=DCMPLX(FR(337),FI(337))
       FR(338) = -3*(Dij1235R(2,1)+Dij1235R(2,2)-Dij1235R(3,1))-4*Dij1
     -   235R(3,2)+7*Dij1235R(6,2)
       FI(338) = -3*(Dij1235I(2,1)+Dij1235I(2,2)-Dij1235I(3,1))-4*Dij1
     -   235I(3,2)+7*Dij1235I(6,2)
       F(338)=DCMPLX(FR(338),FI(338))
       FR(339) = -3*(Dij1235R(2,1)+Dij1235R(2,2))+2*Dij1235R(3,1)+4*Di
     -   j1235R(6,2)
       FI(339) = -3*(Dij1235I(2,1)+Dij1235I(2,2))+2*Dij1235I(3,1)+4*Di
     -   j1235I(6,2)
       F(339)=DCMPLX(FR(339),FI(339))
       FR(340) = (Is34*(-5*C0235R+14*Dij1235R(7,2)+Dij1235R(4,2)*P(151
     -   )-Dij1235R(6,2)*P(169)+2*(p2sq*Dij1235R(2,2)+(D01235R+Dij123
     -   5R(1,1))*P(177))+Dij1235R(2,1)*P(179)+Dij1235R(3,1)*P(180)))
     -   /2.
       FI(340) = (Is34*(-5*C0235I+14*Dij1235I(7,2)+Dij1235I(4,2)*P(151
     -   )-Dij1235I(6,2)*P(169)+2*(p2sq*Dij1235I(2,2)+(D01235I+Dij123
     -   5I(1,1))*P(177))+Dij1235I(2,1)*P(179)+Dij1235I(3,1)*P(180)))
     -   /2.
       F(340)=DCMPLX(FR(340),FI(340))
       FR(341) = 2*Is12*(C0345R+s12*Dij1345R(2,2)-3*s12*Dij1345R(4,2)-
     -   4*s12*Dij1345R(4,3)+s34*Dij1345R(5,2)-s34*Dij1345R(6,2)+p3sq
     -   *(Dij1345R(1,2)-Dij1345R(4,2)-Dij1345R(5,2)+Dij1345R(6,2))-s
     -   45*(D01345R+Dij1345R(1,2)-Dij1345R(3,1)+Dij1345R(4,2)+2*(Dij
     -   1345R(1,1)-Dij1345R(2,1)-Dij1345R(2,2)-Dij1345R(5,2)+Dij1345
     -   R(6,2)))-p5sq*(Dij1345R(3,1)-2*(Dij1345R(1,2)-Dij1345R(4,2)-
     -   Dij1345R(5,2)+Dij1345R(6,2)))+8*(Dij1345R(11,3)-Dij1345R(12,
     -   3))+(Dij1345R(1,1)-Dij1345R(2,1))*P(23)+2*(s12*Dij1345R(1,2)
     -   -s34*Dij1345R(1,2)+s12*Dij1345R(1,3)-s34*Dij1345R(2,3)+s34*D
     -   ij1345R(4,2)-s34*Dij1345R(4,3)+Dij1345R(7,2)+Dij1345R(6,3)*P
     -   (24)))
       FI(341) = 2*Is12*(C0345I+s12*Dij1345I(2,2)-3*s12*Dij1345I(4,2)-
     -   4*s12*Dij1345I(4,3)+s34*Dij1345I(5,2)-s34*Dij1345I(6,2)+p3sq
     -   *(Dij1345I(1,2)-Dij1345I(4,2)-Dij1345I(5,2)+Dij1345I(6,2))-s
     -   45*(D01345I+Dij1345I(1,2)-Dij1345I(3,1)+Dij1345I(4,2)+2*(Dij
     -   1345I(1,1)-Dij1345I(2,1)-Dij1345I(2,2)-Dij1345I(5,2)+Dij1345
     -   I(6,2)))-p5sq*(Dij1345I(3,1)-2*(Dij1345I(1,2)-Dij1345I(4,2)-
     -   Dij1345I(5,2)+Dij1345I(6,2)))+8*(Dij1345I(11,3)-Dij1345I(12,
     -   3))+(Dij1345I(1,1)-Dij1345I(2,1))*P(23)+2*(s12*Dij1345I(1,2)
     -   -s34*Dij1345I(1,2)+s12*Dij1345I(1,3)-s34*Dij1345I(2,3)+s34*D
     -   ij1345I(4,2)-s34*Dij1345I(4,3)+Dij1345I(7,2)+Dij1345I(6,3)*P
     -   (24)))
       F(341)=DCMPLX(FR(341),FI(341))
       FR(342) = 2*Is12*(-(s45*Dij1345R(1,2))+p3sq*(Dij1345R(1,2)-Dij1
     -   345R(3,1)+Dij1345R(3,2)-2*Dij1345R(5,2))-7*s12*Dij1345R(5,2)
     -   -5*(s12*Dij1345R(3,1)-s34*Dij1345R(3,1)-s34*Dij1345R(5,2))+3
     -   *(s45*Dij1345R(3,1)+s12*Dij1345R(3,2)-s34*Dij1345R(3,2)+s45*
     -   Dij1345R(5,2))+2*(C0345R-s34*Dij1345R(1,2)+s12*Dij1345R(1,3)
     -   -s34*Dij1345R(2,1)-s45*Dij1345R(3,2)+s12*Dij1345R(4,2)-s45*D
     -   ij1345R(4,2)-s34*Dij1345R(4,3)+p5sq*(Dij1345R(1,2)-Dij1345R(
     -   3,1)+Dij1345R(3,2)-2*Dij1345R(5,2))-s12*Dij1345R(6,2)+s45*Di
     -   j1345R(6,2)+s12*Dij1345R(7,3)-s34*Dij1345R(9,3))+4*(s12*Dij1
     -   345R(1,2)-s34*Dij1345R(4,2)-s12*Dij1345R(5,3)+s34*Dij1345R(6
     -   ,2)+Dij1345R(7,2)+s34*Dij1345R(10,3))+8*(Dij1345R(11,3)-Dij1
     -   345R(13,3))+D01345R*P(25)+Dij1345R(1,1)*P(26))
       FI(342) = 2*Is12*(-(s45*Dij1345I(1,2))+p3sq*(Dij1345I(1,2)-Dij1
     -   345I(3,1)+Dij1345I(3,2)-2*Dij1345I(5,2))-7*s12*Dij1345I(5,2)
     -   -5*(s12*Dij1345I(3,1)-s34*Dij1345I(3,1)-s34*Dij1345I(5,2))+3
     -   *(s45*Dij1345I(3,1)+s12*Dij1345I(3,2)-s34*Dij1345I(3,2)+s45*
     -   Dij1345I(5,2))+2*(C0345I-s34*Dij1345I(1,2)+s12*Dij1345I(1,3)
     -   -s34*Dij1345I(2,1)-s45*Dij1345I(3,2)+s12*Dij1345I(4,2)-s45*D
     -   ij1345I(4,2)-s34*Dij1345I(4,3)+p5sq*(Dij1345I(1,2)-Dij1345I(
     -   3,1)+Dij1345I(3,2)-2*Dij1345I(5,2))-s12*Dij1345I(6,2)+s45*Di
     -   j1345I(6,2)+s12*Dij1345I(7,3)-s34*Dij1345I(9,3))+4*(s12*Dij1
     -   345I(1,2)-s34*Dij1345I(4,2)-s12*Dij1345I(5,3)+s34*Dij1345I(6
     -   ,2)+Dij1345I(7,2)+s34*Dij1345I(10,3))+8*(Dij1345I(11,3)-Dij1
     -   345I(13,3))+D01345I*P(25)+Dij1345I(1,1)*P(26))
       F(342)=DCMPLX(FR(342),FI(342))
       FR(343) = 2*Is12*(C0345R+Cij345R(1,1)-Cij345R(2,1)-s12*Dij1345R
     -   (2,2)+p3sq*(Dij1345R(1,2)+Dij1345R(2,2)-2*Dij1345R(4,2))-3*s
     -   12*Dij1345R(4,2)+s34*Dij1345R(5,2)+s45*(-D01345R-Dij1345R(1,
     -   2)+Dij1345R(2,1)+Dij1345R(2,2)-2*(Dij1345R(1,1)-Dij1345R(3,1
     -   ))-Dij1345R(3,2)+Dij1345R(5,2))-s34*Dij1345R(6,2)+p5sq*(2*Di
     -   j1345R(1,2)-3*Dij1345R(3,1)+Dij1345R(3,2)-4*Dij1345R(5,2)+Di
     -   j1345R(6,2))-2*(s34*Dij1345R(1,2)-s12*Dij1345R(1,3)-s34*Dij1
     -   345R(2,2)+s12*Dij1345R(4,3)+s34*Dij1345R(4,3)+s12*Dij1345R(5
     -   ,2)+s12*Dij1345R(5,3)-s12*Dij1345R(6,2)-s34*Dij1345R(6,3)+s3
     -   4*Dij1345R(8,3)-(s12+s34)*Dij1345R(10,3))+8*(Dij1345R(7,2)+D
     -   ij1345R(11,3))+4*(s12*Dij1345R(1,2)-Dij1345R(12,3)-Dij1345R(
     -   13,3))+Dij1345R(1,1)*P(27)-Dij1345R(2,1)*P(28))
       FI(343) = 2*Is12*(C0345I+Cij345I(1,1)-Cij345I(2,1)-s12*Dij1345I
     -   (2,2)+p3sq*(Dij1345I(1,2)+Dij1345I(2,2)-2*Dij1345I(4,2))-3*s
     -   12*Dij1345I(4,2)+s34*Dij1345I(5,2)+s45*(-D01345I-Dij1345I(1,
     -   2)+Dij1345I(2,1)+Dij1345I(2,2)-2*(Dij1345I(1,1)-Dij1345I(3,1
     -   ))-Dij1345I(3,2)+Dij1345I(5,2))-s34*Dij1345I(6,2)+p5sq*(2*Di
     -   j1345I(1,2)-3*Dij1345I(3,1)+Dij1345I(3,2)-4*Dij1345I(5,2)+Di
     -   j1345I(6,2))-2*(s34*Dij1345I(1,2)-s12*Dij1345I(1,3)-s34*Dij1
     -   345I(2,2)+s12*Dij1345I(4,3)+s34*Dij1345I(4,3)+s12*Dij1345I(5
     -   ,2)+s12*Dij1345I(5,3)-s12*Dij1345I(6,2)-s34*Dij1345I(6,3)+s3
     -   4*Dij1345I(8,3)-(s12+s34)*Dij1345I(10,3))+8*(Dij1345I(7,2)+D
     -   ij1345I(11,3))+4*(s12*Dij1345I(1,2)-Dij1345I(12,3)-Dij1345I(
     -   13,3))+Dij1345I(1,1)*P(27)-Dij1345I(2,1)*P(28))
       F(343)=DCMPLX(FR(343),FI(343))
       FR(344) = 2*Is12*(-Cij345R(1,1)+Cij345R(2,1)-s45*Dij1345R(1,2)-
     -   s12*Dij1345R(2,1)+s45*Dij1345R(2,1)-s12*Dij1345R(2,2)+s45*Di
     -   j1345R(2,2)-s34*Dij1345R(3,2)-s45*Dij1345R(3,2)+5*s34*Dij134
     -   5R(5,2)+4*s45*Dij1345R(5,2)+3*(s34*Dij1345R(3,1)-s45*Dij1345
     -   R(4,2)-s12*Dij1345R(5,2))+p5sq*(Dij1345R(3,2)+2*(Dij1345R(1,
     -   2)-Dij1345R(2,1)-Dij1345R(4,2)-Dij1345R(5,2))+Dij1345R(6,2))
     -   +p3sq*(Dij1345R(1,2)-Dij1345R(2,2)-Dij1345R(3,1)+Dij1345R(3,
     -   2)+Dij1345R(4,2)-3*Dij1345R(5,2)+Dij1345R(6,2))+2*(C0345R+s1
     -   2*Dij1345R(1,2)-s34*Dij1345R(1,2)+s12*Dij1345R(1,3)+s34*Dij1
     -   345R(2,2)-s12*Dij1345R(3,1)+s45*Dij1345R(3,1)-s34*Dij1345R(4
     -   ,2)-s12*Dij1345R(4,3)-s34*Dij1345R(4,3)-s12*Dij1345R(5,3)+s1
     -   2*Dij1345R(6,2)-s34*Dij1345R(6,2)+s34*Dij1345R(6,3)-Dij1345R
     -   (7,2)-s34*Dij1345R(8,3)+(s12+s34)*Dij1345R(10,3)-2*(-2*Dij13
     -   45R(11,3)+Dij1345R(12,3)+Dij1345R(13,3)))+D01345R*P(25)+Dij1
     -   345R(1,1)*P(29))
       FI(344) = 2*Is12*(-Cij345I(1,1)+Cij345I(2,1)-s45*Dij1345I(1,2)-
     -   s12*Dij1345I(2,1)+s45*Dij1345I(2,1)-s12*Dij1345I(2,2)+s45*Di
     -   j1345I(2,2)-s34*Dij1345I(3,2)-s45*Dij1345I(3,2)+5*s34*Dij134
     -   5I(5,2)+4*s45*Dij1345I(5,2)+3*(s34*Dij1345I(3,1)-s45*Dij1345
     -   I(4,2)-s12*Dij1345I(5,2))+p5sq*(Dij1345I(3,2)+2*(Dij1345I(1,
     -   2)-Dij1345I(2,1)-Dij1345I(4,2)-Dij1345I(5,2))+Dij1345I(6,2))
     -   +p3sq*(Dij1345I(1,2)-Dij1345I(2,2)-Dij1345I(3,1)+Dij1345I(3,
     -   2)+Dij1345I(4,2)-3*Dij1345I(5,2)+Dij1345I(6,2))+2*(C0345I+s1
     -   2*Dij1345I(1,2)-s34*Dij1345I(1,2)+s12*Dij1345I(1,3)+s34*Dij1
     -   345I(2,2)-s12*Dij1345I(3,1)+s45*Dij1345I(3,1)-s34*Dij1345I(4
     -   ,2)-s12*Dij1345I(4,3)-s34*Dij1345I(4,3)-s12*Dij1345I(5,3)+s1
     -   2*Dij1345I(6,2)-s34*Dij1345I(6,2)+s34*Dij1345I(6,3)-Dij1345I
     -   (7,2)-s34*Dij1345I(8,3)+(s12+s34)*Dij1345I(10,3)-2*(-2*Dij13
     -   45I(11,3)+Dij1345I(12,3)+Dij1345I(13,3)))+D01345I*P(25)+Dij1
     -   345I(1,1)*P(29))
       F(344)=DCMPLX(FR(344),FI(344))
       FR(345) = -2*s12*Dij1345R(1,2)+C0345R*P(30)+Dij1345R(4,2)*P(41)
     -   +Dij1345R(5,2)*P(42)+Is12*(-2*D01345R*P(31)-Dij1345R(1,1)*P(
     -   33)+Dij1345R(2,1)*P(35)+Dij1345R(3,1)*P(39)+s34*(Cij345R(1,1
     -   )-4*Dij1345R(12,3)-Dij1345R(2,2)*P(40)+Dij1345R(6,2)*P(43)))
     -   +4*(Dij1345R(11,3)+Dij1345R(7,2)*P(44))
       FI(345) = -2*s12*Dij1345I(1,2)+C0345I*P(30)+Dij1345I(4,2)*P(41)
     -   +Dij1345I(5,2)*P(42)+Is12*(-2*D01345I*P(31)-Dij1345I(1,1)*P(
     -   33)+Dij1345I(2,1)*P(35)+Dij1345I(3,1)*P(39)+s34*(Cij345I(1,1
     -   )-4*Dij1345I(12,3)-Dij1345I(2,2)*P(40)+Dij1345I(6,2)*P(43)))
     -   +4*(Dij1345I(11,3)+Dij1345I(7,2)*P(44))
       F(345)=DCMPLX(FR(345),FI(345))
       FR(346) = -2*Is12*(C0345R-Cij345R(1,1)+s12*Dij1345R(2,2)+3*s12*
     -   Dij1345R(4,2)-s34*Dij1345R(5,2)+s34*Dij1345R(6,2)+s45*(D0134
     -   5R+Dij1345R(1,2)+Dij1345R(2,1)-Dij1345R(2,2)-3*Dij1345R(3,1)
     -   +2*(Dij1345R(1,1)-Dij1345R(5,2))+Dij1345R(6,2))-p5sq*(-3*Dij
     -   1345R(3,1)+2*(Dij1345R(1,2)-Dij1345R(5,2))+Dij1345R(6,2))-p3
     -   sq*(Dij1345R(1,2)+Dij1345R(2,2)-2*Dij1345R(4,2)-Dij1345R(5,2
     -   )+Dij1345R(6,2))+2*(s34*Dij1345R(1,2)-s12*Dij1345R(1,3)-s34*
     -   Dij1345R(2,2)+s12*Dij1345R(4,3)+s34*Dij1345R(4,3)-s34*Dij134
     -   5R(6,3))-12*Dij1345R(7,2)-8*Dij1345R(11,3)-4*(s12*Dij1345R(1
     -   ,2)-Dij1345R(12,3))-Dij1345R(1,1)*P(23)+Dij1345R(2,1)*P(45))
       FI(346) = -2*Is12*(C0345I-Cij345I(1,1)+s12*Dij1345I(2,2)+3*s12*
     -   Dij1345I(4,2)-s34*Dij1345I(5,2)+s34*Dij1345I(6,2)+s45*(D0134
     -   5I+Dij1345I(1,2)+Dij1345I(2,1)-Dij1345I(2,2)-3*Dij1345I(3,1)
     -   +2*(Dij1345I(1,1)-Dij1345I(5,2))+Dij1345I(6,2))-p5sq*(-3*Dij
     -   1345I(3,1)+2*(Dij1345I(1,2)-Dij1345I(5,2))+Dij1345I(6,2))-p3
     -   sq*(Dij1345I(1,2)+Dij1345I(2,2)-2*Dij1345I(4,2)-Dij1345I(5,2
     -   )+Dij1345I(6,2))+2*(s34*Dij1345I(1,2)-s12*Dij1345I(1,3)-s34*
     -   Dij1345I(2,2)+s12*Dij1345I(4,3)+s34*Dij1345I(4,3)-s34*Dij134
     -   5I(6,3))-12*Dij1345I(7,2)-8*Dij1345I(11,3)-4*(s12*Dij1345I(1
     -   ,2)-Dij1345I(12,3))-Dij1345I(1,1)*P(23)+Dij1345I(2,1)*P(45))
       F(346)=DCMPLX(FR(346),FI(346))
       FR(347) = 2*Is12*(Cij345R(2,1)-s45*Dij1345R(1,2)-s34*Dij1345R(3
     -   ,2)-s45*Dij1345R(3,2)+p5sq*(Dij1345R(3,2)+2*(Dij1345R(1,2)-D
     -   ij1345R(3,1)-Dij1345R(5,2)))+5*(s34*Dij1345R(3,1)-s12*Dij134
     -   5R(5,2)+s34*Dij1345R(5,2))-s12*Dij1345R(6,2)+s45*Dij1345R(6,
     -   2)+p3sq*(Dij1345R(1,1)+Dij1345R(1,2)-Dij1345R(3,1)+Dij1345R(
     -   3,2)-3*Dij1345R(5,2)+Dij1345R(6,2))+8*(Dij1345R(7,2)+Dij1345
     -   R(11,3))+4*(s12*Dij1345R(1,2)-s12*Dij1345R(3,1)+s45*Dij1345R
     -   (3,1)-s34*Dij1345R(4,2)+s45*Dij1345R(5,2)-Dij1345R(13,3))+D0
     -   1345R*P(46)-2*(s34*Dij1345R(1,2)-s12*Dij1345R(1,3)-s12*Dij13
     -   45R(2,1)+s34*Dij1345R(2,1)+s45*Dij1345R(2,1)-s12*Dij1345R(4,
     -   2)+s45*Dij1345R(4,2)+s34*Dij1345R(4,3)+s12*Dij1345R(5,3)-s34
     -   *Dij1345R(6,2)-s34*Dij1345R(10,3)-Dij1345R(1,1)*P(47)))
       FI(347) = 2*Is12*(Cij345I(2,1)-s45*Dij1345I(1,2)-s34*Dij1345I(3
     -   ,2)-s45*Dij1345I(3,2)+p5sq*(Dij1345I(3,2)+2*(Dij1345I(1,2)-D
     -   ij1345I(3,1)-Dij1345I(5,2)))+5*(s34*Dij1345I(3,1)-s12*Dij134
     -   5I(5,2)+s34*Dij1345I(5,2))-s12*Dij1345I(6,2)+s45*Dij1345I(6,
     -   2)+p3sq*(Dij1345I(1,1)+Dij1345I(1,2)-Dij1345I(3,1)+Dij1345I(
     -   3,2)-3*Dij1345I(5,2)+Dij1345I(6,2))+8*(Dij1345I(7,2)+Dij1345
     -   I(11,3))+4*(s12*Dij1345I(1,2)-s12*Dij1345I(3,1)+s45*Dij1345I
     -   (3,1)-s34*Dij1345I(4,2)+s45*Dij1345I(5,2)-Dij1345I(13,3))+D0
     -   1345I*P(46)-2*(s34*Dij1345I(1,2)-s12*Dij1345I(1,3)-s12*Dij13
     -   45I(2,1)+s34*Dij1345I(2,1)+s45*Dij1345I(2,1)-s12*Dij1345I(4,
     -   2)+s45*Dij1345I(4,2)+s34*Dij1345I(4,3)+s12*Dij1345I(5,3)-s34
     -   *Dij1345I(6,2)-s34*Dij1345I(10,3)-Dij1345I(1,1)*P(47)))
       F(347)=DCMPLX(FR(347),FI(347))
       FR(348) = 2*Is12*(-Cij345R(1,1)+s12*Dij1345R(1,2)-s45*Dij1345R(
     -   1,2)-s12*Dij1345R(2,1)-s12*Dij1345R(2,2)+s34*Dij1345R(2,2)+s
     -   45*Dij1345R(2,2)-s45*Dij1345R(3,1)+s12*Dij1345R(4,2)-s34*Dij
     -   1345R(4,2)+3*(C0345R-s45*Dij1345R(4,2))+p3sq*(Dij1345R(1,2)-
     -   Dij1345R(2,2)+Dij1345R(4,2)-Dij1345R(5,2))+s34*Dij1345R(5,2)
     -   -s45*Dij1345R(6,2)+p5sq*(Dij1345R(3,1)+2*(Dij1345R(1,2)-Dij1
     -   345R(2,1)-Dij1345R(4,2)-Dij1345R(5,2))+Dij1345R(6,2))-2*(s34
     -   *Dij1345R(1,2)-s12*Dij1345R(1,3)-s45*Dij1345R(2,1)+s12*Dij13
     -   45R(4,3)+s34*Dij1345R(4,3)-s45*Dij1345R(5,2)-s34*Dij1345R(6,
     -   3))-4*(Dij1345R(7,2)-2*Dij1345R(11,3)+Dij1345R(12,3))+D01345
     -   R*P(49)+Dij1345R(1,1)*P(51))
       FI(348) = 2*Is12*(-Cij345I(1,1)+s12*Dij1345I(1,2)-s45*Dij1345I(
     -   1,2)-s12*Dij1345I(2,1)-s12*Dij1345I(2,2)+s34*Dij1345I(2,2)+s
     -   45*Dij1345I(2,2)-s45*Dij1345I(3,1)+s12*Dij1345I(4,2)-s34*Dij
     -   1345I(4,2)+3*(C0345I-s45*Dij1345I(4,2))+p3sq*(Dij1345I(1,2)-
     -   Dij1345I(2,2)+Dij1345I(4,2)-Dij1345I(5,2))+s34*Dij1345I(5,2)
     -   -s45*Dij1345I(6,2)+p5sq*(Dij1345I(3,1)+2*(Dij1345I(1,2)-Dij1
     -   345I(2,1)-Dij1345I(4,2)-Dij1345I(5,2))+Dij1345I(6,2))-2*(s34
     -   *Dij1345I(1,2)-s12*Dij1345I(1,3)-s45*Dij1345I(2,1)+s12*Dij13
     -   45I(4,3)+s34*Dij1345I(4,3)-s45*Dij1345I(5,2)-s34*Dij1345I(6,
     -   3))-4*(Dij1345I(7,2)-2*Dij1345I(11,3)+Dij1345I(12,3))+D01345
     -   I*P(49)+Dij1345I(1,1)*P(51))
       F(348)=DCMPLX(FR(348),FI(348))
       FR(349) = 2*Is12*(-Cij345R(2,1)+p3sq*Dij1345R(1,2)-s45*Dij1345R
     -   (1,2)-s12*Dij1345R(2,1)-s34*Dij1345R(2,1)+s45*Dij1345R(2,1)-
     -   p5sq*Dij1345R(3,1)+p5sq*Dij1345R(3,2)-s45*Dij1345R(3,2)+3*(C
     -   0345R+s12*Dij1345R(1,2)-s34*Dij1345R(4,2))-s12*Dij1345R(5,2)
     -   +s34*Dij1345R(5,2)+s45*Dij1345R(5,2)-p3sq*Dij1345R(6,2)-s12*
     -   Dij1345R(6,2)+s34*Dij1345R(6,2)+s45*Dij1345R(6,2)+2*(p5sq*Di
     -   j1345R(1,2)-s34*Dij1345R(1,2)+s12*Dij1345R(1,3)+s12*Dij1345R
     -   (4,2)-s45*Dij1345R(4,2)-s34*Dij1345R(4,3)-s12*Dij1345R(5,3)+
     -   Dij1345R(7,2)+s34*Dij1345R(10,3))+8*Dij1345R(11,3)-4*(p5sq*D
     -   ij1345R(5,2)+Dij1345R(13,3))+D01345R*P(53)+Dij1345R(1,1)*P(5
     -   4))
       FI(349) = 2*Is12*(-Cij345I(2,1)+p3sq*Dij1345I(1,2)-s45*Dij1345I
     -   (1,2)-s12*Dij1345I(2,1)-s34*Dij1345I(2,1)+s45*Dij1345I(2,1)-
     -   p5sq*Dij1345I(3,1)+p5sq*Dij1345I(3,2)-s45*Dij1345I(3,2)+3*(C
     -   0345I+s12*Dij1345I(1,2)-s34*Dij1345I(4,2))-s12*Dij1345I(5,2)
     -   +s34*Dij1345I(5,2)+s45*Dij1345I(5,2)-p3sq*Dij1345I(6,2)-s12*
     -   Dij1345I(6,2)+s34*Dij1345I(6,2)+s45*Dij1345I(6,2)+2*(p5sq*Di
     -   j1345I(1,2)-s34*Dij1345I(1,2)+s12*Dij1345I(1,3)+s12*Dij1345I
     -   (4,2)-s45*Dij1345I(4,2)-s34*Dij1345I(4,3)-s12*Dij1345I(5,3)+
     -   Dij1345I(7,2)+s34*Dij1345I(10,3))+8*Dij1345I(11,3)-4*(p5sq*D
     -   ij1345I(5,2)+Dij1345I(13,3))+D01345I*P(53)+Dij1345I(1,1)*P(5
     -   4))
       F(349)=DCMPLX(FR(349),FI(349))
       FR(350) = 2*Is12*(C0345R+p3sq*Dij1345R(1,2)-s45*Dij1345R(1,2)+s
     -   12*Dij1345R(2,1)-s34*Dij1345R(2,1)-s45*Dij1345R(2,1)-p5sq*Di
     -   j1345R(3,1)+s45*Dij1345R(3,1)+3*(s12*Dij1345R(1,2)-s34*Dij13
     -   45R(4,2))+2*(p5sq*Dij1345R(1,2)-s34*Dij1345R(1,2)+s12*Dij134
     -   5R(1,3)+s12*Dij1345R(4,2)-s45*Dij1345R(4,2)-s34*Dij1345R(4,3
     -   ))+(s34+2*s45)*Dij1345R(5,2)+6*Dij1345R(7,2)+8*Dij1345R(11,3
     -   )+D01345R*P(55)+Dij1345R(1,1)*P(57)-Dij1345R(5,2)*P(58))
       FI(350) = 2*Is12*(C0345I+p3sq*Dij1345I(1,2)-s45*Dij1345I(1,2)+s
     -   12*Dij1345I(2,1)-s34*Dij1345I(2,1)-s45*Dij1345I(2,1)-p5sq*Di
     -   j1345I(3,1)+s45*Dij1345I(3,1)+3*(s12*Dij1345I(1,2)-s34*Dij13
     -   45I(4,2))+2*(p5sq*Dij1345I(1,2)-s34*Dij1345I(1,2)+s12*Dij134
     -   5I(1,3)+s12*Dij1345I(4,2)-s45*Dij1345I(4,2)-s34*Dij1345I(4,3
     -   ))+(s34+2*s45)*Dij1345I(5,2)+6*Dij1345I(7,2)+8*Dij1345I(11,3
     -   )+D01345I*P(55)+Dij1345I(1,1)*P(57)-Dij1345I(5,2)*P(58))
       F(350)=DCMPLX(FR(350),FI(350))
       FR(351) = Dij1345R(1,2)*P(40)+C0345R*P(59)+D01345R*P(61)+Dij134
     -   5R(2,1)*P(65)+Dij1345R(5,2)*P(69)+Is12*(2*p3sq*s34*Dij1345R(
     -   2,2)+Dij1345R(1,1)*P(62)+Dij1345R(3,1)*P(16)*P(66)-Dij1345R(
     -   4,2)*P(68)+s34*Dij1345R(6,2)*P(70))-2*Dij1345R(7,2)*P(71)
       FI(351) = Dij1345I(1,2)*P(40)+C0345I*P(59)+D01345I*P(61)+Dij134
     -   5I(2,1)*P(65)+Dij1345I(5,2)*P(69)+Is12*(2*p3sq*s34*Dij1345I(
     -   2,2)+Dij1345I(1,1)*P(62)+Dij1345I(3,1)*P(16)*P(66)-Dij1345I(
     -   4,2)*P(68)+s34*Dij1345I(6,2)*P(70))-2*Dij1345I(7,2)*P(71)
       F(351)=DCMPLX(FR(351),FI(351))
       FR(352) = -(Dij1345R(1,2)*P(40))-Dij1345R(5,2)*P(69)+2*Dij1345R
     -   (7,2)*P(71)+C0345R*P(72)+D01345R*P(74)+Is12*(-2*p3sq*s34*Dij
     -   1345R(2,2)-Dij1345R(1,1)*P(62)+Dij1345R(4,2)*P(68)+s34*Dij13
     -   45R(6,2)*P(69)+Dij1345R(2,1)*P(77)+Dij1345R(3,1)*P(16)*P(78)
     -   )
       FI(352) = -(Dij1345I(1,2)*P(40))-Dij1345I(5,2)*P(69)+2*Dij1345I
     -   (7,2)*P(71)+C0345I*P(72)+D01345I*P(74)+Is12*(-2*p3sq*s34*Dij
     -   1345I(2,2)-Dij1345I(1,1)*P(62)+Dij1345I(4,2)*P(68)+s34*Dij13
     -   45I(6,2)*P(69)+Dij1345I(2,1)*P(77)+Dij1345I(3,1)*P(16)*P(78)
     -   )
       F(352)=DCMPLX(FR(352),FI(352))
       FR(353) = 2*Is12*(C0345R+Cij345R(1,1)-6*Dij1345R(7,2)-4*(Dij134
     -   5R(11,3)-Dij1345R(12,3))+2*(-Dij1345R(1,1)+Dij1345R(2,1))*P(
     -   8)-Dij1345R(1,2)*P(79)-Dij1345R(2,2)*P(80)+Dij1345R(4,2)*P(8
     -   1)+(-Dij1345R(5,2)+Dij1345R(6,2))*P(82))
       FI(353) = 2*Is12*(C0345I+Cij345I(1,1)-6*Dij1345I(7,2)-4*(Dij134
     -   5I(11,3)-Dij1345I(12,3))+2*(-Dij1345I(1,1)+Dij1345I(2,1))*P(
     -   8)-Dij1345I(1,2)*P(79)-Dij1345I(2,2)*P(80)+Dij1345I(4,2)*P(8
     -   1)+(-Dij1345I(5,2)+Dij1345I(6,2))*P(82))
       F(353)=DCMPLX(FR(353),FI(353))
       FR(354) = -2*Is12*(2*C0345R+Cij345R(1,1)+p3sq*Dij1345R(1,2)-3*s
     -   12*Dij1345R(1,2)-p3sq*Dij1345R(2,2)-s12*Dij1345R(2,2)-p5sq*D
     -   ij1345R(3,1)-(p5sq+2*(p3sq-s34))*Dij1345R(5,2)+s45*(-D01345R
     -   -Dij1345R(1,2)-2*(Dij1345R(1,1)-Dij1345R(2,1))+Dij1345R(2,2)
     -   +Dij1345R(3,1)+Dij1345R(5,2)-Dij1345R(6,2))+4*(s12*Dij1345R(
     -   4,2)-Dij1345R(7,2)+Dij1345R(11,3)-Dij1345R(12,3))-Dij1345R(1
     -   ,1)*P(83)-Dij1345R(2,1)*P(84)+Dij1345R(6,2)*P(85))
       FI(354) = -2*Is12*(2*C0345I+Cij345I(1,1)+p3sq*Dij1345I(1,2)-3*s
     -   12*Dij1345I(1,2)-p3sq*Dij1345I(2,2)-s12*Dij1345I(2,2)-p5sq*D
     -   ij1345I(3,1)-(p5sq+2*(p3sq-s34))*Dij1345I(5,2)+s45*(-D01345I
     -   -Dij1345I(1,2)-2*(Dij1345I(1,1)-Dij1345I(2,1))+Dij1345I(2,2)
     -   +Dij1345I(3,1)+Dij1345I(5,2)-Dij1345I(6,2))+4*(s12*Dij1345I(
     -   4,2)-Dij1345I(7,2)+Dij1345I(11,3)-Dij1345I(12,3))-Dij1345I(1
     -   ,1)*P(83)-Dij1345I(2,1)*P(84)+Dij1345I(6,2)*P(85))
       F(354)=DCMPLX(FR(354),FI(354))
       FR(355) = -2*Is12*(3*C0345R+Cij345R(2,1)-(3*s12+s45)*Dij1345R(1
     -   ,2)+p3sq*(Dij1345R(1,2)-2*Dij1345R(2,2)-Dij1345R(3,1))+(s34+
     -   2*s45)*Dij1345R(3,1)-4*(2*Dij1345R(7,2)-Dij1345R(11,3)+Dij13
     -   45R(13,3))+2*Dij1345R(4,2)*P(40)+Dij1345R(3,2)*P(43)-D01345R
     -   *P(86)-Dij1345R(5,2)*P(87)-Dij1345R(6,2)*P(88)+Dij1345R(1,1)
     -   *P(90)+Dij1345R(2,1)*P(91))
       FI(355) = -2*Is12*(3*C0345I+Cij345I(2,1)-(3*s12+s45)*Dij1345I(1
     -   ,2)+p3sq*(Dij1345I(1,2)-2*Dij1345I(2,2)-Dij1345I(3,1))+(s34+
     -   2*s45)*Dij1345I(3,1)-4*(2*Dij1345I(7,2)-Dij1345I(11,3)+Dij13
     -   45I(13,3))+2*Dij1345I(4,2)*P(40)+Dij1345I(3,2)*P(43)-D01345I
     -   *P(86)-Dij1345I(5,2)*P(87)-Dij1345I(6,2)*P(88)+Dij1345I(1,1)
     -   *P(90)+Dij1345I(2,1)*P(91))
       F(355)=DCMPLX(FR(355),FI(355))
       FR(356) = 2*Is12*(C0345R+Cij345R(2,1)-12*Dij1345R(7,2)-4*(Dij13
     -   45R(11,3)-Dij1345R(13,3))-2*(p3sq*Dij1345R(2,2)+Dij1345R(1,1
     -   )*P(8))+Dij1345R(3,1)*P(17)-Dij1345R(3,2)*P(52)-Dij1345R(1,2
     -   )*P(79)+Dij1345R(2,1)*P(92)+Dij1345R(4,2)*P(93)-Dij1345R(5,2
     -   )*P(94)+Dij1345R(6,2)*P(95))
       FI(356) = 2*Is12*(C0345I+Cij345I(2,1)-12*Dij1345I(7,2)-4*(Dij13
     -   45I(11,3)-Dij1345I(13,3))-2*(p3sq*Dij1345I(2,2)+Dij1345I(1,1
     -   )*P(8))+Dij1345I(3,1)*P(17)-Dij1345I(3,2)*P(52)-Dij1345I(1,2
     -   )*P(79)+Dij1345I(2,1)*P(92)+Dij1345I(4,2)*P(93)-Dij1345I(5,2
     -   )*P(94)+Dij1345I(6,2)*P(95))
       F(356)=DCMPLX(FR(356),FI(356))
       FR(357) = 2*Is12*(3*C0345R-2*(p3sq*Dij1345R(2,2)+8*Dij1345R(7,2
     -   )+2*Dij1345R(11,3)+Dij1345R(1,1)*P(10)-Dij1345R(3,1)*P(52))+
     -   Dij1345R(6,2)*P(69)-Dij1345R(1,2)*P(79)-Dij1345R(5,2)*P(82)+
     -   Dij1345R(4,2)*P(93)-Dij1345R(2,1)*P(96))
       FI(357) = 2*Is12*(3*C0345I-2*(p3sq*Dij1345I(2,2)+8*Dij1345I(7,2
     -   )+2*Dij1345I(11,3)+Dij1345I(1,1)*P(10)-Dij1345I(3,1)*P(52))+
     -   Dij1345I(6,2)*P(69)-Dij1345I(1,2)*P(79)-Dij1345I(5,2)*P(82)+
     -   Dij1345I(4,2)*P(93)-Dij1345I(2,1)*P(96))
       F(357)=DCMPLX(FR(357),FI(357))
       FR(358) = -2*Is12*(-(s12*(3*Dij1345R(1,2)+2*Dij1345R(2,1)))+p3s
     -   q*(Dij1345R(1,2)-2*Dij1345R(2,2))+p5sq*Dij1345R(3,1)-s45*(Di
     -   j1345R(1,2)-2*Dij1345R(2,1)+Dij1345R(3,1))-10*Dij1345R(7,2)+
     -   4*(C0345R+Dij1345R(11,3))+2*Dij1345R(4,2)*P(40)+Dij1345R(6,2
     -   )*P(69)-Dij1345R(5,2)*P(97)+D01345R*P(100)+Dij1345R(1,1)*P(1
     -   02))
       FI(358) = -2*Is12*(-(s12*(3*Dij1345I(1,2)+2*Dij1345I(2,1)))+p3s
     -   q*(Dij1345I(1,2)-2*Dij1345I(2,2))+p5sq*Dij1345I(3,1)-s45*(Di
     -   j1345I(1,2)-2*Dij1345I(2,1)+Dij1345I(3,1))-10*Dij1345I(7,2)+
     -   4*(C0345I+Dij1345I(11,3))+2*Dij1345I(4,2)*P(40)+Dij1345I(6,2
     -   )*P(69)-Dij1345I(5,2)*P(97)+D01345I*P(100)+Dij1345I(1,1)*P(1
     -   02))
       F(358)=DCMPLX(FR(358),FI(358))
       FR(359) = -2*Is12*(2*(C0345R+Dij1345R(7,2)+2*Dij1345R(11,3)+D01
     -   345R*P(10))-Dij1345R(1,2)*P(13)+Dij1345R(3,1)*P(43)+Dij1345R
     -   (5,2)*P(82)+Dij1345R(1,1)*P(103)-Dij1345R(2,1)*P(104)-Dij134
     -   5R(4,2)*P(105))
       FI(359) = -2*Is12*(2*(C0345I+Dij1345I(7,2)+2*Dij1345I(11,3)+D01
     -   345I*P(10))-Dij1345I(1,2)*P(13)+Dij1345I(3,1)*P(43)+Dij1345I
     -   (5,2)*P(82)+Dij1345I(1,1)*P(103)-Dij1345I(2,1)*P(104)-Dij134
     -   5I(4,2)*P(105))
       F(359)=DCMPLX(FR(359),FI(359))
       FR(360) = 2*Is12*(3*C0345R-Cij345R(1,1)+s12*(Dij1345R(1,2)-3*Di
     -   j1345R(2,1)+Dij1345R(2,2)-Dij1345R(4,2))+p5sq*(4*Dij1345R(3,
     -   1)+Dij1345R(5,2))-s45*(-Dij1345R(1,2)-2*Dij1345R(2,1)+Dij134
     -   5R(2,2)+4*Dij1345R(3,1)+Dij1345R(4,2)+Dij1345R(5,2))-s34*(Di
     -   j1345R(3,1)+2*Dij1345R(5,2))-p3sq*(Dij1345R(1,1)+Dij1345R(1,
     -   2)+Dij1345R(2,2)-Dij1345R(3,1)-Dij1345R(4,2)-2*Dij1345R(5,2)
     -   +Dij1345R(6,2))-10*Dij1345R(7,2)-4*(Dij1345R(11,3)-Dij1345R(
     -   12,3))+D01345R*P(106)+2*Dij1345R(1,1)*P(107)+Dij1345R(6,2)*P
     -   (108))
       FI(360) = 2*Is12*(3*C0345I-Cij345I(1,1)+s12*(Dij1345I(1,2)-3*Di
     -   j1345I(2,1)+Dij1345I(2,2)-Dij1345I(4,2))+p5sq*(4*Dij1345I(3,
     -   1)+Dij1345I(5,2))-s45*(-Dij1345I(1,2)-2*Dij1345I(2,1)+Dij134
     -   5I(2,2)+4*Dij1345I(3,1)+Dij1345I(4,2)+Dij1345I(5,2))-s34*(Di
     -   j1345I(3,1)+2*Dij1345I(5,2))-p3sq*(Dij1345I(1,1)+Dij1345I(1,
     -   2)+Dij1345I(2,2)-Dij1345I(3,1)-Dij1345I(4,2)-2*Dij1345I(5,2)
     -   +Dij1345I(6,2))-10*Dij1345I(7,2)-4*(Dij1345I(11,3)-Dij1345I(
     -   12,3))+D01345I*P(106)+2*Dij1345I(1,1)*P(107)+Dij1345I(6,2)*P
     -   (108))
       F(360)=DCMPLX(FR(360),FI(360))
       FR(361) = 2*Is12*(-Cij345R(2,1)+3*p5sq*Dij1345R(3,1)+2*(C0345R-
     -   s34*Dij1345R(3,1))-p5sq*Dij1345R(3,2)+s34*Dij1345R(3,2)+s12*
     -   (Dij1345R(1,2)-3*Dij1345R(2,1)+Dij1345R(3,1)-Dij1345R(4,2))+
     -   s45*(Dij1345R(1,2)-5*Dij1345R(3,1)+Dij1345R(3,2)+Dij1345R(4,
     -   2)+3*(Dij1345R(2,1)-Dij1345R(5,2)))-p3sq*(Dij1345R(1,2)+Dij1
     -   345R(2,1)-2*Dij1345R(3,1)+Dij1345R(3,2)+Dij1345R(4,2)-4*Dij1
     -   345R(5,2))+(p5sq-2*s34)*Dij1345R(5,2)-6*Dij1345R(7,2)-4*(Dij
     -   1345R(11,3)-Dij1345R(13,3))+D01345R*P(109)-Dij1345R(1,1)*P(1
     -   10)-Dij1345R(6,2)*P(111))
       FI(361) = 2*Is12*(-Cij345I(2,1)+3*p5sq*Dij1345I(3,1)+2*(C0345I-
     -   s34*Dij1345I(3,1))-p5sq*Dij1345I(3,2)+s34*Dij1345I(3,2)+s12*
     -   (Dij1345I(1,2)-3*Dij1345I(2,1)+Dij1345I(3,1)-Dij1345I(4,2))+
     -   s45*(Dij1345I(1,2)-5*Dij1345I(3,1)+Dij1345I(3,2)+Dij1345I(4,
     -   2)+3*(Dij1345I(2,1)-Dij1345I(5,2)))-p3sq*(Dij1345I(1,2)+Dij1
     -   345I(2,1)-2*Dij1345I(3,1)+Dij1345I(3,2)+Dij1345I(4,2)-4*Dij1
     -   345I(5,2))+(p5sq-2*s34)*Dij1345I(5,2)-6*Dij1345I(7,2)-4*(Dij
     -   1345I(11,3)-Dij1345I(13,3))+D01345I*P(109)-Dij1345I(1,1)*P(1
     -   10)-Dij1345I(6,2)*P(111))
       F(361)=DCMPLX(FR(361),FI(361))
       FR(362) = 2*Is12*(C0345R-4*(Dij1345R(7,2)+Dij1345R(11,3))-Dij13
     -   45R(4,2)*P(40)+Dij1345R(5,2)*P(97)-D01345R*P(112)+Dij1345R(1
     -   ,1)*P(113)-Dij1345R(2,1)*P(114)+Dij1345R(3,1)*P(115)+Dij1345
     -   R(1,2)*P(116))
       FI(362) = 2*Is12*(C0345I-4*(Dij1345I(7,2)+Dij1345I(11,3))-Dij13
     -   45I(4,2)*P(40)+Dij1345I(5,2)*P(97)-D01345I*P(112)+Dij1345I(1
     -   ,1)*P(113)-Dij1345I(2,1)*P(114)+Dij1345I(3,1)*P(115)+Dij1345
     -   I(1,2)*P(116))
       F(362)=DCMPLX(FR(362),FI(362))
       FR(363) = -2*Is12*(-Cij345R(1,1)-8*Dij1345R(7,2)+4*(C0345R+Dij1
     -   345R(11,3)-Dij1345R(12,3))+2*D01345R*P(8)-Dij1345R(1,2)*P(13
     -   )-Dij1345R(2,2)*P(41)+Dij1345R(6,2)*P(52)+Dij1345R(3,1)*P(66
     -   )+Dij1345R(5,2)*P(82)+Dij1345R(1,1)*P(117)-Dij1345R(2,1)*P(1
     -   19)-Dij1345R(4,2)*P(121))
       FI(363) = -2*Is12*(-Cij345I(1,1)-8*Dij1345I(7,2)+4*(C0345I+Dij1
     -   345I(11,3)-Dij1345I(12,3))+2*D01345I*P(8)-Dij1345I(1,2)*P(13
     -   )-Dij1345I(2,2)*P(41)+Dij1345I(6,2)*P(52)+Dij1345I(3,1)*P(66
     -   )+Dij1345I(5,2)*P(82)+Dij1345I(1,1)*P(117)-Dij1345I(2,1)*P(1
     -   19)-Dij1345I(4,2)*P(121))
       F(363)=DCMPLX(FR(363),FI(363))
       FR(364) = -2*Is12*(-Cij345R(2,1)+4*(C0345R+Dij1345R(11,3)-Dij13
     -   45R(13,3))-2*(Dij1345R(7,2)-D01345R*P(8))-Dij1345R(1,2)*P(13
     -   )-Dij1345R(6,2)*P(41)+Dij1345R(3,2)*P(52)-Dij1345R(4,2)*P(10
     -   5)+Dij1345R(1,1)*P(117)-Dij1345R(2,1)*P(122)+Dij1345R(3,1)*P
     -   (123)+Dij1345R(5,2)*P(125))
       FI(364) = -2*Is12*(-Cij345I(2,1)+4*(C0345I+Dij1345I(11,3)-Dij13
     -   45I(13,3))-2*(Dij1345I(7,2)-D01345I*P(8))-Dij1345I(1,2)*P(13
     -   )-Dij1345I(6,2)*P(41)+Dij1345I(3,2)*P(52)-Dij1345I(4,2)*P(10
     -   5)+Dij1345I(1,1)*P(117)-Dij1345I(2,1)*P(122)+Dij1345I(3,1)*P
     -   (123)+Dij1345I(5,2)*P(125))
       F(364)=DCMPLX(FR(364),FI(364))
       FR(365) = (Is12*(2*p3sq*s34*Dij1345R(2,2)+Dij1345R(1,1)*P(62)+D
     -   ij1345R(3,1)*P(16)*P(66)-Dij1345R(4,2)*P(68)+s34*Dij1345R(6,
     -   2)*P(70)+C0345R*P(126)+D01345R*P(128)+Dij1345R(2,1)*P(129)))
     -   /2.+(Dij1345R(1,2)*P(40)+Dij1345R(5,2)*P(69)+2*Dij1345R(7,2)
     -   *P(130))/2.
       FI(365) = (Is12*(2*p3sq*s34*Dij1345I(2,2)+Dij1345I(1,1)*P(62)+D
     -   ij1345I(3,1)*P(16)*P(66)-Dij1345I(4,2)*P(68)+s34*Dij1345I(6,
     -   2)*P(70)+C0345I*P(126)+D01345I*P(128)+Dij1345I(2,1)*P(129)))
     -   /2.+(Dij1345I(1,2)*P(40)+Dij1345I(5,2)*P(69)+2*Dij1345I(7,2)
     -   *P(130))/2.
       F(365)=DCMPLX(FR(365),FI(365))
       FR(366) = Is12*(2*C0345R+Cij345R(1,1)+p3sq*Dij1345R(1,2)-3*s12*
     -   Dij1345R(1,2)-p3sq*Dij1345R(2,2)-s12*Dij1345R(2,2)-p5sq*Dij1
     -   345R(3,1)-(p5sq+2*(p3sq-s34))*Dij1345R(5,2)+s45*(-D01345R-Di
     -   j1345R(1,2)-2*(Dij1345R(1,1)-Dij1345R(2,1))+Dij1345R(2,2)+Di
     -   j1345R(3,1)+Dij1345R(5,2)-Dij1345R(6,2))+4*(s12*Dij1345R(4,2
     -   )-Dij1345R(7,2)+Dij1345R(11,3)-Dij1345R(12,3))+Dij1345R(2,1)
     -   *P(83)+Dij1345R(1,1)*P(84)+Dij1345R(6,2)*P(85))
       FI(366) = Is12*(2*C0345I+Cij345I(1,1)+p3sq*Dij1345I(1,2)-3*s12*
     -   Dij1345I(1,2)-p3sq*Dij1345I(2,2)-s12*Dij1345I(2,2)-p5sq*Dij1
     -   345I(3,1)-(p5sq+2*(p3sq-s34))*Dij1345I(5,2)+s45*(-D01345I-Di
     -   j1345I(1,2)-2*(Dij1345I(1,1)-Dij1345I(2,1))+Dij1345I(2,2)+Di
     -   j1345I(3,1)+Dij1345I(5,2)-Dij1345I(6,2))+4*(s12*Dij1345I(4,2
     -   )-Dij1345I(7,2)+Dij1345I(11,3)-Dij1345I(12,3))+Dij1345I(2,1)
     -   *P(83)+Dij1345I(1,1)*P(84)+Dij1345I(6,2)*P(85))
       F(366)=DCMPLX(FR(366),FI(366))
       FR(367) = Is12*(3*C0345R+Cij345R(2,1)+p5sq*Dij1345R(3,2)-s45*(D
     -   ij1345R(1,2)-Dij1345R(2,1)+Dij1345R(3,2)-2*(Dij1345R(3,1)-Di
     -   j1345R(4,2))-3*Dij1345R(5,2))-p5sq*Dij1345R(5,2)+s34*(Dij134
     -   5R(3,1)-Dij1345R(3,2)+2*Dij1345R(5,2))-s12*(3*Dij1345R(1,2)+
     -   Dij1345R(2,1)-2*(Dij1345R(4,2)+Dij1345R(5,2)))+p3sq*(Dij1345
     -   R(1,1)+Dij1345R(1,2)-Dij1345R(3,1)+Dij1345R(3,2)-2*(Dij1345R
     -   (2,2)-Dij1345R(4,2)+2*Dij1345R(5,2)-Dij1345R(6,2)))-8*Dij134
     -   5R(7,2)+4*(Dij1345R(11,3)-Dij1345R(13,3))+D01345R*P(25)+2*Di
     -   j1345R(1,1)*P(89)-Dij1345R(6,2)*P(131))
       FI(367) = Is12*(3*C0345I+Cij345I(2,1)+p5sq*Dij1345I(3,2)-s45*(D
     -   ij1345I(1,2)-Dij1345I(2,1)+Dij1345I(3,2)-2*(Dij1345I(3,1)-Di
     -   j1345I(4,2))-3*Dij1345I(5,2))-p5sq*Dij1345I(5,2)+s34*(Dij134
     -   5I(3,1)-Dij1345I(3,2)+2*Dij1345I(5,2))-s12*(3*Dij1345I(1,2)+
     -   Dij1345I(2,1)-2*(Dij1345I(4,2)+Dij1345I(5,2)))+p3sq*(Dij1345
     -   I(1,1)+Dij1345I(1,2)-Dij1345I(3,1)+Dij1345I(3,2)-2*(Dij1345I
     -   (2,2)-Dij1345I(4,2)+2*Dij1345I(5,2)-Dij1345I(6,2)))-8*Dij134
     -   5I(7,2)+4*(Dij1345I(11,3)-Dij1345I(13,3))+D01345I*P(25)+2*Di
     -   j1345I(1,1)*P(89)-Dij1345I(6,2)*P(131))
       F(367)=DCMPLX(FR(367),FI(367))
       FR(368) = Is12*(p5sq*Dij1345R(3,1)+s12*(-3*Dij1345R(1,2)-2*(Dij
     -   1345R(2,1)-Dij1345R(4,2)))-s45*(Dij1345R(1,2)+Dij1345R(3,1)-
     -   2*(Dij1345R(2,1)-Dij1345R(4,2))-Dij1345R(5,2))-(p5sq-2*s34)*
     -   Dij1345R(5,2)-s34*Dij1345R(6,2)+p3sq*(Dij1345R(1,2)-2*(Dij13
     -   45R(2,2)-Dij1345R(4,2)+Dij1345R(5,2))+Dij1345R(6,2))-10*Dij1
     -   345R(7,2)+4*(C0345R+Dij1345R(11,3))+D01345R*P(100)+Dij1345R(
     -   1,1)*P(102))
       FI(368) = Is12*(p5sq*Dij1345I(3,1)+s12*(-3*Dij1345I(1,2)-2*(Dij
     -   1345I(2,1)-Dij1345I(4,2)))-s45*(Dij1345I(1,2)+Dij1345I(3,1)-
     -   2*(Dij1345I(2,1)-Dij1345I(4,2))-Dij1345I(5,2))-(p5sq-2*s34)*
     -   Dij1345I(5,2)-s34*Dij1345I(6,2)+p3sq*(Dij1345I(1,2)-2*(Dij13
     -   45I(2,2)-Dij1345I(4,2)+Dij1345I(5,2))+Dij1345I(6,2))-10*Dij1
     -   345I(7,2)+4*(C0345I+Dij1345I(11,3))+D01345I*P(100)+Dij1345I(
     -   1,1)*P(102))
       F(368)=DCMPLX(FR(368),FI(368))
       FR(369) = Is12*(-C0345R-Cij345R(1,1)+6*Dij1345R(7,2)+4*(Dij1345
     -   R(11,3)-Dij1345R(12,3))+2*(Dij1345R(1,1)-Dij1345R(2,1))*P(8)
     -   +Dij1345R(1,2)*P(79)+Dij1345R(2,2)*P(80)-Dij1345R(4,2)*P(81)
     -   +Dij1345R(5,2)*P(82)+Dij1345R(6,2)*P(132))
       FI(369) = Is12*(-C0345I-Cij345I(1,1)+6*Dij1345I(7,2)+4*(Dij1345
     -   I(11,3)-Dij1345I(12,3))+2*(Dij1345I(1,1)-Dij1345I(2,1))*P(8)
     -   +Dij1345I(1,2)*P(79)+Dij1345I(2,2)*P(80)-Dij1345I(4,2)*P(81)
     -   +Dij1345I(5,2)*P(82)+Dij1345I(6,2)*P(132))
       F(369)=DCMPLX(FR(369),FI(369))
       FR(370) = Is12*(-C0345R-Cij345R(2,1)+12*Dij1345R(7,2)+4*(Dij134
     -   5R(11,3)-Dij1345R(13,3))+2*(p3sq*Dij1345R(2,2)+Dij1345R(1,1)
     -   *P(8))+Dij1345R(3,1)*P(19)+Dij1345R(3,2)*P(52)+Dij1345R(1,2)
     -   *P(79)-Dij1345R(2,1)*P(92)-Dij1345R(4,2)*P(93)+Dij1345R(5,2)
     -   *P(94)+Dij1345R(6,2)*P(133))
       FI(370) = Is12*(-C0345I-Cij345I(2,1)+12*Dij1345I(7,2)+4*(Dij134
     -   5I(11,3)-Dij1345I(13,3))+2*(p3sq*Dij1345I(2,2)+Dij1345I(1,1)
     -   *P(8))+Dij1345I(3,1)*P(19)+Dij1345I(3,2)*P(52)+Dij1345I(1,2)
     -   *P(79)-Dij1345I(2,1)*P(92)-Dij1345I(4,2)*P(93)+Dij1345I(5,2)
     -   *P(94)+Dij1345I(6,2)*P(133))
       F(370)=DCMPLX(FR(370),FI(370))
       FR(371) = Is12*(-3*C0345R+16*Dij1345R(7,2)+4*Dij1345R(11,3)+2*(
     -   p3sq*Dij1345R(2,2)+Dij1345R(1,1)*P(10)-Dij1345R(3,1)*P(52))+
     -   Dij1345R(6,2)*P(70)+Dij1345R(1,2)*P(79)+Dij1345R(5,2)*P(82)-
     -   Dij1345R(4,2)*P(93)+Dij1345R(2,1)*P(96))
       FI(371) = Is12*(-3*C0345I+16*Dij1345I(7,2)+4*Dij1345I(11,3)+2*(
     -   p3sq*Dij1345I(2,2)+Dij1345I(1,1)*P(10)-Dij1345I(3,1)*P(52))+
     -   Dij1345I(6,2)*P(70)+Dij1345I(1,2)*P(79)+Dij1345I(5,2)*P(82)-
     -   Dij1345I(4,2)*P(93)+Dij1345I(2,1)*P(96))
       F(371)=DCMPLX(FR(371),FI(371))
       FR(372) = 4*Is12*(Dij1345R(2,2)-Dij1345R(2,3)-Dij1345R(4,2)-Dij
     -   1345R(4,3)+2*Dij1345R(6,3))
       FI(372) = 4*Is12*(Dij1345I(2,2)-Dij1345I(2,3)-Dij1345I(4,2)-Dij
     -   1345I(4,3)+2*Dij1345I(6,3))
       F(372)=DCMPLX(FR(372),FI(372))
       FR(373) = 4*Is12*(Dij1345R(2,2)-Dij1345R(4,2)-Dij1345R(4,3)-Dij
     -   1345R(5,2)+Dij1345R(6,2)+Dij1345R(6,3)-Dij1345R(8,3)+Dij1345
     -   R(10,3))
       FI(373) = 4*Is12*(Dij1345I(2,2)-Dij1345I(4,2)-Dij1345I(4,3)-Dij
     -   1345I(5,2)+Dij1345I(6,2)+Dij1345I(6,3)-Dij1345I(8,3)+Dij1345
     -   I(10,3))
       F(373)=DCMPLX(FR(373),FI(373))
       FR(374) = 4*Is12*(Dij1345R(2,2)-Dij1345R(4,2)-Dij1345R(4,3)+Dij
     -   1345R(6,3))
       FI(374) = 4*Is12*(Dij1345I(2,2)-Dij1345I(4,2)-Dij1345I(4,3)+Dij
     -   1345I(6,3))
       F(374)=DCMPLX(FR(374),FI(374))
       FR(375) = -4*Is12*(Dij1345R(2,1)-Dij1345R(3,1)+3*Dij1345R(4,2)+
     -   Dij1345R(4,3)-2*(Dij1345R(2,2)+Dij1345R(5,2))+Dij1345R(6,2)-
     -   Dij1345R(6,3)+Dij1345R(8,3)-Dij1345R(10,3))
       FI(375) = -4*Is12*(Dij1345I(2,1)-Dij1345I(3,1)+3*Dij1345I(4,2)+
     -   Dij1345I(4,3)-2*(Dij1345I(2,2)+Dij1345I(5,2))+Dij1345I(6,2)-
     -   Dij1345I(6,3)+Dij1345I(8,3)-Dij1345I(10,3))
       F(375)=DCMPLX(FR(375),FI(375))
       FR(376) = -4*Is12*(Dij1345R(3,2)+Dij1345R(4,3)-Dij1345R(5,2)+3*
     -   (Dij1345R(4,2)-Dij1345R(6,2))+Dij1345R(9,3)+2*(Dij1345R(2,1)
     -   -Dij1345R(3,1)-Dij1345R(10,3)))
       FI(376) = -4*Is12*(Dij1345I(3,2)+Dij1345I(4,3)-Dij1345I(5,2)+3*
     -   (Dij1345I(4,2)-Dij1345I(6,2))+Dij1345I(9,3)+2*(Dij1345I(2,1)
     -   -Dij1345I(3,1)-Dij1345I(10,3)))
       F(376)=DCMPLX(FR(376),FI(376))
       FR(377) = 4*Is12*(-3*Dij1345R(4,2)-Dij1345R(4,3)-2*(Dij1345R(2,
     -   1)-Dij1345R(3,1)-Dij1345R(5,2))+Dij1345R(6,2)+Dij1345R(10,3)
     -   )
       FI(377) = 4*Is12*(-3*Dij1345I(4,2)-Dij1345I(4,3)-2*(Dij1345I(2,
     -   1)-Dij1345I(3,1)-Dij1345I(5,2))+Dij1345I(6,2)+Dij1345I(10,3)
     -   )
       F(377)=DCMPLX(FR(377),FI(377))
       FR(378) = -2*Is12*(D01345R+Dij1345R(1,1)+5*Dij1345R(4,2)+3*(Dij
     -   1345R(2,1)-Dij1345R(6,2))+2*(Dij1345R(4,3)+Dij1345R(5,2)-Dij
     -   1345R(10,3)))
       FI(378) = -2*Is12*(D01345I+Dij1345I(1,1)+5*Dij1345I(4,2)+3*(Dij
     -   1345I(2,1)-Dij1345I(6,2))+2*(Dij1345I(4,3)+Dij1345I(5,2)-Dij
     -   1345I(10,3)))
       F(378)=DCMPLX(FR(378),FI(378))
       FR(379) = -2*Is12*(D01345R+Dij1345R(1,1)-3*Dij1345R(2,2)+5*Dij1
     -   345R(4,2)+2*(Dij1345R(2,1)+Dij1345R(4,3)-Dij1345R(6,3)))
       FI(379) = -2*Is12*(D01345I+Dij1345I(1,1)-3*Dij1345I(2,2)+5*Dij1
     -   345I(4,2)+2*(Dij1345I(2,1)+Dij1345I(4,3)-Dij1345I(6,3)))
       F(379)=DCMPLX(FR(379),FI(379))
       FR(380) = Is12*(3*C0345R+Cij345R(1,1)+s45*(3*Dij1345R(2,1)+Dij1
     -   345R(2,2))-p3sq*(Dij1345R(2,1)+Dij1345R(2,2)-Dij1345R(3,1))-
     -   (s34-3*(p5sq-s45))*Dij1345R(3,1)-s12*(Dij1345R(2,2)+2*(Dij13
     -   45R(2,1)-Dij1345R(4,2)))-4*(2*Dij1345R(7,2)+Dij1345R(12,3))+
     -   2*D01345R*P(8)+Dij1345R(6,2)*P(43)+Dij1345R(1,1)*P(134))
       FI(380) = Is12*(3*C0345I+Cij345I(1,1)+s45*(3*Dij1345I(2,1)+Dij1
     -   345I(2,2))-p3sq*(Dij1345I(2,1)+Dij1345I(2,2)-Dij1345I(3,1))-
     -   (s34-3*(p5sq-s45))*Dij1345I(3,1)-s12*(Dij1345I(2,2)+2*(Dij13
     -   45I(2,1)-Dij1345I(4,2)))-4*(2*Dij1345I(7,2)+Dij1345I(12,3))+
     -   2*D01345I*P(8)+Dij1345I(6,2)*P(43)+Dij1345I(1,1)*P(134))
       F(380)=DCMPLX(FR(380),FI(380))
       FR(381) = -2*Is12*(D01345R+Dij1345R(1,1)+3*Dij1345R(2,1)+5*Dij1
     -   345R(4,2)+2*Dij1345R(4,3))
       FI(381) = -2*Is12*(D01345I+Dij1345I(1,1)+3*Dij1345I(2,1)+5*Dij1
     -   345I(4,2)+2*Dij1345I(4,3))
       F(381)=DCMPLX(FR(381),FI(381))
       FR(382) = Is12*(-5*C0345R+14*Dij1345R(7,2)+2*(p3sq*Dij1345R(2,2
     -   )-D01345R*P(8))-Dij1345R(4,2)*P(40)+Dij1345R(6,2)*P(70)+Dij1
     -   345R(3,1)*P(78)-Dij1345R(1,1)*P(135)+Dij1345R(2,1)*P(136))
       FI(382) = Is12*(-5*C0345I+14*Dij1345I(7,2)+2*(p3sq*Dij1345I(2,2
     -   )-D01345I*P(8))-Dij1345I(4,2)*P(40)+Dij1345I(6,2)*P(70)+Dij1
     -   345I(3,1)*P(78)-Dij1345I(1,1)*P(135)+Dij1345I(2,1)*P(136))
       F(382)=DCMPLX(FR(382),FI(382))
       FR(383) = Is12*(5*C0345R-14*Dij1345R(7,2)-2*(p3sq*Dij1345R(2,2)
     -   -D01345R*P(8))+Dij1345R(4,2)*P(40)+Dij1345R(3,1)*P(66)+Dij13
     -   45R(6,2)*P(69)+Dij1345R(1,1)*P(135)-Dij1345R(2,1)*P(136))
       FI(383) = Is12*(5*C0345I-14*Dij1345I(7,2)-2*(p3sq*Dij1345I(2,2)
     -   -D01345I*P(8))+Dij1345I(4,2)*P(40)+Dij1345I(3,1)*P(66)+Dij13
     -   45I(6,2)*P(69)+Dij1345I(1,1)*P(135)-Dij1345I(2,1)*P(136))
       F(383)=DCMPLX(FR(383),FI(383))
       FR(384) = Is12*(-C0345R+4*(Dij1345R(7,2)+Dij1345R(11,3))+Dij134
     -   5R(4,2)*P(40)-Dij1345R(5,2)*P(97)+D01345R*P(112)+Dij1345R(2,
     -   1)*P(114)-Dij1345R(1,2)*P(116)+Dij1345R(1,1)*P(137)+Dij1345R
     -   (3,1)*P(138))
       FI(384) = Is12*(-C0345I+4*(Dij1345I(7,2)+Dij1345I(11,3))+Dij134
     -   5I(4,2)*P(40)-Dij1345I(5,2)*P(97)+D01345I*P(112)+Dij1345I(2,
     -   1)*P(114)-Dij1345I(1,2)*P(116)+Dij1345I(1,1)*P(137)+Dij1345I
     -   (3,1)*P(138))
       F(384)=DCMPLX(FR(384),FI(384))
       FR(385) = Is12*(-3*C0345R+Cij345R(1,1)+s34*Dij1345R(3,1)-s12*(D
     -   ij1345R(1,2)-3*Dij1345R(2,1)+Dij1345R(2,2)-Dij1345R(4,2))-p5
     -   sq*(4*Dij1345R(3,1)+Dij1345R(5,2))+s45*(-Dij1345R(1,2)-2*Dij
     -   1345R(2,1)+Dij1345R(2,2)+4*Dij1345R(3,1)+Dij1345R(4,2)+Dij13
     -   45R(5,2))+p3sq*(Dij1345R(1,1)+Dij1345R(1,2)+Dij1345R(2,2)-Di
     -   j1345R(3,1)-Dij1345R(4,2)-2*Dij1345R(5,2)+Dij1345R(6,2))+10*
     -   Dij1345R(7,2)+4*(Dij1345R(11,3)-Dij1345R(12,3))-D01345R*P(10
     -   6)+2*(s34*Dij1345R(5,2)-Dij1345R(1,1)*P(107))-Dij1345R(6,2)*
     -   P(108))
       FI(385) = Is12*(-3*C0345I+Cij345I(1,1)+s34*Dij1345I(3,1)-s12*(D
     -   ij1345I(1,2)-3*Dij1345I(2,1)+Dij1345I(2,2)-Dij1345I(4,2))-p5
     -   sq*(4*Dij1345I(3,1)+Dij1345I(5,2))+s45*(-Dij1345I(1,2)-2*Dij
     -   1345I(2,1)+Dij1345I(2,2)+4*Dij1345I(3,1)+Dij1345I(4,2)+Dij13
     -   45I(5,2))+p3sq*(Dij1345I(1,1)+Dij1345I(1,2)+Dij1345I(2,2)-Di
     -   j1345I(3,1)-Dij1345I(4,2)-2*Dij1345I(5,2)+Dij1345I(6,2))+10*
     -   Dij1345I(7,2)+4*(Dij1345I(11,3)-Dij1345I(12,3))-D01345I*P(10
     -   6)+2*(s34*Dij1345I(5,2)-Dij1345I(1,1)*P(107))-Dij1345I(6,2)*
     -   P(108))
       F(385)=DCMPLX(FR(385),FI(385))
       FR(386) = Is12*(-2*C0345R+Cij345R(2,1)+6*Dij1345R(7,2)+4*(Dij13
     -   45R(11,3)-Dij1345R(13,3))+Dij1345R(4,2)*P(40)+Dij1345R(3,2)*
     -   P(43)-D01345R*P(109)+Dij1345R(1,1)*P(110)+Dij1345R(6,2)*P(11
     -   1)-Dij1345R(1,2)*P(116)+Dij1345R(2,1)*P(139)-Dij1345R(3,1)*P
     -   (140)-Dij1345R(5,2)*P(141))
       FI(386) = Is12*(-2*C0345I+Cij345I(2,1)+6*Dij1345I(7,2)+4*(Dij13
     -   45I(11,3)-Dij1345I(13,3))+Dij1345I(4,2)*P(40)+Dij1345I(3,2)*
     -   P(43)-D01345I*P(109)+Dij1345I(1,1)*P(110)+Dij1345I(6,2)*P(11
     -   1)-Dij1345I(1,2)*P(116)+Dij1345I(2,1)*P(139)-Dij1345I(3,1)*P
     -   (140)-Dij1345I(5,2)*P(141))
       F(386)=DCMPLX(FR(386),FI(386))
       FR(387) = Is12*(2*(C0345R+Dij1345R(7,2)+2*Dij1345R(11,3)+D01345
     -   R*P(10))-Dij1345R(1,2)*P(13)+Dij1345R(3,1)*P(43)+Dij1345R(5,
     -   2)*P(82)+Dij1345R(1,1)*P(103)-Dij1345R(2,1)*P(104)-Dij1345R(
     -   4,2)*P(105))
       FI(387) = Is12*(2*(C0345I+Dij1345I(7,2)+2*Dij1345I(11,3)+D01345
     -   I*P(10))-Dij1345I(1,2)*P(13)+Dij1345I(3,1)*P(43)+Dij1345I(5,
     -   2)*P(82)+Dij1345I(1,1)*P(103)-Dij1345I(2,1)*P(104)-Dij1345I(
     -   4,2)*P(105))
       F(387)=DCMPLX(FR(387),FI(387))
       FR(388) = Is12*(-Cij345R(1,1)-8*Dij1345R(7,2)+4*(C0345R+Dij1345
     -   R(11,3)-Dij1345R(12,3))+2*D01345R*P(8)-Dij1345R(1,2)*P(13)-D
     -   ij1345R(2,2)*P(41)+Dij1345R(6,2)*P(52)+Dij1345R(3,1)*P(66)+D
     -   ij1345R(5,2)*P(82)+Dij1345R(1,1)*P(117)-Dij1345R(2,1)*P(119)
     -   +Dij1345R(4,2)*P(142))
       FI(388) = Is12*(-Cij345I(1,1)-8*Dij1345I(7,2)+4*(C0345I+Dij1345
     -   I(11,3)-Dij1345I(12,3))+2*D01345I*P(8)-Dij1345I(1,2)*P(13)-D
     -   ij1345I(2,2)*P(41)+Dij1345I(6,2)*P(52)+Dij1345I(3,1)*P(66)+D
     -   ij1345I(5,2)*P(82)+Dij1345I(1,1)*P(117)-Dij1345I(2,1)*P(119)
     -   +Dij1345I(4,2)*P(142))
       F(388)=DCMPLX(FR(388),FI(388))
       FR(389) = Is12*(-Cij345R(2,1)+4*(C0345R+Dij1345R(11,3)-Dij1345R
     -   (13,3))-2*(Dij1345R(7,2)-D01345R*P(8))-Dij1345R(1,2)*P(13)-D
     -   ij1345R(6,2)*P(41)+Dij1345R(3,2)*P(52)-Dij1345R(4,2)*P(105)+
     -   Dij1345R(1,1)*P(117)-Dij1345R(2,1)*P(122)+Dij1345R(3,1)*P(12
     -   3)+Dij1345R(5,2)*P(125))
       FI(389) = Is12*(-Cij345I(2,1)+4*(C0345I+Dij1345I(11,3)-Dij1345I
     -   (13,3))-2*(Dij1345I(7,2)-D01345I*P(8))-Dij1345I(1,2)*P(13)-D
     -   ij1345I(6,2)*P(41)+Dij1345I(3,2)*P(52)-Dij1345I(4,2)*P(105)+
     -   Dij1345I(1,1)*P(117)-Dij1345I(2,1)*P(122)+Dij1345I(3,1)*P(12
     -   3)+Dij1345I(5,2)*P(125))
       F(389)=DCMPLX(FR(389),FI(389))
       FR(390) = (Is12*(-5*C0345R+14*Dij1345R(7,2)+2*(p3sq*Dij1345R(2,
     -   2)-D01345R*P(8))-Dij1345R(4,2)*P(40)+Dij1345R(6,2)*P(70)+Dij
     -   1345R(3,1)*P(78)-Dij1345R(1,1)*P(135)+Dij1345R(2,1)*P(136)))
     -   /2.
       FI(390) = (Is12*(-5*C0345I+14*Dij1345I(7,2)+2*(p3sq*Dij1345I(2,
     -   2)-D01345I*P(8))-Dij1345I(4,2)*P(40)+Dij1345I(6,2)*P(70)+Dij
     -   1345I(3,1)*P(78)-Dij1345I(1,1)*P(135)+Dij1345I(2,1)*P(136)))
     -   /2.
       F(390)=DCMPLX(FR(390),FI(390))
       FR(391) = 2*Is15s34*(-2*B025R+s15*(Cij125R(2,1)-2*(Cij125R(2,2)
     -   -Cij125R(3,2)))-4*Cij125R(4,2)+(C0125R+Cij125R(1,1))*P(20))
       FI(391) = 2*Is15s34*(-2*B025I+s15*(Cij125I(2,1)-2*(Cij125I(2,2)
     -   -Cij125I(3,2)))-4*Cij125I(4,2)+(C0125I+Cij125I(1,1))*P(20))
       F(391)=DCMPLX(FR(391),FI(391))
       FR(392) = -2*Is34*(Cij125R(2,1)+2*Cij125R(2,2))
       FI(392) = -2*Is34*(Cij125I(2,1)+2*Cij125I(2,2))
       F(392)=DCMPLX(FR(392),FI(392))
       FR(393) = -2*Is15s34*(s15*Cij125R(2,1)+4*Cij125R(4,2)+2*(B025R+
     -   (C0125R+Cij125R(1,1))*P(21)))
       FI(393) = -2*Is15s34*(s15*Cij125I(2,1)+4*Cij125I(4,2)+2*(B025I+
     -   (C0125I+Cij125I(1,1))*P(21)))
       F(393)=DCMPLX(FR(393),FI(393))
       FR(394) = 2*Is15s34*(2*B025R-s15*(Cij125R(2,1)-2*(Cij125R(2,2)-
     -   Cij125R(3,2)))+4*Cij125R(4,2)+(C0125R+Cij125R(1,1))*P(22))
       FI(394) = 2*Is15s34*(2*B025I-s15*(Cij125I(2,1)-2*(Cij125I(2,2)-
     -   Cij125I(3,2)))+4*Cij125I(4,2)+(C0125I+Cij125I(1,1))*P(22))
       F(394)=DCMPLX(FR(394),FI(394))
       FR(395) = 2*Is34*(Cij125R(2,1)+2*Cij125R(2,2))
       FI(395) = 2*Is34*(Cij125I(2,1)+2*Cij125I(2,2))
       F(395)=DCMPLX(FR(395),FI(395))
       FR(396) = Is15s34*(-2*B025R+s15*(Cij125R(2,1)-2*(Cij125R(2,2)-C
     -   ij125R(3,2)))-4*Cij125R(4,2)+(C0125R+Cij125R(1,1))*P(20))
       FI(396) = Is15s34*(-2*B025I+s15*(Cij125I(2,1)-2*(Cij125I(2,2)-C
     -   ij125I(3,2)))-4*Cij125I(4,2)+(C0125I+Cij125I(1,1))*P(20))
       F(396)=DCMPLX(FR(396),FI(396))
       FR(397) = Is34*(-Cij125R(2,1)-2*Cij125R(2,2))
       FI(397) = Is34*(-Cij125I(2,1)-2*Cij125I(2,2))
       F(397)=DCMPLX(FR(397),FI(397))
       FR(398) = Is15s34*(-(s15*Cij125R(2,1))-4*Cij125R(4,2)-2*(B025R+
     -   (C0125R+Cij125R(1,1))*P(21)))
       FI(398) = Is15s34*(-(s15*Cij125I(2,1))-4*Cij125I(4,2)-2*(B025I+
     -   (C0125I+Cij125I(1,1))*P(21)))
       F(398)=DCMPLX(FR(398),FI(398))
       FR(399) = Is15s34*(s15*Cij125R(2,1)+4*Cij125R(4,2)+2*(B025R+(C0
     -   125R+Cij125R(1,1))*P(21)))
       FI(399) = Is15s34*(s15*Cij125I(2,1)+4*Cij125I(4,2)+2*(B025I+(C0
     -   125I+Cij125I(1,1))*P(21)))
       F(399)=DCMPLX(FR(399),FI(399))
       FR(400) = (Is15s34*(-(s15*Cij125R(2,1))-4*Cij125R(4,2)-2*(B025R
     -   +(C0125R+Cij125R(1,1))*P(21))))/2.
       FI(400) = (Is15s34*(-(s15*Cij125I(2,1))-4*Cij125I(4,2)-2*(B025I
     -   +(C0125I+Cij125I(1,1))*P(21))))/2.
       F(400)=DCMPLX(FR(400),FI(400))
       FR(401) = 2*Is12*(C0135R+Cij135R(1,1)-2*(Cij135R(2,1)-Cij135R(2
     -   ,2)+Cij135R(3,2)))
       FI(401) = 2*Is12*(C0135I+Cij135I(1,1)-2*(Cij135I(2,1)-Cij135I(2
     -   ,2)+Cij135I(3,2)))
       F(401)=DCMPLX(FR(401),FI(401))
       FR(402) = 2*Is12*Is34*(-3*s12*Cij135R(2,1)-s34*Cij135R(2,1)+4*C
     -   ij135R(4,2)+C0135R*P(1)+Cij135R(1,1)*P(2)+2*(B035R+s12*Cij13
     -   5R(1,2)+s34*Cij135R(2,2)-Cij135R(3,2)*P(3)))
       FI(402) = 2*Is12*Is34*(-3*s12*Cij135I(2,1)-s34*Cij135I(2,1)+4*C
     -   ij135I(4,2)+C0135I*P(1)+Cij135I(1,1)*P(2)+2*(B035I+s12*Cij13
     -   5I(1,2)+s34*Cij135I(2,2)-Cij135I(3,2)*P(3)))
       F(402)=DCMPLX(FR(402),FI(402))
       FR(403) = 2*Is12*(C0135R+Cij135R(1,1)-3*Cij135R(2,1)-2*Cij135R(
     -   3,2))
       FI(403) = 2*Is12*(C0135I+Cij135I(1,1)-3*Cij135I(2,1)-2*Cij135I(
     -   3,2))
       F(403)=DCMPLX(FR(403),FI(403))
       FR(404) = 2*Is12*Is34*(2*(B035R+s12*Cij135R(1,2)-s34*(Cij135R(2
     -   ,1)+Cij135R(3,2))+2*Cij135R(4,2))+C0135R*P(4)+Cij135R(1,1)*P
     -   (5))
       FI(404) = 2*Is12*Is34*(2*(B035I+s12*Cij135I(1,2)-s34*(Cij135I(2
     -   ,1)+Cij135I(3,2))+2*Cij135I(4,2))+C0135I*P(4)+Cij135I(1,1)*P
     -   (5))
       F(404)=DCMPLX(FR(404),FI(404))
       FR(405) = Is12*(Cij135R(2,1)*P(7)-2*(B035R+2*Cij135R(4,2)+Cij13
     -   5R(1,1)*P(6)+C0135R*P(8)))
       FI(405) = Is12*(Cij135I(2,1)*P(7)-2*(B035I+2*Cij135I(4,2)+Cij13
     -   5I(1,1)*P(6)+C0135I*P(8)))
       F(405)=DCMPLX(FR(405),FI(405))
       FR(406) = Is12*(2*(B035R+2*Cij135R(4,2)+Cij135R(1,1)*P(6)+C0135
     -   R*P(8))+Cij135R(2,1)*P(9))
       FI(406) = Is12*(2*(B035I+2*Cij135I(4,2)+Cij135I(1,1)*P(6)+C0135
     -   I*P(8))+Cij135I(2,1)*P(9))
       F(406)=DCMPLX(FR(406),FI(406))
       FR(407) = -2*Is12*Is34*(-3*s12*Cij135R(2,1)-s34*Cij135R(2,1)+4*
     -   Cij135R(4,2)+C0135R*P(1)+Cij135R(1,1)*P(2)+2*(B035R+s12*Cij1
     -   35R(1,2)+s34*Cij135R(2,2)-Cij135R(3,2)*P(3)))
       FI(407) = -2*Is12*Is34*(-3*s12*Cij135I(2,1)-s34*Cij135I(2,1)+4*
     -   Cij135I(4,2)+C0135I*P(1)+Cij135I(1,1)*P(2)+2*(B035I+s12*Cij1
     -   35I(1,2)+s34*Cij135I(2,2)-Cij135I(3,2)*P(3)))
       F(407)=DCMPLX(FR(407),FI(407))
       FR(408) = 2*Is12*Is34*(2*B035R+s34*Cij135R(2,1)+4*Cij135R(4,2)+
     -   C0135R*P(11)+Cij135R(1,1)*P(12))
       FI(408) = 2*Is12*Is34*(2*B035I+s34*Cij135I(2,1)+4*Cij135I(4,2)+
     -   C0135I*P(11)+Cij135I(1,1)*P(12))
       F(408)=DCMPLX(FR(408),FI(408))
       FR(409) = -2*Is12*Is34*(2*(B035R+s12*Cij135R(1,2)-s34*(Cij135R(
     -   2,1)+Cij135R(3,2))+2*Cij135R(4,2))+C0135R*P(4)+Cij135R(1,1)*
     -   P(5))
       FI(409) = -2*Is12*Is34*(2*(B035I+s12*Cij135I(1,2)-s34*(Cij135I(
     -   2,1)+Cij135I(3,2))+2*Cij135I(4,2))+C0135I*P(4)+Cij135I(1,1)*
     -   P(5))
       F(409)=DCMPLX(FR(409),FI(409))
       FR(410) = -2*Is12*(C0135R+Cij135R(1,1)-2*(Cij135R(2,1)-Cij135R(
     -   2,2)+Cij135R(3,2)))
       FI(410) = -2*Is12*(C0135I+Cij135I(1,1)-2*(Cij135I(2,1)-Cij135I(
     -   2,2)+Cij135I(3,2)))
       F(410)=DCMPLX(FR(410),FI(410))
       FR(411) = -2*Is12*Is34*(2*B035R+s34*Cij135R(2,1)+4*Cij135R(4,2)
     -   +C0135R*P(11)+Cij135R(1,1)*P(12))
       FI(411) = -2*Is12*Is34*(2*B035I+s34*Cij135I(2,1)+4*Cij135I(4,2)
     -   +C0135I*P(11)+Cij135I(1,1)*P(12))
       F(411)=DCMPLX(FR(411),FI(411))
       FR(412) = -2*Is12*(C0135R+Cij135R(1,1)-3*Cij135R(2,1)-2*Cij135R
     -   (3,2))
       FI(412) = -2*Is12*(C0135I+Cij135I(1,1)-3*Cij135I(2,1)-2*Cij135I
     -   (3,2))
       F(412)=DCMPLX(FR(412),FI(412))
       FR(413) = (Is12*(Cij135R(2,1)*P(7)-2*(B035R+2*Cij135R(4,2)+Cij1
     -   35R(1,1)*P(6)+C0135R*P(8))))/2.
       FI(413) = (Is12*(Cij135I(2,1)*P(7)-2*(B035I+2*Cij135I(4,2)+Cij1
     -   35I(1,1)*P(6)+C0135I*P(8))))/2.
       F(413)=DCMPLX(FR(413),FI(413))
       FR(414) = Is12*Is34*(-3*s12*Cij135R(2,1)-s34*Cij135R(2,1)+4*Cij
     -   135R(4,2)+C0135R*P(1)+Cij135R(1,1)*P(2)+2*(B035R+s12*Cij135R
     -   (1,2)+s34*Cij135R(2,2)-Cij135R(3,2)*P(3)))
       FI(414) = Is12*Is34*(-3*s12*Cij135I(2,1)-s34*Cij135I(2,1)+4*Cij
     -   135I(4,2)+C0135I*P(1)+Cij135I(1,1)*P(2)+2*(B035I+s12*Cij135I
     -   (1,2)+s34*Cij135I(2,2)-Cij135I(3,2)*P(3)))
       F(414)=DCMPLX(FR(414),FI(414))
       FR(415) = Is12*Is34*(2*(B035R+s12*Cij135R(1,2)-s34*(Cij135R(2,1
     -   )+Cij135R(3,2))+2*Cij135R(4,2))+C0135R*P(4)+Cij135R(1,1)*P(5
     -   ))
       FI(415) = Is12*Is34*(2*(B035I+s12*Cij135I(1,2)-s34*(Cij135I(2,1
     -   )+Cij135I(3,2))+2*Cij135I(4,2))+C0135I*P(4)+Cij135I(1,1)*P(5
     -   ))
       F(415)=DCMPLX(FR(415),FI(415))
       FR(416) = Is12*Is34*(-2*B035R-s34*Cij135R(2,1)-4*Cij135R(4,2)-C
     -   ij135R(1,1)*P(12)+C0135R*P(13))
       FI(416) = Is12*Is34*(-2*B035I-s34*Cij135I(2,1)-4*Cij135I(4,2)-C
     -   ij135I(1,1)*P(12)+C0135I*P(13))
       F(416)=DCMPLX(FR(416),FI(416))
       FR(417) = Is12*Is34*(2*B035R+s34*Cij135R(2,1)+4*Cij135R(4,2)+C0
     -   135R*P(11)+Cij135R(1,1)*P(12))
       FI(417) = Is12*Is34*(2*B035I+s34*Cij135I(2,1)+4*Cij135I(4,2)+C0
     -   135I*P(11)+Cij135I(1,1)*P(12))
       F(417)=DCMPLX(FR(417),FI(417))
       FR(418) = Is12*(C0135R+Cij135R(1,1)-2*(Cij135R(2,1)-Cij135R(2,2
     -   )+Cij135R(3,2)))
       FI(418) = Is12*(C0135I+Cij135I(1,1)-2*(Cij135I(2,1)-Cij135I(2,2
     -   )+Cij135I(3,2)))
       F(418)=DCMPLX(FR(418),FI(418))
       FR(419) = Is12*(C0135R+Cij135R(1,1)-3*Cij135R(2,1)-2*Cij135R(3,
     -   2))
       FI(419) = Is12*(C0135I+Cij135I(1,1)-3*Cij135I(2,1)-2*Cij135I(3,
     -   2))
       F(419)=DCMPLX(FR(419),FI(419))
       FR(420) = (Is12*Is34*(-2*B035R-s34*Cij135R(2,1)-4*Cij135R(4,2)-
     -   Cij135R(1,1)*P(12)+C0135R*P(13)))/2.
       FI(420) = (Is12*Is34*(-2*B035I-s34*Cij135I(2,1)-4*Cij135I(4,2)-
     -   Cij135I(1,1)*P(12)+C0135I*P(13)))/2.
       F(420)=DCMPLX(FR(420),FI(420))
       FR(421) = 2*Is12*Is45*(-3*s45*Cij145R(2,1)+2*(B012R+s45*Cij145R
     -   (1,2)-s45*Cij145R(3,2))+4*Cij145R(4,2)+C0145R*P(14)+Cij145R(
     -   1,1)*P(15))
       FI(421) = 2*Is12*Is45*(-3*s45*Cij145I(2,1)+2*(B012I+s45*Cij145I
     -   (1,2)-s45*Cij145I(3,2))+4*Cij145I(4,2)+C0145I*P(14)+Cij145I(
     -   1,1)*P(15))
       F(421)=DCMPLX(FR(421),FI(421))
       FR(422) = -(Is12*Is45*P(16)*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1
     -   )*P(14)+C0145R*P(17)))
       FI(422) = -(Is12*Is45*P(16)*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1
     -   )*P(14)+C0145I*P(17)))
       F(422)=DCMPLX(FR(422),FI(422))
       FR(423) = 4*Is12*Is45*(B012R+C0145R*p5sq+s45*Cij145R(1,2)+2*Cij
     -   145R(4,2)+Cij145R(1,1)*P(18))
       FI(423) = 4*Is12*Is45*(B012I+C0145I*p5sq+s45*Cij145I(1,2)+2*Cij
     -   145I(4,2)+Cij145I(1,1)*P(18))
       F(423)=DCMPLX(FR(423),FI(423))
       FR(424) = Is12*Is45*P(16)*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)*
     -   P(14)+C0145R*P(17))
       FI(424) = Is12*Is45*P(16)*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)*
     -   P(14)+C0145I*P(17))
       F(424)=DCMPLX(FR(424),FI(424))
       FR(425) = -2*Is12*Is45*(-3*s45*Cij145R(2,1)+2*(B012R+s45*Cij145
     -   R(1,2)-s45*Cij145R(3,2))+4*Cij145R(4,2)+C0145R*P(14)+Cij145R
     -   (1,1)*P(15))
       FI(425) = -2*Is12*Is45*(-3*s45*Cij145I(2,1)+2*(B012I+s45*Cij145
     -   I(1,2)-s45*Cij145I(3,2))+4*Cij145I(4,2)+C0145I*P(14)+Cij145I
     -   (1,1)*P(15))
       F(425)=DCMPLX(FR(425),FI(425))
       FR(426) = -4*Is12*Is45*(B012R+C0145R*p5sq+s45*Cij145R(1,2)+2*Ci
     -   j145R(4,2)+Cij145R(1,1)*P(18))
       FI(426) = -4*Is12*Is45*(B012I+C0145I*p5sq+s45*Cij145I(1,2)+2*Ci
     -   j145I(4,2)+Cij145I(1,1)*P(18))
       F(426)=DCMPLX(FR(426),FI(426))
       FR(427) = Is12*(-8*Cij145R(1,1)+6*Cij145R(2,1)-4*(C0145R+Cij145
     -   R(1,2)-Cij145R(3,2)))
       FI(427) = Is12*(-8*Cij145I(1,1)+6*Cij145I(2,1)-4*(C0145I+Cij145
     -   I(1,2)-Cij145I(3,2)))
       F(427)=DCMPLX(FR(427),FI(427))
       FR(428) = -2*Is12*Is45*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)*P(1
     -   4)-C0145R*P(19))
       FI(428) = -2*Is12*Is45*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)*P(1
     -   4)-C0145I*P(19))
       F(428)=DCMPLX(FR(428),FI(428))
       FR(429) = 2*Is12*Is45*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)*P(14
     -   )+C0145R*P(17))
       FI(429) = 2*Is12*Is45*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)*P(14
     -   )+C0145I*P(17))
       F(429)=DCMPLX(FR(429),FI(429))
       FR(430) = -2*Is12*(C0145R+3*Cij145R(1,1)+2*Cij145R(1,2))
       FI(430) = -2*Is12*(C0145I+3*Cij145I(1,1)+2*Cij145I(1,2))
       F(430)=DCMPLX(FR(430),FI(430))
       FR(431) = (Is12*Is45*P(16)*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)
     -   *P(14)+C0145R*P(17)))/2.
       FI(431) = (Is12*Is45*P(16)*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)
     -   *P(14)+C0145I*P(17)))/2.
       F(431)=DCMPLX(FR(431),FI(431))
       FR(432) = Is12*Is45*(-3*s45*Cij145R(2,1)+2*(B012R+s45*Cij145R(1
     -   ,2)-s45*Cij145R(3,2))+4*Cij145R(4,2)+C0145R*P(14)+Cij145R(1,
     -   1)*P(15))
       FI(432) = Is12*Is45*(-3*s45*Cij145I(2,1)+2*(B012I+s45*Cij145I(1
     -   ,2)-s45*Cij145I(3,2))+4*Cij145I(4,2)+C0145I*P(14)+Cij145I(1,
     -   1)*P(15))
       F(432)=DCMPLX(FR(432),FI(432))
       FR(433) = 2*Is12*Is45*(B012R+C0145R*p5sq+s45*Cij145R(1,2)+2*Cij
     -   145R(4,2)+Cij145R(1,1)*P(18))
       FI(433) = 2*Is12*Is45*(B012I+C0145I*p5sq+s45*Cij145I(1,2)+2*Cij
     -   145I(4,2)+Cij145I(1,1)*P(18))
       F(433)=DCMPLX(FR(433),FI(433))
       FR(434) = Is12*Is45*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)*P(14)+
     -   C0145R*P(17))
       FI(434) = Is12*Is45*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)*P(14)+
     -   C0145I*P(17))
       F(434)=DCMPLX(FR(434),FI(434))
       FR(435) = Is12*Is45*(-2*(B012R+2*Cij145R(4,2))-Cij145R(1,1)*P(1
     -   4)+C0145R*P(19))
       FI(435) = Is12*Is45*(-2*(B012I+2*Cij145I(4,2))-Cij145I(1,1)*P(1
     -   4)+C0145I*P(19))
       F(435)=DCMPLX(FR(435),FI(435))
       FR(436) = Is12*(4*Cij145R(1,1)-3*Cij145R(2,1)+2*(C0145R+Cij145R
     -   (1,2)-Cij145R(3,2)))
       FI(436) = Is12*(4*Cij145I(1,1)-3*Cij145I(2,1)+2*(C0145I+Cij145I
     -   (1,2)-Cij145I(3,2)))
       F(436)=DCMPLX(FR(436),FI(436))
       FR(437) = Is12*(C0145R+3*Cij145R(1,1)+2*Cij145R(1,2))
       FI(437) = Is12*(C0145I+3*Cij145I(1,1)+2*Cij145I(1,2))
       F(437)=DCMPLX(FR(437),FI(437))
       FR(438) = (Is12*Is45*(-2*B012R-4*Cij145R(4,2)-Cij145R(1,1)*P(14
     -   )-C0145R*P(17)))/2.
       FI(438) = (Is12*Is45*(-2*B012I-4*Cij145I(4,2)-Cij145I(1,1)*P(14
     -   )-C0145I*P(17)))/2.
       F(438)=DCMPLX(FR(438),FI(438))
       KI(1) = EijI(2,4)+EijI(17,4)+EijI(18,4)+EijI(20,4)+4*EijI(26,4)
     -   -2*(EijI(8,4)+EijI(9,4)+EijI(23,4)+EijI(29,4))
       KR(1) = EijR(2,4)+EijR(17,4)+EijR(18,4)+EijR(20,4)+4*EijR(26,4)
     -   -2*(EijR(8,4)+EijR(9,4)+EijR(23,4)+EijR(29,4))
       K(1)=DCMPLX(KR(1),KI(1))
       KI(2) = EijI(2,4)-EijI(8,4)+EijI(11,4)-EijI(12,4)-3*(EijI(9,4)-
     -   EijI(20,4)-EijI(26,4)+EijI(29,4))
       KR(2) = EijR(2,4)-EijR(8,4)+EijR(11,4)-EijR(12,4)-3*(EijR(9,4)-
     -   EijR(20,4)-EijR(26,4)+EijR(29,4))
       K(2)=DCMPLX(KR(2),KI(2))
       KI(3) = EijI(2,4)-EijI(8,4)-EijI(10,4)+EijI(20,4)+EijI(27,4)-Ei
     -   jI(29,4)+EijI(30,4)-EijI(31,4)-2*(EijI(9,4)-EijI(26,4)-EijI(
     -   28,4)+EijI(35,4))
       KR(3) = EijR(2,4)-EijR(8,4)-EijR(10,4)+EijR(20,4)+EijR(27,4)-Ei
     -   jR(29,4)+EijR(30,4)-EijR(31,4)-2*(EijR(9,4)-EijR(26,4)-EijR(
     -   28,4)+EijR(35,4))
       K(3)=DCMPLX(KR(3),KI(3))
       KI(4) = EijI(2,4)-2*EijI(8,4)+EijI(17,4)
       KR(4) = EijR(2,4)-2*EijR(8,4)+EijR(17,4)
       K(4)=DCMPLX(KR(4),KI(4))
       KI(5) = -EijI(2,3)-EijI(5,3)+EijI(6,3)+EijI(9,3)+2*(EijI(8,3)-E
     -   ijI(17,3))
       KR(5) = -EijR(2,3)-EijR(5,3)+EijR(6,3)+EijR(9,3)+2*(EijR(8,3)-E
     -   ijR(17,3))
       K(5)=DCMPLX(KR(5),KI(5))
       KI(6) = -EijI(2,3)+EijI(8,3)+EijI(11,3)-EijI(12,3)+2*(EijI(9,3)
     -   -EijI(17,3))
       KR(6) = -EijR(2,3)+EijR(8,3)+EijR(11,3)-EijR(12,3)+2*(EijR(9,3)
     -   -EijR(17,3))
       K(6)=DCMPLX(KR(6),KI(6))
       KI(7) = -Dij1235I(2,3)+Dij1235I(6,3)
       KR(7) = -Dij1235R(2,3)+Dij1235R(6,3)
       K(7)=DCMPLX(KR(7),KI(7))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c Computation of Fa fucntions. Depends on the external currents, through
c the contraction of the moments with the currents
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrc(p1,mup3)
       p1mup5 = dotrc(p1,mup5)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup5 = dotrc(p2,mup5)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup5 = dotrc(p3,mup5)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup5 = dotrc(p4,mup5)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrc(p5,mup3)
       p5mup5 = dotrc(p5,mup5)
       mup2mup3=dotcc(mup2,mup3)
       mup2mup5=dotcc(mup2,mup5)
       mup3mup5=dotcc(mup3,mup5)
c      Print*," p1mup2 ",  p1mup2 
c      Print*," p1mup3 ",  p1mup3 
c      Print*," p1mup5 ",  p1mup5 
c      Print*," p2mup2 ",  p2mup2 
c      Print*," p2mup3 ",  p2mup3 
c      Print*," p2mup5 ",  p2mup5 
c      Print*," p3mup2 ",  p3mup2 
c      Print*," p3mup3 ",  p3mup3 
c      Print*," p3mup5 ",  p3mup5 
c      Print*," p4mup2 ",  p4mup2 
c      Print*," p4mup3 ",  p4mup3 
c      Print*," p4mup5 ",  p4mup5 
c      Print*," p5mup2 ",  p5mup2 
c      Print*," p5mup3 ",  p5mup3 
c      Print*," p5mup5 ",  p5mup5 
c       Print*," mup2mup3", mup2mup3  
c       Print*," mup2mup5",mup2mup5   
c       Print*,"mup3mup5", mup3mup5  
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call FaFunctionNoAbe_pmmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup
     -   3,p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5m
     -   up2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
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
       SMB(1) = SC1c(barpsi_p4,mup2,psi_p1,alpha)
       SMB(2) = SC1c(barpsi_p4,mup3,psi_p1,alpha)
       SMB(3) = SC1c(barpsi_p4,mup5,psi_p1,alpha)
       SMB(4) = SC1r(barpsi_p4,p3,psi_p1,alpha)
       SMB(5) = SC1r(barpsi_p4,p5,psi_p1,alpha)
       SMB(6) = SC3ccc(barpsi_p4,mup2,mup3,mup5,psi_p1,alpha)
       SMB(7) = SC3rcc(barpsi_p4,p3,mup2,mup3,psi_p1,alpha)
       SMB(8) = SC3rcc(barpsi_p4,p3,mup2,mup5,psi_p1,alpha)
       SMB(9) = SC3rcc(barpsi_p4,p3,mup3,mup5,psi_p1,alpha)
       SMB(10) = SC3rrc(barpsi_p4,p3,p5,mup2,psi_p1,alpha)
       SMB(11) = SC3rrc(barpsi_p4,p3,p5,mup3,psi_p1,alpha)
       SMB(12) = SC3rrc(barpsi_p4,p3,p5,mup5,psi_p1,alpha)
       SMB(13) = SC3rcc(barpsi_p4,p5,mup2,mup3,psi_p1,alpha)
       SMB(14) = SC3rcc(barpsi_p4,p5,mup2,mup5,psi_p1,alpha)
       SMB(15) = SC3rcc(barpsi_p4,p5,mup3,mup5,psi_p1,alpha)
       SMB(16) = SC5rrccc(barpsi_p4,p3,p5,mup2,mup3,mup5,psi_p1,alpha)
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
c************************************************************************************
c************************************************************************************
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
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c result(1) is the pentagon
c result(2,3) are boxes with boson 2,3 within the loop
c result(4) is triangle on q1, result(6) triangle on q6

       result(1) = Fa(1)*SMB(1)+Fa(2)*SMB(2)+Fa(3)*SMB(3)+Fa(4)*SMB(4)
     -   +Fa(5)*SMB(5)+F(194)*SMB(6)+Fa(6)*SMB(7)+Fa(7)*SMB(8)+Fa(8)*
     -   SMB(9)+Fa(9)*SMB(10)+Fa(10)*SMB(11)+Fa(11)*SMB(12)+Fa(12)*SM
     -   B(13)+Fa(13)*SMB(14)+Fa(14)*SMB(15)+F(268)*SMB(16)
       result(2) = Fa(15)*SMB(1)+Fa(16)*SMB(2)+Fa(17)*SMB(3)+Fa(18)*SM
     -   B(4)+Fa(19)*SMB(5)+F(316)*SMB(6)+Fa(20)*SMB(7)+p5mup3*F(320)
     -   *SMB(8)+Fa(21)*SMB(9)+Fa(22)*SMB(11)+mup2mup3*F(333)*SMB(12)
     -   +Fa(23)*SMB(13)+Fa(24)*SMB(14)+Fa(25)*SMB(15)+F(340)*SMB(16)
       result(3) = Fa(26)*SMB(1)+mup2mup5*F(351)*SMB(2)+mup2mup3*F(352
     -   )*SMB(3)+Fa(27)*SMB(4)+Fa(28)*SMB(5)+F(365)*SMB(6)+Fa(29)*SM
     -   B(7)+Fa(30)*SMB(8)+Fa(31)*SMB(10)+mup2mup5*F(382)*SMB(11)+mu
     -   p2mup3*F(383)*SMB(12)+Fa(32)*SMB(13)+Fa(33)*SMB(14)+F(390)*S
     -   MB(16)
       result(4) = Fa(34)*SMB(1)+Fa(35)*SMB(3)+Fa(36)*SMB(4)+Fa(37)*SM
     -   B(7)+p5mup3*F(398)*SMB(8)+p5mup2*F(399)*SMB(9)+mup2mup3*F(39
     -   9)*SMB(12)+Fa(38)*SMB(14)+F(400)*SMB(16)
       result(5) = Fa(39)*SMB(1)+mup2mup5*F(405)*SMB(2)+mup2mup3*F(406
     -   )*SMB(3)+Fa(40)*SMB(4)+Fa(41)*SMB(5)+F(413)*SMB(6)+Fa(42)*SM
     -   B(7)+p5mup3*F(416)*SMB(8)+mup2mup5*F(416)*SMB(11)+mup2mup3*F
     -   (417)*SMB(12)+Fa(43)*SMB(13)+Fa(44)*SMB(14)+F(420)*SMB(16)
       result(6) = Fa(45)*SMB(1)+mup2mup5*F(424)*SMB(2)+mup2mup3*F(422
     -   )*SMB(3)+Fa(46)*SMB(4)+Fa(47)*SMB(5)+F(431)*SMB(6)+Fa(48)*SM
     -   B(7)+mup3mup5*F(434)*SMB(10)+mup2mup5*F(435)*SMB(11)+mup2mup
     -   3*F(434)*SMB(12)+Fa(49)*SMB(13)+Fa(50)*SMB(14)+F(438)*SMB(16
     -   )
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       result(1) =  -result(1)
       result(2) =  -result(2)
       result(3) =  -result(3)
       result(4) =  -result(4)
       result(5) =  -result(5)
       result(6) =  -result(6)
       If (ngluon.eq.1) then
       result(1)=(result(1)+result(2)+result(3)+result(4)+result(5)+re
     -   sult(6))
       elseif (ngluon.eq.2) then
c  One gluon case is:  result(i)*CA. True color factor is f^abc T^c T^b= (-I)/2 CA.
c  To restore I should multiply for 2*I
c The color factor should be provided from outside
c result 1: posglu inside loop
c result 2: posglu left of loop
c result 3: posglu right of loop
       if (posgluon.eq.2) then
       result(1)=(0d0,2d0)*(result(1)+result(2))
       result(2)=(0d0,2d0)*(result(3)+result(5)+result(6))
       result(3)=(0d0,2d0)*(result(4))
       elseif (posgluon.eq.3) then
       result(1)=(0d0,2d0)*(result(1)+result(3))
       result(3)=(0d0,2d0)*(result(2)+result(4)+result(5))
       result(2)=(0d0,2d0)*(result(6))
       endif
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     -   . Look to the heading for explanation" 
       endif
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      If(gaugetest.eq.1) then 
cc      The mu_p2 is replaced for the incoming moment 
c      mup2->p2 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2=dotrr(p1,p2)
       p2mup2=dotrr(p2,p2)
       p3mup2=dotrr(p3,p2)
       p4mup2=dotrr(p4,p2)
       p5mup2=dotrr(p5,p2)
       mup2mup3=dotrc(p2,mup3)
       mup2mup5=dotrc(p2,mup5)
c       Print*," mup2mup3", mup2mup3  
c       Print*," mup2mup5",mup2mup5   
c       Print*,"mup3mup5", mup3mup5  
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FafunctionNoAbe_pmmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup
     -   3,p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5m
     -   up2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2)
c************************************************************************************
c************************************************************************************
       SMG(1) = -SMB(4)-SMB(5)
       SMG(2) = SMB(2)
       SMG(3) = SMB(3)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = 2*p1mup5*SMB(2)-2*p1mup3*SMB(3)-SMB(9)-SMB(15)
       SMG(7) = -(s34*SMB(2))-2*p1mup3*SMB(4)-SMB(11)
       SMG(8) = -(s34*SMB(3))-2*p1mup5*SMB(4)-SMB(12)
       SMG(9) = SMB(9)
       SMG(10) = (-s12+s34)*SMB(4)+s34*SMB(5)
       SMG(11) = SMB(11)
       SMG(12) = SMB(12)
       SMG(13) = (-s12+s34)*SMB(2)-2*p1mup3*SMB(5)+SMB(11)
       SMG(14) = (-s12+s34)*SMB(3)-2*p1mup5*SMB(5)+SMB(12)
       SMG(15) = SMB(15)
       SMG(16) = (-s12+s34)*SMB(9)+2*p1mup5*SMB(11)-2*p1mup3*SMB(12)+s
     -   34*SMB(15)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c       Print*," SMG(5) ", SMG(5)
c       Print*," SMG(6) ", SMG(6)
c       Print*," SMG(7) ", SMG(7)
c       Print*," SMG(8) ", SMG(8)
c       Print*," SMG(9) ", SMG(9)
c       Print*," SMG(10) ", SMG(10)
c       Print*," SMG(11) ", SMG(11)
c       Print*," SMG(12) ", SMG(12)
c       Print*," SMG(13) ", SMG(13)
c       Print*," SMG(14) ", SMG(14)
c       Print*," SMG(15) ", SMG(15)
c       Print*," SMG(16) ", SMG(16)
c************************************************************************************
c************************************************************************************
c    This should be equal to the boxline with momenta p1,p2+p3,p4,p5:
c************************************************************************************
c************************************************************************************
       resultgauge(1) = Fa(1)*SMG(1)+Fa(15)*SMG(1)+Fa(26)*SMG(1)+Fa(34
     -   )*SMG(1)+Fa(39)*SMG(1)+Fa(45)*SMG(1)+mup2mup5*F(351)*SMG(2)+
     -   mup2mup5*F(405)*SMG(2)+mup2mup5*F(424)*SMG(2)+Fa(2)*SMG(2)+F
     -   a(16)*SMG(2)+mup2mup3*F(352)*SMG(3)+mup2mup3*F(406)*SMG(3)+m
     -   up2mup3*F(422)*SMG(3)+Fa(3)*SMG(3)+Fa(17)*SMG(3)+Fa(35)*SMG(
     -   3)+Fa(4)*SMG(4)+Fa(18)*SMG(4)+Fa(27)*SMG(4)+Fa(36)*SMG(4)+Fa
     -   (40)*SMG(4)+Fa(46)*SMG(4)+Fa(5)*SMG(5)+Fa(19)*SMG(5)+Fa(28)*
     -   SMG(5)+Fa(41)*SMG(5)+Fa(47)*SMG(5)+F(194)*SMG(6)+F(316)*SMG(
     -   6)+F(365)*SMG(6)+F(413)*SMG(6)+F(431)*SMG(6)+Fa(6)*SMG(7)+Fa
     -   (20)*SMG(7)+Fa(29)*SMG(7)+Fa(37)*SMG(7)+Fa(42)*SMG(7)+Fa(48)
     -   *SMG(7)+p5mup3*F(320)*SMG(8)+p5mup3*F(398)*SMG(8)+p5mup3*F(4
     -   16)*SMG(8)+Fa(7)*SMG(8)+Fa(30)*SMG(8)+p5mup2*F(399)*SMG(9)+F
     -   a(8)*SMG(9)+Fa(21)*SMG(9)+mup3mup5*F(434)*SMG(10)+Fa(9)*SMG(
     -   10)+Fa(31)*SMG(10)+mup2mup5*F(382)*SMG(11)+mup2mup5*F(416)*S
     -   MG(11)+mup2mup5*F(435)*SMG(11)+Fa(10)*SMG(11)+Fa(22)*SMG(11)
     -   +mup2mup3*F(333)*SMG(12)+mup2mup3*F(383)*SMG(12)+mup2mup3*F(
     -   399)*SMG(12)+mup2mup3*F(417)*SMG(12)+mup2mup3*F(434)*SMG(12)
     -   +Fa(11)*SMG(12)+Fa(12)*SMG(13)+Fa(23)*SMG(13)+Fa(32)*SMG(13)
     -   +Fa(43)*SMG(13)+Fa(49)*SMG(13)+Fa(13)*SMG(14)+Fa(24)*SMG(14)
     -   +Fa(33)*SMG(14)+Fa(38)*SMG(14)+Fa(44)*SMG(14)+Fa(50)*SMG(14)
     -   +Fa(14)*SMG(15)+Fa(25)*SMG(15)+F(268)*SMG(16)+F(340)*SMG(16)
     -   +F(390)*SMG(16)+F(400)*SMG(16)+F(420)*SMG(16)+F(438)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
           resultgauge(2)=(0d0,0d0)      
           resultgauge(3)=(0d0,0d0)      
           return      
           endif      
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      If(gaugetest.eq.2) then 
cc      The mu_p2 is replaced for the incoming moment 
c      mup2->p2 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup3=dotrr(p1,p3)
       p2mup3=dotrr(p2,p3)
       p3mup3=dotrr(p3,p3)
       p4mup3=dotrr(p4,p3)
       p5mup3=dotrr(p5,p3)
       mup2mup3=dotrc(p3,mup2)
       mup3mup5=dotrc(p3,mup5)
c       Print*," mup2mup3", mup2mup3  
c       Print*," mup2mup5",mup2mup5   
c       Print*,"mup3mup5", mup3mup5  
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FafunctionNoAbe_pmmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup
     -   3,p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5m
     -   up2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = SMB(4)
       SMG(3) = SMB(3)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = 2*p3mup2*SMB(3)-SMB(8)
       SMG(7) = -(p3sq*SMB(1))+2*p3mup2*SMB(4)
       SMG(8) = SMB(8)
       SMG(9) = p3sq*SMB(3)
       SMG(10) = SMB(10)
       SMG(11) = (s12-s34-s45)*SMB(4)-p3sq*SMB(5)
       SMG(12) = SMB(12)
       SMG(13) = (-s12+s34+s45)*SMB(1)+2*p3mup2*SMB(5)+SMB(10)
       SMG(14) = SMB(14)
       SMG(15) = (s12-s34-s45)*SMB(3)-SMB(12)
       SMG(16) = (-s12+s34+s45)*SMB(8)+2*p3mup2*SMB(12)+p3sq*SMB(14)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c       Print*," SMG(5) ", SMG(5)
c       Print*," SMG(6) ", SMG(6)
c       Print*," SMG(7) ", SMG(7)
c       Print*," SMG(8) ", SMG(8)
c       Print*," SMG(9) ", SMG(9)
c       Print*," SMG(10) ", SMG(10)
c       Print*," SMG(11) ", SMG(11)
c       Print*," SMG(12) ", SMG(12)
c       Print*," SMG(13) ", SMG(13)
c       Print*," SMG(14) ", SMG(14)
c       Print*," SMG(15) ", SMG(15)
c       Print*," SMG(16) ", SMG(16)
c************************************************************************************
c************************************************************************************
c    This should be equal to the boxline with momenta p1,p2+p3,p4,p5:
c************************************************************************************
c************************************************************************************
       resultgauge(1) = Fa(1)*SMG(1)+Fa(15)*SMG(1)+Fa(26)*SMG(1)+Fa(34
     -   )*SMG(1)+Fa(39)*SMG(1)+Fa(45)*SMG(1)+mup2mup5*F(351)*SMG(2)+
     -   mup2mup5*F(405)*SMG(2)+mup2mup5*F(424)*SMG(2)+Fa(2)*SMG(2)+F
     -   a(16)*SMG(2)+mup2mup3*F(352)*SMG(3)+mup2mup3*F(406)*SMG(3)+m
     -   up2mup3*F(422)*SMG(3)+Fa(3)*SMG(3)+Fa(17)*SMG(3)+Fa(35)*SMG(
     -   3)+Fa(4)*SMG(4)+Fa(18)*SMG(4)+Fa(27)*SMG(4)+Fa(36)*SMG(4)+Fa
     -   (40)*SMG(4)+Fa(46)*SMG(4)+Fa(5)*SMG(5)+Fa(19)*SMG(5)+Fa(28)*
     -   SMG(5)+Fa(41)*SMG(5)+Fa(47)*SMG(5)+F(194)*SMG(6)+F(316)*SMG(
     -   6)+F(365)*SMG(6)+F(413)*SMG(6)+F(431)*SMG(6)+Fa(6)*SMG(7)+Fa
     -   (20)*SMG(7)+Fa(29)*SMG(7)+Fa(37)*SMG(7)+Fa(42)*SMG(7)+Fa(48)
     -   *SMG(7)+p5mup3*F(320)*SMG(8)+p5mup3*F(398)*SMG(8)+p5mup3*F(4
     -   16)*SMG(8)+Fa(7)*SMG(8)+Fa(30)*SMG(8)+p5mup2*F(399)*SMG(9)+F
     -   a(8)*SMG(9)+Fa(21)*SMG(9)+mup3mup5*F(434)*SMG(10)+Fa(9)*SMG(
     -   10)+Fa(31)*SMG(10)+mup2mup5*F(382)*SMG(11)+mup2mup5*F(416)*S
     -   MG(11)+mup2mup5*F(435)*SMG(11)+Fa(10)*SMG(11)+Fa(22)*SMG(11)
     -   +mup2mup3*F(333)*SMG(12)+mup2mup3*F(383)*SMG(12)+mup2mup3*F(
     -   399)*SMG(12)+mup2mup3*F(417)*SMG(12)+mup2mup3*F(434)*SMG(12)
     -   +Fa(11)*SMG(12)+Fa(12)*SMG(13)+Fa(23)*SMG(13)+Fa(32)*SMG(13)
     -   +Fa(43)*SMG(13)+Fa(49)*SMG(13)+Fa(13)*SMG(14)+Fa(24)*SMG(14)
     -   +Fa(33)*SMG(14)+Fa(38)*SMG(14)+Fa(44)*SMG(14)+Fa(50)*SMG(14)
     -   +Fa(14)*SMG(15)+Fa(25)*SMG(15)+F(268)*SMG(16)+F(340)*SMG(16)
     -   +F(390)*SMG(16)+F(400)*SMG(16)+F(420)*SMG(16)+F(438)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
           resultgauge(2)=(0d0,0d0)      
           resultgauge(3)=(0d0,0d0)      
           return      
           endif      
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      If(gaugetest.eq.3) then 
cc      The mu_p2 is replaced for the incoming moment 
c      mup2->p2 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup5=dotrr(p1,p5)
       p2mup5=dotrr(p2,p5)
       p3mup5=dotrr(p3,p5)
       p4mup5=dotrr(p4,p5)
       p5mup5=dotrr(p5,p5)
       mup2mup5=dotrc(p5,mup2)
       mup3mup5=dotrc(p5,mup3)
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FafunctionNoAbe_pmmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup
     -   3,p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5m
     -   up2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = SMB(2)
       SMG(3) = SMB(5)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = 2*p5mup3*SMB(1)-2*p5mup2*SMB(2)+SMB(13)
       SMG(7) = SMB(7)
       SMG(8) = 2*p5mup2*SMB(4)-SMB(10)
       SMG(9) = 2*p5mup3*SMB(4)-SMB(11)
       SMG(10) = SMB(10)
       SMG(11) = SMB(11)
       SMG(12) = p5sq*SMB(4)
       SMG(13) = SMB(13)
       SMG(14) = -(p5sq*SMB(1))+2*p5mup2*SMB(5)
       SMG(15) = -(p5sq*SMB(2))+2*p5mup3*SMB(5)
       SMG(16) = p5sq*SMB(7)+2*p5mup3*SMB(10)-2*p5mup2*SMB(11)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c       Print*," SMG(5) ", SMG(5)
c       Print*," SMG(6) ", SMG(6)
c       Print*," SMG(7) ", SMG(7)
c       Print*," SMG(8) ", SMG(8)
c       Print*," SMG(9) ", SMG(9)
c       Print*," SMG(10) ", SMG(10)
c       Print*," SMG(11) ", SMG(11)
c       Print*," SMG(12) ", SMG(12)
c       Print*," SMG(13) ", SMG(13)
c       Print*," SMG(14) ", SMG(14)
c       Print*," SMG(15) ", SMG(15)
c       Print*," SMG(16) ", SMG(16)
c************************************************************************************
c************************************************************************************
c    This should be equal to the boxline with momenta p1,p2+p3,p4,p5:
c************************************************************************************
c************************************************************************************
       resultgauge(1) = Fa(1)*SMG(1)+Fa(15)*SMG(1)+Fa(26)*SMG(1)+Fa(34
     -   )*SMG(1)+Fa(39)*SMG(1)+Fa(45)*SMG(1)+mup2mup5*F(351)*SMG(2)+
     -   mup2mup5*F(405)*SMG(2)+mup2mup5*F(424)*SMG(2)+Fa(2)*SMG(2)+F
     -   a(16)*SMG(2)+mup2mup3*F(352)*SMG(3)+mup2mup3*F(406)*SMG(3)+m
     -   up2mup3*F(422)*SMG(3)+Fa(3)*SMG(3)+Fa(17)*SMG(3)+Fa(35)*SMG(
     -   3)+Fa(4)*SMG(4)+Fa(18)*SMG(4)+Fa(27)*SMG(4)+Fa(36)*SMG(4)+Fa
     -   (40)*SMG(4)+Fa(46)*SMG(4)+Fa(5)*SMG(5)+Fa(19)*SMG(5)+Fa(28)*
     -   SMG(5)+Fa(41)*SMG(5)+Fa(47)*SMG(5)+F(194)*SMG(6)+F(316)*SMG(
     -   6)+F(365)*SMG(6)+F(413)*SMG(6)+F(431)*SMG(6)+Fa(6)*SMG(7)+Fa
     -   (20)*SMG(7)+Fa(29)*SMG(7)+Fa(37)*SMG(7)+Fa(42)*SMG(7)+Fa(48)
     -   *SMG(7)+p5mup3*F(320)*SMG(8)+p5mup3*F(398)*SMG(8)+p5mup3*F(4
     -   16)*SMG(8)+Fa(7)*SMG(8)+Fa(30)*SMG(8)+p5mup2*F(399)*SMG(9)+F
     -   a(8)*SMG(9)+Fa(21)*SMG(9)+mup3mup5*F(434)*SMG(10)+Fa(9)*SMG(
     -   10)+Fa(31)*SMG(10)+mup2mup5*F(382)*SMG(11)+mup2mup5*F(416)*S
     -   MG(11)+mup2mup5*F(435)*SMG(11)+Fa(10)*SMG(11)+Fa(22)*SMG(11)
     -   +mup2mup3*F(333)*SMG(12)+mup2mup3*F(383)*SMG(12)+mup2mup3*F(
     -   399)*SMG(12)+mup2mup3*F(417)*SMG(12)+mup2mup3*F(434)*SMG(12)
     -   +Fa(11)*SMG(12)+Fa(12)*SMG(13)+Fa(23)*SMG(13)+Fa(32)*SMG(13)
     -   +Fa(43)*SMG(13)+Fa(49)*SMG(13)+Fa(13)*SMG(14)+Fa(24)*SMG(14)
     -   +Fa(33)*SMG(14)+Fa(38)*SMG(14)+Fa(44)*SMG(14)+Fa(50)*SMG(14)
     -   +Fa(14)*SMG(15)+Fa(25)*SMG(15)+F(268)*SMG(16)+F(340)*SMG(16)
     -   +F(390)*SMG(16)+F(400)*SMG(16)+F(420)*SMG(16)+F(438)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
           resultgauge(2)=(0d0,0d0)      
           resultgauge(3)=(0d0,0d0)      
           return      
           endif      
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      If(gaugetest.eq.4) then 
cc      The mu_p2 is replaced for the incoming moment 
c      mup2->p2 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2=dotrr(p1,p2)
       p2mup2=dotrr(p2,p2)
       p3mup2=dotrr(p3,p2)
       p4mup2=dotrr(p4,p2)
       p5mup2=dotrr(p5,p2)
       mup2mup3=dotrc(p2,mup3)
       mup2mup5=dotrc(p2,mup5)
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FafunctionNoAbe_pmmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup
     -   3,p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5m
     -   up2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2)
c************************************************************************************
c************************************************************************************
       SMG(1) = -SMB(4)-SMB(5)
       SMG(2) = SMB(2)
       SMG(3) = SMB(3)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = 2*p1mup5*SMB(2)-2*p1mup3*SMB(3)-SMB(9)-SMB(15)
       SMG(7) = -(s34*SMB(2))-2*p1mup3*SMB(4)-SMB(11)
       SMG(8) = -(s34*SMB(3))-2*p1mup5*SMB(4)-SMB(12)
       SMG(9) = SMB(9)
       SMG(10) = (-s12+s34)*SMB(4)+s34*SMB(5)
       SMG(11) = SMB(11)
       SMG(12) = SMB(12)
       SMG(13) = (-s12+s34)*SMB(2)-2*p1mup3*SMB(5)+SMB(11)
       SMG(14) = (-s12+s34)*SMB(3)-2*p1mup5*SMB(5)+SMB(12)
       SMG(15) = SMB(15)
       SMG(16) = (-s12+s34)*SMB(9)+2*p1mup5*SMB(11)-2*p1mup3*SMB(12)+s
     -   34*SMB(15)
c************************************************************************************
c************************************************************************************
c    This should be equal to the boxline with momenta p1,p2+p3,p4,p5:
c************************************************************************************
c************************************************************************************
       resultgauge(1) = Fa(1)*SMG(1)+Fa(15)*SMG(1)+Fa(26)*SMG(1)+Fa(34
     -   )*SMG(1)+Fa(39)*SMG(1)+Fa(45)*SMG(1)+mup2mup5*F(351)*SMG(2)+
     -   mup2mup5*F(405)*SMG(2)+mup2mup5*F(424)*SMG(2)+Fa(2)*SMG(2)+F
     -   a(16)*SMG(2)+mup2mup3*F(352)*SMG(3)+mup2mup3*F(406)*SMG(3)+m
     -   up2mup3*F(422)*SMG(3)+Fa(3)*SMG(3)+Fa(17)*SMG(3)+Fa(35)*SMG(
     -   3)+Fa(4)*SMG(4)+Fa(18)*SMG(4)+Fa(27)*SMG(4)+Fa(36)*SMG(4)+Fa
     -   (40)*SMG(4)+Fa(46)*SMG(4)+Fa(5)*SMG(5)+Fa(19)*SMG(5)+Fa(28)*
     -   SMG(5)+Fa(41)*SMG(5)+Fa(47)*SMG(5)+F(194)*SMG(6)+F(316)*SMG(
     -   6)+F(365)*SMG(6)+F(413)*SMG(6)+F(431)*SMG(6)+Fa(6)*SMG(7)+Fa
     -   (20)*SMG(7)+Fa(29)*SMG(7)+Fa(37)*SMG(7)+Fa(42)*SMG(7)+Fa(48)
     -   *SMG(7)+p5mup3*F(320)*SMG(8)+p5mup3*F(398)*SMG(8)+p5mup3*F(4
     -   16)*SMG(8)+Fa(7)*SMG(8)+Fa(30)*SMG(8)+p5mup2*F(399)*SMG(9)+F
     -   a(8)*SMG(9)+Fa(21)*SMG(9)+mup3mup5*F(434)*SMG(10)+Fa(9)*SMG(
     -   10)+Fa(31)*SMG(10)+mup2mup5*F(382)*SMG(11)+mup2mup5*F(416)*S
     -   MG(11)+mup2mup5*F(435)*SMG(11)+Fa(10)*SMG(11)+Fa(22)*SMG(11)
     -   +mup2mup3*F(333)*SMG(12)+mup2mup3*F(383)*SMG(12)+mup2mup3*F(
     -   399)*SMG(12)+mup2mup3*F(417)*SMG(12)+mup2mup3*F(434)*SMG(12)
     -   +Fa(11)*SMG(12)+Fa(12)*SMG(13)+Fa(23)*SMG(13)+Fa(32)*SMG(13)
     -   +Fa(43)*SMG(13)+Fa(49)*SMG(13)+Fa(13)*SMG(14)+Fa(24)*SMG(14)
     -   +Fa(33)*SMG(14)+Fa(38)*SMG(14)+Fa(44)*SMG(14)+Fa(50)*SMG(14)
     -   +Fa(14)*SMG(15)+Fa(25)*SMG(15)+F(268)*SMG(16)+F(340)*SMG(16)
     -   +F(390)*SMG(16)+F(400)*SMG(16)+F(420)*SMG(16)+F(438)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c      The mu_p3 is replaced for the incoming moment 
c      mup3->p3 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrr(p1,p3)
       p1mup5 = dotrc(p1,mup5)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrr(p2,p3)
       p2mup5 = dotrc(p2,mup5)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrr(p3,p3)
       p3mup5 = dotrc(p3,mup5)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrr(p4,p3)
       p4mup5 = dotrc(p4,mup5)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrr(p5,p3)
       p5mup5 = dotrc(p5,mup5)
       mup2mup3=dotrc(p3,mup2)
       mup2mup5=dotcc(mup2,mup5)
       mup3mup5=dotrc(p3,mup5)
c       Print*," p1mup2 ",  p1mup2 
c       Print*," p1mup3 ",  p1mup3 
c       Print*," p1mup4 ",  p1mup4 
c       Print*," p2mup2 ",  p2mup2 
c       Print*," p2mup3 ",  p2mup3 
c       Print*," p2mup4 ",  p2mup4 
c       Print*," p3mup2 ",  p3mup2 
c       Print*," p3mup3 ",  p3mup3 
c       Print*," p3mup4 ",  p3mup4 
c       Print*," p4mup2 ",  p4mup2 
c       Print*," p4mup3 ",  p4mup3 
c       Print*," p4mup4 ",  p4mup4 
c       Print*," p5mup2 ",  p5mup2 
c       Print*," p5mup3 ",  p5mup3 
c       Print*," p5mup4 ",  p5mup4 
c       Print*," mup2mup3", mup2mup3  
c       Print*," mup2mup4",mup2mup4   
c       Print*,"mup3mup4", mup3mup4  
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FafunctionNoAbe_pmmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup
     -   3,p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5m
     -   up2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:3)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = SMB(4)
       SMG(3) = SMB(3)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = 2*p3mup2*SMB(3)-SMB(8)
       SMG(7) = -(p3sq*SMB(1))+2*p3mup2*SMB(4)
       SMG(8) = SMB(8)
       SMG(9) = p3sq*SMB(3)
       SMG(10) = SMB(10)
       SMG(11) = (s12-s34-s45)*SMB(4)-p3sq*SMB(5)
       SMG(12) = SMB(12)
       SMG(13) = (-s12+s34+s45)*SMB(1)+2*p3mup2*SMB(5)+SMB(10)
       SMG(14) = SMB(14)
       SMG(15) = (s12-s34-s45)*SMB(3)-SMB(12)
       SMG(16) = (-s12+s34+s45)*SMB(8)+2*p3mup2*SMB(12)+p3sq*SMB(14)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c       Print*," SMG(5) ", SMG(5)
c       Print*," SMG(6) ", SMG(6)
c       Print*," SMG(7) ", SMG(7)
c       Print*," SMG(8) ", SMG(8)
c       Print*," SMG(9) ", SMG(9)
c       Print*," SMG(10) ", SMG(10)
c       Print*," SMG(11) ", SMG(11)
c       Print*," SMG(12) ", SMG(12)
c       Print*," SMG(13) ", SMG(13)
c       Print*," SMG(14) ", SMG(14)
c       Print*," SMG(15) ", SMG(15)
c       Print*," SMG(16) ", SMG(16)
c************************************************************************************
c************************************************************************************
c    This should be equal to the subsctraction of boxline with momenta (p1,p2,p3+p4,p5) -(p1,p2+p3,p4,p5):
c************************************************************************************
c************************************************************************************
       resultgauge(2) = Fa(1)*SMG(1)+Fa(15)*SMG(1)+Fa(26)*SMG(1)+Fa(34
     -   )*SMG(1)+Fa(39)*SMG(1)+Fa(45)*SMG(1)+mup2mup5*F(351)*SMG(2)+
     -   mup2mup5*F(405)*SMG(2)+mup2mup5*F(424)*SMG(2)+Fa(2)*SMG(2)+F
     -   a(16)*SMG(2)+mup2mup3*F(352)*SMG(3)+mup2mup3*F(406)*SMG(3)+m
     -   up2mup3*F(422)*SMG(3)+Fa(3)*SMG(3)+Fa(17)*SMG(3)+Fa(35)*SMG(
     -   3)+Fa(4)*SMG(4)+Fa(18)*SMG(4)+Fa(27)*SMG(4)+Fa(36)*SMG(4)+Fa
     -   (40)*SMG(4)+Fa(46)*SMG(4)+Fa(5)*SMG(5)+Fa(19)*SMG(5)+Fa(28)*
     -   SMG(5)+Fa(41)*SMG(5)+Fa(47)*SMG(5)+F(194)*SMG(6)+F(316)*SMG(
     -   6)+F(365)*SMG(6)+F(413)*SMG(6)+F(431)*SMG(6)+Fa(6)*SMG(7)+Fa
     -   (20)*SMG(7)+Fa(29)*SMG(7)+Fa(37)*SMG(7)+Fa(42)*SMG(7)+Fa(48)
     -   *SMG(7)+p5mup3*F(320)*SMG(8)+p5mup3*F(398)*SMG(8)+p5mup3*F(4
     -   16)*SMG(8)+Fa(7)*SMG(8)+Fa(30)*SMG(8)+p5mup2*F(399)*SMG(9)+F
     -   a(8)*SMG(9)+Fa(21)*SMG(9)+mup3mup5*F(434)*SMG(10)+Fa(9)*SMG(
     -   10)+Fa(31)*SMG(10)+mup2mup5*F(382)*SMG(11)+mup2mup5*F(416)*S
     -   MG(11)+mup2mup5*F(435)*SMG(11)+Fa(10)*SMG(11)+Fa(22)*SMG(11)
     -   +mup2mup3*F(333)*SMG(12)+mup2mup3*F(383)*SMG(12)+mup2mup3*F(
     -   399)*SMG(12)+mup2mup3*F(417)*SMG(12)+mup2mup3*F(434)*SMG(12)
     -   +Fa(11)*SMG(12)+Fa(12)*SMG(13)+Fa(23)*SMG(13)+Fa(32)*SMG(13)
     -   +Fa(43)*SMG(13)+Fa(49)*SMG(13)+Fa(13)*SMG(14)+Fa(24)*SMG(14)
     -   +Fa(33)*SMG(14)+Fa(38)*SMG(14)+Fa(44)*SMG(14)+Fa(50)*SMG(14)
     -   +Fa(14)*SMG(15)+Fa(25)*SMG(15)+F(268)*SMG(16)+F(340)*SMG(16)
     -   +F(390)*SMG(16)+F(400)*SMG(16)+F(420)*SMG(16)+F(438)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(2) =-resultgauge(2)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c      The mu_p5 is replaced for the incoming moment 
c      mup5->p5 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrc(p1,mup3)
       p1mup5 = dotrr(p1,p5)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup5 = dotrr(p2,p5)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup5 = dotrr(p3,p5)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup5 = dotrr(p4,p5)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrc(p5,mup3)
       p5mup5 = dotrr(p5,p5)
       mup2mup3=dotcc(mup3,mup2)
       mup2mup5=dotrc(p5,mup2)
       mup3mup5=dotrc(p5,mup3)
c       Print*," p1mup2 ",  p1mup2 
c       Print*," p1mup3 ",  p1mup3 
c       Print*," p1mup4 ",  p1mup4 
c       Print*," p2mup2 ",  p2mup2 
c       Print*," p2mup3 ",  p2mup3 
c       Print*," p2mup4 ",  p2mup4 
c       Print*," p3mup2 ",  p3mup2 
c       Print*," p3mup3 ",  p3mup3 
c       Print*," p3mup4 ",  p3mup4 
c       Print*," p4mup2 ",  p4mup2 
c       Print*," p4mup3 ",  p4mup3 
c       Print*," p4mup4 ",  p4mup4 
c       Print*," p5mup2 ",  p5mup2 
c       Print*," p5mup3 ",  p5mup3 
c       Print*," p5mup4 ",  p5mup4 
c       Print*," mup2mup3", mup2mup3  
c       Print*," mup2mup4",mup2mup4   
c       Print*,"mup3mup4", mup3mup4  
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FafunctionNoAbe_pmmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup
     -   3,p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5m
     -   up2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:3)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = SMB(2)
       SMG(3) = SMB(5)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = 2*p5mup3*SMB(1)-2*p5mup2*SMB(2)+SMB(13)
       SMG(7) = SMB(7)
       SMG(8) = 2*p5mup2*SMB(4)-SMB(10)
       SMG(9) = 2*p5mup3*SMB(4)-SMB(11)
       SMG(10) = SMB(10)
       SMG(11) = SMB(11)
       SMG(12) = p5sq*SMB(4)
       SMG(13) = SMB(13)
       SMG(14) = -(p5sq*SMB(1))+2*p5mup2*SMB(5)
       SMG(15) = -(p5sq*SMB(2))+2*p5mup3*SMB(5)
       SMG(16) = p5sq*SMB(7)+2*p5mup3*SMB(10)-2*p5mup2*SMB(11)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c       Print*," SMG(5) ", SMG(5)
c       Print*," SMG(6) ", SMG(6)
c       Print*," SMG(7) ", SMG(7)
c       Print*," SMG(8) ", SMG(8)
c       Print*," SMG(9) ", SMG(9)
c       Print*," SMG(10) ", SMG(10)
c       Print*," SMG(11) ", SMG(11)
c       Print*," SMG(12) ", SMG(12)
c       Print*," SMG(13) ", SMG(13)
c       Print*," SMG(14) ", SMG(14)
c       Print*," SMG(15) ", SMG(15)
c       Print*," SMG(16) ", SMG(16)
c************************************************************************************
c************************************************************************************
c    This should be equal to the subsctraction of boxline with momenta -(p1,p2,p3+p4,p5) :
c************************************************************************************
c************************************************************************************
       resultgauge(3) = Fa(1)*SMG(1)+Fa(15)*SMG(1)+Fa(26)*SMG(1)+Fa(34
     -   )*SMG(1)+Fa(39)*SMG(1)+Fa(45)*SMG(1)+mup2mup5*F(351)*SMG(2)+
     -   mup2mup5*F(405)*SMG(2)+mup2mup5*F(424)*SMG(2)+Fa(2)*SMG(2)+F
     -   a(16)*SMG(2)+mup2mup3*F(352)*SMG(3)+mup2mup3*F(406)*SMG(3)+m
     -   up2mup3*F(422)*SMG(3)+Fa(3)*SMG(3)+Fa(17)*SMG(3)+Fa(35)*SMG(
     -   3)+Fa(4)*SMG(4)+Fa(18)*SMG(4)+Fa(27)*SMG(4)+Fa(36)*SMG(4)+Fa
     -   (40)*SMG(4)+Fa(46)*SMG(4)+Fa(5)*SMG(5)+Fa(19)*SMG(5)+Fa(28)*
     -   SMG(5)+Fa(41)*SMG(5)+Fa(47)*SMG(5)+F(194)*SMG(6)+F(316)*SMG(
     -   6)+F(365)*SMG(6)+F(413)*SMG(6)+F(431)*SMG(6)+Fa(6)*SMG(7)+Fa
     -   (20)*SMG(7)+Fa(29)*SMG(7)+Fa(37)*SMG(7)+Fa(42)*SMG(7)+Fa(48)
     -   *SMG(7)+p5mup3*F(320)*SMG(8)+p5mup3*F(398)*SMG(8)+p5mup3*F(4
     -   16)*SMG(8)+Fa(7)*SMG(8)+Fa(30)*SMG(8)+p5mup2*F(399)*SMG(9)+F
     -   a(8)*SMG(9)+Fa(21)*SMG(9)+mup3mup5*F(434)*SMG(10)+Fa(9)*SMG(
     -   10)+Fa(31)*SMG(10)+mup2mup5*F(382)*SMG(11)+mup2mup5*F(416)*S
     -   MG(11)+mup2mup5*F(435)*SMG(11)+Fa(10)*SMG(11)+Fa(22)*SMG(11)
     -   +mup2mup3*F(333)*SMG(12)+mup2mup3*F(383)*SMG(12)+mup2mup3*F(
     -   399)*SMG(12)+mup2mup3*F(417)*SMG(12)+mup2mup3*F(434)*SMG(12)
     -   +Fa(11)*SMG(12)+Fa(12)*SMG(13)+Fa(23)*SMG(13)+Fa(32)*SMG(13)
     -   +Fa(43)*SMG(13)+Fa(49)*SMG(13)+Fa(13)*SMG(14)+Fa(24)*SMG(14)
     -   +Fa(33)*SMG(14)+Fa(38)*SMG(14)+Fa(44)*SMG(14)+Fa(50)*SMG(14)
     -   +Fa(14)*SMG(15)+Fa(25)*SMG(15)+F(268)*SMG(16)+F(340)*SMG(16)
     -   +F(390)*SMG(16)+F(400)*SMG(16)+F(420)*SMG(16)+F(438)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(3) =-resultgauge(3)
           endif      
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End

