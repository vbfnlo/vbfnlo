       subroutine penlineNoAbemmm(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mup2
     -   ,mup3,mup5,alpha,musqIn,gaugetest,comp,resultgauge,result)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 3/7/2008
c Modified:4/8/2008
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
       Complex*16 barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3),mup5(0:3) 
       Real*8   p1sq, p1p2, p1p3, p1p4, p1p5, p2sq, p2p3, p2p4, 
     -          p2p5, p3sq, p3p4, p3p5, p4sq, p4p5, p5sq
       Real*8   s12, s13, s14, s15, s23, s24, s25, s34, s35, s45
       Complex*16   p1mup2, p1mup3, p1mup5, p2mup2, p2mup3, p2mup5, 
     -          p3mup2, p3mup3, p3mup5, p4mup2, p4mup3, p4mup5, 
     -          p5mup2, p5mup3, p5mup5
       Complex*16 mup2mup3,mup2mup5,mup3mup5
       Real*8 dotrr
       Complex*16 B0fin,C0fin,D0fin,E0fin
       EXTERNAL dotrr,B0fin,C0fin,D0fin,E0fin
       Complex*16 B012,B013,B014,B015,B023,B024,B025,B034,B035,B045
      Real*8 B012R,B013R,B014R,B015R,B023R,B024R,B025R,B034R,B035R,B045R
      Real*8 B012I,B013I,B014I,B015I,B023I,B024I,B025I,B034I,B035I,B045I
       Complex*16 C0123,C0124,C0125,C0134,C0135,C0145,C0234,C0235,C024
     -   5,C0345
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
       Complex*16 D01234,D01235,D01245,D01345,D02345,EE0
        Real*8 D01234R,D01235R,D01245R,D01345R,D02345R,EE0R
        Real*8 D01234I,D01235I,D01245I,D01345I,D02345I,EE0I
        Real*8 Dij1234R(13,3),Dij1235R(13,3),Dij1245R(13,3),Dij1345R(13
     -   ,3),Dij2345R(13,3),EijR(46,4)
        Real*8 Dij1234I(13,3),Dij1235I(13,3),Dij1245I(13,3),Dij1345I(13
     -   ,3),Dij2345I(13,3),EijI(46,4)
       Complex*16 SMB(16),SMG(16) ,Fa(14),F(280)
      Real*8 FI(280),FR(280)
      Complex*16 SC1c,SC3ccc,SC3rcc,SC3rrc,SC5rrccc,SC1r,dotrc,dotcc,result
       Complex*16 resultgauge(3) 
       Real*8 musq,p4t,Is12,Is34,Is15s34,Is45,P(685) 
       EXTERNAL   dotrc,dotcc,SC1c,SC1r,SC3ccc,SC3rcc,SC3rrc,SC5rrccc
       Integer alpha,comp,gaugetest
       common/FfunctionsNOABEmmm/F,P
       SAVE/FfunctionsNOABEmmm/
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
       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)
       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)
       p5sq = dotrr(p5,p5)
       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s15 = (p1sq +p5sq+ 2*p1p5) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s25 = (p2sq +p5sq+ 2*p2p5) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s35 = (p3sq +p5sq+ 2*p3p5) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
c       Write(*,'(a5,F20.10)')," p1sq ", p1sq 
c       Write(*,'(a5,F20.10)')," p1p2 ", p1p2
c       Write(*,'(a5,F20.10)')," p1p3 ", p1p3
c       Write(*,'(a5,F20.10)')," p1p4 ", p1p4
c       Write(*,'(a5,F20.10)')," p1p5 ", p1p5
c       Write(*,'(a5,F20.10)')," p2sq ", p2sq 
c       Write(*,'(a5,F20.10)')," p2p3 ", p2p3
c       Write(*,'(a5,F20.10)')," p2p4 ", p2p4
c       Write(*,'(a5,F20.10)')," p2p5 ", p2p5
c       Write(*,'(a5,F20.10)')," p3sq ", p3sq 
c       Write(*,'(a5,F20.10)')," p3p4 ", p3p4
c       Write(*,'(a5,F20.10)')," p3p5 ", p3p5
c       Write(*,'(a5,F20.10)')," p4sq ", p4sq 
c       Write(*,'(a5,F20.10)')," p4p5 ", p4p5
c       Write(*,'(a5,F20.10)')," p5sq ", p5sq 
       p4t=p1sq+p2sq+p3sq+p5sq+2d0*(p1p2+p1p3+p1p5+p2p3+p2p5+p3p5)
       If(abs(p4t).gt.1d-07) then
       If(p5sq.gt.1d-05)then
       p5(0)=-(p1(0)+p2(0)+p3(0)+p4(0))
       p5(1)=-(p1(1)+p2(1)+p3(1)+p4(1))
       p5(2)=-(p1(2)+p2(2)+p3(2)+p4(2))
       p5(3)=-(p1(3)+p2(3)+p3(3)+p4(3))
       write(*,*)"Momemtum is not conserved: p4sq=0. here=",p4t
       write(*,*)"Corrected by p5=-(p1+p2+p3+p4)"
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p1p5 = dotrr(p1,p5)
       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)
       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)
       p5sq = dotrr(p5,p5)
       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s15 = (p1sq +p5sq+ 2*p1p5) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s25 = (p2sq +p5sq+ 2*p2p5) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s35 = (p3sq +p5sq+ 2*p3p5) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
       elseIf(p2sq.gt.1d-05)then
       p2(0)=-(p1(0)+p4(0)+p3(0)+p5(0))
       p2(1)=-(p1(1)+p4(1)+p3(1)+p5(1))
       p2(2)=-(p1(2)+p4(2)+p3(2)+p5(2))
       p2(3)=-(p1(3)+p4(3)+p3(3)+p5(3))
       write(*,*)"Momemtum is not conserved: p4sq=0. here=",p4t
       write(*,*)"Corrected by p2=-(p1+p3+p4+p5)"
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p1p5 = dotrr(p1,p5)
       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)
       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)
       p5sq = dotrr(p5,p5)
       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s15 = (p1sq +p5sq+ 2*p1p5) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s25 = (p2sq +p5sq+ 2*p2p5) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s35 = (p3sq +p5sq+ 2*p3p5) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
       elseIf(p3sq.gt.1d-05)then
       p3(0)=-(p1(0)+p4(0)+p2(0)+p5(0))
       p3(1)=-(p1(1)+p4(1)+p2(1)+p5(1))
       p3(2)=-(p1(2)+p4(2)+p2(2)+p5(2))
       p3(3)=-(p1(3)+p4(3)+p2(3)+p5(3))
       write(*,*)"Momemtum is not conserved: p4sq=0. here=",p4t
       write(*,*) "Corrected by p3=-(p1+p2+p4+p5)"
       p1sq = dotrr(p1,p1)
       p1p2 = dotrr(p1,p2)
       p1p3 = dotrr(p1,p3)
       p1p4 = dotrr(p1,p4)
       p1p5 = dotrr(p1,p5)
       p2sq = dotrr(p2,p2)
       p2p3 = dotrr(p2,p3)
       p2p4 = dotrr(p2,p4)
       p2p5 = dotrr(p2,p5)
       p3sq = dotrr(p3,p3)
       p3p4 = dotrr(p3,p4)
       p3p5 = dotrr(p3,p5)
       p4sq = dotrr(p4,p4)
       p4p5 = dotrr(p4,p5)
       p5sq = dotrr(p5,p5)
       s12 = (p1sq +p2sq+ 2*p1p2) 
       s13 = (p1sq +p3sq+ 2*p1p3) 
       s14 = (p1sq +p4sq+ 2*p1p4) 
       s15 = (p1sq +p5sq+ 2*p1p5) 
       s23 = (p2sq +p3sq+ 2*p2p3) 
       s24 = (p2sq +p4sq+ 2*p2p4) 
       s25 = (p2sq +p5sq+ 2*p2p5) 
       s34 = (p3sq +p4sq+ 2*p3p4) 
       s35 = (p3sq +p5sq+ 2*p3p5) 
       s45 = (p4sq +p5sq+ 2*p4p5) 
       else
       write(*,*)"Momemtum is not conserved: p4sq=0 here=",p4t
       endif
       endif
       Is12=1d0/s12
       Is15s34=1d0/(s15*s34)
       Is34=1d0/s34
       Is45=1d0/s45
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
        if (comp.gt.0) then    
c    Calling C_ij,D_ij Functions    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       B012=B0fin(p1sq,musq)             
       B013=B0fin(s12,musq)   
       B014=B0fin(s45,musq)   
       B015=B0fin(p5sq,musq)              
       B023=B0fin(p2sq,musq)  
       B024=B0fin(s23,musq)   
       B025=B0fin(s15,musq)   
       B034=B0fin(p3sq,musq)  
       B035=B0fin(s34,musq)  
       B045=B0fin(p4sq,musq)  
       B012R=Dble(B012)          
       B013R=Dble(B013) 
       B014R=Dble(B014)   
       B015R=Dble(B015)           
       B023R=Dble(B023 ) 
       B024R=Dble( B024)  
       B025R=Dble(B025)   
       B034R=Dble(B034)  
       B035R=Dble(B035)  
       B045R=Dble(B045)  
       B012I=DIMAG(B012)          
       B013I=DIMAG(B013) 
       B014I=DIMAG(B014)   
       B015I=DIMAG(B015)           
       B023I=DIMAG(B023) 
       B024I=DIMAG(B024)  
       B025I=DIMAG(B025)   
       B034I=DIMAG(B034)  
       B035I=DIMAG(B035)  
       B045I=DIMAG(B045)  
c************************************************************************************
c************************************************************************************
c************************************************************************************
       C0145=C0fin(s45,p4sq,p5sq,musq)  
       C0125=C0fin(p1sq,s15,p5sq,musq)  
       C0123=C0fin(p1sq,p2sq,s12,musq)  
       C0124=C0fin(p1sq,s23,s45,musq)  
       C0135=C0fin(s12,s34,p5sq,musq)  
       C0134=C0fin(s12,p3sq,s45,musq)  
       C0234=C0fin(p2sq,p3sq,s23,musq)  
       C0235=C0fin(p2sq,s34,s15,musq)  
       C0245=C0fin(s23,p4sq,s15,musq)  
       C0345=C0fin(p3sq,p4sq,s34,musq)  
c************************************************************************************
c************************************************************************************
       call tens_red3_new_Re_Com(s45,p4sq,p5sq,B045,B015,B014,C0145,C0
     -   145R,C0145I,Cij145R,Cij145I) 
       call tens_red3_new_Re_Com(p1sq,s15,p5sq,B025,B015,B012,C0125,C0
     -   125R,C0125I,Cij125R,Cij125I) 
       call tens_red3_new_Re_Com(p1sq,p2sq,s12,B023,B013,B012,C0123,C0
     -   123R,C0123I,Cij123R,Cij123I) 
       call tens_red3_new_Re_Com(p1sq,s23,s45,B024,B014,B012,C0124,C01
     -   24R,C0124I,Cij124R,Cij124I) 
       call tens_red3_new_Re_Com(s12,s34,p5sq,B035,B015,B013,C0135,C01
     -   35R,C0135I,Cij135R,Cij135I) 
       call tens_red3_new_Re_Com(s12,p3sq,s45,B034,B014,B013,C0134,C01
     -   34R,C0134I,Cij134R,Cij134I) 
       call tens_red3_new_Re_Com(p2sq,p3sq,s23,B034,B024,B023,C0234,C0
     -   234R,C0234I,Cij234R,Cij234I) 
       call tens_red3_new_Re_Com(p2sq,s34,s15,B035,B025,B023,C0235,C02
     -   35R,C0235I,Cij235R,Cij235I) 
       call tens_red3_new_Re_Com(s23,p4sq,s15,B045,B025,B024,C0245,C02
     -   45R,C0245I,Cij245R,Cij245I) 
       call tens_red3_new_Re_Com(p3sq,p4sq,s34,B045,B035,B034,C0345,C0
     -   345R,C0345I,Cij345R,Cij345I) 
c************************************************************************************
c************************************************************************************
       D02345=D0fin(s23,s34,p2sq,p3sq,p4sq,s15,musq)
       D01345=D0fin(s45,s34,s12,p3sq,p4sq,p5sq,musq)
       D01245=D0fin(s45,s15,p1sq,s23,p4sq,p5sq,musq)
       D01235=D0fin(s12,s15,p1sq,p2sq,s34,p5sq,musq)
       D01234=D0fin(s12,s23,p1sq,p2sq,p3sq,s45,musq)
c************************************************************************************
c************************************************************************************
       call tens_red4_new_Re_Com(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
     -   C0245R,C0235R,C0234R,Cij345R,Cij245R,Cij235R,Cij234R,C0345I,
     -   C0245I,C0235I,C0234I,Cij345I,Cij245I,Cij235I,Cij234I,D02345,
     -   D02345R,D02345I,Dij2345R,Dij2345I)
        call tens_red4_new_Re_Com(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
     -   4,C0345R,C0145R,C0135R,C0134R,Cij345R,Cij145R,Cij135R,Cij134
     -   R,C0345I,C0145I,C0135I,C0134I,Cij345I,Cij145I,Cij135I,Cij134
     -   I,D01345,D01345R,D01345I,Dij1345R,Dij1345I)
        call tens_red4_new_Re_Com(0d0,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p4
     -   ,C0245R,C0145R,C0125R,C0124R,Cij245R,Cij145R,Cij125R,Cij124R
     -   ,C0245I,C0145I,C0125I,C0124I,Cij245I,Cij145I,Cij125I,Cij124I
     -   ,D01245,D01245R,D01245I,Dij1245R,Dij1245I)
        call tens_red4_new_Re_Com(0d0,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p4
     -   ,C0235R,C0135R,C0125R,C0123R,Cij235R,Cij135R,Cij125R,Cij123R
     -   ,C0235I,C0135I,C0125I,C0123I,Cij235I,Cij135I,Cij125I,Cij123I
     -   ,D01235,D01235R,D01235I,Dij1235R,Dij1235I)
        call tens_red4_new_Re_Com(0d0,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,C
     -   0134R,C0124R,C0123R,Cij234R,Cij134R,Cij124R,Cij123R,C0234I,C
     -   0134I,C0124I,C0123I,Cij234I,Cij134I,Cij124I,Cij123I,D01234,D
     -   01234R,D01234I,Dij1234R,Dij1234I)
c************************************************************************************
c************************************************************************************
       EE0=E0fin(p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
c       p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
     -   ,D02345
     -   ,D01345,D01245,D01235,D01234)
       EE0R=Dble(EE0) 
       EE0I=DIMAG(EE0) 
c************************************************************************************
c************************************************************************************
       call tens_red5_new_Re_Com(0d0,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2p
     -   3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,Dij2345R
     -   ,Dij1345R,Dij1245R,Dij1235R,Dij1234R,D02345I,D01345I,D01245I
     -   ,D01235I,D01234I,Dij2345I,Dij1345I,Dij1245I,Dij1235I,Dij1234
     -   I,EijR,EijI)
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
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       P(1) = p2sq+p3sq-s23-s34
       P(2) = p3sq+s15-s23-s34
       FR(1) = -4*(D02345R+Dij2345R(1,1)+s12*(EijR(2,2)-EijR(5,2))+s34
     -   *(EijR(6,2)-EijR(8,2))-p2sq*(EijR(2,3)+EijR(5,3)-2*EijR(8,3)
     -   )+4*EijR(21,3)-2*(EijR(36,4)+EijR(37,4)+2*(EijR(22,3)-EijR(4
     -   0,4)))+(EijR(6,3)+EijR(9,3)-2*EijR(17,3))*P(1)-(EijR(7,3)+Ei
     -   jR(10,3)-2*EijR(18,3))*P(2))
       FI(1) = -4*(D02345I+Dij2345I(1,1)+s12*(EijI(2,2)-EijI(5,2))+s34
     -   *(EijI(6,2)-EijI(8,2))-p2sq*(EijI(2,3)+EijI(5,3)-2*EijI(8,3)
     -   )+4*EijI(21,3)-2*(EijI(36,4)+EijI(37,4)+2*(EijI(22,3)-EijI(4
     -   0,4)))+(EijI(6,3)+EijI(9,3)-2*EijI(17,3))*P(1)-(EijI(7,3)+Ei
     -   jI(10,3)-2*EijI(18,3))*P(2))
       F(1)=DCMPLX(FR(1),FI(1))
       P(3) = p5sq-s15
       P(4) = p5sq-s12-s15+s34
       P(5) = p2sq+s12
       P(6) = s12-s34
       P(7) = p5sq-s45
       P(8) = p5sq-s15+s34
       P(9) = p2sq-3*s12+2*P(8)
       P(10) = p2sq-s12-s23+s45
       P(11) = p5sq+s15+s23-s45-2*P(6)
       P(12) = p2sq+s12-s23-s45
       P(13) = p5sq-2*s12-s15-s23-s45
       P(14) = p3sq-s34
       P(15) = p2sq+s12-s23-s45+2*P(14)
       P(16) = s12-s45
       P(17) = p2sq-p5sq+s12-s15+s23+s45-2*P(14)
       P(18) = p5sq+s15-s23-s45+2*P(14)
       FR(2) = 2*(D02345R+Dij2345R(1,2)+Dij2345R(3,1)-Dij2345R(5,2)-6*
     -   EijR(11,2)+4*(EijR(37,4)-EijR(40,4)+EijR(42,4)-EijR(44,4))-(
     -   EijR(2,2)-EijR(2,3)+EijR(8,3))*P(5)+EijR(5,2)*P(9)-EijR(6,2)
     -   *P(10)-EijR(7,2)*P(11)+EijR(8,2)*P(12)-EijR(9,2)*P(13)+(-Eij
     -   R(9,3)+EijR(17,3)-EijR(19,3)+EijR(20,3))*P(15)+2*(EijR(22,3)
     -   -EijR(24,3)-(EijR(1,1)+EijR(1,2))*P(3)+EijR(2,1)*P(4)+EijR(4
     -   ,1)*P(6)+EijR(4,2)*P(7)-EijR(10,2)*P(16))+(-EijR(10,3)+EijR(
     -   18,3))*P(17)+(EijR(14,3)-EijR(15,3))*P(18))
       FI(2) = 2*(D02345I+Dij2345I(1,2)+Dij2345I(3,1)-Dij2345I(5,2)-6*
     -   EijI(11,2)+4*(EijI(37,4)-EijI(40,4)+EijI(42,4)-EijI(44,4))-(
     -   EijI(2,2)-EijI(2,3)+EijI(8,3))*P(5)+EijI(5,2)*P(9)-EijI(6,2)
     -   *P(10)-EijI(7,2)*P(11)+EijI(8,2)*P(12)-EijI(9,2)*P(13)+(-Eij
     -   I(9,3)+EijI(17,3)-EijI(19,3)+EijI(20,3))*P(15)+2*(EijI(22,3)
     -   -EijI(24,3)-(EijI(1,1)+EijI(1,2))*P(3)+EijI(2,1)*P(4)+EijI(4
     -   ,1)*P(6)+EijI(4,2)*P(7)-EijI(10,2)*P(16))+(-EijI(10,3)+EijI(
     -   18,3))*P(17)+(EijI(14,3)-EijI(15,3))*P(18))
       F(2)=DCMPLX(FR(2),FI(2))
       P(19) = s15+s34
       P(20) = p2sq-3*s12-s23+s45+2*P(19)
       P(21) = p5sq-s15+s23-s45
       P(22) = p2sq+3*s12+2*s15-s23-s45
       P(23) = p5sq+s15-s23-s45
       P(24) = p2sq+p3sq+s12-s34
       P(25) = -s23-s45+2*P(24)
       FR(3) = 2*(D02345R+Dij2345R(1,2)+Dij2345R(2,1)-Dij2345R(4,2)-6*
     -   EijR(11,2)+4*(EijR(37,4)-EijR(40,4)+EijR(41,4)-EijR(43,4))-(
     -   EijR(2,2)-EijR(2,3)+EijR(8,3))*P(5)+EijR(5,2)*P(9)+(-EijR(11
     -   ,3)+EijR(12,3))*P(15)+2*(EijR(22,3)-EijR(23,3)-(EijR(1,1)+Ei
     -   jR(1,2))*P(3)+EijR(2,1)*P(4)+EijR(3,1)*P(6)+EijR(10,2)*P(7)-
     -   EijR(3,2)*P(16))+(EijR(10,3)-EijR(18,3)+EijR(19,3)-EijR(20,3
     -   ))*P(18)-EijR(6,2)*P(20)-EijR(7,2)*P(21)+EijR(8,2)*P(22)-Eij
     -   R(9,2)*P(23)+(-EijR(9,3)+EijR(17,3))*P(25))
       FI(3) = 2*(D02345I+Dij2345I(1,2)+Dij2345I(2,1)-Dij2345I(4,2)-6*
     -   EijI(11,2)+4*(EijI(37,4)-EijI(40,4)+EijI(41,4)-EijI(43,4))-(
     -   EijI(2,2)-EijI(2,3)+EijI(8,3))*P(5)+EijI(5,2)*P(9)+(-EijI(11
     -   ,3)+EijI(12,3))*P(15)+2*(EijI(22,3)-EijI(23,3)-(EijI(1,1)+Ei
     -   jI(1,2))*P(3)+EijI(2,1)*P(4)+EijI(3,1)*P(6)+EijI(10,2)*P(7)-
     -   EijI(3,2)*P(16))+(EijI(10,3)-EijI(18,3)+EijI(19,3)-EijI(20,3
     -   ))*P(18)-EijI(6,2)*P(20)-EijI(7,2)*P(21)+EijI(8,2)*P(22)-Eij
     -   I(9,2)*P(23)+(-EijI(9,3)+EijI(17,3))*P(25))
       F(3)=DCMPLX(FR(3),FI(3))
       P(26) = 2*p5sq+s12-s34
       P(27) = 2*p5sq+5*s12-s34
       P(28) = 3*s12+s34
       P(29) = s12+s34
       P(30) = s34+s45
       P(31) = p5sq-P(30)
       P(32) = p3sq+s12+2*P(31)
       P(33) = p2sq-s12
       P(34) = p2sq-p5sq-s12-s15
       P(35) = p5sq+s12-s34
       P(36) = p3sq-s45+2*P(35)
       P(37) = s12+2*s45
       P(38) = p2sq+p5sq-s12-s15-2*s34
       P(39) = p5sq-s34
       P(40) = p3sq+3*s12+s45+2*P(39)
       P(41) = p3sq-s34+2*P(7)
       P(42) = -p5sq+s15+2*P(29)
       P(43) = 2*s12+s34
       P(44) = 1+2*Is12*s34
       P(45) = p2sq-s12-2*s15+s23-s45
       P(46) = p2sq-s12-s23-2*s34+s45
       P(47) = p2sq-s12-s15
       P(48) = s23-s45+2*P(47)
       P(49) = p2sq+3*s12+s23+2*s34-s45
       P(50) = s23-s45
       P(51) = p2sq-3*s12-4*s34-2*P(50)
       FR(4) = 2*(Dij2345R(1,1)-Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(2
     -   ,2)+4*(EijR(37,4)+EijR(38,4))-8*(EijR(22,3)-EijR(23,3)+EijR(
     -   43,4))-EijR(2,3)*P(5)+(-EijR(10,3)-EijR(13,3)+2*(EijR(2,1)-E
     -   ijR(3,1)+EijR(5,2)-EijR(6,2)+EijR(20,3)))*P(21)+(EijR(9,2)-E
     -   ijR(10,2))*P(23)+(EijR(2,2)+2*(EijR(8,3)+EijR(11,3)-2*EijR(1
     -   7,3)))*P(33)+Is34*(C0235R-Cij235R(1,1)+Cij235R(2,1)-4*(Dij12
     -   35R(7,2)-Dij1235R(12,3)+Dij1235R(13,3))+2*(Cij135R(1,2)+(D01
     -   235R+Dij1235R(1,1))*P(3))-Dij1235R(2,2)*P(5)+Is12*(4*Cij135R
     -   (4,2)+C0135R*P(26)+Cij135R(1,1)*P(27)-Cij135R(2,1)*P(28)+2*(
     -   B035R-Cij135R(3,2)*P(29)))+(Dij1235R(2,1)+2*(Dij1235R(4,2)-D
     -   ij1235R(5,2)))*P(33)-Dij1235R(3,1)*P(34)+Dij1235R(3,2)*P(38)
     -   +Dij1235R(6,2)*P(42))+Is12*(C0345R-D01345R*s45-Dij1345R(3,1)
     -   *P(7)+(Dij1345R(1,1)-Dij1345R(2,1))*P(32)+Dij1345R(1,2)*P(36
     -   )+Dij1345R(2,2)*P(37)-Dij1345R(4,2)*P(40)+(-Dij1345R(5,2)+Di
     -   j1345R(6,2))*P(41)+2*(Cij135R(2,2)-s34*Dij1345R(2,3)+Dij1345
     -   R(7,2)+4*(Dij1345R(11,3)-Dij1345R(12,3))-Dij1345R(4,3)*P(43)
     -   ))+2*(Dij1345R(1,3)+Dij2345R(4,2)+Dij1345R(6,3)*P(44))+EijR(
     -   3,2)*P(45)-EijR(3,3)*P(46)-EijR(8,2)*P(48)+EijR(9,3)*P(49)+E
     -   ijR(12,3)*P(51))
       FI(4) = 2*(Dij2345I(1,1)-Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(2
     -   ,2)+4*(EijI(37,4)+EijI(38,4))-8*(EijI(22,3)-EijI(23,3)+EijI(
     -   43,4))-EijI(2,3)*P(5)+(-EijI(10,3)-EijI(13,3)+2*(EijI(2,1)-E
     -   ijI(3,1)+EijI(5,2)-EijI(6,2)+EijI(20,3)))*P(21)+(EijI(9,2)-E
     -   ijI(10,2))*P(23)+(EijI(2,2)+2*(EijI(8,3)+EijI(11,3)-2*EijI(1
     -   7,3)))*P(33)+Is34*(C0235I-Cij235I(1,1)+Cij235I(2,1)-4*(Dij12
     -   35I(7,2)-Dij1235I(12,3)+Dij1235I(13,3))+2*(Cij135I(1,2)+(D01
     -   235I+Dij1235I(1,1))*P(3))-Dij1235I(2,2)*P(5)+Is12*(4*Cij135I
     -   (4,2)+C0135I*P(26)+Cij135I(1,1)*P(27)-Cij135I(2,1)*P(28)+2*(
     -   B035I-Cij135I(3,2)*P(29)))+(Dij1235I(2,1)+2*(Dij1235I(4,2)-D
     -   ij1235I(5,2)))*P(33)-Dij1235I(3,1)*P(34)+Dij1235I(3,2)*P(38)
     -   +Dij1235I(6,2)*P(42))+Is12*(C0345I-D01345I*s45-Dij1345I(3,1)
     -   *P(7)+(Dij1345I(1,1)-Dij1345I(2,1))*P(32)+Dij1345I(1,2)*P(36
     -   )+Dij1345I(2,2)*P(37)-Dij1345I(4,2)*P(40)+(-Dij1345I(5,2)+Di
     -   j1345I(6,2))*P(41)+2*(Cij135I(2,2)-s34*Dij1345I(2,3)+Dij1345
     -   I(7,2)+4*(Dij1345I(11,3)-Dij1345I(12,3))-Dij1345I(4,3)*P(43)
     -   ))+2*(Dij1345I(1,3)+Dij2345I(4,2)+Dij1345I(6,3)*P(44))+EijI(
     -   3,2)*P(45)-EijI(3,3)*P(46)-EijI(8,2)*P(48)+EijI(9,3)*P(49)+E
     -   ijI(12,3)*P(51))
       F(4)=DCMPLX(FR(4),FI(4))
       P(52) = p3sq+3*s12+2*P(31)
       P(53) = p3sq+3*s12-2*s34-s45
       P(54) = 3*p5sq-2*s45
       P(55) = p3sq+4*s12-s45+2*P(39)
       P(56) = p3sq-s12+2*s34+s45
       P(57) = 2*p3sq+3*s12
       P(58) = 4*p5sq+2*s12-s34-s45
       P(59) = p5sq+2*s12-s34
       P(60) = p2sq-3*s12
       P(61) = s15-s23+s45
       P(62) = p2sq-s15-s34
       P(63) = 4*s12+s23-s45-2*P(62)
       P(64) = p5sq-s15+3*s23-s45-2*P(6)
       P(65) = s23-s45+2*P(29)
       P(66) = s15-s23
       P(67) = p5sq-s45-2*P(6)-3*P(66)
       P(68) = p2sq-p5sq+s12+s15-s23+s45
       P(69) = p5sq-s15-2*P(29)
       FR(5) = 2*(-Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(2,1)+Dij2345R(
     -   4,2)+Dij2345R(5,2)-Dij2345R(6,2)-4*(EijR(22,3)-EijR(24,3)-Ei
     -   jR(37,4)+EijR(43,4)+EijR(44,4)-EijR(45,4))-EijR(2,3)*P(5)+(2
     -   *(EijR(2,1)+EijR(5,2))+EijR(15,3)-EijR(16,3))*P(21)+Is34*(C0
     -   235R-Cij235R(1,1)+Cij235R(2,1)-4*(Dij1235R(7,2)-Dij1235R(12,
     -   3)+Dij1235R(13,3))+2*(Cij135R(1,2)+(D01235R+Dij1235R(1,1))*P
     -   (3))-Dij1235R(2,2)*P(5)+Is12*(4*Cij135R(4,2)+C0135R*P(26)+Ci
     -   j135R(1,1)*P(27)-Cij135R(2,1)*P(28)+2*(B035R-Cij135R(3,2)*P(
     -   29)))+(Dij1235R(2,1)+2*(Dij1235R(4,2)-Dij1235R(5,2)))*P(33)-
     -   Dij1235R(3,1)*P(34)+Dij1235R(3,2)*P(38)+Dij1235R(6,2)*P(42))
     -   +(EijR(3,2)+EijR(12,3)-EijR(13,3))*P(46)+Is12*(C0345R-D01345
     -   R*s45+Cij345R(1,1)-Cij345R(2,1)+8*(Dij1345R(7,2)+Dij1345R(11
     -   ,3))-4*(Dij1345R(12,3)+Dij1345R(13,3))+Dij1345R(3,2)*P(7)+2*
     -   (Cij135R(2,2)+s34*(Dij1345R(6,3)-Dij1345R(8,3))+(-Dij1345R(4
     -   ,3)+Dij1345R(10,3))*P(29))+Dij1345R(1,1)*P(52)-Dij1345R(2,1)
     -   *P(53)-Dij1345R(3,1)*P(54)+Dij1345R(1,2)*P(55)+Dij1345R(2,2)
     -   *P(56)-Dij1345R(4,2)*P(57)-Dij1345R(5,2)*P(58)+Dij1345R(6,2)
     -   *P(59))+EijR(2,2)*P(60)+2*(Dij1345R(1,3)-Dij1345R(5,3)-p5sq*
     -   EijR(4,1)-EijR(7,2)*P(3)+(EijR(8,3)-EijR(17,3)-EijR(18,3)+Ei
     -   jR(19,3))*P(33)-EijR(6,2)*P(50)+EijR(3,1)*P(61))+EijR(8,2)*P
     -   (63)-EijR(9,2)*P(64)+EijR(9,3)*P(65)+EijR(10,2)*P(67)+EijR(1
     -   0,3)*P(68)+EijR(20,3)*P(69))
       FI(5) = 2*(-Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(2,1)+Dij2345I(
     -   4,2)+Dij2345I(5,2)-Dij2345I(6,2)-4*(EijI(22,3)-EijI(24,3)-Ei
     -   jI(37,4)+EijI(43,4)+EijI(44,4)-EijI(45,4))-EijI(2,3)*P(5)+(2
     -   *(EijI(2,1)+EijI(5,2))+EijI(15,3)-EijI(16,3))*P(21)+Is34*(C0
     -   235I-Cij235I(1,1)+Cij235I(2,1)-4*(Dij1235I(7,2)-Dij1235I(12,
     -   3)+Dij1235I(13,3))+2*(Cij135I(1,2)+(D01235I+Dij1235I(1,1))*P
     -   (3))-Dij1235I(2,2)*P(5)+Is12*(4*Cij135I(4,2)+C0135I*P(26)+Ci
     -   j135I(1,1)*P(27)-Cij135I(2,1)*P(28)+2*(B035I-Cij135I(3,2)*P(
     -   29)))+(Dij1235I(2,1)+2*(Dij1235I(4,2)-Dij1235I(5,2)))*P(33)-
     -   Dij1235I(3,1)*P(34)+Dij1235I(3,2)*P(38)+Dij1235I(6,2)*P(42))
     -   +(EijI(3,2)+EijI(12,3)-EijI(13,3))*P(46)+Is12*(C0345I-D01345
     -   I*s45+Cij345I(1,1)-Cij345I(2,1)+8*(Dij1345I(7,2)+Dij1345I(11
     -   ,3))-4*(Dij1345I(12,3)+Dij1345I(13,3))+Dij1345I(3,2)*P(7)+2*
     -   (Cij135I(2,2)+s34*(Dij1345I(6,3)-Dij1345I(8,3))+(-Dij1345I(4
     -   ,3)+Dij1345I(10,3))*P(29))+Dij1345I(1,1)*P(52)-Dij1345I(2,1)
     -   *P(53)-Dij1345I(3,1)*P(54)+Dij1345I(1,2)*P(55)+Dij1345I(2,2)
     -   *P(56)-Dij1345I(4,2)*P(57)-Dij1345I(5,2)*P(58)+Dij1345I(6,2)
     -   *P(59))+EijI(2,2)*P(60)+2*(Dij1345I(1,3)-Dij1345I(5,3)-p5sq*
     -   EijI(4,1)-EijI(7,2)*P(3)+(EijI(8,3)-EijI(17,3)-EijI(18,3)+Ei
     -   jI(19,3))*P(33)-EijI(6,2)*P(50)+EijI(3,1)*P(61))+EijI(8,2)*P
     -   (63)-EijI(9,2)*P(64)+EijI(9,3)*P(65)+EijI(10,2)*P(67)+EijI(1
     -   0,3)*P(68)+EijI(20,3)*P(69))
       F(5)=DCMPLX(FR(5),FI(5))
       P(70) = -2*p5sq+s15
       P(71) = p2sq-2*s12+s34
       P(72) = p2sq-2*s12
       P(73) = p5sq-s12+2*s34
       P(74) = p2sq-s15
       P(75) = p5sq-s12-s34+2*P(74)
       P(76) = -p5sq+s12+s34-2*P(74)
       P(77) = s23-2*s45
       P(78) = p3sq+s12+s23-2*P(30)
       P(79) = p3sq+s12-2*s34-s45
       P(80) = p5sq-2*s45
       P(81) = p2sq+s34+s45
       P(82) = -p3sq-5*s12-s23+2*P(81)
       P(83) = s23+s34
       P(84) = p2sq-P(83)
       P(85) = -p3sq-s12+s45+2*P(84)
       P(86) = p5sq-2*s15+3*s23-4*P(6)
       P(87) = p2sq-s23-2*s34
       P(88) = s12-s23
       P(89) = s15-s34
       P(90) = p3sq-3*P(88)-2*P(89)
       P(91) = p2sq-s12-s23-s34+s45
       P(92) = p5sq-s12-s15+s23-s45
       P(93) = p2sq-p5sq+s15-s23+s45
       FR(6) = Is15s34*(-4*(B025R+2*Cij125R(4,2))+2*(C0125R+Cij125R(1,
     -   1))*P(70))+2*(EE0R*s45-4*Dij2345R(1,1)+3*Dij2345R(3,1)-EijR(
     -   4,2)*P(7)-EijR(1,2)*P(50)+Is34*(Cij125R(2,1)-2*(C0235R+Cij12
     -   5R(2,2)-Cij125R(3,2))-Cij235R(1,1)+4*(Dij1235R(7,2)-Dij1235R
     -   (11,3)+Dij1235R(12,3))+D01235R*P(6)+(-Dij1235R(1,2)+Dij1235R
     -   (2,2))*P(33)-Dij1235R(1,1)*P(71)+Dij1235R(2,1)*P(72)+Dij1235
     -   R(3,1)*P(73)+Dij1235R(5,2)*P(75)+Dij1235R(6,2)*P(76))-EijR(1
     -   ,1)*P(77)+EijR(2,1)*P(78)-(EijR(3,1)-EijR(10,2))*P(79)+EijR(
     -   4,1)*P(80)-EijR(5,2)*P(82)+EijR(6,2)*P(85)+EijR(7,2)*P(86)-E
     -   ijR(9,2)*P(90)-2*(D02345R+Dij2345R(1,2)-Dij2345R(5,2)+s12*Ei
     -   jR(2,3)-p2sq*EijR(8,3)-EijR(11,2)-2*(EijR(22,3)-EijR(24,3)+E
     -   ijR(37,4)-EijR(40,4)+EijR(42,4)-EijR(44,4))+(EijR(14,3)-EijR
     -   (15,3))*P(21)+(EijR(5,3)-EijR(7,3))*P(33)-EijR(2,2)*P(72)+Ei
     -   jR(8,2)*P(87)+(EijR(9,3)-EijR(17,3)+EijR(19,3)-EijR(20,3))*P
     -   (91)+EijR(10,3)*P(92)+EijR(18,3)*P(93)))
       FI(6) = Is15s34*(-4*(B025I+2*Cij125I(4,2))+2*(C0125I+Cij125I(1,
     -   1))*P(70))+2*(EE0I*s45-4*Dij2345I(1,1)+3*Dij2345I(3,1)-EijI(
     -   4,2)*P(7)-EijI(1,2)*P(50)+Is34*(Cij125I(2,1)-2*(C0235I+Cij12
     -   5I(2,2)-Cij125I(3,2))-Cij235I(1,1)+4*(Dij1235I(7,2)-Dij1235I
     -   (11,3)+Dij1235I(12,3))+D01235I*P(6)+(-Dij1235I(1,2)+Dij1235I
     -   (2,2))*P(33)-Dij1235I(1,1)*P(71)+Dij1235I(2,1)*P(72)+Dij1235
     -   I(3,1)*P(73)+Dij1235I(5,2)*P(75)+Dij1235I(6,2)*P(76))-EijI(1
     -   ,1)*P(77)+EijI(2,1)*P(78)-(EijI(3,1)-EijI(10,2))*P(79)+EijI(
     -   4,1)*P(80)-EijI(5,2)*P(82)+EijI(6,2)*P(85)+EijI(7,2)*P(86)-E
     -   ijI(9,2)*P(90)-2*(D02345I+Dij2345I(1,2)-Dij2345I(5,2)+s12*Ei
     -   jI(2,3)-p2sq*EijI(8,3)-EijI(11,2)-2*(EijI(22,3)-EijI(24,3)+E
     -   ijI(37,4)-EijI(40,4)+EijI(42,4)-EijI(44,4))+(EijI(14,3)-EijI
     -   (15,3))*P(21)+(EijI(5,3)-EijI(7,3))*P(33)-EijI(2,2)*P(72)+Ei
     -   jI(8,2)*P(87)+(EijI(9,3)-EijI(17,3)+EijI(19,3)-EijI(20,3))*P
     -   (91)+EijI(10,3)*P(92)+EijI(18,3)*P(93)))
       F(6)=DCMPLX(FR(6),FI(6))
       P(94) = p5sq-s12+2*P(74)
       FR(7) = -2*Is15s34*(2*(B025R+C0235R*s15+p5sq*Cij125R(1,1))+4*Ci
     -   j125R(4,2)-s15*(Cij125R(1,1)+Cij125R(2,1)-2*(Cij125R(2,2)-Ci
     -   j125R(3,2))-Cij235R(1,1)+s12*(Dij1235R(1,2)-2*Dij1235R(2,1)-
     -   Dij1235R(2,2)-Dij1235R(3,1))+p5sq*Dij1235R(3,1)+(p5sq-s12-2*
     -   s15)*Dij1235R(5,2)+p2sq*(-Dij1235R(1,2)+Dij1235R(2,1)+Dij123
     -   5R(2,2)+2*Dij1235R(5,2))+4*(Dij1235R(7,2)-Dij1235R(11,3)+Dij
     -   1235R(12,3))+s34*(-Dij1235R(5,2)+Dij1235R(6,2)+3*Dij2345R(2,
     -   1)+p3sq*(EijR(2,1)-EijR(3,1)+EijR(3,2)+EijR(5,2)-EijR(6,2)-E
     -   ijR(8,2))+s12*(EijR(2,1)-EijR(3,1)+EijR(3,2)+3*(EijR(5,2)-Ei
     -   jR(6,2))+EijR(8,2)+4*EijR(9,3)-2*(EijR(2,2)+EijR(2,3)-EijR(5
     -   ,3)+EijR(6,3)-EijR(11,3)+EijR(12,3)+EijR(17,3)))+2*(Dij1235R
     -   (3,1)-Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(4,2)+p2sq*(EijR(2
     -   ,2)-EijR(5,2)-EijR(5,3)+EijR(6,2)+EijR(6,3)-EijR(8,2)+EijR(8
     -   ,3)-EijR(9,3)-EijR(11,3)+EijR(12,3))+s34*(-EijR(2,1)+EijR(3,
     -   1)-EijR(3,2)-EijR(5,2)+EijR(6,2)+EijR(8,2)+EijR(9,3)+EijR(11
     -   ,3)-EijR(12,3)-EijR(17,3))+s15*(-EijR(7,2)+EijR(9,2)+EijR(10
     -   ,3)-EijR(18,3)+EijR(19,3)-EijR(20,3)))+s45*(EE0R+EijR(1,2)-E
     -   ijR(3,2)-EijR(4,1)+EijR(7,2)+EijR(10,2)+2*(EijR(1,1)-EijR(2,
     -   1)-EijR(5,2)+EijR(8,2)-EijR(9,2)-EijR(9,3)+EijR(10,3)-EijR(1
     -   1,3)+EijR(12,3)+EijR(17,3)-EijR(18,3)+EijR(19,3)-EijR(20,3))
     -   )+p5sq*(EijR(4,1)-EijR(7,2)-EijR(10,2)+2*(EijR(9,2)-EijR(10,
     -   3)+EijR(18,3)-EijR(19,3)+EijR(20,3)))+s23*(-EijR(1,1)-EijR(1
     -   ,2)+EijR(2,1)+EijR(5,2)-EijR(6,2)+EijR(8,2)+2*(EijR(7,2)-Eij
     -   R(9,2)+EijR(9,3)-EijR(10,3)+EijR(11,3)-EijR(12,3)-EijR(17,3)
     -   +EijR(18,3)-EijR(19,3)+EijR(20,3)))+4*(EijR(11,2)+EijR(21,3)
     -   -EijR(23,3)+EijR(37,4)-EijR(40,4)+EijR(41,4)-EijR(43,4))))-D
     -   01235R*s15*P(6)-C0125R*P(70)+s15*Dij1235R(1,1)*P(71)+s15*Dij
     -   1235R(6,2)*P(94))
       FI(7) = -2*Is15s34*(2*(B025I+C0235I*s15+p5sq*Cij125I(1,1))+4*Ci
     -   j125I(4,2)-s15*(Cij125I(1,1)+Cij125I(2,1)-2*(Cij125I(2,2)-Ci
     -   j125I(3,2))-Cij235I(1,1)+s12*(Dij1235I(1,2)-2*Dij1235I(2,1)-
     -   Dij1235I(2,2)-Dij1235I(3,1))+p5sq*Dij1235I(3,1)+(p5sq-s12-2*
     -   s15)*Dij1235I(5,2)+p2sq*(-Dij1235I(1,2)+Dij1235I(2,1)+Dij123
     -   5I(2,2)+2*Dij1235I(5,2))+4*(Dij1235I(7,2)-Dij1235I(11,3)+Dij
     -   1235I(12,3))+s34*(-Dij1235I(5,2)+Dij1235I(6,2)+3*Dij2345I(2,
     -   1)+p3sq*(EijI(2,1)-EijI(3,1)+EijI(3,2)+EijI(5,2)-EijI(6,2)-E
     -   ijI(8,2))+s12*(EijI(2,1)-EijI(3,1)+EijI(3,2)+3*(EijI(5,2)-Ei
     -   jI(6,2))+EijI(8,2)+4*EijI(9,3)-2*(EijI(2,2)+EijI(2,3)-EijI(5
     -   ,3)+EijI(6,3)-EijI(11,3)+EijI(12,3)+EijI(17,3)))+2*(Dij1235I
     -   (3,1)-Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(4,2)+p2sq*(EijI(2
     -   ,2)-EijI(5,2)-EijI(5,3)+EijI(6,2)+EijI(6,3)-EijI(8,2)+EijI(8
     -   ,3)-EijI(9,3)-EijI(11,3)+EijI(12,3))+s34*(-EijI(2,1)+EijI(3,
     -   1)-EijI(3,2)-EijI(5,2)+EijI(6,2)+EijI(8,2)+EijI(9,3)+EijI(11
     -   ,3)-EijI(12,3)-EijI(17,3))+s15*(-EijI(7,2)+EijI(9,2)+EijI(10
     -   ,3)-EijI(18,3)+EijI(19,3)-EijI(20,3)))+s45*(EE0I+EijI(1,2)-E
     -   ijI(3,2)-EijI(4,1)+EijI(7,2)+EijI(10,2)+2*(EijI(1,1)-EijI(2,
     -   1)-EijI(5,2)+EijI(8,2)-EijI(9,2)-EijI(9,3)+EijI(10,3)-EijI(1
     -   1,3)+EijI(12,3)+EijI(17,3)-EijI(18,3)+EijI(19,3)-EijI(20,3))
     -   )+p5sq*(EijI(4,1)-EijI(7,2)-EijI(10,2)+2*(EijI(9,2)-EijI(10,
     -   3)+EijI(18,3)-EijI(19,3)+EijI(20,3)))+s23*(-EijI(1,1)-EijI(1
     -   ,2)+EijI(2,1)+EijI(5,2)-EijI(6,2)+EijI(8,2)+2*(EijI(7,2)-Eij
     -   I(9,2)+EijI(9,3)-EijI(10,3)+EijI(11,3)-EijI(12,3)-EijI(17,3)
     -   +EijI(18,3)-EijI(19,3)+EijI(20,3)))+4*(EijI(11,2)+EijI(21,3)
     -   -EijI(23,3)+EijI(37,4)-EijI(40,4)+EijI(41,4)-EijI(43,4))))-D
     -   01235I*s15*P(6)-C0125I*P(70)+s15*Dij1235I(1,1)*P(71)+s15*Dij
     -   1235I(6,2)*P(94))
       F(7)=DCMPLX(FR(7),FI(7))
       P(95) = 2*p5sq+s45
       P(96) = 2*p5sq+5*s45
       P(97) = -s45+2*P(35)
       P(98) = p3sq-2*s45+4*P(35)
       P(99) = 2*p5sq+s12-s45
       P(100) = p3sq-3*s34+2*P(16)
       P(101) = p3sq+p5sq-s34-s45
       P(102) = p5sq+s34
       P(103) = p3sq-3*s45-2*P(102)
       P(104) = p3sq+s12
       P(105) = 2*p5sq-5*s34-4*s45+3*P(104)
       P(106) = p3sq+p5sq+2*P(6)
       P(107) = p5sq-P(61)
       P(108) = p2sq-s12-2*P(107)
       P(109) = p2sq-s12-s23+s45-2*P(3)
       P(110) = p5sq+s15-s23+s45
       P(111) = p2sq-s12-2*P(50)
       P(112) = p5sq+s15-s45
       P(113) = p2sq-s12-2*P(112)
       P(114) = -p2sq+p5sq+s12+s15
       P(115) = -p5sq+s15-s23+s45
       FR(8) = 2*(D02345R-Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(3,1)+Di
     -   j2345R(4,2)+Dij2345R(5,2)-Dij2345R(6,2)-8*EijR(22,3)-4*(EijR
     -   (11,2)-EijR(23,3)-EijR(24,3)-EijR(37,4)+EijR(43,4)+EijR(44,4
     -   )-EijR(45,4))-EijR(2,3)*P(5)-EijR(3,2)*P(10)+(2*EijR(5,2)+Ei
     -   jR(15,3))*P(21)-EijR(4,2)*P(23)+Is34*(C0235R-Cij235R(1,1)+Ci
     -   j235R(2,1)-4*(Dij1235R(7,2)-Dij1235R(12,3)+Dij1235R(13,3))+2
     -   *(Cij135R(1,2)+(D01235R+Dij1235R(1,1))*P(3))-Dij1235R(2,2)*P
     -   (5)+Is12*(4*Cij135R(4,2)+C0135R*P(26)+Cij135R(1,1)*P(27)-Cij
     -   135R(2,1)*P(28)+2*(B035R-Cij135R(3,2)*P(29)))+(Dij1235R(2,1)
     -   +2*(Dij1235R(4,2)-Dij1235R(5,2)))*P(33)-Dij1235R(3,1)*P(34)+
     -   Dij1235R(3,2)*P(38)+Dij1235R(6,2)*P(42))+(EijR(12,3)-EijR(13
     -   ,3))*P(46)+2*(Dij1345R(1,3)-Dij1345R(5,3)+Dij2345R(1,1)+(EE0
     -   R+EijR(1,1)-EijR(6,2))*P(3)+(EijR(8,3)-EijR(17,3)-EijR(18,3)
     -   +EijR(19,3))*P(33)-EijR(7,2)*P(50))+EijR(9,3)*P(65)+EijR(10,
     -   3)*P(68)+EijR(20,3)*P(69)+Is12*(-3*Cij145R(2,1)-Cij345R(1,1)
     -   +Cij345R(2,1)+2*(C0345R+Cij135R(2,2)+Cij145R(1,2)-Cij145R(3,
     -   2)-Dij1345R(7,2)-s34*(-Dij1345R(6,3)+Dij1345R(8,3))-2*(-2*Di
     -   j1345R(11,3)+Dij1345R(12,3)+Dij1345R(13,3))+(-Dij1345R(4,3)+
     -   Dij1345R(10,3))*P(29))+Dij1345R(1,2)*P(36)-Dij1345R(2,2)*P(7
     -   9)+Is45*(2*B012R+4*Cij145R(4,2)+C0145R*P(95)+Cij145R(1,1)*P(
     -   96))+D01345R*P(97)+Dij1345R(1,1)*P(98)-Dij1345R(2,1)*P(99)-D
     -   ij1345R(3,1)*P(100)+Dij1345R(3,2)*P(101)+Dij1345R(4,2)*P(103
     -   )-Dij1345R(5,2)*P(105)+Dij1345R(6,2)*P(106))-EijR(2,1)*P(108
     -   )+EijR(3,1)*P(109)+EijR(4,1)*P(110)+EijR(8,2)*P(111)-EijR(9,
     -   2)*P(113)-EijR(10,2)*P(114)+EijR(16,3)*P(115))
       FI(8) = 2*(D02345I-Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(3,1)+Di
     -   j2345I(4,2)+Dij2345I(5,2)-Dij2345I(6,2)-8*EijI(22,3)-4*(EijI
     -   (11,2)-EijI(23,3)-EijI(24,3)-EijI(37,4)+EijI(43,4)+EijI(44,4
     -   )-EijI(45,4))-EijI(2,3)*P(5)-EijI(3,2)*P(10)+(2*EijI(5,2)+Ei
     -   jI(15,3))*P(21)-EijI(4,2)*P(23)+Is34*(C0235I-Cij235I(1,1)+Ci
     -   j235I(2,1)-4*(Dij1235I(7,2)-Dij1235I(12,3)+Dij1235I(13,3))+2
     -   *(Cij135I(1,2)+(D01235I+Dij1235I(1,1))*P(3))-Dij1235I(2,2)*P
     -   (5)+Is12*(4*Cij135I(4,2)+C0135I*P(26)+Cij135I(1,1)*P(27)-Cij
     -   135I(2,1)*P(28)+2*(B035I-Cij135I(3,2)*P(29)))+(Dij1235I(2,1)
     -   +2*(Dij1235I(4,2)-Dij1235I(5,2)))*P(33)-Dij1235I(3,1)*P(34)+
     -   Dij1235I(3,2)*P(38)+Dij1235I(6,2)*P(42))+(EijI(12,3)-EijI(13
     -   ,3))*P(46)+2*(Dij1345I(1,3)-Dij1345I(5,3)+Dij2345I(1,1)+(EE0
     -   I+EijI(1,1)-EijI(6,2))*P(3)+(EijI(8,3)-EijI(17,3)-EijI(18,3)
     -   +EijI(19,3))*P(33)-EijI(7,2)*P(50))+EijI(9,3)*P(65)+EijI(10,
     -   3)*P(68)+EijI(20,3)*P(69)+Is12*(-3*Cij145I(2,1)-Cij345I(1,1)
     -   +Cij345I(2,1)+2*(C0345I+Cij135I(2,2)+Cij145I(1,2)-Cij145I(3,
     -   2)-Dij1345I(7,2)-s34*(-Dij1345I(6,3)+Dij1345I(8,3))-2*(-2*Di
     -   j1345I(11,3)+Dij1345I(12,3)+Dij1345I(13,3))+(-Dij1345I(4,3)+
     -   Dij1345I(10,3))*P(29))+Dij1345I(1,2)*P(36)-Dij1345I(2,2)*P(7
     -   9)+Is45*(2*B012I+4*Cij145I(4,2)+C0145I*P(95)+Cij145I(1,1)*P(
     -   96))+D01345I*P(97)+Dij1345I(1,1)*P(98)-Dij1345I(2,1)*P(99)-D
     -   ij1345I(3,1)*P(100)+Dij1345I(3,2)*P(101)+Dij1345I(4,2)*P(103
     -   )-Dij1345I(5,2)*P(105)+Dij1345I(6,2)*P(106))-EijI(2,1)*P(108
     -   )+EijI(3,1)*P(109)+EijI(4,1)*P(110)+EijI(8,2)*P(111)-EijI(9,
     -   2)*P(113)-EijI(10,2)*P(114)+EijI(16,3)*P(115))
       F(8)=DCMPLX(FR(8),FI(8))
       P(116) = p3sq+6*s12-2*s45+4*P(39)
       P(117) = p3sq+2*p5sq-3*s45+5*P(6)
       P(118) = p3sq+3*P(6)+2*P(7)
       P(119) = s12-2*s34-s45
       P(120) = 2*p3sq+4*p5sq+7*s12-5*s34-3*s45
       P(121) = -p2sq+s12+2*P(107)
       P(122) = p5sq-3*s15+s23-s45
       P(123) = s12+s15-s23-s34
       P(124) = 2*s12+s15-s23-s34
       P(125) = -p5sq+s15-s23+s45+2*P(5)
       P(126) = p2sq+s12-2*P(107)
       FR(9) = 2*(D02345R-Dij2345R(1,2)-Dij2345R(3,2)-4*(Dij1345R(5,3)
     -   +EijR(11,2)+EijR(22,3)-EijR(24,3)-EijR(37,4)-EijR(39,4))-8*E
     -   ijR(44,4)-EijR(2,3)*P(5)+EijR(3,1)*P(10)+Is34*(C0235R-Cij235
     -   R(1,1)+Cij235R(2,1)-4*(Dij1235R(7,2)-Dij1235R(12,3)+Dij1235R
     -   (13,3))+2*(Cij135R(1,2)+(D01235R+Dij1235R(1,1))*P(3))-Dij123
     -   5R(2,2)*P(5)+Is12*(4*Cij135R(4,2)+C0135R*P(26)+Cij135R(1,1)*
     -   P(27)-Cij135R(2,1)*P(28)+2*(B035R-Cij135R(3,2)*P(29)))+(Dij1
     -   235R(2,1)+2*(Dij1235R(4,2)-Dij1235R(5,2)))*P(33)-Dij1235R(3,
     -   1)*P(34)+Dij1235R(3,2)*P(38)+Dij1235R(6,2)*P(42))-(EijR(9,3)
     -   +EijR(16,3)-2*EijR(20,3))*P(46)+EijR(4,3)*P(115)+Is12*(-3*Ci
     -   j145R(2,1)+Dij1345R(1,2)*P(55)+Is45*(2*B012R+4*Cij145R(4,2)+
     -   C0145R*P(95)+Cij145R(1,1)*P(96))+D01345R*P(97)+Dij1345R(1,1)
     -   *P(116)-Dij1345R(3,1)*P(117)+Dij1345R(3,2)*P(118)+2*(C0345R+
     -   Cij135R(2,2)+Cij145R(1,2)-Cij145R(3,2)-s34*(Dij1345R(2,1)+Di
     -   j1345R(4,3)+Dij1345R(9,3)-2*Dij1345R(10,3))+2*(Dij1345R(7,2)
     -   +2*(Dij1345R(11,3)-Dij1345R(13,3)))+(Dij1345R(4,2)-Dij1345R(
     -   6,2))*P(119))-Dij1345R(5,2)*P(120))+EijR(2,1)*P(121)-EijR(4,
     -   1)*P(122)+2*(Dij1345R(1,3)+Dij1345R(7,3)+Dij2345R(5,2)-s12*E
     -   ijR(2,2)+s34*EijR(8,2)-s34*EijR(10,2)+(EE0R+EijR(1,1))*P(3)+
     -   (EijR(5,2)-EijR(7,2))*P(21)+(EijR(8,3)+EijR(14,3)-2*EijR(18,
     -   3))*P(33)-EijR(4,2)*P(123)+EijR(9,2)*P(124))+EijR(10,3)*P(12
     -   5)-EijR(15,3)*P(126))
       FI(9) = 2*(D02345I-Dij2345I(1,2)-Dij2345I(3,2)-4*(Dij1345I(5,3)
     -   +EijI(11,2)+EijI(22,3)-EijI(24,3)-EijI(37,4)-EijI(39,4))-8*E
     -   ijI(44,4)-EijI(2,3)*P(5)+EijI(3,1)*P(10)+Is34*(C0235I-Cij235
     -   I(1,1)+Cij235I(2,1)-4*(Dij1235I(7,2)-Dij1235I(12,3)+Dij1235I
     -   (13,3))+2*(Cij135I(1,2)+(D01235I+Dij1235I(1,1))*P(3))-Dij123
     -   5I(2,2)*P(5)+Is12*(4*Cij135I(4,2)+C0135I*P(26)+Cij135I(1,1)*
     -   P(27)-Cij135I(2,1)*P(28)+2*(B035I-Cij135I(3,2)*P(29)))+(Dij1
     -   235I(2,1)+2*(Dij1235I(4,2)-Dij1235I(5,2)))*P(33)-Dij1235I(3,
     -   1)*P(34)+Dij1235I(3,2)*P(38)+Dij1235I(6,2)*P(42))-(EijI(9,3)
     -   +EijI(16,3)-2*EijI(20,3))*P(46)+EijI(4,3)*P(115)+Is12*(-3*Ci
     -   j145I(2,1)+Dij1345I(1,2)*P(55)+Is45*(2*B012I+4*Cij145I(4,2)+
     -   C0145I*P(95)+Cij145I(1,1)*P(96))+D01345I*P(97)+Dij1345I(1,1)
     -   *P(116)-Dij1345I(3,1)*P(117)+Dij1345I(3,2)*P(118)+2*(C0345I+
     -   Cij135I(2,2)+Cij145I(1,2)-Cij145I(3,2)-s34*(Dij1345I(2,1)+Di
     -   j1345I(4,3)+Dij1345I(9,3)-2*Dij1345I(10,3))+2*(Dij1345I(7,2)
     -   +2*(Dij1345I(11,3)-Dij1345I(13,3)))+(Dij1345I(4,2)-Dij1345I(
     -   6,2))*P(119))-Dij1345I(5,2)*P(120))+EijI(2,1)*P(121)-EijI(4,
     -   1)*P(122)+2*(Dij1345I(1,3)+Dij1345I(7,3)+Dij2345I(5,2)-s12*E
     -   ijI(2,2)+s34*EijI(8,2)-s34*EijI(10,2)+(EE0I+EijI(1,1))*P(3)+
     -   (EijI(5,2)-EijI(7,2))*P(21)+(EijI(8,3)+EijI(14,3)-2*EijI(18,
     -   3))*P(33)-EijI(4,2)*P(123)+EijI(9,2)*P(124))+EijI(10,3)*P(12
     -   5)-EijI(15,3)*P(126))
       F(9)=DCMPLX(FR(9),FI(9))
       P(127) = s23-s45+2*P(3)
       P(128) = 2*p5sq-s34
       P(129) = p5sq+2*s45
       P(130) = 2*p5sq+4*s12-s34
       P(131) = p5sq+s12-s45
       P(132) = p3sq-3*s34+2*P(131)
       P(133) = p2sq-s12+2*P(3)
       P(134) = p5sq+s12
       P(135) = 2*p3sq-5*s34-3*s45+4*P(134)
       P(136) = s12+2*P(7)
       P(137) = p5sq-s12
       P(138) = s34+2*P(137)
       P(139) = p3sq+s12-s45+2*P(39)
       P(140) = p3sq+s12-s34-s45
       P(141) = p3sq-2*p5sq+s12-s34-3*s45
       P(142) = -p2sq-p5sq+s12+s15+2*s34
       P(143) = 2*p2sq+s12
       P(144) = s23-3*s45+2*P(3)
       P(145) = 3*P(33)+2*P(107)
       P(146) = 2*P(3)+3*P(33)
       P(147) = 2*p2sq+s12-s23+3*s34+s45
       FR(10) = Is12*(4*Is45*(B012R+C0145R*p5sq+2*Cij145R(4,2)+Cij145R
     -   (1,1)*P(129))+2*(3*C0345R-Cij345R(1,1)-4*(Dij1345R(7,2)-2*Di
     -   j1345R(11,3)+Dij1345R(12,3))+(Dij1345R(3,1)+Dij1345R(6,2))*P
     -   (7)-2*(Cij135R(2,1)+Cij135R(3,2)-Cij145R(1,2)-s34*Dij1345R(6
     -   ,3)+Dij1345R(4,3)*P(29))-Dij1345R(5,2)*P(41)+D01345R*P(132)+
     -   Dij1345R(1,1)*P(135)-Dij1345R(2,1)*P(136)+Dij1345R(1,2)*P(13
     -   9)-Dij1345R(2,2)*P(140)+Dij1345R(4,2)*P(141)))+2*(Dij2345R(1
     -   ,1)-Dij2345R(1,2)+Dij2345R(4,2)-14*EijR(11,2)-10*EijR(22,3)+
     -   6*EijR(23,3)+4*(EijR(37,4)-EijR(43,4))-EijR(2,3)*P(5)+3*(D02
     -   345R-s34*EijR(3,2)+EijR(4,1)*P(7))+(2*EijR(2,1)-EijR(10,3)+E
     -   ijR(20,3))*P(21)+EijR(9,2)*P(23)+2*(Dij1345R(1,3)+(EijR(8,3)
     -   -EijR(17,3))*P(33))+EijR(12,3)*P(46)+EijR(9,3)*P(65)+(EE0R+E
     -   ijR(1,1))*P(127)+Is34*(-Cij125R(2,1)-Cij235R(1,1)+p2sq*Dij12
     -   35R(2,1)-6*Dij1235R(7,2)+4*Dij1235R(12,3)-Dij1235R(2,2)*P(5)
     -   +2*(C0235R-Cij125R(2,2)+Cij135R(1,2)+Dij1235R(4,2)*P(33))+Is
     -   12*(2*B035R+4*Cij135R(4,2)+C0135R*P(128)+Cij135R(1,1)*P(130)
     -   )+(D01235R+Dij1235R(1,1))*P(133)+Dij1235R(3,1)*P(138)+Dij123
     -   5R(6,2)*P(142))-EijR(2,2)*P(143)-EijR(3,1)*P(144)+EijR(5,2)*
     -   P(145)-EijR(6,2)*P(146)+EijR(8,2)*P(147))
       FI(10) = Is12*(4*Is45*(B012I+C0145I*p5sq+2*Cij145I(4,2)+Cij145I
     -   (1,1)*P(129))+2*(3*C0345I-Cij345I(1,1)-4*(Dij1345I(7,2)-2*Di
     -   j1345I(11,3)+Dij1345I(12,3))+(Dij1345I(3,1)+Dij1345I(6,2))*P
     -   (7)-2*(Cij135I(2,1)+Cij135I(3,2)-Cij145I(1,2)-s34*Dij1345I(6
     -   ,3)+Dij1345I(4,3)*P(29))-Dij1345I(5,2)*P(41)+D01345I*P(132)+
     -   Dij1345I(1,1)*P(135)-Dij1345I(2,1)*P(136)+Dij1345I(1,2)*P(13
     -   9)-Dij1345I(2,2)*P(140)+Dij1345I(4,2)*P(141)))+2*(Dij2345I(1
     -   ,1)-Dij2345I(1,2)+Dij2345I(4,2)-14*EijI(11,2)-10*EijI(22,3)+
     -   6*EijI(23,3)+4*(EijI(37,4)-EijI(43,4))-EijI(2,3)*P(5)+3*(D02
     -   345I-s34*EijI(3,2)+EijI(4,1)*P(7))+(2*EijI(2,1)-EijI(10,3)+E
     -   ijI(20,3))*P(21)+EijI(9,2)*P(23)+2*(Dij1345I(1,3)+(EijI(8,3)
     -   -EijI(17,3))*P(33))+EijI(12,3)*P(46)+EijI(9,3)*P(65)+(EE0I+E
     -   ijI(1,1))*P(127)+Is34*(-Cij125I(2,1)-Cij235I(1,1)+p2sq*Dij12
     -   35I(2,1)-6*Dij1235I(7,2)+4*Dij1235I(12,3)-Dij1235I(2,2)*P(5)
     -   +2*(C0235I-Cij125I(2,2)+Cij135I(1,2)+Dij1235I(4,2)*P(33))+Is
     -   12*(2*B035I+4*Cij135I(4,2)+C0135I*P(128)+Cij135I(1,1)*P(130)
     -   )+(D01235I+Dij1235I(1,1))*P(133)+Dij1235I(3,1)*P(138)+Dij123
     -   5I(6,2)*P(142))-EijI(2,2)*P(143)-EijI(3,1)*P(144)+EijI(5,2)*
     -   P(145)-EijI(6,2)*P(146)+EijI(8,2)*P(147))
       F(10)=DCMPLX(FR(10),FI(10))
       P(148) = -s12+s34
       P(149) = s12-s34+2*P(3)
       P(150) = -3*s12+2*P(8)
       P(151) = -5*s12+2*P(8)
       P(152) = s12-s45+3*P(14)
       P(153) = p5sq-s45+3*P(14)
       FR(11) = 2*(Dij2345R(1,2)+s12*EijR(2,2)-4*(EijR(11,2)-EijR(37,4
     -   )+EijR(40,4))+2*(D02345R+Dij2345R(1,1)+EijR(21,3)-EijR(1,2)*
     -   P(3))+(EijR(2,3)-EijR(8,3))*P(5)-EijR(4,1)*P(7)+(-EijR(9,3)+
     -   EijR(17,3))*P(15)+EijR(3,1)*P(16)+(EijR(10,3)-EijR(18,3))*P(
     -   18)-3*(EijR(7,2)*P(101)-EijR(6,2)*P(140))+EE0R*P(148)-EijR(1
     -   ,1)*P(149)+EijR(2,1)*P(150)+EijR(5,2)*P(151)-EijR(8,2)*P(152
     -   )+EijR(9,2)*P(153))
       FI(11) = 2*(Dij2345I(1,2)+s12*EijI(2,2)-4*(EijI(11,2)-EijI(37,4
     -   )+EijI(40,4))+2*(D02345I+Dij2345I(1,1)+EijI(21,3)-EijI(1,2)*
     -   P(3))+(EijI(2,3)-EijI(8,3))*P(5)-EijI(4,1)*P(7)+(-EijI(9,3)+
     -   EijI(17,3))*P(15)+EijI(3,1)*P(16)+(EijI(10,3)-EijI(18,3))*P(
     -   18)-3*(EijI(7,2)*P(101)-EijI(6,2)*P(140))+EE0I*P(148)-EijI(1
     -   ,1)*P(149)+EijI(2,1)*P(150)+EijI(5,2)*P(151)-EijI(8,2)*P(152
     -   )+EijI(9,2)*P(153))
       F(11)=DCMPLX(FR(11),FI(11))
       P(154) = -s45+2*P(39)
       P(155) = p3sq+2*P(16)+4*P(39)
       P(156) = s12-s34-s45
       P(157) = p3sq+2*p5sq-5*s34+4*P(16)
       P(158) = 3*p3sq+2*p5sq-4*s45+5*P(6)
       P(159) = p2sq-3*s12-2*P(107)
       P(160) = p2sq-2*p3sq-s12-s23+s45
       P(161) = p3sq-s12
       P(162) = p5sq+s15-s23+s45+2*P(161)
       P(163) = p2sq-p5sq-s12+s15-s23+s45
       P(164) = 3*p2sq-s12-2*P(112)
       FR(12) = 2*(D02345R-Dij2345R(1,2)-Dij2345R(3,1)+Dij2345R(5,2)+4
     -   *(EijR(11,2)-EijR(22,3)+EijR(37,4)-EijR(44,4))-EijR(2,3)*P(5
     -   )+EijR(10,2)*P(15)+EijR(15,3)*P(21)-EijR(4,2)*P(23)+(-EijR(9
     -   ,3)+EijR(20,3))*P(46)+EijR(10,3)*P(68)+Is12*(C0135R+Cij135R(
     -   1,1)-3*Cij145R(2,1)+Cij345R(2,1)+8*(Dij1345R(7,2)+Dij1345R(1
     -   1,3))-4*Dij1345R(13,3)+Dij1345R(1,2)*P(55)+Dij1345R(6,2)*P(5
     -   6)+Is45*(2*B012R+4*Cij145R(4,2)+C0145R*P(95)+Cij145R(1,1)*P(
     -   96))+Dij1345R(3,2)*P(101)+D01345R*P(154)+Dij1345R(1,1)*P(155
     -   )-2*(Cij135R(2,1)-Cij135R(2,2)+Cij135R(3,2)-Cij145R(1,2)+Cij
     -   145R(3,2)-s34*(-Dij1345R(4,3)+Dij1345R(10,3))-Dij1345R(4,2)*
     -   P(119)-Dij1345R(2,1)*P(156))-Dij1345R(3,1)*P(157)-Dij1345R(5
     -   ,2)*P(158))-EijR(2,1)*P(159)+EijR(3,1)*P(160)+EijR(4,1)*P(16
     -   2)-2*(Dij1235R(3,2)-Dij1235R(6,2)-Dij1345R(1,3)+Dij1345R(5,3
     -   )-Dij2345R(1,1)-p2sq*EijR(2,2)-(EE0R+EijR(1,1))*P(3)-EijR(7,
     -   2)*P(10)-(EijR(8,3)-EijR(18,3))*P(33)+EijR(8,2)*P(140)+EijR(
     -   5,2)*P(163))-EijR(9,2)*P(164))
       FI(12) = 2*(D02345I-Dij2345I(1,2)-Dij2345I(3,1)+Dij2345I(5,2)+4
     -   *(EijI(11,2)-EijI(22,3)+EijI(37,4)-EijI(44,4))-EijI(2,3)*P(5
     -   )+EijI(10,2)*P(15)+EijI(15,3)*P(21)-EijI(4,2)*P(23)+(-EijI(9
     -   ,3)+EijI(20,3))*P(46)+EijI(10,3)*P(68)+Is12*(C0135I+Cij135I(
     -   1,1)-3*Cij145I(2,1)+Cij345I(2,1)+8*(Dij1345I(7,2)+Dij1345I(1
     -   1,3))-4*Dij1345I(13,3)+Dij1345I(1,2)*P(55)+Dij1345I(6,2)*P(5
     -   6)+Is45*(2*B012I+4*Cij145I(4,2)+C0145I*P(95)+Cij145I(1,1)*P(
     -   96))+Dij1345I(3,2)*P(101)+D01345I*P(154)+Dij1345I(1,1)*P(155
     -   )-2*(Cij135I(2,1)-Cij135I(2,2)+Cij135I(3,2)-Cij145I(1,2)+Cij
     -   145I(3,2)-s34*(-Dij1345I(4,3)+Dij1345I(10,3))-Dij1345I(4,2)*
     -   P(119)-Dij1345I(2,1)*P(156))-Dij1345I(3,1)*P(157)-Dij1345I(5
     -   ,2)*P(158))-EijI(2,1)*P(159)+EijI(3,1)*P(160)+EijI(4,1)*P(16
     -   2)-2*(Dij1235I(3,2)-Dij1235I(6,2)-Dij1345I(1,3)+Dij1345I(5,3
     -   )-Dij2345I(1,1)-p2sq*EijI(2,2)-(EE0I+EijI(1,1))*P(3)-EijI(7,
     -   2)*P(10)-(EijI(8,3)-EijI(18,3))*P(33)+EijI(8,2)*P(140)+EijI(
     -   5,2)*P(163))-EijI(9,2)*P(164))
       F(12)=DCMPLX(FR(12),FI(12))
       P(165) = p3sq+s12-2*s34+s45
       P(166) = 3*p2sq-s12
       P(167) = p3sq-s15-s34
       P(168) = 4*p2sq-s23-s45+2*P(167)
       FR(13) = 2*(Dij2345R(1,1)-Dij2345R(1,2)-Dij2345R(2,1)+Dij2345R(
     -   4,2)-4*(-3*EijR(11,2)+EijR(22,3)-EijR(37,4)+EijR(43,4))-EijR
     -   (2,3)*P(5)+EijR(3,2)*P(15)-(-2*EijR(2,1)+EijR(10,3)-EijR(20,
     -   3))*P(21)-(-EijR(9,2)+EijR(10,2))*P(23)+EijR(12,3)*P(46)+Eij
     -   R(9,3)*P(65)-2*(D02345R+Dij1235R(3,2)-Dij1235R(6,2)-Dij1345R
     -   (1,3)+EijR(4,1)*P(7)-EijR(6,2)*P(10)+(-EijR(8,3)+EijR(17,3))
     -   *P(33)-EijR(3,1)*P(66)+EijR(5,2)*P(163))+Is12*(C0135R-C0345R
     -   -D01345R*s45+Cij135R(1,1)+Cij345R(1,1)+12*Dij1345R(7,2)+8*Di
     -   j1345R(11,3)-4*Dij1345R(12,3)-3*Dij1345R(3,1)*P(7)-2*(Cij135
     -   R(2,1)-Cij135R(2,2)+Cij135R(3,2)-s34*Dij1345R(6,3)+Dij1345R(
     -   4,3)*P(29))+Dij1345R(1,1)*P(32)-Dij1345R(5,2)*P(41)+Dij1345R
     -   (1,2)*P(55)+Dij1345R(2,2)*P(56)-Dij1345R(4,2)*P(57)+Dij1345R
     -   (6,2)*P(101)-Dij1345R(2,1)*P(165))+EijR(2,2)*P(166)-EijR(8,2
     -   )*P(168))
       FI(13) = 2*(Dij2345I(1,1)-Dij2345I(1,2)-Dij2345I(2,1)+Dij2345I(
     -   4,2)-4*(-3*EijI(11,2)+EijI(22,3)-EijI(37,4)+EijI(43,4))-EijI
     -   (2,3)*P(5)+EijI(3,2)*P(15)-(-2*EijI(2,1)+EijI(10,3)-EijI(20,
     -   3))*P(21)-(-EijI(9,2)+EijI(10,2))*P(23)+EijI(12,3)*P(46)+Eij
     -   I(9,3)*P(65)-2*(D02345I+Dij1235I(3,2)-Dij1235I(6,2)-Dij1345I
     -   (1,3)+EijI(4,1)*P(7)-EijI(6,2)*P(10)+(-EijI(8,3)+EijI(17,3))
     -   *P(33)-EijI(3,1)*P(66)+EijI(5,2)*P(163))+Is12*(C0135I-C0345I
     -   -D01345I*s45+Cij135I(1,1)+Cij345I(1,1)+12*Dij1345I(7,2)+8*Di
     -   j1345I(11,3)-4*Dij1345I(12,3)-3*Dij1345I(3,1)*P(7)-2*(Cij135
     -   I(2,1)-Cij135I(2,2)+Cij135I(3,2)-s34*Dij1345I(6,3)+Dij1345I(
     -   4,3)*P(29))+Dij1345I(1,1)*P(32)-Dij1345I(5,2)*P(41)+Dij1345I
     -   (1,2)*P(55)+Dij1345I(2,2)*P(56)-Dij1345I(4,2)*P(57)+Dij1345I
     -   (6,2)*P(101)-Dij1345I(2,1)*P(165))+EijI(2,2)*P(166)-EijI(8,2
     -   )*P(168))
       F(13)=DCMPLX(FR(13),FI(13))
       P(169) = s15-s23+s34
       P(170) = p3sq-P(169)
       P(171) = p5sq-s45+2*P(170)
       P(172) = p3sq-s15+s23-s34
       FR(14) = -2*(D01235R-EE0R*s45+Dij1235R(1,1)-EijR(4,1)*P(7)+EijR
     -   (1,2)*P(50)+EijR(1,1)*P(77)-EijR(2,1)*P(78)+EijR(3,1)*P(79)+
     -   EijR(5,2)*P(82)-EijR(6,2)*P(85)-EijR(7,2)*P(171)+2*(D02345R-
     -   Dij1235R(3,1)+Dij1235R(5,2)-Dij1235R(6,2)+Dij2345R(1,2)+s12*
     -   EijR(2,3)-p2sq*EijR(8,3)-EijR(11,2)+2*(Dij2345R(1,1)-EijR(22
     -   ,3)-EijR(37,4)+EijR(40,4))-EijR(18,3)*P(21)+EijR(5,3)*P(33)-
     -   EijR(2,2)*P(72)+EijR(8,2)*P(87)+(EijR(9,3)-EijR(17,3))*P(91)
     -   -EijR(10,3)*P(115)+EijR(9,2)*P(172)))
       FI(14) = -2*(D01235I-EE0I*s45+Dij1235I(1,1)-EijI(4,1)*P(7)+EijI
     -   (1,2)*P(50)+EijI(1,1)*P(77)-EijI(2,1)*P(78)+EijI(3,1)*P(79)+
     -   EijI(5,2)*P(82)-EijI(6,2)*P(85)-EijI(7,2)*P(171)+2*(D02345I-
     -   Dij1235I(3,1)+Dij1235I(5,2)-Dij1235I(6,2)+Dij2345I(1,2)+s12*
     -   EijI(2,3)-p2sq*EijI(8,3)-EijI(11,2)+2*(Dij2345I(1,1)-EijI(22
     -   ,3)-EijI(37,4)+EijI(40,4))-EijI(18,3)*P(21)+EijI(5,3)*P(33)-
     -   EijI(2,2)*P(72)+EijI(8,2)*P(87)+(EijI(9,3)-EijI(17,3))*P(91)
     -   -EijI(10,3)*P(115)+EijI(9,2)*P(172)))
       F(14)=DCMPLX(FR(14),FI(14))
       P(173) = -2*p5sq+s34
       P(174) = s12+s34-s45
       P(175) = p3sq+3*s12-s45+2*P(39)
       P(176) = -3*s34+2*P(16)
       P(177) = 4*p5sq+s12-s34-s45
       P(178) = -s12+2*P(107)
       P(179) = 2*p2sq+3*s12
       P(180) = s15-s45
       P(181) = -s23+3*P(6)+2*P(180)
       P(182) = p2sq-s12-s23-5*s34+s45
       P(183) = p2sq-s23
       P(184) = -p5sq+s15+s45+2*P(6)+3*P(183)
       FR(15) = Is34*(2*(-Cij125R(2,1)-Cij235R(1,1)+p2sq*Dij1235R(2,1)
     -   -6*Dij1235R(7,2)+4*Dij1235R(12,3)-Dij1235R(2,2)*P(5)+2*(C023
     -   5R-Cij125R(2,2)+Cij135R(1,2)+Dij1235R(4,2)*P(33))+(D01235R+D
     -   ij1235R(1,1))*P(133)+Dij1235R(3,1)*P(138)+Dij1235R(6,2)*P(14
     -   2))+Is12*(4*B035R+8*Cij135R(4,2)+2*(Cij135R(1,1)*P(130)-C013
     -   5R*P(173))))+Is12*(4*Is45*(B012R+C0145R*p5sq+2*Cij145R(4,2)+
     -   Cij145R(1,1)*P(129))+2*(3*C0345R-Cij345R(2,1)-p5sq*Dij1345R(
     -   3,1)-2*(Cij135R(2,1)+Cij135R(3,2)-Cij145R(1,2)+s34*Dij1345R(
     -   4,3)-Dij1345R(7,2)-s34*Dij1345R(10,3))+8*Dij1345R(11,3)-4*Di
     -   j1345R(13,3)+Dij1345R(3,2)*P(7)+D01345R*P(118)+Dij1345R(1,1)
     -   *P(120)-Dij1345R(6,2)*P(140)-Dij1345R(2,1)*P(174)+Dij1345R(1
     -   ,2)*P(175)+Dij1345R(4,2)*P(176)-Dij1345R(5,2)*P(177)))+2*(-D
     -   ij2345R(1,1)-Dij2345R(1,2)+Dij2345R(5,2)-3*s34*EijR(10,2)-12
     -   *EijR(11,2)-6*(EijR(22,3)-EijR(24,3))+4*(EijR(37,4)-EijR(44,
     -   4))-EijR(2,3)*P(5)+EijR(15,3)*P(21)+2*(D02345R+Dij1345R(1,3)
     -   -Dij1345R(5,3)+(EijR(8,3)-EijR(18,3))*P(33))+(-EijR(9,3)+Eij
     -   R(20,3))*P(46)+EijR(10,3)*P(68)+(EE0R+EijR(1,1))*P(127)+EijR
     -   (5,2)*P(145)-EijR(7,2)*P(146)-EijR(3,1)*P(176)+EijR(2,1)*P(1
     -   78)-EijR(2,2)*P(179)+EijR(4,1)*P(181)-EijR(8,2)*P(182)+EijR(
     -   9,2)*P(184))
       FI(15) = Is34*(2*(-Cij125I(2,1)-Cij235I(1,1)+p2sq*Dij1235I(2,1)
     -   -6*Dij1235I(7,2)+4*Dij1235I(12,3)-Dij1235I(2,2)*P(5)+2*(C023
     -   5I-Cij125I(2,2)+Cij135I(1,2)+Dij1235I(4,2)*P(33))+(D01235I+D
     -   ij1235I(1,1))*P(133)+Dij1235I(3,1)*P(138)+Dij1235I(6,2)*P(14
     -   2))+Is12*(4*B035I+8*Cij135I(4,2)+2*(Cij135I(1,1)*P(130)-C013
     -   5I*P(173))))+Is12*(4*Is45*(B012I+C0145I*p5sq+2*Cij145I(4,2)+
     -   Cij145I(1,1)*P(129))+2*(3*C0345I-Cij345I(2,1)-p5sq*Dij1345I(
     -   3,1)-2*(Cij135I(2,1)+Cij135I(3,2)-Cij145I(1,2)+s34*Dij1345I(
     -   4,3)-Dij1345I(7,2)-s34*Dij1345I(10,3))+8*Dij1345I(11,3)-4*Di
     -   j1345I(13,3)+Dij1345I(3,2)*P(7)+D01345I*P(118)+Dij1345I(1,1)
     -   *P(120)-Dij1345I(6,2)*P(140)-Dij1345I(2,1)*P(174)+Dij1345I(1
     -   ,2)*P(175)+Dij1345I(4,2)*P(176)-Dij1345I(5,2)*P(177)))+2*(-D
     -   ij2345I(1,1)-Dij2345I(1,2)+Dij2345I(5,2)-3*s34*EijI(10,2)-12
     -   *EijI(11,2)-6*(EijI(22,3)-EijI(24,3))+4*(EijI(37,4)-EijI(44,
     -   4))-EijI(2,3)*P(5)+EijI(15,3)*P(21)+2*(D02345I+Dij1345I(1,3)
     -   -Dij1345I(5,3)+(EijI(8,3)-EijI(18,3))*P(33))+(-EijI(9,3)+Eij
     -   I(20,3))*P(46)+EijI(10,3)*P(68)+(EE0I+EijI(1,1))*P(127)+EijI
     -   (5,2)*P(145)-EijI(7,2)*P(146)-EijI(3,1)*P(176)+EijI(2,1)*P(1
     -   78)-EijI(2,2)*P(179)+EijI(4,1)*P(181)-EijI(8,2)*P(182)+EijI(
     -   9,2)*P(184))
       F(15)=DCMPLX(FR(15),FI(15))
       P(185) = 2*p5sq-s45
       P(186) = -2+3*Is12*s34
       P(187) = -(p5sq*s34)+P(6)**2
       P(188) = -p2sq+s34
       P(189) = -p5sq+s34
       P(190) = 5*s12**2-7*s12*s34+2*s34*P(189)
       P(191) = p2sq+s12+2*s15-s34
       P(192) = p3sq+2*s12
       P(193) = -2*s12*s45+3*s34*s45+P(6)*P(192)
       P(194) = p5sq*s12-s12*s45
       P(195) = p5sq*s34-s34*s45
       P(196) = -(s12*s34)+s34**2+p3sq*P(6)+2*P(194)-3*P(195)
       P(197) = -p5sq+s45
       P(198) = p3sq+s12-s45
       P(199) = p3sq+s12+2*s34-s45
       P(200) = s12-2*s34
       P(201) = -s12+s34+2*P(3)
       P(202) = s12*P(26)+p2sq*P(201)
       P(203) = -s12-2*s15+s34
       P(204) = 2*s15+s34
       P(205) = p2sq**2+s12*P(203)-p2sq*P(204)
       P(206) = p3sq*p5sq+p5sq*s12-s12*s15-p5sq*s45+s15*s45
       P(207) = s12**2+s12*s23-s12*s34-s23*s34-s12*s45+s34*s45+p2sq*P(
     -   148)+2*P(206)
       P(208) = s15+s45
       P(209) = s12-3*s34-2*P(208)
       P(210) = p3sq*p5sq+p5sq**2+s15*s45
       P(211) = -(s12*s15)+s12*s23+s15*s34-s23*s34-s12*s45+s34*s45+p5s
     -   q*P(209)+2*P(210)
       P(212) = 2*p5sq+s12
       P(213) = p5sq-s12-2*s15
       P(214) = p2sq**2+s12*P(212)+2*p2sq*P(213)
       P(215) = s12+s23+2*P(3)
       P(216) = s12+2*P(3)
       P(217) = -2*s12-s23+s45
       P(218) = p2sq**2+2*p3sq*P(3)+s12*P(215)-s45*P(216)+p2sq*P(217)
       P(219) = p3sq*p5sq+p5sq**2-p3sq*s15-p5sq*s15-p5sq*s34+s15*s34
       P(220) = -(p5sq*s12)+s12*s15-s12*s23+p2sq*P(21)+s45*P(216)-2*P(
     -   219)
       P(221) = s12-s23-s34
       P(222) = -s23+s34-2*P(180)
       P(223) = s12**2+p3sq*s15-s15*s45
       P(224) = p2sq**2+s23*s34-s34*s45+p2sq*P(221)+s12*P(222)-2*P(223
     -   )
       P(225) = s12+s15
       P(226) = s34-2*P(225)
       P(227) = s15-s23+2*s45
       P(228) = p3sq*s15-s15*s45
       P(229) = s15*s34+s23*s34-s34*s45+p2sq*P(66)+p5sq*P(226)+s12*P(2
     -   27)-2*P(228)
       P(230) = p5sq+s12+s23
       P(231) = p2sq+s15-2*P(230)
       P(232) = p5sq*s12+s45**2
       P(233) = p2sq*p5sq+s12*s15-p5sq*s23-s12*s23-s45*P(231)-2*P(232)
       P(234) = p5sq+s15+s34
       P(235) = 3*p2sq-s12-2*P(234)
       FR(16) = 4*Dij1345R(11,3)+8*EijR(46,4)+(Dij2345R(2,1)+EijR(3,2)
     -   *P(10))*P(16)-EijR(4,2)*P(7)*P(21)-Dij1345R(5,2)*P(101)+C034
     -   5R*P(186)+D02345R*P(188)-Dij2345R(1,1)*P(191)+Dij2345R(3,1)*
     -   P(197)+Dij1345R(4,2)*P(199)+Is12*(s34*Dij1345R(6,2)*P(101)-I
     -   s45*P(6)*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)*P(95)+C0145R*P
     -   (185))-2*D01345R*P(187)-Dij1345R(1,1)*P(190)+Dij1345R(2,1)*P
     -   (193)-Dij1345R(3,1)*P(196)+s34*(Cij345R(1,1)-4*Dij1345R(12,3
     -   )-Dij1345R(2,2)*P(198))+4*Dij1345R(7,2)*P(200))-EijR(2,1)*P(
     -   202)+EijR(2,2)*P(205)+EijR(3,1)*P(207)-EijR(4,1)*P(211)-EijR
     -   (5,2)*P(214)+EijR(6,2)*P(218)+EijR(7,2)*P(220)-EijR(8,2)*P(2
     -   24)+EijR(9,2)*P(229)-EijR(10,2)*P(233)-2*(s12*Dij1345R(1,2)-
     -   Dij2345R(7,2)-EijR(22,3)*P(5)+EijR(23,3)*P(15)-EijR(24,3)*P(
     -   18)-P(3)*(-(EE0R*P(6))+EijR(1,2)*P(33)+EijR(1,1)*P(71))-EijR
     -   (11,2)*P(235))
       FI(16) = 4*Dij1345I(11,3)+8*EijI(46,4)+(Dij2345I(2,1)+EijI(3,2)
     -   *P(10))*P(16)-EijI(4,2)*P(7)*P(21)-Dij1345I(5,2)*P(101)+C034
     -   5I*P(186)+D02345I*P(188)-Dij2345I(1,1)*P(191)+Dij2345I(3,1)*
     -   P(197)+Dij1345I(4,2)*P(199)+Is12*(s34*Dij1345I(6,2)*P(101)-I
     -   s45*P(6)*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)*P(95)+C0145I*P
     -   (185))-2*D01345I*P(187)-Dij1345I(1,1)*P(190)+Dij1345I(2,1)*P
     -   (193)-Dij1345I(3,1)*P(196)+s34*(Cij345I(1,1)-4*Dij1345I(12,3
     -   )-Dij1345I(2,2)*P(198))+4*Dij1345I(7,2)*P(200))-EijI(2,1)*P(
     -   202)+EijI(2,2)*P(205)+EijI(3,1)*P(207)-EijI(4,1)*P(211)-EijI
     -   (5,2)*P(214)+EijI(6,2)*P(218)+EijI(7,2)*P(220)-EijI(8,2)*P(2
     -   24)+EijI(9,2)*P(229)-EijI(10,2)*P(233)-2*(s12*Dij1345I(1,2)-
     -   Dij2345I(7,2)-EijI(22,3)*P(5)+EijI(23,3)*P(15)-EijI(24,3)*P(
     -   18)-P(3)*(-(EE0I*P(6))+EijI(1,2)*P(33)+EijI(1,1)*P(71))-EijI
     -   (11,2)*P(235))
       F(16)=DCMPLX(FR(16),FI(16))
       P(236) = p3sq+s12-3*s34+2*P(7)
       P(237) = 2*p3sq+4*p5sq-5*s34+3*P(16)
       P(238) = s12+2*P(107)
       P(239) = -3*s34+2*P(198)
       P(240) = p2sq-s12+2*P(107)
       P(241) = p2sq+2*p3sq+s12-s23-5*s34-s45
       FR(17) = Is12*(4*Is45*(B012R+C0145R*p5sq+2*Cij145R(4,2)+Cij145R
     -   (1,1)*P(129))+2*(C0135R+C0345R+Cij135R(1,1)-3*Cij135R(2,1)-2
     -   *(Cij135R(3,2)-Cij145R(1,2)+s34*Dij1345R(4,3))+6*Dij1345R(7,
     -   2)+8*Dij1345R(11,3)-Dij1345R(3,1)*P(7)-Dij1345R(5,2)*P(41)+D
     -   ij1345R(2,1)*P(156)+Dij1345R(1,2)*P(175)+Dij1345R(4,2)*P(176
     -   )+D01345R*P(236)+Dij1345R(1,1)*P(237)))+2*(Dij1235R(3,1)+Dij
     -   2345R(1,1)-Dij2345R(1,2)-s12*EijR(2,2)-6*EijR(22,3)-4*(EijR(
     -   11,2)-EijR(37,4))-EijR(2,3)*P(5)+EijR(9,2)*P(23)+2*(D02345R+
     -   Dij1235R(6,2)+Dij1345R(1,3)+EijR(8,3)*P(33))+EijR(4,1)*P(41)
     -   -EijR(9,3)*P(46)+EijR(10,3)*P(115)+(EE0R+EijR(1,1))*P(127)+E
     -   ijR(2,1)*P(238)-EijR(3,1)*P(239)+EijR(5,2)*P(240)-EijR(8,2)*
     -   P(241))
       FI(17) = Is12*(4*Is45*(B012I+C0145I*p5sq+2*Cij145I(4,2)+Cij145I
     -   (1,1)*P(129))+2*(C0135I+C0345I+Cij135I(1,1)-3*Cij135I(2,1)-2
     -   *(Cij135I(3,2)-Cij145I(1,2)+s34*Dij1345I(4,3))+6*Dij1345I(7,
     -   2)+8*Dij1345I(11,3)-Dij1345I(3,1)*P(7)-Dij1345I(5,2)*P(41)+D
     -   ij1345I(2,1)*P(156)+Dij1345I(1,2)*P(175)+Dij1345I(4,2)*P(176
     -   )+D01345I*P(236)+Dij1345I(1,1)*P(237)))+2*(Dij1235I(3,1)+Dij
     -   2345I(1,1)-Dij2345I(1,2)-s12*EijI(2,2)-6*EijI(22,3)-4*(EijI(
     -   11,2)-EijI(37,4))-EijI(2,3)*P(5)+EijI(9,2)*P(23)+2*(D02345I+
     -   Dij1235I(6,2)+Dij1345I(1,3)+EijI(8,3)*P(33))+EijI(4,1)*P(41)
     -   -EijI(9,3)*P(46)+EijI(10,3)*P(115)+(EE0I+EijI(1,1))*P(127)+E
     -   ijI(2,1)*P(238)-EijI(3,1)*P(239)+EijI(5,2)*P(240)-EijI(8,2)*
     -   P(241))
       F(17)=DCMPLX(FR(17),FI(17))
       P(242) = 2*p3sq-3*s45+4*P(6)
       P(243) = p2sq-s34
       P(244) = p3sq+s12-s23+2*P(243)
       P(245) = -s23+s45+2*P(33)
       P(246) = -3*s45+4*P(6)
       P(247) = 4*p2sq+p3sq-s12+s45-2*P(83)
       P(248) = p3sq-p5sq+s12-s23+s45+2*P(243)
       P(249) = s23+s34-s45
       P(250) = 4*p2sq+p3sq-p5sq-s12-2*P(249)
       FR(18) = 2*(Dij1235R(2,2)-Dij1235R(4,2)+Dij1235R(5,2)-Dij1235R(
     -   6,2)+Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(3,1)-Dij2345R(5,2)
     -   +2*(Dij1235R(2,3)-Dij1235R(6,3)-Dij1235R(7,3)+Dij1235R(9,3)+
     -   p3sq*EijR(6,2)-p3sq*EijR(8,2))-6*(EijR(21,3)-EijR(22,3))-4*(
     -   Dij1235R(8,3)-Dij1235R(10,3)-EijR(37,4)+EijR(40,4)-EijR(42,4
     -   )+EijR(44,4))+EijR(14,3)*P(7)+EijR(15,3)*P(197)+(-EijR(9,3)+
     -   EijR(17,3)-EijR(19,3)+EijR(20,3))*P(198)+(EijR(2,2)-EijR(5,2
     -   ))*P(242)+EijR(2,3)*P(244)+(EijR(5,3)-EijR(7,3))*P(245)+(Eij
     -   R(7,2)-EijR(9,2))*P(246)-EijR(8,3)*P(247)-EijR(10,3)*P(248)+
     -   EijR(18,3)*P(250))
       FI(18) = 2*(Dij1235I(2,2)-Dij1235I(4,2)+Dij1235I(5,2)-Dij1235I(
     -   6,2)+Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(3,1)-Dij2345I(5,2)
     -   +2*(Dij1235I(2,3)-Dij1235I(6,3)-Dij1235I(7,3)+Dij1235I(9,3)+
     -   p3sq*EijI(6,2)-p3sq*EijI(8,2))-6*(EijI(21,3)-EijI(22,3))-4*(
     -   Dij1235I(8,3)-Dij1235I(10,3)-EijI(37,4)+EijI(40,4)-EijI(42,4
     -   )+EijI(44,4))+EijI(14,3)*P(7)+EijI(15,3)*P(197)+(-EijI(9,3)+
     -   EijI(17,3)-EijI(19,3)+EijI(20,3))*P(198)+(EijI(2,2)-EijI(5,2
     -   ))*P(242)+EijI(2,3)*P(244)+(EijI(5,3)-EijI(7,3))*P(245)+(Eij
     -   I(7,2)-EijI(9,2))*P(246)-EijI(8,3)*P(247)-EijI(10,3)*P(248)+
     -   EijI(18,3)*P(250))
       F(18)=DCMPLX(FR(18),FI(18))
       P(251) = -p2sq+s12+s23-s45
       P(252) = p3sq+2*P(156)
       FR(19) = 4*(Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(3,2)-Dij1235R(
     -   3,3)-3*(Dij1235R(8,3)-Dij1235R(9,3))-p3sq*(EijR(8,2)+EijR(9,
     -   3)-EijR(10,2)+EijR(16,3)-2*EijR(20,3))+(EijR(2,3)-2*EijR(10,
     -   3)+EijR(15,3))*P(1)-2*(Dij1235R(6,2)-EijR(22,3)+EijR(24,3)-E
     -   ijR(37,4)-EijR(39,4)+2*EijR(44,4)-EijR(18,3)*P(10))+EijR(2,2
     -   )*P(140)+EijR(4,2)*P(156)+(EijR(8,3)+EijR(14,3))*P(251)-EijR
     -   (9,2)*P(252))
       FI(19) = 4*(Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(3,2)-Dij1235I(
     -   3,3)-3*(Dij1235I(8,3)-Dij1235I(9,3))-p3sq*(EijI(8,2)+EijI(9,
     -   3)-EijI(10,2)+EijI(16,3)-2*EijI(20,3))+(EijI(2,3)-2*EijI(10,
     -   3)+EijI(15,3))*P(1)-2*(Dij1235I(6,2)-EijI(22,3)+EijI(24,3)-E
     -   ijI(37,4)-EijI(39,4)+2*EijI(44,4)-EijI(18,3)*P(10))+EijI(2,2
     -   )*P(140)+EijI(4,2)*P(156)+(EijI(8,3)+EijI(14,3))*P(251)-EijI
     -   (9,2)*P(252))
       F(19)=DCMPLX(FR(19),FI(19))
       P(253) = p3sq+s12-2*s34
       P(254) = -3*s45+2*P(253)
       P(255) = -s45+2*P(253)
       P(256) = 2*p2sq+p3sq-s23-s34
       FR(20) = -2*(-Dij1235R(2,2)+Dij1235R(4,2)-Dij1235R(5,2)+Dij1235
     -   R(6,2)+Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(2,1)+Dij2345R(4,
     -   2)+4*(Dij1235R(8,3)-Dij1235R(10,3)+2*(EijR(21,3)-EijR(22,3))
     -   -EijR(37,4)+EijR(40,4)-EijR(41,4)+EijR(43,4))-(2*(EijR(7,2)-
     -   EijR(9,2))+EijR(10,3)-EijR(18,3)+EijR(19,3))*P(7)+EijR(9,3)*
     -   P(25)-EijR(20,3)*P(197)+(EijR(11,3)-EijR(12,3))*P(198)-EijR(
     -   2,3)*P(244)+(-EijR(5,3)+EijR(6,3))*P(245)+EijR(8,3)*P(247)+(
     -   -EijR(2,2)+EijR(5,2))*P(254)+(-EijR(6,2)+EijR(8,2))*P(255)+2
     -   *(D02345R-Dij1235R(2,3)+Dij1235R(6,3)+Dij1235R(7,3)-Dij1235R
     -   (9,3)-EijR(17,3)*P(256)))
       FI(20) = -2*(-Dij1235I(2,2)+Dij1235I(4,2)-Dij1235I(5,2)+Dij1235
     -   I(6,2)+Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(2,1)+Dij2345I(4,
     -   2)+4*(Dij1235I(8,3)-Dij1235I(10,3)+2*(EijI(21,3)-EijI(22,3))
     -   -EijI(37,4)+EijI(40,4)-EijI(41,4)+EijI(43,4))-(2*(EijI(7,2)-
     -   EijI(9,2))+EijI(10,3)-EijI(18,3)+EijI(19,3))*P(7)+EijI(9,3)*
     -   P(25)-EijI(20,3)*P(197)+(EijI(11,3)-EijI(12,3))*P(198)-EijI(
     -   2,3)*P(244)+(-EijI(5,3)+EijI(6,3))*P(245)+EijI(8,3)*P(247)+(
     -   -EijI(2,2)+EijI(5,2))*P(254)+(-EijI(6,2)+EijI(8,2))*P(255)+2
     -   *(D02345I-Dij1235I(2,3)+Dij1235I(6,3)+Dij1235I(7,3)-Dij1235I
     -   (9,3)-EijI(17,3)*P(256)))
       F(20)=DCMPLX(FR(20),FI(20))
       P(257) = s12+s23-s45
       P(258) = p3sq-s12-s23-s34
       P(259) = s12-s23-s34+s45
       P(260) = -s12+s23+s34-s45
       P(261) = p3sq-s12-2*s23-s34+s45
       P(262) = p5sq-s12-s23
       P(263) = p5sq-s12-s23-s45+2*P(14)
       P(264) = p3sq-P(83)
       P(265) = p5sq-s12+2*P(264)
       FR(21) = -4*(Dij1235R(2,1)-Dij1235R(2,3)-Dij1235R(3,1)+Dij1235R
     -   (3,2)+Dij1235R(4,2)-Dij1235R(5,2)-Dij1235R(6,2)+Dij1235R(6,3
     -   )+Dij1235R(7,3)-Dij1235R(9,3)+Dij2345R(1,2)-Dij2345R(5,2)+2*
     -   (Dij1235R(8,3)-Dij1235R(10,3)-EijR(11,2)-EijR(21,3)+2*EijR(2
     -   2,3)-EijR(24,3)-EijR(37,4)+EijR(40,4)-EijR(42,4)+EijR(44,4))
     -   +(-EijR(2,1)+EijR(4,1))*P(6)+EijR(4,2)*P(7)+(EijR(8,2)-EijR(
     -   9,3)-EijR(10,2)+EijR(17,3)-EijR(19,3)+EijR(20,3))*P(16)+(Eij
     -   R(5,3)-EijR(7,3))*P(50)+(EijR(14,3)-EijR(15,3))*P(101)-EijR(
     -   2,2)*P(257)-EijR(2,3)*P(258)-EijR(5,2)*P(259)-EijR(7,2)*P(26
     -   0)+EijR(8,3)*P(261)-EijR(9,2)*P(262)+EijR(10,3)*P(263)-EijR(
     -   18,3)*P(265))
       FI(21) = -4*(Dij1235I(2,1)-Dij1235I(2,3)-Dij1235I(3,1)+Dij1235I
     -   (3,2)+Dij1235I(4,2)-Dij1235I(5,2)-Dij1235I(6,2)+Dij1235I(6,3
     -   )+Dij1235I(7,3)-Dij1235I(9,3)+Dij2345I(1,2)-Dij2345I(5,2)+2*
     -   (Dij1235I(8,3)-Dij1235I(10,3)-EijI(11,2)-EijI(21,3)+2*EijI(2
     -   2,3)-EijI(24,3)-EijI(37,4)+EijI(40,4)-EijI(42,4)+EijI(44,4))
     -   +(-EijI(2,1)+EijI(4,1))*P(6)+EijI(4,2)*P(7)+(EijI(8,2)-EijI(
     -   9,3)-EijI(10,2)+EijI(17,3)-EijI(19,3)+EijI(20,3))*P(16)+(Eij
     -   I(5,3)-EijI(7,3))*P(50)+(EijI(14,3)-EijI(15,3))*P(101)-EijI(
     -   2,2)*P(257)-EijI(2,3)*P(258)-EijI(5,2)*P(259)-EijI(7,2)*P(26
     -   0)+EijI(8,3)*P(261)-EijI(9,2)*P(262)+EijI(10,3)*P(263)-EijI(
     -   18,3)*P(265))
       F(21)=DCMPLX(FR(21),FI(21))
       P(266) = p2sq+s12-2*P(61)
       P(267) = -2*s15+s23+s34-s45
       P(268) = p5sq+P(267)
       P(269) = p2sq-3*s12+2*P(268)
       P(270) = p2sq-3*s12+s23+2*s34-s45
       P(271) = p2sq+s23+3*P(16)
       P(272) = p3sq-2*s12-s23-s34+s45
       P(273) = p3sq-s34-2*P(257)
       FR(22) = -2*(D02345R+4*(Dij1235R(8,3)-Dij1235R(10,3))+Dij2345R(
     -   1,1)-EijR(7,2)*P(21)-EijR(9,2)*P(23)-EijR(2,2)*P(266)+EijR(5
     -   ,2)*P(269)-EijR(6,2)*P(270)+EijR(8,2)*P(271)+2*(Dij1235R(2,1
     -   )-Dij1235R(2,3)-Dij1235R(3,1)+Dij1235R(3,2)+Dij1235R(4,2)-Di
     -   j1235R(5,2)-Dij1235R(6,2)+Dij1235R(6,3)+Dij1235R(7,3)-Dij123
     -   5R(9,3)+Dij2345R(1,2)-Dij2345R(4,2)-3*EijR(11,2)+2*(EijR(22,
     -   3)-EijR(23,3)-EijR(37,4)+EijR(40,4)-EijR(41,4)+EijR(43,4))-(
     -   EijR(1,1)+EijR(1,2))*P(3)+EijR(2,1)*P(4)+EijR(3,1)*P(6)+EijR
     -   (10,2)*P(7)-(EijR(3,2)+EijR(11,3)-EijR(12,3))*P(16)+(EijR(5,
     -   3)-EijR(6,3))*P(50)+(EijR(10,3)-EijR(18,3)+EijR(19,3)-EijR(2
     -   0,3))*P(101)-EijR(2,3)*P(258)+EijR(8,3)*P(261)+EijR(9,3)*P(2
     -   72)-EijR(17,3)*P(273)))
       FI(22) = -2*(D02345I+4*(Dij1235I(8,3)-Dij1235I(10,3))+Dij2345I(
     -   1,1)-EijI(7,2)*P(21)-EijI(9,2)*P(23)-EijI(2,2)*P(266)+EijI(5
     -   ,2)*P(269)-EijI(6,2)*P(270)+EijI(8,2)*P(271)+2*(Dij1235I(2,1
     -   )-Dij1235I(2,3)-Dij1235I(3,1)+Dij1235I(3,2)+Dij1235I(4,2)-Di
     -   j1235I(5,2)-Dij1235I(6,2)+Dij1235I(6,3)+Dij1235I(7,3)-Dij123
     -   5I(9,3)+Dij2345I(1,2)-Dij2345I(4,2)-3*EijI(11,2)+2*(EijI(22,
     -   3)-EijI(23,3)-EijI(37,4)+EijI(40,4)-EijI(41,4)+EijI(43,4))-(
     -   EijI(1,1)+EijI(1,2))*P(3)+EijI(2,1)*P(4)+EijI(3,1)*P(6)+EijI
     -   (10,2)*P(7)-(EijI(3,2)+EijI(11,3)-EijI(12,3))*P(16)+(EijI(5,
     -   3)-EijI(6,3))*P(50)+(EijI(10,3)-EijI(18,3)+EijI(19,3)-EijI(2
     -   0,3))*P(101)-EijI(2,3)*P(258)+EijI(8,3)*P(261)+EijI(9,3)*P(2
     -   72)-EijI(17,3)*P(273)))
       F(22)=DCMPLX(FR(22),FI(22))
       P(274) = p3sq-s34-s45
       P(275) = -p3sq+s34
       P(276) = p3sq+p5sq-s34-2*s45
       P(277) = p2sq+2*p3sq-s23-s34
       FR(23) = 4*(Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(3,2)-Dij1235R(
     -   3,3)-3*(Dij1235R(8,3)-Dij1235R(9,3))-Dij2345R(1,1)+Dij2345R(
     -   3,1)+p3sq*EijR(12,3)-p3sq*EijR(13,3)+4*(EijR(22,3)-EijR(24,3
     -   ))-2*(Dij1235R(6,2)-EijR(37,4)+EijR(43,4)+EijR(44,4)-EijR(45
     -   ,4))+(EijR(2,3)-EijR(10,3))*P(1)+EijR(4,2)*P(7)+(EijR(17,3)+
     -   EijR(18,3))*P(10)+EijR(10,2)*P(14)+(EijR(8,3)+EijR(19,3))*P(
     -   251)+EijR(2,2)*P(274)+EijR(8,2)*P(275)-EijR(9,2)*P(276)+(-Ei
     -   jR(9,3)+EijR(20,3))*P(277))
       FI(23) = 4*(Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(3,2)-Dij1235I(
     -   3,3)-3*(Dij1235I(8,3)-Dij1235I(9,3))-Dij2345I(1,1)+Dij2345I(
     -   3,1)+p3sq*EijI(12,3)-p3sq*EijI(13,3)+4*(EijI(22,3)-EijI(24,3
     -   ))-2*(Dij1235I(6,2)-EijI(37,4)+EijI(43,4)+EijI(44,4)-EijI(45
     -   ,4))+(EijI(2,3)-EijI(10,3))*P(1)+EijI(4,2)*P(7)+(EijI(17,3)+
     -   EijI(18,3))*P(10)+EijI(10,2)*P(14)+(EijI(8,3)+EijI(19,3))*P(
     -   251)+EijI(2,2)*P(274)+EijI(8,2)*P(275)-EijI(9,2)*P(276)+(-Ei
     -   jI(9,3)+EijI(20,3))*P(277))
       F(23)=DCMPLX(FR(23),FI(23))
       P(278) = p3sq-s15-s34-s45
       P(279) = -s15-s45+2*P(14)
       P(280) = 3*p3sq+2*P(84)
       P(281) = p2sq+3*p3sq-s23-s34
       FR(24) = 4*(Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(3,2)-Dij1235R(
     -   3,3)-3*(Dij1235R(8,3)-Dij1235R(9,3))-Dij2345R(1,1)+Dij2345R(
     -   2,1)-p5sq*EijR(2,1)+p5sq*EijR(3,1)-p3sq*EijR(3,3)+6*(EijR(22
     -   ,3)-EijR(23,3))-4*EijR(43,4)+EijR(2,3)*P(1)+(-EijR(5,2)+EijR
     -   (6,2))*P(3)+(-EijR(9,2)+EijR(10,2))*P(7)-2*(Dij1235R(6,2)-Ei
     -   jR(37,4)-EijR(38,4)-EijR(17,3)*P(10))+EijR(3,2)*P(14)+(EijR(
     -   8,3)+EijR(11,3))*P(251)+EijR(2,2)*P(278)-EijR(8,2)*P(279)-Ei
     -   jR(9,3)*P(280)+EijR(12,3)*P(281))
       FI(24) = 4*(Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(3,2)-Dij1235I(
     -   3,3)-3*(Dij1235I(8,3)-Dij1235I(9,3))-Dij2345I(1,1)+Dij2345I(
     -   2,1)-p5sq*EijI(2,1)+p5sq*EijI(3,1)-p3sq*EijI(3,3)+6*(EijI(22
     -   ,3)-EijI(23,3))-4*EijI(43,4)+EijI(2,3)*P(1)+(-EijI(5,2)+EijI
     -   (6,2))*P(3)+(-EijI(9,2)+EijI(10,2))*P(7)-2*(Dij1235I(6,2)-Ei
     -   jI(37,4)-EijI(38,4)-EijI(17,3)*P(10))+EijI(3,2)*P(14)+(EijI(
     -   8,3)+EijI(11,3))*P(251)+EijI(2,2)*P(278)-EijI(8,2)*P(279)-Ei
     -   jI(9,3)*P(280)+EijI(12,3)*P(281))
       F(24)=DCMPLX(FR(24),FI(24))
       P(282) = p3sq+s12-s15-s34-s45
       P(283) = 2*p3sq+s12-s34-s45
       P(284) = s12-s15-s34-s45
       P(285) = -s12+s34+s45
       FR(25) = 4*(Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(3,2)-Dij1235R(
     -   3,3)-3*(Dij1235R(8,3)-Dij1235R(9,3))+p5sq*(-EijR(2,1)+EijR(4
     -   ,1))+p3sq*(EijR(3,2)+EijR(12,3)-EijR(13,3))-2*(Dij1235R(6,2)
     -   -2*EijR(22,3)+EijR(23,3)+EijR(24,3)-EijR(37,4)+EijR(43,4)+Ei
     -   jR(44,4)-EijR(45,4))+(EijR(2,3)-EijR(10,3))*P(1)+(-EijR(5,2)
     -   +EijR(7,2))*P(3)+(EijR(17,3)+EijR(18,3))*P(10)+(EijR(8,3)+Ei
     -   jR(19,3))*P(251)+(-EijR(9,3)+EijR(20,3))*P(277)+EijR(2,2)*P(
     -   282)-EijR(8,2)*P(283)-EijR(9,2)*P(284)-EijR(10,2)*P(285))
       FI(25) = 4*(Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(3,2)-Dij1235I(
     -   3,3)-3*(Dij1235I(8,3)-Dij1235I(9,3))+p5sq*(-EijI(2,1)+EijI(4
     -   ,1))+p3sq*(EijI(3,2)+EijI(12,3)-EijI(13,3))-2*(Dij1235I(6,2)
     -   -2*EijI(22,3)+EijI(23,3)+EijI(24,3)-EijI(37,4)+EijI(43,4)+Ei
     -   jI(44,4)-EijI(45,4))+(EijI(2,3)-EijI(10,3))*P(1)+(-EijI(5,2)
     -   +EijI(7,2))*P(3)+(EijI(17,3)+EijI(18,3))*P(10)+(EijI(8,3)+Ei
     -   jI(19,3))*P(251)+(-EijI(9,3)+EijI(20,3))*P(277)+EijI(2,2)*P(
     -   282)-EijI(8,2)*P(283)-EijI(9,2)*P(284)-EijI(10,2)*P(285))
       F(25)=DCMPLX(FR(25),FI(25))
       P(286) = p3sq+s23-2*s34+3*P(16)
       P(287) = p3sq-s12-s23-2*s34
       P(288) = s23-s34
       P(289) = p3sq+3*s12-5*s45+2*P(288)
       P(290) = s34-s45
       P(291) = p3sq-s12-3*s23-2*P(290)
       P(292) = p5sq-s45+2*P(14)
       P(293) = -3*s23-4*s34+s45+2*P(161)
       P(294) = -p5sq+s45-2*P(14)
       FR(26) = 2*(Dij1235R(1,1)+Dij1235R(1,2)-Dij1235R(2,1)-Dij1235R(
     -   2,2)+2*(Dij1235R(2,3)+Dij1235R(4,3)-Dij1235R(5,2)-Dij1235R(5
     -   ,3)+Dij1235R(6,2)-Dij1235R(8,3))-Dij2345R(1,1)-Dij2345R(1,2)
     -   -s45*EijR(1,1)+s45*EijR(2,1)-4*(Dij1235R(6,3)-Dij1235R(10,3)
     -   -EijR(21,3)+EijR(22,3)-EijR(36,4)-EijR(37,4))-8*EijR(40,4)+(
     -   -EijR(7,2)+EijR(9,2))*P(7)+EijR(1,3)*P(50)+EijR(1,2)*P(77)+(
     -   EijR(6,2)+EijR(6,3)-EijR(8,2)+EijR(9,3)-2*EijR(17,3))*P(198)
     -   +EijR(2,2)*P(286)+EijR(2,3)*P(287)-EijR(5,2)*P(289)+EijR(5,3
     -   )*P(291)-(EijR(7,3)-2*EijR(18,3))*P(292)-EijR(8,3)*P(293)+Ei
     -   jR(10,3)*P(294))
       FI(26) = 2*(Dij1235I(1,1)+Dij1235I(1,2)-Dij1235I(2,1)-Dij1235I(
     -   2,2)+2*(Dij1235I(2,3)+Dij1235I(4,3)-Dij1235I(5,2)-Dij1235I(5
     -   ,3)+Dij1235I(6,2)-Dij1235I(8,3))-Dij2345I(1,1)-Dij2345I(1,2)
     -   -s45*EijI(1,1)+s45*EijI(2,1)-4*(Dij1235I(6,3)-Dij1235I(10,3)
     -   -EijI(21,3)+EijI(22,3)-EijI(36,4)-EijI(37,4))-8*EijI(40,4)+(
     -   -EijI(7,2)+EijI(9,2))*P(7)+EijI(1,3)*P(50)+EijI(1,2)*P(77)+(
     -   EijI(6,2)+EijI(6,3)-EijI(8,2)+EijI(9,3)-2*EijI(17,3))*P(198)
     -   +EijI(2,2)*P(286)+EijI(2,3)*P(287)-EijI(5,2)*P(289)+EijI(5,3
     -   )*P(291)-(EijI(7,3)-2*EijI(18,3))*P(292)-EijI(8,3)*P(293)+Ei
     -   jI(10,3)*P(294))
       F(26)=DCMPLX(FR(26),FI(26))
       FR(27) = 2*(Dij1235R(2,2)-Dij1235R(4,2)+Dij2345R(1,1)+Dij2345R(
     -   1,2)-6*(EijR(21,3)-EijR(22,3))+4*(EijR(37,4)-EijR(40,4))+(Ei
     -   jR(10,3)-EijR(18,3))*P(7)+(-EijR(9,3)+EijR(17,3))*P(198)+(Ei
     -   jR(2,2)-EijR(5,2))*P(242)+EijR(2,3)*P(244)+EijR(5,3)*P(245)-
     -   EijR(8,3)*P(247)+2*(Dij1235R(2,3)+Dij1235R(5,2)-Dij1235R(6,2
     -   )-Dij1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3)+p3sq*EijR(6,2)-
     -   p3sq*EijR(8,2)+EijR(9,2)*P(14)+EijR(7,2)*P(275)))
       FI(27) = 2*(Dij1235I(2,2)-Dij1235I(4,2)+Dij2345I(1,1)+Dij2345I(
     -   1,2)-6*(EijI(21,3)-EijI(22,3))+4*(EijI(37,4)-EijI(40,4))+(Ei
     -   jI(10,3)-EijI(18,3))*P(7)+(-EijI(9,3)+EijI(17,3))*P(198)+(Ei
     -   jI(2,2)-EijI(5,2))*P(242)+EijI(2,3)*P(244)+EijI(5,3)*P(245)-
     -   EijI(8,3)*P(247)+2*(Dij1235I(2,3)+Dij1235I(5,2)-Dij1235I(6,2
     -   )-Dij1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3)+p3sq*EijI(6,2)-
     -   p3sq*EijI(8,2)+EijI(9,2)*P(14)+EijI(7,2)*P(275)))
       F(27)=DCMPLX(FR(27),FI(27))
       P(295) = p2sq-p3sq+s15+s34+s45
       P(296) = p2sq-p5sq-s12+s15
       P(297) = p2sq-p3sq+s12+s34
       FR(28) = 4*(D02345R+Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(3,2)+D
     -   ij1235R(9,3)-Dij2345R(1,1)+Dij2345R(2,1)-p5sq*EijR(2,1)+s45*
     -   EijR(3,1)-6*EijR(11,2)+p3sq*EijR(12,3)+4*EijR(22,3)-2*(Dij12
     -   35R(6,2)+Dij1235R(8,3)+EijR(23,3)-EijR(37,4)+EijR(43,4))+Eij
     -   R(2,3)*P(1)+(EijR(4,1)-EijR(9,2)+EijR(10,2))*P(7)+EijR(17,3)
     -   *P(10)-EijR(3,2)*P(16)-EijR(6,2)*P(33)+EijR(8,3)*P(251)-EijR
     -   (9,3)*P(277)-EijR(2,2)*P(295)+EijR(5,2)*P(296)+EijR(8,2)*P(2
     -   97))
       FI(28) = 4*(D02345I+Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(3,2)+D
     -   ij1235I(9,3)-Dij2345I(1,1)+Dij2345I(2,1)-p5sq*EijI(2,1)+s45*
     -   EijI(3,1)-6*EijI(11,2)+p3sq*EijI(12,3)+4*EijI(22,3)-2*(Dij12
     -   35I(6,2)+Dij1235I(8,3)+EijI(23,3)-EijI(37,4)+EijI(43,4))+Eij
     -   I(2,3)*P(1)+(EijI(4,1)-EijI(9,2)+EijI(10,2))*P(7)+EijI(17,3)
     -   *P(10)-EijI(3,2)*P(16)-EijI(6,2)*P(33)+EijI(8,3)*P(251)-EijI
     -   (9,3)*P(277)-EijI(2,2)*P(295)+EijI(5,2)*P(296)+EijI(8,2)*P(2
     -   97))
       F(28)=DCMPLX(FR(28),FI(28))
       P(298) = p3sq-s45+3*P(6)
       P(299) = p2sq+p3sq
       P(300) = 4*s12-s23-5*s34-2*s45+3*P(299)
       P(301) = s45-3*P(6)
       P(302) = -s23+s45+3*P(33)
       P(303) = 3*p3sq+2*P(16)
       P(304) = 3*p2sq+p3sq-2*p5sq+4*s12-s23-5*s34
       P(305) = p3sq+2*P(16)
       FR(29) = 2*(-Dij1235R(2,1)+Dij1235R(2,2)+Dij1235R(3,1)-Dij1235R
     -   (3,2)-4*Dij1235R(8,3)-p3sq*(EijR(3,1)+2*(EijR(9,3)-EijR(20,3
     -   )))+2*(Dij1235R(2,3)+Dij1235R(9,3)+Dij2345R(1,1)-Dij2345R(3,
     -   1)+3*EijR(11,2)+EijR(22,3)+EijR(24,3)+2*(EijR(37,4)-EijR(44,
     -   4))+(EijR(2,3)-EijR(10,3))*P(1)-EijR(4,2)*P(7)+EijR(18,3)*P(
     -   10)+EijR(8,3)*P(251))+EijR(2,1)*P(298)+EijR(2,2)*P(300)+EijR
     -   (4,1)*P(301)+(-EijR(5,2)+EijR(7,2))*P(302)-EijR(8,2)*P(303)-
     -   EijR(9,2)*P(304)+EijR(10,2)*P(305))
       FI(29) = -2*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(3,1)+Dij1235I
     -   (3,2)+4*Dij1235I(8,3)+p3sq*(EijI(3,1)+2*(EijI(9,3)-EijI(20,3
     -   )))-2*(Dij1235I(2,3)+Dij1235I(9,3)+Dij2345I(1,1)-Dij2345I(3,
     -   1)+3*EijI(11,2)+EijI(22,3)+EijI(24,3)+2*(EijI(37,4)-EijI(44,
     -   4))-(-EijI(2,3)+EijI(10,3))*P(1)-EijI(4,2)*P(7)+EijI(18,3)*P
     -   (10)+EijI(8,3)*P(251))-EijI(2,1)*P(298)-EijI(2,2)*P(300)-Eij
     -   I(4,1)*P(301)+(EijI(5,2)-EijI(7,2))*P(302)+EijI(8,2)*P(303)+
     -   EijI(9,2)*P(304)-EijI(10,2)*P(305))
       F(29)=DCMPLX(FR(29),FI(29))
       P(306) = 3*s12+s23-2*s34-s45
       P(307) = p3sq+s12+s23-3*s34-2*s45
       P(308) = p3sq-3*s12-s34-s45
       P(309) = 2*p3sq+s12-s45
       P(310) = 2*p3sq-s12+s45
       P(311) = -p5sq+s45+3*P(14)
       P(312) = -p3sq-p5sq+s34+s45
       FR(30) = -2*(D01235R+Dij1235R(1,1)+Dij1235R(2,2)-Dij1235R(3,1)+
     -   Dij1235R(4,2)+Dij1235R(5,2)-3*Dij1235R(6,2)+Dij2345R(1,1)-6*
     -   EijR(21,3)+10*EijR(22,3)-4*(EijR(37,4)-EijR(40,4))-EE0R*P(6)
     -   -EijR(4,1)*P(7)+(EijR(3,1)-2*(EijR(9,3)-EijR(17,3)))*P(16)+(
     -   EijR(1,2)+2*EijR(5,3))*P(50)-EijR(7,2)*P(153)-EijR(1,1)*P(25
     -   9)-EijR(2,1)*P(306)-EijR(2,2)*P(307)+EijR(5,2)*P(308)+EijR(6
     -   ,2)*P(309)-EijR(8,2)*P(310)+EijR(9,2)*P(311)+2*(Dij1235R(2,1
     -   )-Dij1235R(2,3)+Dij1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)+D
     -   ij2345R(1,2)-EijR(11,2)-EijR(18,3)*P(101)-EijR(2,3)*P(258)+E
     -   ijR(8,3)*P(261)-EijR(10,3)*P(312)))
       FI(30) = -2*(D01235I+Dij1235I(1,1)+Dij1235I(2,2)-Dij1235I(3,1)+
     -   Dij1235I(4,2)+Dij1235I(5,2)-3*Dij1235I(6,2)+Dij2345I(1,1)-6*
     -   EijI(21,3)+10*EijI(22,3)-4*(EijI(37,4)-EijI(40,4))-EE0I*P(6)
     -   -EijI(4,1)*P(7)+(EijI(3,1)-2*(EijI(9,3)-EijI(17,3)))*P(16)+(
     -   EijI(1,2)+2*EijI(5,3))*P(50)-EijI(7,2)*P(153)-EijI(1,1)*P(25
     -   9)-EijI(2,1)*P(306)-EijI(2,2)*P(307)+EijI(5,2)*P(308)+EijI(6
     -   ,2)*P(309)-EijI(8,2)*P(310)+EijI(9,2)*P(311)+2*(Dij1235I(2,1
     -   )-Dij1235I(2,3)+Dij1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)+D
     -   ij2345I(1,2)-EijI(11,2)-EijI(18,3)*P(101)-EijI(2,3)*P(258)+E
     -   ijI(8,3)*P(261)-EijI(10,3)*P(312)))
       F(30)=DCMPLX(FR(30),FI(30))
       P(313) = p3sq+2*s12-3*s34-s45
       P(314) = -s23-5*s34+2*P(16)+3*P(299)
       P(315) = p3sq+2*s12-3*s34
       P(316) = 3*p2sq-s23-5*s34-2*s45+4*P(104)
       FR(31) = -2*(D02345R+Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(3,1)+
     -   Dij1235R(3,2)+(EijR(4,1)+2*EijR(10,2))*P(7)-2*(Dij1235R(2,3)
     -   +Dij1235R(9,3)-Dij2345R(2,1)+4*EijR(11,2)+p3sq*EijR(12,3)+3*
     -   EijR(22,3)+EijR(23,3)-2*(Dij1235R(8,3)-EijR(37,4)+EijR(43,4)
     -   )+EijR(2,3)*P(1)+EijR(17,3)*P(10)+EijR(8,3)*P(251)-EijR(9,3)
     -   *P(277))+(EijR(5,2)-EijR(6,2))*P(302)-EijR(3,2)*P(305)-EijR(
     -   2,1)*P(313)-EijR(2,2)*P(314)+EijR(3,1)*P(315)+EijR(8,2)*P(31
     -   6))
       FI(31) = -2*(D02345I+Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(3,1)+
     -   Dij1235I(3,2)+(EijI(4,1)+2*EijI(10,2))*P(7)-2*(Dij1235I(2,3)
     -   +Dij1235I(9,3)-Dij2345I(2,1)+4*EijI(11,2)+p3sq*EijI(12,3)+3*
     -   EijI(22,3)+EijI(23,3)-2*(Dij1235I(8,3)-EijI(37,4)+EijI(43,4)
     -   )+EijI(2,3)*P(1)+EijI(17,3)*P(10)+EijI(8,3)*P(251)-EijI(9,3)
     -   *P(277))+(EijI(5,2)-EijI(6,2))*P(302)-EijI(3,2)*P(305)-EijI(
     -   2,1)*P(313)-EijI(2,2)*P(314)+EijI(3,1)*P(315)+EijI(8,2)*P(31
     -   6))
       F(31)=DCMPLX(FR(31),FI(31))
       P(317) = p2sq-p3sq+s34+s45
       P(318) = p2sq-p3sq-p5sq+s34+2*s45
       FR(32) = 4*(Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(3,2)+Dij1235R(
     -   9,3)-Dij2345R(1,1)+Dij2345R(3,1)-s12*EijR(2,1)+p3sq*(EijR(3,
     -   1)-EijR(9,3)+EijR(20,3))-2*(Dij1235R(6,2)+Dij1235R(8,3)+2*Ei
     -   jR(11,2)-EijR(22,3)+EijR(24,3)-EijR(37,4)+EijR(44,4))+(EijR(
     -   2,3)-EijR(10,3))*P(1)+EijR(4,2)*P(7)+EijR(18,3)*P(10)+(EijR(
     -   8,2)-EijR(10,2))*P(16)+(EijR(5,2)-EijR(7,2))*P(33)-EijR(4,1)
     -   *P(161)+EijR(8,3)*P(251)-EijR(2,2)*P(317)+EijR(9,2)*P(318))
       FI(32) = 4*(Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(3,2)+Dij1235I(
     -   9,3)-Dij2345I(1,1)+Dij2345I(3,1)-s12*EijI(2,1)+p3sq*(EijI(3,
     -   1)-EijI(9,3)+EijI(20,3))-2*(Dij1235I(6,2)+Dij1235I(8,3)+2*Ei
     -   jI(11,2)-EijI(22,3)+EijI(24,3)-EijI(37,4)+EijI(44,4))+(EijI(
     -   2,3)-EijI(10,3))*P(1)+EijI(4,2)*P(7)+EijI(18,3)*P(10)+(EijI(
     -   8,2)-EijI(10,2))*P(16)+(EijI(5,2)-EijI(7,2))*P(33)-EijI(4,1)
     -   *P(161)+EijI(8,3)*P(251)-EijI(2,2)*P(317)+EijI(9,2)*P(318))
       F(32)=DCMPLX(FR(32),FI(32))
       P(319) = 4*s12-5*s34
       P(320) = 3*s12-s34
       P(321) = s12+Is12*s34*P(189)
       P(322) = -4*s34+s45+2*P(321)
       P(323) = p3sq*s12+6*s12**2-9*s12*s34+2*s34*P(189)
       P(324) = 2*s12-s34
       P(325) = 3*s12-4*s34
       P(326) = p3sq*P(324)+P(16)*P(325)
       P(327) = -s12+s45
       P(328) = p3sq-s34+3*P(7)
       P(329) = p2sq+p5sq-s12+s15
       P(330) = p2sq+p5sq-s12-s15
       P(331) = s34*P(16)+p3sq*P(43)
       P(332) = p5sq-s15+2*P(33)
       P(333) = 3*s12-7*s34
       P(334) = s15+s23+s34-s45
       P(335) = p5sq+P(334)
       P(336) = 3*p2sq-s12-2*P(335)
       P(337) = p2sq-3*s12+2*P(264)
       P(338) = p2sq-s23-3*P(16)
       P(339) = -s15+s23+3*P(7)+2*P(14)
       FR(33) = -C0235R-Cij235R(1,1)+Cij235R(2,1)-6*Dij2345R(7,2)+8*Ei
     -   jR(46,4)+Dij1345R(5,2)*P(14)+P(7)*(Dij2345R(3,1)+EijR(4,2)*P
     -   (21))+(Dij1235R(2,1)+Dij1235R(2,2))*P(33)+4*(Dij1235R(7,2)+D
     -   ij1235R(12,3)-Dij1235R(13,3)+EijR(21,3)*P(50))+Dij2345R(1,1)
     -   *P(191)+Dij1345R(1,2)*P(198)+EijR(2,1)*P(202)-EijR(2,2)*P(20
     -   5)-EijR(3,1)*P(207)+EijR(4,1)*P(211)+EijR(5,2)*P(214)-EijR(6
     -   ,2)*P(218)-EijR(7,2)*P(220)+EijR(8,2)*P(224)-EijR(9,2)*P(229
     -   )+EijR(10,2)*P(233)+D02345R*P(243)+EijR(3,2)*P(16)*P(251)+D0
     -   1345R*P(322)+Dij2345R(2,1)*P(327)-Dij1235R(3,1)*P(329)+Dij12
     -   35R(3,2)*P(330)-Dij1235R(6,2)*P(332)+Is12*(-4*Cij135R(4,2)+I
     -   s45*P(6)*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)*P(95)+C0145R*P
     -   (185))+s34*Dij1345R(6,2)*P(275)+C0345R*P(319)+Cij135R(2,1)*P
     -   (320)+Dij1345R(1,1)*P(323)-Dij1345R(2,1)*P(326)+Dij1345R(3,1
     -   )*P(6)*P(328)-Dij1345R(4,2)*P(331)-2*(B035R-p3sq*s34*Dij1345
     -   R(2,2)+C0135R*P(35)+Cij135R(1,1)*P(59)+Dij1345R(7,2)*P(333))
     -   )-2*(P(3)*(D01235R+Dij1235R(1,1)-EE0R*P(6)+EijR(1,2)*P(33)+E
     -   ijR(1,1)*P(71))+EijR(11,2)*P(336)-EijR(22,3)*P(337)+EijR(23,
     -   3)*P(338)+EijR(24,3)*P(339))
       FI(33) = -C0235I-Cij235I(1,1)+Cij235I(2,1)-6*Dij2345I(7,2)+8*Ei
     -   jI(46,4)+Dij1345I(5,2)*P(14)+P(7)*(Dij2345I(3,1)+EijI(4,2)*P
     -   (21))+(Dij1235I(2,1)+Dij1235I(2,2))*P(33)+4*(Dij1235I(7,2)+D
     -   ij1235I(12,3)-Dij1235I(13,3)+EijI(21,3)*P(50))+Dij2345I(1,1)
     -   *P(191)+Dij1345I(1,2)*P(198)+EijI(2,1)*P(202)-EijI(2,2)*P(20
     -   5)-EijI(3,1)*P(207)+EijI(4,1)*P(211)+EijI(5,2)*P(214)-EijI(6
     -   ,2)*P(218)-EijI(7,2)*P(220)+EijI(8,2)*P(224)-EijI(9,2)*P(229
     -   )+EijI(10,2)*P(233)+D02345I*P(243)+EijI(3,2)*P(16)*P(251)+D0
     -   1345I*P(322)+Dij2345I(2,1)*P(327)-Dij1235I(3,1)*P(329)+Dij12
     -   35I(3,2)*P(330)-Dij1235I(6,2)*P(332)+Is12*(-4*Cij135I(4,2)+I
     -   s45*P(6)*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)*P(95)+C0145I*P
     -   (185))+s34*Dij1345I(6,2)*P(275)+C0345I*P(319)+Cij135I(2,1)*P
     -   (320)+Dij1345I(1,1)*P(323)-Dij1345I(2,1)*P(326)+Dij1345I(3,1
     -   )*P(6)*P(328)-Dij1345I(4,2)*P(331)-2*(B035I-p3sq*s34*Dij1345
     -   I(2,2)+C0135I*P(35)+Cij135I(1,1)*P(59)+Dij1345I(7,2)*P(333))
     -   )-2*(P(3)*(D01235I+Dij1235I(1,1)-EE0I*P(6)+EijI(1,2)*P(33)+E
     -   ijI(1,1)*P(71))+EijI(11,2)*P(336)-EijI(22,3)*P(337)+EijI(23,
     -   3)*P(338)+EijI(24,3)*P(339))
       F(33)=DCMPLX(FR(33),FI(33))
       P(340) = p3sq+s12-3*s34-s45
       P(341) = p2sq+3*p3sq-s23-5*s34+2*P(16)
       FR(34) = 2*(-Dij1235R(2,1)+Dij1235R(2,2)+Dij1235R(3,1)-Dij1235R
     -   (6,2)+p3sq*EijR(3,1)-p3sq*EijR(8,2)+4*EijR(37,4)-EijR(5,2)*P
     -   (10)-EijR(4,1)*P(14)+2*(Dij1235R(2,3)-Dij1235R(8,3)-p3sq*Eij
     -   R(9,3)-EijR(11,2)+EijR(22,3)+EijR(2,3)*P(1)+EijR(8,3)*P(251)
     -   )+EijR(2,1)*P(340)+EijR(2,2)*P(341))
       FI(34) = 2*(-Dij1235I(2,1)+Dij1235I(2,2)+Dij1235I(3,1)-Dij1235I
     -   (6,2)+p3sq*EijI(3,1)-p3sq*EijI(8,2)+4*EijI(37,4)-EijI(5,2)*P
     -   (10)-EijI(4,1)*P(14)+2*(Dij1235I(2,3)-Dij1235I(8,3)-p3sq*Eij
     -   I(9,3)-EijI(11,2)+EijI(22,3)+EijI(2,3)*P(1)+EijI(8,3)*P(251)
     -   )+EijI(2,1)*P(340)+EijI(2,2)*P(341))
       F(34)=DCMPLX(FR(34),FI(34))
       P(342) = p5sq-s12+s34
       P(343) = -s12+s15+s34
       P(344) = p5sq+s12-2*s15-s34
       P(345) = s15+s45-2*P(6)
       P(346) = s12-P(19)
       P(347) = p5sq-s45+2*P(346)
       FR(35) = 4*(-2*(EijR(21,3)-EijR(22,3)-EijR(37,4)+EijR(40,4)-Eij
     -   R(41,4)+EijR(43,4))+(-EijR(5,3)+EijR(6,3))*P(3)+(-EijR(10,3)
     -   +EijR(18,3)-EijR(19,3)+EijR(20,3))*P(7)+(EijR(11,3)-EijR(12,
     -   3))*P(156)+(EijR(2,2)-EijR(5,2)+EijR(6,2)-EijR(8,2))*P(342)+
     -   EijR(2,3)*P(343)+EijR(8,3)*P(344)-EijR(9,3)*P(345)-EijR(17,3
     -   )*P(347))
       FI(35) = 4*(-2*(EijI(21,3)-EijI(22,3)-EijI(37,4)+EijI(40,4)-Eij
     -   I(41,4)+EijI(43,4))+(-EijI(5,3)+EijI(6,3))*P(3)+(-EijI(10,3)
     -   +EijI(18,3)-EijI(19,3)+EijI(20,3))*P(7)+(EijI(11,3)-EijI(12,
     -   3))*P(156)+(EijI(2,2)-EijI(5,2)+EijI(6,2)-EijI(8,2))*P(342)+
     -   EijI(2,3)*P(343)+EijI(8,3)*P(344)-EijI(9,3)*P(345)-EijI(17,3
     -   )*P(347))
       F(35)=DCMPLX(FR(35),FI(35))
       P(348) = p5sq-s12+s15+s34-s45
       P(349) = s12-2*s15-s34+s45
       FR(36) = 4*(-2*(EijR(21,3)-EijR(22,3)-EijR(37,4)+EijR(40,4)-Eij
     -   R(42,4)+EijR(44,4))+(-EijR(5,3)+EijR(7,3))*P(3)+(-EijR(14,3)
     -   +EijR(15,3))*P(7)+(EijR(9,3)-EijR(17,3)+EijR(19,3)-EijR(20,3
     -   ))*P(156)+(EijR(2,2)-EijR(5,2)+EijR(7,2)-EijR(9,2))*P(342)+E
     -   ijR(2,3)*P(343)+EijR(8,3)*P(344)-EijR(10,3)*P(348)-EijR(18,3
     -   )*P(349))
       FI(36) = 4*(-2*(EijI(21,3)-EijI(22,3)-EijI(37,4)+EijI(40,4)-Eij
     -   I(42,4)+EijI(44,4))+(-EijI(5,3)+EijI(7,3))*P(3)+(-EijI(14,3)
     -   +EijI(15,3))*P(7)+(EijI(9,3)-EijI(17,3)+EijI(19,3)-EijI(20,3
     -   ))*P(156)+(EijI(2,2)-EijI(5,2)+EijI(7,2)-EijI(9,2))*P(342)+E
     -   ijI(2,3)*P(343)+EijI(8,3)*P(344)-EijI(10,3)*P(348)-EijI(18,3
     -   )*P(349))
       F(36)=DCMPLX(FR(36),FI(36))
       P(350) = p5sq-2*s12+s34
       FR(37) = -4*(D02345R+Dij2345R(1,1)+4*(EijR(21,3)-EijR(22,3))-2*
     -   (EijR(37,4)-EijR(40,4))+EijR(5,3)*P(3)+(EijR(10,3)-EijR(18,3
     -   ))*P(7)+(-EijR(7,2)+EijR(9,2))*P(101)+(EijR(6,2)-EijR(8,2))*
     -   P(140)+(-EijR(9,3)+EijR(17,3))*P(156)-EijR(2,3)*P(343)-EijR(
     -   8,3)*P(344)+(-EijR(2,2)+EijR(5,2))*P(350))
       FI(37) = -4*(D02345I+Dij2345I(1,1)+4*(EijI(21,3)-EijI(22,3))-2*
     -   (EijI(37,4)-EijI(40,4))+EijI(5,3)*P(3)+(EijI(10,3)-EijI(18,3
     -   ))*P(7)+(-EijI(7,2)+EijI(9,2))*P(101)+(EijI(6,2)-EijI(8,2))*
     -   P(140)+(-EijI(9,3)+EijI(17,3))*P(156)-EijI(2,3)*P(343)-EijI(
     -   8,3)*P(344)+(-EijI(2,2)+EijI(5,2))*P(350))
       F(37)=DCMPLX(FR(37),FI(37))
       P(351) = p2sq+s12-2*s15-s34
       P(352) = p2sq+2*p5sq+s12-4*s15-s34
       P(353) = -p2sq+p5sq-s12+s15+s34
       P(354) = p2sq-p5sq+s12-s15-s34
       P(355) = p2sq-2*p5sq-s12+s45
       P(356) = p2sq+p3sq-s23-2*P(19)
       P(357) = p2sq-s12-2*s15+s45
       P(358) = -s23+s45+2*P(3)
       P(359) = -p2sq+s12+s23+2*P(3)
       P(360) = p2sq-s12-s23-2*P(3)
       P(361) = 2*s15+s23+s34
       P(362) = p5sq-P(361)
       P(363) = p2sq+p3sq+s45+2*P(362)
       P(364) = p2sq+p3sq-s15-s23-2*s34
       P(365) = p2sq+p3sq-s23-2*s34
       P(366) = p2sq+p3sq+p5sq
       P(367) = -3*s23+s45-4*P(19)+2*P(366)
       FR(38) = 2*(D02345R+3*Dij2345R(1,1)-6*(EijR(11,2)+EijR(22,3)-Ei
     -   jR(23,3))+4*(EijR(37,4)-EijR(40,4)+EijR(41,4)-EijR(43,4))+(E
     -   ijR(7,2)-EijR(9,2))*P(23)+(-EijR(10,3)+EijR(18,3)-EijR(19,3)
     -   +EijR(20,3))*P(66)+Is34*(C0235R+Cij235R(1,1)-6*Dij1235R(7,2)
     -   -4*(Dij1235R(11,3)-Dij1235R(12,3))-2*(Dij1235R(1,1)+Dij1235R
     -   (1,2)-Dij1235R(2,1))*P(3)-Dij1235R(2,2)*P(351)+Dij1235R(4,2)
     -   *P(352)+Dij1235R(5,2)*P(353)+Dij1235R(6,2)*P(354))-EijR(2,2)
     -   *P(355)-EijR(2,3)*P(356)+EijR(5,2)*P(357)+(-EijR(5,3)+EijR(6
     -   ,3))*P(358)+EijR(6,2)*P(359)+EijR(8,2)*P(360)+EijR(8,3)*P(36
     -   3)+2*(Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(4,2)-(EijR(1,1)+E
     -   ijR(1,2)-EijR(2,1))*P(3)+EijR(9,3)*P(364))+(EijR(11,3)-EijR(
     -   12,3))*P(365)-EijR(17,3)*P(367))
       FI(38) = 2*(D02345I+3*Dij2345I(1,1)-6*(EijI(11,2)+EijI(22,3)-Ei
     -   jI(23,3))+4*(EijI(37,4)-EijI(40,4)+EijI(41,4)-EijI(43,4))+(E
     -   ijI(7,2)-EijI(9,2))*P(23)+(-EijI(10,3)+EijI(18,3)-EijI(19,3)
     -   +EijI(20,3))*P(66)+Is34*(C0235I+Cij235I(1,1)-6*Dij1235I(7,2)
     -   -4*(Dij1235I(11,3)-Dij1235I(12,3))-2*(Dij1235I(1,1)+Dij1235I
     -   (1,2)-Dij1235I(2,1))*P(3)-Dij1235I(2,2)*P(351)+Dij1235I(4,2)
     -   *P(352)+Dij1235I(5,2)*P(353)+Dij1235I(6,2)*P(354))-EijI(2,2)
     -   *P(355)-EijI(2,3)*P(356)+EijI(5,2)*P(357)+(-EijI(5,3)+EijI(6
     -   ,3))*P(358)+EijI(6,2)*P(359)+EijI(8,2)*P(360)+EijI(8,3)*P(36
     -   3)+2*(Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(4,2)-(EijI(1,1)+E
     -   ijI(1,2)-EijI(2,1))*P(3)+EijI(9,3)*P(364))+(EijI(11,3)-EijI(
     -   12,3))*P(365)-EijI(17,3)*P(367))
       F(38)=DCMPLX(FR(38),FI(38))
       P(368) = p2sq+s12-s23-2*s34-s45
       P(369) = 3*p5sq-s15-s23+2*P(156)
       P(370) = p2sq+p3sq-3*s15-2*s34
       P(371) = p2sq+p3sq-5*s15-s23+s45+2*P(39)
       FR(39) = 2*(D02345R+3*Dij2345R(1,1)-6*(EijR(11,2)+EijR(21,3)-Ei
     -   jR(24,3))+4*(EijR(37,4)-EijR(40,4)+EijR(42,4)-EijR(44,4))+2*
     -   (Dij2345R(1,2)-Dij2345R(3,1)-Dij2345R(5,2)-(EijR(1,1)+EijR(1
     -   ,2)-EijR(2,1))*P(3))+(-EijR(14,3)+EijR(15,3))*P(66)+Is34*(C0
     -   235R+Cij235R(1,1)-6*Dij1235R(7,2)-4*(Dij1235R(11,3)-Dij1235R
     -   (12,3))-2*(Dij1235R(1,1)+Dij1235R(1,2)-Dij1235R(2,1))*P(3)-D
     -   ij1235R(2,2)*P(351)+Dij1235R(4,2)*P(352)+Dij1235R(5,2)*P(353
     -   )+Dij1235R(6,2)*P(354))-EijR(2,2)*P(355)-EijR(2,3)*P(356)+Ei
     -   jR(5,2)*P(357)+(-EijR(5,3)+EijR(7,3))*P(358)+EijR(8,3)*P(363
     -   )+(EijR(9,3)-EijR(17,3)+EijR(19,3)-EijR(20,3))*P(365)+(-EijR
     -   (6,2)+EijR(8,2))*P(368)+(EijR(7,2)-EijR(9,2))*P(369)+EijR(10
     -   ,3)*P(370)-EijR(18,3)*P(371))
       FI(39) = 2*(D02345I+3*Dij2345I(1,1)-6*(EijI(11,2)+EijI(21,3)-Ei
     -   jI(24,3))+4*(EijI(37,4)-EijI(40,4)+EijI(42,4)-EijI(44,4))+2*
     -   (Dij2345I(1,2)-Dij2345I(3,1)-Dij2345I(5,2)-(EijI(1,1)+EijI(1
     -   ,2)-EijI(2,1))*P(3))+(-EijI(14,3)+EijI(15,3))*P(66)+Is34*(C0
     -   235I+Cij235I(1,1)-6*Dij1235I(7,2)-4*(Dij1235I(11,3)-Dij1235I
     -   (12,3))-2*(Dij1235I(1,1)+Dij1235I(1,2)-Dij1235I(2,1))*P(3)-D
     -   ij1235I(2,2)*P(351)+Dij1235I(4,2)*P(352)+Dij1235I(5,2)*P(353
     -   )+Dij1235I(6,2)*P(354))-EijI(2,2)*P(355)-EijI(2,3)*P(356)+Ei
     -   jI(5,2)*P(357)+(-EijI(5,3)+EijI(7,3))*P(358)+EijI(8,3)*P(363
     -   )+(EijI(9,3)-EijI(17,3)+EijI(19,3)-EijI(20,3))*P(365)+(-EijI
     -   (6,2)+EijI(8,2))*P(368)+(EijI(7,2)-EijI(9,2))*P(369)+EijI(10
     -   ,3)*P(370)-EijI(18,3)*P(371))
       F(39)=DCMPLX(FR(39),FI(39))
       P(372) = 2*p5sq-s12+s34
       P(373) = p2sq-s12+2*s15+s34
       P(374) = -p2sq+s12+2*P(3)
       P(375) = p2sq-s12-2*P(3)
       P(376) = s12-s15
       P(377) = p2sq+p5sq+3*s34-2*P(376)
       P(378) = p3sq+s12-s23+s45-2*P(234)
       P(379) = p3sq+s12-s23-2*P(19)
       P(380) = p3sq+s12-2*P(102)
       P(381) = s15+s34-s45
       P(382) = p3sq-3*p5sq+s12-s23-2*P(381)
       P(383) = p3sq+s12-s15
       P(384) = -s23-4*s34-s45+2*P(383)
       P(385) = p3sq-p5sq+s12-s23+s45-2*P(19)
       P(386) = p5sq+s23+4*s34-2*P(383)
       FR(40) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(3,1)-Dij2345R(
     -   4,2)-Dij2345R(5,2)+Dij2345R(6,2)+6*(EijR(11,2)+EijR(22,3)-Ei
     -   jR(23,3))+4*(EijR(37,4)-EijR(43,4)-EijR(44,4)+EijR(45,4))+(E
     -   ijR(4,2)+EijR(15,3)-EijR(16,3))*P(7)+(-EijR(12,3)+EijR(13,3)
     -   )*P(79)+(EijR(2,1)-EijR(4,1))*P(185)+(EijR(5,2)-EijR(7,2)+Ei
     -   jR(8,3)-EijR(17,3)-EijR(18,3)+EijR(19,3))*P(358)+Is34*(Cij23
     -   5R(1,1)-Cij235R(2,1)+6*Dij1235R(7,2)+4*(Dij1235R(12,3)-Dij12
     -   35R(13,3))+Dij1235R(3,2)*P(73)+(Dij1235R(2,1)-Dij1235R(3,1))
     -   *P(372)+Dij1235R(2,2)*P(373)+Dij1235R(4,2)*P(374)+Dij1235R(5
     -   ,2)*P(375)-Dij1235R(6,2)*P(377))-EijR(2,2)*P(378)-EijR(2,3)*
     -   P(379)+(EijR(8,2)-EijR(10,2))*P(380)+EijR(9,2)*P(382)+EijR(9
     -   ,3)*P(384)+EijR(10,3)*P(385)+EijR(20,3)*P(386))
       FI(40) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(3,1)-Dij2345I(
     -   4,2)-Dij2345I(5,2)+Dij2345I(6,2)+6*(EijI(11,2)+EijI(22,3)-Ei
     -   jI(23,3))+4*(EijI(37,4)-EijI(43,4)-EijI(44,4)+EijI(45,4))+(E
     -   ijI(4,2)+EijI(15,3)-EijI(16,3))*P(7)+(-EijI(12,3)+EijI(13,3)
     -   )*P(79)+(EijI(2,1)-EijI(4,1))*P(185)+(EijI(5,2)-EijI(7,2)+Ei
     -   jI(8,3)-EijI(17,3)-EijI(18,3)+EijI(19,3))*P(358)+Is34*(Cij23
     -   5I(1,1)-Cij235I(2,1)+6*Dij1235I(7,2)+4*(Dij1235I(12,3)-Dij12
     -   35I(13,3))+Dij1235I(3,2)*P(73)+(Dij1235I(2,1)-Dij1235I(3,1))
     -   *P(372)+Dij1235I(2,2)*P(373)+Dij1235I(4,2)*P(374)+Dij1235I(5
     -   ,2)*P(375)-Dij1235I(6,2)*P(377))-EijI(2,2)*P(378)-EijI(2,3)*
     -   P(379)+(EijI(8,2)-EijI(10,2))*P(380)+EijI(9,2)*P(382)+EijI(9
     -   ,3)*P(384)+EijI(10,3)*P(385)+EijI(20,3)*P(386))
       F(40)=DCMPLX(FR(40),FI(40))
       P(387) = -2*p5sq+s45
       P(388) = -s23+s45-4*P(102)+2*P(383)
       P(389) = -4*s15-2*s23-6*s34-s45+3*P(104)
       P(390) = -s23-6*s34+3*P(104)-2*P(208)
       FR(41) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(2,1)+Dij2345R(
     -   2,2)-2*Dij2345R(4,2)+4*(EijR(22,3)-EijR(23,3)+EijR(37,4)+Eij
     -   R(38,4)-2*EijR(43,4))-(EijR(9,2)-EijR(10,2)+EijR(10,3)+EijR(
     -   13,3)-2*EijR(20,3))*P(7)+EijR(3,3)*P(79)+EijR(2,1)*P(185)+(E
     -   ijR(5,2)-EijR(6,2)+EijR(8,3)+EijR(11,3)-2*EijR(17,3))*P(358)
     -   +Is34*(Cij235R(1,1)-Cij235R(2,1)+6*Dij1235R(7,2)+4*(Dij1235R
     -   (12,3)-Dij1235R(13,3))+Dij1235R(3,2)*P(73)+(Dij1235R(2,1)-Di
     -   j1235R(3,1))*P(372)+Dij1235R(2,2)*P(373)+Dij1235R(4,2)*P(374
     -   )+Dij1235R(5,2)*P(375)-Dij1235R(6,2)*P(377))-EijR(2,2)*P(378
     -   )-EijR(2,3)*P(379)-EijR(3,2)*P(380)+EijR(3,1)*P(387)+EijR(8,
     -   2)*P(388)+EijR(9,3)*P(389)-EijR(12,3)*P(390))
       FI(41) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(2,1)+Dij2345I(
     -   2,2)-2*Dij2345I(4,2)+4*(EijI(22,3)-EijI(23,3)+EijI(37,4)+Eij
     -   I(38,4)-2*EijI(43,4))-(EijI(9,2)-EijI(10,2)+EijI(10,3)+EijI(
     -   13,3)-2*EijI(20,3))*P(7)+EijI(3,3)*P(79)+EijI(2,1)*P(185)+(E
     -   ijI(5,2)-EijI(6,2)+EijI(8,3)+EijI(11,3)-2*EijI(17,3))*P(358)
     -   +Is34*(Cij235I(1,1)-Cij235I(2,1)+6*Dij1235I(7,2)+4*(Dij1235I
     -   (12,3)-Dij1235I(13,3))+Dij1235I(3,2)*P(73)+(Dij1235I(2,1)-Di
     -   j1235I(3,1))*P(372)+Dij1235I(2,2)*P(373)+Dij1235I(4,2)*P(374
     -   )+Dij1235I(5,2)*P(375)-Dij1235I(6,2)*P(377))-EijI(2,2)*P(378
     -   )-EijI(2,3)*P(379)-EijI(3,2)*P(380)+EijI(3,1)*P(387)+EijI(8,
     -   2)*P(388)+EijI(9,3)*P(389)-EijI(12,3)*P(390))
       F(41)=DCMPLX(FR(41),FI(41))
       P(391) = p3sq-5*p5sq+s12-s23+3*s45-2*P(19)
       P(392) = s23+2*P(19)
       P(393) = p3sq+s12-P(392)
       P(394) = -p5sq+s45+2*P(393)
       P(395) = p5sq+s15+s34-s45
       P(396) = p3sq+s12-s23-2*P(395)
       FR(42) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(3,1)+Dij2345R(
     -   3,2)-2*Dij2345R(5,2)+6*EijR(11,2)+12*(EijR(22,3)-EijR(24,3))
     -   +4*(EijR(37,4)+EijR(39,4))-8*EijR(44,4)-EijR(4,3)*P(7)+EijR(
     -   4,2)*P(54)+(EijR(8,2)+EijR(9,3)-EijR(10,2)+EijR(16,3)-2*EijR
     -   (20,3))*P(79)+EijR(2,1)*P(185)+(EijR(5,2)-EijR(7,2)+EijR(8,3
     -   )+EijR(14,3)-2*EijR(18,3))*P(358)+Is34*(Cij235R(1,1)-Cij235R
     -   (2,1)+6*Dij1235R(7,2)+4*(Dij1235R(12,3)-Dij1235R(13,3))+Dij1
     -   235R(3,2)*P(73)+(Dij1235R(2,1)-Dij1235R(3,1))*P(372)+Dij1235
     -   R(2,2)*P(373)+Dij1235R(4,2)*P(374)+Dij1235R(5,2)*P(375)-Dij1
     -   235R(6,2)*P(377))-EijR(2,2)*P(378)-EijR(2,3)*P(379)+EijR(4,1
     -   )*P(387)+EijR(9,2)*P(391)+EijR(10,3)*P(394)-EijR(15,3)*P(396
     -   ))
       FI(42) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(3,1)+Dij2345I(
     -   3,2)-2*Dij2345I(5,2)+6*EijI(11,2)+12*(EijI(22,3)-EijI(24,3))
     -   +4*(EijI(37,4)+EijI(39,4))-8*EijI(44,4)-EijI(4,3)*P(7)+EijI(
     -   4,2)*P(54)+(EijI(8,2)+EijI(9,3)-EijI(10,2)+EijI(16,3)-2*EijI
     -   (20,3))*P(79)+EijI(2,1)*P(185)+(EijI(5,2)-EijI(7,2)+EijI(8,3
     -   )+EijI(14,3)-2*EijI(18,3))*P(358)+Is34*(Cij235I(1,1)-Cij235I
     -   (2,1)+6*Dij1235I(7,2)+4*(Dij1235I(12,3)-Dij1235I(13,3))+Dij1
     -   235I(3,2)*P(73)+(Dij1235I(2,1)-Dij1235I(3,1))*P(372)+Dij1235
     -   I(2,2)*P(373)+Dij1235I(4,2)*P(374)+Dij1235I(5,2)*P(375)-Dij1
     -   235I(6,2)*P(377))-EijI(2,2)*P(378)-EijI(2,3)*P(379)+EijI(4,1
     -   )*P(387)+EijI(9,2)*P(391)+EijI(10,3)*P(394)-EijI(15,3)*P(396
     -   ))
       F(42)=DCMPLX(FR(42),FI(42))
       P(397) = p3sq-p5sq+s12-s15
       P(398) = -s23-4*s34+2*P(397)
       FR(43) = 2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(
     -   4,2)-Dij2345R(5,2)+Dij2345R(6,2)+10*EijR(22,3)-6*EijR(23,3)-
     -   4*(EijR(24,3)-EijR(37,4)+EijR(43,4)+EijR(44,4)-EijR(45,4))+(
     -   EijR(15,3)-EijR(16,3))*P(7)+(-EijR(9,2)+EijR(10,2))*P(54)-(E
     -   ijR(3,2)+EijR(12,3)-EijR(13,3))*P(79)+EijR(2,1)*P(185)+(EijR
     -   (5,2)-EijR(6,2)+EijR(8,3)-EijR(17,3)-EijR(18,3)+EijR(19,3))*
     -   P(358)+Is34*(Cij235R(1,1)-Cij235R(2,1)+6*Dij1235R(7,2)+4*(Di
     -   j1235R(12,3)-Dij1235R(13,3))+Dij1235R(3,2)*P(73)+(Dij1235R(2
     -   ,1)-Dij1235R(3,1))*P(372)+Dij1235R(2,2)*P(373)+Dij1235R(4,2)
     -   *P(374)+Dij1235R(5,2)*P(375)-Dij1235R(6,2)*P(377))-EijR(2,2)
     -   *P(378)-EijR(2,3)*P(379)+EijR(9,3)*P(384)+EijR(10,3)*P(385)+
     -   EijR(20,3)*P(386)+EijR(3,1)*P(387)+EijR(8,2)*P(398))
       FI(43) = 2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(
     -   4,2)-Dij2345I(5,2)+Dij2345I(6,2)+10*EijI(22,3)-6*EijI(23,3)-
     -   4*(EijI(24,3)-EijI(37,4)+EijI(43,4)+EijI(44,4)-EijI(45,4))+(
     -   EijI(15,3)-EijI(16,3))*P(7)+(-EijI(9,2)+EijI(10,2))*P(54)-(E
     -   ijI(3,2)+EijI(12,3)-EijI(13,3))*P(79)+EijI(2,1)*P(185)+(EijI
     -   (5,2)-EijI(6,2)+EijI(8,3)-EijI(17,3)-EijI(18,3)+EijI(19,3))*
     -   P(358)+Is34*(Cij235I(1,1)-Cij235I(2,1)+6*Dij1235I(7,2)+4*(Di
     -   j1235I(12,3)-Dij1235I(13,3))+Dij1235I(3,2)*P(73)+(Dij1235I(2
     -   ,1)-Dij1235I(3,1))*P(372)+Dij1235I(2,2)*P(373)+Dij1235I(4,2)
     -   *P(374)+Dij1235I(5,2)*P(375)-Dij1235I(6,2)*P(377))-EijI(2,2)
     -   *P(378)-EijI(2,3)*P(379)+EijI(9,3)*P(384)+EijI(10,3)*P(385)+
     -   EijI(20,3)*P(386)+EijI(3,1)*P(387)+EijI(8,2)*P(398))
       F(43)=DCMPLX(FR(43),FI(43))
       P(399) = p2sq+s12-2*P(19)
       P(400) = p2sq+4*p5sq+s12-6*s15-2*s34
       P(401) = p2sq+p5sq+s12-3*s15-2*s34
       FR(44) = 2*(D02345R+Dij2345R(1,2)+4*(EijR(36,4)+EijR(37,4)+2*(E
     -   ijR(21,3)-EijR(22,3)-EijR(40,4)))-(EijR(7,3)+EijR(10,3)-2*Ei
     -   jR(18,3))*P(23)+(EijR(6,3)+EijR(9,3)-2*EijR(17,3))*P(368)-Ei
     -   jR(2,3)*P(399)-EijR(5,3)*P(400)+2*(Dij2345R(1,1)+(EijR(1,2)+
     -   EijR(1,3)+EijR(2,2)-2*EijR(5,2))*P(3)+EijR(8,3)*P(401)))
       FI(44) = 2*(D02345I+Dij2345I(1,2)+4*(EijI(36,4)+EijI(37,4)+2*(E
     -   ijI(21,3)-EijI(22,3)-EijI(40,4)))-(EijI(7,3)+EijI(10,3)-2*Ei
     -   jI(18,3))*P(23)+(EijI(6,3)+EijI(9,3)-2*EijI(17,3))*P(368)-Ei
     -   jI(2,3)*P(399)-EijI(5,3)*P(400)+2*(Dij2345I(1,1)+(EijI(1,2)+
     -   EijI(1,3)+EijI(2,2)-2*EijI(5,2))*P(3)+EijI(8,3)*P(401)))
       F(44)=DCMPLX(FR(44),FI(44))
       P(402) = p2sq-s12+s45-2*P(102)
       P(403) = p2sq-s12+s45-2*P(19)
       P(404) = p3sq+s34
       P(405) = p2sq+s12-s23-s45-2*P(404)
       P(406) = -p5sq-s15+s23+s45+2*P(14)
       FR(45) = 2*(D02345R+4*(Dij1235R(5,2)-Dij1235R(6,2))+3*(Dij1235R
     -   (2,2)-Dij1235R(4,2)+Dij2345R(1,1))+2*(Dij2345R(1,2)-3*EijR(1
     -   1,2)-5*EijR(21,3)+2*(EijR(22,3)+EijR(37,4)-EijR(40,4))-(EijR
     -   (1,1)+EijR(1,2)-EijR(2,1))*P(3))+(-EijR(10,3)+EijR(18,3))*P(
     -   66)-EijR(2,3)*P(356)-EijR(5,3)*P(358)+EijR(8,3)*P(363)+(EijR
     -   (9,3)-EijR(17,3))*P(365)-EijR(2,2)*P(402)+EijR(5,2)*P(403)+(
     -   -EijR(6,2)+EijR(8,2))*P(405)+(-EijR(7,2)+EijR(9,2))*P(406))
       FI(45) = 2*(D02345I+4*(Dij1235I(5,2)-Dij1235I(6,2))+3*(Dij1235I
     -   (2,2)-Dij1235I(4,2)+Dij2345I(1,1))+2*(Dij2345I(1,2)-3*EijI(1
     -   1,2)-5*EijI(21,3)+2*(EijI(22,3)+EijI(37,4)-EijI(40,4))-(EijI
     -   (1,1)+EijI(1,2)-EijI(2,1))*P(3))+(-EijI(10,3)+EijI(18,3))*P(
     -   66)-EijI(2,3)*P(356)-EijI(5,3)*P(358)+EijI(8,3)*P(363)+(EijI
     -   (9,3)-EijI(17,3))*P(365)-EijI(2,2)*P(402)+EijI(5,2)*P(403)+(
     -   -EijI(6,2)+EijI(8,2))*P(405)+(-EijI(7,2)+EijI(9,2))*P(406))
       F(45)=DCMPLX(FR(45),FI(45))
       P(407) = -p5sq+s15
       P(408) = 2*p5sq-3*s12+s34
       P(409) = p2sq+p5sq+s15+s34
       P(410) = -p3sq-s12+s23-s45+2*P(409)
       P(411) = p5sq+s12-s34-s45
       P(412) = -s23+s45-2*P(296)
       P(413) = p3sq+s12-4*s34-s45
       P(414) = p2sq-s45
       P(415) = 3*p5sq+2*P(414)
       FR(46) = -4*Is15s34*(B025R+2*Cij125R(4,2)+C0125R*P(3)-Cij125R(1
     -   ,1)*P(407))+2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(5,2)+Eij
     -   R(15,3)*P(7)+(EijR(9,3)-EijR(20,3))*P(79)+EijR(2,1)*P(185)+(
     -   EijR(8,3)-EijR(18,3))*P(358)-EijR(2,3)*P(379)+EijR(10,3)*P(3
     -   85)+Is34*(-Cij125R(2,1)+Cij235R(1,1)+10*Dij1235R(7,2)+4*Dij1
     -   235R(12,3)-Dij1235R(6,2)*P(73)-2*(C0235R+Dij1235R(3,1)*P(137
     -   ))+Dij1235R(2,2)*P(373)+Dij1235R(4,2)*P(374)+Dij1235R(2,1)*P
     -   (408))+EijR(2,2)*P(410)-2*(D02345R-s34*EijR(10,2)-7*EijR(11,
     -   2)-8*EijR(22,3)+5*EijR(24,3)-2*(EijR(37,4)-EijR(44,4))-EijR(
     -   7,2)*P(33)-EijR(3,1)*P(156)+EijR(4,1)*P(411))+EijR(5,2)*P(41
     -   2)+EijR(8,2)*P(413)-EijR(9,2)*P(415))
       FI(46) = -4*Is15s34*(B025I+2*Cij125I(4,2)+C0125I*P(3)-Cij125I(1
     -   ,1)*P(407))+2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(5,2)+Eij
     -   I(15,3)*P(7)+(EijI(9,3)-EijI(20,3))*P(79)+EijI(2,1)*P(185)+(
     -   EijI(8,3)-EijI(18,3))*P(358)-EijI(2,3)*P(379)+EijI(10,3)*P(3
     -   85)+Is34*(-Cij125I(2,1)+Cij235I(1,1)+10*Dij1235I(7,2)+4*Dij1
     -   235I(12,3)-Dij1235I(6,2)*P(73)-2*(C0235I+Dij1235I(3,1)*P(137
     -   ))+Dij1235I(2,2)*P(373)+Dij1235I(4,2)*P(374)+Dij1235I(2,1)*P
     -   (408))+EijI(2,2)*P(410)-2*(D02345I-s34*EijI(10,2)-7*EijI(11,
     -   2)-8*EijI(22,3)+5*EijI(24,3)-2*(EijI(37,4)-EijI(44,4))-EijI(
     -   7,2)*P(33)-EijI(3,1)*P(156)+EijI(4,1)*P(411))+EijI(5,2)*P(41
     -   2)+EijI(8,2)*P(413)-EijI(9,2)*P(415))
       F(46)=DCMPLX(FR(46),FI(46))
       P(416) = p2sq+p5sq
       P(417) = p3sq+s12-4*s34-2*P(416)
       FR(47) = -4*Is15s34*(B025R+2*Cij125R(4,2)+C0125R*P(3)-Cij125R(1
     -   ,1)*P(407))+2*(Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(4,2)-(2*
     -   EijR(4,1)+EijR(9,2)+EijR(10,3)-EijR(20,3))*P(7)-2*(D02345R+s
     -   45*EijR(3,1)-s34*EijR(3,2)-7*EijR(11,2)-5*(EijR(22,3)-EijR(2
     -   3,3))-2*(EijR(37,4)-EijR(43,4))-EijR(6,2)*P(33))-EijR(12,3)*
     -   P(79)+EijR(2,1)*P(185)-(-EijR(8,3)+EijR(17,3))*P(358)-EijR(2
     -   ,3)*P(379)+EijR(9,3)*P(384)+Is34*(-Cij125R(2,1)+Cij235R(1,1)
     -   +10*Dij1235R(7,2)+4*Dij1235R(12,3)-Dij1235R(6,2)*P(73)-2*(C0
     -   235R+Dij1235R(3,1)*P(137))+Dij1235R(2,2)*P(373)+Dij1235R(4,2
     -   )*P(374)+Dij1235R(2,1)*P(408))+EijR(2,2)*P(410)+EijR(5,2)*P(
     -   412)+EijR(8,2)*P(417))
       FI(47) = -4*Is15s34*(B025I+2*Cij125I(4,2)+C0125I*P(3)-Cij125I(1
     -   ,1)*P(407))+2*(Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(4,2)-(2*
     -   EijI(4,1)+EijI(9,2)+EijI(10,3)-EijI(20,3))*P(7)-2*(D02345I+s
     -   45*EijI(3,1)-s34*EijI(3,2)-7*EijI(11,2)-5*(EijI(22,3)-EijI(2
     -   3,3))-2*(EijI(37,4)-EijI(43,4))-EijI(6,2)*P(33))-EijI(12,3)*
     -   P(79)+EijI(2,1)*P(185)-(-EijI(8,3)+EijI(17,3))*P(358)-EijI(2
     -   ,3)*P(379)+EijI(9,3)*P(384)+Is34*(-Cij125I(2,1)+Cij235I(1,1)
     -   +10*Dij1235I(7,2)+4*Dij1235I(12,3)-Dij1235I(6,2)*P(73)-2*(C0
     -   235I+Dij1235I(3,1)*P(137))+Dij1235I(2,2)*P(373)+Dij1235I(4,2
     -   )*P(374)+Dij1235I(2,1)*P(408))+EijI(2,2)*P(410)+EijI(5,2)*P(
     -   412)+EijI(8,2)*P(417))
       F(47)=DCMPLX(FR(47),FI(47))
       P(418) = -s45+2*P(342)
       P(419) = p2sq-p5sq-s15
       P(420) = p3sq+3*s12-s23-4*s34+s45+2*P(419)
       P(421) = p3sq-2*s34+3*P(16)
       P(422) = -s23+s45+2*P(330)
       P(423) = p2sq+p3sq-p5sq-s15-s45
       P(424) = -s23+6*P(6)+2*P(423)
       FR(48) = 2*(3*(Dij1235R(2,1)+Dij1235R(2,2)-Dij1235R(3,1))-7*Dij
     -   1235R(6,2)-Dij2345R(1,1)+Dij2345R(1,2)+Dij2345R(2,1)-Dij2345
     -   R(4,2)+14*EijR(22,3)-10*EijR(23,3)+4*(Dij1235R(3,2)-EijR(11,
     -   2)+EijR(37,4)-EijR(43,4))+(-3*(EijR(9,2)-EijR(10,2))-EijR(10
     -   ,3)+EijR(20,3))*P(7)-EijR(12,3)*P(79)+(EijR(8,3)-EijR(17,3))
     -   *P(358)-EijR(2,3)*P(379)+EijR(9,3)*P(384)+(EijR(2,1)-EijR(3,
     -   1))*P(418)-EijR(2,2)*P(420)-EijR(3,2)*P(421)+(EijR(5,2)-EijR
     -   (6,2))*P(422)+EijR(8,2)*P(424))
       FI(48) = 2*(3*(Dij1235I(2,1)+Dij1235I(2,2)-Dij1235I(3,1))-7*Dij
     -   1235I(6,2)-Dij2345I(1,1)+Dij2345I(1,2)+Dij2345I(2,1)-Dij2345
     -   I(4,2)+14*EijI(22,3)-10*EijI(23,3)+4*(Dij1235I(3,2)-EijI(11,
     -   2)+EijI(37,4)-EijI(43,4))+(-3*(EijI(9,2)-EijI(10,2))-EijI(10
     -   ,3)+EijI(20,3))*P(7)-EijI(12,3)*P(79)+(EijI(8,3)-EijI(17,3))
     -   *P(358)-EijI(2,3)*P(379)+EijI(9,3)*P(384)+(EijI(2,1)-EijI(3,
     -   1))*P(418)-EijI(2,2)*P(420)-EijI(3,2)*P(421)+(EijI(5,2)-EijI
     -   (6,2))*P(422)+EijI(8,2)*P(424))
       F(48)=DCMPLX(FR(48),FI(48))
       P(425) = p3sq-5*p5sq+3*s12-s23+2*P(74)-4*P(290)
       FR(49) = 2*(3*(Dij1235R(2,1)+Dij1235R(2,2)-Dij1235R(3,1))+4*Dij
     -   1235R(3,2)-7*Dij1235R(6,2)-Dij2345R(1,1)+Dij2345R(1,2)+Dij23
     -   45R(3,1)-Dij2345R(5,2)+2*(EijR(11,2)+8*EijR(22,3)-5*EijR(24,
     -   3)+2*(EijR(37,4)-EijR(44,4)))+(3*EijR(4,2)+EijR(15,3))*P(7)+
     -   (EijR(9,3)-EijR(20,3))*P(79)+(EijR(8,3)-EijR(18,3))*P(358)-E
     -   ijR(2,3)*P(379)+EijR(10,3)*P(385)+(EijR(2,1)-EijR(4,1))*P(41
     -   8)-EijR(2,2)*P(420)+(EijR(8,2)-EijR(10,2))*P(421)+(EijR(5,2)
     -   -EijR(7,2))*P(422)+EijR(9,2)*P(425))
       FI(49) = 2*(3*(Dij1235I(2,1)+Dij1235I(2,2)-Dij1235I(3,1))+4*Dij
     -   1235I(3,2)-7*Dij1235I(6,2)-Dij2345I(1,1)+Dij2345I(1,2)+Dij23
     -   45I(3,1)-Dij2345I(5,2)+2*(EijI(11,2)+8*EijI(22,3)-5*EijI(24,
     -   3)+2*(EijI(37,4)-EijI(44,4)))+(3*EijI(4,2)+EijI(15,3))*P(7)+
     -   (EijI(9,3)-EijI(20,3))*P(79)+(EijI(8,3)-EijI(18,3))*P(358)-E
     -   ijI(2,3)*P(379)+EijI(10,3)*P(385)+(EijI(2,1)-EijI(4,1))*P(41
     -   8)-EijI(2,2)*P(420)+(EijI(8,2)-EijI(10,2))*P(421)+(EijI(5,2)
     -   -EijI(7,2))*P(422)+EijI(9,2)*P(425))
       F(49)=DCMPLX(FR(49),FI(49))
       P(426) = -3*s12+s34
       P(427) = p5sq*s34-P(6)**2
       P(428) = s12*s45-2*P(427)
       P(429) = p2sq+s23-3*s34+2*P(16)
       P(430) = p2sq+s15
       P(431) = -p3sq+s23+s34+2*P(430)
       P(432) = p2sq-s23+s34
       P(433) = -p3sq+s34-3*P(7)
       P(434) = p2sq-s12+s15
       P(435) = -p3sq-s15+s23+s34
       P(436) = p2sq-s15+s34
       P(437) = 3*p2sq-s15+s34
       P(438) = 2*p5sq+s12-s34-s45
       P(439) = p2sq*P(201)+s12*P(438)
       P(440) = -3*s34+2*P(376)
       P(441) = p3sq-2*s15-s23-s34
       P(442) = p2sq**2+p2sq*P(440)+s12*P(441)
       P(443) = p5sq+s12-s15-s34
       P(444) = p3sq+2*P(443)
       P(445) = p3sq*p5sq+p5sq*s12-s12*s15
       P(446) = -(p2sq*s12)+s12**2+s12*s23+p2sq*s34-s12*s34-s23*s34+s4
     -   5**2-s45*P(444)+2*P(445)
       P(447) = p2sq-s23-s34
       P(448) = p3sq**2+p3sq*P(176)+P(16)*P(447)
       P(449) = p3sq+3*p5sq+s12-2*P(19)
       P(450) = p3sq*p5sq+p5sq**2-p5sq*s15
       P(451) = p5sq*s12-s12*s15+s12*s23-3*p5sq*s34+s15*s34-s23*s34+s4
     -   5**2-s45*P(449)+2*P(450)
       P(452) = p5sq-2*s15-s34
       P(453) = -s12-s23+s45+2*P(102)
       P(454) = p2sq**2+2*p2sq*P(452)+s12*P(453)
       P(455) = -p5sq+s12+s15
       P(456) = -s23+s45
       P(457) = s23+2*s34+s45
       P(458) = p3sq*p5sq+p5sq*s12+s12*s34-s15*P(104)+s45*P(455)
       P(459) = p2sq**2-s12**2+P(30)*P(50)+p3sq*P(456)-p2sq*P(457)+2*P
     -   (458)
       P(460) = p3sq-p5sq+s12+2*s15+s23-s34
       P(461) = -(p5sq*s12)-s12*s15-p3sq*s23-p5sq*s23+s12*s23+s23*s34-
     -   s45**2+p2sq*P(23)+2*P(219)+s45*P(460)
       P(462) = s15+s23+s34
       P(463) = -(s23*s34)+s12*P(462)
       P(464) = -s23-2*s45+3*P(6)
       P(465) = -2*s15-s23-s45+3*P(6)
       P(466) = p2sq**2+p3sq**2+s45*P(361)-2*P(463)+p2sq*P(464)+p3sq*P
     -   (465)
       P(467) = s15-s23+2*P(7)
       P(468) = p5sq+s12-2*s15-s23-3*s34-s45
       P(469) = s23*s34+s34**2
       P(470) = p3sq**2+s12*s15-s12*s23-s12*s34+s15*s34-p5sq*P(361)+s4
     -   5*P(361)+p2sq*P(467)+p3sq*P(468)+2*P(469)
       P(471) = s15-s23-s34
       P(472) = s15-2*s23
       P(473) = 2*p5sq+s12-3*P(30)
       P(474) = p3sq**2+p2sq*P(7)+2*s34*P(30)-p5sq*P(83)+s12*P(471)-s4
     -   5*P(472)+p3sq*P(473)
       P(475) = p2sq+p3sq+s12
       P(476) = s15+s23+s45
       P(477) = -9*s34+3*P(475)-2*P(476)
       FR(50) = 4*Dij2345R(7,2)-8*(Dij1235R(7,2)-EijR(46,4))+EijR(4,2)
     -   *P(2)*P(7)-Dij1345R(5,2)*P(14)-(Dij1235R(2,1)-Dij1235R(4,2)+
     -   Dij1235R(5,2))*P(33)-Dij1345R(1,2)*P(198)-D02345R*P(429)-Dij
     -   2345R(1,1)*P(431)+Dij2345R(2,1)*P(432)+Is12*(4*Cij135R(4,2)+
     -   s34*Dij1345R(6,2)*P(14)-Is45*P(6)*(2*B012R+4*Cij145R(4,2)+Ci
     -   j145R(1,1)*P(95)+C0145R*P(185))-Dij1345R(1,1)*P(323)+Dij1345
     -   R(2,1)*P(326)+Dij1345R(4,2)*P(331)+2*(B035R-p3sq*s34*Dij1345
     -   R(2,2)+C0135R*P(35)+Cij135R(1,1)*P(59)+Dij1345R(7,2)*P(333))
     -   +Cij135R(2,1)*P(426)-D01345R*P(428)+P(6)*(-5*C0345R+Dij1345R
     -   (3,1)*P(433)))+Dij1235R(3,1)*P(434)+Dij2345R(3,1)*P(435)-Dij
     -   1235R(3,2)*P(436)+Dij1235R(6,2)*P(437)-EijR(2,1)*P(439)+EijR
     -   (2,2)*P(442)+EijR(3,1)*P(446)+EijR(3,2)*P(448)-EijR(4,1)*P(4
     -   51)-EijR(5,2)*P(454)+EijR(6,2)*P(459)-EijR(7,2)*P(461)-EijR(
     -   8,2)*P(466)+EijR(9,2)*P(470)-EijR(10,2)*P(474)+2*(C0235R-p2s
     -   q*Dij1235R(2,2)-EijR(24,3)*P(66)+P(3)*(D01235R+Dij1235R(1,1)
     -   -EE0R*P(6)+EijR(1,2)*P(33)+EijR(1,1)*P(71))-EijR(22,3)*P(356
     -   )+EijR(21,3)*P(358)+EijR(23,3)*P(365)+EijR(11,2)*P(477))
       FI(50) = 4*Dij2345I(7,2)-8*(Dij1235I(7,2)-EijI(46,4))+EijI(4,2)
     -   *P(2)*P(7)-Dij1345I(5,2)*P(14)-(Dij1235I(2,1)-Dij1235I(4,2)+
     -   Dij1235I(5,2))*P(33)-Dij1345I(1,2)*P(198)-D02345I*P(429)-Dij
     -   2345I(1,1)*P(431)+Dij2345I(2,1)*P(432)+Is12*(4*Cij135I(4,2)+
     -   s34*Dij1345I(6,2)*P(14)-Is45*P(6)*(2*B012I+4*Cij145I(4,2)+Ci
     -   j145I(1,1)*P(95)+C0145I*P(185))-Dij1345I(1,1)*P(323)+Dij1345
     -   I(2,1)*P(326)+Dij1345I(4,2)*P(331)+2*(B035I-p3sq*s34*Dij1345
     -   I(2,2)+C0135I*P(35)+Cij135I(1,1)*P(59)+Dij1345I(7,2)*P(333))
     -   +Cij135I(2,1)*P(426)-D01345I*P(428)+P(6)*(-5*C0345I+Dij1345I
     -   (3,1)*P(433)))+Dij1235I(3,1)*P(434)+Dij2345I(3,1)*P(435)-Dij
     -   1235I(3,2)*P(436)+Dij1235I(6,2)*P(437)-EijI(2,1)*P(439)+EijI
     -   (2,2)*P(442)+EijI(3,1)*P(446)+EijI(3,2)*P(448)-EijI(4,1)*P(4
     -   51)-EijI(5,2)*P(454)+EijI(6,2)*P(459)-EijI(7,2)*P(461)-EijI(
     -   8,2)*P(466)+EijI(9,2)*P(470)-EijI(10,2)*P(474)+2*(C0235I-p2s
     -   q*Dij1235I(2,2)-EijI(24,3)*P(66)+P(3)*(D01235I+Dij1235I(1,1)
     -   -EE0I*P(6)+EijI(1,2)*P(33)+EijI(1,1)*P(71))-EijI(22,3)*P(356
     -   )+EijI(21,3)*P(358)+EijI(23,3)*P(365)+EijI(11,2)*P(477))
       F(50)=DCMPLX(FR(50),FI(50))
       P(478) = p5sq+s15
       P(479) = p3sq+3*s12-s23-4*s34+s45-2*P(478)
       P(480) = p3sq-4*s34+3*P(16)
       FR(51) = 2*(3*(Dij1235R(2,1)+Dij1235R(2,2))-Dij2345R(1,1)+Dij23
     -   45R(1,2)+10*EijR(11,2)+20*EijR(22,3)-4*(Dij1235R(6,2)-EijR(3
     -   7,4))-(2*EijR(4,1)+3*EijR(9,2)+EijR(10,3))*P(7)+EijR(9,3)*P(
     -   79)-2*(D02345R+Dij1235R(3,1)-EijR(3,1)*P(156))+(EijR(5,2)+Ei
     -   jR(8,3))*P(358)-EijR(2,3)*P(379)+EijR(2,1)*P(418)-EijR(2,2)*
     -   P(479)+EijR(8,2)*P(480))
       FI(51) = 2*(3*(Dij1235I(2,1)+Dij1235I(2,2))-Dij2345I(1,1)+Dij23
     -   45I(1,2)+10*EijI(11,2)+20*EijI(22,3)-4*(Dij1235I(6,2)-EijI(3
     -   7,4))-(2*EijI(4,1)+3*EijI(9,2)+EijI(10,3))*P(7)+EijI(9,3)*P(
     -   79)-2*(D02345I+Dij1235I(3,1)-EijI(3,1)*P(156))+(EijI(5,2)+Ei
     -   jI(8,3))*P(358)-EijI(2,3)*P(379)+EijI(2,1)*P(418)-EijI(2,2)*
     -   P(479)+EijI(8,2)*P(480))
       F(51)=DCMPLX(FR(51),FI(51))
       FR(52) = -8*(EijR(2,3)-EijR(2,4)+EijR(5,3)+EijR(5,4)-EijR(6,3)-
     -   EijR(6,4)-EijR(9,3)+EijR(9,4)-2*(EijR(8,3)-EijR(17,3))+3*(Ei
     -   jR(8,4)-EijR(17,4)+EijR(23,4)-EijR(26,4)))
       FI(52) = -8*(EijI(2,3)-EijI(2,4)+EijI(5,3)+EijI(5,4)-EijI(6,3)-
     -   EijI(6,4)-EijI(9,3)+EijI(9,4)-2*(EijI(8,3)-EijI(17,3))+3*(Ei
     -   jI(8,4)-EijI(17,4)+EijI(23,4)-EijI(26,4)))
       F(52)=DCMPLX(FR(52),FI(52))
       FR(53) = 8*(EijR(2,4)+EijR(17,4)+EijR(18,4)+EijR(20,4)-2*(EijR(
     -   8,4)+EijR(9,4)+EijR(23,4)-2*EijR(26,4)+EijR(29,4)))
       FI(53) = 8*(EijI(2,4)+EijI(17,4)+EijI(18,4)+EijI(20,4)-2*(EijI(
     -   8,4)+EijI(9,4)+EijI(23,4)-2*EijI(26,4)+EijI(29,4)))
       F(53)=DCMPLX(FR(53),FI(53))
       FR(54) = 8*(EijR(2,4)-EijR(6,3)+EijR(7,3)-EijR(9,3)-EijR(9,4)+E
     -   ijR(10,3)-EijR(10,4)+EijR(17,4)-EijR(23,4)-EijR(24,4)+EijR(2
     -   5,4)+EijR(28,4)-2*(EijR(8,4)-EijR(17,3)+EijR(18,3)-EijR(26,4
     -   )-EijR(27,4)+EijR(35,4)))
       FI(54) = 8*(EijI(2,4)-EijI(6,3)+EijI(7,3)-EijI(9,3)-EijI(9,4)+E
     -   ijI(10,3)-EijI(10,4)+EijI(17,4)-EijI(23,4)-EijI(24,4)+EijI(2
     -   5,4)+EijI(28,4)-2*(EijI(8,4)-EijI(17,3)+EijI(18,3)-EijI(26,4
     -   )-EijI(27,4)+EijI(35,4)))
       F(54)=DCMPLX(FR(54),FI(54))
       FR(55) = 8*(EijR(2,4)-EijR(6,3)-EijR(9,3)-EijR(9,4)+EijR(17,4)-
     -   EijR(23,4)-2*(EijR(8,4)-EijR(17,3)-EijR(26,4)))
       FI(55) = 8*(EijI(2,4)-EijI(6,3)-EijI(9,3)-EijI(9,4)+EijI(17,4)-
     -   EijI(23,4)-2*(EijI(8,4)-EijI(17,3)-EijI(26,4)))
       F(55)=DCMPLX(FR(55),FI(55))
       FR(56) = -8*(EijR(2,2)-EijR(2,4)-EijR(5,2)-EijR(5,3)+EijR(6,2)+
     -   EijR(6,3)-EijR(8,2)+EijR(8,3)-EijR(9,3)-EijR(11,3)+EijR(12,3
     -   )-EijR(17,4)-EijR(18,4)-EijR(20,4)-4*EijR(26,4)+2*(EijR(8,4)
     -   +EijR(9,4)+EijR(23,4)+EijR(29,4)))
       FI(56) = -8*(EijI(2,2)-EijI(2,4)-EijI(5,2)-EijI(5,3)+EijI(6,2)+
     -   EijI(6,3)-EijI(8,2)+EijI(8,3)-EijI(9,3)-EijI(11,3)+EijI(12,3
     -   )-EijI(17,4)-EijI(18,4)-EijI(20,4)-4*EijI(26,4)+2*(EijI(8,4)
     -   +EijI(9,4)+EijI(23,4)+EijI(29,4)))
       F(56)=DCMPLX(FR(56),FI(56))
       FR(57) = 8*(EijR(2,3)+EijR(2,4)-EijR(8,3)-EijR(8,4)-EijR(11,3)+
     -   EijR(11,4)+EijR(12,3)-EijR(12,4)-2*(EijR(9,3)-EijR(17,3))-3*
     -   (EijR(9,4)-EijR(20,4)-EijR(26,4)+EijR(29,4)))
       FI(57) = 8*(EijI(2,3)+EijI(2,4)-EijI(8,3)-EijI(8,4)-EijI(11,3)+
     -   EijI(11,4)+EijI(12,3)-EijI(12,4)-2*(EijI(9,3)-EijI(17,3))-3*
     -   (EijI(9,4)-EijI(20,4)-EijI(26,4)+EijI(29,4)))
       F(57)=DCMPLX(FR(57),FI(57))
       FR(58) = 8*(EijR(2,3)+EijR(2,4)-EijR(8,3)-EijR(8,4)-EijR(10,4)-
     -   EijR(11,3)+EijR(12,3)+EijR(20,4)+EijR(27,4)-EijR(29,4)+EijR(
     -   30,4)-EijR(31,4)-2*(EijR(9,3)+EijR(9,4)-EijR(17,3)-EijR(26,4
     -   )-EijR(28,4)+EijR(35,4)))
       FI(58) = 8*(EijI(2,3)+EijI(2,4)-EijI(8,3)-EijI(8,4)-EijI(10,4)-
     -   EijI(11,3)+EijI(12,3)+EijI(20,4)+EijI(27,4)-EijI(29,4)+EijI(
     -   30,4)-EijI(31,4)-2*(EijI(9,3)+EijI(9,4)-EijI(17,3)-EijI(26,4
     -   )-EijI(28,4)+EijI(35,4)))
       F(58)=DCMPLX(FR(58),FI(58))
       FR(59) = 8*(EijR(2,3)+EijR(2,4)-EijR(8,3)-EijR(8,4)-EijR(11,3)+
     -   EijR(12,3)+EijR(20,4)-2*(EijR(9,3)+EijR(9,4)-EijR(17,3)-EijR
     -   (26,4))-EijR(29,4))
       FI(59) = 8*(EijI(2,3)+EijI(2,4)-EijI(8,3)-EijI(8,4)-EijI(11,3)+
     -   EijI(12,3)+EijI(20,4)-2*(EijI(9,3)+EijI(9,4)-EijI(17,3)-EijI
     -   (26,4))-EijI(29,4))
       F(59)=DCMPLX(FR(59),FI(59))
       FR(60) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)-EijR(5,2)+EijR(6,2)+
     -   EijR(6,3)+EijR(7,3)-EijR(8,2)+3*EijR(8,3)+EijR(9,4)+EijR(10,
     -   4)-EijR(17,3)-EijR(17,4)-EijR(18,3)-EijR(19,3)+EijR(20,3)+Ei
     -   jR(23,4)+EijR(24,4)-EijR(25,4)-EijR(28,4)-2*(EijR(5,3)-EijR(
     -   8,4)+EijR(26,4)+EijR(27,4)-EijR(35,4)))
       FI(60) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)-EijI(5,2)+EijI(6,2)+
     -   EijI(6,3)+EijI(7,3)-EijI(8,2)+3*EijI(8,3)+EijI(9,4)+EijI(10,
     -   4)-EijI(17,3)-EijI(17,4)-EijI(18,3)-EijI(19,3)+EijI(20,3)+Ei
     -   jI(23,4)+EijI(24,4)-EijI(25,4)-EijI(28,4)-2*(EijI(5,3)-EijI(
     -   8,4)+EijI(26,4)+EijI(27,4)-EijI(35,4)))
       F(60)=DCMPLX(FR(60),FI(60))
       FR(61) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-E
     -   ijR(8,4)-EijR(10,3)-EijR(10,4)-EijR(11,3)+EijR(12,3)-3*(EijR
     -   (9,3)-EijR(17,3))+EijR(18,3)-EijR(19,3)+EijR(20,3)+EijR(20,4
     -   )+EijR(27,4)-EijR(29,4)+EijR(30,4)-EijR(31,4)+2*(EijR(2,3)-E
     -   ijR(8,3)-EijR(9,4)+EijR(26,4)+EijR(28,4)-EijR(35,4)))
       FI(61) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-E
     -   ijI(8,4)-EijI(10,3)-EijI(10,4)-EijI(11,3)+EijI(12,3)-3*(EijI
     -   (9,3)-EijI(17,3))+EijI(18,3)-EijI(19,3)+EijI(20,3)+EijI(20,4
     -   )+EijI(27,4)-EijI(29,4)+EijI(30,4)-EijI(31,4)+2*(EijI(2,3)-E
     -   ijI(8,3)-EijI(9,4)+EijI(26,4)+EijI(28,4)-EijI(35,4)))
       F(61)=DCMPLX(FR(61),FI(61))
       FR(62) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-E
     -   ijR(8,4)-EijR(9,4)+EijR(21,4)+EijR(26,4)-EijR(32,4)+EijR(33,
     -   4)-EijR(34,4)+2*(EijR(2,3)-EijR(8,3)-EijR(9,3)-EijR(10,3)-Ei
     -   jR(10,4)+EijR(17,3)+EijR(18,3)-EijR(19,3)+EijR(20,3)+EijR(27
     -   ,4)+EijR(28,4)-EijR(35,4)))
       FI(62) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-E
     -   ijI(8,4)-EijI(9,4)+EijI(21,4)+EijI(26,4)-EijI(32,4)+EijI(33,
     -   4)-EijI(34,4)+2*(EijI(2,3)-EijI(8,3)-EijI(9,3)-EijI(10,3)-Ei
     -   jI(10,4)+EijI(17,3)+EijI(18,3)-EijI(19,3)+EijI(20,3)+EijI(27
     -   ,4)+EijI(28,4)-EijI(35,4)))
       F(62)=DCMPLX(FR(62),FI(62))
       FR(63) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(7,2)-E
     -   ijR(8,2)-EijR(8,4)+EijR(9,2)-EijR(9,4)-EijR(10,3)-EijR(10,4)
     -   +2*(EijR(2,3)-EijR(8,3)-EijR(9,3)+EijR(17,3))+EijR(18,3)-Eij
     -   R(19,3)+EijR(20,3)+EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,
     -   4))
       FI(63) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(7,2)-E
     -   ijI(8,2)-EijI(8,4)+EijI(9,2)-EijI(9,4)-EijI(10,3)-EijI(10,4)
     -   +2*(EijI(2,3)-EijI(8,3)-EijI(9,3)+EijI(17,3))+EijI(18,3)-Eij
     -   I(19,3)+EijI(20,3)+EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,
     -   4))
       F(63)=DCMPLX(FR(63),FI(63))
       FR(64) = -8*Is34*(-Dij1235R(2,2)+Dij1235R(2,3)+Dij1235R(4,2)+Di
     -   j1235R(4,3)-2*Dij1235R(6,3)+s34*(EijR(2,2)-EijR(2,3)-EijR(2,
     -   4)-EijR(5,2)+EijR(6,2)+EijR(6,3)-EijR(8,2)+3*EijR(8,3)+EijR(
     -   9,4)-EijR(17,3)-EijR(17,4)+EijR(23,4)-2*(EijR(5,3)-EijR(8,4)
     -   +EijR(26,4))))
       FI(64) = -8*Is34*(-Dij1235I(2,2)+Dij1235I(2,3)+Dij1235I(4,2)+Di
     -   j1235I(4,3)-2*Dij1235I(6,3)+s34*(EijI(2,2)-EijI(2,3)-EijI(2,
     -   4)-EijI(5,2)+EijI(6,2)+EijI(6,3)-EijI(8,2)+3*EijI(8,3)+EijI(
     -   9,4)-EijI(17,3)-EijI(17,4)+EijI(23,4)-2*(EijI(5,3)-EijI(8,4)
     -   +EijI(26,4))))
       F(64)=DCMPLX(FR(64),FI(64))
       FR(65) = 8*Is34*(-Dij1235R(2,3)+Dij1235R(5,2)-Dij1235R(6,2)+Dij
     -   1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)+s34*(EijR(2,2)+EijR(
     -   2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-EijR(8,4)-EijR(11,3)+EijR
     -   (12,3)-3*(EijR(9,3)-EijR(17,3))+EijR(20,4)+2*(EijR(2,3)-EijR
     -   (8,3)-EijR(9,4)+EijR(26,4))-EijR(29,4)))
       FI(65) = 8*Is34*(-Dij1235I(2,3)+Dij1235I(5,2)-Dij1235I(6,2)+Dij
     -   1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)+s34*(EijI(2,2)+EijI(
     -   2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-EijI(8,4)-EijI(11,3)+EijI
     -   (12,3)-3*(EijI(9,3)-EijI(17,3))+EijI(20,4)+2*(EijI(2,3)-EijI
     -   (8,3)-EijI(9,4)+EijI(26,4))-EijI(29,4)))
       F(65)=DCMPLX(FR(65),FI(65))
       FR(66) = 8*Is34*(-Dij1235R(2,3)+Dij1235R(5,2)-Dij1235R(6,2)+Dij
     -   1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)+s34*(EijR(2,2)+EijR(
     -   2,4)-EijR(5,2)+EijR(6,2)+EijR(7,2)-EijR(8,2)-EijR(8,4)-EijR(
     -   9,2)-EijR(9,4)-EijR(10,3)-EijR(10,4)+2*(EijR(2,3)-EijR(8,3)-
     -   EijR(9,3)+EijR(17,3))+EijR(18,3)-EijR(19,3)+EijR(20,3)+EijR(
     -   26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4)))
       FI(66) = 8*Is34*(-Dij1235I(2,3)+Dij1235I(5,2)-Dij1235I(6,2)+Dij
     -   1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)+s34*(EijI(2,2)+EijI(
     -   2,4)-EijI(5,2)+EijI(6,2)+EijI(7,2)-EijI(8,2)-EijI(8,4)-EijI(
     -   9,2)-EijI(9,4)-EijI(10,3)-EijI(10,4)+2*(EijI(2,3)-EijI(8,3)-
     -   EijI(9,3)+EijI(17,3))+EijI(18,3)-EijI(19,3)+EijI(20,3)+EijI(
     -   26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4)))
       F(66)=DCMPLX(FR(66),FI(66))
       FR(67) = 8*Is34*(-Dij1235R(2,3)+Dij1235R(6,3)+s34*(EijR(2,2)+Ei
     -   jR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-EijR(8,4)-EijR(9,4)+2*
     -   (EijR(2,3)-EijR(8,3)-EijR(9,3)+EijR(17,3))+EijR(26,4)))
       FI(67) = 8*Is34*(-Dij1235I(2,3)+Dij1235I(6,3)+s34*(EijI(2,2)+Ei
     -   jI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-EijI(8,4)-EijI(9,4)+2*
     -   (EijI(2,3)-EijI(8,3)-EijI(9,3)+EijI(17,3))+EijI(26,4)))
       F(67)=DCMPLX(FR(67),FI(67))
       FR(68) = -8*(EijR(2,3)-EijR(2,4)-EijR(8,3)-EijR(11,3)+EijR(12,3
     -   )-EijR(17,4)-EijR(18,4)-EijR(20,4)-4*EijR(26,4)+2*(EijR(8,4)
     -   -EijR(9,3)+EijR(9,4)+EijR(17,3)+EijR(23,4)+EijR(29,4)))
       FI(68) = -8*(EijI(2,3)-EijI(2,4)-EijI(8,3)-EijI(11,3)+EijI(12,3
     -   )-EijI(17,4)-EijI(18,4)-EijI(20,4)-4*EijI(26,4)+2*(EijI(8,4)
     -   -EijI(9,3)+EijI(9,4)+EijI(17,3)+EijI(23,4)+EijI(29,4)))
       F(68)=DCMPLX(FR(68),FI(68))
       FR(69) = 8*(EijR(2,4)-EijR(8,4)+EijR(11,4)-EijR(12,4)-3*(EijR(9
     -   ,4)-EijR(20,4)-EijR(26,4)+EijR(29,4)))
       FI(69) = 8*(EijI(2,4)-EijI(8,4)+EijI(11,4)-EijI(12,4)-3*(EijI(9
     -   ,4)-EijI(20,4)-EijI(26,4)+EijI(29,4)))
       F(69)=DCMPLX(FR(69),FI(69))
       FR(70) = 8*(EijR(2,4)-EijR(8,4)-EijR(9,3)+EijR(10,3)-EijR(10,4)
     -   -EijR(11,3)+EijR(12,3)+EijR(17,3)-EijR(18,3)+EijR(19,3)-EijR
     -   (20,3)+EijR(20,4)+EijR(27,4)-EijR(29,4)+EijR(30,4)-EijR(31,4
     -   )-2*(EijR(9,4)-EijR(26,4)-EijR(28,4)+EijR(35,4)))
       FI(70) = 8*(EijI(2,4)-EijI(8,4)-EijI(9,3)+EijI(10,3)-EijI(10,4)
     -   -EijI(11,3)+EijI(12,3)+EijI(17,3)-EijI(18,3)+EijI(19,3)-EijI
     -   (20,3)+EijI(20,4)+EijI(27,4)-EijI(29,4)+EijI(30,4)-EijI(31,4
     -   )-2*(EijI(9,4)-EijI(26,4)-EijI(28,4)+EijI(35,4)))
       F(70)=DCMPLX(FR(70),FI(70))
       FR(71) = 8*(EijR(2,4)-EijR(8,4)-EijR(9,3)-EijR(11,3)+EijR(12,3)
     -   +EijR(17,3)+EijR(20,4)-2*(EijR(9,4)-EijR(26,4))-EijR(29,4))
       FI(71) = 8*(EijI(2,4)-EijI(8,4)-EijI(9,3)-EijI(11,3)+EijI(12,3)
     -   +EijI(17,3)+EijI(20,4)-2*(EijI(9,4)-EijI(26,4))-EijI(29,4))
       F(71)=DCMPLX(FR(71),FI(71))
       FR(72) = -8*(EijR(2,2)-EijR(2,4)+EijR(3,2)-EijR(3,3)+EijR(8,3)+
     -   EijR(8,4)-EijR(9,3)+EijR(11,3)-EijR(11,4)+EijR(12,4)-2*(EijR
     -   (8,2)-EijR(12,3)+EijR(17,3))+3*(EijR(9,4)-EijR(20,4)-EijR(26
     -   ,4)+EijR(29,4)))
       FI(72) = -8*(EijI(2,2)-EijI(2,4)+EijI(3,2)-EijI(3,3)+EijI(8,3)+
     -   EijI(8,4)-EijI(9,3)+EijI(11,3)-EijI(11,4)+EijI(12,4)-2*(EijI
     -   (8,2)-EijI(12,3)+EijI(17,3))+3*(EijI(9,4)-EijI(20,4)-EijI(26
     -   ,4)+EijI(29,4)))
       F(72)=DCMPLX(FR(72),FI(72))
       FR(73) = 8*(EijR(2,3)+EijR(2,4)-EijR(3,3)+EijR(3,4)-3*(EijR(9,3
     -   )-EijR(12,3))-4*(EijR(9,4)+EijR(12,4))+6*EijR(20,4))
       FI(73) = 8*(EijI(2,3)+EijI(2,4)-EijI(3,3)+EijI(3,4)-3*(EijI(9,3
     -   )-EijI(12,3))-4*(EijI(9,4)+EijI(12,4))+6*EijI(20,4))
       F(73)=DCMPLX(FR(73),FI(73))
       FR(74) = 8*(EijR(2,3)+EijR(2,4)-EijR(3,3)-EijR(10,4)-EijR(12,4)
     -   +EijR(13,4)-3*(EijR(9,3)+EijR(9,4)-EijR(12,3)-EijR(20,4)-Eij
     -   R(28,4)+EijR(31,4)))
       FI(74) = 8*(EijI(2,3)+EijI(2,4)-EijI(3,3)-EijI(10,4)-EijI(12,4)
     -   +EijI(13,4)-3*(EijI(9,3)+EijI(9,4)-EijI(12,3)-EijI(20,4)-Eij
     -   I(28,4)+EijI(31,4)))
       F(74)=DCMPLX(FR(74),FI(74))
       FR(75) = 8*(EijR(2,3)+EijR(2,4)-EijR(3,3)-EijR(12,4)-3*(EijR(9,
     -   3)+EijR(9,4)-EijR(12,3)-EijR(20,4)))
       FI(75) = 8*(EijI(2,3)+EijI(2,4)-EijI(3,3)-EijI(12,4)-3*(EijI(9,
     -   3)+EijI(9,4)-EijI(12,3)-EijI(20,4)))
       F(75)=DCMPLX(FR(75),FI(75))
       FR(76) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)+EijR(3,2)+EijR(8,4)+
     -   EijR(9,3)+EijR(10,4)+EijR(11,3)-EijR(13,3)-3*EijR(17,3)-EijR
     -   (18,3)+EijR(19,3)+EijR(20,3)-EijR(20,4)-EijR(27,4)+EijR(29,4
     -   )-EijR(30,4)+EijR(31,4)-2*(EijR(8,2)-EijR(8,3)-EijR(9,4)+Eij
     -   R(26,4)+EijR(28,4)-EijR(35,4)))
       FI(76) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)+EijI(3,2)+EijI(8,4)+
     -   EijI(9,3)+EijI(10,4)+EijI(11,3)-EijI(13,3)-3*EijI(17,3)-EijI
     -   (18,3)+EijI(19,3)+EijI(20,3)-EijI(20,4)-EijI(27,4)+EijI(29,4
     -   )-EijI(30,4)+EijI(31,4)-2*(EijI(8,2)-EijI(8,3)-EijI(9,4)+Eij
     -   I(26,4)+EijI(28,4)-EijI(35,4)))
       F(76)=DCMPLX(FR(76),FI(76))
       FR(77) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,2)-EijR(3,3)-5*EijR(9,3)
     -   -EijR(10,3)-EijR(10,4)+4*EijR(12,3)-EijR(12,4)-EijR(13,3)+Ei
     -   jR(13,4)+2*(EijR(2,3)-EijR(8,2)+EijR(20,3))-3*(EijR(9,4)-Eij
     -   R(20,4)-EijR(28,4)+EijR(31,4)))
       FI(77) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,2)-EijI(3,3)-5*EijI(9,3)
     -   -EijI(10,3)-EijI(10,4)+4*EijI(12,3)-EijI(12,4)-EijI(13,3)+Ei
     -   jI(13,4)+2*(EijI(2,3)-EijI(8,2)+EijI(20,3))-3*(EijI(9,4)-Eij
     -   I(20,4)-EijI(28,4)+EijI(31,4)))
       F(77)=DCMPLX(FR(77),FI(77))
       FR(78) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,2)+EijR(20,4)+EijR(21,4)
     -   +EijR(22,4)-4*(EijR(9,3)-EijR(20,3)-EijR(28,4))+2*(EijR(2,3)
     -   -EijR(8,2)-EijR(9,4)-EijR(10,3)-EijR(10,4)+EijR(12,3)-EijR(1
     -   3,3)-EijR(31,4)-EijR(34,4)))
       FI(78) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,2)+EijI(20,4)+EijI(21,4)
     -   +EijI(22,4)-4*(EijI(9,3)-EijI(20,3)-EijI(28,4))+2*(EijI(2,3)
     -   -EijI(8,2)-EijI(9,4)-EijI(10,3)-EijI(10,4)+EijI(12,3)-EijI(1
     -   3,3)-EijI(31,4)-EijI(34,4)))
       F(78)=DCMPLX(FR(78),FI(78))
       FR(79) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,2)+EijR(9,2)-4*EijR(9,3)
     -   -EijR(10,2)-EijR(10,3)-EijR(10,4)-EijR(13,3)+EijR(20,4)+2*(E
     -   ijR(2,3)-EijR(8,2)-EijR(9,4)+EijR(12,3)+EijR(20,3)+EijR(28,4
     -   ))-EijR(31,4))
       FI(79) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,2)+EijI(9,2)-4*EijI(9,3)
     -   -EijI(10,2)-EijI(10,3)-EijI(10,4)-EijI(13,3)+EijI(20,4)+2*(E
     -   ijI(2,3)-EijI(8,2)-EijI(9,4)+EijI(12,3)+EijI(20,3)+EijI(28,4
     -   ))-EijI(31,4))
       F(79)=DCMPLX(FR(79),FI(79))
       FR(80) = -8*Is34*(-Dij1235R(2,1)+Dij1235R(2,2)+Dij1235R(2,3)+Di
     -   j1235R(3,1)-2*(Dij1235R(4,2)-Dij1235R(5,2))-Dij1235R(6,2)-Di
     -   j1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3)+s34*(EijR(2,2)-EijR
     -   (2,3)-EijR(2,4)+EijR(3,2)+EijR(8,4)+EijR(9,3)+EijR(11,3)-3*E
     -   ijR(17,3)-EijR(20,4)-2*(EijR(8,2)-EijR(8,3)-EijR(9,4)+EijR(2
     -   6,4))+EijR(29,4)))
       FI(80) = -8*Is34*(-Dij1235I(2,1)+Dij1235I(2,2)+Dij1235I(2,3)+Di
     -   j1235I(3,1)-2*(Dij1235I(4,2)-Dij1235I(5,2))-Dij1235I(6,2)-Di
     -   j1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3)+s34*(EijI(2,2)-EijI
     -   (2,3)-EijI(2,4)+EijI(3,2)+EijI(8,4)+EijI(9,3)+EijI(11,3)-3*E
     -   ijI(17,3)-EijI(20,4)-2*(EijI(8,2)-EijI(8,3)-EijI(9,4)+EijI(2
     -   6,4))+EijI(29,4)))
       F(80)=DCMPLX(FR(80),FI(80))
       FR(81) = -8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)-s34*(EijR(2,2)+EijR(2,4)+EijR(3,2)-EijR(3,3)+2*
     -   (EijR(2,3)-EijR(8,2))-5*EijR(9,3)+4*EijR(12,3)-EijR(12,4)-3*
     -   (EijR(9,4)-EijR(20,4))))
       FI(81) = -8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)-s34*(EijI(2,2)+EijI(2,4)+EijI(3,2)-EijI(3,3)+2*
     -   (EijI(2,3)-EijI(8,2))-5*EijI(9,3)+4*EijI(12,3)-EijI(12,4)-3*
     -   (EijI(9,4)-EijI(20,4))))
       F(81)=DCMPLX(FR(81),FI(81))
       FR(82) = -8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(-EijR(2,2)-EijR(2,4)-EijR(3,2)+EijR(9,2)+4
     -   *EijR(9,3)-EijR(10,2)+EijR(10,3)+EijR(10,4)+EijR(13,3)-EijR(
     -   20,4)-2*(EijR(2,3)-EijR(8,2)-EijR(9,4)+EijR(12,3)+EijR(20,3)
     -   +EijR(28,4))+EijR(31,4)))
       FI(82) = -8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(-EijI(2,2)-EijI(2,4)-EijI(3,2)+EijI(9,2)+4
     -   *EijI(9,3)-EijI(10,2)+EijI(10,3)+EijI(10,4)+EijI(13,3)-EijI(
     -   20,4)-2*(EijI(2,3)-EijI(8,2)-EijI(9,4)+EijI(12,3)+EijI(20,3)
     -   +EijI(28,4))+EijI(31,4)))
       F(82)=DCMPLX(FR(82),FI(82))
       FR(83) = 8*Is34*(-Dij1235R(2,1)-Dij1235R(2,3)+Dij1235R(3,1)-2*(
     -   Dij1235R(2,2)-Dij1235R(6,2))+Dij1235R(8,3)+s34*(EijR(2,2)+Ei
     -   jR(2,4)+EijR(3,2)+2*(EijR(2,3)-EijR(8,2)-2*EijR(9,3)-EijR(9,
     -   4)+EijR(12,3))+EijR(20,4)))
       FI(83) = 8*Is34*(-Dij1235I(2,1)-Dij1235I(2,3)+Dij1235I(3,1)-2*(
     -   Dij1235I(2,2)-Dij1235I(6,2))+Dij1235I(8,3)+s34*(EijI(2,2)+Ei
     -   jI(2,4)+EijI(3,2)+2*(EijI(2,3)-EijI(8,2)-2*EijI(9,3)-EijI(9,
     -   4)+EijI(12,3))+EijI(20,4)))
       F(83)=DCMPLX(FR(83),FI(83))
       FR(84) = -8*(EijR(2,3)-EijR(2,4)-EijR(8,3)-EijR(9,3)+EijR(9,4)-
     -   EijR(10,3)+EijR(10,4)+EijR(17,3)-EijR(17,4)+EijR(18,3)-EijR(
     -   19,3)+EijR(20,3)+EijR(23,4)+EijR(24,4)-EijR(25,4)-EijR(28,4)
     -   +2*(EijR(8,4)-EijR(26,4)-EijR(27,4)+EijR(35,4)))
       FI(84) = -8*(EijI(2,3)-EijI(2,4)-EijI(8,3)-EijI(9,3)+EijI(9,4)-
     -   EijI(10,3)+EijI(10,4)+EijI(17,3)-EijI(17,4)+EijI(18,3)-EijI(
     -   19,3)+EijI(20,3)+EijI(23,4)+EijI(24,4)-EijI(25,4)-EijI(28,4)
     -   +2*(EijI(8,4)-EijI(26,4)-EijI(27,4)+EijI(35,4)))
       F(84)=DCMPLX(FR(84),FI(84))
       FR(85) = 8*(EijR(2,4)-EijR(8,4)-EijR(10,4)+EijR(20,4)+EijR(27,4
     -   )-EijR(29,4)+EijR(30,4)-EijR(31,4)-2*(EijR(9,4)-EijR(26,4)-E
     -   ijR(28,4)+EijR(35,4)))
       FI(85) = 8*(EijI(2,4)-EijI(8,4)-EijI(10,4)+EijI(20,4)+EijI(27,4
     -   )-EijI(29,4)+EijI(30,4)-EijI(31,4)-2*(EijI(9,4)-EijI(26,4)-E
     -   ijI(28,4)+EijI(35,4)))
       F(85)=DCMPLX(FR(85),FI(85))
       FR(86) = 8*(EijR(2,4)-EijR(8,4)-EijR(9,3)-EijR(9,4)+EijR(10,3)+
     -   EijR(14,3)-EijR(15,3)+EijR(17,3)-EijR(18,3)-EijR(19,3)+EijR(
     -   20,3)+EijR(21,4)+EijR(26,4)-EijR(32,4)+EijR(33,4)-EijR(34,4)
     -   -2*(EijR(10,4)-EijR(27,4)-EijR(28,4)+EijR(35,4)))
       FI(86) = 8*(EijI(2,4)-EijI(8,4)-EijI(9,3)-EijI(9,4)+EijI(10,3)+
     -   EijI(14,3)-EijI(15,3)+EijI(17,3)-EijI(18,3)-EijI(19,3)+EijI(
     -   20,3)+EijI(21,4)+EijI(26,4)-EijI(32,4)+EijI(33,4)-EijI(34,4)
     -   -2*(EijI(10,4)-EijI(27,4)-EijI(28,4)+EijI(35,4)))
       F(86)=DCMPLX(FR(86),FI(86))
       FR(87) = 8*(EijR(2,4)-EijR(6,2)+EijR(7,2)+EijR(8,2)-EijR(8,4)-E
     -   ijR(9,2)-EijR(9,3)-EijR(9,4)-EijR(10,4)+EijR(17,3)-EijR(19,3
     -   )+EijR(20,3)+EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4))
       FI(87) = 8*(EijI(2,4)-EijI(6,2)+EijI(7,2)+EijI(8,2)-EijI(8,4)-E
     -   ijI(9,2)-EijI(9,3)-EijI(9,4)-EijI(10,4)+EijI(17,3)-EijI(19,3
     -   )+EijI(20,3)+EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4))
       F(87)=DCMPLX(FR(87),FI(87))
       FR(88) = -8*(EijR(2,2)-EijR(2,4)-EijR(8,2)+EijR(8,3)+EijR(8,4)-
     -   EijR(9,2)-EijR(9,3)+EijR(10,2)+EijR(10,4)+EijR(12,3)-EijR(13
     -   ,3)-EijR(17,3)-EijR(18,3)+EijR(19,3)+EijR(20,3)-EijR(20,4)-E
     -   ijR(27,4)+EijR(29,4)-EijR(30,4)+EijR(31,4)+2*(EijR(9,4)-EijR
     -   (26,4)-EijR(28,4)+EijR(35,4)))
       FI(88) = -8*(EijI(2,2)-EijI(2,4)-EijI(8,2)+EijI(8,3)+EijI(8,4)-
     -   EijI(9,2)-EijI(9,3)+EijI(10,2)+EijI(10,4)+EijI(12,3)-EijI(13
     -   ,3)-EijI(17,3)-EijI(18,3)+EijI(19,3)+EijI(20,3)-EijI(20,4)-E
     -   ijI(27,4)+EijI(29,4)-EijI(30,4)+EijI(31,4)+2*(EijI(9,4)-EijI
     -   (26,4)-EijI(28,4)+EijI(35,4)))
       F(88)=DCMPLX(FR(88),FI(88))
       FR(89) = 8*(EijR(2,3)+EijR(2,4)-EijR(10,3)-EijR(10,4)+EijR(12,3
     -   )-EijR(12,4)-EijR(13,3)+EijR(13,4)-2*(EijR(9,3)-EijR(20,3))-
     -   3*(EijR(9,4)-EijR(20,4)-EijR(28,4)+EijR(31,4)))
       FI(89) = 8*(EijI(2,3)+EijI(2,4)-EijI(10,3)-EijI(10,4)+EijI(12,3
     -   )-EijI(12,4)-EijI(13,3)+EijI(13,4)-2*(EijI(9,3)-EijI(20,3))-
     -   3*(EijI(9,4)-EijI(20,4)-EijI(28,4)+EijI(31,4)))
       F(89)=DCMPLX(FR(89),FI(89))
       FR(90) = 8*(EijR(2,3)+EijR(2,4)-EijR(10,3)+EijR(12,3)-EijR(13,3
     -   )+EijR(20,4)+EijR(21,4)+EijR(22,4)+4*EijR(28,4)-2*(EijR(9,3)
     -   +EijR(9,4)+EijR(10,4)-EijR(20,3)+EijR(31,4)+EijR(34,4)))
       FI(90) = 8*(EijI(2,3)+EijI(2,4)-EijI(10,3)+EijI(12,3)-EijI(13,3
     -   )+EijI(20,4)+EijI(21,4)+EijI(22,4)+4*EijI(28,4)-2*(EijI(9,3)
     -   +EijI(9,4)+EijI(10,4)-EijI(20,3)+EijI(31,4)+EijI(34,4)))
       F(90)=DCMPLX(FR(90),FI(90))
       FR(91) = 8*(EijR(2,3)+EijR(2,4)-EijR(3,2)+EijR(8,2)-EijR(9,2)+E
     -   ijR(10,2)-EijR(10,3)-EijR(10,4)+EijR(12,3)-EijR(13,3)+EijR(2
     -   0,4)-2*(EijR(9,3)+EijR(9,4)-EijR(20,3)-EijR(28,4))-EijR(31,4
     -   ))
       FI(91) = 8*(EijI(2,3)+EijI(2,4)-EijI(3,2)+EijI(8,2)-EijI(9,2)+E
     -   ijI(10,2)-EijI(10,3)-EijI(10,4)+EijI(12,3)-EijI(13,3)+EijI(2
     -   0,4)-2*(EijI(9,3)+EijI(9,4)-EijI(20,3)-EijI(28,4))-EijI(31,4
     -   ))
       F(91)=DCMPLX(FR(91),FI(91))
       FR(92) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)-EijR(8,2)+EijR(8,4)-
     -   EijR(9,2)+EijR(9,4)+EijR(10,2)+EijR(10,3)+EijR(14,3)-EijR(16
     -   ,3)-EijR(17,3)-3*EijR(18,3)+EijR(19,3)+EijR(20,3)-EijR(21,4)
     -   -EijR(26,4)+EijR(32,4)-EijR(33,4)+EijR(34,4)+2*(EijR(8,3)+Ei
     -   jR(10,4)-EijR(27,4)-EijR(28,4)+EijR(35,4)))
       FI(92) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)-EijI(8,2)+EijI(8,4)-
     -   EijI(9,2)+EijI(9,4)+EijI(10,2)+EijI(10,3)+EijI(14,3)-EijI(16
     -   ,3)-EijI(17,3)-3*EijI(18,3)+EijI(19,3)+EijI(20,3)-EijI(21,4)
     -   -EijI(26,4)+EijI(32,4)-EijI(33,4)+EijI(34,4)+2*(EijI(8,3)+Ei
     -   jI(10,4)-EijI(27,4)-EijI(28,4)+EijI(35,4)))
       F(92)=DCMPLX(FR(92),FI(92))
       FR(93) = 8*(EijR(2,2)+EijR(2,4)-EijR(8,2)-EijR(9,2)+EijR(10,2)-
     -   3*(EijR(9,3)+EijR(10,3))+EijR(12,3)-EijR(13,3)+EijR(15,3)-Ei
     -   jR(16,3)+EijR(20,4)+EijR(21,4)+EijR(22,4)+4*(EijR(20,3)+EijR
     -   (28,4))+2*(EijR(2,3)-EijR(9,4)-EijR(10,4)-EijR(31,4)-EijR(34
     -   ,4)))
       FI(93) = 8*(EijI(2,2)+EijI(2,4)-EijI(8,2)-EijI(9,2)+EijI(10,2)-
     -   3*(EijI(9,3)+EijI(10,3))+EijI(12,3)-EijI(13,3)+EijI(15,3)-Ei
     -   jI(16,3)+EijI(20,4)+EijI(21,4)+EijI(22,4)+4*(EijI(20,3)+EijI
     -   (28,4))+2*(EijI(2,3)-EijI(9,4)-EijI(10,4)-EijI(31,4)-EijI(34
     -   ,4)))
       F(93)=DCMPLX(FR(93),FI(93))
       FR(94) = 8*(EijR(2,2)+EijR(2,4)-EijR(8,2)-EijR(9,2)-EijR(9,4)+E
     -   ijR(10,2)-EijR(15,4)+2*(EijR(2,3)-EijR(9,3)+EijR(15,3)-EijR(
     -   16,3))+EijR(16,4)-4*(EijR(10,3)-EijR(20,3))-3*(EijR(10,4)-Ei
     -   jR(21,4)-EijR(28,4)+EijR(34,4)))
       FI(94) = 8*(EijI(2,2)+EijI(2,4)-EijI(8,2)-EijI(9,2)-EijI(9,4)+E
     -   ijI(10,2)-EijI(15,4)+2*(EijI(2,3)-EijI(9,3)+EijI(15,3)-EijI(
     -   16,3))+EijI(16,4)-4*(EijI(10,3)-EijI(20,3))-3*(EijI(10,4)-Ei
     -   jI(21,4)-EijI(28,4)+EijI(34,4)))
       F(94)=DCMPLX(FR(94),FI(94))
       FR(95) = 8*(EijR(2,2)+EijR(2,4)+EijR(3,1)-EijR(4,1)-EijR(9,2)-E
     -   ijR(9,4)+EijR(15,3)-EijR(16,3)-3*(EijR(10,3)-EijR(20,3))+Eij
     -   R(21,4)+2*(EijR(2,3)-EijR(9,3)-EijR(10,4)+EijR(28,4))-EijR(3
     -   4,4))
       FI(95) = 8*(EijI(2,2)+EijI(2,4)+EijI(3,1)-EijI(4,1)-EijI(9,2)-E
     -   ijI(9,4)+EijI(15,3)-EijI(16,3)-3*(EijI(10,3)-EijI(20,3))+Eij
     -   I(21,4)+2*(EijI(2,3)-EijI(9,3)-EijI(10,4)+EijI(28,4))-EijI(3
     -   4,4))
       F(95)=DCMPLX(FR(95),FI(95))
       FR(96) = 8*Is34*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(2,3)-Dij1
     -   235R(3,1)+2*(Dij1235R(4,2)-Dij1235R(5,2))+Dij1235R(6,2)+Dij1
     -   235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)+s34*(EijR(2,3)+EijR(2
     -   ,4)+EijR(5,2)-EijR(7,2)+EijR(8,2)-EijR(8,4)-EijR(9,4)-EijR(1
     -   0,2)-EijR(10,3)-EijR(10,4)+EijR(17,3)-2*(EijR(2,2)+EijR(8,3)
     -   -EijR(9,2)-EijR(18,3))-EijR(19,3)+EijR(26,4)+EijR(27,4)+EijR
     -   (28,4)-EijR(35,4)))
       FI(96) = 8*Is34*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(2,3)-Dij1
     -   235I(3,1)+2*(Dij1235I(4,2)-Dij1235I(5,2))+Dij1235I(6,2)+Dij1
     -   235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)+s34*(EijI(2,3)+EijI(2
     -   ,4)+EijI(5,2)-EijI(7,2)+EijI(8,2)-EijI(8,4)-EijI(9,4)-EijI(1
     -   0,2)-EijI(10,3)-EijI(10,4)+EijI(17,3)-2*(EijI(2,2)+EijI(8,3)
     -   -EijI(9,2)-EijI(18,3))-EijI(19,3)+EijI(26,4)+EijI(27,4)+EijI
     -   (28,4)-EijI(35,4)))
       F(96)=DCMPLX(FR(96),FI(96))
       FR(97) = -8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,1)-EijR(2,4)-EijR(3,1)+EijR(10,4)-E
     -   ijR(12,3)+EijR(13,3)+3*(EijR(9,3)-EijR(20,3))-EijR(20,4)-2*(
     -   EijR(2,3)-EijR(9,4)-EijR(10,3)+EijR(28,4))+EijR(31,4)))
       FI(97) = -8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,1)-EijI(2,4)-EijI(3,1)+EijI(10,4)-E
     -   ijI(12,3)+EijI(13,3)+3*(EijI(9,3)-EijI(20,3))-EijI(20,4)-2*(
     -   EijI(2,3)-EijI(9,4)-EijI(10,3)+EijI(28,4))+EijI(31,4)))
       F(97)=DCMPLX(FR(97),FI(97))
       FR(98) = -8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,1)-EijR(2,4)-EijR(4,1)+EijR(8,2)+Ei
     -   jR(9,4)-EijR(10,2)-EijR(15,3)+EijR(16,3)+3*(EijR(10,3)-EijR(
     -   20,3))-EijR(21,4)-2*(EijR(2,3)-EijR(9,3)-EijR(10,4)+EijR(28,
     -   4))+EijR(34,4)))
       FI(98) = -8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,1)-EijI(2,4)-EijI(4,1)+EijI(8,2)+Ei
     -   jI(9,4)-EijI(10,2)-EijI(15,3)+EijI(16,3)+3*(EijI(10,3)-EijI(
     -   20,3))-EijI(21,4)-2*(EijI(2,3)-EijI(9,3)-EijI(10,4)+EijI(28,
     -   4))+EijI(34,4)))
       F(98)=DCMPLX(FR(98),FI(98))
       FR(99) = -8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+2*(
     -   Dij1235R(2,2)-Dij1235R(6,2))-Dij1235R(8,3)+s34*(EijR(2,1)-Ei
     -   jR(2,4)-EijR(3,1)+EijR(9,2)+EijR(9,4)-EijR(10,2)+EijR(10,4)-
     -   2*(EijR(2,3)-EijR(9,3)-EijR(10,3)+EijR(20,3))-EijR(28,4)))
       FI(99) = -8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+2*(
     -   Dij1235I(2,2)-Dij1235I(6,2))-Dij1235I(8,3)+s34*(EijI(2,1)-Ei
     -   jI(2,4)-EijI(3,1)+EijI(9,2)+EijI(9,4)-EijI(10,2)+EijI(10,4)-
     -   2*(EijI(2,3)-EijI(9,3)-EijI(10,3)+EijI(20,3))-EijI(28,4)))
       F(99)=DCMPLX(FR(99),FI(99))
       FR(100) = -4*(EijR(2,2)-EijR(5,2)+EijR(5,3)+EijR(6,2)-EijR(6,3)
     -   -EijR(8,2)-4*(EijR(8,3)-EijR(8,4))+3*(EijR(2,3)-EijR(9,3))-2
     -   *(EijR(2,4)-EijR(9,4)+EijR(17,4)-EijR(23,4)-2*(EijR(17,3)-Ei
     -   jR(26,4))))
       FI(100) = -4*(EijI(2,2)-EijI(5,2)+EijI(5,3)+EijI(6,2)-EijI(6,3)
     -   -EijI(8,2)-4*(EijI(8,3)-EijI(8,4))+3*(EijI(2,3)-EijI(9,3))-2
     -   *(EijI(2,4)-EijI(9,4)+EijI(17,4)-EijI(23,4)-2*(EijI(17,3)-Ei
     -   jI(26,4))))
       F(100)=DCMPLX(FR(100),FI(100))
       FR(101) = -4*(EijR(2,3)-EijR(8,3)-EijR(11,3)+EijR(12,3)+4*(EijR
     -   (9,4)-EijR(26,4))+2*(EijR(2,2)-EijR(2,4)-EijR(5,2)+EijR(6,2)
     -   -EijR(8,2)+EijR(8,4)-EijR(9,3)+EijR(17,3)-EijR(20,4)+EijR(29
     -   ,4)))
       FI(101) = -4*(EijI(2,3)-EijI(8,3)-EijI(11,3)+EijI(12,3)+4*(EijI
     -   (9,4)-EijI(26,4))+2*(EijI(2,2)-EijI(2,4)-EijI(5,2)+EijI(6,2)
     -   -EijI(8,2)+EijI(8,4)-EijI(9,3)+EijI(17,3)-EijI(20,4)+EijI(29
     -   ,4)))
       F(101)=DCMPLX(FR(101),FI(101))
       FR(102) = -4*(EijR(2,3)-EijR(6,2)+EijR(8,2)-EijR(8,3)+EijR(9,3)
     -   -EijR(17,3)+3*(EijR(7,2)-EijR(9,2)-EijR(10,3)+EijR(18,3))-Ei
     -   jR(19,3)+EijR(20,3)+2*(EijR(2,2)-EijR(2,4)-EijR(5,2)+EijR(8,
     -   4)+EijR(9,4)+EijR(10,4)-EijR(26,4)-EijR(27,4)-EijR(28,4)+Eij
     -   R(35,4)))
       FI(102) = -4*(EijI(2,3)-EijI(6,2)+EijI(8,2)-EijI(8,3)+EijI(9,3)
     -   -EijI(17,3)+3*(EijI(7,2)-EijI(9,2)-EijI(10,3)+EijI(18,3))-Ei
     -   jI(19,3)+EijI(20,3)+2*(EijI(2,2)-EijI(2,4)-EijI(5,2)+EijI(8,
     -   4)+EijI(9,4)+EijI(10,4)-EijI(26,4)-EijI(27,4)-EijI(28,4)+Eij
     -   I(35,4)))
       F(102)=DCMPLX(FR(102),FI(102))
       FR(103) = -4*(EijR(2,3)+EijR(6,2)-EijR(8,2)-EijR(8,3)+EijR(9,3)
     -   -EijR(17,3)+2*(EijR(2,2)-EijR(2,4)-EijR(5,2)+EijR(8,4)+EijR(
     -   9,4)-EijR(26,4)))
       FI(103) = -4*(EijI(2,3)+EijI(6,2)-EijI(8,2)-EijI(8,3)+EijI(9,3)
     -   -EijI(17,3)+2*(EijI(2,2)-EijI(2,4)-EijI(5,2)+EijI(8,4)+EijI(
     -   9,4)-EijI(26,4)))
       F(103)=DCMPLX(FR(103),FI(103))
       FR(104) = -4*(EijR(2,1)+EijR(2,3)-EijR(3,1)+EijR(3,2)+EijR(5,2)
     -   -EijR(6,2)+EijR(8,3)-EijR(11,3)-3*(EijR(8,2)-EijR(12,3))-4*(
     -   EijR(9,3)-EijR(9,4)+EijR(26,4))+2*(EijR(2,2)-EijR(2,4)+EijR(
     -   8,4)-EijR(20,4)+EijR(29,4)))
       FI(104) = -4*(EijI(2,1)+EijI(2,3)-EijI(3,1)+EijI(3,2)+EijI(5,2)
     -   -EijI(6,2)+EijI(8,3)-EijI(11,3)-3*(EijI(8,2)-EijI(12,3))-4*(
     -   EijI(9,3)-EijI(9,4)+EijI(26,4))+2*(EijI(2,2)-EijI(2,4)+EijI(
     -   8,4)-EijI(20,4)+EijI(29,4)))
       F(104)=DCMPLX(FR(104),FI(104))
       FR(105) = -4*(EijR(2,2)-EijR(2,3)+EijR(3,2)-EijR(3,3)+EijR(9,3)
     -   +EijR(12,3)-2*(EijR(2,4)+EijR(8,2)-EijR(12,4))+6*(EijR(9,4)-
     -   EijR(20,4)))
       FI(105) = -4*(EijI(2,2)-EijI(2,3)+EijI(3,2)-EijI(3,3)+EijI(9,3)
     -   +EijI(12,3)-2*(EijI(2,4)+EijI(8,2)-EijI(12,4))+6*(EijI(9,4)-
     -   EijI(20,4)))
       F(105)=DCMPLX(FR(105),FI(105))
       FR(106) = 4*(-EijR(2,2)+EijR(2,3)+EijR(3,2)+EijR(10,3)+EijR(12,
     -   3)+EijR(13,3)-4*(EijR(9,4)-EijR(28,4))+2*(EijR(2,4)+EijR(9,2
     -   )-EijR(9,3)-EijR(10,2)-EijR(10,4)-EijR(20,3)+EijR(20,4)-EijR
     -   (31,4)))
       FI(106) = 4*(-EijI(2,2)+EijI(2,3)+EijI(3,2)+EijI(10,3)+EijI(12,
     -   3)+EijI(13,3)-4*(EijI(9,4)-EijI(28,4))+2*(EijI(2,4)+EijI(9,2
     -   )-EijI(9,3)-EijI(10,2)-EijI(10,4)-EijI(20,3)+EijI(20,4)-EijI
     -   (31,4)))
       F(106)=DCMPLX(FR(106),FI(106))
       FR(107) = -4*(EijR(2,2)-EijR(2,3)+EijR(3,2)+4*EijR(9,4)-EijR(12
     -   ,3)-2*(EijR(2,4)+EijR(8,2)-EijR(9,3)+EijR(20,4)))
       FI(107) = -4*(EijI(2,2)-EijI(2,3)+EijI(3,2)+4*EijI(9,4)-EijI(12
     -   ,3)-2*(EijI(2,4)+EijI(8,2)-EijI(9,3)+EijI(20,4)))
       F(107)=DCMPLX(FR(107),FI(107))
       FR(108) = -4*(EijR(2,1)+EijR(2,2)-EijR(2,3)-EijR(3,1)+EijR(6,2)
     -   -4*EijR(8,2)-EijR(9,3)+EijR(10,2)-EijR(10,3)-EijR(17,3)-EijR
     -   (18,3)-EijR(19,3)-3*(EijR(7,2)-EijR(8,3)-EijR(20,3))-2*(EijR
     -   (2,4)-EijR(5,2)-EijR(8,4)-EijR(9,2)-EijR(9,4)-EijR(10,4)+Eij
     -   R(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4)))
       FI(108) = -4*(EijI(2,1)+EijI(2,2)-EijI(2,3)-EijI(3,1)+EijI(6,2)
     -   -4*EijI(8,2)-EijI(9,3)+EijI(10,2)-EijI(10,3)-EijI(17,3)-EijI
     -   (18,3)-EijI(19,3)-3*(EijI(7,2)-EijI(8,3)-EijI(20,3))-2*(EijI
     -   (2,4)-EijI(5,2)-EijI(8,4)-EijI(9,2)-EijI(9,4)-EijI(10,4)+Eij
     -   I(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4)))
       F(108)=DCMPLX(FR(108),FI(108))
       FR(109) = 4*(EijR(2,1)+3*EijR(2,3)-EijR(3,1)-EijR(3,2)-EijR(8,2
     -   )-EijR(9,2)+EijR(10,2)-EijR(10,3)+EijR(12,3)+EijR(13,3)-4*(E
     -   ijR(9,3)+EijR(9,4)-EijR(28,4))+2*(EijR(2,2)+EijR(2,4)-EijR(1
     -   0,4)+EijR(20,4)-EijR(31,4)))
       FI(109) = 4*(EijI(2,1)+3*EijI(2,3)-EijI(3,1)-EijI(3,2)-EijI(8,2
     -   )-EijI(9,2)+EijI(10,2)-EijI(10,3)+EijI(12,3)+EijI(13,3)-4*(E
     -   ijI(9,3)+EijI(9,4)-EijI(28,4))+2*(EijI(2,2)+EijI(2,4)-EijI(1
     -   0,4)+EijI(20,4)-EijI(31,4)))
       F(109)=DCMPLX(FR(109),FI(109))
       FR(110) = 4*(EijR(2,1)-EijR(3,1)+3*(EijR(2,3)-EijR(9,3))-EijR(1
     -   5,3)+EijR(16,3)-4*(EijR(10,4)-EijR(28,4))+2*(EijR(2,2)+EijR(
     -   2,4)-EijR(8,2)-EijR(9,4)-EijR(10,3)+EijR(20,3)+EijR(21,4)-Ei
     -   jR(34,4)))
       FI(110) = 4*(EijI(2,1)-EijI(3,1)+3*(EijI(2,3)-EijI(9,3))-EijI(1
     -   5,3)+EijI(16,3)-4*(EijI(10,4)-EijI(28,4))+2*(EijI(2,2)+EijI(
     -   2,4)-EijI(8,2)-EijI(9,4)-EijI(10,3)+EijI(20,3)+EijI(21,4)-Ei
     -   jI(34,4)))
       F(110)=DCMPLX(FR(110),FI(110))
       P(481) = s45+2*P(6)
       P(482) = 3*s12+2*P(39)
       P(483) = 3*s12+2*P(14)
       P(484) = 3*p5sq-s34-2*P(161)
       P(485) = s34-s45+2*P(161)
       P(486) = -s23+s45+2*P(434)
       P(487) = p5sq+s15+s23-s45
       P(488) = p2sq+p5sq-s12-s15+s23-s45
       FR(111) = 2*(-Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(2,1)+Dij2345
     -   R(4,2)+Dij2345R(5,2)-Dij2345R(6,2)-4*(-2*EijR(22,3)+EijR(23,
     -   3)+EijR(24,3)-EijR(37,4)+EijR(43,4)+EijR(44,4)-EijR(45,4))-2
     -   *(s15*EijR(3,1)-p5sq*EijR(4,1)+(EijR(2,1)+EijR(5,2)-EijR(7,2
     -   ))*P(3))+(EijR(3,2)+EijR(12,3)-EijR(13,3))*P(10)-(EijR(9,2)-
     -   EijR(15,3)+EijR(16,3))*P(21)+(EijR(2,2)+EijR(2,3))*P(33)-Eij
     -   R(9,3)*P(245)+EijR(20,3)*P(332)+Is12*(C0345R+Cij345R(2,1)-12
     -   *Dij1345R(7,2)-4*(Dij1345R(11,3)-Dij1345R(13,3))-Dij1345R(3,
     -   2)*P(7)-2*(p3sq*Dij1345R(2,2)+Dij1345R(1,1)*P(35))+Dij1345R(
     -   3,1)*P(185)+Dij1345R(2,1)*P(481)-Dij1345R(1,2)*P(482)+Dij134
     -   5R(4,2)*P(483)+Dij1345R(5,2)*P(484)+Dij1345R(6,2)*P(485))-Ei
     -   jR(8,2)*P(486)+EijR(10,2)*P(487)-EijR(10,3)*P(488))
       FI(111) = 2*(-Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(2,1)+Dij2345
     -   I(4,2)+Dij2345I(5,2)-Dij2345I(6,2)-4*(-2*EijI(22,3)+EijI(23,
     -   3)+EijI(24,3)-EijI(37,4)+EijI(43,4)+EijI(44,4)-EijI(45,4))-2
     -   *(s15*EijI(3,1)-p5sq*EijI(4,1)+(EijI(2,1)+EijI(5,2)-EijI(7,2
     -   ))*P(3))+(EijI(3,2)+EijI(12,3)-EijI(13,3))*P(10)-(EijI(9,2)-
     -   EijI(15,3)+EijI(16,3))*P(21)+(EijI(2,2)+EijI(2,3))*P(33)-Eij
     -   I(9,3)*P(245)+EijI(20,3)*P(332)+Is12*(C0345I+Cij345I(2,1)-12
     -   *Dij1345I(7,2)-4*(Dij1345I(11,3)-Dij1345I(13,3))-Dij1345I(3,
     -   2)*P(7)-2*(p3sq*Dij1345I(2,2)+Dij1345I(1,1)*P(35))+Dij1345I(
     -   3,1)*P(185)+Dij1345I(2,1)*P(481)-Dij1345I(1,2)*P(482)+Dij134
     -   5I(4,2)*P(483)+Dij1345I(5,2)*P(484)+Dij1345I(6,2)*P(485))-Ei
     -   jI(8,2)*P(486)+EijI(10,2)*P(487)-EijI(10,3)*P(488))
       F(111)=DCMPLX(FR(111),FI(111))
       P(489) = p3sq+s45+2*P(6)
       P(490) = p3sq+2*p5sq+5*s12-4*s34+s45
       P(491) = p3sq-p5sq-s34+s45
       P(492) = p2sq-s12+2*s15-s23+s45
       P(493) = 3*P(33)-2*P(50)
       FR(112) = 2*(-Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(2,1)-Dij2345
     -   R(2,2)+2*(Dij2345R(4,2)+2*(EijR(22,3)-EijR(23,3)+EijR(37,4)+
     -   EijR(38,4)-2*EijR(43,4))-(EijR(2,1)-EijR(3,1)+EijR(5,2)-EijR
     -   (6,2))*P(3))-EijR(3,3)*P(10)-(EijR(9,2)-EijR(10,2)+EijR(10,3
     -   )+EijR(13,3)-2*EijR(20,3))*P(21)+(EijR(2,2)+EijR(2,3))*P(33)
     -   -EijR(9,3)*P(302)-EijR(8,2)*P(486)+Is12*(C0345R+Cij345R(1,1)
     -   -6*Dij1345R(7,2)-4*(Dij1345R(11,3)-Dij1345R(12,3))+2*(-Dij13
     -   45R(1,1)+Dij1345R(2,1))*P(35)-Dij1345R(1,2)*P(482)-Dij1345R(
     -   2,2)*P(489)+Dij1345R(4,2)*P(490)+(-Dij1345R(5,2)+Dij1345R(6,
     -   2))*P(491))+EijR(3,2)*P(492)+EijR(12,3)*P(493))
       FI(112) = 2*(-Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(2,1)-Dij2345
     -   I(2,2)+2*(Dij2345I(4,2)+2*(EijI(22,3)-EijI(23,3)+EijI(37,4)+
     -   EijI(38,4)-2*EijI(43,4))-(EijI(2,1)-EijI(3,1)+EijI(5,2)-EijI
     -   (6,2))*P(3))-EijI(3,3)*P(10)-(EijI(9,2)-EijI(10,2)+EijI(10,3
     -   )+EijI(13,3)-2*EijI(20,3))*P(21)+(EijI(2,2)+EijI(2,3))*P(33)
     -   -EijI(9,3)*P(302)-EijI(8,2)*P(486)+Is12*(C0345I+Cij345I(1,1)
     -   -6*Dij1345I(7,2)-4*(Dij1345I(11,3)-Dij1345I(12,3))+2*(-Dij13
     -   45I(1,1)+Dij1345I(2,1))*P(35)-Dij1345I(1,2)*P(482)-Dij1345I(
     -   2,2)*P(489)+Dij1345I(4,2)*P(490)+(-Dij1345I(5,2)+Dij1345I(6,
     -   2))*P(491))+EijI(3,2)*P(492)+EijI(12,3)*P(493))
       F(112)=DCMPLX(FR(112),FI(112))
       FR(113) = -2*(D02345R+Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(4,2)
     -   +p5sq*(-2*(EijR(1,1)+EijR(1,2)-EijR(2,1)-EijR(5,2))-EijR(7,2
     -   )+EijR(9,2)+EijR(10,3)-EijR(18,3)+EijR(19,3)-EijR(20,3))+s15
     -   *(EijR(7,2)+2*(EijR(1,1)+EijR(1,2)-EijR(2,1)-EijR(5,2)-EijR(
     -   6,2)+EijR(8,2))-EijR(9,2)-EijR(10,3)+EijR(18,3)-EijR(19,3)+E
     -   ijR(20,3))+2*(Dij2345R(1,1)-3*EijR(11,2)-EijR(22,3)+EijR(23,
     -   3)-2*(EijR(37,4)-EijR(40,4)+EijR(41,4)-EijR(43,4)))+(-EijR(2
     -   ,2)-EijR(2,3)+EijR(5,2)-EijR(6,2)+EijR(8,2)+EijR(8,3)+EijR(1
     -   1,3)-EijR(12,3)+2*(EijR(9,3)-EijR(17,3)))*P(33)+(EijR(6,2)-E
     -   ijR(7,2)-EijR(8,2)+EijR(9,2)-EijR(9,3)+EijR(10,3)-EijR(11,3)
     -   +EijR(12,3)+EijR(17,3)-EijR(18,3)+EijR(19,3)-EijR(20,3))*P(5
     -   0))
       FI(113) = -2*(D02345I+Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(4,2)
     -   +p5sq*(-2*(EijI(1,1)+EijI(1,2)-EijI(2,1)-EijI(5,2))-EijI(7,2
     -   )+EijI(9,2)+EijI(10,3)-EijI(18,3)+EijI(19,3)-EijI(20,3))+s15
     -   *(EijI(7,2)+2*(EijI(1,1)+EijI(1,2)-EijI(2,1)-EijI(5,2)-EijI(
     -   6,2)+EijI(8,2))-EijI(9,2)-EijI(10,3)+EijI(18,3)-EijI(19,3)+E
     -   ijI(20,3))+2*(Dij2345I(1,1)-3*EijI(11,2)-EijI(22,3)+EijI(23,
     -   3)-2*(EijI(37,4)-EijI(40,4)+EijI(41,4)-EijI(43,4)))+(-EijI(2
     -   ,2)-EijI(2,3)+EijI(5,2)-EijI(6,2)+EijI(8,2)+EijI(8,3)+EijI(1
     -   1,3)-EijI(12,3)+2*(EijI(9,3)-EijI(17,3)))*P(33)+(EijI(6,2)-E
     -   ijI(7,2)-EijI(8,2)+EijI(9,2)-EijI(9,3)+EijI(10,3)-EijI(11,3)
     -   +EijI(12,3)+EijI(17,3)-EijI(18,3)+EijI(19,3)-EijI(20,3))*P(5
     -   0))
       F(113)=DCMPLX(FR(113),FI(113))
       FR(114) = 2*(Dij1345R(1,1)+3*Dij1345R(1,2)-Dij1345R(2,1)+Dij134
     -   5R(2,2)+2*(Dij1235R(3,2)-Dij1235R(6,2)+s23*(-EijR(2,1)+EijR(
     -   3,1)-EijR(5,2)+EijR(6,2)))-4*(Dij1345R(4,2)-EijR(22,3)+EijR(
     -   23,3)-EijR(37,4)-EijR(38,4)+2*EijR(43,4))-Is12*Is34*(2*B035R
     -   +C0235R*s12+4*Cij135R(4,2)+p5sq*(2*Cij135R(1,1)+s12*(Dij1235
     -   R(3,1)+Dij1235R(3,2)-Dij1235R(6,2)))-s12*(-5*Cij135R(1,1)+3*
     -   Cij135R(2,1)-2*(Cij135R(1,2)-Cij135R(3,2))+Cij235R(1,1)-Cij2
     -   35R(2,1)+p2sq*(Dij1235R(2,2)+Dij1235R(3,1)-Dij1235R(3,2)-2*(
     -   Dij1235R(4,2)-Dij1235R(5,2)))+s12*(Dij1235R(2,2)-Dij1235R(3,
     -   1)+Dij1235R(3,2)+2*(Dij1235R(4,2)-Dij1235R(5,2)-Dij1235R(6,2
     -   )))-s15*(Dij1235R(3,1)-Dij1235R(3,2)+Dij1235R(6,2))+4*(Dij12
     -   35R(7,2)-Dij1235R(12,3)+Dij1235R(13,3)))+s34*(-C0135R-D01345
     -   R*s45-Cij135R(1,1)-Cij135R(2,1)+Cij345R(1,1)-p5sq*Dij1345R(3
     -   ,1)+p3sq*(Dij1345R(1,1)+Dij1345R(1,2)-Dij1345R(2,1)-Dij1345R
     -   (2,2)-2*(Dij1345R(5,2)-Dij1345R(6,2)))-(p5sq-2*s34)*(Dij1345
     -   R(5,2)-Dij1345R(6,2))-4*(Dij1345R(7,2)-Dij1345R(11,3)+Dij134
     -   5R(12,3))+s45*(-Dij1345R(1,2)+Dij1345R(2,2)+Dij1345R(3,1)+Di
     -   j1345R(5,2)-Dij1345R(6,2)-2*(Dij1345R(1,1)-Dij1345R(2,1)-s12
     -   *(-EijR(2,1)+EijR(3,1)-EijR(5,2)+EijR(6,2))))+2*(C0345R+Cij1
     -   35R(2,2)-Cij135R(3,2)+s12*(s15-s23)*EijR(9,2)))+s12*(Dij1235
     -   R(2,1)*P(33)+2*((D01235R+Dij1235R(1,1))*P(3)-s34*EijR(10,2)*
     -   P(66)))+C0135R*P(212)))
       FI(114) = 2*(Dij1345I(1,1)+3*Dij1345I(1,2)-Dij1345I(2,1)+Dij134
     -   5I(2,2)+2*(Dij1235I(3,2)-Dij1235I(6,2)+s23*(-EijI(2,1)+EijI(
     -   3,1)-EijI(5,2)+EijI(6,2)))-4*(Dij1345I(4,2)-EijI(22,3)+EijI(
     -   23,3)-EijI(37,4)-EijI(38,4)+2*EijI(43,4))-Is12*Is34*(2*B035I
     -   +C0235I*s12+4*Cij135I(4,2)+p5sq*(2*Cij135I(1,1)+s12*(Dij1235
     -   I(3,1)+Dij1235I(3,2)-Dij1235I(6,2)))-s12*(-5*Cij135I(1,1)+3*
     -   Cij135I(2,1)-2*(Cij135I(1,2)-Cij135I(3,2))+Cij235I(1,1)-Cij2
     -   35I(2,1)+p2sq*(Dij1235I(2,2)+Dij1235I(3,1)-Dij1235I(3,2)-2*(
     -   Dij1235I(4,2)-Dij1235I(5,2)))+s12*(Dij1235I(2,2)-Dij1235I(3,
     -   1)+Dij1235I(3,2)+2*(Dij1235I(4,2)-Dij1235I(5,2)-Dij1235I(6,2
     -   )))-s15*(Dij1235I(3,1)-Dij1235I(3,2)+Dij1235I(6,2))+4*(Dij12
     -   35I(7,2)-Dij1235I(12,3)+Dij1235I(13,3)))+s34*(-C0135I-D01345
     -   I*s45-Cij135I(1,1)-Cij135I(2,1)+Cij345I(1,1)-p5sq*Dij1345I(3
     -   ,1)+p3sq*(Dij1345I(1,1)+Dij1345I(1,2)-Dij1345I(2,1)-Dij1345I
     -   (2,2)-2*(Dij1345I(5,2)-Dij1345I(6,2)))-(p5sq-2*s34)*(Dij1345
     -   I(5,2)-Dij1345I(6,2))-4*(Dij1345I(7,2)-Dij1345I(11,3)+Dij134
     -   5I(12,3))+s45*(-Dij1345I(1,2)+Dij1345I(2,2)+Dij1345I(3,1)+Di
     -   j1345I(5,2)-Dij1345I(6,2)-2*(Dij1345I(1,1)-Dij1345I(2,1)-s12
     -   *(-EijI(2,1)+EijI(3,1)-EijI(5,2)+EijI(6,2))))+2*(C0345I+Cij1
     -   35I(2,2)-Cij135I(3,2)+s12*(s15-s23)*EijI(9,2)))+s12*(Dij1235
     -   I(2,1)*P(33)+2*((D01235I+Dij1235I(1,1))*P(3)-s34*EijI(10,2)*
     -   P(66)))+C0135I*P(212)))
       F(114)=DCMPLX(FR(114),FI(114))
       P(494) = p2sq-s12-2*s15
       P(495) = p2sq+2*p5sq-s12-4*s15
       FR(115) = 2*(D02345R+Dij2345R(1,1)-6*EijR(11,2)-4*(EijR(22,3)-E
     -   ijR(23,3)-EijR(37,4)+EijR(40,4)-EijR(41,4)+EijR(43,4))-2*(Ei
     -   jR(1,1)+EijR(1,2)-EijR(2,1))*P(3)+(-EijR(6,2)+EijR(8,2))*P(1
     -   0)+(-EijR(7,2)+EijR(9,2))*P(21)-EijR(2,2)*P(494)+EijR(5,2)*P
     -   (495))
       FI(115) = 2*(D02345I+Dij2345I(1,1)-6*EijI(11,2)-4*(EijI(22,3)-E
     -   ijI(23,3)-EijI(37,4)+EijI(40,4)-EijI(41,4)+EijI(43,4))-2*(Ei
     -   jI(1,1)+EijI(1,2)-EijI(2,1))*P(3)+(-EijI(6,2)+EijI(8,2))*P(1
     -   0)+(-EijI(7,2)+EijI(9,2))*P(21)-EijI(2,2)*P(494)+EijI(5,2)*P
     -   (495))
       F(115)=DCMPLX(FR(115),FI(115))
       P(496) = 2*p5sq-s15
       FR(116) = 2*(D01235R+Dij1235R(1,1)-2*Dij1235R(3,1)+Dij1235R(5,2
     -   )-Dij1235R(6,2)-Dij2345R(2,1)+s23*(EijR(1,1)+EijR(1,2))-s45*
     -   (EE0R+2*EijR(1,1)+EijR(1,2))+p3sq*(-EijR(2,1)+EijR(3,1)-EijR
     -   (3,2)-EijR(5,2)+EijR(6,2)+EijR(8,2))-4*(EijR(11,2)+EijR(21,3
     -   )-EijR(22,3)-EijR(37,4)+EijR(40,4)-EijR(41,4)+EijR(43,4)))+I
     -   s15s34*(4*(B025R+C0235R*s15+p5sq*Cij125R(1,1))+8*Cij125R(4,2
     -   )-2*(s15*(D01235R*s12+Cij125R(1,1)+Cij125R(2,1)-Cij235R(1,1)
     -   +p2sq*(-Dij1235R(1,1)-Dij1235R(1,2)+Dij1235R(2,1)+Dij1235R(2
     -   ,2)+2*(Dij1235R(5,2)-Dij1235R(6,2)))+p5sq*(Dij1235R(3,1)+Dij
     -   1235R(5,2)-Dij1235R(6,2))+4*(Dij1235R(7,2)-Dij1235R(11,3)+Di
     -   j1235R(12,3))+s34*(-(s45*(-EijR(3,2)+EijR(4,1)+2*(EijR(2,1)+
     -   EijR(5,2))))+s23*(EijR(2,1)+EijR(5,2)-EijR(6,2))+p5sq*(EijR(
     -   4,1)-EijR(7,2)))+s34*s45*EijR(7,2)+s23*s34*EijR(8,2)+s12*(Di
     -   j1235R(1,2)+2*(Dij1235R(1,1)-Dij1235R(2,1))-Dij1235R(2,2)-Di
     -   j1235R(3,1)-Dij1235R(5,2)+Dij1235R(6,2)+s34*(-EijR(2,1)+EijR
     -   (3,1)-EijR(3,2)-EijR(5,2)+EijR(6,2)+EijR(8,2)))-2*(Cij125R(2
     -   ,2)-Cij125R(3,2)+s15*Dij1235R(5,2)-s15*Dij1235R(6,2)+s15*s34
     -   *EijR(7,2)-s23*s34*EijR(7,2)-(s15-s23)*s34*EijR(9,2)))+s15*s
     -   34*EijR(10,2)*P(7)-C0125R*P(496)))
       FI(116) = 2*(D01235I+Dij1235I(1,1)-2*Dij1235I(3,1)+Dij1235I(5,2
     -   )-Dij1235I(6,2)-Dij2345I(2,1)+s23*(EijI(1,1)+EijI(1,2))-s45*
     -   (EE0I+2*EijI(1,1)+EijI(1,2))+p3sq*(-EijI(2,1)+EijI(3,1)-EijI
     -   (3,2)-EijI(5,2)+EijI(6,2)+EijI(8,2))-4*(EijI(11,2)+EijI(21,3
     -   )-EijI(22,3)-EijI(37,4)+EijI(40,4)-EijI(41,4)+EijI(43,4)))+I
     -   s15s34*(4*(B025I+C0235I*s15+p5sq*Cij125I(1,1))+8*Cij125I(4,2
     -   )-2*(s15*(D01235I*s12+Cij125I(1,1)+Cij125I(2,1)-Cij235I(1,1)
     -   +p2sq*(-Dij1235I(1,1)-Dij1235I(1,2)+Dij1235I(2,1)+Dij1235I(2
     -   ,2)+2*(Dij1235I(5,2)-Dij1235I(6,2)))+p5sq*(Dij1235I(3,1)+Dij
     -   1235I(5,2)-Dij1235I(6,2))+4*(Dij1235I(7,2)-Dij1235I(11,3)+Di
     -   j1235I(12,3))+s34*(-(s45*(-EijI(3,2)+EijI(4,1)+2*(EijI(2,1)+
     -   EijI(5,2))))+s23*(EijI(2,1)+EijI(5,2)-EijI(6,2))+p5sq*(EijI(
     -   4,1)-EijI(7,2)))+s34*s45*EijI(7,2)+s23*s34*EijI(8,2)+s12*(Di
     -   j1235I(1,2)+2*(Dij1235I(1,1)-Dij1235I(2,1))-Dij1235I(2,2)-Di
     -   j1235I(3,1)-Dij1235I(5,2)+Dij1235I(6,2)+s34*(-EijI(2,1)+EijI
     -   (3,1)-EijI(3,2)-EijI(5,2)+EijI(6,2)+EijI(8,2)))-2*(Cij125I(2
     -   ,2)-Cij125I(3,2)+s15*Dij1235I(5,2)-s15*Dij1235I(6,2)+s15*s34
     -   *EijI(7,2)-s23*s34*EijI(7,2)-(s15-s23)*s34*EijI(9,2)))+s15*s
     -   34*EijI(10,2)*P(7)-C0125I*P(496)))
       F(116)=DCMPLX(FR(116),FI(116))
       FR(117) = 4*(p5sq*(EijR(2,1)-EijR(3,1)+EijR(5,2)-EijR(6,2))+s15
     -   *(EijR(2,2)-EijR(5,2)+EijR(6,2)-EijR(8,2))+2*(EijR(37,4)+Eij
     -   R(38,4)-2*EijR(43,4)))
       FI(117) = 4*(p5sq*(EijI(2,1)-EijI(3,1)+EijI(5,2)-EijI(6,2))+s15
     -   *(EijI(2,2)-EijI(5,2)+EijI(6,2)-EijI(8,2))+2*(EijI(37,4)+Eij
     -   I(38,4)-2*EijI(43,4)))
       F(117)=DCMPLX(FR(117),FI(117))
       FR(118) = 4*(p5sq*(EijR(2,1)-EijR(4,1)+EijR(5,2)-EijR(7,2))+s15
     -   *(EijR(2,2)-EijR(5,2)+EijR(7,2)-EijR(9,2))-2*(EijR(23,3)-Eij
     -   R(24,3)-EijR(37,4)+EijR(43,4)+EijR(44,4)-EijR(45,4)))
       FI(118) = 4*(p5sq*(EijI(2,1)-EijI(4,1)+EijI(5,2)-EijI(7,2))+s15
     -   *(EijI(2,2)-EijI(5,2)+EijI(7,2)-EijI(9,2))-2*(EijI(23,3)-Eij
     -   I(24,3)-EijI(37,4)+EijI(43,4)+EijI(44,4)-EijI(45,4)))
       F(118)=DCMPLX(FR(118),FI(118))
       P(497) = p3sq+2*P(411)
       P(498) = p3sq-s34-2*s45
       P(499) = -p3sq+3*s12+s45
       P(500) = 4*p3sq+p5sq-3*s45-2*P(29)
       P(501) = -2*p3sq+s12+s34-s45
       P(502) = p2sq+p5sq-s12+s15-s23-s45
       FR(119) = -2*(D02345R+Dij2345R(1,1)-Dij2345R(3,1)-4*(EijR(11,2)
     -   +EijR(22,3)-EijR(24,3)+EijR(37,4)-EijR(43,4)-EijR(44,4)+EijR
     -   (45,4))+(EijR(3,1)+EijR(8,2)-EijR(10,2))*P(10)-EijR(4,2)*P(2
     -   3)-EijR(2,2)*P(33)+2*((EE0R+EijR(1,1))*P(3)+(EijR(5,2)-EijR(
     -   7,2))*P(50))+EijR(4,1)*P(110)-EijR(2,1)*P(111)+Is34*(C0235R-
     -   Cij235R(1,1)+Cij235R(2,1)-4*(Dij1235R(7,2)-Dij1235R(12,3)+Di
     -   j1235R(13,3))+2*(Cij135R(1,2)+(D01235R+Dij1235R(1,1))*P(3))-
     -   Dij1235R(2,2)*P(5)+Is12*(4*Cij135R(4,2)+C0135R*P(26)+Cij135R
     -   (1,1)*P(27)-Cij135R(2,1)*P(28)+2*(B035R-Cij135R(3,2)*P(29)))
     -   +(Dij1235R(2,1)+2*(Dij1235R(4,2)-Dij1235R(5,2)))*P(33)-Dij12
     -   35R(3,1)*P(34)+Dij1235R(6,2)*P(42)-Dij1235R(3,2)*P(142))+Is1
     -   2*(3*(C0345R-Cij145R(2,1))+Cij345R(2,1)-8*Dij1345R(7,2)+4*(D
     -   ij1345R(11,3)-Dij1345R(13,3))-Dij1345R(2,1)*P(16)+Is45*(2*B0
     -   12R+4*Cij145R(4,2)+C0145R*P(95)+Cij145R(1,1)*P(96))+D01345R*
     -   P(97)+Dij1345R(3,2)*P(101)+2*(Cij135R(2,2)+Cij145R(1,2)-Cij1
     -   45R(3,2)-p3sq*Dij1345R(2,2)+Dij1345R(4,2)*P(198))+Dij1345R(1
     -   ,1)*P(497)-Dij1345R(3,1)*P(498)-Dij1345R(1,2)*P(499)-Dij1345
     -   R(5,2)*P(500)-Dij1345R(6,2)*P(501))+EijR(9,2)*P(502))
       FI(119) = -2*(D02345I+Dij2345I(1,1)-Dij2345I(3,1)-4*(EijI(11,2)
     -   +EijI(22,3)-EijI(24,3)+EijI(37,4)-EijI(43,4)-EijI(44,4)+EijI
     -   (45,4))+(EijI(3,1)+EijI(8,2)-EijI(10,2))*P(10)-EijI(4,2)*P(2
     -   3)-EijI(2,2)*P(33)+2*((EE0I+EijI(1,1))*P(3)+(EijI(5,2)-EijI(
     -   7,2))*P(50))+EijI(4,1)*P(110)-EijI(2,1)*P(111)+Is34*(C0235I-
     -   Cij235I(1,1)+Cij235I(2,1)-4*(Dij1235I(7,2)-Dij1235I(12,3)+Di
     -   j1235I(13,3))+2*(Cij135I(1,2)+(D01235I+Dij1235I(1,1))*P(3))-
     -   Dij1235I(2,2)*P(5)+Is12*(4*Cij135I(4,2)+C0135I*P(26)+Cij135I
     -   (1,1)*P(27)-Cij135I(2,1)*P(28)+2*(B035I-Cij135I(3,2)*P(29)))
     -   +(Dij1235I(2,1)+2*(Dij1235I(4,2)-Dij1235I(5,2)))*P(33)-Dij12
     -   35I(3,1)*P(34)+Dij1235I(6,2)*P(42)-Dij1235I(3,2)*P(142))+Is1
     -   2*(3*(C0345I-Cij145I(2,1))+Cij345I(2,1)-8*Dij1345I(7,2)+4*(D
     -   ij1345I(11,3)-Dij1345I(13,3))-Dij1345I(2,1)*P(16)+Is45*(2*B0
     -   12I+4*Cij145I(4,2)+C0145I*P(95)+Cij145I(1,1)*P(96))+D01345I*
     -   P(97)+Dij1345I(3,2)*P(101)+2*(Cij135I(2,2)+Cij145I(1,2)-Cij1
     -   45I(3,2)-p3sq*Dij1345I(2,2)+Dij1345I(4,2)*P(198))+Dij1345I(1
     -   ,1)*P(497)-Dij1345I(3,1)*P(498)-Dij1345I(1,2)*P(499)-Dij1345
     -   I(5,2)*P(500)-Dij1345I(6,2)*P(501))+EijI(9,2)*P(502))
       F(119)=DCMPLX(FR(119),FI(119))
       FR(120) = -4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+
     -   2*(Dij1235R(2,3)-Dij1235R(8,3))-s34*(-EijR(2,1)+3*EijR(2,3)+
     -   EijR(3,1)-EijR(3,2)+EijR(8,2)-4*(EijR(9,3)+EijR(9,4))+EijR(1
     -   2,3)+2*(EijR(2,4)+EijR(20,4))))
       FI(120) = -4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+
     -   2*(Dij1235I(2,3)-Dij1235I(8,3))-s34*(-EijI(2,1)+3*EijI(2,3)+
     -   EijI(3,1)-EijI(3,2)+EijI(8,2)-4*(EijI(9,3)+EijI(9,4))+EijI(1
     -   2,3)+2*(EijI(2,4)+EijI(20,4))))
       F(120)=DCMPLX(FR(120),FI(120))
       FR(121) = 4*Is34*(D01235R+Dij1235R(1,1)-Dij1235R(2,2)+3*Dij1235
     -   R(4,2)+2*(Dij1235R(2,1)-Dij1235R(2,3)+Dij1235R(6,3))+s34*(-E
     -   ijR(2,1)+EijR(2,3)+EijR(3,1)-EijR(6,2)+4*EijR(8,2)-3*(EijR(2
     -   ,2)+EijR(8,3))+EijR(9,3)+EijR(17,3)+2*(EijR(2,4)-EijR(8,4)-E
     -   ijR(9,4)+EijR(26,4))))
       FI(121) = 4*Is34*(D01235I+Dij1235I(1,1)-Dij1235I(2,2)+3*Dij1235
     -   I(4,2)+2*(Dij1235I(2,1)-Dij1235I(2,3)+Dij1235I(6,3))+s34*(-E
     -   ijI(2,1)+EijI(2,3)+EijI(3,1)-EijI(6,2)+4*EijI(8,2)-3*(EijI(2
     -   ,2)+EijI(8,3))+EijI(9,3)+EijI(17,3)+2*(EijI(2,4)-EijI(8,4)-E
     -   ijI(9,4)+EijI(26,4))))
       F(121)=DCMPLX(FR(121),FI(121))
       FR(122) = -4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+
     -   2*(Dij1235R(2,3)-Dij1235R(8,3))+s34*(EijR(2,1)+EijR(3,1)-Eij
     -   R(4,1)-EijR(9,2)-3*(EijR(2,3)-EijR(9,3))+EijR(10,2)+EijR(10,
     -   3)-EijR(20,3)-2*(EijR(2,4)-EijR(8,2)-EijR(9,4)-EijR(10,4)+Ei
     -   jR(28,4))))
       FI(122) = -4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+
     -   2*(Dij1235I(2,3)-Dij1235I(8,3))+s34*(EijI(2,1)+EijI(3,1)-Eij
     -   I(4,1)-EijI(9,2)-3*(EijI(2,3)-EijI(9,3))+EijI(10,2)+EijI(10,
     -   3)-EijI(20,3)-2*(EijI(2,4)-EijI(8,2)-EijI(9,4)-EijI(10,4)+Ei
     -   jI(28,4))))
       F(122)=DCMPLX(FR(122),FI(122))
       FR(123) = 4*(EijR(2,1)+EijR(3,1)-EijR(4,1)+EijR(9,2)+3*(EijR(2,
     -   3)-EijR(9,3))-EijR(10,2)-EijR(10,3)+EijR(20,3)+2*(EijR(2,2)+
     -   EijR(2,4)-EijR(9,4)-EijR(10,4)+EijR(28,4)))
       FI(123) = 4*(EijI(2,1)+EijI(3,1)-EijI(4,1)+EijI(9,2)+3*(EijI(2,
     -   3)-EijI(9,3))-EijI(10,2)-EijI(10,3)+EijI(20,3)+2*(EijI(2,2)+
     -   EijI(2,4)-EijI(9,4)-EijI(10,4)+EijI(28,4)))
       F(123)=DCMPLX(FR(123),FI(123))
       P(503) = p2sq+s12+s15
       P(504) = -p2sq+p5sq+s12-s15
       P(505) = p2sq+p3sq+2*s12-s45
       FR(124) = -4*(D02345R-Dij2345R(1,1)+Dij2345R(2,1)-p5sq*EijR(2,1
     -   )+s45*EijR(3,1)+2*(-3*EijR(11,2)+EijR(23,3)-EijR(37,4)+EijR(
     -   43,4))+EijR(4,1)*P(7)-EijR(6,2)*P(33)+(-EijR(9,2)+EijR(10,2)
     -   )*P(101)-EijR(3,2)*P(198)-EijR(2,2)*P(503)-EijR(5,2)*P(504)+
     -   EijR(8,2)*P(505))
       FI(124) = -4*(D02345I-Dij2345I(1,1)+Dij2345I(2,1)-p5sq*EijI(2,1
     -   )+s45*EijI(3,1)+2*(-3*EijI(11,2)+EijI(23,3)-EijI(37,4)+EijI(
     -   43,4))+EijI(4,1)*P(7)-EijI(6,2)*P(33)+(-EijI(9,2)+EijI(10,2)
     -   )*P(101)-EijI(3,2)*P(198)-EijI(2,2)*P(503)-EijI(5,2)*P(504)+
     -   EijI(8,2)*P(505))
       F(124)=DCMPLX(FR(124),FI(124))
       P(506) = p3sq+p5sq+s12-s34
       P(507) = -3*s45+2*P(506)
       P(508) = 3*p2sq+s12
       P(509) = p2sq-s12-s23+3*s45
       P(510) = p2sq-s12+s23-s45
       P(511) = 2*p3sq-s23-s45+3*P(5)
       P(512) = p5sq-s15+s23-s45+2*P(14)
       FR(125) = Is12*(-4*Is45*(B012R+C0145R*p5sq+2*Cij145R(4,2)+Cij14
     -   5R(1,1)*P(129))-2*(-10*Dij1345R(7,2)+4*(C0345R+Dij1345R(11,3
     -   ))+Dij1345R(3,1)*P(7)+Dij1345R(6,2)*P(14)+D01345R*P(32)-2*(C
     -   ij135R(2,1)+Cij135R(3,2)-Cij145R(1,2)+p3sq*Dij1345R(2,2)+Dij
     -   1345R(2,1)*P(16)-Dij1345R(4,2)*P(198))-Dij1345R(5,2)*P(292)-
     -   Dij1345R(1,2)*P(499)+Dij1345R(1,1)*P(507)))-2*(-Dij2345R(1,1
     -   )+3*(D02345R+EijR(4,1)*P(7))-EijR(2,1)*P(111)+(EE0R+EijR(1,1
     -   ))*P(127)+Is34*(-Cij125R(2,1)-Cij235R(1,1)+p2sq*Dij1235R(2,1
     -   )-6*Dij1235R(7,2)+4*Dij1235R(12,3)-Dij1235R(2,2)*P(5)+2*(C02
     -   35R-Cij125R(2,2)+Cij135R(1,2)+Dij1235R(4,2)*P(33))-Dij1235R(
     -   6,2)*P(38)+(D01235R+Dij1235R(1,1))*P(133)+Dij1235R(3,1)*P(13
     -   8)+Is12*(2*B035R+4*Cij135R(4,2)+Cij135R(1,1)*P(130)-C0135R*P
     -   (173)))-EijR(2,2)*P(508)+EijR(3,1)*P(509)-2*(-Dij2345R(2,1)+
     -   7*EijR(11,2)+EijR(22,3)+EijR(23,3)+2*(EijR(37,4)-EijR(43,4))
     -   +EijR(6,2)*P(33)+EijR(3,2)*P(198)+EijR(10,2)*P(312)-EijR(5,2
     -   )*P(510))+EijR(8,2)*P(511)-EijR(9,2)*P(512))
       FI(125) = Is12*(-4*Is45*(B012I+C0145I*p5sq+2*Cij145I(4,2)+Cij14
     -   5I(1,1)*P(129))-2*(-10*Dij1345I(7,2)+4*(C0345I+Dij1345I(11,3
     -   ))+Dij1345I(3,1)*P(7)+Dij1345I(6,2)*P(14)+D01345I*P(32)-2*(C
     -   ij135I(2,1)+Cij135I(3,2)-Cij145I(1,2)+p3sq*Dij1345I(2,2)+Dij
     -   1345I(2,1)*P(16)-Dij1345I(4,2)*P(198))-Dij1345I(5,2)*P(292)-
     -   Dij1345I(1,2)*P(499)+Dij1345I(1,1)*P(507)))-2*(-Dij2345I(1,1
     -   )+3*(D02345I+EijI(4,1)*P(7))-EijI(2,1)*P(111)+(EE0I+EijI(1,1
     -   ))*P(127)+Is34*(-Cij125I(2,1)-Cij235I(1,1)+p2sq*Dij1235I(2,1
     -   )-6*Dij1235I(7,2)+4*Dij1235I(12,3)-Dij1235I(2,2)*P(5)+2*(C02
     -   35I-Cij125I(2,2)+Cij135I(1,2)+Dij1235I(4,2)*P(33))-Dij1235I(
     -   6,2)*P(38)+(D01235I+Dij1235I(1,1))*P(133)+Dij1235I(3,1)*P(13
     -   8)+Is12*(2*B035I+4*Cij135I(4,2)+Cij135I(1,1)*P(130)-C0135I*P
     -   (173)))-EijI(2,2)*P(508)+EijI(3,1)*P(509)-2*(-Dij2345I(2,1)+
     -   7*EijI(11,2)+EijI(22,3)+EijI(23,3)+2*(EijI(37,4)-EijI(43,4))
     -   +EijI(6,2)*P(33)+EijI(3,2)*P(198)+EijI(10,2)*P(312)-EijI(5,2
     -   )*P(510))+EijI(8,2)*P(511)-EijI(9,2)*P(512))
       F(125)=DCMPLX(FR(125),FI(125))
       P(513) = s12-2*P(39)
       P(514) = s12+2*P(39)
       P(515) = 2*s34-3*s45
       P(516) = p2sq-s15+3*P(137)
       P(517) = p2sq+3*s12
       P(518) = p2sq-2*p3sq-s23-3*P(16)
       P(519) = p3sq-s15
       P(520) = 6*s12+s23-3*s45+2*P(519)
       FR(126) = Is34*(Is12*(4*B035R+8*Cij135R(4,2)-2*(C0135R*P(513)-C
     -   ij135R(1,1)*P(514)))+2*(3*C0235R+Cij235R(1,1)-8*Dij1235R(7,2
     -   )-4*Dij1235R(12,3)+2*(D01235R+Dij1235R(1,1))*P(3)-Dij1235R(2
     -   ,2)*P(33)-Dij1235R(2,1)*P(60)+Dij1235R(6,2)*P(330)+Dij1235R(
     -   3,1)*P(516)))+2*(-Dij2345R(1,2)-3*(Dij2345R(1,1)-Dij2345R(2,
     -   1))+Dij2345R(4,2)-4*(3*EijR(11,2)-2*EijR(22,3)+EijR(23,3)-Ei
     -   jR(37,4)+EijR(43,4))+EijR(12,3)*P(10)+(-EijR(10,3)+EijR(20,3
     -   ))*P(21)+(EijR(2,3)-2*EijR(6,2))*P(33)-EijR(9,3)*P(245)+2*(D
     -   02345R-EijR(2,1)*P(3)+EijR(4,1)*P(7)-EijR(3,1)*P(180)+EijR(5
     -   ,2)*P(296))+(-EijR(9,2)+EijR(10,2))*P(339)+Is12*(3*C0345R+Ci
     -   j135R(2,1)+Dij1345R(6,2)*P(14)-2*(p3sq*Dij1345R(2,2)+8*Dij13
     -   45R(7,2)+2*Dij1345R(11,3)-Dij1345R(3,1)*P(7)+Dij1345R(1,1)*P
     -   (39))-Dij1345R(1,2)*P(482)+Dij1345R(4,2)*P(483)-Dij1345R(5,2
     -   )*P(491)-Dij1345R(2,1)*P(515))-EijR(2,2)*P(517)+EijR(3,2)*P(
     -   518)+EijR(8,2)*P(520))
       FI(126) = Is34*(Is12*(4*B035I+8*Cij135I(4,2)-2*(C0135I*P(513)-C
     -   ij135I(1,1)*P(514)))+2*(3*C0235I+Cij235I(1,1)-8*Dij1235I(7,2
     -   )-4*Dij1235I(12,3)+2*(D01235I+Dij1235I(1,1))*P(3)-Dij1235I(2
     -   ,2)*P(33)-Dij1235I(2,1)*P(60)+Dij1235I(6,2)*P(330)+Dij1235I(
     -   3,1)*P(516)))+2*(-Dij2345I(1,2)-3*(Dij2345I(1,1)-Dij2345I(2,
     -   1))+Dij2345I(4,2)-4*(3*EijI(11,2)-2*EijI(22,3)+EijI(23,3)-Ei
     -   jI(37,4)+EijI(43,4))+EijI(12,3)*P(10)+(-EijI(10,3)+EijI(20,3
     -   ))*P(21)+(EijI(2,3)-2*EijI(6,2))*P(33)-EijI(9,3)*P(245)+2*(D
     -   02345I-EijI(2,1)*P(3)+EijI(4,1)*P(7)-EijI(3,1)*P(180)+EijI(5
     -   ,2)*P(296))+(-EijI(9,2)+EijI(10,2))*P(339)+Is12*(3*C0345I+Ci
     -   j135I(2,1)+Dij1345I(6,2)*P(14)-2*(p3sq*Dij1345I(2,2)+8*Dij13
     -   45I(7,2)+2*Dij1345I(11,3)-Dij1345I(3,1)*P(7)+Dij1345I(1,1)*P
     -   (39))-Dij1345I(1,2)*P(482)+Dij1345I(4,2)*P(483)-Dij1345I(5,2
     -   )*P(491)-Dij1345I(2,1)*P(515))-EijI(2,2)*P(517)+EijI(3,2)*P(
     -   518)+EijI(8,2)*P(520))
       F(126)=DCMPLX(FR(126),FI(126))
       FR(127) = -4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)+2*Dij1235R(2,3
     -   )+s34*(EijR(2,1)-EijR(3,1)-3*(EijR(2,3)-EijR(9,3))-2*(EijR(2
     -   ,4)-EijR(9,4))))
       FI(127) = -4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)+2*Dij1235I(2,3
     -   )+s34*(EijI(2,1)-EijI(3,1)-3*(EijI(2,3)-EijI(9,3))-2*(EijI(2
     -   ,4)-EijI(9,4))))
       F(127)=DCMPLX(FR(127),FI(127))
       P(521) = p2sq+s12-s34
       P(522) = p2sq+2*s12-s34-s45
       FR(128) = -4*(Dij2345R(1,1)-Dij2345R(2,1)+2*(EijR(11,2)-2*EijR(
     -   22,3)+EijR(23,3)-EijR(37,4)+EijR(43,4))+(EijR(2,1)-EijR(3,1)
     -   )*P(6)+(EijR(9,2)-EijR(10,2))*P(7)+EijR(3,2)*P(16)+(-EijR(5,
     -   2)+EijR(6,2))*P(33)+EijR(2,2)*P(521)-EijR(8,2)*P(522))
       FI(128) = -4*(Dij2345I(1,1)-Dij2345I(2,1)+2*(EijI(11,2)-2*EijI(
     -   22,3)+EijI(23,3)-EijI(37,4)+EijI(43,4))+(EijI(2,1)-EijI(3,1)
     -   )*P(6)+(EijI(9,2)-EijI(10,2))*P(7)+EijI(3,2)*P(16)+(-EijI(5,
     -   2)+EijI(6,2))*P(33)+EijI(2,2)*P(521)-EijI(8,2)*P(522))
       F(128)=DCMPLX(FR(128),FI(128))
       FR(129) = 4*(s34*(EijR(2,1)+EijR(2,2)-EijR(4,1)-EijR(9,2))+p3sq
     -   *(-EijR(3,1)+EijR(4,1)-EijR(8,2)+EijR(9,2))+2*(EijR(11,2)+2*
     -   EijR(22,3)+EijR(37,4)))
       FI(129) = 4*(s34*(EijI(2,1)+EijI(2,2)-EijI(4,1)-EijI(9,2))+p3sq
     -   *(-EijI(3,1)+EijI(4,1)-EijI(8,2)+EijI(9,2))+2*(EijI(11,2)+2*
     -   EijI(22,3)+EijI(37,4)))
       F(129)=DCMPLX(FR(129),FI(129))
       FR(130) = -8*(EijR(2,3)-EijR(2,4)+EijR(5,3)+EijR(5,4)-2*EijR(8,
     -   3)+3*(EijR(8,4)-EijR(17,4)))
       FI(130) = -8*(EijI(2,3)-EijI(2,4)+EijI(5,3)+EijI(5,4)-2*EijI(8,
     -   3)+3*(EijI(8,4)-EijI(17,4)))
       F(130)=DCMPLX(FR(130),FI(130))
       FR(131) = 8*(EijR(2,4)+EijR(6,3)+EijR(9,3)-EijR(9,4)+EijR(17,4)
     -   -EijR(23,4)-2*(EijR(8,4)+EijR(17,3)-EijR(26,4)))
       FI(131) = 8*(EijI(2,4)+EijI(6,3)+EijI(9,3)-EijI(9,4)+EijI(17,4)
     -   -EijI(23,4)-2*(EijI(8,4)+EijI(17,3)-EijI(26,4)))
       F(131)=DCMPLX(FR(131),FI(131))
       FR(132) = 8*(EijR(2,4)+EijR(7,3)+EijR(10,3)-EijR(10,4)+EijR(17,
     -   4)-EijR(24,4)-2*(EijR(8,4)+EijR(18,3)-EijR(27,4)))
       FI(132) = 8*(EijI(2,4)+EijI(7,3)+EijI(10,3)-EijI(10,4)+EijI(17,
     -   4)-EijI(24,4)-2*(EijI(8,4)+EijI(18,3)-EijI(27,4)))
       F(132)=DCMPLX(FR(132),FI(132))
       FR(133) = 8*(EijR(2,4)-2*EijR(8,4)+EijR(17,4))
       FI(133) = 8*(EijI(2,4)-2*EijI(8,4)+EijI(17,4))
       F(133)=DCMPLX(FR(133),FI(133))
       FR(134) = -8*Is34*(Dij1235R(2,2)-Dij1235R(2,3)-Dij1235R(4,2)-Di
     -   j1235R(4,3)+2*Dij1235R(6,3)+s34*(EijR(2,2)+EijR(2,3)-EijR(2,
     -   4)-EijR(5,2)-EijR(8,3)-EijR(9,3)+EijR(9,4)+EijR(17,3)-EijR(1
     -   7,4)+EijR(23,4)+2*(EijR(8,4)-EijR(26,4))))
       FI(134) = -8*Is34*(Dij1235I(2,2)-Dij1235I(2,3)-Dij1235I(4,2)-Di
     -   j1235I(4,3)+2*Dij1235I(6,3)+s34*(EijI(2,2)+EijI(2,3)-EijI(2,
     -   4)-EijI(5,2)-EijI(8,3)-EijI(9,3)+EijI(9,4)+EijI(17,3)-EijI(1
     -   7,4)+EijI(23,4)+2*(EijI(8,4)-EijI(26,4))))
       F(134)=DCMPLX(FR(134),FI(134))
       FR(135) = 8*Is34*(Dij1235R(2,3)-Dij1235R(5,2)+Dij1235R(6,2)-Dij
     -   1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3)+s34*(-EijR(2,2)+EijR
     -   (2,4)+EijR(5,2)-EijR(6,2)+EijR(8,2)-EijR(8,4)+EijR(9,3)+EijR
     -   (11,3)-EijR(12,3)-EijR(17,3)+EijR(20,4)-2*(EijR(9,4)-EijR(26
     -   ,4))-EijR(29,4)))
       FI(135) = 8*Is34*(Dij1235I(2,3)-Dij1235I(5,2)+Dij1235I(6,2)-Dij
     -   1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3)+s34*(-EijI(2,2)+EijI
     -   (2,4)+EijI(5,2)-EijI(6,2)+EijI(8,2)-EijI(8,4)+EijI(9,3)+EijI
     -   (11,3)-EijI(12,3)-EijI(17,3)+EijI(20,4)-2*(EijI(9,4)-EijI(26
     -   ,4))-EijI(29,4)))
       F(135)=DCMPLX(FR(135),FI(135))
       FR(136) = 8*Is34*(Dij1235R(2,3)-Dij1235R(5,2)+Dij1235R(6,2)-Dij
     -   1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3)+s34*(-EijR(2,2)+EijR
     -   (2,4)+EijR(5,2)-EijR(8,4)-2*(EijR(7,2)-EijR(9,2))-EijR(9,4)+
     -   EijR(10,3)-EijR(10,4)-EijR(18,3)+EijR(19,3)-EijR(20,3)+EijR(
     -   26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4)))
       FI(136) = 8*Is34*(Dij1235I(2,3)-Dij1235I(5,2)+Dij1235I(6,2)-Dij
     -   1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3)+s34*(-EijI(2,2)+EijI
     -   (2,4)+EijI(5,2)-EijI(8,4)-2*(EijI(7,2)-EijI(9,2))-EijI(9,4)+
     -   EijI(10,3)-EijI(10,4)-EijI(18,3)+EijI(19,3)-EijI(20,3)+EijI(
     -   26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4)))
       F(136)=DCMPLX(FR(136),FI(136))
       FR(137) = -8*Is34*(-Dij1235R(2,3)+Dij1235R(6,3)+s34*(EijR(2,2)-
     -   EijR(2,4)-EijR(5,2)+EijR(8,4)+EijR(9,4)-EijR(26,4)))
       FI(137) = -8*Is34*(-Dij1235I(2,3)+Dij1235I(6,3)+s34*(EijI(2,2)-
     -   EijI(2,4)-EijI(5,2)+EijI(8,4)+EijI(9,4)-EijI(26,4)))
       F(137)=DCMPLX(FR(137),FI(137))
       FR(138) = -8*Is34*(Dij1235R(2,2)-Dij1235R(2,3)-Dij1235R(4,2)-Di
     -   j1235R(4,3)+2*Dij1235R(6,3)+s34*(EijR(2,2)-EijR(2,4)-EijR(5,
     -   2)-EijR(5,3)+EijR(8,3)-EijR(10,3)+EijR(10,4)-EijR(17,4)+EijR
     -   (18,3)+EijR(24,4)+2*(EijR(8,4)-EijR(27,4))))
       FI(138) = -8*Is34*(Dij1235I(2,2)-Dij1235I(2,3)-Dij1235I(4,2)-Di
     -   j1235I(4,3)+2*Dij1235I(6,3)+s34*(EijI(2,2)-EijI(2,4)-EijI(5,
     -   2)-EijI(5,3)+EijI(8,3)-EijI(10,3)+EijI(10,4)-EijI(17,4)+EijI
     -   (18,3)+EijI(24,4)+2*(EijI(8,4)-EijI(27,4))))
       F(138)=DCMPLX(FR(138),FI(138))
       FR(139) = 8*Is34*(Dij1235R(2,3)-Dij1235R(5,2)+Dij1235R(6,2)-Dij
     -   1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3)+s34*(EijR(2,3)+EijR(
     -   2,4)-EijR(6,2)+EijR(7,2)+EijR(8,2)-EijR(8,3)-EijR(8,4)-EijR(
     -   9,2)-EijR(9,4)-EijR(10,4)+EijR(19,3)-EijR(20,3)+EijR(26,4)+E
     -   ijR(27,4)+EijR(28,4)-EijR(35,4)))
       FI(139) = 8*Is34*(Dij1235I(2,3)-Dij1235I(5,2)+Dij1235I(6,2)-Dij
     -   1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3)+s34*(EijI(2,3)+EijI(
     -   2,4)-EijI(6,2)+EijI(7,2)+EijI(8,2)-EijI(8,3)-EijI(8,4)-EijI(
     -   9,2)-EijI(9,4)-EijI(10,4)+EijI(19,3)-EijI(20,3)+EijI(26,4)+E
     -   ijI(27,4)+EijI(28,4)-EijI(35,4)))
       F(139)=DCMPLX(FR(139),FI(139))
       FR(140) = 8*Is34*(Dij1235R(2,3)-Dij1235R(5,2)+Dij1235R(6,2)-Dij
     -   1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3)+s34*(EijR(2,3)+EijR(
     -   2,4)-EijR(7,2)-EijR(8,3)-EijR(8,4)+EijR(9,2)+EijR(14,3)-EijR
     -   (15,3)+EijR(21,4)-2*(EijR(10,4)-EijR(27,4))-EijR(32,4)))
       FI(140) = 8*Is34*(Dij1235I(2,3)-Dij1235I(5,2)+Dij1235I(6,2)-Dij
     -   1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3)+s34*(EijI(2,3)+EijI(
     -   2,4)-EijI(7,2)-EijI(8,3)-EijI(8,4)+EijI(9,2)+EijI(14,3)-EijI
     -   (15,3)+EijI(21,4)-2*(EijI(10,4)-EijI(27,4))-EijI(32,4)))
       F(140)=DCMPLX(FR(140),FI(140))
       FR(141) = 8*Is34*(Dij1235R(2,3)-Dij1235R(6,3)+s34*(EijR(2,3)+Ei
     -   jR(2,4)-EijR(8,3)-EijR(8,4)-EijR(10,4)+EijR(27,4)))
       FI(141) = 8*Is34*(Dij1235I(2,3)-Dij1235I(6,3)+s34*(EijI(2,3)+Ei
     -   jI(2,4)-EijI(8,3)-EijI(8,4)-EijI(10,4)+EijI(27,4)))
       F(141)=DCMPLX(FR(141),FI(141))
       FR(142) = -8*(EijR(2,2)-EijR(2,4)-EijR(5,2)-EijR(5,3)+EijR(8,3)
     -   +2*EijR(8,4)-EijR(17,4))
       FI(142) = -8*(EijI(2,2)-EijI(2,4)-EijI(5,2)-EijI(5,3)+EijI(8,3)
     -   +2*EijI(8,4)-EijI(17,4))
       F(142)=DCMPLX(FR(142),FI(142))
       FR(143) = 8*(EijR(2,3)+EijR(2,4)-EijR(6,2)+EijR(8,2)-EijR(8,3)-
     -   EijR(8,4)-EijR(9,4)+EijR(26,4))
       FI(143) = 8*(EijI(2,3)+EijI(2,4)-EijI(6,2)+EijI(8,2)-EijI(8,3)-
     -   EijI(8,4)-EijI(9,4)+EijI(26,4))
       F(143)=DCMPLX(FR(143),FI(143))
       FR(144) = 8*(EijR(2,3)+EijR(2,4)-EijR(7,2)-EijR(8,3)-EijR(8,4)+
     -   EijR(9,2)-EijR(10,4)+EijR(27,4))
       FI(144) = 8*(EijI(2,3)+EijI(2,4)-EijI(7,2)-EijI(8,3)-EijI(8,4)+
     -   EijI(9,2)-EijI(10,4)+EijI(27,4))
       F(144)=DCMPLX(FR(144),FI(144))
       FR(145) = 8*(EijR(2,3)+EijR(2,4)-EijR(8,3)-EijR(8,4))
       FI(145) = 8*(EijI(2,3)+EijI(2,4)-EijI(8,3)-EijI(8,4))
       F(145)=DCMPLX(FR(145),FI(145))
       FR(146) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)-EijR(5,2)+EijR(6,2)
     -   -EijR(8,2)+EijR(9,3)+EijR(9,4)+3*(EijR(8,3)-EijR(17,3))-EijR
     -   (17,4)+EijR(23,4)-2*(EijR(5,3)-EijR(6,3)-EijR(8,4)+EijR(26,4
     -   )))
       FI(146) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)-EijI(5,2)+EijI(6,2)
     -   -EijI(8,2)+EijI(9,3)+EijI(9,4)+3*(EijI(8,3)-EijI(17,3))-EijI
     -   (17,4)+EijI(23,4)-2*(EijI(5,3)-EijI(6,3)-EijI(8,4)+EijI(26,4
     -   )))
       F(146)=DCMPLX(FR(146),FI(146))
       FR(147) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-
     -   EijR(8,4)-EijR(11,3)+EijR(12,3)-3*(EijR(9,3)-EijR(17,3))+Eij
     -   R(20,4)+2*(EijR(2,3)-EijR(8,3)-EijR(9,4)+EijR(26,4))-EijR(29
     -   ,4))
       FI(147) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-
     -   EijI(8,4)-EijI(11,3)+EijI(12,3)-3*(EijI(9,3)-EijI(17,3))+Eij
     -   I(20,4)+2*(EijI(2,3)-EijI(8,3)-EijI(9,4)+EijI(26,4))-EijI(29
     -   ,4))
       F(147)=DCMPLX(FR(147),FI(147))
       FR(148) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(7,2)-EijR(8,4)-
     -   EijR(9,2)-EijR(9,4)-EijR(10,3)-EijR(10,4)+2*(EijR(2,3)-EijR(
     -   8,3)-EijR(9,3)+EijR(17,3))+EijR(18,3)-EijR(19,3)+EijR(20,3)+
     -   EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4))
       FI(148) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(7,2)-EijI(8,4)-
     -   EijI(9,2)-EijI(9,4)-EijI(10,3)-EijI(10,4)+2*(EijI(2,3)-EijI(
     -   8,3)-EijI(9,3)+EijI(17,3))+EijI(18,3)-EijI(19,3)+EijI(20,3)+
     -   EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4))
       F(148)=DCMPLX(FR(148),FI(148))
       FR(149) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-
     -   EijR(8,4)-EijR(9,4)+2*(EijR(2,3)-EijR(8,3)-EijR(9,3)+EijR(17
     -   ,3))+EijR(26,4))
       FI(149) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-
     -   EijI(8,4)-EijI(9,4)+2*(EijI(2,3)-EijI(8,3)-EijI(9,3)+EijI(17
     -   ,3))+EijI(26,4))
       F(149)=DCMPLX(FR(149),FI(149))
       FR(150) = -8*Is34*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(2,3)-Di
     -   j1235R(3,1)+2*(Dij1235R(4,2)-Dij1235R(5,2))+Dij1235R(6,2)+Di
     -   j1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)+s34*(EijR(2,1)+EijR
     -   (2,2)-EijR(2,3)-EijR(2,4)-EijR(3,1)+EijR(3,2)+EijR(5,2)-EijR
     -   (6,2)+EijR(8,4)-EijR(12,3)-4*EijR(17,3)-EijR(20,4)-2*(EijR(8
     -   ,2)-EijR(8,3)-EijR(9,3)-EijR(9,4)-EijR(11,3)+EijR(26,4))+Eij
     -   R(29,4)))
       FI(150) = -8*Is34*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(2,3)-Di
     -   j1235I(3,1)+2*(Dij1235I(4,2)-Dij1235I(5,2))+Dij1235I(6,2)+Di
     -   j1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)+s34*(EijI(2,1)+EijI
     -   (2,2)-EijI(2,3)-EijI(2,4)-EijI(3,1)+EijI(3,2)+EijI(5,2)-EijI
     -   (6,2)+EijI(8,4)-EijI(12,3)-4*EijI(17,3)-EijI(20,4)-2*(EijI(8
     -   ,2)-EijI(8,3)-EijI(9,3)-EijI(9,4)-EijI(11,3)+EijI(26,4))+Eij
     -   I(29,4)))
       F(150)=DCMPLX(FR(150),FI(150))
       FR(151) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,2)+EijR(2,4)+EijR(3,2)-EijR(3,3)+2*
     -   (EijR(2,3)-EijR(8,2))-5*EijR(9,3)+4*EijR(12,3)-EijR(12,4)-3*
     -   (EijR(9,4)-EijR(20,4))))
       FI(151) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,2)+EijI(2,4)+EijI(3,2)-EijI(3,3)+2*
     -   (EijI(2,3)-EijI(8,2))-5*EijI(9,3)+4*EijI(12,3)-EijI(12,4)-3*
     -   (EijI(9,4)-EijI(20,4))))
       F(151)=DCMPLX(FR(151),FI(151))
       FR(152) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,2)+EijR(2,4)-EijR(8,2)-4*EijR(9,3)-
     -   EijR(10,3)-EijR(10,4)-EijR(13,3)+EijR(20,4)+2*(EijR(2,3)-Eij
     -   R(9,4)+EijR(12,3)+EijR(20,3)+EijR(28,4))-EijR(31,4)))
       FI(152) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,2)+EijI(2,4)-EijI(8,2)-4*EijI(9,3)-
     -   EijI(10,3)-EijI(10,4)-EijI(13,3)+EijI(20,4)+2*(EijI(2,3)-Eij
     -   I(9,4)+EijI(12,3)+EijI(20,3)+EijI(28,4))-EijI(31,4)))
       F(152)=DCMPLX(FR(152),FI(152))
       FR(153) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+2*(
     -   Dij1235R(2,2)-Dij1235R(6,2))-Dij1235R(8,3)+s34*(EijR(2,2)+Ei
     -   jR(2,4)+EijR(3,2)+2*(EijR(2,3)-EijR(8,2)-2*EijR(9,3)-EijR(9,
     -   4)+EijR(12,3))+EijR(20,4)))
       FI(153) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+2*(
     -   Dij1235I(2,2)-Dij1235I(6,2))-Dij1235I(8,3)+s34*(EijI(2,2)+Ei
     -   jI(2,4)+EijI(3,2)+2*(EijI(2,3)-EijI(8,2)-2*EijI(9,3)-EijI(9,
     -   4)+EijI(12,3))+EijI(20,4)))
       F(153)=DCMPLX(FR(153),FI(153))
       FR(154) = -8*Is34*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(2,3)-Di
     -   j1235R(3,1)+2*(Dij1235R(4,2)-Dij1235R(5,2))+Dij1235R(6,2)+Di
     -   j1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)+s34*(EijR(2,1)+EijR
     -   (2,2)-EijR(2,4)-EijR(3,1)+EijR(5,2)+EijR(7,2)+EijR(8,4)+EijR
     -   (9,4)+EijR(10,2)+EijR(10,3)+EijR(10,4)+3*(EijR(8,3)-EijR(17,
     -   3))-2*(EijR(2,3)+EijR(6,2)+EijR(9,2)-EijR(9,3)+EijR(18,3)-Ei
     -   jR(19,3))-EijR(20,3)-EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(3
     -   5,4)))
       FI(154) = -8*Is34*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(2,3)-Di
     -   j1235I(3,1)+2*(Dij1235I(4,2)-Dij1235I(5,2))+Dij1235I(6,2)+Di
     -   j1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)+s34*(EijI(2,1)+EijI
     -   (2,2)-EijI(2,4)-EijI(3,1)+EijI(5,2)+EijI(7,2)+EijI(8,4)+EijI
     -   (9,4)+EijI(10,2)+EijI(10,3)+EijI(10,4)+3*(EijI(8,3)-EijI(17,
     -   3))-2*(EijI(2,3)+EijI(6,2)+EijI(9,2)-EijI(9,3)+EijI(18,3)-Ei
     -   jI(19,3))-EijI(20,3)-EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(3
     -   5,4)))
       F(154)=DCMPLX(FR(154),FI(154))
       FR(155) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,4)+EijR(3,2)-EijR(9,2)-5*EijR(9,3)+
     -   EijR(10,2)-EijR(10,4)-EijR(13,3)+3*(EijR(2,3)-EijR(8,2)+EijR
     -   (20,3))+EijR(20,4)+2*(EijR(2,2)-EijR(9,4)-EijR(10,3)+EijR(12
     -   ,3)+EijR(28,4))-EijR(31,4)))
       FI(155) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,4)+EijI(3,2)-EijI(9,2)-5*EijI(9,3)+
     -   EijI(10,2)-EijI(10,4)-EijI(13,3)+3*(EijI(2,3)-EijI(8,2)+EijI
     -   (20,3))+EijI(20,4)+2*(EijI(2,2)-EijI(9,4)-EijI(10,3)+EijI(12
     -   ,3)+EijI(28,4))-EijI(31,4)))
       F(155)=DCMPLX(FR(155),FI(155))
       FR(156) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,4)-EijR(9,2)+3*(EijR(2,3)-EijR(9,3)
     -   )-EijR(9,4)+EijR(10,2)+EijR(15,3)-EijR(16,3)-4*(EijR(10,3)-E
     -   ijR(20,3))+EijR(21,4)+2*(EijR(2,2)-EijR(8,2)-EijR(10,4)+EijR
     -   (28,4))-EijR(34,4)))
       FI(156) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,4)-EijI(9,2)+3*(EijI(2,3)-EijI(9,3)
     -   )-EijI(9,4)+EijI(10,2)+EijI(15,3)-EijI(16,3)-4*(EijI(10,3)-E
     -   ijI(20,3))+EijI(21,4)+2*(EijI(2,2)-EijI(8,2)-EijI(10,4)+EijI
     -   (28,4))-EijI(34,4)))
       F(156)=DCMPLX(FR(156),FI(156))
       FR(157) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+2*(
     -   Dij1235R(2,2)-Dij1235R(6,2))-Dij1235R(8,3)+s34*(EijR(2,4)-Ei
     -   jR(3,1)+EijR(4,1)+3*(EijR(2,3)-EijR(8,2)-EijR(9,3))-EijR(9,4
     -   )+EijR(10,2)-EijR(10,4)+2*(EijR(2,2)-EijR(10,3)+EijR(20,3))+
     -   EijR(28,4)))
       FI(157) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+2*(
     -   Dij1235I(2,2)-Dij1235I(6,2))-Dij1235I(8,3)+s34*(EijI(2,4)-Ei
     -   jI(3,1)+EijI(4,1)+3*(EijI(2,3)-EijI(8,2)-EijI(9,3))-EijI(9,4
     -   )+EijI(10,2)-EijI(10,4)+2*(EijI(2,2)-EijI(10,3)+EijI(20,3))+
     -   EijI(28,4)))
       F(157)=DCMPLX(FR(157),FI(157))
       FR(158) = -8*(EijR(2,1)-EijR(2,4)-EijR(3,1)+EijR(8,4)-2*(EijR(2
     -   ,3)-EijR(5,2)+EijR(6,2)-EijR(9,3))+EijR(9,4)+3*(EijR(8,3)-Ei
     -   jR(17,3))-EijR(26,4))
       FI(158) = -8*(EijI(2,1)-EijI(2,4)-EijI(3,1)+EijI(8,4)-2*(EijI(2
     -   ,3)-EijI(5,2)+EijI(6,2)-EijI(9,3))+EijI(9,4)+3*(EijI(8,3)-Ei
     -   jI(17,3))-EijI(26,4))
       F(158)=DCMPLX(FR(158),FI(158))
       FR(159) = 8*(EijR(2,1)+3*(EijR(2,2)+EijR(2,3))+EijR(2,4)-EijR(3
     -   ,1)+EijR(3,2)-4*EijR(8,2)-5*EijR(9,3)-2*(EijR(9,4)-EijR(12,3
     -   ))+EijR(20,4))
       FI(159) = 8*(EijI(2,1)+3*(EijI(2,2)+EijI(2,3))+EijI(2,4)-EijI(3
     -   ,1)+EijI(3,2)-4*EijI(8,2)-5*EijI(9,3)-2*(EijI(9,4)-EijI(12,3
     -   ))+EijI(20,4))
       F(159)=DCMPLX(FR(159),FI(159))
       FR(160) = 8*(EijR(2,1)+EijR(2,4)-EijR(4,1)+3*(EijR(2,2)+EijR(2,
     -   3)-EijR(9,3))-EijR(9,4)+EijR(10,2)-EijR(10,4)-2*(EijR(8,2)+E
     -   ijR(9,2)+EijR(10,3)-EijR(20,3))+EijR(28,4))
       FI(160) = 8*(EijI(2,1)+EijI(2,4)-EijI(4,1)+3*(EijI(2,2)+EijI(2,
     -   3)-EijI(9,3))-EijI(9,4)+EijI(10,2)-EijI(10,4)-2*(EijI(8,2)+E
     -   ijI(9,2)+EijI(10,3)-EijI(20,3))+EijI(28,4))
       F(160)=DCMPLX(FR(160),FI(160))
       FR(161) = 8*(EijR(2,1)+EijR(2,4)-EijR(3,1)+3*(EijR(2,2)+EijR(2,
     -   3)-EijR(8,2)-EijR(9,3))-EijR(9,4))
       FI(161) = 8*(EijI(2,1)+EijI(2,4)-EijI(3,1)+3*(EijI(2,2)+EijI(2,
     -   3)-EijI(8,2)-EijI(9,3))-EijI(9,4))
       F(161)=DCMPLX(FR(161),FI(161))
       FR(162) = -8*(EijR(2,2)-EijR(2,3)-EijR(2,4)-EijR(5,2)+EijR(7,2)
     -   -EijR(9,2)+EijR(10,3)+EijR(10,4)-EijR(17,4)+3*(EijR(8,3)-Eij
     -   R(18,3))+EijR(24,4)-2*(EijR(5,3)-EijR(7,3)-EijR(8,4)+EijR(27
     -   ,4)))
       FI(162) = -8*(EijI(2,2)-EijI(2,3)-EijI(2,4)-EijI(5,2)+EijI(7,2)
     -   -EijI(9,2)+EijI(10,3)+EijI(10,4)-EijI(17,4)+3*(EijI(8,3)-Eij
     -   I(18,3))+EijI(24,4)-2*(EijI(5,3)-EijI(7,3)-EijI(8,4)+EijI(27
     -   ,4)))
       F(162)=DCMPLX(FR(162),FI(162))
       FR(163) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(6,2)-EijR(8,2)-
     -   EijR(8,4)-EijR(9,3)-EijR(9,4)-EijR(10,4)+EijR(17,3)+2*(EijR(
     -   2,3)-EijR(8,3)-EijR(10,3)+EijR(18,3))-EijR(19,3)+EijR(20,3)+
     -   EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(35,4))
       FI(163) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(6,2)-EijI(8,2)-
     -   EijI(8,4)-EijI(9,3)-EijI(9,4)-EijI(10,4)+EijI(17,3)+2*(EijI(
     -   2,3)-EijI(8,3)-EijI(10,3)+EijI(18,3))-EijI(19,3)+EijI(20,3)+
     -   EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(35,4))
       F(163)=DCMPLX(FR(163),FI(163))
       FR(164) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(7,2)-EijR(8,4)-
     -   EijR(9,2)-EijR(14,3)+EijR(15,3)-3*(EijR(10,3)-EijR(18,3))+Ei
     -   jR(21,4)+2*(EijR(2,3)-EijR(8,3)-EijR(10,4)+EijR(27,4))-EijR(
     -   32,4))
       FI(164) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(7,2)-EijI(8,4)-
     -   EijI(9,2)-EijI(14,3)+EijI(15,3)-3*(EijI(10,3)-EijI(18,3))+Ei
     -   jI(21,4)+2*(EijI(2,3)-EijI(8,3)-EijI(10,4)+EijI(27,4))-EijI(
     -   32,4))
       F(164)=DCMPLX(FR(164),FI(164))
       FR(165) = 8*(EijR(2,2)+EijR(2,4)-EijR(5,2)+EijR(7,2)-EijR(8,4)-
     -   EijR(9,2)-EijR(10,4)+2*(EijR(2,3)-EijR(8,3)-EijR(10,3)+EijR(
     -   18,3))+EijR(27,4))
       FI(165) = 8*(EijI(2,2)+EijI(2,4)-EijI(5,2)+EijI(7,2)-EijI(8,4)-
     -   EijI(9,2)-EijI(10,4)+2*(EijI(2,3)-EijI(8,3)-EijI(10,3)+EijI(
     -   18,3))+EijI(27,4))
       F(165)=DCMPLX(FR(165),FI(165))
       FR(166) = -8*Is34*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(2,3)-Di
     -   j1235R(3,1)+2*(Dij1235R(4,2)-Dij1235R(5,2))+Dij1235R(6,2)+Di
     -   j1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)+s34*(EijR(2,1)-EijR
     -   (2,3)-EijR(2,4)-EijR(4,1)-EijR(8,2)+EijR(8,4)+EijR(9,3)+EijR
     -   (9,4)+EijR(10,2)+EijR(10,3)+EijR(10,4)+2*(EijR(5,2)-EijR(7,2
     -   )+EijR(8,3)-EijR(17,3)-EijR(18,3)+EijR(19,3))-EijR(20,3)-Eij
     -   R(26,4)-EijR(27,4)-EijR(28,4)+EijR(35,4)))
       FI(166) = -8*Is34*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(2,3)-Di
     -   j1235I(3,1)+2*(Dij1235I(4,2)-Dij1235I(5,2))+Dij1235I(6,2)+Di
     -   j1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)+s34*(EijI(2,1)-EijI
     -   (2,3)-EijI(2,4)-EijI(4,1)-EijI(8,2)+EijI(8,4)+EijI(9,3)+EijI
     -   (9,4)+EijI(10,2)+EijI(10,3)+EijI(10,4)+2*(EijI(5,2)-EijI(7,2
     -   )+EijI(8,3)-EijI(17,3)-EijI(18,3)+EijI(19,3))-EijI(20,3)-Eij
     -   I(26,4)-EijI(27,4)-EijI(28,4)+EijI(35,4)))
       F(166)=DCMPLX(FR(166),FI(166))
       FR(167) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,1)+EijR(2,4)-EijR(3,1)+EijR(3,2)-Ei
     -   jR(9,2)+EijR(10,2)-EijR(10,4)+EijR(12,3)-EijR(13,3)-3*(EijR(
     -   8,2)+EijR(9,3)-EijR(20,3))+EijR(20,4)+2*(EijR(2,2)+EijR(2,3)
     -   -EijR(9,4)-EijR(10,3)+EijR(28,4))-EijR(31,4)))
       FI(167) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,1)+EijI(2,4)-EijI(3,1)+EijI(3,2)-Ei
     -   jI(9,2)+EijI(10,2)-EijI(10,4)+EijI(12,3)-EijI(13,3)-3*(EijI(
     -   8,2)+EijI(9,3)-EijI(20,3))+EijI(20,4)+2*(EijI(2,2)+EijI(2,3)
     -   -EijI(9,4)-EijI(10,3)+EijI(28,4))-EijI(31,4)))
       F(167)=DCMPLX(FR(167),FI(167))
       FR(168) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,1)+EijR(2,4)-EijR(4,1)-EijR(8,2)-Ei
     -   jR(9,4)+EijR(10,2)+EijR(15,3)-EijR(16,3)-3*(EijR(10,3)-EijR(
     -   20,3))+EijR(21,4)+2*(EijR(2,2)+EijR(2,3)-EijR(9,2)-EijR(9,3)
     -   -EijR(10,4)+EijR(28,4))-EijR(34,4)))
       FI(168) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,1)+EijI(2,4)-EijI(4,1)-EijI(8,2)-Ei
     -   jI(9,4)+EijI(10,2)+EijI(15,3)-EijI(16,3)-3*(EijI(10,3)-EijI(
     -   20,3))+EijI(21,4)+2*(EijI(2,2)+EijI(2,3)-EijI(9,2)-EijI(9,3)
     -   -EijI(10,4)+EijI(28,4))-EijI(34,4)))
       F(168)=DCMPLX(FR(168),FI(168))
       FR(169) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+2*(
     -   Dij1235R(2,2)-Dij1235R(6,2))-Dij1235R(8,3)+s34*(EijR(2,1)+Ei
     -   jR(2,4)-EijR(4,1)-EijR(8,2)-EijR(9,4)+EijR(10,2)-EijR(10,4)+
     -   2*(EijR(2,2)+EijR(2,3)-EijR(9,2)-EijR(9,3)-EijR(10,3)+EijR(2
     -   0,3))+EijR(28,4)))
       FI(169) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+2*(
     -   Dij1235I(2,2)-Dij1235I(6,2))-Dij1235I(8,3)+s34*(EijI(2,1)+Ei
     -   jI(2,4)-EijI(4,1)-EijI(8,2)-EijI(9,4)+EijI(10,2)-EijI(10,4)+
     -   2*(EijI(2,2)+EijI(2,3)-EijI(9,2)-EijI(9,3)-EijI(10,3)+EijI(2
     -   0,3))+EijI(28,4)))
       F(169)=DCMPLX(FR(169),FI(169))
       FR(170) = -8*Is34*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(2,3)-Di
     -   j1235R(3,1)+2*(Dij1235R(4,2)-Dij1235R(5,2))+Dij1235R(6,2)+Di
     -   j1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)+s34*(EijR(2,1)-EijR
     -   (2,4)-EijR(4,1)+EijR(4,2)+EijR(8,4)-EijR(9,2)+3*(EijR(8,3)+E
     -   ijR(10,3))-EijR(15,3)-5*EijR(18,3)-EijR(21,4)-2*(EijR(2,3)-E
     -   ijR(5,2)+EijR(7,2)-EijR(10,4)-EijR(14,3)+EijR(27,4))+EijR(32
     -   ,4)))
       FI(170) = -8*Is34*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(2,3)-Di
     -   j1235I(3,1)+2*(Dij1235I(4,2)-Dij1235I(5,2))+Dij1235I(6,2)+Di
     -   j1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)+s34*(EijI(2,1)-EijI
     -   (2,4)-EijI(4,1)+EijI(4,2)+EijI(8,4)-EijI(9,2)+3*(EijI(8,3)+E
     -   ijI(10,3))-EijI(15,3)-5*EijI(18,3)-EijI(21,4)-2*(EijI(2,3)-E
     -   ijI(5,2)+EijI(7,2)-EijI(10,4)-EijI(14,3)+EijI(27,4))+EijI(32
     -   ,4)))
       F(170)=DCMPLX(FR(170),FI(170))
       FR(171) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,1)+EijR(2,4)-EijR(3,1)+EijR(4,2)-4*
     -   EijR(9,2)-EijR(9,4)-5*EijR(10,3)-EijR(16,3)+3*(EijR(2,2)+Eij
     -   R(2,3)+EijR(20,3))+EijR(21,4)-2*(EijR(8,2)+EijR(9,3)-EijR(10
     -   ,2)+EijR(10,4)-EijR(15,3)-EijR(28,4))-EijR(34,4)))
       FI(171) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,1)+EijI(2,4)-EijI(3,1)+EijI(4,2)-4*
     -   EijI(9,2)-EijI(9,4)-5*EijI(10,3)-EijI(16,3)+3*(EijI(2,2)+Eij
     -   I(2,3)+EijI(20,3))+EijI(21,4)-2*(EijI(8,2)+EijI(9,3)-EijI(10
     -   ,2)+EijI(10,4)-EijI(15,3)-EijI(28,4))-EijI(34,4)))
       F(171)=DCMPLX(FR(171),FI(171))
       FR(172) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,1)+EijR(2,4)-EijR(4,1)+2*EijR(4,2)-
     -   EijR(4,3)-7*EijR(10,3)-5*(EijR(9,2)-EijR(15,3))-EijR(15,4)+3
     -   *(EijR(2,2)+EijR(2,3)-EijR(10,4)+EijR(21,4))))
       FI(172) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,1)+EijI(2,4)-EijI(4,1)+2*EijI(4,2)-
     -   EijI(4,3)-7*EijI(10,3)-5*(EijI(9,2)-EijI(15,3))-EijI(15,4)+3
     -   *(EijI(2,2)+EijI(2,3)-EijI(10,4)+EijI(21,4))))
       F(172)=DCMPLX(FR(172),FI(172))
       FR(173) = 8*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+2*(
     -   Dij1235R(2,2)-Dij1235R(6,2))-Dij1235R(8,3)+s34*(EijR(2,1)+3*
     -   (EijR(2,2)+EijR(2,3))+EijR(2,4)-EijR(4,1)+EijR(4,2)-4*EijR(9
     -   ,2)-5*EijR(10,3)-2*(EijR(10,4)-EijR(15,3))+EijR(21,4)))
       FI(173) = 8*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+2*(
     -   Dij1235I(2,2)-Dij1235I(6,2))-Dij1235I(8,3)+s34*(EijI(2,1)+3*
     -   (EijI(2,2)+EijI(2,3))+EijI(2,4)-EijI(4,1)+EijI(4,2)-4*EijI(9
     -   ,2)-5*EijI(10,3)-2*(EijI(10,4)-EijI(15,3))+EijI(21,4)))
       F(173)=DCMPLX(FR(173),FI(173))
       FR(174) = -8*(EijR(2,1)-EijR(2,4)-EijR(4,1)+EijR(8,4)-2*(EijR(2
     -   ,3)-EijR(5,2)+EijR(7,2)-EijR(10,3))+EijR(10,4)+3*(EijR(8,3)-
     -   EijR(18,3))-EijR(27,4))
       FI(174) = -8*(EijI(2,1)-EijI(2,4)-EijI(4,1)+EijI(8,4)-2*(EijI(2
     -   ,3)-EijI(5,2)+EijI(7,2)-EijI(10,3))+EijI(10,4)+3*(EijI(8,3)-
     -   EijI(18,3))-EijI(27,4))
       F(174)=DCMPLX(FR(174),FI(174))
       FR(175) = 8*(EijR(2,1)+EijR(2,4)-EijR(3,1)-EijR(9,4)+EijR(10,2)
     -   +3*(EijR(2,2)+EijR(2,3)-EijR(10,3))-EijR(10,4)-2*(EijR(8,2)+
     -   EijR(9,2)+EijR(9,3)-EijR(20,3))+EijR(28,4))
       FI(175) = 8*(EijI(2,1)+EijI(2,4)-EijI(3,1)-EijI(9,4)+EijI(10,2)
     -   +3*(EijI(2,2)+EijI(2,3)-EijI(10,3))-EijI(10,4)-2*(EijI(8,2)+
     -   EijI(9,2)+EijI(9,3)-EijI(20,3))+EijI(28,4))
       F(175)=DCMPLX(FR(175),FI(175))
       FR(176) = 8*(EijR(2,1)+3*(EijR(2,2)+EijR(2,3))+EijR(2,4)-EijR(4
     -   ,1)+EijR(4,2)-4*EijR(9,2)-5*EijR(10,3)-2*(EijR(10,4)-EijR(15
     -   ,3))+EijR(21,4))
       FI(176) = 8*(EijI(2,1)+3*(EijI(2,2)+EijI(2,3))+EijI(2,4)-EijI(4
     -   ,1)+EijI(4,2)-4*EijI(9,2)-5*EijI(10,3)-2*(EijI(10,4)-EijI(15
     -   ,3))+EijI(21,4))
       F(176)=DCMPLX(FR(176),FI(176))
       FR(177) = 8*(EijR(2,1)+EijR(2,4)-EijR(4,1)+3*(EijR(2,2)+EijR(2,
     -   3)-EijR(9,2)-EijR(10,3))-EijR(10,4))
       FI(177) = 8*(EijI(2,1)+EijI(2,4)-EijI(4,1)+3*(EijI(2,2)+EijI(2,
     -   3)-EijI(9,2)-EijI(10,3))-EijI(10,4))
       F(177)=DCMPLX(FR(177),FI(177))
       FR(178) = 4*(EijR(1,1)+EijR(1,2)-EijR(2,1)+EijR(2,3)+EijR(5,2)+
     -   3*EijR(5,3)-4*(EijR(8,3)+EijR(8,4))-2*(EijR(2,2)-EijR(2,4)-E
     -   ijR(17,4)))
       FI(178) = 4*(EijI(1,1)+EijI(1,2)-EijI(2,1)+EijI(2,3)+EijI(5,2)+
     -   3*EijI(5,3)-4*(EijI(8,3)+EijI(8,4))-2*(EijI(2,2)-EijI(2,4)-E
     -   ijI(17,4)))
       F(178)=DCMPLX(FR(178),FI(178))
       FR(179) = 4*(EijR(2,2)-EijR(5,2)+3*(EijR(2,3)-EijR(8,3))-EijR(9
     -   ,3)+EijR(17,3)+2*(EijR(2,4)-EijR(8,4)-EijR(9,4)+EijR(26,4)))
       FI(179) = 4*(EijI(2,2)-EijI(5,2)+3*(EijI(2,3)-EijI(8,3))-EijI(9
     -   ,3)+EijI(17,3)+2*(EijI(2,4)-EijI(8,4)-EijI(9,4)+EijI(26,4)))
       F(179)=DCMPLX(FR(179),FI(179))
       FR(180) = 4*(EijR(2,2)-EijR(5,2)+3*(EijR(2,3)-EijR(8,3))-EijR(1
     -   0,3)+EijR(18,3)+2*(EijR(2,4)-EijR(8,4)-EijR(10,4)+EijR(27,4)
     -   ))
       FI(180) = 4*(EijI(2,2)-EijI(5,2)+3*(EijI(2,3)-EijI(8,3))-EijI(1
     -   0,3)+EijI(18,3)+2*(EijI(2,4)-EijI(8,4)-EijI(10,4)+EijI(27,4)
     -   ))
       F(180)=DCMPLX(FR(180),FI(180))
       FR(181) = 4*(EijR(2,2)-EijR(5,2)+3*(EijR(2,3)-EijR(8,3))+2*(Eij
     -   R(2,4)-EijR(8,4)))
       FI(181) = 4*(EijI(2,2)-EijI(5,2)+3*(EijI(2,3)-EijI(8,3))+2*(Eij
     -   I(2,4)-EijI(8,4)))
       F(181)=DCMPLX(FR(181),FI(181))
       FR(182) = 4*(s34*(EijR(2,2)-EijR(5,2)+EijR(7,2)-EijR(9,2))+p3sq
     -   *(EijR(6,2)-EijR(7,2)-EijR(8,2)+EijR(9,2))-2*(EijR(11,2)+2*E
     -   ijR(21,3)-EijR(22,3)-EijR(37,4)+EijR(40,4)))
       FI(182) = 4*(s34*(EijI(2,2)-EijI(5,2)+EijI(7,2)-EijI(9,2))+p3sq
     -   *(EijI(6,2)-EijI(7,2)-EijI(8,2)+EijI(9,2))-2*(EijI(11,2)+2*E
     -   ijI(21,3)-EijI(22,3)-EijI(37,4)+EijI(40,4)))
       F(182)=DCMPLX(FR(182),FI(182))
       FR(183) = 4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+2
     -   *(Dij1235R(2,3)-Dij1235R(8,3))+s34*(EijR(2,2)+3*EijR(2,3)-Ei
     -   jR(8,2)-4*(EijR(9,3)+EijR(9,4))+EijR(12,3)+2*(EijR(2,4)+EijR
     -   (20,4))))
       FI(183) = 4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+2
     -   *(Dij1235I(2,3)-Dij1235I(8,3))+s34*(EijI(2,2)+3*EijI(2,3)-Ei
     -   jI(8,2)-4*(EijI(9,3)+EijI(9,4))+EijI(12,3)+2*(EijI(2,4)+EijI
     -   (20,4))))
       F(183)=DCMPLX(FR(183),FI(183))
       FR(184) = 4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+2
     -   *(Dij1235R(2,3)-Dij1235R(8,3))+s34*(EijR(2,1)+4*EijR(2,2)+5*
     -   EijR(2,3)-EijR(4,1)-EijR(8,2)-3*(EijR(9,2)+EijR(9,3)+EijR(10
     -   ,3))+EijR(20,3)+2*(EijR(2,4)-EijR(9,4)-EijR(10,4)+EijR(28,4)
     -   )))
       FI(184) = 4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+2
     -   *(Dij1235I(2,3)-Dij1235I(8,3))+s34*(EijI(2,1)+4*EijI(2,2)+5*
     -   EijI(2,3)-EijI(4,1)-EijI(8,2)-3*(EijI(9,2)+EijI(9,3)+EijI(10
     -   ,3))+EijI(20,3)+2*(EijI(2,4)-EijI(9,4)-EijI(10,4)+EijI(28,4)
     -   )))
       F(184)=DCMPLX(FR(184),FI(184))
       P(523) = p2sq-3*s12+2*P(3)
       P(524) = p2sq-s23-2*P(14)-3*P(16)
       FR(185) = -2*(3*D02345R+Dij2345R(1,2)+8*EijR(21,3)-10*EijR(22,3
     -   )+4*(Dij2345R(1,1)-EijR(37,4)+EijR(40,4))-2*(EijR(11,2)+(Eij
     -   R(1,1)+EijR(1,2)-EijR(2,1))*P(3))+(EijR(9,3)-EijR(17,3))*P(1
     -   0)+(EijR(10,3)-EijR(18,3))*P(21)+(-EijR(2,3)+EijR(8,3))*P(33
     -   )-EijR(2,2)*P(60)+(-EijR(7,2)+EijR(9,2))*P(339)+EijR(5,2)*P(
     -   523)+(-EijR(6,2)+EijR(8,2))*P(524))
       FI(185) = -2*(3*D02345I+Dij2345I(1,2)+8*EijI(21,3)-10*EijI(22,3
     -   )+4*(Dij2345I(1,1)-EijI(37,4)+EijI(40,4))-2*(EijI(11,2)+(Eij
     -   I(1,1)+EijI(1,2)-EijI(2,1))*P(3))+(EijI(9,3)-EijI(17,3))*P(1
     -   0)+(EijI(10,3)-EijI(18,3))*P(21)+(-EijI(2,3)+EijI(8,3))*P(33
     -   )-EijI(2,2)*P(60)+(-EijI(7,2)+EijI(9,2))*P(339)+EijI(5,2)*P(
     -   523)+(-EijI(6,2)+EijI(8,2))*P(524))
       F(185)=DCMPLX(FR(185),FI(185))
       FR(186) = -4*Is34*(D01235R+Dij1235R(1,1)-Dij1235R(2,2)+3*Dij123
     -   5R(4,2)+2*(Dij1235R(2,1)-Dij1235R(2,3)+Dij1235R(6,3))+s34*(E
     -   E0R+EijR(1,1)+EijR(2,2)-EijR(4,1)+4*EijR(5,2)-EijR(7,2)+5*Ei
     -   jR(8,3)+EijR(10,3)+3*(EijR(2,1)-EijR(2,3)-EijR(18,3))-2*(Eij
     -   R(2,4)-EijR(8,4)+EijR(9,2)-EijR(10,4)+EijR(27,4))))
       FI(186) = -4*Is34*(D01235I+Dij1235I(1,1)-Dij1235I(2,2)+3*Dij123
     -   5I(4,2)+2*(Dij1235I(2,1)-Dij1235I(2,3)+Dij1235I(6,3))+s34*(E
     -   E0I+EijI(1,1)+EijI(2,2)-EijI(4,1)+4*EijI(5,2)-EijI(7,2)+5*Ei
     -   jI(8,3)+EijI(10,3)+3*(EijI(2,1)-EijI(2,3)-EijI(18,3))-2*(Eij
     -   I(2,4)-EijI(8,4)+EijI(9,2)-EijI(10,4)+EijI(27,4))))
       F(186)=DCMPLX(FR(186),FI(186))
       FR(187) = -4*Is34*(D01235R+Dij1235R(1,1)-Dij1235R(2,2)+3*Dij123
     -   5R(4,2)+2*(Dij1235R(2,1)-Dij1235R(2,3)+Dij1235R(6,3))+s34*(E
     -   E0R+EijR(1,1)-EijR(2,3)-EijR(3,1)-EijR(6,2)+EijR(9,3)+3*(Eij
     -   R(2,1)+EijR(5,2)+EijR(8,3)-EijR(17,3))+2*(EijR(2,2)-EijR(2,4
     -   )-EijR(8,2)+EijR(8,4)+EijR(9,4)-EijR(26,4))))
       FI(187) = -4*Is34*(D01235I+Dij1235I(1,1)-Dij1235I(2,2)+3*Dij123
     -   5I(4,2)+2*(Dij1235I(2,1)-Dij1235I(2,3)+Dij1235I(6,3))+s34*(E
     -   E0I+EijI(1,1)-EijI(2,3)-EijI(3,1)-EijI(6,2)+EijI(9,3)+3*(Eij
     -   I(2,1)+EijI(5,2)+EijI(8,3)-EijI(17,3))+2*(EijI(2,2)-EijI(2,4
     -   )-EijI(8,2)+EijI(8,4)+EijI(9,4)-EijI(26,4))))
       F(187)=DCMPLX(FR(187),FI(187))
       FR(188) = 4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+2
     -   *(Dij1235R(2,3)-Dij1235R(8,3))+s34*(EijR(2,1)+5*EijR(2,3)-6*
     -   EijR(10,3)+4*(EijR(2,2)-EijR(10,4))+EijR(15,3)+2*(EijR(2,4)-
     -   EijR(9,2)+EijR(21,4))))
       FI(188) = 4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+2
     -   *(Dij1235I(2,3)-Dij1235I(8,3))+s34*(EijI(2,1)+5*EijI(2,3)-6*
     -   EijI(10,3)+4*(EijI(2,2)-EijI(10,4))+EijI(15,3)+2*(EijI(2,4)-
     -   EijI(9,2)+EijI(21,4))))
       F(188)=DCMPLX(FR(188),FI(188))
       FR(189) = 4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+2
     -   *(Dij1235R(2,3)-Dij1235R(8,3))+s34*(EijR(2,2)+EijR(4,1)-EijR
     -   (8,2)+3*(EijR(2,3)-EijR(9,3))-EijR(10,3)+EijR(20,3)+2*(EijR(
     -   2,4)+EijR(9,2)-EijR(9,4)-EijR(10,4)+EijR(28,4))))
       FI(189) = 4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+2
     -   *(Dij1235I(2,3)-Dij1235I(8,3))+s34*(EijI(2,2)+EijI(4,1)-EijI
     -   (8,2)+3*(EijI(2,3)-EijI(9,3))-EijI(10,3)+EijI(20,3)+2*(EijI(
     -   2,4)+EijI(9,2)-EijI(9,4)-EijI(10,4)+EijI(28,4))))
       F(189)=DCMPLX(FR(189),FI(189))
       P(525) = s45+2*P(35)
       P(526) = p5sq+3*s12-s34+s45
       P(527) = p3sq-2*P(526)
       P(528) = -3*s12+2*s45
       P(529) = p3sq-s34+4*P(7)
       P(530) = -p3sq+s12+s45
       P(531) = p3sq-s12+s45
       P(532) = -p2sq+s12+s23+s45
       P(533) = -s15+s23+3*P(7)
       FR(190) = 2*(3*D02345R-Dij2345R(1,1)-12*EijR(11,2)+8*(EijR(22,3
     -   )-EijR(23,3))+4*(EijR(37,4)-EijR(43,4))-EijR(9,2)*P(23)+EijR
     -   (2,1)*P(111)+2*(Dij1235R(3,2)-Dij1235R(6,2)+Dij2345R(2,1)+(E
     -   E0R+EijR(1,1))*P(3)+EijR(10,2)*P(7)-EijR(3,2)*P(198)+(EijR(5
     -   ,2)-EijR(6,2))*P(245))-EijR(2,2)*P(508)+EijR(8,2)*P(511)+Is1
     -   2*(-C0135R+3*C0345R-Cij135R(1,1)+2*(Cij135R(2,1)-Cij135R(2,2
     -   )+Cij135R(3,2))-Cij345R(1,1)-10*Dij1345R(7,2)-4*(Dij1345R(11
     -   ,3)-Dij1345R(12,3))-Dij1345R(6,2)*P(101)+Is45*(2*B012R+4*Cij
     -   145R(4,2)+Cij145R(1,1)*P(95)+C0145R*P(185))+Dij1345R(5,2)*P(
     -   292)+D01345R*P(525)-Dij1345R(1,1)*P(527)+Dij1345R(2,1)*P(528
     -   )+Dij1345R(3,1)*P(529)-(-Dij1345R(1,2)+Dij1345R(4,2))*P(530)
     -   -Dij1345R(2,2)*P(531))+EijR(3,1)*P(532)+EijR(4,1)*P(533))
       FI(190) = 2*(3*D02345I-Dij2345I(1,1)-12*EijI(11,2)+8*(EijI(22,3
     -   )-EijI(23,3))+4*(EijI(37,4)-EijI(43,4))-EijI(9,2)*P(23)+EijI
     -   (2,1)*P(111)+2*(Dij1235I(3,2)-Dij1235I(6,2)+Dij2345I(2,1)+(E
     -   E0I+EijI(1,1))*P(3)+EijI(10,2)*P(7)-EijI(3,2)*P(198)+(EijI(5
     -   ,2)-EijI(6,2))*P(245))-EijI(2,2)*P(508)+EijI(8,2)*P(511)+Is1
     -   2*(-C0135I+3*C0345I-Cij135I(1,1)+2*(Cij135I(2,1)-Cij135I(2,2
     -   )+Cij135I(3,2))-Cij345I(1,1)-10*Dij1345I(7,2)-4*(Dij1345I(11
     -   ,3)-Dij1345I(12,3))-Dij1345I(6,2)*P(101)+Is45*(2*B012I+4*Cij
     -   145I(4,2)+Cij145I(1,1)*P(95)+C0145I*P(185))+Dij1345I(5,2)*P(
     -   292)+D01345I*P(525)-Dij1345I(1,1)*P(527)+Dij1345I(2,1)*P(528
     -   )+Dij1345I(3,1)*P(529)-(-Dij1345I(1,2)+Dij1345I(4,2))*P(530)
     -   -Dij1345I(2,2)*P(531))+EijI(3,1)*P(532)+EijI(4,1)*P(533))
       F(190)=DCMPLX(FR(190),FI(190))
       P(534) = 5*s12+4*P(39)
       P(535) = p5sq+s12-s34-2*s45
       P(536) = p3sq+2*P(535)
       P(537) = p5sq+s45
       P(538) = -3*s12-4*s34+2*P(537)
       P(539) = 2*p2sq-s12
       P(540) = p2sq-s12-s23+3*s45-2*P(3)
       P(541) = p2sq-s12-s23+2*s34+s45
       P(542) = -s23+s34+s45
       P(543) = 5*p2sq-3*s12+2*P(542)
       FR(191) = -2*(3*D02345R+Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(4,
     -   2)-12*EijR(11,2)-8*(EijR(22,3)-EijR(23,3))-4*(EijR(37,4)-Eij
     -   R(43,4))-EijR(12,3)*P(10)+(2*EijR(9,2)-EijR(10,2)+EijR(10,3)
     -   -EijR(20,3))*P(21)-EijR(2,3)*P(33)+EijR(9,3)*P(245)-EijR(2,1
     -   )*P(375)+Is34*(3*C0235R+Cij235R(1,1)-8*Dij1235R(7,2)-4*Dij12
     -   35R(12,3)+2*(D01235R+Dij1235R(1,1))*P(3)-Dij1235R(2,2)*P(33)
     -   -Dij1235R(2,1)*P(60)+Dij1235R(6,2)*P(330)+Is12*(2*B035R+4*Ci
     -   j135R(4,2)-C0135R*P(513)+Cij135R(1,1)*P(514))+Dij1235R(3,1)*
     -   P(516))+EijR(4,1)*P(533)+Is12*(Cij135R(2,1)-Cij345R(1,1)-8*D
     -   ij1345R(7,2)+4*(C0345R+Dij1345R(11,3)-Dij1345R(12,3))+Dij134
     -   5R(6,2)*P(7)+2*D01345R*P(35)-Dij1345R(2,2)*P(199)+Dij1345R(3
     -   ,1)*P(328)+Is45*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)*P(95)-C
     -   0145R*P(387))+Dij1345R(5,2)*P(491)-Dij1345R(1,2)*P(513)+Dij1
     -   345R(1,1)*P(534)-Dij1345R(2,1)*P(536)-Dij1345R(4,2)*P(538))+
     -   2*(Dij2345R(1,1)+(EE0R+EijR(1,1))*P(3)+(EijR(5,2)-EijR(6,2))
     -   *P(330)-EijR(2,2)*P(539))+EijR(3,1)*P(540)-EijR(3,2)*P(541)+
     -   EijR(8,2)*P(543))
       FI(191) = -2*(3*D02345I+Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(4,
     -   2)-12*EijI(11,2)-8*(EijI(22,3)-EijI(23,3))-4*(EijI(37,4)-Eij
     -   I(43,4))-EijI(12,3)*P(10)+(2*EijI(9,2)-EijI(10,2)+EijI(10,3)
     -   -EijI(20,3))*P(21)-EijI(2,3)*P(33)+EijI(9,3)*P(245)-EijI(2,1
     -   )*P(375)+Is34*(3*C0235I+Cij235I(1,1)-8*Dij1235I(7,2)-4*Dij12
     -   35I(12,3)+2*(D01235I+Dij1235I(1,1))*P(3)-Dij1235I(2,2)*P(33)
     -   -Dij1235I(2,1)*P(60)+Dij1235I(6,2)*P(330)+Is12*(2*B035I+4*Ci
     -   j135I(4,2)-C0135I*P(513)+Cij135I(1,1)*P(514))+Dij1235I(3,1)*
     -   P(516))+EijI(4,1)*P(533)+Is12*(Cij135I(2,1)-Cij345I(1,1)-8*D
     -   ij1345I(7,2)+4*(C0345I+Dij1345I(11,3)-Dij1345I(12,3))+Dij134
     -   5I(6,2)*P(7)+2*D01345I*P(35)-Dij1345I(2,2)*P(199)+Dij1345I(3
     -   ,1)*P(328)+Is45*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)*P(95)-C
     -   0145I*P(387))+Dij1345I(5,2)*P(491)-Dij1345I(1,2)*P(513)+Dij1
     -   345I(1,1)*P(534)-Dij1345I(2,1)*P(536)-Dij1345I(4,2)*P(538))+
     -   2*(Dij2345I(1,1)+(EE0I+EijI(1,1))*P(3)+(EijI(5,2)-EijI(6,2))
     -   *P(330)-EijI(2,2)*P(539))+EijI(3,1)*P(540)-EijI(3,2)*P(541)+
     -   EijI(8,2)*P(543))
       F(191)=DCMPLX(FR(191),FI(191))
       P(544) = p3sq-2*s45
       P(545) = 2*p3sq-3*s45
       P(546) = s12+s23
       FR(192) = 2*(-Dij1235R(3,1)+Dij1345R(1,2)-Dij1345R(4,2)+p2sq*(E
     -   ijR(1,1)+2*(EijR(2,1)+EijR(5,2)))-2*(Dij1235R(6,2)-Dij1345R(
     -   1,1)+Dij1345R(2,1)-EijR(11,2)-3*EijR(22,3)-2*EijR(37,4))+EE0
     -   R*P(10)-Is12*(C0135R+C0145R-C0345R+Cij135R(1,1)-3*(Cij135R(2
     -   ,1)-Cij145R(1,1))-p5sq*(2*Dij1345R(3,1)+Dij1345R(5,2))+s34*(
     -   Dij1345R(3,1)+2*Dij1345R(5,2))+4*(Dij1345R(7,2)+Dij1345R(11,
     -   3))+p3sq*(Dij1345R(1,2)+Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R
     -   (4,2)-2*Dij1345R(5,2)-s12*EijR(4,1))+s12*(s34*EijR(4,1)+2*(s
     -   12+s23)*(EijR(2,1)+EijR(5,2)))-s45*(Dij1345R(1,2)+2*(Dij1345
     -   R(2,1)-Dij1345R(3,1))+Dij1345R(4,2)-Dij1345R(5,2)+s12*(EijR(
     -   1,1)+2*(EijR(2,1)+EijR(5,2))))-2*(Cij135R(3,2)-Cij145R(1,2)+
     -   s12*EijR(9,2)*P(14))+D01345R*P(544)+Dij1345R(1,1)*P(545)+s12
     -   *EijR(1,1)*P(546)))
       FI(192) = 2*(-Dij1235I(3,1)+Dij1345I(1,2)-Dij1345I(4,2)+p2sq*(E
     -   ijI(1,1)+2*(EijI(2,1)+EijI(5,2)))-2*(Dij1235I(6,2)-Dij1345I(
     -   1,1)+Dij1345I(2,1)-EijI(11,2)-3*EijI(22,3)-2*EijI(37,4))+EE0
     -   I*P(10)-Is12*(C0135I+C0145I-C0345I+Cij135I(1,1)-3*(Cij135I(2
     -   ,1)-Cij145I(1,1))-p5sq*(2*Dij1345I(3,1)+Dij1345I(5,2))+s34*(
     -   Dij1345I(3,1)+2*Dij1345I(5,2))+4*(Dij1345I(7,2)+Dij1345I(11,
     -   3))+p3sq*(Dij1345I(1,2)+Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I
     -   (4,2)-2*Dij1345I(5,2)-s12*EijI(4,1))+s12*(s34*EijI(4,1)+2*(s
     -   12+s23)*(EijI(2,1)+EijI(5,2)))-s45*(Dij1345I(1,2)+2*(Dij1345
     -   I(2,1)-Dij1345I(3,1))+Dij1345I(4,2)-Dij1345I(5,2)+s12*(EijI(
     -   1,1)+2*(EijI(2,1)+EijI(5,2))))-2*(Cij135I(3,2)-Cij145I(1,2)+
     -   s12*EijI(9,2)*P(14))+D01345I*P(544)+Dij1345I(1,1)*P(545)+s12
     -   *EijI(1,1)*P(546)))
       F(192)=DCMPLX(FR(192),FI(192))
       P(547) = p3sq+s12+s23-2*s45
       P(548) = p3sq+s12+s23+2*P(414)
       P(549) = 3*p3sq+s12-s45
       P(550) = p5sq-s45+4*P(14)
       FR(193) = 2*(D01235R-EE0R*s45+Dij1235R(1,1)-EijR(4,1)*P(7)+EijR
     -   (1,2)*P(50)+EijR(1,1)*P(77)+(EijR(3,1)-2*EijR(8,2))*P(198)+2
     -   *(-Dij1235R(3,1)+Dij1235R(5,2)-Dij1235R(6,2)+Dij2345R(1,1)+3
     -   *EijR(11,2)+2*(EijR(37,4)-EijR(40,4))+EijR(2,2)*P(5)+EijR(9,
     -   2)*P(292))-EijR(2,1)*P(547)-EijR(5,2)*P(548)+EijR(6,2)*P(549
     -   )-EijR(7,2)*P(550))
       FI(193) = 2*(D01235I-EE0I*s45+Dij1235I(1,1)-EijI(4,1)*P(7)+EijI
     -   (1,2)*P(50)+EijI(1,1)*P(77)+(EijI(3,1)-2*EijI(8,2))*P(198)+2
     -   *(-Dij1235I(3,1)+Dij1235I(5,2)-Dij1235I(6,2)+Dij2345I(1,1)+3
     -   *EijI(11,2)+2*(EijI(37,4)-EijI(40,4))+EijI(2,2)*P(5)+EijI(9,
     -   2)*P(292))-EijI(2,1)*P(547)-EijI(5,2)*P(548)+EijI(6,2)*P(549
     -   )-EijI(7,2)*P(550))
       F(193)=DCMPLX(FR(193),FI(193))
       P(551) = p3sq+p5sq+s12-s34-2*s45
       P(552) = p3sq-s12-2*s34+s45
       P(553) = -3*p5sq+s34+2*P(104)
       P(554) = p2sq-s23+2*s34-3*P(16)
       P(555) = p5sq+s15+s23-3*s45+2*P(6)
       P(556) = p2sq-s12-s23+s34+s45
       P(557) = 3*p2sq-s12+2*P(107)
       FR(194) = -2*(3*D02345R+Dij2345R(1,2)-Dij2345R(3,1)-Dij2345R(5,
     -   2)-12*(EijR(11,2)+EijR(22,3))+8*EijR(24,3)-4*(EijR(37,4)-Eij
     -   R(44,4))+(EijR(9,3)-EijR(20,3))*P(10)-(EijR(4,2)+EijR(15,3))
     -   *P(21)-EijR(2,3)*P(33)-EijR(2,1)*P(375)+EijR(10,3)*P(488)+Is
     -   34*(3*C0235R+Cij235R(1,1)-8*Dij1235R(7,2)-4*Dij1235R(12,3)+2
     -   *(D01235R+Dij1235R(1,1))*P(3)-Dij1235R(2,2)*P(33)-Dij1235R(2
     -   ,1)*P(60)+Dij1235R(6,2)*P(330)+Is12*(2*B035R+4*Cij135R(4,2)-
     -   C0135R*P(513)+Cij135R(1,1)*P(514))+Dij1235R(3,1)*P(516))-Eij
     -   R(10,2)*P(541)+Is12*(Cij135R(2,1)-Cij345R(2,1)+4*(C0345R+Dij
     -   1345R(11,3)-Dij1345R(13,3))+Dij1345R(3,2)*P(7)-2*(Dij1345R(7
     -   ,2)-D01345R*P(35))-Dij1345R(6,2)*P(199)+Is45*(2*B012R+4*Cij1
     -   45R(4,2)+Cij145R(1,1)*P(95)-C0145R*P(387))-Dij1345R(2,1)*P(4
     -   21)-Dij1345R(1,2)*P(513)+Dij1345R(1,1)*P(534)+Dij1345R(3,1)*
     -   P(551)-Dij1345R(4,2)*P(552)+Dij1345R(5,2)*P(553))+EijR(3,1)*
     -   P(554)+EijR(4,1)*P(555)+2*(Dij2345R(1,1)+(EE0R+EijR(1,1))*P(
     -   3)+(EijR(5,2)-EijR(7,2))*P(330)-EijR(2,2)*P(539)+EijR(8,2)*P
     -   (556))+EijR(9,2)*P(557))
       FI(194) = -2*(3*D02345I+Dij2345I(1,2)-Dij2345I(3,1)-Dij2345I(5,
     -   2)-12*(EijI(11,2)+EijI(22,3))+8*EijI(24,3)-4*(EijI(37,4)-Eij
     -   I(44,4))+(EijI(9,3)-EijI(20,3))*P(10)-(EijI(4,2)+EijI(15,3))
     -   *P(21)-EijI(2,3)*P(33)-EijI(2,1)*P(375)+EijI(10,3)*P(488)+Is
     -   34*(3*C0235I+Cij235I(1,1)-8*Dij1235I(7,2)-4*Dij1235I(12,3)+2
     -   *(D01235I+Dij1235I(1,1))*P(3)-Dij1235I(2,2)*P(33)-Dij1235I(2
     -   ,1)*P(60)+Dij1235I(6,2)*P(330)+Is12*(2*B035I+4*Cij135I(4,2)-
     -   C0135I*P(513)+Cij135I(1,1)*P(514))+Dij1235I(3,1)*P(516))-Eij
     -   I(10,2)*P(541)+Is12*(Cij135I(2,1)-Cij345I(2,1)+4*(C0345I+Dij
     -   1345I(11,3)-Dij1345I(13,3))+Dij1345I(3,2)*P(7)-2*(Dij1345I(7
     -   ,2)-D01345I*P(35))-Dij1345I(6,2)*P(199)+Is45*(2*B012I+4*Cij1
     -   45I(4,2)+Cij145I(1,1)*P(95)-C0145I*P(387))-Dij1345I(2,1)*P(4
     -   21)-Dij1345I(1,2)*P(513)+Dij1345I(1,1)*P(534)+Dij1345I(3,1)*
     -   P(551)-Dij1345I(4,2)*P(552)+Dij1345I(5,2)*P(553))+EijI(3,1)*
     -   P(554)+EijI(4,1)*P(555)+2*(Dij2345I(1,1)+(EE0I+EijI(1,1))*P(
     -   3)+(EijI(5,2)-EijI(7,2))*P(330)-EijI(2,2)*P(539)+EijI(8,2)*P
     -   (556))+EijI(9,2)*P(557))
       F(194)=DCMPLX(FR(194),FI(194))
       P(558) = p2sq-p5sq+s12-s34+s45
       FR(195) = -4*(Dij2345R(1,1)-Dij2345R(3,1)+2*(EijR(11,2)-2*EijR(
     -   22,3)+EijR(24,3)-EijR(37,4)+EijR(44,4))+(EijR(2,1)-EijR(4,1)
     -   )*P(6)-EijR(4,2)*P(7)+(-EijR(8,2)+EijR(10,2))*P(16)+(-EijR(5
     -   ,2)+EijR(7,2))*P(33)+EijR(2,2)*P(521)-EijR(9,2)*P(558))
       FI(195) = -4*(Dij2345I(1,1)-Dij2345I(3,1)+2*(EijI(11,2)-2*EijI(
     -   22,3)+EijI(24,3)-EijI(37,4)+EijI(44,4))+(EijI(2,1)-EijI(4,1)
     -   )*P(6)-EijI(4,2)*P(7)+(-EijI(8,2)+EijI(10,2))*P(16)+(-EijI(5
     -   ,2)+EijI(7,2))*P(33)+EijI(2,2)*P(521)-EijI(9,2)*P(558))
       F(195)=DCMPLX(FR(195),FI(195))
       P(559) = 2+Is12*s45
       P(560) = -p3sq+5*s12+2*s45
       P(561) = p3sq+3*P(16)
       P(562) = 3*p5sq+s12-5*s45+2*P(14)
       P(563) = 4*p3sq+p5sq-2*s34-3*s45
       P(564) = p2sq-2*s12-s23+s45
       P(565) = p2sq+p3sq-2*s12-s23+s45
       P(566) = p2sq-p5sq+s12+s45
       FR(196) = -2*(-(D01345R*P(559))+Is12*(C0135R+Cij135R(1,1)-3*Cij
     -   145R(2,1)+2*(C0145R-C0345R-Cij135R(2,1)+Cij135R(2,2)-Cij135R
     -   (3,2)+Cij145R(1,2)-Cij145R(3,2))+Cij345R(2,1)+6*Dij1345R(7,2
     -   )+4*(Cij145R(1,1)+Dij1345R(11,3)-Dij1345R(13,3))+Dij1345R(3,
     -   2)*P(101)+Dij1345R(4,2)*P(198)-Dij1345R(1,2)*P(530)+Dij1345R
     -   (6,2)*P(531)-Dij1345R(1,1)*P(560)+Dij1345R(2,1)*P(561)-Dij13
     -   45R(3,1)*P(562)-Dij1345R(5,2)*P(563))+2*(-Dij1235R(3,2)+Dij1
     -   235R(6,2)+Dij2345R(1,1)-Dij2345R(3,1)-p3sq*EijR(3,1)+2*(EijR
     -   (11,2)-2*(EijR(22,3)-EijR(24,3))-EijR(37,4)+EijR(44,4))+EijR
     -   (2,2)*P(5)-EijR(4,2)*P(7)-(EijR(8,2)-EijR(10,2))*P(198)-(Eij
     -   R(5,2)-EijR(7,2))*P(245)-EijR(2,1)*P(564)+EijR(4,1)*P(565)-E
     -   ijR(9,2)*P(566)))
       FI(196) = -2*(-(D01345I*P(559))+Is12*(C0135I+Cij135I(1,1)-3*Cij
     -   145I(2,1)+2*(C0145I-C0345I-Cij135I(2,1)+Cij135I(2,2)-Cij135I
     -   (3,2)+Cij145I(1,2)-Cij145I(3,2))+Cij345I(2,1)+6*Dij1345I(7,2
     -   )+4*(Cij145I(1,1)+Dij1345I(11,3)-Dij1345I(13,3))+Dij1345I(3,
     -   2)*P(101)+Dij1345I(4,2)*P(198)-Dij1345I(1,2)*P(530)+Dij1345I
     -   (6,2)*P(531)-Dij1345I(1,1)*P(560)+Dij1345I(2,1)*P(561)-Dij13
     -   45I(3,1)*P(562)-Dij1345I(5,2)*P(563))+2*(-Dij1235I(3,2)+Dij1
     -   235I(6,2)+Dij2345I(1,1)-Dij2345I(3,1)-p3sq*EijI(3,1)+2*(EijI
     -   (11,2)-2*(EijI(22,3)-EijI(24,3))-EijI(37,4)+EijI(44,4))+EijI
     -   (2,2)*P(5)-EijI(4,2)*P(7)-(EijI(8,2)-EijI(10,2))*P(198)-(Eij
     -   I(5,2)-EijI(7,2))*P(245)-EijI(2,1)*P(564)+EijI(4,1)*P(565)-E
     -   ijI(9,2)*P(566)))
       F(196)=DCMPLX(FR(196),FI(196))
       FR(197) = 4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)+2*Dij1235R(2,3)
     -   +s34*(EijR(2,2)-EijR(8,2)+3*(EijR(2,3)-EijR(9,3))+2*(EijR(2,
     -   4)-EijR(9,4))))
       FI(197) = 4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)+2*Dij1235I(2,3)
     -   +s34*(EijI(2,2)-EijI(8,2)+3*(EijI(2,3)-EijI(9,3))+2*(EijI(2,
     -   4)-EijI(9,4))))
       F(197)=DCMPLX(FR(197),FI(197))
       FR(198) = -4*(EE0R+EijR(1,1)+EijR(2,2)+3*(EijR(2,1)-EijR(2,3))+
     -   4*EijR(5,2)+5*EijR(8,3)-2*(EijR(2,4)-EijR(8,4)))
       FI(198) = -4*(EE0I+EijI(1,1)+EijI(2,2)+3*(EijI(2,1)-EijI(2,3))+
     -   4*EijI(5,2)+5*EijI(8,3)-2*(EijI(2,4)-EijI(8,4)))
       F(198)=DCMPLX(FR(198),FI(198))
       FR(199) = 4*(EijR(2,1)+4*EijR(2,2)+5*EijR(2,3)-EijR(9,2)-3*EijR
     -   (10,3)+2*(EijR(2,4)-EijR(10,4)))
       FI(199) = 4*(EijI(2,1)+4*EijI(2,2)+5*EijI(2,3)-EijI(9,2)-3*EijI
     -   (10,3)+2*(EijI(2,4)-EijI(10,4)))
       F(199)=DCMPLX(FR(199),FI(199))
       FR(200) = 4*(EijR(2,1)+4*EijR(2,2)+5*EijR(2,3)-EijR(8,2)-3*EijR
     -   (9,3)+2*(EijR(2,4)-EijR(9,4)))
       FI(200) = 4*(EijI(2,1)+4*EijI(2,2)+5*EijI(2,3)-EijI(8,2)-3*EijI
     -   (9,3)+2*(EijI(2,4)-EijI(9,4)))
       F(200)=DCMPLX(FR(200),FI(200))
       FR(201) = 4*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)+2*Dij1235R(2,3)
     -   +s34*(EijR(2,1)+4*EijR(2,2)+5*EijR(2,3)-EijR(9,2)-3*EijR(10,
     -   3)+2*(EijR(2,4)-EijR(10,4))))
       FI(201) = 4*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)+2*Dij1235I(2,3)
     -   +s34*(EijI(2,1)+4*EijI(2,2)+5*EijI(2,3)-EijI(9,2)-3*EijI(10,
     -   3)+2*(EijI(2,4)-EijI(10,4))))
       F(201)=DCMPLX(FR(201),FI(201))
       P(567) = s12+4*P(39)
       P(568) = p2sq-3*s12-2*P(3)
       P(569) = p2sq-p3sq-s23+s34-2*P(16)
       P(570) = p3sq-s15+s23-s34+2*P(7)
       FR(202) = -2*(3*D02345R+Dij2345R(1,2)+4*(Dij2345R(1,1)-EijR(11,
     -   2)-3*EijR(22,3)-EijR(37,4))+EijR(9,3)*P(10)-EijR(2,3)*P(33)-
     -   EijR(10,3)*P(115)+EijR(4,1)*P(339)+EijR(3,1)*P(524)+Is12*(2*
     -   (C0345R+Dij1345R(7,2)+2*Dij1345R(11,3)+D01345R*P(39))-Dij134
     -   5R(2,1)*P(79)+Dij1345R(3,1)*P(101)+Is45*(2*B012R+4*Cij145R(4
     -   ,2)+Cij145R(1,1)*P(95)-C0145R*P(387))+Dij1345R(5,2)*P(491)-D
     -   ij1345R(1,2)*P(513)-Dij1345R(4,2)*P(552)+Dij1345R(1,1)*P(567
     -   ))-EijR(2,1)*P(568)+2*((EE0R+EijR(1,1)+EijR(5,2))*P(3)-EijR(
     -   2,2)*P(72)+EijR(8,2)*P(569)+EijR(9,2)*P(570)))
       FI(202) = -2*(3*D02345I+Dij2345I(1,2)+4*(Dij2345I(1,1)-EijI(11,
     -   2)-3*EijI(22,3)-EijI(37,4))+EijI(9,3)*P(10)-EijI(2,3)*P(33)-
     -   EijI(10,3)*P(115)+EijI(4,1)*P(339)+EijI(3,1)*P(524)+Is12*(2*
     -   (C0345I+Dij1345I(7,2)+2*Dij1345I(11,3)+D01345I*P(39))-Dij134
     -   5I(2,1)*P(79)+Dij1345I(3,1)*P(101)+Is45*(2*B012I+4*Cij145I(4
     -   ,2)+Cij145I(1,1)*P(95)-C0145I*P(387))+Dij1345I(5,2)*P(491)-D
     -   ij1345I(1,2)*P(513)-Dij1345I(4,2)*P(552)+Dij1345I(1,1)*P(567
     -   ))-EijI(2,1)*P(568)+2*((EE0I+EijI(1,1)+EijI(5,2))*P(3)-EijI(
     -   2,2)*P(72)+EijI(8,2)*P(569)+EijI(9,2)*P(570)))
       F(202)=DCMPLX(FR(202),FI(202))
       FR(203) = 4*(EijR(2,1)+4*EijR(2,2)+5*EijR(2,3)+2*EijR(2,4))
       FI(203) = 4*(EijI(2,1)+4*EijI(2,2)+5*EijI(2,3)+2*EijI(2,4))
       F(203)=DCMPLX(FR(203),FI(203))
       P(571) = p5sq+s15+s23+s45
       P(572) = 7*s12-13*s34+3*P(299)-2*P(571)
       P(573) = s23-s45+2*P(33)
       P(574) = p2sq-p3sq+4*s12+s23
       P(575) = p2sq-p3sq-s23-4*P(16)
       P(576) = -s15+s23+4*P(7)+2*P(14)
       FR(204) = Is12*((Is45*P(6)*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)
     -   *P(95)+C0145R*P(185)))/2.d0+(-4*Cij135R(4,2)+s34*Dij1345R(6,2)
     -   *P(275)+Cij135R(2,1)*P(320)+Dij1345R(1,1)*P(323)-Dij1345R(2,
     -   1)*P(326)+P(6)*(5*C0345R+Dij1345R(3,1)*P(328))-Dij1345R(4,2)
     -   *P(331)-2*(B035R-p3sq*s34*Dij1345R(2,2)+C0135R*P(35)+Cij135R
     -   (1,1)*P(59)+Dij1345R(7,2)*P(333))+D01345R*P(428))/2.d0)+(P(2)*
     -   (Dij2345R(3,1)-EijR(4,2)*P(7))+Dij1345R(5,2)*P(14)+(Dij1235R
     -   (2,1)-Dij1235R(4,2)+Dij1235R(5,2))*P(33)+Dij1345R(1,2)*P(198
     -   )+D02345R*P(429)+Dij2345R(1,1)*P(431)-Dij2345R(2,1)*P(432)-D
     -   ij1235R(3,1)*P(434)+Dij1235R(3,2)*P(436)-Dij1235R(6,2)*P(437
     -   )+EijR(2,1)*P(439)-EijR(2,2)*P(442)-EijR(3,1)*P(446)-EijR(3,
     -   2)*P(448)+EijR(4,1)*P(451)+EijR(5,2)*P(454)-EijR(6,2)*P(459)
     -   +EijR(7,2)*P(461)+EijR(8,2)*P(466)-EijR(9,2)*P(470)+EijR(10,
     -   2)*P(474)-2*(C0235R-p2sq*Dij1235R(2,2)-4*(Dij1235R(7,2)-Dij2
     -   345R(7,2))+P(3)*(D01235R+Dij1235R(1,1)-EE0R*P(6)+EijR(1,2)*P
     -   (33)+EijR(1,1)*P(71))+EijR(11,2)*P(572)-EijR(21,3)*P(573)+Ei
     -   jR(22,3)*P(574)+EijR(23,3)*P(575)+EijR(24,3)*P(576)))/2.d0
       FI(204) = Is12*((Is45*P(6)*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)
     -   *P(95)+C0145I*P(185)))/2.d0+(-4*Cij135I(4,2)+s34*Dij1345I(6,2)
     -   *P(275)+Cij135I(2,1)*P(320)+Dij1345I(1,1)*P(323)-Dij1345I(2,
     -   1)*P(326)+P(6)*(5*C0345I+Dij1345I(3,1)*P(328))-Dij1345I(4,2)
     -   *P(331)-2*(B035I-p3sq*s34*Dij1345I(2,2)+C0135I*P(35)+Cij135I
     -   (1,1)*P(59)+Dij1345I(7,2)*P(333))+D01345I*P(428))/2.d0)+(P(2)*
     -   (Dij2345I(3,1)-EijI(4,2)*P(7))+Dij1345I(5,2)*P(14)+(Dij1235I
     -   (2,1)-Dij1235I(4,2)+Dij1235I(5,2))*P(33)+Dij1345I(1,2)*P(198
     -   )+D02345I*P(429)+Dij2345I(1,1)*P(431)-Dij2345I(2,1)*P(432)-D
     -   ij1235I(3,1)*P(434)+Dij1235I(3,2)*P(436)-Dij1235I(6,2)*P(437
     -   )+EijI(2,1)*P(439)-EijI(2,2)*P(442)-EijI(3,1)*P(446)-EijI(3,
     -   2)*P(448)+EijI(4,1)*P(451)+EijI(5,2)*P(454)-EijI(6,2)*P(459)
     -   +EijI(7,2)*P(461)+EijI(8,2)*P(466)-EijI(9,2)*P(470)+EijI(10,
     -   2)*P(474)-2*(C0235I-p2sq*Dij1235I(2,2)-4*(Dij1235I(7,2)-Dij2
     -   345I(7,2))+P(3)*(D01235I+Dij1235I(1,1)-EE0I*P(6)+EijI(1,2)*P
     -   (33)+EijI(1,1)*P(71))+EijI(11,2)*P(572)-EijI(21,3)*P(573)+Ei
     -   jI(22,3)*P(574)+EijI(23,3)*P(575)+EijI(24,3)*P(576)))/2.d0
       F(204)=DCMPLX(FR(204),FI(204))
       P(577) = p3sq-s12-2*s45
       P(578) = -p3sq+s12+2*s45
       P(579) = p3sq-3*s12-s45
       P(580) = 2*p2sq+p3sq+3*s12-s45
       P(581) = p3sq-s45
       P(582) = p2sq+3*s12+2*P(581)
       FR(205) = Is34*(C0235R-Cij235R(1,1)+Cij235R(2,1)-4*(Dij1235R(7,
     -   2)-Dij1235R(12,3)+Dij1235R(13,3))+2*(Cij135R(1,2)+(D01235R+D
     -   ij1235R(1,1))*P(3))-Dij1235R(2,2)*P(5)+Is12*(4*Cij135R(4,2)+
     -   C0135R*P(26)+Cij135R(1,1)*P(27)-Cij135R(2,1)*P(28)+2*(B035R-
     -   Cij135R(3,2)*P(29)))+(Dij1235R(2,1)+2*(Dij1235R(4,2)-Dij1235
     -   R(5,2)))*P(33)-Dij1235R(3,1)*P(34)+Dij1235R(3,2)*P(38)+Dij12
     -   35R(6,2)*P(42))+Is12*(-(D01345R*s45)+2*(C0345R+Cij135R(2,2))
     -   +Cij345R(1,1)-4*(Dij1345R(7,2)-Dij1345R(11,3)+Dij1345R(12,3)
     -   )-Dij1345R(3,1)*P(7)-Dij1345R(2,2)*P(198)+(-Dij1345R(5,2)+Di
     -   j1345R(6,2))*P(292)+Dij1345R(1,1)*P(577)+Dij1345R(2,1)*P(578
     -   )+Dij1345R(1,2)*P(579))-2*(-Dij2345R(1,1)+Dij2345R(1,2)+Dij2
     -   345R(2,1)+Dij2345R(2,2)-2*(Dij1345R(4,2)+Dij2345R(4,2))+s45*
     -   EijR(3,2)-s45*EijR(8,2)+8*(EijR(22,3)-EijR(23,3))+EijR(2,3)*
     -   P(5)+(-EijR(9,2)+EijR(10,2))*P(23)-(EijR(8,3)+EijR(11,3)-2*E
     -   ijR(17,3))*P(33)+(-EijR(2,1)+EijR(3,1)-EijR(5,2)+EijR(6,2))*
     -   P(50)+(EijR(10,3)-2*EijR(20,3))*P(101)-EijR(3,3)*P(198)-EijR
     -   (13,3)*P(312)-EijR(9,3)*P(580)+EijR(12,3)*P(582))
       FI(205) = Is34*(C0235I-Cij235I(1,1)+Cij235I(2,1)-4*(Dij1235I(7,
     -   2)-Dij1235I(12,3)+Dij1235I(13,3))+2*(Cij135I(1,2)+(D01235I+D
     -   ij1235I(1,1))*P(3))-Dij1235I(2,2)*P(5)+Is12*(4*Cij135I(4,2)+
     -   C0135I*P(26)+Cij135I(1,1)*P(27)-Cij135I(2,1)*P(28)+2*(B035I-
     -   Cij135I(3,2)*P(29)))+(Dij1235I(2,1)+2*(Dij1235I(4,2)-Dij1235
     -   I(5,2)))*P(33)-Dij1235I(3,1)*P(34)+Dij1235I(3,2)*P(38)+Dij12
     -   35I(6,2)*P(42))+Is12*(-(D01345I*s45)+2*(C0345I+Cij135I(2,2))
     -   +Cij345I(1,1)-4*(Dij1345I(7,2)-Dij1345I(11,3)+Dij1345I(12,3)
     -   )-Dij1345I(3,1)*P(7)-Dij1345I(2,2)*P(198)+(-Dij1345I(5,2)+Di
     -   j1345I(6,2))*P(292)+Dij1345I(1,1)*P(577)+Dij1345I(2,1)*P(578
     -   )+Dij1345I(1,2)*P(579))-2*(-Dij2345I(1,1)+Dij2345I(1,2)+Dij2
     -   345I(2,1)+Dij2345I(2,2)-2*(Dij1345I(4,2)+Dij2345I(4,2))+s45*
     -   EijI(3,2)-s45*EijI(8,2)+8*(EijI(22,3)-EijI(23,3))+EijI(2,3)*
     -   P(5)+(-EijI(9,2)+EijI(10,2))*P(23)-(EijI(8,3)+EijI(11,3)-2*E
     -   ijI(17,3))*P(33)+(-EijI(2,1)+EijI(3,1)-EijI(5,2)+EijI(6,2))*
     -   P(50)+(EijI(10,3)-2*EijI(20,3))*P(101)-EijI(3,3)*P(198)-EijI
     -   (13,3)*P(312)-EijI(9,3)*P(580)+EijI(12,3)*P(582))
       F(205)=DCMPLX(FR(205),FI(205))
       P(583) = -1+Is34*s12
       P(584) = p3sq-s12+s23-2*s45
       P(585) = -p3sq+s12
       P(586) = p3sq-s12+s23+2*s45
       P(587) = 3*P(7)+2*P(66)
       P(588) = -p3sq+s12+s23+2*s45
       P(589) = 2*p2sq+p3sq+s12-s45
       FR(206) = EE0R*s45+3*Dij2345R(2,1)+(EijR(4,1)+EijR(10,2))*P(7)-
     -   EijR(1,2)*P(50)+Is15s34*(-2*(B025R+2*Cij125R(4,2))+(C0125R+C
     -   ij125R(1,1))*P(70))+Is34*(Cij125R(2,1)-2*(C0235R+Cij125R(2,2
     -   )-Cij125R(3,2))-Cij235R(1,1)+4*(Dij1235R(7,2)-Dij1235R(11,3)
     -   +Dij1235R(12,3))+(-Dij1235R(1,2)+Dij1235R(2,2))*P(33)-Dij123
     -   5R(1,1)*P(71)+Dij1235R(2,1)*P(72)+Dij1235R(3,1)*P(73)+Dij123
     -   5R(5,2)*P(75)+Dij1235R(6,2)*P(76))-EijR(1,1)*P(77)+EijR(3,2)
     -   *P(531)+D01235R*P(583)+(EijR(2,1)+EijR(5,2))*P(584)+EijR(3,1
     -   )*P(585)-EijR(6,2)*P(586)-EijR(7,2)*P(587)+EijR(8,2)*P(588)+
     -   2*(D02345R-Dij2345R(1,2)+Dij2345R(4,2)+2*(p2sq*EijR(8,3)+Eij
     -   R(11,2))+6*(EijR(21,3)-EijR(22,3))-EijR(2,3)*P(5)+EijR(9,2)*
     -   P(23)+(-EijR(5,3)+EijR(6,3))*P(33)+(-EijR(10,3)+EijR(18,3)-E
     -   ijR(19,3)+EijR(20,3))*P(101)+(EijR(11,3)-EijR(12,3))*P(198)+
     -   EijR(9,3)*P(505)-EijR(17,3)*P(589))
       FI(206) = EE0I*s45+3*Dij2345I(2,1)+(EijI(4,1)+EijI(10,2))*P(7)-
     -   EijI(1,2)*P(50)+Is15s34*(-2*(B025I+2*Cij125I(4,2))+(C0125I+C
     -   ij125I(1,1))*P(70))+Is34*(Cij125I(2,1)-2*(C0235I+Cij125I(2,2
     -   )-Cij125I(3,2))-Cij235I(1,1)+4*(Dij1235I(7,2)-Dij1235I(11,3)
     -   +Dij1235I(12,3))+(-Dij1235I(1,2)+Dij1235I(2,2))*P(33)-Dij123
     -   5I(1,1)*P(71)+Dij1235I(2,1)*P(72)+Dij1235I(3,1)*P(73)+Dij123
     -   5I(5,2)*P(75)+Dij1235I(6,2)*P(76))-EijI(1,1)*P(77)+EijI(3,2)
     -   *P(531)+D01235I*P(583)+(EijI(2,1)+EijI(5,2))*P(584)+EijI(3,1
     -   )*P(585)-EijI(6,2)*P(586)-EijI(7,2)*P(587)+EijI(8,2)*P(588)+
     -   2*(D02345I-Dij2345I(1,2)+Dij2345I(4,2)+2*(p2sq*EijI(8,3)+Eij
     -   I(11,2))+6*(EijI(21,3)-EijI(22,3))-EijI(2,3)*P(5)+EijI(9,2)*
     -   P(23)+(-EijI(5,3)+EijI(6,3))*P(33)+(-EijI(10,3)+EijI(18,3)-E
     -   ijI(19,3)+EijI(20,3))*P(101)+(EijI(11,3)-EijI(12,3))*P(198)+
     -   EijI(9,3)*P(505)-EijI(17,3)*P(589))
       F(206)=DCMPLX(FR(206),FI(206))
       P(590) = -4*p3sq-p5sq+3*s45+2*P(29)
       P(591) = s15-s23+3*P(7)
       P(592) = p2sq-s12+s15-s23+3*P(7)
       P(593) = p2sq-p3sq-p5sq+s12+s34+s45
       P(594) = p2sq-p5sq+2*s12+s34
       FR(207) = D02345R+3*(Dij2345R(1,1)-Dij2345R(3,1))+EijR(3,1)*P(1
     -   0)-(EijR(2,2)-2*(EijR(8,3)-EijR(17,3)-EijR(18,3)+EijR(19,3))
     -   )*P(33)+Is34*(C0235R-Cij235R(1,1)+Cij235R(2,1)-4*(Dij1235R(7
     -   ,2)-Dij1235R(12,3)+Dij1235R(13,3))+2*(Cij135R(1,2)+(D01235R+
     -   Dij1235R(1,1))*P(3))-Dij1235R(2,2)*P(5)+Is12*(4*Cij135R(4,2)
     -   +C0135R*P(26)+Cij135R(1,1)*P(27)-Cij135R(2,1)*P(28)+2*(B035R
     -   -Cij135R(3,2)*P(29)))+(Dij1235R(2,1)+2*(Dij1235R(4,2)-Dij123
     -   5R(5,2)))*P(33)-Dij1235R(3,1)*P(34)+Dij1235R(3,2)*P(38)+Dij1
     -   235R(6,2)*P(42))+EijR(4,1)*P(110)-EijR(2,1)*P(111)+(EijR(8,2
     -   )-EijR(10,2))*P(509)+Is12*(3*(C0345R-Cij145R(2,1))+Cij345R(2
     -   ,1)-8*Dij1345R(7,2)+4*(Dij1345R(11,3)-Dij1345R(13,3))-Dij134
     -   5R(2,1)*P(16)+Is45*(2*B012R+4*Cij145R(4,2)+C0145R*P(95)+Cij1
     -   45R(1,1)*P(96))+D01345R*P(97)+Dij1345R(3,2)*P(101)+2*(Cij135
     -   R(2,2)+Cij145R(1,2)-Cij145R(3,2)-p3sq*Dij1345R(2,2)+Dij1345R
     -   (4,2)*P(198))+Dij1345R(1,1)*P(497)-Dij1345R(3,1)*P(498)-Dij1
     -   345R(6,2)*P(501)+Dij1345R(1,2)*P(579)+Dij1345R(5,2)*P(590))-
     -   EijR(4,2)*P(591)+EijR(9,2)*P(592)-2*(Dij2345R(1,2)-Dij2345R(
     -   4,2)-Dij2345R(5,2)+Dij2345R(6,2)+8*EijR(22,3)+2*(EijR(11,2)-
     -   EijR(23,3))-6*EijR(24,3)-(EE0R+EijR(1,1))*P(3)+EijR(2,3)*P(5
     -   )-(EijR(5,2)-EijR(7,2))*P(50)-EijR(15,3)*P(101)-(-EijR(12,3)
     -   +EijR(13,3))*P(198)-EijR(16,3)*P(312)-EijR(9,3)*P(505)-EijR(
     -   10,3)*P(593)+EijR(20,3)*P(594))
       FI(207) = D02345I+3*(Dij2345I(1,1)-Dij2345I(3,1))+EijI(3,1)*P(1
     -   0)-(EijI(2,2)-2*(EijI(8,3)-EijI(17,3)-EijI(18,3)+EijI(19,3))
     -   )*P(33)+Is34*(C0235I-Cij235I(1,1)+Cij235I(2,1)-4*(Dij1235I(7
     -   ,2)-Dij1235I(12,3)+Dij1235I(13,3))+2*(Cij135I(1,2)+(D01235I+
     -   Dij1235I(1,1))*P(3))-Dij1235I(2,2)*P(5)+Is12*(4*Cij135I(4,2)
     -   +C0135I*P(26)+Cij135I(1,1)*P(27)-Cij135I(2,1)*P(28)+2*(B035I
     -   -Cij135I(3,2)*P(29)))+(Dij1235I(2,1)+2*(Dij1235I(4,2)-Dij123
     -   5I(5,2)))*P(33)-Dij1235I(3,1)*P(34)+Dij1235I(3,2)*P(38)+Dij1
     -   235I(6,2)*P(42))+EijI(4,1)*P(110)-EijI(2,1)*P(111)+(EijI(8,2
     -   )-EijI(10,2))*P(509)+Is12*(3*(C0345I-Cij145I(2,1))+Cij345I(2
     -   ,1)-8*Dij1345I(7,2)+4*(Dij1345I(11,3)-Dij1345I(13,3))-Dij134
     -   5I(2,1)*P(16)+Is45*(2*B012I+4*Cij145I(4,2)+C0145I*P(95)+Cij1
     -   45I(1,1)*P(96))+D01345I*P(97)+Dij1345I(3,2)*P(101)+2*(Cij135
     -   I(2,2)+Cij145I(1,2)-Cij145I(3,2)-p3sq*Dij1345I(2,2)+Dij1345I
     -   (4,2)*P(198))+Dij1345I(1,1)*P(497)-Dij1345I(3,1)*P(498)-Dij1
     -   345I(6,2)*P(501)+Dij1345I(1,2)*P(579)+Dij1345I(5,2)*P(590))-
     -   EijI(4,2)*P(591)+EijI(9,2)*P(592)-2*(Dij2345I(1,2)-Dij2345I(
     -   4,2)-Dij2345I(5,2)+Dij2345I(6,2)+8*EijI(22,3)+2*(EijI(11,2)-
     -   EijI(23,3))-6*EijI(24,3)-(EE0I+EijI(1,1))*P(3)+EijI(2,3)*P(5
     -   )-(EijI(5,2)-EijI(7,2))*P(50)-EijI(15,3)*P(101)-(-EijI(12,3)
     -   +EijI(13,3))*P(198)-EijI(16,3)*P(312)-EijI(9,3)*P(505)-EijI(
     -   10,3)*P(593)+EijI(20,3)*P(594))
       F(207)=DCMPLX(FR(207),FI(207))
       P(595) = p2sq-s12-s23+4*s45
       P(596) = 3*P(33)+2*P(50)
       P(597) = 4*p2sq+3*p3sq+5*s12-s23
       P(598) = s15-s23-3*P(14)
       FR(208) = 3*Dij2345R(2,1)-22*EijR(11,2)-16*EijR(22,3)+4*(D02345
     -   R+EijR(23,3)+EijR(4,1)*P(7))+(-3*EijR(6,2)+2*(EijR(8,3)-EijR
     -   (17,3)))*P(33)+(3*EijR(10,2)-2*(EijR(10,3)-EijR(20,3)))*P(10
     -   1)-EijR(2,1)*P(111)+(EE0R+EijR(1,1))*P(127)+Is34*(-Cij125R(2
     -   ,1)-Cij235R(1,1)+p2sq*Dij1235R(2,1)-6*Dij1235R(7,2)+4*Dij123
     -   5R(12,3)-Dij1235R(2,2)*P(5)+2*(C0235R-Cij125R(2,2)+Cij135R(1
     -   ,2)+Dij1235R(4,2)*P(33))+Is12*(2*B035R+4*Cij135R(4,2)+C0135R
     -   *P(128)+Cij135R(1,1)*P(130))+(D01235R+Dij1235R(1,1))*P(133)+
     -   Dij1235R(3,1)*P(138)+Dij1235R(6,2)*P(142))-(3*EijR(3,2)+2*Ei
     -   jR(12,3))*P(198)-2*(Dij2345R(1,2)-Dij2345R(4,2)+EijR(2,3)*P(
     -   5)+EijR(2,2)*P(143)-EijR(9,3)*P(505))+Is12*(-10*Dij1345R(7,2
     -   )+4*(C0345R+Dij1345R(11,3))+Dij1345R(3,1)*P(7)+Dij1345R(6,2)
     -   *P(14)+D01345R*P(32)-2*(Cij135R(2,1)+Cij135R(3,2)-Cij145R(1,
     -   2)+p3sq*Dij1345R(2,2)+Dij1345R(2,1)*P(16)-Is45*(B012R+C0145R
     -   *p5sq+2*Cij145R(4,2)+Cij145R(1,1)*P(129))-Dij1345R(4,2)*P(19
     -   8))-Dij1345R(5,2)*P(292)+Dij1345R(1,1)*P(507)+Dij1345R(1,2)*
     -   P(579))+EijR(3,1)*P(595)+EijR(5,2)*P(596)+EijR(8,2)*P(597)+E
     -   ijR(9,2)*P(598)
       FI(208) = 3*Dij2345I(2,1)-22*EijI(11,2)-16*EijI(22,3)+4*(D02345
     -   I+EijI(23,3)+EijI(4,1)*P(7))+(-3*EijI(6,2)+2*(EijI(8,3)-EijI
     -   (17,3)))*P(33)+(3*EijI(10,2)-2*(EijI(10,3)-EijI(20,3)))*P(10
     -   1)-EijI(2,1)*P(111)+(EE0I+EijI(1,1))*P(127)+Is34*(-Cij125I(2
     -   ,1)-Cij235I(1,1)+p2sq*Dij1235I(2,1)-6*Dij1235I(7,2)+4*Dij123
     -   5I(12,3)-Dij1235I(2,2)*P(5)+2*(C0235I-Cij125I(2,2)+Cij135I(1
     -   ,2)+Dij1235I(4,2)*P(33))+Is12*(2*B035I+4*Cij135I(4,2)+C0135I
     -   *P(128)+Cij135I(1,1)*P(130))+(D01235I+Dij1235I(1,1))*P(133)+
     -   Dij1235I(3,1)*P(138)+Dij1235I(6,2)*P(142))-(3*EijI(3,2)+2*Ei
     -   jI(12,3))*P(198)-2*(Dij2345I(1,2)-Dij2345I(4,2)+EijI(2,3)*P(
     -   5)+EijI(2,2)*P(143)-EijI(9,3)*P(505))+Is12*(-10*Dij1345I(7,2
     -   )+4*(C0345I+Dij1345I(11,3))+Dij1345I(3,1)*P(7)+Dij1345I(6,2)
     -   *P(14)+D01345I*P(32)-2*(Cij135I(2,1)+Cij135I(3,2)-Cij145I(1,
     -   2)+p3sq*Dij1345I(2,2)+Dij1345I(2,1)*P(16)-Is45*(B012I+C0145I
     -   *p5sq+2*Cij145I(4,2)+Cij145I(1,1)*P(129))-Dij1345I(4,2)*P(19
     -   8))-Dij1345I(5,2)*P(292)+Dij1345I(1,1)*P(507)+Dij1345I(1,2)*
     -   P(579))+EijI(3,1)*P(595)+EijI(5,2)*P(596)+EijI(8,2)*P(597)+E
     -   ijI(9,2)*P(598)
       F(208)=DCMPLX(FR(208),FI(208))
       P(599) = -3*p5sq+s34+2*P(161)
       P(600) = -s34+s45-2*P(161)
       P(601) = -p2sq+s12
       P(602) = 2*s15-s45
       P(603) = -p2sq+p3sq+s23
       P(604) = -p3sq-s12-s45+2*P(430)
       P(605) = p3sq-p5sq-s12-s15+s45
       P(606) = p3sq-s12+s15+s23-s45
       FR(209) = Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(3,1)-Dij2345R(4,
     -   2)-Dij2345R(5,2)+Dij2345R(6,2)+6*EijR(11,2)-4*(EijR(22,3)-Ei
     -   jR(24,3))+2*(EijR(2,1)+EijR(5,2))*P(3)+EijR(4,2)*P(7)+EijR(1
     -   3,3)*P(10)+(-EijR(15,3)+EijR(16,3))*P(21)-EijR(2,2)*P(33)-Ei
     -   jR(6,2)*P(50)+EijR(9,3)*P(245)+EijR(12,3)*P(251)-EijR(20,3)*
     -   P(332)-EijR(7,2)*P(358)+EijR(4,1)*P(387)+EijR(10,3)*P(488)+I
     -   s12*(-C0345R-Cij345R(2,1)+12*Dij1345R(7,2)+4*(Dij1345R(11,3)
     -   -Dij1345R(13,3))+Dij1345R(3,2)*P(7)+2*(p3sq*Dij1345R(2,2)+Di
     -   j1345R(1,1)*P(35))+Dij1345R(3,1)*P(387)-Dij1345R(2,1)*P(481)
     -   +Dij1345R(1,2)*P(482)-Dij1345R(4,2)*P(483)+Dij1345R(5,2)*P(5
     -   99)+Dij1345R(6,2)*P(600))+EijR(2,3)*P(601)+EijR(3,1)*P(602)+
     -   EijR(3,2)*P(603)+EijR(8,2)*P(604)+EijR(9,2)*P(605)-EijR(10,2
     -   )*P(606)
       FI(209) = Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(3,1)-Dij2345I(4,
     -   2)-Dij2345I(5,2)+Dij2345I(6,2)+6*EijI(11,2)-4*(EijI(22,3)-Ei
     -   jI(24,3))+2*(EijI(2,1)+EijI(5,2))*P(3)+EijI(4,2)*P(7)+EijI(1
     -   3,3)*P(10)+(-EijI(15,3)+EijI(16,3))*P(21)-EijI(2,2)*P(33)-Ei
     -   jI(6,2)*P(50)+EijI(9,3)*P(245)+EijI(12,3)*P(251)-EijI(20,3)*
     -   P(332)-EijI(7,2)*P(358)+EijI(4,1)*P(387)+EijI(10,3)*P(488)+I
     -   s12*(-C0345I-Cij345I(2,1)+12*Dij1345I(7,2)+4*(Dij1345I(11,3)
     -   -Dij1345I(13,3))+Dij1345I(3,2)*P(7)+2*(p3sq*Dij1345I(2,2)+Di
     -   j1345I(1,1)*P(35))+Dij1345I(3,1)*P(387)-Dij1345I(2,1)*P(481)
     -   +Dij1345I(1,2)*P(482)-Dij1345I(4,2)*P(483)+Dij1345I(5,2)*P(5
     -   99)+Dij1345I(6,2)*P(600))+EijI(2,3)*P(601)+EijI(3,1)*P(602)+
     -   EijI(3,2)*P(603)+EijI(8,2)*P(604)+EijI(9,2)*P(605)-EijI(10,2
     -   )*P(606)
       F(209)=DCMPLX(FR(209),FI(209))
       P(607) = p2sq-s12+2*s15-s23-s45
       P(608) = -s23-s45+2*P(434)
       P(609) = -3*P(33)+2*P(50)
       FR(210) = Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(2,1)+Dij2345R(2,
     -   2)-4*(EijR(22,3)-EijR(23,3))-2*(Dij2345R(4,2)-(EijR(2,1)-Eij
     -   R(3,1)+EijR(5,2)-EijR(6,2))*P(3))+EijR(3,3)*P(10)+(EijR(10,3
     -   )+EijR(13,3)-2*EijR(20,3))*P(21)+(-EijR(9,2)+EijR(10,2))*P(2
     -   3)-EijR(2,2)*P(33)+EijR(9,3)*P(302)+Is12*(-C0345R-Cij345R(1,
     -   1)+6*Dij1345R(7,2)+4*(Dij1345R(11,3)-Dij1345R(12,3))+2*(Dij1
     -   345R(1,1)-Dij1345R(2,1))*P(35)+Dij1345R(1,2)*P(482)+Dij1345R
     -   (2,2)*P(489)-Dij1345R(4,2)*P(490)+(Dij1345R(5,2)-Dij1345R(6,
     -   2))*P(491))+EijR(2,3)*P(601)-EijR(3,2)*P(607)+EijR(8,2)*P(60
     -   8)+EijR(12,3)*P(609)
       FI(210) = Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(2,1)+Dij2345I(2,
     -   2)-4*(EijI(22,3)-EijI(23,3))-2*(Dij2345I(4,2)-(EijI(2,1)-Eij
     -   I(3,1)+EijI(5,2)-EijI(6,2))*P(3))+EijI(3,3)*P(10)+(EijI(10,3
     -   )+EijI(13,3)-2*EijI(20,3))*P(21)+(-EijI(9,2)+EijI(10,2))*P(2
     -   3)-EijI(2,2)*P(33)+EijI(9,3)*P(302)+Is12*(-C0345I-Cij345I(1,
     -   1)+6*Dij1345I(7,2)+4*(Dij1345I(11,3)-Dij1345I(12,3))+2*(Dij1
     -   345I(1,1)-Dij1345I(2,1))*P(35)+Dij1345I(1,2)*P(482)+Dij1345I
     -   (2,2)*P(489)-Dij1345I(4,2)*P(490)+(Dij1345I(5,2)-Dij1345I(6,
     -   2))*P(491))+EijI(2,3)*P(601)-EijI(3,2)*P(607)+EijI(8,2)*P(60
     -   8)+EijI(12,3)*P(609)
       F(210)=DCMPLX(FR(210),FI(210))
       P(610) = s23-s45-2*P(33)
       FR(211) = D02345R+Dij2345R(1,2)-Dij2345R(2,1)-Dij2345R(4,2)+4*E
     -   ijR(21,3)-6*(EijR(11,2)+EijR(22,3))+2*(Dij2345R(1,1)+EijR(23
     -   ,3)-(EijR(1,1)+EijR(1,2)-EijR(2,1))*P(3))+EijR(11,3)*P(10)+(
     -   EijR(10,3)-EijR(18,3)+EijR(19,3)-EijR(20,3))*P(21)+(EijR(7,2
     -   )-EijR(9,2))*P(23)+(-EijR(2,2)+EijR(8,3))*P(33)+EijR(5,2)*P(
     -   133)+EijR(9,3)*P(245)+EijR(12,3)*P(251)+EijR(2,3)*P(601)+(-E
     -   ijR(6,2)+EijR(8,2))*P(607)+EijR(17,3)*P(610)
       FI(211) = D02345I+Dij2345I(1,2)-Dij2345I(2,1)-Dij2345I(4,2)+4*E
     -   ijI(21,3)-6*(EijI(11,2)+EijI(22,3))+2*(Dij2345I(1,1)+EijI(23
     -   ,3)-(EijI(1,1)+EijI(1,2)-EijI(2,1))*P(3))+EijI(11,3)*P(10)+(
     -   EijI(10,3)-EijI(18,3)+EijI(19,3)-EijI(20,3))*P(21)+(EijI(7,2
     -   )-EijI(9,2))*P(23)+(-EijI(2,2)+EijI(8,3))*P(33)+EijI(5,2)*P(
     -   133)+EijI(9,3)*P(245)+EijI(12,3)*P(251)+EijI(2,3)*P(601)+(-E
     -   ijI(6,2)+EijI(8,2))*P(607)+EijI(17,3)*P(610)
       F(211)=DCMPLX(FR(211),FI(211))
       P(611) = p2sq-4*s12
       P(612) = -p2sq+s15+s34-3*P(137)
       P(613) = 2*s15-3*s45
       P(614) = p2sq-3*p3sq-s23-2*P(16)
       P(615) = 3*p3sq+5*s12-2*s15-s45
       P(616) = p3sq+p5sq-P(30)
       P(617) = -s15+s23+2*P(616)
       FR(212) = 3*Dij2345R(1,1)+Dij2345R(1,2)-Dij2345R(4,2)+18*EijR(1
     -   1,2)-4*EijR(22,3)+EijR(10,3)*P(21)+EijR(20,3)*P(115)+(EijR(6
     -   ,2)+EijR(9,3))*P(245)+EijR(12,3)*P(251)-2*(D02345R+Dij2345R(
     -   2,1)-EijR(2,1)*P(3)+EijR(4,1)*P(7)+EijR(5,2)*P(296)+Is15s34*
     -   (B025R+2*Cij125R(4,2)+C0125R*P(3)-Cij125R(1,1)*P(407)))+EijR
     -   (9,2)*P(512)+Is12*(-3*C0345R-Cij135R(2,1)+16*Dij1345R(7,2)+4
     -   *Dij1345R(11,3)-Dij1345R(6,2)*P(14)+2*(p3sq*Dij1345R(2,2)-Di
     -   j1345R(3,1)*P(7)+Dij1345R(1,1)*P(39))+Dij1345R(1,2)*P(482)-D
     -   ij1345R(4,2)*P(483)+Dij1345R(5,2)*P(491)+Dij1345R(2,1)*P(515
     -   ))+EijR(2,2)*P(517)+EijR(2,3)*P(601)+Is34*(-5*C0235R-Cij125R
     -   (2,1)+14*Dij1235R(7,2)+2*(p2sq*Dij1235R(2,2)-(D01235R+Dij123
     -   5R(1,1))*P(3))-Dij1235R(4,2)*P(33)-Dij1235R(6,2)*P(436)+Is12
     -   *(-2*B035R-4*Cij135R(4,2)+C0135R*P(513)-Cij135R(1,1)*P(514))
     -   +Dij1235R(2,1)*P(611)+Dij1235R(3,1)*P(612))+EijR(3,1)*P(613)
     -   -EijR(3,2)*P(614)-EijR(8,2)*P(615)-EijR(10,2)*P(617)
       FI(212) = 3*Dij2345I(1,1)+Dij2345I(1,2)-Dij2345I(4,2)+18*EijI(1
     -   1,2)-4*EijI(22,3)+EijI(10,3)*P(21)+EijI(20,3)*P(115)+(EijI(6
     -   ,2)+EijI(9,3))*P(245)+EijI(12,3)*P(251)-2*(D02345I+Dij2345I(
     -   2,1)-EijI(2,1)*P(3)+EijI(4,1)*P(7)+EijI(5,2)*P(296)+Is15s34*
     -   (B025I+2*Cij125I(4,2)+C0125I*P(3)-Cij125I(1,1)*P(407)))+EijI
     -   (9,2)*P(512)+Is12*(-3*C0345I-Cij135I(2,1)+16*Dij1345I(7,2)+4
     -   *Dij1345I(11,3)-Dij1345I(6,2)*P(14)+2*(p3sq*Dij1345I(2,2)-Di
     -   j1345I(3,1)*P(7)+Dij1345I(1,1)*P(39))+Dij1345I(1,2)*P(482)-D
     -   ij1345I(4,2)*P(483)+Dij1345I(5,2)*P(491)+Dij1345I(2,1)*P(515
     -   ))+EijI(2,2)*P(517)+EijI(2,3)*P(601)+Is34*(-5*C0235I-Cij125I
     -   (2,1)+14*Dij1235I(7,2)+2*(p2sq*Dij1235I(2,2)-(D01235I+Dij123
     -   5I(1,1))*P(3))-Dij1235I(4,2)*P(33)-Dij1235I(6,2)*P(436)+Is12
     -   *(-2*B035I-4*Cij135I(4,2)+C0135I*P(513)-Cij135I(1,1)*P(514))
     -   +Dij1235I(2,1)*P(611)+Dij1235I(3,1)*P(612))+EijI(3,1)*P(613)
     -   -EijI(3,2)*P(614)-EijI(8,2)*P(615)-EijI(10,2)*P(617)
       F(212)=DCMPLX(FR(212),FI(212))
       P(618) = -p3sq+s12+s23
       P(619) = p2sq+2*p5sq-s12-4*s15+s45
       P(620) = p2sq-s12-s23
       P(621) = -p2sq+s12+s23
       P(622) = p3sq-s12-2*s23+s45
       P(623) = -s23+s45+2*P(161)
       P(624) = p3sq-s12-s23+s45
       FR(213) = -D02345R-Dij2345R(1,2)+Dij2345R(2,1)+Dij2345R(4,2)+6*
     -   EijR(11,2)-4*(EijR(21,3)-EijR(22,3))+(-EijR(10,3)+EijR(18,3)
     -   -EijR(19,3)+EijR(20,3))*P(7)+(-EijR(7,2)+EijR(9,2))*P(23)+(-
     -   EijR(5,3)+EijR(6,3))*P(50)+Is34*(-C0235R-Cij235R(1,1)+6*Dij1
     -   235R(7,2)+4*(Dij1235R(11,3)-Dij1235R(12,3))+2*(Dij1235R(1,1)
     -   +Dij1235R(1,2)-Dij1235R(2,1))*P(3)+Dij1235R(2,2)*P(351)-Dij1
     -   235R(4,2)*P(352)-(Dij1235R(5,2)-Dij1235R(6,2))*P(353))+EijR(
     -   2,2)*P(357)+(-EijR(11,3)+EijR(12,3))*P(531)-EijR(2,3)*P(618)
     -   -EijR(5,2)*P(619)+EijR(6,2)*P(620)+EijR(8,2)*P(621)-EijR(8,3
     -   )*P(622)-EijR(9,3)*P(623)-2*(Dij2345R(1,1)-(EijR(1,1)+EijR(1
     -   ,2)-EijR(2,1))*P(3)-EijR(17,3)*P(624))
       FI(213) = -D02345I-Dij2345I(1,2)+Dij2345I(2,1)+Dij2345I(4,2)+6*
     -   EijI(11,2)-4*(EijI(21,3)-EijI(22,3))+(-EijI(10,3)+EijI(18,3)
     -   -EijI(19,3)+EijI(20,3))*P(7)+(-EijI(7,2)+EijI(9,2))*P(23)+(-
     -   EijI(5,3)+EijI(6,3))*P(50)+Is34*(-C0235I-Cij235I(1,1)+6*Dij1
     -   235I(7,2)+4*(Dij1235I(11,3)-Dij1235I(12,3))+2*(Dij1235I(1,1)
     -   +Dij1235I(1,2)-Dij1235I(2,1))*P(3)+Dij1235I(2,2)*P(351)-Dij1
     -   235I(4,2)*P(352)-(Dij1235I(5,2)-Dij1235I(6,2))*P(353))+EijI(
     -   2,2)*P(357)+(-EijI(11,3)+EijI(12,3))*P(531)-EijI(2,3)*P(618)
     -   -EijI(5,2)*P(619)+EijI(6,2)*P(620)+EijI(8,2)*P(621)-EijI(8,3
     -   )*P(622)-EijI(9,3)*P(623)-2*(Dij2345I(1,1)-(EijI(1,1)+EijI(1
     -   ,2)-EijI(2,1))*P(3)-EijI(17,3)*P(624))
       F(213)=DCMPLX(FR(213),FI(213))
       P(625) = p3sq-s12-2*s15-s23+s45
       P(626) = p3sq-p5sq-s12-s23-2*P(180)
       P(627) = p3sq+p5sq-s12-s23-s45
       P(628) = p5sq-s23+2*P(161)
       FR(214) = -Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(3,1)+Dij2345R(4
     -   ,2)+Dij2345R(5,2)-Dij2345R(6,2)-6*EijR(11,2)+2*(EijR(22,3)+E
     -   ijR(23,3))-4*EijR(24,3)-(EijR(4,2)-EijR(15,3)+EijR(16,3))*P(
     -   7)+(EijR(8,3)-EijR(17,3)-EijR(18,3)+EijR(19,3))*P(50)+(-EijR
     -   (8,2)+EijR(10,2))*P(161)+EijR(4,1)*P(185)+(-EijR(5,2)+EijR(7
     -   ,2))*P(358)+Is34*(-Cij235R(1,1)+Cij235R(2,1)-6*Dij1235R(7,2)
     -   -4*(Dij1235R(12,3)-Dij1235R(13,3))-Dij1235R(3,2)*P(73)+(-Dij
     -   1235R(2,1)+Dij1235R(3,1))*P(372)-Dij1235R(2,2)*P(373)+Dij123
     -   5R(5,2)*P(374)+Dij1235R(4,2)*P(375)+Dij1235R(6,2)*P(377))+Ei
     -   jR(2,1)*P(387)+(EijR(12,3)-EijR(13,3))*P(531)-EijR(2,3)*P(61
     -   8)-EijR(9,3)*P(623)+EijR(2,2)*P(625)-EijR(9,2)*P(626)-EijR(1
     -   0,3)*P(627)+EijR(20,3)*P(628)
       FI(214) = -Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(3,1)+Dij2345I(4
     -   ,2)+Dij2345I(5,2)-Dij2345I(6,2)-6*EijI(11,2)+2*(EijI(22,3)+E
     -   ijI(23,3))-4*EijI(24,3)-(EijI(4,2)-EijI(15,3)+EijI(16,3))*P(
     -   7)+(EijI(8,3)-EijI(17,3)-EijI(18,3)+EijI(19,3))*P(50)+(-EijI
     -   (8,2)+EijI(10,2))*P(161)+EijI(4,1)*P(185)+(-EijI(5,2)+EijI(7
     -   ,2))*P(358)+Is34*(-Cij235I(1,1)+Cij235I(2,1)-6*Dij1235I(7,2)
     -   -4*(Dij1235I(12,3)-Dij1235I(13,3))-Dij1235I(3,2)*P(73)+(-Dij
     -   1235I(2,1)+Dij1235I(3,1))*P(372)-Dij1235I(2,2)*P(373)+Dij123
     -   5I(5,2)*P(374)+Dij1235I(4,2)*P(375)+Dij1235I(6,2)*P(377))+Ei
     -   jI(2,1)*P(387)+(EijI(12,3)-EijI(13,3))*P(531)-EijI(2,3)*P(61
     -   8)-EijI(9,3)*P(623)+EijI(2,2)*P(625)-EijI(9,2)*P(626)-EijI(1
     -   0,3)*P(627)+EijI(20,3)*P(628)
       F(214)=DCMPLX(FR(214),FI(214))
       P(629) = p3sq-s12-s15
       P(630) = -s23+s45+2*P(629)
       P(631) = -2*s23+s45+3*P(161)
       P(632) = -s23+2*s45+3*P(161)
       FR(215) = -Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(2,1)-Dij2345R(2
     -   ,2)+2*Dij2345R(4,2)+4*(EijR(22,3)-EijR(23,3))-(-EijR(9,2)+Ei
     -   jR(10,2)+EijR(10,3)+EijR(13,3)-2*EijR(20,3))*P(7)+(EijR(8,3)
     -   +EijR(11,3)-2*EijR(17,3))*P(50)+EijR(3,2)*P(161)+EijR(3,1)*P
     -   (185)-EijR(5,2)*P(358)+EijR(6,2)*P(358)+Is34*(-Cij235R(1,1)+
     -   Cij235R(2,1)-6*Dij1235R(7,2)-4*(Dij1235R(12,3)-Dij1235R(13,3
     -   ))-Dij1235R(3,2)*P(73)+(-Dij1235R(2,1)+Dij1235R(3,1))*P(372)
     -   -Dij1235R(2,2)*P(373)+Dij1235R(5,2)*P(374)+Dij1235R(4,2)*P(3
     -   75)+Dij1235R(6,2)*P(377))+EijR(2,1)*P(387)-EijR(3,3)*P(531)-
     -   EijR(2,3)*P(618)+EijR(2,2)*P(625)-EijR(8,2)*P(630)-EijR(9,3)
     -   *P(631)+EijR(12,3)*P(632)
       FI(215) = -Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(2,1)-Dij2345I(2
     -   ,2)+2*Dij2345I(4,2)+4*(EijI(22,3)-EijI(23,3))-(-EijI(9,2)+Ei
     -   jI(10,2)+EijI(10,3)+EijI(13,3)-2*EijI(20,3))*P(7)+(EijI(8,3)
     -   +EijI(11,3)-2*EijI(17,3))*P(50)+EijI(3,2)*P(161)+EijI(3,1)*P
     -   (185)-EijI(5,2)*P(358)+EijI(6,2)*P(358)+Is34*(-Cij235I(1,1)+
     -   Cij235I(2,1)-6*Dij1235I(7,2)-4*(Dij1235I(12,3)-Dij1235I(13,3
     -   ))-Dij1235I(3,2)*P(73)+(-Dij1235I(2,1)+Dij1235I(3,1))*P(372)
     -   -Dij1235I(2,2)*P(373)+Dij1235I(5,2)*P(374)+Dij1235I(4,2)*P(3
     -   75)+Dij1235I(6,2)*P(377))+EijI(2,1)*P(387)-EijI(3,3)*P(531)-
     -   EijI(2,3)*P(618)+EijI(2,2)*P(625)-EijI(8,2)*P(630)-EijI(9,3)
     -   *P(631)+EijI(12,3)*P(632)
       F(215)=DCMPLX(FR(215),FI(215))
       P(633) = -p3sq+3*s12+s23-s45+2*P(430)
       P(634) = s23-s45+2*P(296)
       P(635) = p3sq+5*s12+2*P(414)
       FR(216) = -3*Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(4,2)+(2*EijR(
     -   4,1)-EijR(10,3)+EijR(20,3))*P(7)+(EijR(8,3)-EijR(17,3))*P(50
     -   )+2*(D02345R+Dij2345R(2,1)+s45*EijR(3,1)-9*EijR(11,2)+EijR(2
     -   2,3)+EijR(23,3)+Is15s34*(B025R+2*Cij125R(4,2)+(C0125R+Cij125
     -   R(1,1))*P(3))-EijR(6,2)*P(33)+EijR(10,2)*P(101)-EijR(3,2)*P(
     -   198))-EijR(9,2)*P(292)+EijR(2,1)*P(387)+Is34*(Cij125R(2,1)-C
     -   ij235R(1,1)-10*Dij1235R(7,2)-4*Dij1235R(12,3)+Dij1235R(6,2)*
     -   P(73)+2*(C0235R+Dij1235R(3,1)*P(137))-Dij1235R(2,2)*P(373)+D
     -   ij1235R(4,2)*P(375)-Dij1235R(2,1)*P(408))+EijR(12,3)*P(531)-
     -   EijR(2,3)*P(618)-EijR(9,3)*P(623)-EijR(2,2)*P(633)+EijR(5,2)
     -   *P(634)+EijR(8,2)*P(635)
       FI(216) = -3*Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(4,2)+(2*EijI(
     -   4,1)-EijI(10,3)+EijI(20,3))*P(7)+(EijI(8,3)-EijI(17,3))*P(50
     -   )+2*(D02345I+Dij2345I(2,1)+s45*EijI(3,1)-9*EijI(11,2)+EijI(2
     -   2,3)+EijI(23,3)+Is15s34*(B025I+2*Cij125I(4,2)+(C0125I+Cij125
     -   I(1,1))*P(3))-EijI(6,2)*P(33)+EijI(10,2)*P(101)-EijI(3,2)*P(
     -   198))-EijI(9,2)*P(292)+EijI(2,1)*P(387)+Is34*(Cij125I(2,1)-C
     -   ij235I(1,1)-10*Dij1235I(7,2)-4*Dij1235I(12,3)+Dij1235I(6,2)*
     -   P(73)+2*(C0235I+Dij1235I(3,1)*P(137))-Dij1235I(2,2)*P(373)+D
     -   ij1235I(4,2)*P(375)-Dij1235I(2,1)*P(408))+EijI(12,3)*P(531)-
     -   EijI(2,3)*P(618)-EijI(9,3)*P(623)-EijI(2,2)*P(633)+EijI(5,2)
     -   *P(634)+EijI(8,2)*P(635)
       F(216)=DCMPLX(FR(216),FI(216))
       FR(217) = 4*(-EijR(6,2)+EijR(6,3)+EijR(8,2)+EijR(9,3)-2*EijR(17
     -   ,3))
       FI(217) = 4*(-EijI(6,2)+EijI(6,3)+EijI(8,2)+EijI(9,3)-2*EijI(17
     -   ,3))
       F(217)=DCMPLX(FR(217),FI(217))
       FR(218) = -4*(EijR(2,1)-EijR(3,1)+EijR(3,2)+EijR(5,2)-EijR(6,2)
     -   -EijR(8,2)-EijR(9,3)-EijR(11,3)+EijR(12,3)+EijR(17,3))
       FI(218) = -4*(EijI(2,1)-EijI(3,1)+EijI(3,2)+EijI(5,2)-EijI(6,2)
     -   -EijI(8,2)-EijI(9,3)-EijI(11,3)+EijI(12,3)+EijI(17,3))
       F(218)=DCMPLX(FR(218),FI(218))
       FR(219) = -4*(EijR(2,1)-EijR(3,1)+EijR(5,2)+EijR(6,2)-3*EijR(7,
     -   2)-2*(EijR(8,2)-EijR(9,2))-EijR(9,3)+EijR(10,2)+EijR(17,3)-E
     -   ijR(19,3)+EijR(20,3))
       FI(219) = -4*(EijI(2,1)-EijI(3,1)+EijI(5,2)+EijI(6,2)-3*EijI(7,
     -   2)-2*(EijI(8,2)-EijI(9,2))-EijI(9,3)+EijI(10,2)+EijI(17,3)-E
     -   ijI(19,3)+EijI(20,3))
       F(219)=DCMPLX(FR(219),FI(219))
       FR(220) = -4*(EijR(2,1)-EijR(3,1)+EijR(5,2)+EijR(6,2)-2*EijR(8,
     -   2)-EijR(9,3)+EijR(17,3))
       FI(220) = -4*(EijI(2,1)-EijI(3,1)+EijI(5,2)+EijI(6,2)-2*EijI(8,
     -   2)-EijI(9,3)+EijI(17,3))
       F(220)=DCMPLX(FR(220),FI(220))
       FR(221) = 4*(EijR(2,1)-EijR(3,1)+EijR(5,2)-EijR(6,2)+EijR(9,3)+
     -   EijR(11,3)-EijR(12,3)-EijR(17,3))
       FI(221) = 4*(EijI(2,1)-EijI(3,1)+EijI(5,2)-EijI(6,2)+EijI(9,3)+
     -   EijI(11,3)-EijI(12,3)-EijI(17,3))
       F(221)=DCMPLX(FR(221),FI(221))
       FR(222) = 4*Is12*(Dij1345R(2,2)-Dij1345R(2,3)-Dij1345R(4,2)-Dij
     -   1345R(4,3)+2*Dij1345R(6,3)+s12*(EijR(3,3)+EijR(9,3)-2*EijR(1
     -   2,3)))
       FI(222) = 4*Is12*(Dij1345I(2,2)-Dij1345I(2,3)-Dij1345I(4,2)-Dij
     -   1345I(4,3)+2*Dij1345I(6,3)+s12*(EijI(3,3)+EijI(9,3)-2*EijI(1
     -   2,3)))
       F(222)=DCMPLX(FR(222),FI(222))
       FR(223) = 4*Is12*(Dij1345R(2,2)-Dij1345R(4,2)-Dij1345R(4,3)-Dij
     -   1345R(5,2)+Dij1345R(6,2)+Dij1345R(6,3)-Dij1345R(8,3)+Dij1345
     -   R(10,3)+s12*(-EijR(3,2)+EijR(8,2)+EijR(9,3)-2*(EijR(9,2)-Eij
     -   R(10,2))-EijR(12,3)+EijR(13,3)-EijR(20,3)))
       FI(223) = 4*Is12*(Dij1345I(2,2)-Dij1345I(4,2)-Dij1345I(4,3)-Dij
     -   1345I(5,2)+Dij1345I(6,2)+Dij1345I(6,3)-Dij1345I(8,3)+Dij1345
     -   I(10,3)+s12*(-EijI(3,2)+EijI(8,2)+EijI(9,3)-2*(EijI(9,2)-Eij
     -   I(10,2))-EijI(12,3)+EijI(13,3)-EijI(20,3)))
       F(223)=DCMPLX(FR(223),FI(223))
       FR(224) = 4*Is12*(Dij1345R(2,2)-Dij1345R(4,2)-Dij1345R(4,3)+Dij
     -   1345R(6,3)+s12*(-EijR(3,2)+EijR(8,2)+EijR(9,3)-EijR(12,3)))
       FI(224) = 4*Is12*(Dij1345I(2,2)-Dij1345I(4,2)-Dij1345I(4,3)+Dij
     -   1345I(6,3)+s12*(-EijI(3,2)+EijI(8,2)+EijI(9,3)-EijI(12,3)))
       F(224)=DCMPLX(FR(224),FI(224))
       FR(225) = 4*(EijR(2,1)-EijR(4,1)+EijR(5,2)+EijR(6,2)-2*EijR(7,2
     -   )-EijR(8,2)+EijR(9,2)+EijR(9,3)-EijR(17,3)+EijR(19,3)-EijR(2
     -   0,3))
       FI(225) = 4*(EijI(2,1)-EijI(4,1)+EijI(5,2)+EijI(6,2)-2*EijI(7,2
     -   )-EijI(8,2)+EijI(9,2)+EijI(9,3)-EijI(17,3)+EijI(19,3)-EijI(2
     -   0,3))
       F(225)=DCMPLX(FR(225),FI(225))
       FR(226) = -4*Is12*(Dij1345R(2,1)-Dij1345R(3,1)+3*Dij1345R(4,2)+
     -   Dij1345R(4,3)-2*(Dij1345R(2,2)+Dij1345R(5,2))+Dij1345R(6,2)-
     -   Dij1345R(6,3)+Dij1345R(8,3)-Dij1345R(10,3)+s12*(-EijR(3,2)+E
     -   ijR(8,2)-EijR(9,2)-EijR(9,3)+EijR(10,2)+EijR(12,3)-EijR(13,3
     -   )+EijR(20,3)))
       FI(226) = -4*Is12*(Dij1345I(2,1)-Dij1345I(3,1)+3*Dij1345I(4,2)+
     -   Dij1345I(4,3)-2*(Dij1345I(2,2)+Dij1345I(5,2))+Dij1345I(6,2)-
     -   Dij1345I(6,3)+Dij1345I(8,3)-Dij1345I(10,3)+s12*(-EijI(3,2)+E
     -   ijI(8,2)-EijI(9,2)-EijI(9,3)+EijI(10,2)+EijI(12,3)-EijI(13,3
     -   )+EijI(20,3)))
       F(226)=DCMPLX(FR(226),FI(226))
       FR(227) = -4*Is12*(Dij1345R(3,2)+Dij1345R(4,3)-Dij1345R(5,2)+3*
     -   (Dij1345R(4,2)-Dij1345R(6,2))+Dij1345R(9,3)+2*(Dij1345R(2,1)
     -   -Dij1345R(3,1)-Dij1345R(10,3))-s12*(-EijR(3,1)+EijR(4,1)+Eij
     -   R(4,2)-EijR(9,2)+EijR(9,3)+EijR(16,3)-2*EijR(20,3)))
       FI(227) = -4*Is12*(Dij1345I(3,2)+Dij1345I(4,3)-Dij1345I(5,2)+3*
     -   (Dij1345I(4,2)-Dij1345I(6,2))+Dij1345I(9,3)+2*(Dij1345I(2,1)
     -   -Dij1345I(3,1)-Dij1345I(10,3))-s12*(-EijI(3,1)+EijI(4,1)+Eij
     -   I(4,2)-EijI(9,2)+EijI(9,3)+EijI(16,3)-2*EijI(20,3)))
       F(227)=DCMPLX(FR(227),FI(227))
       FR(228) = -4*Is12*(3*Dij1345R(4,2)+Dij1345R(4,3)+2*(Dij1345R(2,
     -   1)-Dij1345R(3,1)-Dij1345R(5,2))-Dij1345R(6,2)-Dij1345R(10,3)
     -   +s12*(EijR(3,1)-EijR(4,1)-EijR(9,2)-EijR(9,3)+EijR(10,2)+Eij
     -   R(20,3)))
       FI(228) = -4*Is12*(3*Dij1345I(4,2)+Dij1345I(4,3)+2*(Dij1345I(2,
     -   1)-Dij1345I(3,1)-Dij1345I(5,2))-Dij1345I(6,2)-Dij1345I(10,3)
     -   +s12*(EijI(3,1)-EijI(4,1)-EijI(9,2)-EijI(9,3)+EijI(10,2)+Eij
     -   I(20,3)))
       F(228)=DCMPLX(FR(228),FI(228))
       FR(229) = -2*Is12*(D01345R+Dij1345R(1,1)-3*Dij1345R(2,2)+5*Dij1
     -   345R(4,2)+2*(Dij1345R(2,1)+Dij1345R(4,3)-Dij1345R(6,3))+s12*
     -   (EijR(3,2)-EijR(8,2)-2*(EijR(9,3)-EijR(12,3))))
       FI(229) = -2*Is12*(D01345I+Dij1345I(1,1)-3*Dij1345I(2,2)+5*Dij1
     -   345I(4,2)+2*(Dij1345I(2,1)+Dij1345I(4,3)-Dij1345I(6,3))+s12*
     -   (EijI(3,2)-EijI(8,2)-2*(EijI(9,3)-EijI(12,3))))
       F(229)=DCMPLX(FR(229),FI(229))
       FR(230) = -2*Is12*(D01345R+Dij1345R(1,1)+5*Dij1345R(4,2)+3*(Dij
     -   1345R(2,1)-Dij1345R(6,2))+2*(Dij1345R(4,3)+Dij1345R(5,2)-Dij
     -   1345R(10,3))+s12*(-EijR(3,1)-3*EijR(8,2)+4*EijR(9,2)+EijR(10
     -   ,2)+2*(EijR(4,1)-EijR(9,3)+EijR(20,3))))
       FI(230) = -2*Is12*(D01345I+Dij1345I(1,1)+5*Dij1345I(4,2)+3*(Dij
     -   1345I(2,1)-Dij1345I(6,2))+2*(Dij1345I(4,3)+Dij1345I(5,2)-Dij
     -   1345I(10,3))+s12*(-EijI(3,1)-3*EijI(8,2)+4*EijI(9,2)+EijI(10
     -   ,2)+2*(EijI(4,1)-EijI(9,3)+EijI(20,3))))
       F(230)=DCMPLX(FR(230),FI(230))
       FR(231) = 2*(EE0R+EijR(1,1)-EijR(6,2)+EijR(8,2)+2*(EijR(2,1)+Ei
     -   jR(5,2)+EijR(9,3)-EijR(17,3)))
       FI(231) = 2*(EE0I+EijI(1,1)-EijI(6,2)+EijI(8,2)+2*(EijI(2,1)+Ei
     -   jI(5,2)+EijI(9,3)-EijI(17,3)))
       F(231)=DCMPLX(FR(231),FI(231))
       P(636) = 5*s12+2*P(39)
       P(637) = p3sq+2*s12-3*s45
       FR(232) = 3*D02345R-Dij2345R(1,1)-16*EijR(11,2)+4*EijR(23,3)-Ei
     -   jR(9,2)*P(18)-(EijR(2,1)-2*(EijR(5,2)-EijR(6,2)))*P(33)+2*(D
     -   ij1345R(4,2)+Dij2345R(2,1)+(EE0R+EijR(1,1))*P(3)+EijR(10,2)*
     -   P(101)-EijR(3,2)*P(198))-EijR(2,2)*P(508)+EijR(3,1)*P(509)+E
     -   ijR(8,2)*P(511)+EijR(4,1)*P(533)+Is12*(3*C0345R+Cij345R(1,1)
     -   -8*Dij1345R(7,2)-4*Dij1345R(12,3)+2*D01345R*P(35)+Dij1345R(6
     -   ,2)*P(101)+Is45*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)*P(95)+C
     -   0145R*P(185))-Dij1345R(2,2)*P(198)+Dij1345R(3,1)*P(328)+Dij1
     -   345R(1,1)*P(636)-Dij1345R(2,1)*P(637))
       FI(232) = 3*D02345I-Dij2345I(1,1)-16*EijI(11,2)+4*EijI(23,3)-Ei
     -   jI(9,2)*P(18)-(EijI(2,1)-2*(EijI(5,2)-EijI(6,2)))*P(33)+2*(D
     -   ij1345I(4,2)+Dij2345I(2,1)+(EE0I+EijI(1,1))*P(3)+EijI(10,2)*
     -   P(101)-EijI(3,2)*P(198))-EijI(2,2)*P(508)+EijI(3,1)*P(509)+E
     -   ijI(8,2)*P(511)+EijI(4,1)*P(533)+Is12*(3*C0345I+Cij345I(1,1)
     -   -8*Dij1345I(7,2)-4*Dij1345I(12,3)+2*D01345I*P(35)+Dij1345I(6
     -   ,2)*P(101)+Is45*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)*P(95)+C
     -   0145I*P(185))-Dij1345I(2,2)*P(198)+Dij1345I(3,1)*P(328)+Dij1
     -   345I(1,1)*P(636)-Dij1345I(2,1)*P(637))
       F(232)=DCMPLX(FR(232),FI(232))
       FR(233) = -2*Is12*(D01345R+Dij1345R(1,1)+3*Dij1345R(2,1)+5*Dij1
     -   345R(4,2)+2*Dij1345R(4,3)-s12*(EijR(3,1)+3*EijR(8,2)+2*EijR(
     -   9,3)))
       FI(233) = -2*Is12*(D01345I+Dij1345I(1,1)+3*Dij1345I(2,1)+5*Dij1
     -   345I(4,2)+2*Dij1345I(4,3)-s12*(EijI(3,1)+3*EijI(8,2)+2*EijI(
     -   9,3)))
       F(233)=DCMPLX(FR(233),FI(233))
       FR(234) = -4*Is34*(Dij1235R(2,2)-Dij1235R(2,3)-Dij1235R(4,2)-Di
     -   j1235R(4,3)+2*Dij1235R(6,3)+s34*(EijR(2,2)+EijR(2,3)-EijR(5,
     -   2)+EijR(5,3)-2*EijR(8,3)))
       FI(234) = -4*Is34*(Dij1235I(2,2)-Dij1235I(2,3)-Dij1235I(4,2)-Di
     -   j1235I(4,3)+2*Dij1235I(6,3)+s34*(EijI(2,2)+EijI(2,3)-EijI(5,
     -   2)+EijI(5,3)-2*EijI(8,3)))
       F(234)=DCMPLX(FR(234),FI(234))
       FR(235) = 4*Is34*(Dij1235R(2,3)-Dij1235R(5,2)+Dij1235R(6,2)-Dij
     -   1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3)+s34*(-EijR(2,3)-2*(E
     -   ijR(2,2)-EijR(5,2)+EijR(6,2)-EijR(8,2))+EijR(8,3)+EijR(9,3)-
     -   EijR(17,3)))
       FI(235) = 4*Is34*(Dij1235I(2,3)-Dij1235I(5,2)+Dij1235I(6,2)-Dij
     -   1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3)+s34*(-EijI(2,3)-2*(E
     -   ijI(2,2)-EijI(5,2)+EijI(6,2)-EijI(8,2))+EijI(8,3)+EijI(9,3)-
     -   EijI(17,3)))
       F(235)=DCMPLX(FR(235),FI(235))
       FR(236) = 4*Is34*(Dij1235R(2,3)-Dij1235R(5,2)+Dij1235R(6,2)-Dij
     -   1235R(6,3)-Dij1235R(8,3)+Dij1235R(10,3)+s34*(-EijR(2,3)-2*(E
     -   ijR(2,2)-EijR(5,2))+EijR(8,3)-3*(EijR(7,2)-EijR(9,2))+EijR(1
     -   0,3)-EijR(18,3)))
       FI(236) = 4*Is34*(Dij1235I(2,3)-Dij1235I(5,2)+Dij1235I(6,2)-Dij
     -   1235I(6,3)-Dij1235I(8,3)+Dij1235I(10,3)+s34*(-EijI(2,3)-2*(E
     -   ijI(2,2)-EijI(5,2))+EijI(8,3)-3*(EijI(7,2)-EijI(9,2))+EijI(1
     -   0,3)-EijI(18,3)))
       F(236)=DCMPLX(FR(236),FI(236))
       FR(237) = 4*Is34*(Dij1235R(2,3)-Dij1235R(6,3)+s34*(-EijR(2,3)-2
     -   *(EijR(2,2)-EijR(5,2))+EijR(8,3)))
       FI(237) = -4*Is34*(-Dij1235I(2,3)+Dij1235I(6,3)+s34*(EijI(2,3)+
     -   2*(EijI(2,2)-EijI(5,2))-EijI(8,3)))
       F(237)=DCMPLX(FR(237),FI(237))
       FR(238) = -4*Is34*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(2,3)-Di
     -   j1235R(3,1)+2*(Dij1235R(4,2)-Dij1235R(5,2))+Dij1235R(6,2)+Di
     -   j1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)+s34*(EijR(2,1)+EijR
     -   (2,3)-EijR(3,1)+EijR(5,2)-EijR(6,2)-EijR(8,3)-EijR(9,3)+EijR
     -   (17,3)))
       FI(238) = -4*Is34*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(2,3)-Di
     -   j1235I(3,1)+2*(Dij1235I(4,2)-Dij1235I(5,2))+Dij1235I(6,2)+Di
     -   j1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)+s34*(EijI(2,1)+EijI
     -   (2,3)-EijI(3,1)+EijI(5,2)-EijI(6,2)-EijI(8,3)-EijI(9,3)+EijI
     -   (17,3)))
       F(238)=DCMPLX(FR(238),FI(238))
       FR(239) = 4*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)-s34*(EijR(2,2)+EijR(2,3)+EijR(3,2)-2*(EijR(8,2)
     -   +EijR(9,3))+EijR(12,3)))
       FI(239) = 4*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)-s34*(EijI(2,2)+EijI(2,3)+EijI(3,2)-2*(EijI(8,2)
     -   +EijI(9,3))+EijI(12,3)))
       F(239)=DCMPLX(FR(239),FI(239))
       FR(240) = 4*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(-EijR(2,2)-EijR(2,3)+EijR(8,2)+EijR(9,3)+2
     -   *(EijR(9,2)-EijR(10,2))+EijR(10,3)-EijR(20,3)))
       FI(240) = 4*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(-EijI(2,2)-EijI(2,3)+EijI(8,2)+EijI(9,3)+2
     -   *(EijI(9,2)-EijI(10,2))+EijI(10,3)-EijI(20,3)))
       F(240)=DCMPLX(FR(240),FI(240))
       FR(241) = 4*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+2*(
     -   Dij1235R(2,2)-Dij1235R(6,2))-Dij1235R(8,3)-s34*EijR(2,2)-s34
     -   *EijR(2,3)+s34*EijR(8,2)+s34*EijR(9,3))
       FI(241) = 4*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+2*(
     -   Dij1235I(2,2)-Dij1235I(6,2))-Dij1235I(8,3)-s34*EijI(2,2)-s34
     -   *EijI(2,3)+s34*EijI(8,2)+s34*EijI(9,3))
       F(241)=DCMPLX(FR(241),FI(241))
       FR(242) = -4*Is34*(Dij1235R(2,1)-Dij1235R(2,2)-Dij1235R(2,3)-Di
     -   j1235R(3,1)+2*(Dij1235R(4,2)-Dij1235R(5,2))+Dij1235R(6,2)+Di
     -   j1235R(6,3)+Dij1235R(8,3)-Dij1235R(10,3)+s34*(EijR(2,1)-EijR
     -   (2,2)+EijR(2,3)-EijR(4,1)+2*(EijR(5,2)-EijR(7,2))-EijR(8,3)+
     -   EijR(9,2)-EijR(10,3)+EijR(18,3)))
       FI(242) = -4*Is34*(Dij1235I(2,1)-Dij1235I(2,2)-Dij1235I(2,3)-Di
     -   j1235I(3,1)+2*(Dij1235I(4,2)-Dij1235I(5,2))+Dij1235I(6,2)+Di
     -   j1235I(6,3)+Dij1235I(8,3)-Dij1235I(10,3)+s34*(EijI(2,1)-EijI
     -   (2,2)+EijI(2,3)-EijI(4,1)+2*(EijI(5,2)-EijI(7,2))-EijI(8,3)+
     -   EijI(9,2)-EijI(10,3)+EijI(18,3)))
       F(242)=DCMPLX(FR(242),FI(242))
       FR(243) = 4*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,1)-EijR(2,3)-EijR(3,1)+EijR(9,3)+Ei
     -   jR(10,3)-EijR(20,3)))
       FI(243) = 4*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,1)-EijI(2,3)-EijI(3,1)+EijI(9,3)+Ei
     -   jI(10,3)-EijI(20,3)))
       F(243)=DCMPLX(FR(243),FI(243))
       FR(244) = 4*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+Dij
     -   1235R(3,2)-3*Dij1235R(6,2)+2*(Dij1235R(2,2)-Dij1235R(8,3))+D
     -   ij1235R(9,3)+s34*(EijR(2,1)-EijR(2,3)-EijR(4,1)-EijR(4,2)+Ei
     -   jR(9,2)+2*EijR(10,3)-EijR(15,3)))
       FI(244) = 4*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+Dij
     -   1235I(3,2)-3*Dij1235I(6,2)+2*(Dij1235I(2,2)-Dij1235I(8,3))+D
     -   ij1235I(9,3)+s34*(EijI(2,1)-EijI(2,3)-EijI(4,1)-EijI(4,2)+Ei
     -   jI(9,2)+2*EijI(10,3)-EijI(15,3)))
       F(244)=DCMPLX(FR(244),FI(244))
       FR(245) = 4*Is34*(Dij1235R(2,1)+Dij1235R(2,3)-Dij1235R(3,1)+2*(
     -   Dij1235R(2,2)-Dij1235R(6,2))-Dij1235R(8,3)+s34*EijR(2,1)-s34
     -   *EijR(2,3)-s34*EijR(4,1)+s34*EijR(10,3))
       FI(245) = 4*Is34*(Dij1235I(2,1)+Dij1235I(2,3)-Dij1235I(3,1)+2*(
     -   Dij1235I(2,2)-Dij1235I(6,2))-Dij1235I(8,3)+s34*EijI(2,1)-s34
     -   *EijI(2,3)-s34*EijI(4,1)+s34*EijI(10,3))
       F(245)=DCMPLX(FR(245),FI(245))
       FR(246) = 2*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+2
     -   *(Dij1235R(2,3)-Dij1235R(8,3))+s34*(-EijR(2,1)+EijR(3,1)-3*(
     -   EijR(2,2)-EijR(8,2))-2*(EijR(2,3)-EijR(9,3))))
       FI(246) = 2*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+2
     -   *(Dij1235I(2,3)-Dij1235I(8,3))+s34*(-EijI(2,1)+EijI(3,1)-3*(
     -   EijI(2,2)-EijI(8,2))-2*(EijI(2,3)-EijI(9,3))))
       F(246)=DCMPLX(FR(246),FI(246))
       FR(247) = -2*Is34*(D01235R+Dij1235R(1,1)-Dij1235R(2,2)+3*Dij123
     -   5R(4,2)+2*(Dij1235R(2,1)-Dij1235R(2,3)+Dij1235R(6,3))+s34*(E
     -   E0R+EijR(1,1)+EijR(2,2)+EijR(5,2)+2*(EijR(2,1)+EijR(2,3)-Eij
     -   R(8,3))))
       FI(247) = -2*Is34*(D01235I+Dij1235I(1,1)-Dij1235I(2,2)+3*Dij123
     -   5I(4,2)+2*(Dij1235I(2,1)-Dij1235I(2,3)+Dij1235I(6,3))+s34*(E
     -   E0I+EijI(1,1)+EijI(2,2)+EijI(5,2)+2*(EijI(2,1)+EijI(2,3)-Eij
     -   I(8,3))))
       F(247)=DCMPLX(FR(247),FI(247))
       FR(248) = 2*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)-Dij1235R(6,2)+2
     -   *(Dij1235R(2,3)-Dij1235R(8,3))-s34*(EijR(2,1)+3*EijR(2,2)-5*
     -   EijR(9,2)+2*(EijR(2,3)-EijR(4,1)-EijR(10,3))))
       FI(248) = 2*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)-Dij1235I(6,2)+2
     -   *(Dij1235I(2,3)-Dij1235I(8,3))-s34*(EijI(2,1)+3*EijI(2,2)-5*
     -   EijI(9,2)+2*(EijI(2,3)-EijI(4,1)-EijI(10,3))))
       F(248)=DCMPLX(FR(248),FI(248))
       P(638) = 7*s12+2*P(39)
       P(639) = p3sq+4*P(16)
       FR(249) = -3*D02345R+Dij2345R(1,1)+16*EijR(11,2)-4*EijR(22,3)+E
     -   ijR(9,2)*P(18)+(EijR(2,1)-2*(EijR(5,2)-EijR(6,2)))*P(33)-2*(
     -   Dij2345R(2,1)+(EE0R+EijR(1,1))*P(3)-EijR(3,2)*P(198)-EijR(10
     -   ,2)*P(312))+EijR(2,2)*P(508)-EijR(3,1)*P(509)-EijR(8,2)*P(51
     -   1)+Is34*(-3*C0235R-Cij235R(1,1)+8*Dij1235R(7,2)+4*Dij1235R(1
     -   2,3)-2*(D01235R+Dij1235R(1,1))*P(3)+Dij1235R(2,2)*P(33)+Dij1
     -   235R(2,1)*P(60)-Dij1235R(6,2)*P(330)+Is12*(-2*B035R-4*Cij135
     -   R(4,2)+C0135R*P(513)-Cij135R(1,1)*P(514))-Dij1235R(3,1)*P(51
     -   6))-EijR(4,1)*P(533)+Is12*(-5*C0345R-Cij135R(2,1)+14*Dij1345
     -   R(7,2)-Dij1345R(6,2)*P(14)+2*(p3sq*Dij1345R(2,2)-D01345R*P(3
     -   5))-Dij1345R(4,2)*P(198)-Dij1345R(3,1)*P(328)+Is45*(-2*(B012
     -   R+2*Cij145R(4,2))-Cij145R(1,1)*P(95)+C0145R*P(387))-Dij1345R
     -   (1,1)*P(638)+Dij1345R(2,1)*P(639))
       FI(249) = -3*D02345I+Dij2345I(1,1)+16*EijI(11,2)-4*EijI(22,3)+E
     -   ijI(9,2)*P(18)+(EijI(2,1)-2*(EijI(5,2)-EijI(6,2)))*P(33)-2*(
     -   Dij2345I(2,1)+(EE0I+EijI(1,1))*P(3)-EijI(3,2)*P(198)-EijI(10
     -   ,2)*P(312))+EijI(2,2)*P(508)-EijI(3,1)*P(509)-EijI(8,2)*P(51
     -   1)+Is34*(-3*C0235I-Cij235I(1,1)+8*Dij1235I(7,2)+4*Dij1235I(1
     -   2,3)-2*(D01235I+Dij1235I(1,1))*P(3)+Dij1235I(2,2)*P(33)+Dij1
     -   235I(2,1)*P(60)-Dij1235I(6,2)*P(330)+Is12*(-2*B035I-4*Cij135
     -   I(4,2)+C0135I*P(513)-Cij135I(1,1)*P(514))-Dij1235I(3,1)*P(51
     -   6))-EijI(4,1)*P(533)+Is12*(-5*C0345I-Cij135I(2,1)+14*Dij1345
     -   I(7,2)-Dij1345I(6,2)*P(14)+2*(p3sq*Dij1345I(2,2)-D01345I*P(3
     -   5))-Dij1345I(4,2)*P(198)-Dij1345I(3,1)*P(328)+Is45*(-2*(B012
     -   I+2*Cij145I(4,2))-Cij145I(1,1)*P(95)+C0145I*P(387))-Dij1345I
     -   (1,1)*P(638)+Dij1345I(2,1)*P(639))
       F(249)=DCMPLX(FR(249),FI(249))
       FR(250) = 2*Is34*(Dij1235R(2,1)+3*Dij1235R(2,2)+2*Dij1235R(2,3)
     -   -s34*(EijR(2,1)+3*EijR(2,2)+2*EijR(2,3)))
       FI(250) = 2*Is34*(Dij1235I(2,1)+3*Dij1235I(2,2)+2*Dij1235I(2,3)
     -   -s34*(EijI(2,1)+3*EijI(2,2)+2*EijI(2,3)))
       F(250)=DCMPLX(FR(250),FI(250))
       FR(251) = 8*(EijR(2,3)+EijR(5,3)-EijR(6,3)-EijR(9,3)-2*(EijR(8,
     -   3)-EijR(17,3)))
       FI(251) = 8*(EijI(2,3)+EijI(5,3)-EijI(6,3)-EijI(9,3)-2*(EijI(8,
     -   3)-EijI(17,3)))
       F(251)=DCMPLX(FR(251),FI(251))
       FR(252) = 8*(EijR(2,2)+EijR(2,3)-EijR(5,2)+EijR(6,2)-EijR(8,2)-
     -   EijR(8,3)-EijR(11,3)+EijR(12,3)-2*(EijR(9,3)-EijR(17,3)))
       FI(252) = 8*(EijI(2,2)+EijI(2,3)-EijI(5,2)+EijI(6,2)-EijI(8,2)-
     -   EijI(8,3)-EijI(11,3)+EijI(12,3)-2*(EijI(9,3)-EijI(17,3)))
       F(252)=DCMPLX(FR(252),FI(252))
       FR(253) = 8*(EijR(2,2)+EijR(2,3)-EijR(5,2)+EijR(7,2)-EijR(8,3)-
     -   EijR(9,2)-EijR(9,3)-EijR(10,3)+EijR(17,3)+EijR(18,3)-EijR(19
     -   ,3)+EijR(20,3))
       FI(253) = 8*(EijI(2,2)+EijI(2,3)-EijI(5,2)+EijI(7,2)-EijI(8,3)-
     -   EijI(9,2)-EijI(9,3)-EijI(10,3)+EijI(17,3)+EijI(18,3)-EijI(19
     -   ,3)+EijI(20,3))
       F(253)=DCMPLX(FR(253),FI(253))
       FR(254) = 4*(EijR(6,2)-EijR(8,2)+2*(EijR(2,2)+EijR(2,3)-EijR(5,
     -   2)-EijR(8,3)-EijR(9,3)+EijR(17,3)))
       FI(254) = 4*(EijI(6,2)-EijI(8,2)+2*(EijI(2,2)+EijI(2,3)-EijI(5,
     -   2)-EijI(8,3)-EijI(9,3)+EijI(17,3)))
       F(254)=DCMPLX(FR(254),FI(254))
       FR(255) = 8*(EijR(2,3)-EijR(8,3)-EijR(11,3)+EijR(12,3)-2*(EijR(
     -   9,3)-EijR(17,3)))
       FI(255) = 8*(EijI(2,3)-EijI(8,3)-EijI(11,3)+EijI(12,3)-2*(EijI(
     -   9,3)-EijI(17,3)))
       F(255)=DCMPLX(FR(255),FI(255))
       FR(256) = 8*(EijR(2,2)+EijR(2,3)+EijR(3,2)-EijR(3,3)-2*EijR(8,2
     -   )-3*(EijR(9,3)-EijR(12,3)))
       FI(256) = 8*(EijI(2,2)+EijI(2,3)+EijI(3,2)-EijI(3,3)-2*EijI(8,2
     -   )-3*(EijI(9,3)-EijI(12,3)))
       F(256)=DCMPLX(FR(256),FI(256))
       FR(257) = 8*(EijR(2,2)+EijR(2,3)-EijR(8,2)-EijR(9,2)+EijR(10,2)
     -   -EijR(10,3)+EijR(12,3)-EijR(13,3)-2*(EijR(9,3)-EijR(20,3)))
       FI(257) = 8*(EijI(2,2)+EijI(2,3)-EijI(8,2)-EijI(9,2)+EijI(10,2)
     -   -EijI(10,3)+EijI(12,3)-EijI(13,3)-2*(EijI(9,3)-EijI(20,3)))
       F(257)=DCMPLX(FR(257),FI(257))
       FR(258) = 4*(EijR(3,2)-3*EijR(8,2)-4*EijR(9,3)+2*(EijR(2,2)+Eij
     -   R(2,3)+EijR(12,3)))
       FI(258) = 4*(EijI(3,2)-3*EijI(8,2)-4*EijI(9,3)+2*(EijI(2,2)+Eij
     -   I(2,3)+EijI(12,3)))
       F(258)=DCMPLX(FR(258),FI(258))
       FR(259) = 8*(EijR(2,3)+EijR(6,2)-EijR(7,2)-EijR(8,2)-EijR(8,3)+
     -   EijR(9,2)-EijR(9,3)-EijR(10,3)+EijR(17,3)+EijR(18,3)-EijR(19
     -   ,3)+EijR(20,3))
       FI(259) = 8*(EijI(2,3)+EijI(6,2)-EijI(7,2)-EijI(8,2)-EijI(8,3)+
     -   EijI(9,2)-EijI(9,3)-EijI(10,3)+EijI(17,3)+EijI(18,3)-EijI(19
     -   ,3)+EijI(20,3))
       F(259)=DCMPLX(FR(259),FI(259))
       FR(260) = 8*(EijR(2,2)+EijR(2,3)+EijR(3,2)-EijR(10,3)+EijR(12,3
     -   )-EijR(13,3)-2*(EijR(8,2)+EijR(9,3)-EijR(20,3)))
       FI(260) = 8*(EijI(2,2)+EijI(2,3)+EijI(3,2)-EijI(10,3)+EijI(12,3
     -   )-EijI(13,3)-2*(EijI(8,2)+EijI(9,3)-EijI(20,3)))
       F(260)=DCMPLX(FR(260),FI(260))
       FR(261) = 8*(EijR(2,2)+EijR(2,3)-EijR(8,2)-EijR(9,2)-EijR(9,3)+
     -   EijR(10,2)+EijR(15,3)-EijR(16,3)-2*(EijR(10,3)-EijR(20,3)))
       FI(261) = 8*(EijI(2,2)+EijI(2,3)-EijI(8,2)-EijI(9,2)-EijI(9,3)+
     -   EijI(10,2)+EijI(15,3)-EijI(16,3)-2*(EijI(10,3)-EijI(20,3)))
       F(261)=DCMPLX(FR(261),FI(261))
       FR(262) = 4*(-EijR(3,1)+EijR(4,1)-3*EijR(8,2)+EijR(10,2)+2*(Eij
     -   R(2,2)+EijR(2,3)-EijR(9,3)-EijR(10,3)+EijR(20,3)))
       FI(262) = 4*(-EijI(3,1)+EijI(4,1)-3*EijI(8,2)+EijI(10,2)+2*(Eij
     -   I(2,2)+EijI(2,3)-EijI(9,3)-EijI(10,3)+EijI(20,3)))
       F(262)=DCMPLX(FR(262),FI(262))
       FR(263) = 4*(EijR(2,2)-EijR(5,2)+2*(EijR(2,3)+EijR(6,2)-EijR(8,
     -   2)-EijR(8,3)-EijR(9,3)+EijR(17,3)))
       FI(263) = 4*(EijI(2,2)-EijI(5,2)+2*(EijI(2,3)+EijI(6,2)-EijI(8,
     -   2)-EijI(8,3)-EijI(9,3)+EijI(17,3)))
       F(263)=DCMPLX(FR(263),FI(263))
       FR(264) = 4*(EijR(2,1)+3*EijR(2,2)-EijR(3,1)-5*EijR(8,2)-4*EijR
     -   (9,3)+2*(EijR(2,3)+EijR(3,2)+EijR(12,3)))
       FI(264) = 4*(EijI(2,1)+3*EijI(2,2)-EijI(3,1)-5*EijI(8,2)-4*EijI
     -   (9,3)+2*(EijI(2,3)+EijI(3,2)+EijI(12,3)))
       F(264)=DCMPLX(FR(264),FI(264))
       FR(265) = 4*(EijR(2,1)-EijR(4,1)+3*(EijR(2,2)-EijR(9,2))+2*(Eij
     -   R(2,3)-EijR(8,2)-EijR(9,3)+EijR(10,2)-EijR(10,3)+EijR(20,3))
     -   )
       FI(265) = 4*(EijI(2,1)-EijI(4,1)+3*(EijI(2,2)-EijI(9,2))+2*(Eij
     -   I(2,3)-EijI(8,2)-EijI(9,3)+EijI(10,2)-EijI(10,3)+EijI(20,3))
     -   )
       F(265)=DCMPLX(FR(265),FI(265))
       FR(266) = 4*(EijR(2,1)-EijR(3,1)+3*(EijR(2,2)-EijR(8,2))+2*(Eij
     -   R(2,3)-EijR(9,3)))
       FI(266) = 4*(EijI(2,1)-EijI(3,1)+3*(EijI(2,2)-EijI(8,2))+2*(Eij
     -   I(2,3)-EijI(9,3)))
       F(266)=DCMPLX(FR(266),FI(266))
       P(640) = p2sq-s15-s34+3*P(137)
       P(641) = 5*p2sq+3*s12
       P(642) = 5*p3sq+3*P(16)
       P(643) = -s23+s45+4*P(33)
       P(644) = 6*s12-2*s23-3*s45+5*P(299)
       P(645) = s15-s23+3*P(7)+4*P(14)
       P(646) = 3*P(7)+4*P(14)
       FR(267) = 5*D02345R-3*(Dij2345R(1,1)-Dij2345R(2,1))-26*EijR(11,
     -   2)+8*(EijR(22,3)-EijR(23,3))+2*((EE0R+EijR(1,1))*P(3)+Is15s3
     -   4*(B025R+2*Cij125R(4,2)+(C0125R+Cij125R(1,1))*P(3)))-(EijR(2
     -   ,1)-4*EijR(5,2))*P(33)+EijR(4,1)*P(533)+EijR(3,1)*P(595)+Is1
     -   2*(5*C0345R+Cij135R(2,1)-14*Dij1345R(7,2)+Dij1345R(6,2)*P(14
     -   )-2*(p3sq*Dij1345R(2,2)-D01345R*P(35))+Is45*(2*B012R+4*Cij14
     -   5R(4,2)+Cij145R(1,1)*P(95)+C0145R*P(185))+Dij1345R(4,2)*P(19
     -   8)+Dij1345R(3,1)*P(328)+Dij1345R(1,1)*P(638)-Dij1345R(2,1)*P
     -   (639))+Is34*(5*C0235R+Cij125R(2,1)-14*Dij1235R(7,2)-2*(p2sq*
     -   Dij1235R(2,2)-(D01235R+Dij1235R(1,1))*P(3))+Dij1235R(4,2)*P(
     -   33)+Dij1235R(6,2)*P(436)+Is12*(2*B035R+4*Cij135R(4,2)-C0135R
     -   *P(513)+Cij135R(1,1)*P(514))-Dij1235R(2,1)*P(611)+Dij1235R(3
     -   ,1)*P(640))-EijR(2,2)*P(641)-EijR(3,2)*P(642)-EijR(6,2)*P(64
     -   3)+EijR(8,2)*P(644)-EijR(9,2)*P(645)+EijR(10,2)*P(646)
       FI(267) = 5*D02345I-3*(Dij2345I(1,1)-Dij2345I(2,1))-26*EijI(11,
     -   2)+8*(EijI(22,3)-EijI(23,3))+2*((EE0I+EijI(1,1))*P(3)+Is15s3
     -   4*(B025I+2*Cij125I(4,2)+(C0125I+Cij125I(1,1))*P(3)))-(EijI(2
     -   ,1)-4*EijI(5,2))*P(33)+EijI(4,1)*P(533)+EijI(3,1)*P(595)+Is1
     -   2*(5*C0345I+Cij135I(2,1)-14*Dij1345I(7,2)+Dij1345I(6,2)*P(14
     -   )-2*(p3sq*Dij1345I(2,2)-D01345I*P(35))+Is45*(2*B012I+4*Cij14
     -   5I(4,2)+Cij145I(1,1)*P(95)+C0145I*P(185))+Dij1345I(4,2)*P(19
     -   8)+Dij1345I(3,1)*P(328)+Dij1345I(1,1)*P(638)-Dij1345I(2,1)*P
     -   (639))+Is34*(5*C0235I+Cij125I(2,1)-14*Dij1235I(7,2)-2*(p2sq*
     -   Dij1235I(2,2)-(D01235I+Dij1235I(1,1))*P(3))+Dij1235I(4,2)*P(
     -   33)+Dij1235I(6,2)*P(436)+Is12*(2*B035I+4*Cij135I(4,2)-C0135I
     -   *P(513)+Cij135I(1,1)*P(514))-Dij1235I(2,1)*P(611)+Dij1235I(3
     -   ,1)*P(640))-EijI(2,2)*P(641)-EijI(3,2)*P(642)-EijI(6,2)*P(64
     -   3)+EijI(8,2)*P(644)-EijI(9,2)*P(645)+EijI(10,2)*P(646)
       F(267)=DCMPLX(FR(267),FI(267))
       P(647) = 3*s12-2*s45
       P(648) = 3*p2sq+s12-s23-s45
       FR(268) = -3*D02345R+Dij2345R(1,1)-EijR(2,1)*P(111)+EijR(2,2)*P
     -   (166)-EijR(9,2)*P(406)-2*(Dij1235R(3,2)-Dij1235R(6,2)+Dij234
     -   5R(1,2)+Dij2345R(2,1)-Dij2345R(4,2)-8*EijR(11,2)+6*EijR(22,3
     -   )-4*EijR(23,3)+(EE0R+EijR(1,1))*P(3)+EijR(2,3)*P(5)-EijR(3,2
     -   )*P(16)-(EijR(8,3)-EijR(17,3))*P(33)-(-EijR(10,3)+EijR(20,3)
     -   )*P(101)+EijR(12,3)*P(198)-(-EijR(5,2)+EijR(6,2))*P(245)-Eij
     -   R(10,2)*P(491)-EijR(9,3)*P(505))-EijR(3,1)*P(532)-EijR(4,1)*
     -   P(533)+Is12*(C0135R-3*C0345R+Cij135R(1,1)-2*(Cij135R(2,1)-Ci
     -   j135R(2,2)+Cij135R(3,2))+Cij345R(1,1)+10*Dij1345R(7,2)+4*(Di
     -   j1345R(11,3)-Dij1345R(12,3))+Dij1345R(6,2)*P(101)-Dij1345R(5
     -   ,2)*P(292)+Is45*(-2*(B012R+2*Cij145R(4,2))-Cij145R(1,1)*P(95
     -   )+C0145R*P(387))-D01345R*P(525)+Dij1345R(1,1)*P(527)-Dij1345
     -   R(3,1)*P(529)+(-Dij1345R(1,2)+Dij1345R(4,2))*P(530)+Dij1345R
     -   (2,2)*P(531)+Dij1345R(2,1)*P(647))-EijR(8,2)*P(648)
       FI(268) = -3*D02345I+Dij2345I(1,1)-EijI(2,1)*P(111)+EijI(2,2)*P
     -   (166)-EijI(9,2)*P(406)-2*(Dij1235I(3,2)-Dij1235I(6,2)+Dij234
     -   5I(1,2)+Dij2345I(2,1)-Dij2345I(4,2)-8*EijI(11,2)+6*EijI(22,3
     -   )-4*EijI(23,3)+(EE0I+EijI(1,1))*P(3)+EijI(2,3)*P(5)-EijI(3,2
     -   )*P(16)-(EijI(8,3)-EijI(17,3))*P(33)-(-EijI(10,3)+EijI(20,3)
     -   )*P(101)+EijI(12,3)*P(198)-(-EijI(5,2)+EijI(6,2))*P(245)-Eij
     -   I(10,2)*P(491)-EijI(9,3)*P(505))-EijI(3,1)*P(532)-EijI(4,1)*
     -   P(533)+Is12*(C0135I-3*C0345I+Cij135I(1,1)-2*(Cij135I(2,1)-Ci
     -   j135I(2,2)+Cij135I(3,2))+Cij345I(1,1)+10*Dij1345I(7,2)+4*(Di
     -   j1345I(11,3)-Dij1345I(12,3))+Dij1345I(6,2)*P(101)-Dij1345I(5
     -   ,2)*P(292)+Is45*(-2*(B012I+2*Cij145I(4,2))-Cij145I(1,1)*P(95
     -   )+C0145I*P(387))-D01345I*P(525)+Dij1345I(1,1)*P(527)-Dij1345
     -   I(3,1)*P(529)+(-Dij1345I(1,2)+Dij1345I(4,2))*P(530)+Dij1345I
     -   (2,2)*P(531)+Dij1345I(2,1)*P(647))-EijI(8,2)*P(648)
       F(268)=DCMPLX(FR(268),FI(268))
       P(649) = 2*s12+s45
       P(650) = p3sq-5*s12-2*s45
       P(651) = p2sq+p3sq-p5sq-s34+s45
       FR(269) = Is12*(C0135R+Cij135R(1,1)-3*Cij145R(2,1)+2*(C0145R-C0
     -   345R-Cij135R(2,1)+Cij135R(2,2)-Cij135R(3,2)+Cij145R(1,2)-Cij
     -   145R(3,2))+Cij345R(2,1)+6*Dij1345R(7,2)+4*(Cij145R(1,1)+Dij1
     -   345R(11,3)-Dij1345R(13,3))+Dij1345R(3,2)*P(101)+Dij1345R(4,2
     -   )*P(198)-Dij1345R(1,2)*P(530)+Dij1345R(6,2)*P(531)+Dij1345R(
     -   2,1)*P(561)-Dij1345R(3,1)*P(562)-Dij1345R(5,2)*P(563)-D01345
     -   R*P(649)+Dij1345R(1,1)*P(650))+2*(-Dij1235R(3,2)+Dij1235R(6,
     -   2)+Dij2345R(1,1)-Dij2345R(1,2)-Dij2345R(3,1)+Dij2345R(5,2)+p
     -   2sq*EijR(2,2)-p3sq*EijR(3,1)-6*EijR(22,3)+4*(EijR(11,2)+EijR
     -   (24,3))-EijR(2,3)*P(5)+(-EijR(8,2)+EijR(10,2))*P(16)+(EijR(8
     -   ,3)-EijR(18,3))*P(33)+EijR(15,3)*P(101)+(EijR(9,3)-EijR(20,3
     -   ))*P(198)+(-EijR(5,2)+EijR(7,2))*P(245)+EijR(4,2)*P(491)-Eij
     -   R(2,1)*P(564)+EijR(4,1)*P(565)+EijR(10,3)*P(593)-EijR(9,2)*P
     -   (651))
       FI(269) = Is12*(C0135I+Cij135I(1,1)-3*Cij145I(2,1)+2*(C0145I-C0
     -   345I-Cij135I(2,1)+Cij135I(2,2)-Cij135I(3,2)+Cij145I(1,2)-Cij
     -   145I(3,2))+Cij345I(2,1)+6*Dij1345I(7,2)+4*(Cij145I(1,1)+Dij1
     -   345I(11,3)-Dij1345I(13,3))+Dij1345I(3,2)*P(101)+Dij1345I(4,2
     -   )*P(198)-Dij1345I(1,2)*P(530)+Dij1345I(6,2)*P(531)+Dij1345I(
     -   2,1)*P(561)-Dij1345I(3,1)*P(562)-Dij1345I(5,2)*P(563)-D01345
     -   I*P(649)+Dij1345I(1,1)*P(650))+2*(-Dij1235I(3,2)+Dij1235I(6,
     -   2)+Dij2345I(1,1)-Dij2345I(1,2)-Dij2345I(3,1)+Dij2345I(5,2)+p
     -   2sq*EijI(2,2)-p3sq*EijI(3,1)-6*EijI(22,3)+4*(EijI(11,2)+EijI
     -   (24,3))-EijI(2,3)*P(5)+(-EijI(8,2)+EijI(10,2))*P(16)+(EijI(8
     -   ,3)-EijI(18,3))*P(33)+EijI(15,3)*P(101)+(EijI(9,3)-EijI(20,3
     -   ))*P(198)+(-EijI(5,2)+EijI(7,2))*P(245)+EijI(4,2)*P(491)-Eij
     -   I(2,1)*P(564)+EijI(4,1)*P(565)+EijI(10,3)*P(593)-EijI(9,2)*P
     -   (651))
       F(269)=DCMPLX(FR(269),FI(269))
       P(652) = p2sq+2*s12
       P(653) = p3sq+3*s12+s23+2*P(414)
       P(654) = 5*p3sq+s12-s45
       P(655) = p5sq-s45+6*P(14)
       FR(270) = -D01235R+EE0R*s45-Dij1235R(1,1)+EijR(4,1)*P(7)-EijR(1
     -   ,2)*P(50)-EijR(1,1)*P(77)-(EijR(3,1)-2*(EijR(9,3)-EijR(17,3)
     -   ))*P(198)+EijR(2,1)*P(547)+2*(Dij1235R(3,1)-Dij1235R(5,2)+Di
     -   j1235R(6,2)-Dij2345R(1,2)-2*(Dij2345R(1,1)-p2sq*EijR(8,3))-3
     -   *EijR(11,2)+4*(EijR(21,3)-EijR(22,3))-EijR(2,3)*P(5)-EijR(5,
     -   3)*P(33)+EijR(18,3)*P(101)-EijR(9,2)*P(153)+EijR(8,2)*P(309)
     -   +EijR(10,3)*P(312)-EijR(2,2)*P(652))+EijR(5,2)*P(653)-EijR(6
     -   ,2)*P(654)+EijR(7,2)*P(655)
       FI(270) = -D01235I+EE0I*s45-Dij1235I(1,1)+EijI(4,1)*P(7)-EijI(1
     -   ,2)*P(50)-EijI(1,1)*P(77)-(EijI(3,1)-2*(EijI(9,3)-EijI(17,3)
     -   ))*P(198)+EijI(2,1)*P(547)+2*(Dij1235I(3,1)-Dij1235I(5,2)+Di
     -   j1235I(6,2)-Dij2345I(1,2)-2*(Dij2345I(1,1)-p2sq*EijI(8,3))-3
     -   *EijI(11,2)+4*(EijI(21,3)-EijI(22,3))-EijI(2,3)*P(5)-EijI(5,
     -   3)*P(33)+EijI(18,3)*P(101)-EijI(9,2)*P(153)+EijI(8,2)*P(309)
     -   +EijI(10,3)*P(312)-EijI(2,2)*P(652))+EijI(5,2)*P(653)-EijI(6
     -   ,2)*P(654)+EijI(7,2)*P(655)
       F(270)=DCMPLX(FR(270),FI(270))
       P(656) = -3*s45+2*P(161)
       P(657) = p2sq-s23+s45
       P(658) = s12-2*P(657)
       P(659) = p5sq-s45+5*P(14)
       FR(271) = Dij1235R(3,1)-Dij2345R(1,1)+p3sq*EijR(3,1)-4*EijR(11,
     -   2)-12*EijR(22,3)-EijR(5,2)*P(111)+2*(Dij1235R(6,2)-Dij2345R(
     -   1,2)-EijR(2,3)*P(5)-EijR(4,1)*P(14)+EijR(8,3)*P(33)-EijR(10,
     -   3)*P(101)+EijR(9,3)*P(198))+(EE0R+EijR(1,1))*P(251)-EijR(2,2
     -   )*P(517)+EijR(8,2)*P(549)+Is12*(C0135R+C0145R-C0345R+Cij135R
     -   (1,1)-3*(Cij135R(2,1)-Cij145R(1,1))-2*(Cij135R(3,2)-Cij145R(
     -   1,2))+4*(Dij1345R(7,2)+Dij1345R(11,3))-Dij1345R(3,1)*P(41)+D
     -   ij1345R(4,2)*P(198)-Dij1345R(5,2)*P(292)+Dij1345R(2,1)*P(305
     -   )-Dij1345R(1,2)*P(530)+D01345R*P(544)+Dij1345R(1,1)*P(656))+
     -   EijR(2,1)*P(658)-EijR(9,2)*P(659)
       FI(271) = Dij1235I(3,1)-Dij2345I(1,1)+p3sq*EijI(3,1)-4*EijI(11,
     -   2)-12*EijI(22,3)-EijI(5,2)*P(111)+2*(Dij1235I(6,2)-Dij2345I(
     -   1,2)-EijI(2,3)*P(5)-EijI(4,1)*P(14)+EijI(8,3)*P(33)-EijI(10,
     -   3)*P(101)+EijI(9,3)*P(198))+(EE0I+EijI(1,1))*P(251)-EijI(2,2
     -   )*P(517)+EijI(8,2)*P(549)+Is12*(C0135I+C0145I-C0345I+Cij135I
     -   (1,1)-3*(Cij135I(2,1)-Cij145I(1,1))-2*(Cij135I(3,2)-Cij145I(
     -   1,2))+4*(Dij1345I(7,2)+Dij1345I(11,3))-Dij1345I(3,1)*P(41)+D
     -   ij1345I(4,2)*P(198)-Dij1345I(5,2)*P(292)+Dij1345I(2,1)*P(305
     -   )-Dij1345I(1,2)*P(530)+D01345I*P(544)+Dij1345I(1,1)*P(656))+
     -   EijI(2,1)*P(658)-EijI(9,2)*P(659)
       F(271)=DCMPLX(FR(271),FI(271))
       P(660) = p2sq-5*s12
       P(661) = p2sq-5*s12+2*P(3)
       P(662) = p2sq-4*p3sq-s23+2*s34-3*P(16)
       P(663) = -s15+s23+3*P(7)+4*P(14)
       FR(272) = 5*D02345R+6*Dij2345R(1,1)+Dij2345R(1,2)-2*(EijR(11,2)
     -   -6*EijR(21,3)+7*EijR(22,3)+(EijR(1,1)+EijR(1,2)-EijR(2,1))*P
     -   (3))+EijR(9,3)*P(10)+(EijR(10,3)-EijR(18,3))*P(21)+EijR(8,3)
     -   *P(33)+EijR(17,3)*P(251)+EijR(2,3)*P(601)-EijR(2,2)*P(660)+E
     -   ijR(5,2)*P(661)+(-EijR(6,2)+EijR(8,2))*P(662)+(-EijR(7,2)+Ei
     -   jR(9,2))*P(663)
       FI(272) = 5*D02345I+6*Dij2345I(1,1)+Dij2345I(1,2)-2*(EijI(11,2)
     -   -6*EijI(21,3)+7*EijI(22,3)+(EijI(1,1)+EijI(1,2)-EijI(2,1))*P
     -   (3))+EijI(9,3)*P(10)+(EijI(10,3)-EijI(18,3))*P(21)+EijI(8,3)
     -   *P(33)+EijI(17,3)*P(251)+EijI(2,3)*P(601)-EijI(2,2)*P(660)+E
     -   ijI(5,2)*P(661)+(-EijI(6,2)+EijI(8,2))*P(662)+(-EijI(7,2)+Ei
     -   jI(9,2))*P(663)
       F(272)=DCMPLX(FR(272),FI(272))
       P(664) = 3*s12+4*s34-2*P(537)
       P(665) = p2sq-s12-s23+4*s45-2*P(3)
       P(666) = p2sq+p3sq-s23+2*s34
       P(667) = -s23+s45+2*P(3)+4*P(33)
       P(668) = 7*p2sq+p3sq-3*s23-2*P(6)
       FR(273) = 5*D02345R+Dij2345R(1,2)-Dij2345R(4,2)-26*EijR(11,2)-1
     -   2*(EijR(22,3)-EijR(23,3))+EijR(10,3)*P(21)+(-2*EijR(9,2)+Eij
     -   R(10,2))*P(66)+EijR(20,3)*P(115)+EijR(9,3)*P(245)+EijR(12,3)
     -   *P(251)+2*(Dij2345R(1,1)-Dij2345R(2,1)+(EE0R+EijR(1,1))*P(3)
     -   +Is15s34*(B025R+2*Cij125R(4,2)+(C0125R+Cij125R(1,1))*P(3))-E
     -   ijR(2,2)*P(166)+EijR(5,2)*P(332))-EijR(2,1)*P(375)+EijR(4,1)
     -   *P(533)+EijR(2,3)*P(601)+Is34*(5*C0235R+Cij125R(2,1)-14*Dij1
     -   235R(7,2)-2*(p2sq*Dij1235R(2,2)-(D01235R+Dij1235R(1,1))*P(3)
     -   )+Dij1235R(4,2)*P(33)+Dij1235R(6,2)*P(436)+Is12*(2*B035R+4*C
     -   ij135R(4,2)-C0135R*P(513)+Cij135R(1,1)*P(514))-Dij1235R(2,1)
     -   *P(611)+Dij1235R(3,1)*P(640))+Is12*(Cij135R(2,1)-Cij345R(1,1
     -   )-8*Dij1345R(7,2)+4*(C0345R+Dij1345R(11,3)-Dij1345R(12,3))+D
     -   ij1345R(6,2)*P(7)+2*D01345R*P(35)+Is45*(2*B012R+4*Cij145R(4,
     -   2)+Cij145R(1,1)*P(95)+C0145R*P(185))-Dij1345R(2,2)*P(199)+Di
     -   j1345R(3,1)*P(328)+Dij1345R(5,2)*P(491)-Dij1345R(1,2)*P(513)
     -   +Dij1345R(1,1)*P(534)-Dij1345R(2,1)*P(536)+Dij1345R(4,2)*P(6
     -   64))+EijR(3,1)*P(665)-EijR(3,2)*P(666)-EijR(6,2)*P(667)+EijR
     -   (8,2)*P(668)
       FI(273) = 5*D02345I+Dij2345I(1,2)-Dij2345I(4,2)-26*EijI(11,2)-1
     -   2*(EijI(22,3)-EijI(23,3))+EijI(10,3)*P(21)+(-2*EijI(9,2)+Eij
     -   I(10,2))*P(66)+EijI(20,3)*P(115)+EijI(9,3)*P(245)+EijI(12,3)
     -   *P(251)+2*(Dij2345I(1,1)-Dij2345I(2,1)+(EE0I+EijI(1,1))*P(3)
     -   +Is15s34*(B025I+2*Cij125I(4,2)+(C0125I+Cij125I(1,1))*P(3))-E
     -   ijI(2,2)*P(166)+EijI(5,2)*P(332))-EijI(2,1)*P(375)+EijI(4,1)
     -   *P(533)+EijI(2,3)*P(601)+Is34*(5*C0235I+Cij125I(2,1)-14*Dij1
     -   235I(7,2)-2*(p2sq*Dij1235I(2,2)-(D01235I+Dij1235I(1,1))*P(3)
     -   )+Dij1235I(4,2)*P(33)+Dij1235I(6,2)*P(436)+Is12*(2*B035I+4*C
     -   ij135I(4,2)-C0135I*P(513)+Cij135I(1,1)*P(514))-Dij1235I(2,1)
     -   *P(611)+Dij1235I(3,1)*P(640))+Is12*(Cij135I(2,1)-Cij345I(1,1
     -   )-8*Dij1345I(7,2)+4*(C0345I+Dij1345I(11,3)-Dij1345I(12,3))+D
     -   ij1345I(6,2)*P(7)+2*D01345I*P(35)+Is45*(2*B012I+4*Cij145I(4,
     -   2)+Cij145I(1,1)*P(95)+C0145I*P(185))-Dij1345I(2,2)*P(199)+Di
     -   j1345I(3,1)*P(328)+Dij1345I(5,2)*P(491)-Dij1345I(1,2)*P(513)
     -   +Dij1345I(1,1)*P(534)-Dij1345I(2,1)*P(536)+Dij1345I(4,2)*P(6
     -   64))+EijI(3,1)*P(665)-EijI(3,2)*P(666)-EijI(6,2)*P(667)+EijI
     -   (8,2)*P(668)
       F(273)=DCMPLX(FR(273),FI(273))
       P(669) = -p3sq+s12+2*s34-s45
       P(670) = s12-P(30)
       P(671) = p5sq+s15+s23+2*P(670)
       P(672) = 5*p2sq+p3sq+s23-2*P(225)
       FR(274) = 5*D02345R+Dij2345R(1,2)-Dij2345R(5,2)-20*EijR(11,2)-1
     -   2*(EijR(22,3)-EijR(24,3))+(EijR(9,3)-EijR(20,3))*P(10)+EijR(
     -   4,2)*P(66)+EijR(15,3)*P(115)-EijR(2,1)*P(375)+2*(Dij2345R(1,
     -   1)-Dij2345R(3,1)+(EE0R+EijR(1,1))*P(3)+Is15s34*(B025R+2*Cij1
     -   25R(4,2)+(C0125R+Cij125R(1,1))*P(3))-EijR(2,2)*P(166)+EijR(5
     -   ,2)*P(332)+EijR(8,2)*P(432))+EijR(10,3)*P(488)+EijR(3,1)*P(5
     -   54)+EijR(2,3)*P(601)+Is34*(5*C0235R+Cij125R(2,1)-14*Dij1235R
     -   (7,2)-2*(p2sq*Dij1235R(2,2)-(D01235R+Dij1235R(1,1))*P(3))+Di
     -   j1235R(4,2)*P(33)+Dij1235R(6,2)*P(436)+Is12*(2*B035R+4*Cij13
     -   5R(4,2)-C0135R*P(513)+Cij135R(1,1)*P(514))-Dij1235R(2,1)*P(6
     -   11)+Dij1235R(3,1)*P(640))-EijR(10,2)*P(666)-EijR(7,2)*P(667)
     -   +Is12*(Cij135R(2,1)-Cij345R(2,1)+4*(C0345R+Dij1345R(11,3)-Di
     -   j1345R(13,3))+Dij1345R(3,2)*P(7)-2*(Dij1345R(7,2)-D01345R*P(
     -   35))+Is45*(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)*P(95)+C0145R*
     -   P(185))-Dij1345R(6,2)*P(199)-Dij1345R(2,1)*P(421)-Dij1345R(1
     -   ,2)*P(513)+Dij1345R(1,1)*P(534)+Dij1345R(3,1)*P(551)+Dij1345
     -   R(5,2)*P(553)+Dij1345R(4,2)*P(669))+EijR(4,1)*P(671)+EijR(9,
     -   2)*P(672)
       FI(274) = 5*D02345I+Dij2345I(1,2)-Dij2345I(5,2)-20*EijI(11,2)-1
     -   2*(EijI(22,3)-EijI(24,3))+(EijI(9,3)-EijI(20,3))*P(10)+EijI(
     -   4,2)*P(66)+EijI(15,3)*P(115)-EijI(2,1)*P(375)+2*(Dij2345I(1,
     -   1)-Dij2345I(3,1)+(EE0I+EijI(1,1))*P(3)+Is15s34*(B025I+2*Cij1
     -   25I(4,2)+(C0125I+Cij125I(1,1))*P(3))-EijI(2,2)*P(166)+EijI(5
     -   ,2)*P(332)+EijI(8,2)*P(432))+EijI(10,3)*P(488)+EijI(3,1)*P(5
     -   54)+EijI(2,3)*P(601)+Is34*(5*C0235I+Cij125I(2,1)-14*Dij1235I
     -   (7,2)-2*(p2sq*Dij1235I(2,2)-(D01235I+Dij1235I(1,1))*P(3))+Di
     -   j1235I(4,2)*P(33)+Dij1235I(6,2)*P(436)+Is12*(2*B035I+4*Cij13
     -   5I(4,2)-C0135I*P(513)+Cij135I(1,1)*P(514))-Dij1235I(2,1)*P(6
     -   11)+Dij1235I(3,1)*P(640))-EijI(10,2)*P(666)-EijI(7,2)*P(667)
     -   +Is12*(Cij135I(2,1)-Cij345I(2,1)+4*(C0345I+Dij1345I(11,3)-Di
     -   j1345I(13,3))+Dij1345I(3,2)*P(7)-2*(Dij1345I(7,2)-D01345I*P(
     -   35))+Is45*(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)*P(95)+C0145I*
     -   P(185))-Dij1345I(6,2)*P(199)-Dij1345I(2,1)*P(421)-Dij1345I(1
     -   ,2)*P(513)+Dij1345I(1,1)*P(534)+Dij1345I(3,1)*P(551)+Dij1345
     -   I(5,2)*P(553)+Dij1345I(4,2)*P(669))+EijI(4,1)*P(671)+EijI(9,
     -   2)*P(672)
       F(274)=DCMPLX(FR(274),FI(274))
       P(673) = p2sq-p3sq-s12-s23+s34+s45
       P(674) = p3sq+p5sq-s15+s23-s34-s45
       FR(275) = 5*D02345R+Dij2345R(1,2)-12*(EijR(11,2)+EijR(22,3))+Ei
     -   jR(9,3)*P(10)+EijR(10,3)*P(21)+4*(Dij2345R(1,1)-EijR(2,2)*P(
     -   33))+EijR(4,1)*P(339)+EijR(3,1)*P(524)-EijR(2,1)*P(568)+EijR
     -   (2,3)*P(601)+Is12*(2*(C0345R+Dij1345R(7,2)+2*Dij1345R(11,3)+
     -   D01345R*P(39))-Dij1345R(2,1)*P(79)+Dij1345R(3,1)*P(101)+Is45
     -   *(2*B012R+4*Cij145R(4,2)+Cij145R(1,1)*P(95)+C0145R*P(185))+D
     -   ij1345R(5,2)*P(491)-Dij1345R(1,2)*P(513)+Dij1345R(1,1)*P(567
     -   )+Dij1345R(4,2)*P(669))+2*((EE0R+EijR(1,1))*P(3)+EijR(5,2)*P
     -   (330)+EijR(8,2)*P(673)+EijR(9,2)*P(674))
       FI(275) = 5*D02345I+Dij2345I(1,2)-12*(EijI(11,2)+EijI(22,3))+Ei
     -   jI(9,3)*P(10)+EijI(10,3)*P(21)+4*(Dij2345I(1,1)-EijI(2,2)*P(
     -   33))+EijI(4,1)*P(339)+EijI(3,1)*P(524)-EijI(2,1)*P(568)+EijI
     -   (2,3)*P(601)+Is12*(2*(C0345I+Dij1345I(7,2)+2*Dij1345I(11,3)+
     -   D01345I*P(39))-Dij1345I(2,1)*P(79)+Dij1345I(3,1)*P(101)+Is45
     -   *(2*B012I+4*Cij145I(4,2)+Cij145I(1,1)*P(95)+C0145I*P(185))+D
     -   ij1345I(5,2)*P(491)-Dij1345I(1,2)*P(513)+Dij1345I(1,1)*P(567
     -   )+Dij1345I(4,2)*P(669))+2*((EE0I+EijI(1,1))*P(3)+EijI(5,2)*P
     -   (330)+EijI(8,2)*P(673)+EijI(9,2)*P(674))
       F(275)=DCMPLX(FR(275),FI(275))
       P(675) = 4*p2sq+p3sq+s12-s23-2*s34+s45
       P(676) = p3sq-P(30)
       P(677) = -s23+4*P(5)+2*P(676)
       FR(276) = -4*Dij1235R(3,2)+7*Dij1235R(6,2)-Dij2345R(1,2)-3*(Dij
     -   1235R(2,1)+Dij1235R(2,2)-Dij1235R(3,1)-Dij2345R(1,1)+Dij2345
     -   R(2,1))+Dij2345R(4,2)+12*EijR(11,2)-10*(EijR(22,3)-EijR(23,3
     -   ))+(3*(EijR(9,2)-EijR(10,2))-EijR(10,3)+EijR(20,3))*P(7)+(Ei
     -   jR(8,3)-EijR(17,3))*P(50)+(EijR(2,1)-EijR(3,1))*P(481)+EijR(
     -   12,3)*P(531)+EijR(3,2)*P(561)-EijR(2,3)*P(618)-EijR(9,3)*P(6
     -   23)+(-EijR(5,2)+EijR(6,2))*P(643)+EijR(2,2)*P(675)-EijR(8,2)
     -   *P(677)
       FI(276) = -4*Dij1235I(3,2)+7*Dij1235I(6,2)-Dij2345I(1,2)-3*(Dij
     -   1235I(2,1)+Dij1235I(2,2)-Dij1235I(3,1)-Dij2345I(1,1)+Dij2345
     -   I(2,1))+Dij2345I(4,2)+12*EijI(11,2)-10*(EijI(22,3)-EijI(23,3
     -   ))+(3*(EijI(9,2)-EijI(10,2))-EijI(10,3)+EijI(20,3))*P(7)+(Ei
     -   jI(8,3)-EijI(17,3))*P(50)+(EijI(2,1)-EijI(3,1))*P(481)+EijI(
     -   12,3)*P(531)+EijI(3,2)*P(561)-EijI(2,3)*P(618)-EijI(9,3)*P(6
     -   23)+(-EijI(5,2)+EijI(6,2))*P(643)+EijI(2,2)*P(675)-EijI(8,2)
     -   *P(677)
       F(276)=DCMPLX(FR(276),FI(276))
       P(678) = s45-2*P(29)
       P(679) = -s45+2*P(29)
       FR(277) = -3*(Dij1235R(2,2)-Dij1235R(4,2))-Dij2345R(1,1)-Dij234
     -   5R(1,2)+10*(EijR(21,3)-EijR(22,3))-EijR(10,3)*P(7)+EijR(18,3
     -   )*P(7)-4*(Dij1235R(5,2)-Dij1235R(6,2)+p3sq*EijR(6,2)-p3sq*Ei
     -   jR(8,2)-EijR(7,2)*P(14)+EijR(9,2)*P(14))-EijR(5,3)*P(50)-Eij
     -   R(9,3)*P(531)+EijR(17,3)*P(531)-EijR(2,3)*P(618)-EijR(8,3)*P
     -   (622)+EijR(2,2)*P(678)+EijR(5,2)*P(679)
       FI(277) = -3*(Dij1235I(2,2)-Dij1235I(4,2))-Dij2345I(1,1)-Dij234
     -   5I(1,2)+10*(EijI(21,3)-EijI(22,3))-EijI(10,3)*P(7)+EijI(18,3
     -   )*P(7)-4*(Dij1235I(5,2)-Dij1235I(6,2)+p3sq*EijI(6,2)-p3sq*Ei
     -   jI(8,2)-EijI(7,2)*P(14)+EijI(9,2)*P(14))-EijI(5,3)*P(50)-Eij
     -   I(9,3)*P(531)+EijI(17,3)*P(531)-EijI(2,3)*P(618)-EijI(8,3)*P
     -   (622)+EijI(2,2)*P(678)+EijI(5,2)*P(679)
       F(277)=DCMPLX(FR(277),FI(277))
       P(680) = p2sq+s45
       P(681) = p3sq-3*p5sq+s12-s23-2*s34+4*P(680)
       FR(278) = -4*Dij1235R(3,2)+7*Dij1235R(6,2)-Dij2345R(1,2)-3*(Dij
     -   1235R(2,1)+Dij1235R(2,2)-Dij1235R(3,1)-Dij2345R(1,1)+Dij2345
     -   R(3,1))+Dij2345R(5,2)+6*EijR(11,2)-12*EijR(22,3)+10*EijR(24,
     -   3)+(-3*EijR(4,2)+EijR(15,3))*P(7)+(EijR(8,3)-EijR(18,3))*P(5
     -   0)+(EijR(2,1)-EijR(4,1))*P(481)+(-EijR(9,3)+EijR(20,3))*P(53
     -   1)+(-EijR(8,2)+EijR(10,2))*P(561)-EijR(2,3)*P(618)-EijR(10,3
     -   )*P(627)+(-EijR(5,2)+EijR(7,2))*P(643)+EijR(2,2)*P(675)-EijR
     -   (9,2)*P(681)
       FI(278) = -4*Dij1235I(3,2)+7*Dij1235I(6,2)-Dij2345I(1,2)-3*(Dij
     -   1235I(2,1)+Dij1235I(2,2)-Dij1235I(3,1)-Dij2345I(1,1)+Dij2345
     -   I(3,1))+Dij2345I(5,2)+6*EijI(11,2)-12*EijI(22,3)+10*EijI(24,
     -   3)+(-3*EijI(4,2)+EijI(15,3))*P(7)+(EijI(8,3)-EijI(18,3))*P(5
     -   0)+(EijI(2,1)-EijI(4,1))*P(481)+(-EijI(9,3)+EijI(20,3))*P(53
     -   1)+(-EijI(8,2)+EijI(10,2))*P(561)-EijI(2,3)*P(618)-EijI(10,3
     -   )*P(627)+(-EijI(5,2)+EijI(7,2))*P(643)+EijI(2,2)*P(675)-EijI
     -   (9,2)*P(681)
       F(278)=DCMPLX(FR(278),FI(278))
       P(682) = -2*s34+s45
       P(683) = p3sq-s12-s23+s45+2*P(243)
       P(684) = -p5sq+s45+2*P(14)
       FR(279) = -3*(Dij1235R(2,1)+Dij1235R(2,2))+4*Dij1235R(6,2)+Dij2
     -   345R(1,1)-Dij2345R(1,2)-12*EijR(22,3)-EijR(10,3)*P(7)+2*(Dij
     -   1235R(3,1)+p3sq*EijR(3,1)-EijR(11,2)-EijR(4,1)*P(14))+EijR(8
     -   ,3)*P(50)-EijR(5,2)*P(245)+(EijR(8,2)-EijR(9,3))*P(531)-EijR
     -   (2,3)*P(618)+EijR(2,1)*P(682)+EijR(2,2)*P(683)-EijR(9,2)*P(6
     -   84)
       FI(279) = -3*(Dij1235I(2,1)+Dij1235I(2,2))+4*Dij1235I(6,2)+Dij2
     -   345I(1,1)-Dij2345I(1,2)-12*EijI(22,3)-EijI(10,3)*P(7)+2*(Dij
     -   1235I(3,1)+p3sq*EijI(3,1)-EijI(11,2)-EijI(4,1)*P(14))+EijI(8
     -   ,3)*P(50)-EijI(5,2)*P(245)+(EijI(8,2)-EijI(9,3))*P(531)-EijI
     -   (2,3)*P(618)+EijI(2,1)*P(682)+EijI(2,2)*P(683)-EijI(9,2)*P(6
     -   84)
       F(279)=DCMPLX(FR(279),FI(279))
       P(685) = -3*P(7)-4*P(14)
       FR(280) = Is15s34*(-B025R-2*Cij125R(4,2)-C0125R*P(3)+Cij125R(1,
     -   1)*P(407))+Is34*((Is12*(-2*B035R-4*Cij135R(4,2)+C0135R*P(513
     -   )-Cij135R(1,1)*P(514)))/2.d0+(-5*C0235R-Cij125R(2,1)+14*Dij123
     -   5R(7,2)+2*(p2sq*Dij1235R(2,2)-(D01235R+Dij1235R(1,1))*P(3))-
     -   Dij1235R(4,2)*P(33)-Dij1235R(6,2)*P(436)+Dij1235R(2,1)*P(611
     -   )+Dij1235R(3,1)*P(612))/2.d0)+Is12*((Is45*(-2*(B012R+2*Cij145R
     -   (4,2))-Cij145R(1,1)*P(95)+C0145R*P(387)))/2.d0+(-5*C0345R-Cij1
     -   35R(2,1)+14*Dij1345R(7,2)-Dij1345R(6,2)*P(14)+2*(p3sq*Dij134
     -   5R(2,2)-D01345R*P(35))-Dij1345R(4,2)*P(198)-Dij1345R(3,1)*P(
     -   328)-Dij1345R(1,1)*P(638)+Dij1345R(2,1)*P(639))/2.d0)+(-5*D023
     -   45R+3*(Dij2345R(1,1)-Dij2345R(2,1))+38*EijR(11,2)+(EijR(2,1)
     -   -4*EijR(5,2))*P(33)-2*(EijR(1,1)*P(3)-EE0R*P(407))-EijR(4,1)
     -   *P(533)-EijR(3,1)*P(595)+EijR(2,2)*P(641)+EijR(3,2)*P(642)+E
     -   ijR(6,2)*P(643)-EijR(8,2)*P(644)+EijR(9,2)*P(645)+EijR(10,2)
     -   *P(685))/2.d0
       FI(280) = Is15s34*(-B025I-2*Cij125I(4,2)-C0125I*P(3)+Cij125I(1,
     -   1)*P(407))+Is34*((Is12*(-2*B035I-4*Cij135I(4,2)+C0135I*P(513
     -   )-Cij135I(1,1)*P(514)))/2.d0+(-5*C0235I-Cij125I(2,1)+14*Dij123
     -   5I(7,2)+2*(p2sq*Dij1235I(2,2)-(D01235I+Dij1235I(1,1))*P(3))-
     -   Dij1235I(4,2)*P(33)-Dij1235I(6,2)*P(436)+Dij1235I(2,1)*P(611
     -   )+Dij1235I(3,1)*P(612))/2.d0)+Is12*((Is45*(-2*(B012I+2*Cij145I
     -   (4,2))-Cij145I(1,1)*P(95)+C0145I*P(387)))/2.d0+(-5*C0345I-Cij1
     -   35I(2,1)+14*Dij1345I(7,2)-Dij1345I(6,2)*P(14)+2*(p3sq*Dij134
     -   5I(2,2)-D01345I*P(35))-Dij1345I(4,2)*P(198)-Dij1345I(3,1)*P(
     -   328)-Dij1345I(1,1)*P(638)+Dij1345I(2,1)*P(639))/2.d0)+(-5*D023
     -   45I+3*(Dij2345I(1,1)-Dij2345I(2,1))+38*EijI(11,2)+(EijI(2,1)
     -   -4*EijI(5,2))*P(33)-2*(EijI(1,1)*P(3)-EE0I*P(407))-EijI(4,1)
     -   *P(533)-EijI(3,1)*P(595)+EijI(2,2)*P(641)+EijI(3,2)*P(642)+E
     -   ijI(6,2)*P(643)-EijI(8,2)*P(644)+EijI(9,2)*P(645)+EijI(10,2)
     -   *P(685))/2.d0
       F(280)=DCMPLX(FR(280),FI(280))
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
       call FaFunctionNoAbemmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup3,
     -   p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5mup
     -   2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
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


       result = Fa(1)*SMB(1)+Fa(2)*SMB(2)+Fa(3)*SMB(3)+Fa(4)*SMB(4)+Fa
     -   (5)*SMB(5)+F(204)*SMB(6)+Fa(6)*SMB(7)+Fa(7)*SMB(8)+Fa(8)*SMB
     -   (9)+Fa(9)*SMB(10)+Fa(10)*SMB(11)+Fa(11)*SMB(12)+Fa(12)*SMB(1
     -   3)+Fa(13)*SMB(14)+Fa(14)*SMB(15)+F(280)*SMB(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       result =-result
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
       call FafunctionNoAbemmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup3,
     -   p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5mup
     -   2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
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
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*S
     -   MG(4)+Fa(5)*SMG(5)+F(204)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(280)*SMG(16)
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
       call FafunctionNoAbemmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup3,
     -   p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5mup
     -   2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
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
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*S
     -   MG(4)+Fa(5)*SMG(5)+F(204)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(280)*SMG(16)
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
       call FafunctionNoAbemmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup3,
     -   p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5mup
     -   2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
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
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*S
     -   MG(4)+Fa(5)*SMG(5)+F(204)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(280)*SMG(16)
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
       call FafunctionNoAbemmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup3,
     -   p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5mup
     -   2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
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
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*S
     -   MG(4)+Fa(5)*SMG(5)+F(204)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(280)*SMG(16)
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
       call FafunctionNoAbemmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup3,
     -   p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5mup
     -   2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
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
       resultgauge(2) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*S
     -   MG(4)+Fa(5)*SMG(5)+F(204)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(280)*SMG(16)
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
       call FafunctionNoAbemmm(p1mup2, p1mup3, p1mup5, p2mup2, p2mup3,
     -   p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5,p5mup
     -   2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
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
       resultgauge(3) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*S
     -   MG(4)+Fa(5)*SMG(5)+F(204)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(280)*SMG(16)
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
c************** Calling the Fa functions**********************************************************************
       subroutine FaFunctionNoAbemmm(p1mup2, p1mup3, p1mup5, p2mup2, p
     -   2mup3,p2mup5,p3mup2, p3mup3, p3mup5 , p4mup2, p4mup3, p4mup5
     -   ,p5mup2,p5mup3, p5mup5, mup2mup3, mup2mup5, mup3mup5,Fa)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup5, p2mup2, p2mup3, p2mup5, 
     -          p3mup2, p3mup3, p3mup5, p4mup2, p4mup3, p4mup5, 
     -          p5mup2, p5mup3, p5mup5
       Complex*16 mup2mup3,mup2mup5,mup3mup5
        common/FfunctionsNOABEmmm/F,P
      Complex*16 Fa(14),F(280)
      Real*8 P(685) 
        Fa(1) = p1mup3*(p1mup5*F(1)+p4mup5*F(2)+p3mup5*F(3))+p4mup3*(p3
     -   mup5*F(5)+p1mup5*F(6)+p4mup5*F(9))+p3mup3*(p3mup5*F(4)+p1mup
     -   5*F(7)+p4mup5*F(8)+p5mup5*F(10))+p5mup3*(p4mup5*F(12)+p3mup5
     -   *F(13)+p1mup5*F(14))+mup3mup5*F(16)+p5mup5*(p1mup3*F(11)+p4m
     -   up3*F(15)+p5mup3*F(17))
       Fa(2) = p3mup2*(p1mup5*F(20)+p4mup5*F(23)+p3mup5*F(24))+p4mup2*
     -   (p1mup5*F(18)+p4mup5*F(19)+p3mup5*F(25))+p1mup2*(p4mup5*F(21
     -   )+p3mup5*F(22)+p1mup5*F(26))+p5mup2*(p1mup5*F(27)+p3mup5*F(2
     -   8)+p4mup5*F(32))+mup2mup5*F(33)+p5mup5*(p4mup2*F(29)+p1mup2*
     -   F(30)+p3mup2*F(31)+p5mup2*F(34))
       Fa(3) = p1mup2*(p4mup3*F(39)+p1mup3*F(44)+p5mup3*F(45))+p3mup3*
     -   (p1mup2*F(38)+p4mup2*F(40)+p3mup2*F(41)+p5mup2*F(47))+p3mup2
     -   *(p1mup3*F(35)+p4mup3*F(43)+p5mup3*F(48))+p4mup2*(p1mup3*F(3
     -   6)+p4mup3*F(42)+p5mup3*F(49))+mup2mup3*F(50)+p5mup2*(p1mup3*
     -   F(37)+p4mup3*F(46)+p5mup3*F(51))
       Fa(4) = p1mup2*(p1mup3*(p1mup5*F(52)+p3mup5*F(68)+p4mup5*F(84))
     -   +p4mup3*(p1mup5*F(60)+p3mup5*F(76)+p4mup5*F(92))+p5mup3*(p1m
     -   up5*F(64)+p3mup5*F(80)+p4mup5*F(96)))+p3mup2*(p1mup3*(p1mup5
     -   *F(53)+p3mup5*F(69)+p4mup5*F(85))+p4mup3*(p1mup5*F(61)+p3mup
     -   5*F(77)+p4mup5*F(93))+p5mup3*(p1mup5*F(65)+p3mup5*F(81)+p4mu
     -   p5*F(97)))+p4mup2*(p1mup3*(p1mup5*F(54)+p3mup5*F(70)+p4mup5*
     -   F(86))+p4mup3*(p1mup5*F(62)+p3mup5*F(78)+p4mup5*F(94))+p5mup
     -   3*(p1mup5*F(66)+p3mup5*F(82)+p4mup5*F(98)))+p5mup2*(p1mup3*(
     -   p1mup5*F(55)+p3mup5*F(71)+p4mup5*F(87))+p4mup3*(p1mup5*F(63)
     -   +p3mup5*F(79)+p4mup5*F(95))+p5mup3*(p1mup5*F(67)+p3mup5*F(83
     -   )+p4mup5*F(99)))+p3mup3*(p1mup2*(p1mup5*F(56)+p3mup5*F(72)+p
     -   4mup5*F(88))+p3mup2*(p1mup5*F(57)+p3mup5*F(73)+p4mup5*F(89))
     -   +p4mup2*(p1mup5*F(58)+p3mup5*F(74)+p4mup5*F(90))+p5mup2*(p1m
     -   up5*F(59)+p3mup5*F(75)+p4mup5*F(91))+p5mup5*(p1mup2*F(104)+p
     -   3mup2*F(105)+p4mup2*F(106)+p5mup2*F(107))+mup2mup5*F(112))+m
     -   up2mup3*(p3mup5*F(114)+p1mup5*F(116)+p4mup5*F(119))+mup3mup5
     -   *(p1mup2*F(115)+p3mup2*F(117)+p4mup2*F(118)+p5mup2*F(124))+m
     -   up2mup5*(p4mup3*F(111)+p1mup3*F(113)+p5mup3*F(126))+p5mup5*(
     -   p3mup2*(p1mup3*F(101)+p4mup3*F(109)+p5mup3*F(120))+p1mup2*(p
     -   1mup3*F(100)+p4mup3*F(108)+p5mup3*F(121))+p4mup2*(p1mup3*F(1
     -   02)+p4mup3*F(110)+p5mup3*F(122))+mup2mup3*F(125)+p5mup2*(p1m
     -   up3*F(103)+p4mup3*F(123)+p5mup3*F(127)))
       Fa(5) = p1mup2*(p1mup3*(p1mup5*F(130)+p3mup5*F(146)+p4mup5*F(16
     -   2))+p4mup3*(p1mup5*F(138)+p3mup5*F(154)+p4mup5*F(170))+p5mup
     -   3*(p1mup5*F(142)+p3mup5*F(158)+p4mup5*F(174)))+p3mup2*(p1mup
     -   3*(p1mup5*F(131)+p3mup5*F(147)+p4mup5*F(163))+p4mup3*(p1mup5
     -   *F(139)+p3mup5*F(155)+p4mup5*F(171))+p5mup3*(p1mup5*F(143)+p
     -   3mup5*F(159)+p4mup5*F(175)))+p4mup2*(p1mup3*(p1mup5*F(132)+p
     -   3mup5*F(148)+p4mup5*F(164))+p4mup3*(p1mup5*F(140)+p3mup5*F(1
     -   56)+p4mup5*F(172))+p5mup3*(p1mup5*F(144)+p3mup5*F(160)+p4mup
     -   5*F(176)))+p5mup2*(p1mup3*(p1mup5*F(133)+p3mup5*F(149)+p4mup
     -   5*F(165))+p4mup3*(p1mup5*F(141)+p3mup5*F(157)+p4mup5*F(173))
     -   +p5mup3*(p1mup5*F(145)+p3mup5*F(161)+p4mup5*F(177)))+mup3mup
     -   5*(p3mup2*F(128)+p5mup2*F(129)+p1mup2*F(182)+p4mup2*F(195))+
     -   mup2mup3*(p3mup5*F(190)+p1mup5*F(193)+p4mup5*F(196))+p3mup3*
     -   (p1mup2*(p1mup5*F(134)+p3mup5*F(150)+p4mup5*F(166))+p3mup2*(
     -   p1mup5*F(135)+p3mup5*F(151)+p4mup5*F(167))+p4mup2*(p1mup5*F(
     -   136)+p3mup5*F(152)+p4mup5*F(168))+p5mup2*(p1mup5*F(137)+p3mu
     -   p5*F(153)+p4mup5*F(169))+mup2mup5*F(191)+p5mup5*(p3mup2*F(18
     -   3)+p1mup2*F(187)+p4mup2*F(189)+p5mup2*F(197)))+mup2mup5*(p1m
     -   up3*F(185)+p4mup3*F(194)+p5mup3*F(202))+p5mup5*(mup2mup3*F(1
     -   92)+p1mup2*(p1mup3*F(178)+p4mup3*F(186)+p5mup3*F(198))+p4mup
     -   2*(p1mup3*F(180)+p4mup3*F(188)+p5mup3*F(199))+p3mup2*(p1mup3
     -   *F(179)+p4mup3*F(184)+p5mup3*F(200))+p5mup2*(p1mup3*F(181)+p
     -   4mup3*F(201)+p5mup3*F(203)))
       Fa(6) = p3mup5*F(205)+p1mup5*F(206)+p4mup5*F(207)+p5mup5*F(208)
       Fa(7) = p4mup3*F(209)+p3mup3*F(210)+p1mup3*F(211)+p5mup3*F(212)
       Fa(8) = p1mup2*F(213)+p4mup2*F(214)+p3mup2*F(215)+p5mup2*F(216)
       Fa(9) = p1mup3*(p1mup5*F(217)+p3mup5*F(221)+p4mup5*F(225))+p4mu
     -   p3*(p1mup5*F(219)+p3mup5*F(223)+p4mup5*F(227))+p5mup3*(p1mup
     -   5*F(220)+p3mup5*F(224)+p4mup5*F(228))+p3mup3*(p1mup5*F(218)+
     -   p3mup5*F(222)+p4mup5*F(226)+p5mup5*F(229))+mup3mup5*F(232)+p
     -   5mup5*(p4mup3*F(230)+p1mup3*F(231)+p5mup3*F(233))
       Fa(10) = p1mup2*(p1mup5*F(234)+p3mup5*F(238)+p4mup5*F(242))+p3m
     -   up2*(p1mup5*F(235)+p3mup5*F(239)+p4mup5*F(243))+p4mup2*(p1mu
     -   p5*F(236)+p3mup5*F(240)+p4mup5*F(244))+p5mup2*(p1mup5*F(237)
     -   +p3mup5*F(241)+p4mup5*F(245))+mup2mup5*F(249)+p5mup5*(p3mup2
     -   *F(246)+p1mup2*F(247)+p4mup2*F(248)+p5mup2*F(250))
       Fa(11) = p3mup3*(p1mup2*F(255)+p3mup2*F(256)+p4mup2*F(257)+p5mu
     -   p2*F(258))+p1mup2*(p1mup3*F(251)+p4mup3*F(259)+p5mup3*F(263)
     -   )+p3mup2*(p1mup3*F(252)+p4mup3*F(260)+p5mup3*F(264))+p4mup2*
     -   (p1mup3*F(253)+p4mup3*F(261)+p5mup3*F(265))+p5mup2*(p1mup3*F
     -   (254)+p4mup3*F(262)+p5mup3*F(266))+mup2mup3*F(267)
       Fa(12) = p3mup5*F(268)+p4mup5*F(269)+p1mup5*F(270)+p5mup5*F(271
     -   )
       Fa(13) = p1mup3*F(272)+p3mup3*F(273)+p4mup3*F(274)+p5mup3*F(275
     -   )
       Fa(14) = p3mup2*F(276)+p1mup2*F(277)+p4mup2*F(278)+p5mup2*F(279
     -   )
       Return
       End
