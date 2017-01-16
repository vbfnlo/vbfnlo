       subroutine penlinemmm(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,mup3
     -   ,mup4,alpha,musqIn,gaugetest,comp,resultgauge,result,resultgau
     -   geb,resultb)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 13/2/2008
c Modified:30/7/2008
c ************************************************************************************
c determine the  finite virtual corrections along the quark line i.e 
c the sum of all virtual corrections to the Born diagram 
c psi(p1)    ---->------->------->------->---   bar_psi(p5)
c                   $       $       $              
c                   $       $       $              
c                   $       $       $              
c                   V       V       V              
c                   $       $       $              
c                   $       $       $              
c                  p2       p3      p4             
c                 mu_p2    mu_p3   mu_p4           
c Note: The vertices are just Gamma^(mu_p2)..terms. So the correct
c factor should be added by hand for each boson. i.e, Without adding
c anything it represents pp->3 photons up to the electromagnetic charge 
c (for uu->3 photons is missing (2/3)^3) from the coupling. For W and Z,
c we have to account for the Diracgamma_5 that are missing in the vertex
c Note: To make it shorter in the promgram: mu_p2,...->mup2,... 
c Notation of External momenta: p1+p2+p3+p4+p5=0 
c mu_p2,mu_p3,mu_p4, should be think as external current 
c alpha is the helicity of the initial spinor 
c musq is the renormalization scale energy  
c gaugetest,integer value.Different gauge test: 
c gaugetest=0 should give zero 
c gaugetest=1 give the result after replacing mu_p2 with p2 
c gaugetest=2 give the result after replacing mu_p3 with p3 
c gaugetest=3 give the result after replacing mu_p4 with p4 
c gaugetest=4 calculates the three different gaugue tests 
c gaugetest>5 no calculation of gauge test in case you
c have already done it before(it safes times)
c comp: integer value.The first time called with p1...p5, comp=1
c ATTENTION: ONLY!!!If you have to call the subroutine consecutively with the same arguments
c(p1,p2,p3,p4,p5). Then, comp=-1 (it safes 4000 lines of code) 
c This applies when you have for examples the same diagram for an off-shell photon
c and a Z boson. The differences are in the coupling and  the part that depends on the
c polarization vector that are calculated at the end of this program.
c resultgauge and resultgaugeb are arrays of dimension three. 
c In case you use gaugetest=(0,3). The result is given in the first argument:
c resultgauge(1) and resultgaugeb(1).The argument 2,and 3 is set to zero 
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
       Real*8   p1sq, p1p2, p1p3, p1p4, p1p5, p2sq, p2p3, p2p4, 
     -          p2p5, p3sq, p3p4, p3p5, p4sq, p4p5, p5sq
       Real*8   s12, s13, s14, s15, s23, s24, s25, s34, s35, s45
       Complex*16   p1mup2, p1mup3, p1mup4, p2mup2, p2mup3, p2mup4, 
     -          p3mup2, p3mup3, p3mup4, p4mup2, p4mup3, p4mup4, 
     -          p5mup2, p5mup3, p5mup4
       Complex*16 mup2mup3,mup2mup4,mup3mup4
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
       Complex*16 SMB(16),SMG(16) ,Fa(14),F(271),K(3)
       Real*8 FI(271),FR(271),KI(3),  KR(3)
       Complex*16 barpsi_p5(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3) 
       Complex*16 SC1c,SC3ccc,SC3rcc,SC3rrc,SC5rrccc,SC1r,dotrc,dotcc,result
       Complex*16 resultgauge(3),resultgaugeb(3),resultb 
       Real*8 musq,p5t,P(592),Is12,Is45,Is12s45 ,musqIn
       EXTERNAL   dotrc,dotcc,SC1c,SC1r,SC3ccc,SC3rcc,SC3rrc,SC5rrccc
       Integer alpha,comp,gaugetest
       common/Ffunctionsmmm/F,P
       common/Kfunctionsmmm/K
       SAVE/Ffunctionsmmm/
       SAVE/Kfunctionsmmm/
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
       Is12=1.d0/s12
       Is45=1.d0/s45
       Is12s45=1.d0/(s12*s45)
       p5t=p1sq+p2sq+p3sq+p4sq+2d0*(p1p2+p1p3+p1p4+p2p3+p2p4+p3p4)
       If(abs(p5t).gt.1d-07) then
       If(p4sq.gt.1d-05)then
       p4(0)=-(p1(0)+p2(0)+p3(0)+p5(0))
       p4(1)=-(p1(1)+p2(1)+p3(1)+p5(1))
       p4(2)=-(p1(2)+p2(2)+p3(2)+p5(2))
       p4(3)=-(p1(3)+p2(3)+p3(3)+p5(3))
       write(*,*)"Momemtum is not conserved: p5sq=0. here=",p5t
       write(*,*)"Corrected by p4=-(p1+p2+p3+p5)"
       elseIf(p2sq.gt.1d-05)then
       p2(0)=-(p1(0)+p4(0)+p3(0)+p5(0))
       p2(1)=-(p1(1)+p4(1)+p3(1)+p5(1))
       p2(2)=-(p1(2)+p4(2)+p3(2)+p5(2))
       p2(3)=-(p1(3)+p4(3)+p3(3)+p5(3))
       write(*,*)"Momemtum is not conserved: p5sq=0. here=",p5t
       write(*,*)"Corrected by p2=-(p1+p3+p4+p5)"
       elseIf(p3sq.gt.1d-05)then
       p3(0)=-(p1(0)+p4(0)+p2(0)+p5(0))
       p3(1)=-(p1(1)+p4(1)+p2(1)+p5(1))
       p3(2)=-(p1(2)+p4(2)+p2(2)+p5(2))
       p3(3)=-(p1(3)+p4(3)+p2(3)+p5(3))
       write(*,*)"Momemtum is not conserved: p5sq=0. here=",p5t
       write(*,*) "Corrected by p3=-(p1+p2+p4+p5)"
       else
       write(*,*)"Momemtum is not conserved: p5sq=0 here=",p5t
       endif
       endif
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
       B012=0d0              
       B013=B0fin(s12,musq)   
       B014=B0fin(s45,musq)   
       B015=0d0              
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
       C0145=C0fin(s45,p4sq,0d0,musq) 
       C0125=C0fin(0d0,s15,0d0,musq) 
       C0123=C0fin(0d0,p2sq,s12,musq) 
       C0124=C0fin(0d0,s23,s45,musq) 
       C0135=C0fin(s12,s34,0d0,musq) 
       C0134=C0fin(s12,p3sq,s45,musq) 
       C0234=C0fin(p2sq,p3sq,s23,musq) 
       C0235=C0fin(p2sq,s34,s15,musq) 
       C0245=C0fin(s23,p4sq,s15,musq) 
       C0345=C0fin(p3sq,p4sq,s34,musq) 
c************************************************************************************
c************************************************************************************
       call tens_red3_new_Re_Com(s45,p4sq,0d0,B045,B015,B014,C0145,C01
     -   45R,C0145I,Cij145R,Cij145I) 
       call tens_red3_new_Re_Com(0d0,s15,0d0,B025,B015,B012,C0125,C012
     -   5R,C0125I,Cij125R,Cij125I) 
       call tens_red3_new_Re_Com(0d0,p2sq,s12,B023,B013,B012,C0123,C01
     -   23R,C0123I,Cij123R,Cij123I) 
       call tens_red3_new_Re_Com(0d0,s23,s45,B024,B014,B012,C0124,C012
     -   4R,C0124I,Cij124R,Cij124I) 
       call tens_red3_new_Re_Com(s12,s34,0d0,B035,B015,B013,C0135,C013
     -   5R,C0135I,Cij135R,Cij135I) 
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
       EE0=E0fin(0d0,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
c,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
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
       P(1) = p2sq+p3sq-s23
       P(2) = 2-Is12*P(1)
       P(3) = 3-Is12*P(1)
       P(4) = 3-2*Is12*P(1)
       P(5) = -p2sq-p3sq+s12+s23
       P(6) = p3sq-s23
       P(7) = p2sq-p3sq-s12-s15+s23+s34
       P(8) = p2sq-s12-s15+s34
       P(9) = p2sq+p3sq-s12-s23
       FR(1) = 16*Dij1345R(5,2)*P(2)+8*(-Dij1345R(1,1)-Dij1345R(1,3)-p
     -   2sq*(EijR(2,1)+EijR(5,3)+EijR(15,3)+2*(EijR(5,2)-EijR(9,2)-E
     -   ijR(18,3)))-2*(Dij1345R(1,2)+EijR(11,2)+EijR(36,4)+EijR(39,4
     -   )+2*(EijR(21,3)-EijR(24,3)-EijR(42,4)))+(EijR(3,1)+EijR(6,2)
     -   -EijR(10,2))*P(1)+Dij1345R(5,3)*P(3)-Dij1345R(7,3)*P(4)-EijR
     -   (4,1)*P(6)-(EijR(4,2)-EijR(7,2))*P(7)+(EijR(4,3)+EijR(7,3)-2
     -   *EijR(14,3))*P(8)+Is12*((Dij1345R(2,1)+Dij1345R(4,3)+Dij1345
     -   R(9,3)+2*(Dij1345R(4,2)-Dij1345R(6,2)-Dij1345R(10,3)))*P(1)+
     -   (Dij1345R(3,1)+Dij1345R(3,3))*P(5)+2*Dij1345R(3,2)*P(9)))
       FI(1) = 16*Dij1345I(5,2)*P(2)+8*(-Dij1345I(1,1)-Dij1345I(1,3)-p
     -   2sq*(EijI(2,1)+EijI(5,3)+EijI(15,3)+2*(EijI(5,2)-EijI(9,2)-E
     -   ijI(18,3)))-2*(Dij1345I(1,2)+EijI(11,2)+EijI(36,4)+EijI(39,4
     -   )+2*(EijI(21,3)-EijI(24,3)-EijI(42,4)))+(EijI(3,1)+EijI(6,2)
     -   -EijI(10,2))*P(1)+Dij1345I(5,3)*P(3)-Dij1345I(7,3)*P(4)-EijI
     -   (4,1)*P(6)-(EijI(4,2)-EijI(7,2))*P(7)+(EijI(4,3)+EijI(7,3)-2
     -   *EijI(14,3))*P(8)+Is12*((Dij1345I(2,1)+Dij1345I(4,3)+Dij1345
     -   I(9,3)+2*(Dij1345I(4,2)-Dij1345I(6,2)-Dij1345I(10,3)))*P(1)+
     -   (Dij1345I(3,1)+Dij1345I(3,3))*P(5)+2*Dij1345I(3,2)*P(9)))
       F(1)=DCMPLX(FR(1),FI(1))
       P(10) = p2sq-s23
       P(11) = -p2sq+s12+s23
       P(12) = p3sq-p4sq-3*s12+s34+s45+2*P(10)
       P(13) = p3sq+2*P(10)
       P(14) = p3sq-p4sq-5*s12+s34+s45+2*P(10)
       P(15) = p4sq-s45
       P(16) = p2sq+p3sq-s15+s23-s34+2*P(15)
       P(17) = -p2sq-p3sq+s12+s15+s45
       P(18) = p2sq-s12
       P(19) = p3sq+p4sq-s45
       P(20) = p2sq-s15-s34+2*P(19)
       P(21) = 2*p2sq-p3sq-s12+s23
       P(22) = -p2sq+p3sq+s23
       P(23) = p3sq-s15-s45+2*P(18)
       P(24) = p2sq+p3sq-2*s12-s15-s45
       P(25) = p3sq+s12-s45
       P(26) = 2*p2sq+p3sq-3*s12-s15-s45
       FR(2) = -8*(Dij1345R(1,1)+Dij1345R(1,3)+Dij2345R(1,1)-Dij2345R(
     -   3,1)-Dij2345R(3,2)+Dij2345R(5,2)+p2sq*(EijR(2,1)+EijR(2,2))-
     -   s23*EijR(4,1)+s12*(-EijR(8,3)+EijR(10,3))+2*(Dij1345R(1,2)+p
     -   3sq*EijR(6,2)+2*EijR(11,2)+EijR(39,4)+EijR(40,4)-EijR(42,4)-
     -   EijR(44,4))-EijR(8,2)*P(1)-Dij1345R(5,3)*P(3)+Dij1345R(7,3)*
     -   P(4)-EijR(3,1)*P(10)+Is12*(Cij145R(1,1)+Cij145R(1,2)-Cij145R
     -   (2,1)+Cij145R(2,2)-2*(Cij145R(3,2)-Dij1345R(7,2)-Dij1345R(11
     -   ,3)+Dij1345R(13,3))-(Dij1345R(4,3)+Dij1345R(9,3)-2*Dij1345R(
     -   10,3))*P(1)-Dij1345R(3,3)*P(5)-Dij1345R(2,1)*P(10)-Dij1345R(
     -   3,1)*P(11)-Dij1345R(3,2)*P(12)+(-Dij1345R(4,2)+Dij1345R(6,2)
     -   )*P(13)+Dij1345R(5,2)*P(14))+EijR(4,2)*P(16)+EijR(4,3)*P(17)
     -   +(EijR(5,2)+EijR(5,3)-EijR(7,3))*P(18)-EijR(7,2)*P(20)-EijR(
     -   9,2)*P(21)-EijR(10,2)*P(22)+EijR(14,3)*P(23)+EijR(15,3)*P(24
     -   )+(EijR(16,3)+EijR(17,3)-EijR(19,3)-EijR(20,3))*P(25)-EijR(1
     -   8,3)*P(26))
       FI(2) = -8*(Dij1345I(1,1)+Dij1345I(1,3)+Dij2345I(1,1)-Dij2345I(
     -   3,1)-Dij2345I(3,2)+Dij2345I(5,2)+p2sq*(EijI(2,1)+EijI(2,2))-
     -   s23*EijI(4,1)+s12*(-EijI(8,3)+EijI(10,3))+2*(Dij1345I(1,2)+p
     -   3sq*EijI(6,2)+2*EijI(11,2)+EijI(39,4)+EijI(40,4)-EijI(42,4)-
     -   EijI(44,4))-EijI(8,2)*P(1)-Dij1345I(5,3)*P(3)+Dij1345I(7,3)*
     -   P(4)-EijI(3,1)*P(10)+Is12*(Cij145I(1,1)+Cij145I(1,2)-Cij145I
     -   (2,1)+Cij145I(2,2)-2*(Cij145I(3,2)-Dij1345I(7,2)-Dij1345I(11
     -   ,3)+Dij1345I(13,3))-(Dij1345I(4,3)+Dij1345I(9,3)-2*Dij1345I(
     -   10,3))*P(1)-Dij1345I(3,3)*P(5)-Dij1345I(2,1)*P(10)-Dij1345I(
     -   3,1)*P(11)-Dij1345I(3,2)*P(12)+(-Dij1345I(4,2)+Dij1345I(6,2)
     -   )*P(13)+Dij1345I(5,2)*P(14))+EijI(4,2)*P(16)+EijI(4,3)*P(17)
     -   +(EijI(5,2)+EijI(5,3)-EijI(7,3))*P(18)-EijI(7,2)*P(20)-EijI(
     -   9,2)*P(21)-EijI(10,2)*P(22)+EijI(14,3)*P(23)+EijI(15,3)*P(24
     -   )+(EijI(16,3)+EijI(17,3)-EijI(19,3)-EijI(20,3))*P(25)-EijI(1
     -   8,3)*P(26))
       F(2)=DCMPLX(FR(2),FI(2))
       P(27) = 1-2*Is12*P(1)
       P(28) = s12+s15-s34
       P(29) = p2sq+p3sq+s15-s23-s34
       P(30) = p3sq-s45
       P(31) = p2sq-s15+2*P(30)
       FR(3) = 8*(-Dij1345R(3,2)+Dij1345R(4,2)-Dij1345R(7,3)+Dij1345R(
     -   10,3)+Dij2345R(3,2)-Dij2345R(6,2)+p2sq*EijR(8,2)+s12*EijR(15
     -   ,3)-s12*EijR(20,3)+2*(EijR(23,3)-EijR(39,4)+EijR(45,4))-EijR
     -   (3,2)*P(1)-(EijR(3,1)-EijR(4,1))*P(9)+Is12*(-((Dij1345R(2,2)
     -   -Dij1345R(6,2)+Dij1345R(8,3))*P(1))+(Dij1345R(2,1)+Dij1345R(
     -   3,3))*P(5)+Dij1345R(3,1)*P(9))-EijR(4,3)*P(17)-(EijR(6,2)+Ei
     -   jR(14,3)-EijR(19,3))*P(18)+EijR(13,3)*P(25)-Dij1345R(9,3)*P(
     -   27)-EijR(4,2)*P(28)+EijR(10,2)*P(29)-EijR(16,3)*P(31))
       FI(3) = 8*(-Dij1345I(3,2)+Dij1345I(4,2)-Dij1345I(7,3)+Dij1345I(
     -   10,3)+Dij2345I(3,2)-Dij2345I(6,2)+p2sq*EijI(8,2)+s12*EijI(15
     -   ,3)-s12*EijI(20,3)+2*(EijI(23,3)-EijI(39,4)+EijI(45,4))-EijI
     -   (3,2)*P(1)-(EijI(3,1)-EijI(4,1))*P(9)+Is12*(-((Dij1345I(2,2)
     -   -Dij1345I(6,2)+Dij1345I(8,3))*P(1))+(Dij1345I(2,1)+Dij1345I(
     -   3,3))*P(5)+Dij1345I(3,1)*P(9))-EijI(4,3)*P(17)-(EijI(6,2)+Ei
     -   jI(14,3)-EijI(19,3))*P(18)+EijI(13,3)*P(25)-Dij1345I(9,3)*P(
     -   27)-EijI(4,2)*P(28)+EijI(10,2)*P(29)-EijI(16,3)*P(31))
       F(3)=DCMPLX(FR(3),FI(3))
       P(32) = p2sq+p3sq-2*s12-s23
       P(33) = p2sq+p3sq-s34-s45
       P(34) = p2sq+2*p3sq-s23-s45
       P(35) = -p2sq-2*p3sq+s23+s34+s45
       P(36) = s23-s45
       P(37) = p2sq+p3sq-s45
       P(38) = -s12-s15+2*P(37)
       FR(4) = 8*(-D02345R-Dij1345R(1,1)-Dij1345R(1,2)+Dij1345R(3,1)-D
     -   ij1345R(4,2)-Dij1345R(4,3)+3*Dij1345R(5,2)+Dij1345R(5,3)-Dij
     -   2345R(2,1)+Dij2345R(3,2)-Dij2345R(6,2)-p2sq*EijR(8,2)+s12*(E
     -   ijR(5,2)+EijR(15,3)+EijR(17,3)-EijR(18,3)-EijR(20,3))+2*(Eij
     -   R(11,2)-EijR(21,3)+EijR(24,3)-EijR(39,4)-EijR(41,4)+EijR(42,
     -   4)+EijR(45,4))+EijR(3,2)*P(1)-EijR(4,3)*P(17)-(EijR(6,3)-Eij
     -   R(7,3)-EijR(9,2))*P(18)-EijR(14,3)*P(23)-(EijR(11,3)-EijR(13
     -   ,3))*P(25)-(-Dij1345R(6,2)+Dij1345R(9,3))*P(27)-EijR(16,3)*P
     -   (31)+Is12*((Dij1345R(2,2)+Dij1345R(6,3)-Dij1345R(8,3))*P(1)+
     -   (Dij1345R(3,3)+2*Dij1345R(10,3))*P(5)+(Dij1345R(3,2)+Dij1345
     -   R(7,3))*P(32))-EijR(4,2)*P(33)-EijR(6,2)*P(34)-EijR(7,2)*P(3
     -   5)+EijR(10,2)*P(36)+EijR(19,3)*P(38))
       FI(4) = 8*(-D02345I-Dij1345I(1,1)-Dij1345I(1,2)+Dij1345I(3,1)-D
     -   ij1345I(4,2)-Dij1345I(4,3)+3*Dij1345I(5,2)+Dij1345I(5,3)-Dij
     -   2345I(2,1)+Dij2345I(3,2)-Dij2345I(6,2)-p2sq*EijI(8,2)+s12*(E
     -   ijI(5,2)+EijI(15,3)+EijI(17,3)-EijI(18,3)-EijI(20,3))+2*(Eij
     -   I(11,2)-EijI(21,3)+EijI(24,3)-EijI(39,4)-EijI(41,4)+EijI(42,
     -   4)+EijI(45,4))+EijI(3,2)*P(1)-EijI(4,3)*P(17)-(EijI(6,3)-Eij
     -   I(7,3)-EijI(9,2))*P(18)-EijI(14,3)*P(23)-(EijI(11,3)-EijI(13
     -   ,3))*P(25)-(-Dij1345I(6,2)+Dij1345I(9,3))*P(27)-EijI(16,3)*P
     -   (31)+Is12*((Dij1345I(2,2)+Dij1345I(6,3)-Dij1345I(8,3))*P(1)+
     -   (Dij1345I(3,3)+2*Dij1345I(10,3))*P(5)+(Dij1345I(3,2)+Dij1345
     -   I(7,3))*P(32))-EijI(4,2)*P(33)-EijI(6,2)*P(34)-EijI(7,2)*P(3
     -   5)+EijI(10,2)*P(36)+EijI(19,3)*P(38))
       F(4)=DCMPLX(FR(4),FI(4))
       P(39) = 1-Is12*s34
       P(40) = 1-Is45*s23
       P(41) = 1-Is45*s12
       P(42) = p3sq-s12
       P(43) = p3sq-s12+s45
       P(44) = Is12-Is45
       P(45) = p2sq-s12-s23+s45
       P(46) = -p2sq+s12+s23-s45
       P(47) = p3sq-3*s12+s34+s45+2*P(10)
       P(48) = p3sq-s12+s45+2*P(10)
       P(49) = p3sq-4*s12+s34+s45+2*P(10)
       P(50) = p2sq-s12+s15
       P(51) = p2sq-2*s12
       P(52) = s12-s45
       P(53) = p2sq-p4sq-s23+s34-2*P(52)
       P(54) = p4sq+s12+s15+s23-s34-2*s45
       P(55) = p2sq+p3sq+p4sq-s15+s23-s34-3*s45
       P(56) = -p2sq+p3sq+p4sq+s23-s34+3*P(52)
       P(57) = p3sq+p4sq-s15+s23-s34-3*s45+2*P(18)
       P(58) = p2sq+p3sq-3*s12-s15-s45
       P(59) = p2sq+p3sq-s15-s45
       P(60) = -3*s12+2*P(59)
       FR(5) = 4*Is12s45*(3*B013R+B014R-4*(Cij123R(4,2)+Cij134R(4,2)+C
     -   ij145R(4,2))+2*(p2sq*Cij123R(2,1)-(C0145R+Cij145R(2,1))*P(15
     -   )-(C0123R+Cij123R(1,1))*P(18)-C0134R*P(42)-Cij134R(2,1)*P(43
     -   )))-8*(Dij1345R(1,2)+Dij1345R(1,3)-Dij2345R(1,1)-Dij2345R(1,
     -   2)+Dij2345R(3,1)-Dij2345R(3,2)+s15*(EE0R+EijR(1,1))-s12*EijR
     -   (2,3)+2*(Dij2345R(5,2)-EijR(11,2)-EijR(22,3)+EijR(24,3)+EijR
     -   (37,4)+EijR(39,4)-2*EijR(44,4))-Dij1345R(5,3)*P(3)+Dij1345R(
     -   7,3)*P(4)+EijR(4,3)*P(17)+(EijR(8,3)+EijR(14,3)-2*EijR(18,3)
     -   )*P(18)+Is45*(Cij134R(1,2)-s12*Dij1234R(2,1)-p2sq*Dij1234R(2
     -   ,2)+2*Dij1234R(12,3)+Dij1234R(4,2)*P(18))+(EijR(9,3)+EijR(16
     -   ,3)-2*EijR(20,3))*P(25)-(D01345R+Dij1345R(1,1))*P(39)-(D0123
     -   4R+Dij1234R(1,1))*P(40)-(Dij1234R(3,1)+Dij1234R(6,2))*P(41)+
     -   Cij134R(3,2)*P(44)+Is12*(Cij145R(1,2)+Cij145R(2,2)-2*(Cij145
     -   R(3,2)-Dij1345R(7,2))+4*(Dij1345R(11,3)-Dij1345R(13,3))-(Dij
     -   1345R(4,3)+Dij1345R(9,3)-2*Dij1345R(10,3))*P(1)-Dij1345R(3,3
     -   )*P(5)-Dij1345R(2,1)*P(45)-Dij1345R(3,1)*P(46)-Dij1345R(3,2)
     -   *P(47)+(-Dij1345R(4,2)+Dij1345R(6,2))*P(48)+Dij1345R(5,2)*P(
     -   49))+(EijR(2,1)+EijR(5,2)-EijR(7,2))*P(50)+EijR(2,2)*P(51)-E
     -   ijR(3,1)*P(53)-EijR(4,1)*P(54)+EijR(4,2)*P(55)-(-EijR(8,2)+E
     -   ijR(10,2))*P(56)-EijR(9,2)*P(57)-EijR(10,3)*P(58)+EijR(15,3)
     -   *P(60))
       FI(5) = 4*Is12s45*(3*B013I+B014I-4*(Cij123I(4,2)+Cij134I(4,2)+C
     -   ij145I(4,2))+2*(p2sq*Cij123I(2,1)-(C0145I+Cij145I(2,1))*P(15
     -   )-(C0123I+Cij123I(1,1))*P(18)-C0134I*P(42)-Cij134I(2,1)*P(43
     -   )))-8*(Dij1345I(1,2)+Dij1345I(1,3)-Dij2345I(1,1)-Dij2345I(1,
     -   2)+Dij2345I(3,1)-Dij2345I(3,2)+s15*(EE0I+EijI(1,1))-s12*EijI
     -   (2,3)+2*(Dij2345I(5,2)-EijI(11,2)-EijI(22,3)+EijI(24,3)+EijI
     -   (37,4)+EijI(39,4)-2*EijI(44,4))-Dij1345I(5,3)*P(3)+Dij1345I(
     -   7,3)*P(4)+EijI(4,3)*P(17)+(EijI(8,3)+EijI(14,3)-2*EijI(18,3)
     -   )*P(18)+Is45*(Cij134I(1,2)-s12*Dij1234I(2,1)-p2sq*Dij1234I(2
     -   ,2)+2*Dij1234I(12,3)+Dij1234I(4,2)*P(18))+(EijI(9,3)+EijI(16
     -   ,3)-2*EijI(20,3))*P(25)-(D01345I+Dij1345I(1,1))*P(39)-(D0123
     -   4I+Dij1234I(1,1))*P(40)-(Dij1234I(3,1)+Dij1234I(6,2))*P(41)+
     -   Cij134I(3,2)*P(44)+Is12*(Cij145I(1,2)+Cij145I(2,2)-2*(Cij145
     -   I(3,2)-Dij1345I(7,2))+4*(Dij1345I(11,3)-Dij1345I(13,3))-(Dij
     -   1345I(4,3)+Dij1345I(9,3)-2*Dij1345I(10,3))*P(1)-Dij1345I(3,3
     -   )*P(5)-Dij1345I(2,1)*P(45)-Dij1345I(3,1)*P(46)-Dij1345I(3,2)
     -   *P(47)+(-Dij1345I(4,2)+Dij1345I(6,2))*P(48)+Dij1345I(5,2)*P(
     -   49))+(EijI(2,1)+EijI(5,2)-EijI(7,2))*P(50)+EijI(2,2)*P(51)-E
     -   ijI(3,1)*P(53)-EijI(4,1)*P(54)+EijI(4,2)*P(55)-(-EijI(8,2)+E
     -   ijI(10,2))*P(56)-EijI(9,2)*P(57)-EijI(10,3)*P(58)+EijI(15,3)
     -   *P(60))
       F(5)=DCMPLX(FR(5),FI(5))
       P(61) = -p3sq+s12+s45
       P(62) = p2sq+p3sq+p4sq-s23-s34
       P(63) = p3sq+p4sq-3*s12+s34+s45+2*P(10)
       P(64) = p2sq-2*s12-s23+s34+s45
       P(65) = p3sq+p4sq
       P(66) = p2sq-s12-s23
       P(67) = s45+2*P(65)+3*P(66)
       P(68) = p2sq+s12
       P(69) = p2sq-2*s12-s15-s23+s34+s45
       P(70) = s12+s15+s23-s34-s45
       P(71) = p2sq+p3sq-s15+s23-s34-4*s45
       P(72) = p2sq-s12-2*P(36)
       P(73) = p2sq-s12+2*P(36)
       P(74) = p2sq+2*p3sq+s12+s15-s23-s34-s45
       P(75) = p2sq+2*p3sq-3*s12-s23-s34-s45
       P(76) = 4*s12+s23-s34-5*s45
       P(77) = p2sq-s12-s15+2*P(30)
       FR(6) = 4*Is12s45*(3*B013R+B014R-4*(Cij123R(4,2)+Cij134R(4,2)+C
     -   ij145R(4,2))+2*(p2sq*Cij123R(2,1)-(C0145R+Cij145R(2,1))*P(15
     -   )-(C0123R+Cij123R(1,1))*P(18)-C0134R*P(42)+Cij134R(2,1)*P(61
     -   )))-8*(D02345R+Dij1345R(4,3)-Dij1345R(5,3)-Dij2345R(1,1)-Dij
     -   2345R(2,1)-Dij2345R(3,2)-Dij2345R(4,2)+Dij2345R(5,2)+Dij2345
     -   R(6,2)+s15*(EE0R+EijR(1,1))-s12*EijR(9,3)+s12*EijR(10,3)-10*
     -   EijR(11,2)+2*(Dij2345R(3,1)+EijR(39,4)+EijR(43,4)-EijR(44,4)
     -   -EijR(45,4))+EijR(4,3)*P(17)+(EijR(2,1)+2*EijR(5,2)+EijR(14,
     -   3)+EijR(17,3)-EijR(18,3)-EijR(19,3))*P(18)+Is45*(Cij134R(1,2
     -   )-Cij134R(3,2)-p2sq*Dij1234R(2,2)+s12*(-Dij1234R(2,1)+Dij123
     -   4R(3,1)+Dij1234R(6,2))+2*Dij1234R(12,3)+Dij1234R(4,2)*P(18))
     -   +EijR(15,3)*P(24)-(-EijR(12,3)+EijR(13,3))*P(25)+Dij1345R(9,
     -   3)*P(27)+EijR(16,3)*P(31)-(D01345R+Dij1345R(1,1))*P(39)-(D01
     -   234R+Dij1234R(1,1))*P(40)-EijR(3,2)*P(43)+Is12*(Cij145R(1,2)
     -   +Cij145R(2,2)-2*(Cij145R(3,2)+Dij1345R(7,2)-Dij1345R(12,3)+D
     -   ij1345R(13,3))+(-Dij1345R(6,3)+Dij1345R(8,3))*P(1)-(Dij1345R
     -   (3,3)+2*Dij1345R(10,3))*P(5)-Dij1345R(7,3)*P(32)-Dij1345R(2,
     -   1)*P(45)-Dij1345R(3,1)*P(46)-Dij1345R(2,2)*P(62)-Dij1345R(3,
     -   2)*P(63)+(-Dij1345R(4,2)+Dij1345R(5,2))*P(64)+Dij1345R(6,2)*
     -   P(67))-EijR(2,2)*P(68)-EijR(3,1)*P(69)-EijR(4,1)*P(70)+EijR(
     -   4,2)*P(71)-EijR(6,2)*P(72)-EijR(7,2)*P(73)+EijR(8,2)*P(74)-E
     -   ijR(9,2)*P(75)-EijR(10,2)*P(76)-EijR(20,3)*P(77))
       FI(6) = 4*Is12s45*(3*B013I+B014I-4*(Cij123I(4,2)+Cij134I(4,2)+C
     -   ij145I(4,2))+2*(p2sq*Cij123I(2,1)-(C0145I+Cij145I(2,1))*P(15
     -   )-(C0123I+Cij123I(1,1))*P(18)-C0134I*P(42)+Cij134I(2,1)*P(61
     -   )))-8*(D02345I+Dij1345I(4,3)-Dij1345I(5,3)-Dij2345I(1,1)-Dij
     -   2345I(2,1)-Dij2345I(3,2)-Dij2345I(4,2)+Dij2345I(5,2)+Dij2345
     -   I(6,2)+s15*(EE0I+EijI(1,1))-s12*EijI(9,3)+s12*EijI(10,3)-10*
     -   EijI(11,2)+2*(Dij2345I(3,1)+EijI(39,4)+EijI(43,4)-EijI(44,4)
     -   -EijI(45,4))+EijI(4,3)*P(17)+(EijI(2,1)+2*EijI(5,2)+EijI(14,
     -   3)+EijI(17,3)-EijI(18,3)-EijI(19,3))*P(18)+Is45*(Cij134I(1,2
     -   )-Cij134I(3,2)-p2sq*Dij1234I(2,2)+s12*(-Dij1234I(2,1)+Dij123
     -   4I(3,1)+Dij1234I(6,2))+2*Dij1234I(12,3)+Dij1234I(4,2)*P(18))
     -   +EijI(15,3)*P(24)-(-EijI(12,3)+EijI(13,3))*P(25)+Dij1345I(9,
     -   3)*P(27)+EijI(16,3)*P(31)-(D01345I+Dij1345I(1,1))*P(39)-(D01
     -   234I+Dij1234I(1,1))*P(40)-EijI(3,2)*P(43)+Is12*(Cij145I(1,2)
     -   +Cij145I(2,2)-2*(Cij145I(3,2)+Dij1345I(7,2)-Dij1345I(12,3)+D
     -   ij1345I(13,3))+(-Dij1345I(6,3)+Dij1345I(8,3))*P(1)-(Dij1345I
     -   (3,3)+2*Dij1345I(10,3))*P(5)-Dij1345I(7,3)*P(32)-Dij1345I(2,
     -   1)*P(45)-Dij1345I(3,1)*P(46)-Dij1345I(2,2)*P(62)-Dij1345I(3,
     -   2)*P(63)+(-Dij1345I(4,2)+Dij1345I(5,2))*P(64)+Dij1345I(6,2)*
     -   P(67))-EijI(2,2)*P(68)-EijI(3,1)*P(69)-EijI(4,1)*P(70)+EijI(
     -   4,2)*P(71)-EijI(6,2)*P(72)-EijI(7,2)*P(73)+EijI(8,2)*P(74)-E
     -   ijI(9,2)*P(75)-EijI(10,2)*P(76)-EijI(20,3)*P(77))
       F(6)=DCMPLX(FR(6),FI(6))
       P(78) = p2sq+p3sq-s12-s23+s45
       P(79) = p4sq+2*P(9)
       P(80) = -s12+s45+2*P(1)
       P(81) = -2*p2sq+s12
       P(82) = p4sq-3*s12+2*P(1)
       P(83) = p2sq+p3sq-p4sq-s23+s34-2*P(52)
       P(84) = -p3sq+p4sq+s12+s15+s23-s34-2*s45
       P(85) = p2sq-p3sq-p4sq-s12-s15+s23+s34-s45
       P(86) = 2*p2sq-s12
       P(87) = p2sq-p4sq-s23+s34+2*P(43)
       P(88) = p2sq-2*p3sq-p4sq-s12+s15+s23+s34
       P(89) = p2sq+p3sq-p4sq-s23+s34-3*P(52)
       P(90) = 2*p2sq-s12-s15+s34
       FR(7) = 4*Is12s45*(3*B013R+B014R-2*(C0134R*p3sq-p2sq*Cij123R(2,
     -   1)+p4sq*Cij145R(2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,
     -   2))+C0145R*P(15)+(C0123R+Cij123R(1,1))*P(18)+Cij134R(2,1)*P(
     -   43)))-8*(Dij1345R(1,2)+Dij1345R(1,3)-Dij2345R(1,1)+Dij2345R(
     -   3,1)+s15*(EE0R+2*EijR(1,1)+EijR(1,2))-s12*EijR(2,2)+p2sq*(Ei
     -   jR(8,3)-EijR(10,3))-2*(2*EijR(11,2)-EijR(22,3)+EijR(24,3)-Ei
     -   jR(39,4)-EijR(40,4)+EijR(42,4)+EijR(44,4))-Dij1345R(5,3)*P(3
     -   )+Dij1345R(7,3)*P(4)+(-EijR(4,3)+EijR(14,3))*P(8)+EijR(2,1)*
     -   P(18)-EijR(9,2)*P(23)+EijR(8,2)*P(25)-(D01345R+Dij1345R(1,1)
     -   )*P(39)+Cij134R(3,2)*P(44)+Is45*(Cij134R(1,2)+s23*(D01234R+D
     -   ij1234R(1,1))-2*(p2sq*Dij1234R(2,2)+Dij1234R(7,2)-Dij1234R(1
     -   1,3))+Dij1234R(6,2)*P(1)+Dij1234R(2,1)*P(18)-Dij1234R(3,1)*P
     -   (45)-Dij1234R(5,2)*P(78)-Dij1234R(4,2)*P(81))-Is12*(Cij145R(
     -   1,1)-2*(Dij1345R(11,3)-Dij1345R(13,3))+(Dij1345R(4,3)+Dij134
     -   5R(9,3)-2*Dij1345R(10,3))*P(1)+Dij1345R(3,3)*P(5)+(Dij1345R(
     -   2,1)-Dij1345R(3,1))*P(78)+Dij1345R(3,2)*P(79)+(Dij1345R(4,2)
     -   -Dij1345R(6,2))*P(80)-Dij1345R(5,2)*P(82))-EijR(3,1)*P(83)-E
     -   ijR(4,1)*P(84)+EijR(4,2)*P(85)+EijR(5,2)*P(86)-EijR(6,2)*P(8
     -   7)-EijR(7,2)*P(88)+EijR(10,2)*P(89)+(EijR(15,3)-EijR(18,3))*
     -   P(90))
       FI(7) = 4*Is12s45*(3*B013I+B014I-2*(C0134I*p3sq-p2sq*Cij123I(2,
     -   1)+p4sq*Cij145I(2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,
     -   2))+C0145I*P(15)+(C0123I+Cij123I(1,1))*P(18)+Cij134I(2,1)*P(
     -   43)))-8*(Dij1345I(1,2)+Dij1345I(1,3)-Dij2345I(1,1)+Dij2345I(
     -   3,1)+s15*(EE0I+2*EijI(1,1)+EijI(1,2))-s12*EijI(2,2)+p2sq*(Ei
     -   jI(8,3)-EijI(10,3))-2*(2*EijI(11,2)-EijI(22,3)+EijI(24,3)-Ei
     -   jI(39,4)-EijI(40,4)+EijI(42,4)+EijI(44,4))-Dij1345I(5,3)*P(3
     -   )+Dij1345I(7,3)*P(4)+(-EijI(4,3)+EijI(14,3))*P(8)+EijI(2,1)*
     -   P(18)-EijI(9,2)*P(23)+EijI(8,2)*P(25)-(D01345I+Dij1345I(1,1)
     -   )*P(39)+Cij134I(3,2)*P(44)+Is45*(Cij134I(1,2)+s23*(D01234I+D
     -   ij1234I(1,1))-2*(p2sq*Dij1234I(2,2)+Dij1234I(7,2)-Dij1234I(1
     -   1,3))+Dij1234I(6,2)*P(1)+Dij1234I(2,1)*P(18)-Dij1234I(3,1)*P
     -   (45)-Dij1234I(5,2)*P(78)-Dij1234I(4,2)*P(81))-Is12*(Cij145I(
     -   1,1)-2*(Dij1345I(11,3)-Dij1345I(13,3))+(Dij1345I(4,3)+Dij134
     -   5I(9,3)-2*Dij1345I(10,3))*P(1)+Dij1345I(3,3)*P(5)+(Dij1345I(
     -   2,1)-Dij1345I(3,1))*P(78)+Dij1345I(3,2)*P(79)+(Dij1345I(4,2)
     -   -Dij1345I(6,2))*P(80)-Dij1345I(5,2)*P(82))-EijI(3,1)*P(83)-E
     -   ijI(4,1)*P(84)+EijI(4,2)*P(85)+EijI(5,2)*P(86)-EijI(6,2)*P(8
     -   7)-EijI(7,2)*P(88)+EijI(10,2)*P(89)+(EijI(15,3)-EijI(18,3))*
     -   P(90))
       F(7)=DCMPLX(FR(7),FI(7))
       P(91) = 1-Is45*p2sq
       P(92) = -1+Is12*s34
       P(93) = p2sq-p4sq-s12-s23+s34+s45
       P(94) = p2sq+p3sq-p4sq-2*s12-s23
       P(95) = p3sq-s34-s45
       P(96) = -p3sq+2*s12+s34+s45
       P(97) = p2sq-p4sq-s12-s23
       P(98) = p3sq+s34+s45+2*P(97)
       P(99) = p4sq+s12+s15-s34-s45
       P(100) = -p2sq+p4sq+s23-s34+2*P(25)
       P(101) = p2sq+p3sq+p4sq-s34
       P(102) = s15+2*P(45)
       P(103) = p2sq-s12+s15-2*P(36)
       P(104) = p3sq+s12-s23
       P(105) = p2sq-p3sq+s23
       P(106) = p4sq-s34-s45
       P(107) = 3*p3sq+s12+s23+2*P(106)
       FR(8) = -8*(-D02345R+Dij1345R(1,2)-Dij1345R(3,1)+Dij1345R(4,3)-
     -   Dij1345R(5,3)-Dij2345R(3,2)-Dij2345R(4,2)+Dij2345R(5,2)+Dij2
     -   345R(6,2)+p2sq*EijR(2,2)-s12*EijR(9,3)+s12*EijR(10,3)+8*EijR
     -   (11,2)+2*(EijR(22,3)-EijR(23,3)+EijR(39,4)+EijR(43,4)-EijR(4
     -   4,4)-EijR(45,4))+EijR(4,3)*P(17)-(EijR(5,2)-EijR(14,3)-EijR(
     -   17,3)+EijR(18,3)+EijR(19,3))*P(18)+Is45*(Cij134R(1,1)+Cij134
     -   R(3,2)-p2sq*(Dij1234R(2,1)+Dij1234R(6,2))+2*(Dij1234R(7,2)+D
     -   ij1234R(13,3))+Dij1234R(5,2)*P(18))+EijR(15,3)*P(24)-(-EijR(
     -   12,3)+EijR(13,3))*P(25)+Dij1345R(9,3)*P(27)+EijR(16,3)*P(31)
     -   -Dij1234R(3,2)*P(41)+(Cij134R(2,1)+Cij134R(2,2))*P(44)-EijR(
     -   20,3)*P(77)-Dij1234R(3,1)*P(91)-Dij1345R(1,1)*P(92)-Is12*(C0
     -   134R-s34*Dij1345R(2,1)-4*Dij1345R(7,2)-2*(Dij1345R(12,3)-Dij
     -   1345R(13,3))+(Dij1345R(6,3)-Dij1345R(8,3))*P(1)+(Dij1345R(3,
     -   3)+2*Dij1345R(10,3))*P(5)+Dij1345R(7,3)*P(32)+Dij1345R(2,2)*
     -   P(93)+Dij1345R(3,2)*P(94)+Dij1345R(4,2)*P(95)+Dij1345R(5,2)*
     -   P(96)-Dij1345R(6,2)*P(98))-(-EijR(3,1)+EijR(4,1))*P(99)+EijR
     -   (3,2)*P(100)+EijR(4,2)*P(101)+EijR(6,2)*P(102)-EijR(7,2)*P(1
     -   03)-EijR(8,2)*P(104)-EijR(9,2)*P(105)-EijR(10,2)*P(107))
       FI(8) = -8*(-D02345I+Dij1345I(1,2)-Dij1345I(3,1)+Dij1345I(4,3)-
     -   Dij1345I(5,3)-Dij2345I(3,2)-Dij2345I(4,2)+Dij2345I(5,2)+Dij2
     -   345I(6,2)+p2sq*EijI(2,2)-s12*EijI(9,3)+s12*EijI(10,3)+8*EijI
     -   (11,2)+2*(EijI(22,3)-EijI(23,3)+EijI(39,4)+EijI(43,4)-EijI(4
     -   4,4)-EijI(45,4))+EijI(4,3)*P(17)-(EijI(5,2)-EijI(14,3)-EijI(
     -   17,3)+EijI(18,3)+EijI(19,3))*P(18)+Is45*(Cij134I(1,1)+Cij134
     -   I(3,2)-p2sq*(Dij1234I(2,1)+Dij1234I(6,2))+2*(Dij1234I(7,2)+D
     -   ij1234I(13,3))+Dij1234I(5,2)*P(18))+EijI(15,3)*P(24)-(-EijI(
     -   12,3)+EijI(13,3))*P(25)+Dij1345I(9,3)*P(27)+EijI(16,3)*P(31)
     -   -Dij1234I(3,2)*P(41)+(Cij134I(2,1)+Cij134I(2,2))*P(44)-EijI(
     -   20,3)*P(77)-Dij1234I(3,1)*P(91)-Dij1345I(1,1)*P(92)-Is12*(C0
     -   134I-s34*Dij1345I(2,1)-4*Dij1345I(7,2)-2*(Dij1345I(12,3)-Dij
     -   1345I(13,3))+(Dij1345I(6,3)-Dij1345I(8,3))*P(1)+(Dij1345I(3,
     -   3)+2*Dij1345I(10,3))*P(5)+Dij1345I(7,3)*P(32)+Dij1345I(2,2)*
     -   P(93)+Dij1345I(3,2)*P(94)+Dij1345I(4,2)*P(95)+Dij1345I(5,2)*
     -   P(96)-Dij1345I(6,2)*P(98))-(-EijI(3,1)+EijI(4,1))*P(99)+EijI
     -   (3,2)*P(100)+EijI(4,2)*P(101)+EijI(6,2)*P(102)-EijI(7,2)*P(1
     -   03)-EijI(8,2)*P(104)-EijI(9,2)*P(105)-EijI(10,2)*P(107))
       F(8)=DCMPLX(FR(8),FI(8))
       P(108) = 1-3*Is12*P(1)
       P(109) = 2-3*Is12*P(1)
       P(110) = 3*s12-2*P(1)
       P(111) = p3sq+2*s12+s15-s34-s45
       P(112) = p3sq+s12-s34-s45
       P(113) = p2sq+s15+2*P(112)
       P(114) = p2sq+s12-s15+3*P(30)
       P(115) = p2sq-s15
       P(116) = -s12+3*P(30)+2*P(115)
       FR(9) = -8*(Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R(4,2)-Dij1345R(
     -   5,2)+Dij1345R(6,3)+Dij1345R(7,3)-Dij2345R(2,1)-Dij2345R(2,2)
     -   +Dij2345R(3,1)-Dij2345R(3,2)-s12*(EijR(8,2)-EijR(9,2)+EijR(1
     -   2,3)+EijR(15,3)-2*EijR(20,3))-2*(Dij1345R(10,3)-Dij2345R(6,2
     -   )-EijR(23,3)+EijR(24,3)-EijR(38,4)-EijR(39,4)+2*EijR(45,4))+
     -   EijR(4,3)*P(17)+(EijR(6,2)-EijR(7,2)+EijR(11,3)+EijR(14,3)-2
     -   *EijR(19,3))*P(18)+Is45*(Cij134R(1,1)-Cij134R(2,1)-Cij134R(2
     -   ,2)+Cij134R(3,2)+s12*Dij1234R(3,2)-p2sq*(Dij1234R(2,1)-Dij12
     -   34R(3,1)+Dij1234R(6,2))+2*(Dij1234R(7,2)+Dij1234R(13,3))+Dij
     -   1234R(5,2)*P(18))+EijR(3,3)*P(25)-(-EijR(3,1)+EijR(4,1))*P(2
     -   8)+EijR(4,2)*P(33)-Dij1345R(8,3)*P(108)+Dij1345R(9,3)*P(109)
     -   -Is12*(Dij1345R(2,3)*P(1)+Dij1345R(3,3)*P(5)+Dij1345R(2,2)*P
     -   (9)+Dij1345R(3,2)*P(32)+Dij1345R(6,2)*P(110))+EijR(3,2)*P(11
     -   1)-EijR(10,2)*P(113)-EijR(13,3)*P(114)+EijR(16,3)*P(116))
       FI(9) = -8*(Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I(4,2)-Dij1345I(
     -   5,2)+Dij1345I(6,3)+Dij1345I(7,3)-Dij2345I(2,1)-Dij2345I(2,2)
     -   +Dij2345I(3,1)-Dij2345I(3,2)-s12*(EijI(8,2)-EijI(9,2)+EijI(1
     -   2,3)+EijI(15,3)-2*EijI(20,3))-2*(Dij1345I(10,3)-Dij2345I(6,2
     -   )-EijI(23,3)+EijI(24,3)-EijI(38,4)-EijI(39,4)+2*EijI(45,4))+
     -   EijI(4,3)*P(17)+(EijI(6,2)-EijI(7,2)+EijI(11,3)+EijI(14,3)-2
     -   *EijI(19,3))*P(18)+Is45*(Cij134I(1,1)-Cij134I(2,1)-Cij134I(2
     -   ,2)+Cij134I(3,2)+s12*Dij1234I(3,2)-p2sq*(Dij1234I(2,1)-Dij12
     -   34I(3,1)+Dij1234I(6,2))+2*(Dij1234I(7,2)+Dij1234I(13,3))+Dij
     -   1234I(5,2)*P(18))+EijI(3,3)*P(25)-(-EijI(3,1)+EijI(4,1))*P(2
     -   8)+EijI(4,2)*P(33)-Dij1345I(8,3)*P(108)+Dij1345I(9,3)*P(109)
     -   -Is12*(Dij1345I(2,3)*P(1)+Dij1345I(3,3)*P(5)+Dij1345I(2,2)*P
     -   (9)+Dij1345I(3,2)*P(32)+Dij1345I(6,2)*P(110))+EijI(3,2)*P(11
     -   1)-EijI(10,2)*P(113)-EijI(13,3)*P(114)+EijI(16,3)*P(116))
       F(9)=DCMPLX(FR(9),FI(9))
       P(117) = 2*s12-3*P(1)
       P(118) = p2sq+p3sq-s12-s15-s23+s34
       P(119) = p3sq-s12-s15-s23+s34
       P(120) = 2*p3sq-s15-s23+s34-s45
       P(121) = p2sq-2*p3sq-s12+s23+s34+s45
       P(122) = s12+s23-s34-s45
       FR(10) = -8*(-D02345R+Dij1345R(4,3)-Dij1345R(5,3)-Dij2345R(2,1)
     -   -s12*EijR(8,2)-p2sq*(-EijR(2,1)+EijR(9,2)-EijR(15,3)-EijR(17
     -   ,3)+EijR(18,3)+EijR(20,3))+2*(2*EijR(11,2)+EijR(23,3)-EijR(2
     -   4,3)+EijR(39,4)+EijR(41,4)-EijR(42,4)-EijR(45,4))+Is45*(Cij1
     -   34R(1,2)-Cij134R(3,2)+p2sq*Dij1234R(4,2)+4*Dij1234R(7,2)+2*(
     -   Cij134R(1,1)-Cij134R(2,1)-p2sq*Dij1234R(6,2)+Dij1234R(11,3))
     -   +Dij1234R(3,2)*P(1)-Dij1234R(5,2)*P(6))+EijR(4,2)*P(7)-(EijR
     -   (4,3)-EijR(14,3)-EijR(16,3)+EijR(19,3))*P(8)+EijR(3,2)*P(25)
     -   +Dij1345R(9,3)*P(27)+EijR(5,2)*P(68)-Is12*((Dij1345R(2,2)+Di
     -   j1345R(6,3)-Dij1345R(8,3))*P(1)+(Dij1345R(3,1)+Dij1345R(3,3)
     -   +2*Dij1345R(10,3))*P(5)+(Dij1345R(2,1)+2*Dij1345R(3,2))*P(9)
     -   +(Dij1345R(4,2)-Dij1345R(5,2)+Dij1345R(7,3))*P(32)+Dij1345R(
     -   6,2)*P(117))-EijR(3,1)*P(118)+EijR(4,1)*P(119)-EijR(6,2)*P(1
     -   20)-EijR(7,2)*P(121)-EijR(10,2)*P(122))
       FI(10) = -8*(-D02345I+Dij1345I(4,3)-Dij1345I(5,3)-Dij2345I(2,1)
     -   -s12*EijI(8,2)-p2sq*(-EijI(2,1)+EijI(9,2)-EijI(15,3)-EijI(17
     -   ,3)+EijI(18,3)+EijI(20,3))+2*(2*EijI(11,2)+EijI(23,3)-EijI(2
     -   4,3)+EijI(39,4)+EijI(41,4)-EijI(42,4)-EijI(45,4))+Is45*(Cij1
     -   34I(1,2)-Cij134I(3,2)+p2sq*Dij1234I(4,2)+4*Dij1234I(7,2)+2*(
     -   Cij134I(1,1)-Cij134I(2,1)-p2sq*Dij1234I(6,2)+Dij1234I(11,3))
     -   +Dij1234I(3,2)*P(1)-Dij1234I(5,2)*P(6))+EijI(4,2)*P(7)-(EijI
     -   (4,3)-EijI(14,3)-EijI(16,3)+EijI(19,3))*P(8)+EijI(3,2)*P(25)
     -   +Dij1345I(9,3)*P(27)+EijI(5,2)*P(68)-Is12*((Dij1345I(2,2)+Di
     -   j1345I(6,3)-Dij1345I(8,3))*P(1)+(Dij1345I(3,1)+Dij1345I(3,3)
     -   +2*Dij1345I(10,3))*P(5)+(Dij1345I(2,1)+2*Dij1345I(3,2))*P(9)
     -   +(Dij1345I(4,2)-Dij1345I(5,2)+Dij1345I(7,3))*P(32)+Dij1345I(
     -   6,2)*P(117))-EijI(3,1)*P(118)+EijI(4,1)*P(119)-EijI(6,2)*P(1
     -   20)-EijI(7,2)*P(121)-EijI(10,2)*P(122))
       F(10)=DCMPLX(FR(10),FI(10))
       P(123) = p2sq+p3sq+s12-s23
       P(124) = 2*p2sq+p3sq-s12-s23
       P(125) = -s12+2*P(1)
       FR(11) = 8*(-Dij1345R(3,2)+Dij1345R(5,3)+Dij2345R(3,1)+Dij2345R
     -   (3,2)+p2sq*EijR(9,2)+s12*EijR(15,3)-s12*EijR(18,3)-2*(EijR(3
     -   9,4)-EijR(42,4))-EijR(10,2)*P(1)-EijR(4,3)*P(17)+EijR(7,3)*P
     -   (18)-EijR(14,3)*P(23)+(-EijR(16,3)+EijR(19,3))*P(25)+Is12*(-
     -   ((Dij1345R(6,2)-Dij1345R(9,3)+Dij1345R(10,3))*P(1))+Dij1345R
     -   (3,3)*P(5)+Dij1345R(7,3)*P(32)+Dij1345R(5,2)*P(123))-EijR(4,
     -   2)*P(124)+EijR(7,2)*P(125))
       FI(11) = 8*(-Dij1345I(3,2)+Dij1345I(5,3)+Dij2345I(3,1)+Dij2345I
     -   (3,2)+p2sq*EijI(9,2)+s12*EijI(15,3)-s12*EijI(18,3)-2*(EijI(3
     -   9,4)-EijI(42,4))-EijI(10,2)*P(1)-EijI(4,3)*P(17)+EijI(7,3)*P
     -   (18)-EijI(14,3)*P(23)+(-EijI(16,3)+EijI(19,3))*P(25)+Is12*(-
     -   ((Dij1345I(6,2)-Dij1345I(9,3)+Dij1345I(10,3))*P(1))+Dij1345I
     -   (3,3)*P(5)+Dij1345I(7,3)*P(32)+Dij1345I(5,2)*P(123))-EijI(4,
     -   2)*P(124)+EijI(7,2)*P(125))
       F(11)=DCMPLX(FR(11),FI(11))
       P(126) = p2sq+p3sq+p4sq-2*s12-s23+s34
       P(127) = p2sq+p4sq-2*s12-s23+s34+s45
       P(128) = p4sq-s34
       P(129) = p2sq+p3sq-s12+s15-s23+s45
       P(130) = s12+s23-s34-3*s45
       P(131) = -p3sq+s12-s15+s23+s34
       P(132) = p2sq+s12-s23-s45
       FR(12) = -8*(D02345R-Dij1345R(1,2)-Dij1345R(5,3)+Dij2345R(3,1)-
     -   Dij2345R(3,2)+Dij2345R(5,2)-p2sq*EijR(2,2)+s12*EijR(10,3)+Ei
     -   jR(8,2)*P(1)+EijR(4,3)*P(17)+(EijR(5,2)+EijR(14,3)-EijR(18,3
     -   ))*P(18)+EijR(15,3)*P(24)+(EijR(16,3)-EijR(20,3))*P(25)-2*(2
     -   *EijR(11,2)+EijR(22,3)-EijR(24,3)-EijR(39,4)+EijR(44,4)+EijR
     -   (7,2)*P(36))-Dij1345R(1,1)*P(39)+EijR(3,1)*P(78)+Is12*(C0134
     -   R-C0145R-Cij145R(1,1)+Cij145R(2,2)-Cij145R(3,2)-4*Dij1345R(7
     -   ,2)-2*Dij1345R(13,3)+(Dij1345R(4,2)-Dij1345R(9,3)+Dij1345R(1
     -   0,3))*P(1)-Dij1345R(3,3)*P(5)-Dij1345R(7,3)*P(32)+Dij1345R(2
     -   ,1)*P(62)-Dij1345R(5,2)*P(95)-Dij1345R(3,1)*P(126)-Dij1345R(
     -   3,2)*P(127)+Dij1345R(6,2)*P(128))-EijR(4,1)*P(129)+EijR(4,2)
     -   *P(130)+EijR(9,2)*P(131)-EijR(10,2)*P(132))
       FI(12) = -8*(D02345I-Dij1345I(1,2)-Dij1345I(5,3)+Dij2345I(3,1)-
     -   Dij2345I(3,2)+Dij2345I(5,2)-p2sq*EijI(2,2)+s12*EijI(10,3)+Ei
     -   jI(8,2)*P(1)+EijI(4,3)*P(17)+(EijI(5,2)+EijI(14,3)-EijI(18,3
     -   ))*P(18)+EijI(15,3)*P(24)+(EijI(16,3)-EijI(20,3))*P(25)-2*(2
     -   *EijI(11,2)+EijI(22,3)-EijI(24,3)-EijI(39,4)+EijI(44,4)+EijI
     -   (7,2)*P(36))-Dij1345I(1,1)*P(39)+EijI(3,1)*P(78)+Is12*(C0134
     -   I-C0145I-Cij145I(1,1)+Cij145I(2,2)-Cij145I(3,2)-4*Dij1345I(7
     -   ,2)-2*Dij1345I(13,3)+(Dij1345I(4,2)-Dij1345I(9,3)+Dij1345I(1
     -   0,3))*P(1)-Dij1345I(3,3)*P(5)-Dij1345I(7,3)*P(32)+Dij1345I(2
     -   ,1)*P(62)-Dij1345I(5,2)*P(95)-Dij1345I(3,1)*P(126)-Dij1345I(
     -   3,2)*P(127)+Dij1345I(6,2)*P(128))-EijI(4,1)*P(129)+EijI(4,2)
     -   *P(130)+EijI(9,2)*P(131)-EijI(10,2)*P(132))
       F(12)=DCMPLX(FR(12),FI(12))
       P(133) = 1-Is12*P(1)
       P(134) = p2sq+s15-s23+s45
       P(135) = -2*p3sq+s15+s23+s45
       FR(13) = 8*(D01345R+Dij1345R(1,2)+Dij1345R(5,3)-s12*(EijR(2,1)+
     -   EijR(5,2))+p2sq*(-EijR(15,3)+EijR(18,3))+2*(Dij1345R(1,1)+Ei
     -   jR(21,3)-EijR(39,4)+EijR(42,4))+EijR(10,2)*P(1)+EijR(4,2)*P(
     -   6)+(EijR(4,3)-EijR(14,3))*P(8)+Is12*(-((Dij1345R(2,1)+Dij134
     -   5R(4,2)-Dij1345R(9,3)+Dij1345R(10,3))*P(1))+Dij1345R(3,3)*P(
     -   5)+Dij1345R(3,2)*P(9)+Dij1345R(7,3)*P(32))-(EijR(3,1)+EijR(6
     -   ,2))*P(45)-Dij1345R(3,1)*P(133)+EijR(4,1)*P(134)+EijR(7,2)*P
     -   (135))
       FI(13) = 8*(D01345I+Dij1345I(1,2)+Dij1345I(5,3)-s12*(EijI(2,1)+
     -   EijI(5,2))+p2sq*(-EijI(15,3)+EijI(18,3))+2*(Dij1345I(1,1)+Ei
     -   jI(21,3)-EijI(39,4)+EijI(42,4))+EijI(10,2)*P(1)+EijI(4,2)*P(
     -   6)+(EijI(4,3)-EijI(14,3))*P(8)+Is12*(-((Dij1345I(2,1)+Dij134
     -   5I(4,2)-Dij1345I(9,3)+Dij1345I(10,3))*P(1))+Dij1345I(3,3)*P(
     -   5)+Dij1345I(3,2)*P(9)+Dij1345I(7,3)*P(32))-(EijI(3,1)+EijI(6
     -   ,2))*P(45)-Dij1345I(3,1)*P(133)+EijI(4,1)*P(134)+EijI(7,2)*P
     -   (135))
       F(13)=DCMPLX(FR(13),FI(13))
       P(136) = s12+s34
       P(137) = p4sq+s12
       P(138) = p2sq-s23+s34+s45
       P(139) = -p2sq+p4sq+s12+s23-s34-s45
       P(140) = -p3sq+s12+s15
       P(141) = 2*p2sq+p3sq+p4sq-s12-s23+s45
       P(142) = 2*p2sq-p3sq+s34+s45
       P(143) = p2sq-p3sq-p4sq-s12-s23+s34+s45
       FR(14) = 8*(D01345R+Dij1345R(5,3)-Dij2345R(1,1)+Dij2345R(3,1)+D
     -   ij2345R(3,2)-Dij2345R(5,2)+p3sq*EijR(3,1)-s12*(EijR(2,1)+Eij
     -   R(10,3))-2*(2*EijR(11,2)+EijR(24,3)+EijR(39,4)-EijR(44,4))-E
     -   ijR(4,3)*P(17)+(EijR(5,2)-EijR(14,3)+EijR(18,3))*P(18)-EijR(
     -   15,3)*P(24)+(EijR(8,2)-EijR(16,3)+EijR(20,3))*P(25)-EijR(2,2
     -   )*P(68)+EijR(7,2)*P(103)+Is12*(C0134R-p4sq*Dij1345R(3,1)-2*(
     -   Dij1345R(7,2)-Dij1345R(13,3))+(Dij1345R(9,3)-Dij1345R(10,3))
     -   *P(1)+Dij1345R(3,3)*P(5)+Dij1345R(7,3)*P(32)+Dij1345R(2,1)*P
     -   (128)+Dij1345R(1,1)*P(136)-Dij1345R(3,2)*P(137)+Dij1345R(5,2
     -   )*P(138)+Dij1345R(6,2)*P(139))+EijR(4,1)*P(140)-EijR(4,2)*P(
     -   141)+EijR(9,2)*P(142)-EijR(10,2)*P(143))
       FI(14) = 8*(D01345I+Dij1345I(5,3)-Dij2345I(1,1)+Dij2345I(3,1)+D
     -   ij2345I(3,2)-Dij2345I(5,2)+p3sq*EijI(3,1)-s12*(EijI(2,1)+Eij
     -   I(10,3))-2*(2*EijI(11,2)+EijI(24,3)+EijI(39,4)-EijI(44,4))-E
     -   ijI(4,3)*P(17)+(EijI(5,2)-EijI(14,3)+EijI(18,3))*P(18)-EijI(
     -   15,3)*P(24)+(EijI(8,2)-EijI(16,3)+EijI(20,3))*P(25)-EijI(2,2
     -   )*P(68)+EijI(7,2)*P(103)+Is12*(C0134I-p4sq*Dij1345I(3,1)-2*(
     -   Dij1345I(7,2)-Dij1345I(13,3))+(Dij1345I(9,3)-Dij1345I(10,3))
     -   *P(1)+Dij1345I(3,3)*P(5)+Dij1345I(7,3)*P(32)+Dij1345I(2,1)*P
     -   (128)+Dij1345I(1,1)*P(136)-Dij1345I(3,2)*P(137)+Dij1345I(5,2
     -   )*P(138)+Dij1345I(6,2)*P(139))+EijI(4,1)*P(140)-EijI(4,2)*P(
     -   141)+EijI(9,2)*P(142)-EijI(10,2)*P(143))
       F(14)=DCMPLX(FR(14),FI(14))
       P(144) = -p3sq+s12+s15+s45
       P(145) = p2sq+s12+s15-s23+s45
       FR(15) = -8*(-D01345R+D02345R-Dij1345R(1,1)+Dij1345R(3,2)-Dij13
     -   45R(6,2)+Dij1345R(7,3)-Dij1345R(10,3)+Dij2345R(2,1)-Dij2345R
     -   (3,2)+Dij2345R(6,2)+s12*(EijR(2,1)-EijR(15,3)+EijR(20,3))-2*
     -   (EijR(11,2)-EijR(39,4)+EijR(45,4))+Is12*(Dij1345R(8,3)*P(1)-
     -   Dij1345R(3,3)*P(5))+EijR(4,3)*P(17)-(EijR(6,2)+EijR(7,2)-Eij
     -   R(14,3)+EijR(19,3))*P(18)-(EijR(3,2)+EijR(13,3))*P(25)+Dij13
     -   45R(9,3)*P(27)-EijR(3,1)*P(30)+EijR(16,3)*P(31)+EijR(8,2)*P(
     -   68)+EijR(4,2)*P(124)-EijR(4,1)*P(144)-EijR(10,2)*P(145))
       FI(15) = -8*(-D01345I+D02345I-Dij1345I(1,1)+Dij1345I(3,2)-Dij13
     -   45I(6,2)+Dij1345I(7,3)-Dij1345I(10,3)+Dij2345I(2,1)-Dij2345I
     -   (3,2)+Dij2345I(6,2)+s12*(EijI(2,1)-EijI(15,3)+EijI(20,3))-2*
     -   (EijI(11,2)-EijI(39,4)+EijI(45,4))+Is12*(Dij1345I(8,3)*P(1)-
     -   Dij1345I(3,3)*P(5))+EijI(4,3)*P(17)-(EijI(6,2)+EijI(7,2)-Eij
     -   I(14,3)+EijI(19,3))*P(18)-(EijI(3,2)+EijI(13,3))*P(25)+Dij13
     -   45I(9,3)*P(27)-EijI(3,1)*P(30)+EijI(16,3)*P(31)+EijI(8,2)*P(
     -   68)+EijI(4,2)*P(124)-EijI(4,1)*P(144)-EijI(10,2)*P(145))
       F(15)=DCMPLX(FR(15),FI(15))
       P(146) = s45+2*P(9)
       P(147) = p3sq+s12
       P(148) = s12*s45+P(9)*P(147)
       P(149) = p2sq+p3sq-p4sq-s12-s23+s34+s45-Is12*s34*P(1)
       P(150) = -p3sq+s23
       P(151) = s23*P(9)+s45*P(150)
       P(152) = p3sq-p4sq+s34+2*P(45)
       P(153) = -s12+s34
       P(154) = p4sq*s12+s12**2+s12*s23-s12*s34-s23*s34-s12*s45+p2sq*P
     -   (153)+p3sq*P(153)
       P(155) = p3sq-s15-s23+s34+2*P(18)
       P(156) = p2sq*s45+P(9)*P(18)
       P(157) = 3*p2sq+p3sq-s23
       P(158) = 2*p2sq*s45+P(9)*P(157)
       P(159) = p3sq-s23+2*P(18)
       P(160) = -p2sq+s12
       P(161) = -p3sq-p4sq+s12+s23+s45
       P(162) = p2sq*P(144)-s12*P(161)
       P(163) = s12-s15+s23
       P(164) = p4sq+s12-s45
       P(165) = p4sq+s15-s34
       P(166) = -s45**2+p2sq*P(42)+s12*P(163)-p3sq*P(164)+s45*P(165)
       P(167) = -p4sq-3*s23+s34-2*P(52)
       P(168) = p4sq+s23-s34+2*P(52)
       P(169) = 3*p3sq+s12-s45
       P(170) = 2*p3sq**2+p3sq*P(167)-P(52)*P(168)+p2sq*P(169)
       P(171) = s12+s15+s45
       P(172) = -p2sq+s12+s34+s45
       P(173) = p3sq*P(15)-p4sq*P(171)+s45*P(172)
       P(174) = p4sq+s45
       P(175) = p3sq-s15
       P(176) = 5*s12+s23-s34+s45-2*P(175)
       P(177) = s12+s23
       P(178) = 2*p3sq+s34+s45-3*P(177)
       P(179) = -2*p2sq**2+p2sq*P(176)+s12*P(178)
       P(180) = 2*s12+s15
       P(181) = 2*p3sq-p4sq-5*s12-3*s23+s34+4*s45
       P(182) = p2sq**2-s23*s45
       P(183) = p4sq*s12-s12*s15+s23**2-s12*s34-4*s12*s45+s15*s45+s45*
     -   *2+3*s12*P(177)-p3sq*P(180)+p2sq*P(181)+2*P(182)
       P(184) = p4sq-2*s23+3*s45
       P(185) = s12*s23-s23*s34-s23*s45+s34*s45
       P(186) = -(p4sq*s12)-p3sq*s15+s15*s23+s23**2+s15*s34-3*s12*s45+
     -   s45**2+p2sq*P(184)+2*P(185)
       P(187) = 4*p3sq-p4sq+s12-3*s23+s34
       P(188) = p2sq**2+p3sq**2-s12*s23+s12*s34
       P(189) = 4*s12+s15+s23-s34
       P(190) = s12+s15+3*s23-s34+s45
       P(191) = p4sq*s12+4*s12**2+s12*s15-s23**2-p2sq*P(187)-2*P(188)-
     -   s45*P(189)+p3sq*P(190)
       P(192) = 2*p3sq-s12-s23-s34-s45
       P(193) = -5*s12-s15+s34-3*P(36)
       P(194) = s15+2*s34
       P(195) = -s23+s34+s45
       P(196) = 4*p3sq+p4sq+s45-2*P(136)
       P(197) = p3sq**2+s12**2+s23*s34-s45**2
       P(198) = s12*s15+s12*s23+s12*s34-s34**2+p3sq*P(193)-s45*P(194)+
     -   p4sq*P(195)+p2sq*P(196)+2*P(197)
       P(199) = s45+3*P(9)
       P(200) = s34+s45
       P(201) = 8*p2sq+5*p3sq-p4sq-s15-7*P(177)+3*P(200)
       FR(16) = Is12s45*(2*P(9)*(3*B013R+B014R-2*(C0134R*p3sq-p2sq*Cij
     -   123R(2,1)+(C0123R+Cij123R(1,1))*P(18)))+4*(P(5)*(p4sq*Cij145
     -   R(2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2))+C0145R*P(1
     -   5))-Cij134R(2,1)*P(148)))+Is12*(8*(Dij1345R(12,3)*P(1)+Dij13
     -   45R(13,3)*P(5))+4*((Cij145R(1,1)+2*Dij1345R(7,2))*P(9)-Dij13
     -   45R(1,1)*P(154)))-4*(4*EijR(46,4)-EE0R*s15*P(5)+Dij2345R(2,1
     -   )*P(53)-D01345R*P(149)+D02345R*P(152)+s15*EijR(1,1)*P(159)-s
     -   15*EijR(1,2)*P(160)-EijR(2,1)*P(162)-EijR(3,1)*P(166)-EijR(3
     -   ,2)*P(170)-EijR(4,1)*P(173)-EijR(5,2)*P(179)-EijR(6,2)*P(183
     -   )+EijR(7,2)*P(186)-EijR(8,2)*P(191)-P(155)*(Dij2345R(1,1)-Di
     -   j2345R(3,1)+EijR(2,2)*P(68)+EijR(4,2)*P(174)+EijR(9,2)*P(192
     -   ))+EijR(10,2)*P(198)+Is45*(P(9)*(-2*p2sq*Dij1234R(2,2)+Dij12
     -   34R(4,2)*P(18)-Dij1234R(3,1)*P(45))-(Dij1234R(3,2)*P(1)+Dij1
     -   234R(5,2)*P(18))*P(78)-Cij134R(1,1)*P(146)+(D01234R+Dij1234R
     -   (1,1))*P(151)+Dij1234R(2,1)*P(156)+Dij1234R(6,2)*P(158)-2*Di
     -   j1234R(7,2)*P(199))+2*(Dij1345R(11,3)+p2sq*EijR(22,3)-EijR(2
     -   4,3)*P(8)-EijR(11,2)*P(201)))
       FI(16) = Is12s45*(2*P(9)*(3*B013I+B014I-2*(C0134I*p3sq-p2sq*Cij
     -   123I(2,1)+(C0123I+Cij123I(1,1))*P(18)))+4*(P(5)*(p4sq*Cij145
     -   I(2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2))+C0145I*P(1
     -   5))-Cij134I(2,1)*P(148)))+Is12*(8*(Dij1345I(12,3)*P(1)+Dij13
     -   45I(13,3)*P(5))+4*((Cij145I(1,1)+2*Dij1345I(7,2))*P(9)-Dij13
     -   45I(1,1)*P(154)))-4*(4*EijI(46,4)-EE0I*s15*P(5)+Dij2345I(2,1
     -   )*P(53)-D01345I*P(149)+D02345I*P(152)+s15*EijI(1,1)*P(159)-s
     -   15*EijI(1,2)*P(160)-EijI(2,1)*P(162)-EijI(3,1)*P(166)-EijI(3
     -   ,2)*P(170)-EijI(4,1)*P(173)-EijI(5,2)*P(179)-EijI(6,2)*P(183
     -   )+EijI(7,2)*P(186)-EijI(8,2)*P(191)-P(155)*(Dij2345I(1,1)-Di
     -   j2345I(3,1)+EijI(2,2)*P(68)+EijI(4,2)*P(174)+EijI(9,2)*P(192
     -   ))+EijI(10,2)*P(198)+Is45*(P(9)*(-2*p2sq*Dij1234I(2,2)+Dij12
     -   34I(4,2)*P(18)-Dij1234I(3,1)*P(45))-(Dij1234I(3,2)*P(1)+Dij1
     -   234I(5,2)*P(18))*P(78)-Cij134I(1,1)*P(146)+(D01234I+Dij1234I
     -   (1,1))*P(151)+Dij1234I(2,1)*P(156)+Dij1234I(6,2)*P(158)-2*Di
     -   j1234I(7,2)*P(199))+2*(Dij1345I(11,3)+p2sq*EijI(22,3)-EijI(2
     -   4,3)*P(8)-EijI(11,2)*P(201)))
       F(16)=DCMPLX(FR(16),FI(16))
       P(202) = p2sq+p3sq-s12-s15-s45
       FR(17) = 8*(-Dij1345R(3,1)-Dij1345R(5,2)-Dij1345R(7,3)+Dij2345R
     -   (3,1)+Dij2345R(3,2)+s12*EijR(9,2)+s12*EijR(15,3)-2*(EijR(24,
     -   3)+EijR(39,4))+Is12*((-Dij1345R(3,2)+Dij1345R(6,2)+Dij1345R(
     -   9,3))*P(1)+Dij1345R(3,3)*P(5))-EijR(14,3)*P(18)-EijR(16,3)*P
     -   (25)+EijR(10,2)*P(45)-EijR(4,2)*P(134)+EijR(4,3)*P(202))
       FI(17) = 8*(-Dij1345I(3,1)-Dij1345I(5,2)-Dij1345I(7,3)+Dij2345I
     -   (3,1)+Dij2345I(3,2)+s12*EijI(9,2)+s12*EijI(15,3)-2*(EijI(24,
     -   3)+EijI(39,4))+Is12*((-Dij1345I(3,2)+Dij1345I(6,2)+Dij1345I(
     -   9,3))*P(1)+Dij1345I(3,3)*P(5))-EijI(14,3)*P(18)-EijI(16,3)*P
     -   (25)+EijI(10,2)*P(45)-EijI(4,2)*P(134)+EijI(4,3)*P(202))
       F(17)=DCMPLX(FR(17),FI(17))
       P(203) = p3sq+p4sq+s12-s34-s45
       P(204) = -p2sq+p3sq+p4sq+s12+s23-s34-s45
       FR(18) = -8*(p3sq*(EijR(3,1)-EijR(4,1)+EijR(6,3)+EijR(16,3)+2*(
     -   EijR(6,2)-EijR(10,2)-EijR(19,3)))+2*(EijR(11,2)+EijR(36,4)+E
     -   ijR(39,4)+2*(EijR(21,3)-EijR(24,3)-EijR(42,4)))-(EijR(5,2)+E
     -   ijR(5,3)-EijR(9,2)+EijR(15,3)-2*EijR(18,3))*P(1)+(EijR(7,3)-
     -   2*EijR(14,3))*P(93)-EijR(4,3)*P(139)+Is12*(Cij145R(1,1)+Cij1
     -   45R(1,2)-Cij145R(2,1)+Cij145R(2,2)+p3sq*(Dij1345R(2,1)-Dij13
     -   45R(3,1)+Dij1345R(4,2)-Dij1345R(6,2))-2*(Cij145R(3,2)-Dij134
     -   5R(7,2)-Dij1345R(11,3)+Dij1345R(13,3))+(Dij1345R(3,2)-Dij134
     -   5R(5,2))*P(203))+(EijR(4,2)-EijR(7,2))*P(204))
       FI(18) = -8*(p3sq*(EijI(3,1)-EijI(4,1)+EijI(6,3)+EijI(16,3)+2*(
     -   EijI(6,2)-EijI(10,2)-EijI(19,3)))+2*(EijI(11,2)+EijI(36,4)+E
     -   ijI(39,4)+2*(EijI(21,3)-EijI(24,3)-EijI(42,4)))-(EijI(5,2)+E
     -   ijI(5,3)-EijI(9,2)+EijI(15,3)-2*EijI(18,3))*P(1)+(EijI(7,3)-
     -   2*EijI(14,3))*P(93)-EijI(4,3)*P(139)+Is12*(Cij145I(1,1)+Cij1
     -   45I(1,2)-Cij145I(2,1)+Cij145I(2,2)+p3sq*(Dij1345I(2,1)-Dij13
     -   45I(3,1)+Dij1345I(4,2)-Dij1345I(6,2))-2*(Cij145I(3,2)-Dij134
     -   5I(7,2)-Dij1345I(11,3)+Dij1345I(13,3))+(Dij1345I(3,2)-Dij134
     -   5I(5,2))*P(203))+(EijI(4,2)-EijI(7,2))*P(204))
       F(18)=DCMPLX(FR(18),FI(18))
       P(205) = p2sq+p3sq+p4sq-s23-s34-s45
       P(206) = p2sq+p3sq-p4sq-s12-s23
       P(207) = p2sq+2*p3sq+p4sq-s23-s34-s45
       P(208) = -p3sq+s12+s23
       P(209) = p3sq-p4sq-s23+2*P(18)
       P(210) = p3sq-s12-s23
       P(211) = p2sq-p4sq+2*P(210)
       P(212) = -p4sq-3*s12+2*P(1)
       P(213) = -s12+s45
       FR(19) = 8*(Dij2345R(1,1)-Dij2345R(3,1)-Dij2345R(3,2)+Dij2345R(
     -   5,2)+p3sq*EijR(6,2)-p3sq*EijR(8,2)-4*(EijR(21,3)+EijR(22,3))
     -   +8*EijR(24,3)-2*(EijR(39,4)+EijR(40,4)-EijR(42,4)-EijR(44,4)
     -   )+EijR(5,2)*P(9)+EijR(5,3)*P(18)+(EijR(16,3)+EijR(17,3))*P(5
     -   2)-EijR(9,2)*P(66)+EijR(7,3)*P(160)+EijR(4,2)*P(205)-EijR(4,
     -   3)*P(206)-EijR(7,2)*P(207)+(-EijR(8,3)+EijR(10,3))*P(208)+Ei
     -   jR(14,3)*P(209)+EijR(15,3)*P(211)-EijR(18,3)*P(212)+(EijR(19
     -   ,3)+EijR(20,3))*P(213))
       FI(19) = 8*(Dij2345I(1,1)-Dij2345I(3,1)-Dij2345I(3,2)+Dij2345I(
     -   5,2)+p3sq*EijI(6,2)-p3sq*EijI(8,2)-4*(EijI(21,3)+EijI(22,3))
     -   +8*EijI(24,3)-2*(EijI(39,4)+EijI(40,4)-EijI(42,4)-EijI(44,4)
     -   )+EijI(5,2)*P(9)+EijI(5,3)*P(18)+(EijI(16,3)+EijI(17,3))*P(5
     -   2)-EijI(9,2)*P(66)+EijI(7,3)*P(160)+EijI(4,2)*P(205)-EijI(4,
     -   3)*P(206)-EijI(7,2)*P(207)+(-EijI(8,3)+EijI(10,3))*P(208)+Ei
     -   jI(14,3)*P(209)+EijI(15,3)*P(211)-EijI(18,3)*P(212)+(EijI(19
     -   ,3)+EijI(20,3))*P(213))
       F(19)=DCMPLX(FR(19),FI(19))
       P(214) = p2sq-s23+s45
       P(215) = p3sq+p4sq-s12+2*P(214)
       P(216) = p2sq-p3sq+p4sq+2*s12
       FR(20) = 8*(Dij2345R(1,1)-Dij2345R(3,1)+s12*EijR(2,1)-p3sq*(Eij
     -   R(3,1)-EijR(16,3)+EijR(20,3))+2*(2*EijR(11,2)-EijR(39,4)+Eij
     -   R(44,4))-EijR(5,2)*P(18)-EijR(8,2)*P(25)+EijR(4,1)*P(42)-Eij
     -   R(10,2)*P(43)+EijR(14,3)*P(45)+EijR(18,3)*P(46)+EijR(2,2)*P(
     -   68)-EijR(7,2)*P(72)+(-EijR(4,3)+EijR(15,3))*P(78)+EijR(4,2)*
     -   P(215)-EijR(9,2)*P(216))
       FI(20) = 8*(Dij2345I(1,1)-Dij2345I(3,1)+s12*EijI(2,1)-p3sq*(Eij
     -   I(3,1)-EijI(16,3)+EijI(20,3))+2*(2*EijI(11,2)-EijI(39,4)+Eij
     -   I(44,4))-EijI(5,2)*P(18)-EijI(8,2)*P(25)+EijI(4,1)*P(42)-Eij
     -   I(10,2)*P(43)+EijI(14,3)*P(45)+EijI(18,3)*P(46)+EijI(2,2)*P(
     -   68)-EijI(7,2)*P(72)+(-EijI(4,3)+EijI(15,3))*P(78)+EijI(4,2)*
     -   P(215)-EijI(9,2)*P(216))
       F(20)=DCMPLX(FR(20),FI(20))
       P(217) = s12-s34-2*s45
       P(218) = p2sq-s12-s34-s45+2*P(6)
       FR(21) = 8*(-Dij2345R(1,1)+Dij2345R(3,1)+p3sq*(EijR(16,3)-EijR(
     -   20,3))-2*(EijR(11,2)-EijR(22,3)+EijR(39,4)-EijR(44,4))+EijR(
     -   4,1)*P(5)+EijR(2,1)*P(9)+EijR(5,2)*P(18)-EijR(7,2)*P(36)+Eij
     -   R(10,2)*P(43)+EijR(14,3)*P(45)+EijR(18,3)*P(46)+EijR(8,2)*P(
     -   52)+(-EijR(4,3)+EijR(15,3))*P(78)-EijR(2,2)*P(208)+EijR(4,2)
     -   *P(217)-EijR(9,2)*P(218))
       FI(21) = 8*(-Dij2345I(1,1)+Dij2345I(3,1)+p3sq*(EijI(16,3)-EijI(
     -   20,3))-2*(EijI(11,2)-EijI(22,3)+EijI(39,4)-EijI(44,4))+EijI(
     -   4,1)*P(5)+EijI(2,1)*P(9)+EijI(5,2)*P(18)-EijI(7,2)*P(36)+Eij
     -   I(10,2)*P(43)+EijI(14,3)*P(45)+EijI(18,3)*P(46)+EijI(8,2)*P(
     -   52)+(-EijI(4,3)+EijI(15,3))*P(78)-EijI(2,2)*P(208)+EijI(4,2)
     -   *P(217)-EijI(9,2)*P(218))
       F(21)=DCMPLX(FR(21),FI(21))
       P(219) = s12-s34-s45
       P(220) = p2sq+2*p3sq-s12-s23+s45
       FR(22) = -8*(Dij2345R(2,1)-Dij2345R(3,1)+p3sq*EijR(13,3)+4*EijR
     -   (24,3)-2*(EijR(23,3)-EijR(39,4)+EijR(45,4))-EijR(4,1)*P(5)+(
     -   -EijR(3,1)+EijR(9,2))*P(9)-EijR(6,2)*P(18)-EijR(14,3)*P(45)-
     -   EijR(19,3)*P(46)-EijR(3,2)*P(52)+EijR(4,3)*P(78)+EijR(8,2)*P
     -   (208)+EijR(10,2)*P(217)-EijR(4,2)*P(219)-EijR(16,3)*P(220))
       FI(22) = -8*(Dij2345I(2,1)-Dij2345I(3,1)+p3sq*EijI(13,3)+4*EijI
     -   (24,3)-2*(EijI(23,3)-EijI(39,4)+EijI(45,4))-EijI(4,1)*P(5)+(
     -   -EijI(3,1)+EijI(9,2))*P(9)-EijI(6,2)*P(18)-EijI(14,3)*P(45)-
     -   EijI(19,3)*P(46)-EijI(3,2)*P(52)+EijI(4,3)*P(78)+EijI(8,2)*P
     -   (208)+EijI(10,2)*P(217)-EijI(4,2)*P(219)-EijI(16,3)*P(220))
       F(22)=DCMPLX(FR(22),FI(22))
       FR(23) = 8*(p3sq*(EijR(10,2)+EijR(16,3))-2*(EijR(24,3)+EijR(39,
     -   4))-EijR(9,2)*P(1)+EijR(4,2)*P(10)+EijR(14,3)*P(45)-EijR(4,3
     -   )*P(78))
       FI(23) = 8*(p3sq*(EijI(10,2)+EijI(16,3))-2*(EijI(24,3)+EijI(39,
     -   4))-EijI(9,2)*P(1)+EijI(4,2)*P(10)+EijI(14,3)*P(45)-EijI(4,3
     -   )*P(78))
       F(23)=DCMPLX(FR(23),FI(23))
       P(221) = p2sq-s12-s15
       P(222) = p2sq-2*s12-s15+s34
       P(223) = -s12+s34+2*P(214)
       FR(24) = 8*(D02345R+Dij2345R(2,1)+s12*EijR(2,1)-p3sq*EijR(13,3)
     -   -2*(EijR(11,2)+EijR(24,3)+EijR(39,4)-EijR(45,4))-EijR(6,2)*P
     -   (18)-EijR(3,2)*P(25)-EijR(3,1)*P(30)+EijR(14,3)*P(45)+EijR(1
     -   9,3)*P(46)+EijR(4,2)*P(48)-EijR(4,1)*P(61)+EijR(8,2)*P(68)-E
     -   ijR(4,3)*P(78)+EijR(16,3)*P(220)-EijR(7,2)*P(221)+EijR(9,2)*
     -   P(222)-EijR(10,2)*P(223))
       FI(24) = 8*(D02345I+Dij2345I(2,1)+s12*EijI(2,1)-p3sq*EijI(13,3)
     -   -2*(EijI(11,2)+EijI(24,3)+EijI(39,4)-EijI(45,4))-EijI(6,2)*P
     -   (18)-EijI(3,2)*P(25)-EijI(3,1)*P(30)+EijI(14,3)*P(45)+EijI(1
     -   9,3)*P(46)+EijI(4,2)*P(48)-EijI(4,1)*P(61)+EijI(8,2)*P(68)-E
     -   ijI(4,3)*P(78)+EijI(16,3)*P(220)-EijI(7,2)*P(221)+EijI(9,2)*
     -   P(222)-EijI(10,2)*P(223))
       F(24)=DCMPLX(FR(24),FI(24))
       P(224) = p2sq+p3sq-s23-s34-s45
       P(225) = -p2sq-p3sq+p4sq+s12+s23
       P(226) = p2sq-s23-s45
       P(227) = p2sq+p3sq-p4sq-s23-s45
       P(228) = 2*p2sq+p3sq-p4sq-s12-s23-s45
       FR(25) = 8*(D02345R+Dij2345R(2,1)-Dij2345R(3,2)+Dij2345R(6,2)-p
     -   3sq*EijR(3,2)+s12*(-EijR(5,2)+EijR(9,2))-2*(2*(EijR(11,2)+Ei
     -   jR(23,3)-EijR(24,3))+EijR(39,4)+EijR(41,4)-EijR(42,4)-EijR(4
     -   5,4))+EijR(6,3)*P(18)+EijR(6,2)*P(34)+EijR(7,2)*P(35)+EijR(1
     -   1,3)*P(52)+EijR(7,3)*P(160)+(-EijR(15,3)-EijR(17,3)+EijR(18,
     -   3)+EijR(20,3))*P(208)+EijR(14,3)*P(209)+EijR(13,3)*P(213)+Ei
     -   jR(4,2)*P(224)+EijR(4,3)*P(225)-EijR(10,2)*P(226)+EijR(16,3)
     -   *P(227)-EijR(19,3)*P(228))
       FI(25) = 8*(D02345I+Dij2345I(2,1)-Dij2345I(3,2)+Dij2345I(6,2)-p
     -   3sq*EijI(3,2)+s12*(-EijI(5,2)+EijI(9,2))-2*(2*(EijI(11,2)+Ei
     -   jI(23,3)-EijI(24,3))+EijI(39,4)+EijI(41,4)-EijI(42,4)-EijI(4
     -   5,4))+EijI(6,3)*P(18)+EijI(6,2)*P(34)+EijI(7,2)*P(35)+EijI(1
     -   1,3)*P(52)+EijI(7,3)*P(160)+(-EijI(15,3)-EijI(17,3)+EijI(18,
     -   3)+EijI(20,3))*P(208)+EijI(14,3)*P(209)+EijI(13,3)*P(213)+Ei
     -   jI(4,2)*P(224)+EijI(4,3)*P(225)-EijI(10,2)*P(226)+EijI(16,3)
     -   *P(227)-EijI(19,3)*P(228))
       F(25)=DCMPLX(FR(25),FI(25))
       P(229) = 1-Is45*p3sq
       P(230) = s12-s34
       P(231) = -s15+2*P(18)
       P(232) = p2sq-s12-s15+s23-s45
       P(233) = p3sq-s23-s34
       P(234) = s15+2*P(233)
       FR(26) = -8*(-D02345R+Dij2345R(1,1)-Dij2345R(3,1)+p3sq*(EijR(3,
     -   2)-EijR(12,3)+EijR(13,3))+2*(4*EijR(11,2)-EijR(22,3)+EijR(23
     -   ,3)+EijR(39,4)+EijR(43,4)-EijR(44,4)-EijR(45,4)-EijR(8,2)*P(
     -   1))-(-EijR(6,2)+EijR(14,3)+EijR(17,3))*P(45)+Is45*(Cij123R(2
     -   ,1)+Cij123R(2,2)-p3sq*(Dij1234R(3,1)+Dij1234R(6,2))+2*(Dij12
     -   34R(7,2)+Dij1234R(12,3))+Dij1234R(2,2)*P(9)-Dij1234R(4,2)*P(
     -   45))-(EijR(18,3)+EijR(19,3))*P(46)-(-EijR(4,3)+EijR(15,3))*P
     -   (78)+EijR(2,2)*P(90)+EijR(10,2)*P(132)-EijR(4,1)*P(153)-(Eij
     -   R(16,3)-EijR(20,3))*P(220)-EijR(4,2)*P(224)-Dij1234R(2,1)*P(
     -   229)-EijR(2,1)*P(230)-EijR(5,2)*P(231)+EijR(7,2)*P(232)+EijR
     -   (9,2)*P(234))
       FI(26) = -8*(-D02345I+Dij2345I(1,1)-Dij2345I(3,1)+p3sq*(EijI(3,
     -   2)-EijI(12,3)+EijI(13,3))+2*(4*EijI(11,2)-EijI(22,3)+EijI(23
     -   ,3)+EijI(39,4)+EijI(43,4)-EijI(44,4)-EijI(45,4)-EijI(8,2)*P(
     -   1))-(-EijI(6,2)+EijI(14,3)+EijI(17,3))*P(45)+Is45*(Cij123I(2
     -   ,1)+Cij123I(2,2)-p3sq*(Dij1234I(3,1)+Dij1234I(6,2))+2*(Dij12
     -   34I(7,2)+Dij1234I(12,3))+Dij1234I(2,2)*P(9)-Dij1234I(4,2)*P(
     -   45))-(EijI(18,3)+EijI(19,3))*P(46)-(-EijI(4,3)+EijI(15,3))*P
     -   (78)+EijI(2,2)*P(90)+EijI(10,2)*P(132)-EijI(4,1)*P(153)-(Eij
     -   I(16,3)-EijI(20,3))*P(220)-EijI(4,2)*P(224)-Dij1234I(2,1)*P(
     -   229)-EijI(2,1)*P(230)-EijI(5,2)*P(231)+EijI(7,2)*P(232)+EijI
     -   (9,2)*P(234))
       F(26)=DCMPLX(FR(26),FI(26))
       P(235) = p4sq+s12-s34-s45
       P(236) = p2sq+p3sq+s12-s23-s45+2*P(128)
       FR(27) = -8*(-(p3sq*(EijR(8,2)+EijR(9,3)-EijR(10,2)+EijR(16,3)-
     -   2*EijR(20,3)))+2*(EijR(22,3)-EijR(24,3)+EijR(37,4)+EijR(39,4
     -   )-2*EijR(44,4))-(EijR(5,2)-EijR(7,2)+EijR(8,3)+EijR(14,3)-2*
     -   EijR(18,3))*P(45)-EijR(4,2)*P(62)+(EijR(4,3)+EijR(10,3)-2*Ei
     -   jR(15,3))*P(78)+Is45*(Cij123R(2,1)+Cij123R(2,2)+p3sq*(Dij123
     -   4R(2,1)-Dij1234R(3,1)-Dij1234R(6,2))+2*(Dij1234R(7,2)+Dij123
     -   4R(12,3))-Dij1234R(4,2)*P(45)+Dij1234R(2,2)*P(78))-(EijR(2,1
     -   )+EijR(2,2)-EijR(4,1))*P(235)+EijR(9,2)*P(236))
       FI(27) = -8*(-(p3sq*(EijI(8,2)+EijI(9,3)-EijI(10,2)+EijI(16,3)-
     -   2*EijI(20,3)))+2*(EijI(22,3)-EijI(24,3)+EijI(37,4)+EijI(39,4
     -   )-2*EijI(44,4))-(EijI(5,2)-EijI(7,2)+EijI(8,3)+EijI(14,3)-2*
     -   EijI(18,3))*P(45)-EijI(4,2)*P(62)+(EijI(4,3)+EijI(10,3)-2*Ei
     -   jI(15,3))*P(78)+Is45*(Cij123I(2,1)+Cij123I(2,2)+p3sq*(Dij123
     -   4I(2,1)-Dij1234I(3,1)-Dij1234I(6,2))+2*(Dij1234I(7,2)+Dij123
     -   4I(12,3))-Dij1234I(4,2)*P(45)+Dij1234I(2,2)*P(78))-(EijI(2,1
     -   )+EijI(2,2)-EijI(4,1))*P(235)+EijI(9,2)*P(236))
       F(27)=DCMPLX(FR(27),FI(27))
       P(237) = -2+Is45*s23
       P(238) = p2sq-p3sq+s23-s45
       P(239) = 2*p2sq-p3sq+s23-s45
       P(240) = p3sq-p4sq-s12+s34+s45
       P(241) = -p2sq+p3sq+p4sq+s12+s23-s34-2*s45
       P(242) = p2sq-s12-s23+s34+2*s45
       P(243) = p2sq+p4sq+s12-s23-s34
       P(244) = -p2sq+p4sq+s12+s23-s34+2*P(30)
       P(245) = p2sq+p3sq+p4sq+s12-s23-s34
       P(246) = p2sq-s12-2*s23+s34+3*s45
       P(247) = p2sq+s45
       P(248) = -s23+s34+2*P(247)
       P(249) = s34+3*s45+2*P(10)
       FR(28) = -8*(-D01234R+Dij2345R(1,1)-Dij2345R(3,1)-Dij2345R(3,2)
     -   +Dij2345R(5,2)-p3sq*(-EijR(3,1)-2*EijR(6,2)+EijR(8,2)+EijR(1
     -   0,2))+2*(EijR(11,2)-EijR(21,3)+EijR(24,3)+EijR(39,4)+EijR(40
     -   ,4)-EijR(42,4)-EijR(44,4))-(-EijR(5,3)+EijR(7,3))*P(36)-Dij1
     -   234R(1,2)*P(40)-(-EijR(16,3)-EijR(17,3)+EijR(19,3)+EijR(20,3
     -   ))*P(52)-(EijR(8,3)-EijR(10,3))*P(68)+Is12*(Cij145R(1,1)+Cij
     -   145R(1,2)-Cij145R(2,1)+Cij145R(2,2)+p3sq*(Dij1345R(2,1)-Dij1
     -   345R(3,1)+Dij1345R(4,2)-Dij1345R(6,2))-2*(Cij145R(3,2)-Dij13
     -   45R(7,2)-Dij1345R(11,3)+Dij1345R(13,3))+(Dij1345R(3,2)-Dij13
     -   45R(5,2))*P(203))-EijR(2,1)*P(235)+Dij1234R(1,1)*P(237)+Is45
     -   *(C0123R-C0134R+Cij123R(1,1)+Cij123R(2,1)+Cij123R(3,2)+6*Dij
     -   1234R(7,2)+2*(p2sq*Dij1234R(2,2)+Dij1234R(11,3))-Dij1234R(6,
     -   2)*P(1)+Dij1234R(5,2)*P(10)-Dij1234R(3,1)*P(22)-Dij1234R(2,1
     -   )*P(238)-Dij1234R(4,2)*P(239))-EijR(4,1)*P(240)+EijR(4,2)*P(
     -   241)+EijR(4,3)*P(242)-EijR(5,2)*P(243)-EijR(7,2)*P(244)+EijR
     -   (9,2)*P(245)-EijR(14,3)*P(246)-EijR(15,3)*P(248)+EijR(18,3)*
     -   P(249))
       FI(28) = -8*(-D01234I+Dij2345I(1,1)-Dij2345I(3,1)-Dij2345I(3,2)
     -   +Dij2345I(5,2)-p3sq*(-EijI(3,1)-2*EijI(6,2)+EijI(8,2)+EijI(1
     -   0,2))+2*(EijI(11,2)-EijI(21,3)+EijI(24,3)+EijI(39,4)+EijI(40
     -   ,4)-EijI(42,4)-EijI(44,4))-(-EijI(5,3)+EijI(7,3))*P(36)-Dij1
     -   234I(1,2)*P(40)-(-EijI(16,3)-EijI(17,3)+EijI(19,3)+EijI(20,3
     -   ))*P(52)-(EijI(8,3)-EijI(10,3))*P(68)+Is12*(Cij145I(1,1)+Cij
     -   145I(1,2)-Cij145I(2,1)+Cij145I(2,2)+p3sq*(Dij1345I(2,1)-Dij1
     -   345I(3,1)+Dij1345I(4,2)-Dij1345I(6,2))-2*(Cij145I(3,2)-Dij13
     -   45I(7,2)-Dij1345I(11,3)+Dij1345I(13,3))+(Dij1345I(3,2)-Dij13
     -   45I(5,2))*P(203))-EijI(2,1)*P(235)+Dij1234I(1,1)*P(237)+Is45
     -   *(C0123I-C0134I+Cij123I(1,1)+Cij123I(2,1)+Cij123I(3,2)+6*Dij
     -   1234I(7,2)+2*(p2sq*Dij1234I(2,2)+Dij1234I(11,3))-Dij1234I(6,
     -   2)*P(1)+Dij1234I(5,2)*P(10)-Dij1234I(3,1)*P(22)-Dij1234I(2,1
     -   )*P(238)-Dij1234I(4,2)*P(239))-EijI(4,1)*P(240)+EijI(4,2)*P(
     -   241)+EijI(4,3)*P(242)-EijI(5,2)*P(243)-EijI(7,2)*P(244)+EijI
     -   (9,2)*P(245)-EijI(14,3)*P(246)-EijI(15,3)*P(248)+EijI(18,3)*
     -   P(249))
       F(28)=DCMPLX(FR(28),FI(28))
       P(250) = p2sq+p3sq-s23-s34
       P(251) = -p2sq+s12+s15
       P(252) = -s34+2*P(1)
       P(253) = p2sq+3*p3sq-s12-s23+s45
       P(254) = 3*p3sq+2*P(45)
       FR(29) = -8*(Dij2345R(2,1)-Dij2345R(3,1)-p3sq*EijR(3,3)-4*EijR(
     -   23,3)-EijR(3,2)*P(1)-(-EijR(8,2)+EijR(9,2))*P(8)-(EijR(11,3)
     -   +EijR(14,3))*P(45)+2*(EijR(38,4)+EijR(39,4)+2*(EijR(24,3)-Ei
     -   jR(45,4))-EijR(19,3)*P(46))+EijR(4,3)*P(78)-EijR(4,1)*P(153)
     -   +Is45*(C0134R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,2)-2*(Dij12
     -   34R(7,2)-Dij1234R(13,3))+Dij1234R(6,2)*P(9)-Dij1234R(5,2)*P(
     -   45)-Dij1234R(3,1)*P(214))-EijR(3,1)*P(230)-EijR(4,2)*P(250)-
     -   (-EijR(6,2)+EijR(7,2))*P(251)+EijR(10,2)*P(252)+EijR(13,3)*P
     -   (253)-EijR(16,3)*P(254))
       FI(29) = -8*(Dij2345I(2,1)-Dij2345I(3,1)-p3sq*EijI(3,3)-4*EijI(
     -   23,3)-EijI(3,2)*P(1)-(-EijI(8,2)+EijI(9,2))*P(8)-(EijI(11,3)
     -   +EijI(14,3))*P(45)+2*(EijI(38,4)+EijI(39,4)+2*(EijI(24,3)-Ei
     -   jI(45,4))-EijI(19,3)*P(46))+EijI(4,3)*P(78)-EijI(4,1)*P(153)
     -   +Is45*(C0134I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,2)-2*(Dij12
     -   34I(7,2)-Dij1234I(13,3))+Dij1234I(6,2)*P(9)-Dij1234I(5,2)*P(
     -   45)-Dij1234I(3,1)*P(214))-EijI(3,1)*P(230)-EijI(4,2)*P(250)-
     -   (-EijI(6,2)+EijI(7,2))*P(251)+EijI(10,2)*P(252)+EijI(13,3)*P
     -   (253)-EijI(16,3)*P(254))
       F(29)=DCMPLX(FR(29),FI(29))
       P(255) = -p4sq-s12+s34+s45
       P(256) = p2sq+p3sq+p4sq-s23-s34+s45
       P(257) = p2sq+p3sq-p4sq-s12-s23+s34+s45
       P(258) = p3sq-p4sq-2*s12-s23+s34+s45
       P(259) = p2sq+3*p3sq+p4sq-s12-s23-s34+s45
       FR(30) = -8*(D02345R-p2sq*EijR(2,2)+p3sq*(-2*EijR(3,2)-EijR(12,
     -   3)+EijR(13,3))-2*(4*EijR(11,2)+EijR(22,3)-EijR(24,3)-EijR(39
     -   ,4)-EijR(43,4)+EijR(44,4)+EijR(45,4))+EijR(5,2)*P(18)-(2*Eij
     -   R(6,2)+EijR(14,3)+EijR(17,3))*P(45)-(EijR(18,3)+EijR(19,3))*
     -   P(46)+EijR(7,2)*P(72)+(EijR(4,3)-EijR(15,3))*P(78)+Is45*(C01
     -   34R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,2)-2*(Dij1234R(7,2)-D
     -   ij1234R(13,3))-Dij1234R(3,1)*P(10)-Dij1234R(5,2)*P(45)+Dij12
     -   34R(6,2)*P(78))+(-EijR(16,3)+EijR(20,3))*P(220)-EijR(3,1)*P(
     -   235)-EijR(4,1)*P(255)-EijR(4,2)*P(256)+EijR(8,2)*P(257)-EijR
     -   (9,2)*P(258)+EijR(10,2)*P(259))
       FI(30) = -8*(D02345I-p2sq*EijI(2,2)+p3sq*(-2*EijI(3,2)-EijI(12,
     -   3)+EijI(13,3))-2*(4*EijI(11,2)+EijI(22,3)-EijI(24,3)-EijI(39
     -   ,4)-EijI(43,4)+EijI(44,4)+EijI(45,4))+EijI(5,2)*P(18)-(2*Eij
     -   I(6,2)+EijI(14,3)+EijI(17,3))*P(45)-(EijI(18,3)+EijI(19,3))*
     -   P(46)+EijI(7,2)*P(72)+(EijI(4,3)-EijI(15,3))*P(78)+Is45*(C01
     -   34I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,2)-2*(Dij1234I(7,2)-D
     -   ij1234I(13,3))-Dij1234I(3,1)*P(10)-Dij1234I(5,2)*P(45)+Dij12
     -   34I(6,2)*P(78))+(-EijI(16,3)+EijI(20,3))*P(220)-EijI(3,1)*P(
     -   235)-EijI(4,1)*P(255)-EijI(4,2)*P(256)+EijI(8,2)*P(257)-EijI
     -   (9,2)*P(258)+EijI(10,2)*P(259))
       F(30)=DCMPLX(FR(30),FI(30))
       P(260) = p3sq-p4sq+s12-s34-s45
       P(261) = -p3sq+s34+s45
       P(262) = p2sq-p3sq+s12+s23
       P(263) = s12-s23
       P(264) = p3sq-2*p4sq+s12-s45
       P(265) = 2*s12+s15-s34
       P(266) = p3sq-s12+s34-s45
       P(267) = -p2sq+p3sq+s12+s23-s34-3*s45
       P(268) = 3*s12+s15-s34
       P(269) = -p2sq+s23-3*s45+2*P(147)
       P(270) = -p2sq+2*p3sq+s12+s15+s23-s34-3*s45
       P(271) = p2sq-2*s12-s23+3*s45
       P(272) = p2sq-2*s12-s23+s34+3*s45
       P(273) = p2sq+s34+4*s45-2*P(177)
       P(274) = -1+Is45*s23
       FR(31) = 4*Is12s45*(3*B013R+B014R-2*(p2sq*(C0123R+Cij123R(1,1))
     -   +2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2))+(C0145R+Cij145R(
     -   2,1))*P(15)-Cij123R(2,1)*P(18)+Cij134R(2,1)*P(25)+C0134R*P(4
     -   2)))-8*(D02345R+Dij2345R(2,1)-Dij2345R(3,2)+Dij2345R(6,2)+s1
     -   5*(EE0R+2*EijR(1,1)+EijR(1,2))-p3sq*EijR(3,2)+EijR(3,1)*P(25
     -   )-(-EijR(6,3)+EijR(7,3))*P(36)-(D01345R+Dij1345R(1,1))*P(39)
     -   -(D01234R+Dij1234R(1,2))*P(40)-(-EijR(11,3)+EijR(13,3))*P(52
     -   )-(EijR(15,3)+EijR(17,3)-EijR(18,3)-EijR(20,3))*P(68)+EijR(4
     -   ,3)*P(242)-EijR(14,3)*P(246)+Is45*(Cij123R(3,2)-2*(Cij134R(1
     -   ,1)-p2sq*Dij1234R(6,2)-Dij1234R(11,3))-Dij1234R(3,2)*P(1)-Di
     -   j1234R(3,1)*P(42)-Dij1234R(2,1)*P(208)-Dij1234R(4,2)*P(262)+
     -   Dij1234R(5,2)*P(263))+Is12*(Cij145R(1,2)+Cij145R(2,2)-2*(Cij
     -   145R(3,2)+Dij1345R(7,2)-Dij1345R(12,3)+Dij1345R(13,3))+(Dij1
     -   345R(2,1)-Dij1345R(3,1))*P(30)-Dij1345R(2,2)*P(128)+Dij1345R
     -   (3,2)*P(260)+(-Dij1345R(4,2)+Dij1345R(5,2))*P(261)-Dij1345R(
     -   6,2)*P(264))-EijR(2,1)*P(265)-EijR(4,1)*P(266)+EijR(4,2)*P(2
     -   67)-(EijR(5,2)-EijR(9,2))*P(268)+EijR(6,2)*P(269)-EijR(7,2)*
     -   P(270)+EijR(10,2)*P(271)-EijR(16,3)*P(272)+EijR(19,3)*P(273)
     -   -2*(3*EijR(11,2)+EijR(21,3)-EijR(24,3)-EijR(39,4)-EijR(41,4)
     -   +EijR(42,4)+EijR(45,4)-Dij1234R(1,1)*P(274)))
       FI(31) = 4*Is12s45*(3*B013I+B014I-2*(p2sq*(C0123I+Cij123I(1,1))
     -   +2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2))+(C0145I+Cij145I(
     -   2,1))*P(15)-Cij123I(2,1)*P(18)+Cij134I(2,1)*P(25)+C0134I*P(4
     -   2)))-8*(D02345I+Dij2345I(2,1)-Dij2345I(3,2)+Dij2345I(6,2)+s1
     -   5*(EE0I+2*EijI(1,1)+EijI(1,2))-p3sq*EijI(3,2)+EijI(3,1)*P(25
     -   )-(-EijI(6,3)+EijI(7,3))*P(36)-(D01345I+Dij1345I(1,1))*P(39)
     -   -(D01234I+Dij1234I(1,2))*P(40)-(-EijI(11,3)+EijI(13,3))*P(52
     -   )-(EijI(15,3)+EijI(17,3)-EijI(18,3)-EijI(20,3))*P(68)+EijI(4
     -   ,3)*P(242)-EijI(14,3)*P(246)+Is45*(Cij123I(3,2)-2*(Cij134I(1
     -   ,1)-p2sq*Dij1234I(6,2)-Dij1234I(11,3))-Dij1234I(3,2)*P(1)-Di
     -   j1234I(3,1)*P(42)-Dij1234I(2,1)*P(208)-Dij1234I(4,2)*P(262)+
     -   Dij1234I(5,2)*P(263))+Is12*(Cij145I(1,2)+Cij145I(2,2)-2*(Cij
     -   145I(3,2)+Dij1345I(7,2)-Dij1345I(12,3)+Dij1345I(13,3))+(Dij1
     -   345I(2,1)-Dij1345I(3,1))*P(30)-Dij1345I(2,2)*P(128)+Dij1345I
     -   (3,2)*P(260)+(-Dij1345I(4,2)+Dij1345I(5,2))*P(261)-Dij1345I(
     -   6,2)*P(264))-EijI(2,1)*P(265)-EijI(4,1)*P(266)+EijI(4,2)*P(2
     -   67)-(EijI(5,2)-EijI(9,2))*P(268)+EijI(6,2)*P(269)-EijI(7,2)*
     -   P(270)+EijI(10,2)*P(271)-EijI(16,3)*P(272)+EijI(19,3)*P(273)
     -   -2*(3*EijI(11,2)+EijI(21,3)-EijI(24,3)-EijI(39,4)-EijI(41,4)
     -   +EijI(42,4)+EijI(45,4)-Dij1234I(1,1)*P(274)))
       F(31)=DCMPLX(FR(31),FI(31))
       P(275) = p3sq-s12+2*P(10)
       P(276) = s12-2*P(1)
       FR(32) = 8*(-Dij2345R(3,1)-Dij2345R(3,2)+p3sq*EijR(10,2)+4*EijR
     -   (24,3)-2*(EijR(39,4)-EijR(42,4))+EijR(16,3)*P(52)+EijR(7,3)*
     -   P(160)-EijR(4,3)*P(206)+(-EijR(15,3)+EijR(18,3))*P(208)+EijR
     -   (14,3)*P(209)+EijR(19,3)*P(213)+EijR(4,2)*P(275)+EijR(7,2)*P
     -   (276))
       FI(32) = 8*(-Dij2345I(3,1)-Dij2345I(3,2)+p3sq*EijI(10,2)+4*EijI
     -   (24,3)-2*(EijI(39,4)-EijI(42,4))+EijI(16,3)*P(52)+EijI(7,3)*
     -   P(160)-EijI(4,3)*P(206)+(-EijI(15,3)+EijI(18,3))*P(208)+EijI
     -   (14,3)*P(209)+EijI(19,3)*P(213)+EijI(4,2)*P(275)+EijI(7,2)*P
     -   (276))
       F(32)=DCMPLX(FR(32),FI(32))
       P(277) = p4sq-s12+s34
       P(278) = -p2sq+s12+s23-s34-2*s45
       P(279) = 2*p3sq-s45
       FR(33) = -8*(D01345R-Dij2345R(3,1)-Dij2345R(3,2)+p3sq*(EijR(3,1
     -   )+EijR(6,2))-2*(EijR(11,2)+EijR(21,3)-EijR(24,3)-EijR(39,4)+
     -   EijR(42,4))-(EijR(2,1)+EijR(5,2)-EijR(9,2))*P(1)+EijR(4,1)*P
     -   (10)+EijR(4,2)*P(30)-EijR(7,3)*P(36)+(EijR(16,3)-EijR(19,3))
     -   *P(52)+(-EijR(15,3)+EijR(18,3))*P(68)-EijR(14,3)*P(246)+Is12
     -   *(C0134R-C0145R-Cij145R(1,1)+Cij145R(2,2)-Cij145R(3,2)-4*Dij
     -   1345R(7,2)-2*Dij1345R(13,3)+(Dij1345R(2,1)+Dij1345R(6,2))*P(
     -   128)+Dij1345R(1,1)*P(136)+Dij1345R(3,2)*P(260)+Dij1345R(5,2)
     -   *P(261)-Dij1345R(3,1)*P(277))-EijR(4,3)*P(278)-EijR(7,2)*P(2
     -   79))
       FI(33) = -8*(D01345I-Dij2345I(3,1)-Dij2345I(3,2)+p3sq*(EijI(3,1
     -   )+EijI(6,2))-2*(EijI(11,2)+EijI(21,3)-EijI(24,3)-EijI(39,4)+
     -   EijI(42,4))-(EijI(2,1)+EijI(5,2)-EijI(9,2))*P(1)+EijI(4,1)*P
     -   (10)+EijI(4,2)*P(30)-EijI(7,3)*P(36)+(EijI(16,3)-EijI(19,3))
     -   *P(52)+(-EijI(15,3)+EijI(18,3))*P(68)-EijI(14,3)*P(246)+Is12
     -   *(C0134I-C0145I-Cij145I(1,1)+Cij145I(2,2)-Cij145I(3,2)-4*Dij
     -   1345I(7,2)-2*Dij1345I(13,3)+(Dij1345I(2,1)+Dij1345I(6,2))*P(
     -   128)+Dij1345I(1,1)*P(136)+Dij1345I(3,2)*P(260)+Dij1345I(5,2)
     -   *P(261)-Dij1345I(3,1)*P(277))-EijI(4,3)*P(278)-EijI(7,2)*P(2
     -   79))
       F(33)=DCMPLX(FR(33),FI(33))
       P(280) = s45*P(1)+p3sq*P(9)
       P(281) = -p2sq-p3sq+p4sq+s12+s23-s34-s45+Is12*s34*P(1)
       P(282) = -p2sq-2*p3sq+p4sq+s12+s23
       P(283) = s34*P(125)+s12*P(282)
       P(284) = -p4sq+s34+s45
       P(285) = p2sq*P(106)+p3sq*P(235)+s23*P(284)
       P(286) = p4sq*P(1)-s45*P(9)
       P(287) = p3sq+p4sq-s34
       P(288) = s12-s23-s34-s45
       P(289) = p3sq**2+p4sq*s12-s12*s34+s23*s34+s23*s45+p2sq*P(95)+p3
     -   sq*P(288)
       P(290) = s12+s23+s34+s45
       P(291) = p3sq**2+2*p4sq*s12+s12*s34+s23*s34+p2sq*P(95)+s45*P(17
     -   7)-p3sq*P(290)
       P(292) = s12+3*s23
       P(293) = p3sq+3*p4sq-s34
       P(294) = 3*p4sq-s12-s23-s34
       P(295) = p3sq**2+s34*P(177)-p4sq*P(292)+p2sq*P(293)+p3sq*P(294)
       P(296) = p3sq-s12-s15-s45
       P(297) = s12*P(161)+p2sq*P(296)
       P(298) = -p3sq+s12
       P(299) = s45**2-s12*P(163)+p3sq*P(164)-s45*P(165)+p2sq*P(298)
       P(300) = p4sq+3*s23-s34+2*P(52)
       P(301) = -3*p3sq-s12+s45
       P(302) = -2*p3sq**2+P(52)*P(168)+p3sq*P(300)+p2sq*P(301)
       P(303) = -p4sq+s45
       P(304) = s12+s34+s45
       P(305) = p2sq*s45+p4sq*P(171)+p3sq*P(303)-s45*P(304)
       P(306) = -p3sq+s15+s23-s34-2*P(18)
       P(307) = -2*p3sq-s34-s45+3*P(177)
       P(308) = 2*p2sq**2-p2sq*P(176)+s12*P(307)
       P(309) = s12*P(177)
       P(310) = p4sq*s12-s12*s15+s23**2-s12*s34-4*s12*s45+s15*s45+s45*
     -   *2-p3sq*P(180)+p2sq*P(181)+2*P(182)+3*P(309)
       P(311) = -(p4sq*s12)-4*s12**2-s12*s15+s23**2+p2sq*P(187)+2*P(18
     -   8)+s45*P(189)-p3sq*P(190)
       P(312) = 8*p2sq+5*p3sq-p4sq-7*s12-s15-6*s23+3*s34+2*s45
       FR(34) = Is12s45*(2*(3*B013R+B014R)*P(5)+4*(P(9)*(-(p2sq*Cij123
     -   R(2,1))+p4sq*Cij145R(2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij14
     -   5R(4,2))+C0145R*P(15)+(C0123R+Cij123R(1,1))*P(18)+Cij134R(2,
     -   1)*P(25))+C0134R*P(280)))+4*(-4*EijR(46,4)+Dij2345R(2,1)*P(5
     -   3)+D02345R*P(152)+P(68)*(2*EijR(22,3)-EijR(2,2)*P(155))+s15*
     -   (EE0R*P(9)+EijR(1,2)*P(18)+EijR(1,1)*P(159))+EijR(7,2)*P(186
     -   )-P(155)*(Dij2345R(1,1)-Dij2345R(3,1)+EijR(9,2)*P(192))+EijR
     -   (10,2)*P(198)+Is45*(P(9)*(-2*p2sq*Dij1234R(2,2)+Dij1234R(4,2
     -   )*P(18)-Dij1234R(3,1)*P(45))-(Dij1234R(3,2)*P(1)+Dij1234R(5,
     -   2)*P(18))*P(78)+(D01234R+Dij1234R(1,1))*P(151)+Dij1234R(2,1)
     -   *P(156)+Dij1234R(6,2)*P(158)+2*(Cij134R(1,1)*P(5)-Dij1234R(7
     -   ,2)*P(199)))+Dij1345R(1,2)*P(261)+D01345R*P(281)+Is12*((Cij1
     -   45R(1,1)+2*p4sq*Dij1345R(3,2))*P(5)+Dij1345R(1,1)*P(283)+Dij
     -   1345R(2,1)*P(285)-Dij1345R(3,1)*P(286)-Dij1345R(2,2)*P(1)*P(
     -   287)+Dij1345R(4,2)*P(289)-Dij1345R(5,2)*P(291)+Dij1345R(6,2)
     -   *P(295))+EijR(2,1)*P(297)+EijR(3,1)*P(299)+EijR(3,2)*P(302)+
     -   EijR(4,1)*P(305)+EijR(4,2)*P(174)*P(306)+EijR(5,2)*P(308)-Ei
     -   jR(6,2)*P(310)+EijR(8,2)*P(311)+2*(Dij2345R(7,2)-EijR(21,3)*
     -   P(36)-EijR(23,3)*P(52)+Dij1345R(7,2)*P(109)-EijR(24,3)*P(242
     -   )-EijR(11,2)*P(312)))
       FI(34) = Is12s45*(2*(3*B013I+B014I)*P(5)+4*(P(9)*(-(p2sq*Cij123
     -   I(2,1))+p4sq*Cij145I(2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij14
     -   5I(4,2))+C0145I*P(15)+(C0123I+Cij123I(1,1))*P(18)+Cij134I(2,
     -   1)*P(25))+C0134I*P(280)))+4*(-4*EijI(46,4)+Dij2345I(2,1)*P(5
     -   3)+D02345I*P(152)+P(68)*(2*EijI(22,3)-EijI(2,2)*P(155))+s15*
     -   (EE0I*P(9)+EijI(1,2)*P(18)+EijI(1,1)*P(159))+EijI(7,2)*P(186
     -   )-P(155)*(Dij2345I(1,1)-Dij2345I(3,1)+EijI(9,2)*P(192))+EijI
     -   (10,2)*P(198)+Is45*(P(9)*(-2*p2sq*Dij1234I(2,2)+Dij1234I(4,2
     -   )*P(18)-Dij1234I(3,1)*P(45))-(Dij1234I(3,2)*P(1)+Dij1234I(5,
     -   2)*P(18))*P(78)+(D01234I+Dij1234I(1,1))*P(151)+Dij1234I(2,1)
     -   *P(156)+Dij1234I(6,2)*P(158)+2*(Cij134I(1,1)*P(5)-Dij1234I(7
     -   ,2)*P(199)))+Dij1345I(1,2)*P(261)+D01345I*P(281)+Is12*((Cij1
     -   45I(1,1)+2*p4sq*Dij1345I(3,2))*P(5)+Dij1345I(1,1)*P(283)+Dij
     -   1345I(2,1)*P(285)-Dij1345I(3,1)*P(286)-Dij1345I(2,2)*P(1)*P(
     -   287)+Dij1345I(4,2)*P(289)-Dij1345I(5,2)*P(291)+Dij1345I(6,2)
     -   *P(295))+EijI(2,1)*P(297)+EijI(3,1)*P(299)+EijI(3,2)*P(302)+
     -   EijI(4,1)*P(305)+EijI(4,2)*P(174)*P(306)+EijI(5,2)*P(308)-Ei
     -   jI(6,2)*P(310)+EijI(8,2)*P(311)+2*(Dij2345I(7,2)-EijI(21,3)*
     -   P(36)-EijI(23,3)*P(52)+Dij1345I(7,2)*P(109)-EijI(24,3)*P(242
     -   )-EijI(11,2)*P(312)))
       F(34)=DCMPLX(FR(34),FI(34))
       P(313) = s15-2*P(36)
       P(314) = p2sq-p3sq+p4sq-s15+s23+s34-s45
       P(315) = -p4sq-s12+s34+s45+2*P(10)
       P(316) = s15-s23+s45
       P(317) = p2sq+p4sq-s12-s15+s34
       FR(35) = 8*(D02345R-p2sq*EijR(2,2)+p4sq*(-2*EijR(4,2)+EijR(4,3)
     -   )-2*(2*EijR(11,2)+EijR(22,3)-3*EijR(24,3)+EijR(39,4)-EijR(44
     -   ,4))+EijR(8,2)*P(1)+EijR(10,3)*P(8)+EijR(5,2)*P(18)+(EijR(3,
     -   1)-EijR(4,1))*P(78)+(EijR(16,3)-EijR(20,3))*P(93)+EijR(7,2)*
     -   P(313)+EijR(9,2)*P(314)-EijR(10,2)*P(315)+(-EijR(14,3)+EijR(
     -   18,3))*P(316)-EijR(15,3)*P(317))
       FI(35) = 8*(D02345I-p2sq*EijI(2,2)+p4sq*(-2*EijI(4,2)+EijI(4,3)
     -   )-2*(2*EijI(11,2)+EijI(22,3)-3*EijI(24,3)+EijI(39,4)-EijI(44
     -   ,4))+EijI(8,2)*P(1)+EijI(10,3)*P(8)+EijI(5,2)*P(18)+(EijI(3,
     -   1)-EijI(4,1))*P(78)+(EijI(16,3)-EijI(20,3))*P(93)+EijI(7,2)*
     -   P(313)+EijI(9,2)*P(314)-EijI(10,2)*P(315)+(-EijI(14,3)+EijI(
     -   18,3))*P(316)-EijI(15,3)*P(317))
       F(35)=DCMPLX(FR(35),FI(35))
       P(318) = p2sq-2*p4sq-s12-s23+s34+s45
       FR(36) = 8*(Dij2345R(2,1)-Dij2345R(3,1)+s45*EijR(4,2)+p4sq*EijR
     -   (4,3)+4*EijR(24,3)-2*(s45*EijR(10,2)+EijR(23,3)+EijR(39,4)-E
     -   ijR(45,4))+(-EijR(15,3)+EijR(20,3))*P(8)+(-EijR(3,1)+EijR(4,
     -   1)+EijR(9,2))*P(9)-EijR(6,2)*P(18)-EijR(3,2)*P(52)+EijR(13,3
     -   )*P(139)-EijR(8,2)*P(210)+(-EijR(14,3)+EijR(19,3))*P(316)+Ei
     -   jR(16,3)*P(318))
       FI(36) = 8*(Dij2345I(2,1)-Dij2345I(3,1)+s45*EijI(4,2)+p4sq*EijI
     -   (4,3)+4*EijI(24,3)-2*(s45*EijI(10,2)+EijI(23,3)+EijI(39,4)-E
     -   ijI(45,4))+(-EijI(15,3)+EijI(20,3))*P(8)+(-EijI(3,1)+EijI(4,
     -   1)+EijI(9,2))*P(9)-EijI(6,2)*P(18)-EijI(3,2)*P(52)+EijI(13,3
     -   )*P(139)-EijI(8,2)*P(210)+(-EijI(14,3)+EijI(19,3))*P(316)+Ei
     -   jI(16,3)*P(318))
       F(36)=DCMPLX(FR(36),FI(36))
       P(319) = -p2sq+s12+s15-s34
       FR(37) = 8*(-Dij2345R(3,1)+s45*EijR(4,2)+p4sq*EijR(4,3)+4*EijR(
     -   24,3)-2*EijR(39,4)+EijR(9,2)*P(9)-EijR(10,2)*P(78)+EijR(16,3
     -   )*P(93)-EijR(14,3)*P(316)+EijR(15,3)*P(319))
       FI(37) = 8*(-Dij2345I(3,1)+s45*EijI(4,2)+p4sq*EijI(4,3)+4*EijI(
     -   24,3)-2*EijI(39,4)+EijI(9,2)*P(9)-EijI(10,2)*P(78)+EijI(16,3
     -   )*P(93)-EijI(14,3)*P(316)+EijI(15,3)*P(319))
       F(37)=DCMPLX(FR(37),FI(37))
       P(320) = p2sq-s45+2*P(210)
       FR(38) = 8*(Dij2345R(1,1)+p4sq*EijR(4,3)-2*(Dij2345R(3,1)-s45*E
     -   ijR(4,2)-EijR(11,2)+EijR(22,3)-3*EijR(24,3)+EijR(39,4)-EijR(
     -   44,4))+EijR(10,3)*P(8)+(-EijR(2,1)+EijR(4,1))*P(9)-EijR(5,2)
     -   *P(18)+EijR(7,2)*P(36)-EijR(10,2)*P(43)-EijR(8,2)*P(52)+EijR
     -   (16,3)*P(93)+EijR(20,3)*P(139)-EijR(2,2)*P(210)+(-EijR(14,3)
     -   +EijR(18,3))*P(316)-EijR(15,3)*P(317)+EijR(9,2)*P(320))
       FI(38) = 8*(Dij2345I(1,1)+p4sq*EijI(4,3)-2*(Dij2345I(3,1)-s45*E
     -   ijI(4,2)-EijI(11,2)+EijI(22,3)-3*EijI(24,3)+EijI(39,4)-EijI(
     -   44,4))+EijI(10,3)*P(8)+(-EijI(2,1)+EijI(4,1))*P(9)-EijI(5,2)
     -   *P(18)+EijI(7,2)*P(36)-EijI(10,2)*P(43)-EijI(8,2)*P(52)+EijI
     -   (16,3)*P(93)+EijI(20,3)*P(139)-EijI(2,2)*P(210)+(-EijI(14,3)
     -   +EijI(18,3))*P(316)-EijI(15,3)*P(317)+EijI(9,2)*P(320))
       F(38)=DCMPLX(FR(38),FI(38))
       P(321) = p2sq+p3sq-p4sq-s12-s23+s45
       FR(39) = -8*(Dij2345R(3,1)-p4sq*EijR(4,3)+p2sq*EijR(8,2)+2*(Eij
     -   R(23,3)-3*EijR(24,3)+EijR(39,4)-EijR(45,4))-EijR(3,2)*P(1)-E
     -   ijR(4,1)*P(5)+(EijR(15,3)-EijR(20,3))*P(8)-EijR(3,1)*P(9)+Ei
     -   jR(4,2)*P(15)-EijR(6,2)*P(18)-EijR(13,3)*P(139)+(EijR(14,3)-
     -   EijR(19,3))*P(316)-EijR(16,3)*P(318)+EijR(10,2)*P(321))
       FI(39) = -8*(Dij2345I(3,1)-p4sq*EijI(4,3)+p2sq*EijI(8,2)+2*(Eij
     -   I(23,3)-3*EijI(24,3)+EijI(39,4)-EijI(45,4))-EijI(3,2)*P(1)-E
     -   ijI(4,1)*P(5)+(EijI(15,3)-EijI(20,3))*P(8)-EijI(3,1)*P(9)+Ei
     -   jI(4,2)*P(15)-EijI(6,2)*P(18)-EijI(13,3)*P(139)+(EijI(14,3)-
     -   EijI(19,3))*P(316)-EijI(16,3)*P(318)+EijI(10,2)*P(321))
       F(39)=DCMPLX(FR(39),FI(39))
       P(322) = p4sq-2*s45
       P(323) = 2*p3sq+s15-s23-s45
       P(324) = -p4sq+s34+2*P(78)
       P(325) = s15-s45
       P(326) = p2sq+p3sq-2*s12-s15-s23+s34
       FR(40) = 8*(Dij2345R(3,2)-s45*EijR(4,3)+s12*(EijR(2,1)+EijR(5,2
     -   ))+s15*EijR(7,3)-2*(EijR(21,3)-EijR(24,3)+EijR(39,4)-EijR(42
     -   ,4))+(EijR(3,1)+EijR(6,2))*P(45)+(EijR(16,3)-EijR(19,3))*P(8
     -   3)-EijR(4,1)*P(214)+EijR(9,2)*P(222)-EijR(4,2)*P(322)+EijR(7
     -   ,2)*P(323)-EijR(10,2)*P(324)-EijR(14,3)*P(325)+(-EijR(15,3)+
     -   EijR(18,3))*P(326))
       FI(40) = 8*(Dij2345I(3,2)-s45*EijI(4,3)+s12*(EijI(2,1)+EijI(5,2
     -   ))+s15*EijI(7,3)-2*(EijI(21,3)-EijI(24,3)+EijI(39,4)-EijI(42
     -   ,4))+(EijI(3,1)+EijI(6,2))*P(45)+(EijI(16,3)-EijI(19,3))*P(8
     -   3)-EijI(4,1)*P(214)+EijI(9,2)*P(222)-EijI(4,2)*P(322)+EijI(7
     -   ,2)*P(323)-EijI(10,2)*P(324)-EijI(14,3)*P(325)+(-EijI(15,3)+
     -   EijI(18,3))*P(326))
       F(40)=DCMPLX(FR(40),FI(40))
       P(327) = p2sq+s15-s23
       P(328) = -s45+2*P(175)
       P(329) = p3sq+p4sq-s34+2*P(52)
       P(330) = p4sq+s15-s23+s45
       FR(41) = 8*(D01345R-Dij2345R(3,1)+p4sq*EijR(4,3)+p3sq*(EijR(3,1
     -   )+EijR(6,2))-2*(EijR(11,2)+EijR(21,3)-2*EijR(24,3)+EijR(39,4
     -   )-EijR(42,4))-(EijR(2,1)+EijR(5,2))*P(1)+EijR(18,3)*P(8)-Eij
     -   R(4,2)*P(15)+(EijR(16,3)-EijR(19,3))*P(93)+Is12*(C0134R-p4sq
     -   *(Dij1345R(3,1)+Dij1345R(3,2))-2*(Dij1345R(7,2)-Dij1345R(13,
     -   3))-Dij1345R(5,2)*P(112)+Dij1345R(2,1)*P(128)+Dij1345R(1,1)*
     -   P(136)+Dij1345R(6,2)*P(203))+EijR(7,3)*P(316)+EijR(15,3)*P(3
     -   19)+EijR(9,2)*P(326)+EijR(4,1)*P(327)-EijR(7,2)*P(328)+EijR(
     -   10,2)*P(329)-EijR(14,3)*P(330))
       FI(41) = 8*(D01345I-Dij2345I(3,1)+p4sq*EijI(4,3)+p3sq*(EijI(3,1
     -   )+EijI(6,2))-2*(EijI(11,2)+EijI(21,3)-2*EijI(24,3)+EijI(39,4
     -   )-EijI(42,4))-(EijI(2,1)+EijI(5,2))*P(1)+EijI(18,3)*P(8)-Eij
     -   I(4,2)*P(15)+(EijI(16,3)-EijI(19,3))*P(93)+Is12*(C0134I-p4sq
     -   *(Dij1345I(3,1)+Dij1345I(3,2))-2*(Dij1345I(7,2)-Dij1345I(13,
     -   3))-Dij1345I(5,2)*P(112)+Dij1345I(2,1)*P(128)+Dij1345I(1,1)*
     -   P(136)+Dij1345I(6,2)*P(203))+EijI(7,3)*P(316)+EijI(15,3)*P(3
     -   19)+EijI(9,2)*P(326)+EijI(4,1)*P(327)-EijI(7,2)*P(328)+EijI(
     -   10,2)*P(329)-EijI(14,3)*P(330))
       F(41)=DCMPLX(FR(41),FI(41))
       P(331) = s12+s23+s45
       P(332) = p3sq**2+s23*s45+p2sq*P(30)-p3sq*P(331)
       P(333) = -p3sq-3*p4sq+s34
       P(334) = -3*p4sq+s12+s23+s34
       P(335) = -p3sq**2-s34*P(177)+p4sq*P(292)+p2sq*P(333)+p3sq*P(334
     -   )
       P(336) = -(p2sq*P(144))+s12*P(161)
       P(337) = 2*p2sq**2-p2sq*P(176)-s12*P(178)
       P(338) = -2*s12+3*P(1)
       P(339) = 8*p2sq+5*p3sq-p4sq-s15-6*P(177)+2*P(200)
       FR(42) = Is45*(-8*Dij1234R(13,3)*P(78)+4*((Cij134R(1,1)+2*(Dij1
     -   234R(7,2)+Dij1234R(12,3)))*P(9)-(D01234R+Dij1234R(1,1))*P(15
     -   1)))+Is12s45*(2*P(9)*(3*B013R+B014R+2*(p2sq*Cij123R(2,1)-(C0
     -   123R+Cij123R(1,1))*P(18)-C0134R*P(43)))+4*(P(5)*(p4sq*Cij145
     -   R(2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2))+C0145R*P(1
     -   5))-Cij134R(2,1)*P(332)))+4*(-4*EijR(46,4)+EE0R*s15*P(5)-Dij
     -   2345R(2,1)*P(53)+D01345R*P(149)-D02345R*P(152)-s15*(EijR(1,2
     -   )*P(18)+EijR(1,1)*P(159))+EijR(3,1)*P(166)+EijR(4,1)*P(173)-
     -   EijR(7,2)*P(186)+P(155)*(Dij2345R(1,1)-Dij2345R(3,1)+EijR(2,
     -   2)*P(68)+EijR(4,2)*P(174)+EijR(9,2)*P(192))-EijR(10,2)*P(198
     -   )-Dij1345R(1,2)*P(261)-EijR(3,2)*P(302)+EijR(6,2)*P(310)-Eij
     -   R(8,2)*P(311)-EijR(2,1)*P(336)-EijR(5,2)*P(337)+Is12*((Cij14
     -   5R(1,1)+2*p4sq*Dij1345R(3,2))*P(9)-Dij1345R(1,1)*P(283)-Dij1
     -   345R(2,1)*P(285)+Dij1345R(3,1)*P(286)+Dij1345R(2,2)*P(1)*P(2
     -   87)-Dij1345R(4,2)*P(289)+Dij1345R(5,2)*P(291)+Dij1345R(6,2)*
     -   P(335)+2*Dij1345R(7,2)*P(338))+2*(p4sq*EijR(24,3)+EijR(23,3)
     -   *P(93)-EijR(21,3)*P(316)+EijR(22,3)*P(319)+EijR(11,2)*P(339)
     -   ))
       FI(42) = Is45*(-8*Dij1234I(13,3)*P(78)+4*((Cij134I(1,1)+2*(Dij1
     -   234I(7,2)+Dij1234I(12,3)))*P(9)-(D01234I+Dij1234I(1,1))*P(15
     -   1)))+Is12s45*(2*P(9)*(3*B013I+B014I+2*(p2sq*Cij123I(2,1)-(C0
     -   123I+Cij123I(1,1))*P(18)-C0134I*P(43)))+4*(P(5)*(p4sq*Cij145
     -   I(2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2))+C0145I*P(1
     -   5))-Cij134I(2,1)*P(332)))+4*(-4*EijI(46,4)+EE0I*s15*P(5)-Dij
     -   2345I(2,1)*P(53)+D01345I*P(149)-D02345I*P(152)-s15*(EijI(1,2
     -   )*P(18)+EijI(1,1)*P(159))+EijI(3,1)*P(166)+EijI(4,1)*P(173)-
     -   EijI(7,2)*P(186)+P(155)*(Dij2345I(1,1)-Dij2345I(3,1)+EijI(2,
     -   2)*P(68)+EijI(4,2)*P(174)+EijI(9,2)*P(192))-EijI(10,2)*P(198
     -   )-Dij1345I(1,2)*P(261)-EijI(3,2)*P(302)+EijI(6,2)*P(310)-Eij
     -   I(8,2)*P(311)-EijI(2,1)*P(336)-EijI(5,2)*P(337)+Is12*((Cij14
     -   5I(1,1)+2*p4sq*Dij1345I(3,2))*P(9)-Dij1345I(1,1)*P(283)-Dij1
     -   345I(2,1)*P(285)+Dij1345I(3,1)*P(286)+Dij1345I(2,2)*P(1)*P(2
     -   87)-Dij1345I(4,2)*P(289)+Dij1345I(5,2)*P(291)+Dij1345I(6,2)*
     -   P(335)+2*Dij1345I(7,2)*P(338))+2*(p4sq*EijI(24,3)+EijI(23,3)
     -   *P(93)-EijI(21,3)*P(316)+EijI(22,3)*P(319)+EijI(11,2)*P(339)
     -   ))
       F(42)=DCMPLX(FR(42),FI(42))
       P(340) = 1-Is12*s45
       P(341) = -2+Is12*s34
       P(342) = -1+Is12*s45
       P(343) = 3*p2sq+2*P(210)
       P(344) = 3*P(10)+2*P(43)
       P(345) = p4sq-s34+2*P(52)
       P(346) = p4sq+s15-2*s45
       P(347) = p4sq-3*s45
       P(348) = s15-s34
       P(349) = p2sq+p3sq-4*s12-s23-2*P(348)
       P(350) = p2sq+p3sq-s23-4*P(52)-2*P(128)
       P(351) = p4sq+3*P(325)
       P(352) = 2*s15-s45
       P(353) = s15-2*s45
       FR(43) = 4*Is12s45*(3*B013R+B014R-2*(p2sq*(C0123R+Cij123R(1,1))
     -   +p4sq*Cij145R(2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2)
     -   )+C0145R*P(15)-Cij123R(2,1)*P(18)+Cij134R(2,1)*P(43)-C0134R*
     -   P(298)))+8*(Dij1345R(1,2)+Dij2345R(3,2)-s15*(EE0R+3*(EijR(1,
     -   1)+EijR(1,2))+EijR(1,3))-s45*EijR(4,3)+2*(EijR(11,2)-EijR(36
     -   ,4)-EijR(39,4))+4*EijR(42,4)+D01345R*P(39)-Cij134R(3,2)*P(44
     -   )+(EijR(6,3)+EijR(16,3)-2*EijR(19,3))*P(83)+Is12*(Cij145R(1,
     -   1)-s45*Dij1345R(3,1)+p4sq*Dij1345R(3,2)-2*(Dij1345R(11,3)-Di
     -   j1345R(13,3))-Dij1345R(5,2)*P(137))+EijR(2,1)*P(265)-(EijR(5
     -   ,3)+EijR(15,3)-2*EijR(18,3))*P(326)+Dij1345R(6,2)*P(340)-Dij
     -   1345R(1,1)*P(341)+(Dij1345R(2,1)+Dij1345R(4,2))*P(342)-Is45*
     -   (Cij123R(3,2)+Cij134R(1,2)+s23*(D01234R+2*Dij1234R(1,1)+Dij1
     -   234R(1,2))+4*(Dij1234R(7,2)+Dij1234R(11,3))-Dij1234R(4,3)*P(
     -   9)+Dij1234R(3,1)*P(48)+Dij1234R(5,3)*P(78)-Dij1234R(2,1)*P(1
     -   24)-Dij1234R(4,2)*P(343)+Dij1234R(5,2)*P(344))-EijR(3,1)*P(3
     -   45)+EijR(4,1)*P(346)-EijR(4,2)*P(347)+(-EijR(5,2)+EijR(9,2))
     -   *P(349)+(EijR(6,2)-EijR(10,2))*P(350)+EijR(7,2)*P(351)+EijR(
     -   7,3)*P(352)-EijR(14,3)*P(353))
       FI(43) = 4*Is12s45*(3*B013I+B014I-2*(p2sq*(C0123I+Cij123I(1,1))
     -   +p4sq*Cij145I(2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2)
     -   )+C0145I*P(15)-Cij123I(2,1)*P(18)+Cij134I(2,1)*P(43)-C0134I*
     -   P(298)))+8*(Dij1345I(1,2)+Dij2345I(3,2)-s15*(EE0I+3*(EijI(1,
     -   1)+EijI(1,2))+EijI(1,3))-s45*EijI(4,3)+2*(EijI(11,2)-EijI(36
     -   ,4)-EijI(39,4))+4*EijI(42,4)+D01345I*P(39)-Cij134I(3,2)*P(44
     -   )+(EijI(6,3)+EijI(16,3)-2*EijI(19,3))*P(83)+Is12*(Cij145I(1,
     -   1)-s45*Dij1345I(3,1)+p4sq*Dij1345I(3,2)-2*(Dij1345I(11,3)-Di
     -   j1345I(13,3))-Dij1345I(5,2)*P(137))+EijI(2,1)*P(265)-(EijI(5
     -   ,3)+EijI(15,3)-2*EijI(18,3))*P(326)+Dij1345I(6,2)*P(340)-Dij
     -   1345I(1,1)*P(341)+(Dij1345I(2,1)+Dij1345I(4,2))*P(342)-Is45*
     -   (Cij123I(3,2)+Cij134I(1,2)+s23*(D01234I+2*Dij1234I(1,1)+Dij1
     -   234I(1,2))+4*(Dij1234I(7,2)+Dij1234I(11,3))-Dij1234I(4,3)*P(
     -   9)+Dij1234I(3,1)*P(48)+Dij1234I(5,3)*P(78)-Dij1234I(2,1)*P(1
     -   24)-Dij1234I(4,2)*P(343)+Dij1234I(5,2)*P(344))-EijI(3,1)*P(3
     -   45)+EijI(4,1)*P(346)-EijI(4,2)*P(347)+(-EijI(5,2)+EijI(9,2))
     -   *P(349)+(EijI(6,2)-EijI(10,2))*P(350)+EijI(7,2)*P(351)+EijI(
     -   7,3)*P(352)-EijI(14,3)*P(353))
       F(43)=DCMPLX(FR(43),FI(43))
       P(354) = p2sq-s12+2*P(6)
       P(355) = p4sq-s12-s45
       P(356) = 2*p4sq+s12+s15-s34-3*s45
       P(357) = p2sq+p4sq-s12+s15-s23
       P(358) = -2*p3sq+3*s15
       P(359) = p4sq-s34+2*P(25)
       P(360) = p3sq+p4sq-s15-s45
       P(361) = p2sq-3*s12-s23+s34+2*P(360)
       P(362) = -p2sq+s23+4*P(52)+2*P(128)
       P(363) = p2sq+p4sq-s12-s23+s34+s45
       FR(44) = 4*Is12s45*(3*B013R+B014R+2*(p2sq*Cij123R(2,1)-p4sq*Cij
     -   145R(2,1)-2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2))-C0145R*
     -   P(15)-(C0123R+Cij123R(1,1))*P(18)-Cij134R(2,1)*P(43)+C0134R*
     -   P(298)))+8*(Dij1345R(1,2)+Dij2345R(1,1)-s15*(EE0R+2*EijR(1,1
     -   )+EijR(1,2))+p4sq*EijR(4,3)-2*(Dij2345R(3,1)-EijR(11,2)+2*Ei
     -   jR(21,3)+EijR(22,3)-3*EijR(24,3)+EijR(39,4)+EijR(40,4)-EijR(
     -   42,4)-EijR(44,4))+(-EijR(8,3)+EijR(10,3))*P(8)+D01345R*P(39)
     -   -Cij134R(3,2)*P(44)+EijR(6,2)*P(87)+(EijR(16,3)+EijR(17,3))*
     -   P(93)+Is12*(Cij145R(1,1)-s45*Dij1345R(3,1)+p4sq*Dij1345R(3,2
     -   )-2*(Dij1345R(11,3)-Dij1345R(13,3))-Dij1345R(5,2)*P(137))+(E
     -   ijR(19,3)+EijR(20,3))*P(139)+EijR(2,2)*P(265)+(-EijR(5,3)+Ei
     -   jR(7,3))*P(316)-EijR(15,3)*P(317)-EijR(4,2)*P(322)-EijR(14,3
     -   )*P(330)+Dij1345R(6,2)*P(340)-Dij1345R(1,1)*P(341)+(Dij1345R
     -   (2,1)+Dij1345R(4,2))*P(342)-EijR(3,1)*P(345)-Is45*(Cij134R(1
     -   ,2)+s23*(D01234R+Dij1234R(1,1))+2*Dij1234R(12,3)-Dij1234R(6,
     -   3)*P(9)+(Dij1234R(5,2)+Dij1234R(10,3))*P(78)+(Dij1234R(3,1)+
     -   Dij1234R(6,2))*P(80)-Dij1234R(2,2)*P(124)+Dij1234R(2,1)*P(27
     -   6)-Dij1234R(4,2)*P(354))-EijR(2,1)*P(355)+EijR(4,1)*P(356)-E
     -   ijR(5,2)*P(357)+EijR(7,2)*P(358)-EijR(8,2)*P(359)+EijR(9,2)*
     -   P(361)+EijR(10,2)*P(362)+EijR(18,3)*P(363))
       FI(44) = 4*Is12s45*(3*B013I+B014I+2*(p2sq*Cij123I(2,1)-p4sq*Cij
     -   145I(2,1)-2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2))-C0145I*
     -   P(15)-(C0123I+Cij123I(1,1))*P(18)-Cij134I(2,1)*P(43)+C0134I*
     -   P(298)))+8*(Dij1345I(1,2)+Dij2345I(1,1)-s15*(EE0I+2*EijI(1,1
     -   )+EijI(1,2))+p4sq*EijI(4,3)-2*(Dij2345I(3,1)-EijI(11,2)+2*Ei
     -   jI(21,3)+EijI(22,3)-3*EijI(24,3)+EijI(39,4)+EijI(40,4)-EijI(
     -   42,4)-EijI(44,4))+(-EijI(8,3)+EijI(10,3))*P(8)+D01345I*P(39)
     -   -Cij134I(3,2)*P(44)+EijI(6,2)*P(87)+(EijI(16,3)+EijI(17,3))*
     -   P(93)+Is12*(Cij145I(1,1)-s45*Dij1345I(3,1)+p4sq*Dij1345I(3,2
     -   )-2*(Dij1345I(11,3)-Dij1345I(13,3))-Dij1345I(5,2)*P(137))+(E
     -   ijI(19,3)+EijI(20,3))*P(139)+EijI(2,2)*P(265)+(-EijI(5,3)+Ei
     -   jI(7,3))*P(316)-EijI(15,3)*P(317)-EijI(4,2)*P(322)-EijI(14,3
     -   )*P(330)+Dij1345I(6,2)*P(340)-Dij1345I(1,1)*P(341)+(Dij1345I
     -   (2,1)+Dij1345I(4,2))*P(342)-EijI(3,1)*P(345)-Is45*(Cij134I(1
     -   ,2)+s23*(D01234I+Dij1234I(1,1))+2*Dij1234I(12,3)-Dij1234I(6,
     -   3)*P(9)+(Dij1234I(5,2)+Dij1234I(10,3))*P(78)+(Dij1234I(3,1)+
     -   Dij1234I(6,2))*P(80)-Dij1234I(2,2)*P(124)+Dij1234I(2,1)*P(27
     -   6)-Dij1234I(4,2)*P(354))-EijI(2,1)*P(355)+EijI(4,1)*P(356)-E
     -   ijI(5,2)*P(357)+EijI(7,2)*P(358)-EijI(8,2)*P(359)+EijI(9,2)*
     -   P(361)+EijI(10,2)*P(362)+EijI(18,3)*P(363))
       F(44)=DCMPLX(FR(44),FI(44))
       P(364) = p3sq-s12+s15-s23+s34
       P(365) = 2*p3sq+p4sq+s15-s23-s45
       P(366) = p2sq+2*p3sq-p4sq-s12-s23+s34+s45
       P(367) = p2sq-p3sq+p4sq-s15-s45-2*P(230)
       P(368) = -p2sq-p3sq+2*s12+s15+s23-s34+s45
       P(369) = p2sq+p3sq-2*s12-s23+s34-s45
       FR(45) = -8*(Dij2345R(1,1)-Dij2345R(3,2)+Dij2345R(5,2)+s45*EijR
     -   (4,3)+p3sq*EijR(6,2)+s15*(EijR(5,3)-EijR(7,3))+EijR(2,2)*P(8
     -   )+2*(2*EijR(11,2)+EijR(21,3)-EijR(24,3)+EijR(39,4)+EijR(40,4
     -   )-EijR(42,4)-EijR(44,4)+EijR(4,2)*P(15))+Is45*(Cij123R(2,1)+
     -   Cij123R(2,2)+2*(Dij1234R(7,2)+Dij1234R(12,3))-(Dij1234R(2,2)
     -   +Dij1234R(6,3))*P(9)+(Dij1234R(3,1)+Dij1234R(6,2))*P(45)-(Di
     -   j1234R(2,1)+Dij1234R(4,2))*P(66)+(Dij1234R(5,2)+Dij1234R(10,
     -   3))*P(78))+(-EijR(16,3)-EijR(17,3)+EijR(19,3)+EijR(20,3))*P(
     -   83)-EijR(4,1)*P(153)-EijR(2,1)*P(230)+EijR(10,2)*P(257)+EijR
     -   (14,3)*P(325)+(EijR(8,3)-EijR(10,3))*P(326)+EijR(5,2)*P(364)
     -   -EijR(7,2)*P(365)-EijR(8,2)*P(366)-EijR(9,2)*P(367)-EijR(15,
     -   3)*P(368)-EijR(18,3)*P(369))
       FI(45) = -8*(Dij2345I(1,1)-Dij2345I(3,2)+Dij2345I(5,2)+s45*EijI
     -   (4,3)+p3sq*EijI(6,2)+s15*(EijI(5,3)-EijI(7,3))+EijI(2,2)*P(8
     -   )+2*(2*EijI(11,2)+EijI(21,3)-EijI(24,3)+EijI(39,4)+EijI(40,4
     -   )-EijI(42,4)-EijI(44,4)+EijI(4,2)*P(15))+Is45*(Cij123I(2,1)+
     -   Cij123I(2,2)+2*(Dij1234I(7,2)+Dij1234I(12,3))-(Dij1234I(2,2)
     -   +Dij1234I(6,3))*P(9)+(Dij1234I(3,1)+Dij1234I(6,2))*P(45)-(Di
     -   j1234I(2,1)+Dij1234I(4,2))*P(66)+(Dij1234I(5,2)+Dij1234I(10,
     -   3))*P(78))+(-EijI(16,3)-EijI(17,3)+EijI(19,3)+EijI(20,3))*P(
     -   83)-EijI(4,1)*P(153)-EijI(2,1)*P(230)+EijI(10,2)*P(257)+EijI
     -   (14,3)*P(325)+(EijI(8,3)-EijI(10,3))*P(326)+EijI(5,2)*P(364)
     -   -EijI(7,2)*P(365)-EijI(8,2)*P(366)-EijI(9,2)*P(367)-EijI(15,
     -   3)*P(368)-EijI(18,3)*P(369))
       F(45)=DCMPLX(FR(45),FI(45))
       P(370) = p2sq+p4sq-s12-s15+s34-s45
       P(371) = -2*p4sq+s45
       P(372) = p2sq+3*p4sq-s12-s15+s34-2*s45
       P(373) = p4sq+2*P(8)
       P(374) = p2sq+2*p4sq-s12-s15+s34
       FR(46) = 8*(D02345R+Dij2345R(1,1)-Dij2345R(3,1)+p4sq*EijR(4,3)-
     -   2*(4*(EijR(22,3)-EijR(24,3))+EijR(37,4)+EijR(39,4)+2*(EijR(1
     -   1,2)-EijR(44,4)))-EijR(2,3)*P(8)+(-EijR(2,1)+EijR(4,1))*P(15
     -   )+Is45*((Dij1234R(2,1)+2*Dij1234R(2,2)+Dij1234R(2,3))*P(9)-(
     -   Dij1234R(3,1)+2*Dij1234R(6,2)+Dij1234R(8,3))*P(78))+(EijR(8,
     -   2)+EijR(9,3)+EijR(16,3))*P(93)+(EijR(10,2)+2*EijR(20,3))*P(1
     -   39)-(EijR(5,2)-EijR(7,2)+EijR(8,3)+EijR(14,3)-2*EijR(18,3))*
     -   P(316)-EijR(2,2)*P(370)+EijR(4,2)*P(371)+EijR(9,2)*P(372)+Ei
     -   jR(10,3)*P(373)-EijR(15,3)*P(374))
       FI(46) = 8*(D02345I+Dij2345I(1,1)-Dij2345I(3,1)+p4sq*EijI(4,3)-
     -   2*(4*(EijI(22,3)-EijI(24,3))+EijI(37,4)+EijI(39,4)+2*(EijI(1
     -   1,2)-EijI(44,4)))-EijI(2,3)*P(8)+(-EijI(2,1)+EijI(4,1))*P(15
     -   )+Is45*((Dij1234I(2,1)+2*Dij1234I(2,2)+Dij1234I(2,3))*P(9)-(
     -   Dij1234I(3,1)+2*Dij1234I(6,2)+Dij1234I(8,3))*P(78))+(EijI(8,
     -   2)+EijI(9,3)+EijI(16,3))*P(93)+(EijI(10,2)+2*EijI(20,3))*P(1
     -   39)-(EijI(5,2)-EijI(7,2)+EijI(8,3)+EijI(14,3)-2*EijI(18,3))*
     -   P(316)-EijI(2,2)*P(370)+EijI(4,2)*P(371)+EijI(9,2)*P(372)+Ei
     -   jI(10,3)*P(373)-EijI(15,3)*P(374))
       F(46)=DCMPLX(FR(46),FI(46))
       P(375) = 2*p2sq+p3sq-s12-s23+s45
       P(376) = s12-s45-2*P(1)
       P(377) = p2sq-s12+s45
       P(378) = p3sq+2*p4sq+s12-s34-s45
       P(379) = p4sq-s34+s45-2*P(175)
       P(380) = p4sq+s45-2*P(175)
       P(381) = p3sq+s12-s34
       P(382) = 3*P(15)+2*P(381)
       P(383) = p2sq-2*p4sq-s12-s15+s34
       FR(47) = 8*(D02345R+Dij2345R(2,1)-Dij2345R(3,1)+p4sq*EijR(4,3)-
     -   s12*EijR(5,2)-2*(2*EijR(21,3)+EijR(23,3)+3*(EijR(11,2)-EijR(
     -   24,3))+EijR(39,4)+EijR(41,4)-EijR(42,4)-EijR(45,4))+(-EijR(1
     -   5,3)-EijR(17,3)+EijR(18,3)+EijR(20,3))*P(8)-EijR(9,2)*P(28)+
     -   Is12s45*(Cij134R(2,1)+Cij134R(2,2))*P(52)+EijR(11,3)*P(93)+(
     -   -EijR(3,1)+EijR(4,1))*P(99)+EijR(13,3)*P(139)+EijR(8,2)*P(26
     -   5)+(-EijR(6,3)+EijR(7,3))*P(316)+EijR(16,3)*P(318)-EijR(14,3
     -   )*P(330)-EijR(3,2)*P(359)+EijR(4,2)*P(371)-Is45*(Cij134R(1,1
     -   )+Cij134R(3,2)+2*(Dij1234R(7,2)+Dij1234R(13,3))-(Dij1234R(4,
     -   2)+Dij1234R(10,3))*P(9)+Dij1234R(7,3)*P(78)-(Dij1234R(2,1)+D
     -   ij1234R(6,2))*P(124)+Dij1234R(3,1)*P(375)-Dij1234R(3,2)*P(37
     -   6)+Dij1234R(5,2)*P(377))+Is12*(C0134R+s34*Dij1345R(1,1)-s34*
     -   Dij1345R(2,1)-p4sq*Dij1345R(3,2)-4*Dij1345R(7,2)-2*(Dij1345R
     -   (12,3)-Dij1345R(13,3))+(Dij1345R(4,2)-Dij1345R(5,2))*P(112)-
     -   Dij1345R(2,2)*P(203)+Dij1345R(6,2)*P(378))-EijR(6,2)*P(379)+
     -   EijR(7,2)*P(380)+EijR(10,2)*P(382)-EijR(19,3)*P(383))
       FI(47) = 8*(D02345I+Dij2345I(2,1)-Dij2345I(3,1)+p4sq*EijI(4,3)-
     -   s12*EijI(5,2)-2*(2*EijI(21,3)+EijI(23,3)+3*(EijI(11,2)-EijI(
     -   24,3))+EijI(39,4)+EijI(41,4)-EijI(42,4)-EijI(45,4))+(-EijI(1
     -   5,3)-EijI(17,3)+EijI(18,3)+EijI(20,3))*P(8)-EijI(9,2)*P(28)+
     -   Is12s45*(Cij134I(2,1)+Cij134I(2,2))*P(52)+EijI(11,3)*P(93)+(
     -   -EijI(3,1)+EijI(4,1))*P(99)+EijI(13,3)*P(139)+EijI(8,2)*P(26
     -   5)+(-EijI(6,3)+EijI(7,3))*P(316)+EijI(16,3)*P(318)-EijI(14,3
     -   )*P(330)-EijI(3,2)*P(359)+EijI(4,2)*P(371)-Is45*(Cij134I(1,1
     -   )+Cij134I(3,2)+2*(Dij1234I(7,2)+Dij1234I(13,3))-(Dij1234I(4,
     -   2)+Dij1234I(10,3))*P(9)+Dij1234I(7,3)*P(78)-(Dij1234I(2,1)+D
     -   ij1234I(6,2))*P(124)+Dij1234I(3,1)*P(375)-Dij1234I(3,2)*P(37
     -   6)+Dij1234I(5,2)*P(377))+Is12*(C0134I+s34*Dij1345I(1,1)-s34*
     -   Dij1345I(2,1)-p4sq*Dij1345I(3,2)-4*Dij1345I(7,2)-2*(Dij1345I
     -   (12,3)-Dij1345I(13,3))+(Dij1345I(4,2)-Dij1345I(5,2))*P(112)-
     -   Dij1345I(2,2)*P(203)+Dij1345I(6,2)*P(378))-EijI(6,2)*P(379)+
     -   EijI(7,2)*P(380)+EijI(10,2)*P(382)-EijI(19,3)*P(383))
       F(47)=DCMPLX(FR(47),FI(47))
       P(384) = 2*p3sq+s15-s23+s34-s45
       P(385) = p3sq-p4sq
       P(386) = p2sq-s12-s23+3*s45+2*P(385)
       P(387) = p2sq+p3sq-p4sq-2*s12-s23+s34+3*s45
       P(388) = -p2sq-p3sq+p4sq+2*s12+s15+s23-s34-3*s45
       FR(48) = -8*(D02345R+Dij2345R(2,1)-Dij2345R(3,2)+Dij2345R(6,2)+
     -   s45*EijR(4,3)-s12*EijR(5,2)+s15*(EijR(6,3)-EijR(7,3))-2*(2*E
     -   ijR(11,2)-EijR(39,4)-EijR(41,4)+EijR(42,4)+EijR(45,4))+EijR(
     -   8,2)*P(8)+Is45*(C0134R+p2sq*Dij1234R(2,1)+p3sq*Dij1234R(5,2)
     -   -2*(Dij1234R(7,2)-Dij1234R(13,3))-(Dij1234R(6,2)+Dij1234R(10
     -   ,3))*P(9)-Dij1234R(3,1)*P(10)+Dij1234R(3,2)*P(45)+Dij1234R(7
     -   ,3)*P(78))+(-EijR(11,3)+EijR(13,3))*P(83)-EijR(4,1)*P(153)-E
     -   ijR(9,2)*P(222)-EijR(3,1)*P(230)+EijR(4,2)*P(322)-EijR(7,2)*
     -   P(323)+EijR(14,3)*P(325)+(EijR(15,3)+EijR(17,3)-EijR(18,3)-E
     -   ijR(20,3))*P(326)-EijR(3,2)*P(366)+EijR(6,2)*P(384)+EijR(10,
     -   2)*P(386)-EijR(16,3)*P(387)-EijR(19,3)*P(388))
       FI(48) = -8*(D02345I+Dij2345I(2,1)-Dij2345I(3,2)+Dij2345I(6,2)+
     -   s45*EijI(4,3)-s12*EijI(5,2)+s15*(EijI(6,3)-EijI(7,3))-2*(2*E
     -   ijI(11,2)-EijI(39,4)-EijI(41,4)+EijI(42,4)+EijI(45,4))+EijI(
     -   8,2)*P(8)+Is45*(C0134I+p2sq*Dij1234I(2,1)+p3sq*Dij1234I(5,2)
     -   -2*(Dij1234I(7,2)-Dij1234I(13,3))-(Dij1234I(6,2)+Dij1234I(10
     -   ,3))*P(9)-Dij1234I(3,1)*P(10)+Dij1234I(3,2)*P(45)+Dij1234I(7
     -   ,3)*P(78))+(-EijI(11,3)+EijI(13,3))*P(83)-EijI(4,1)*P(153)-E
     -   ijI(9,2)*P(222)-EijI(3,1)*P(230)+EijI(4,2)*P(322)-EijI(7,2)*
     -   P(323)+EijI(14,3)*P(325)+(EijI(15,3)+EijI(17,3)-EijI(18,3)-E
     -   ijI(20,3))*P(326)-EijI(3,2)*P(366)+EijI(6,2)*P(384)+EijI(10,
     -   2)*P(386)-EijI(16,3)*P(387)-EijI(19,3)*P(388))
       F(48)=DCMPLX(FR(48),FI(48))
       P(389) = p2sq+p3sq+p4sq-s23-s45
       P(390) = p3sq-p4sq-s12+2*s45
       P(391) = 2*p4sq-s15+s23-s45
       FR(49) = 8*(-Dij1234R(6,2)+Dij2345R(1,1)+Dij2345R(2,1)+p2sq*Eij
     -   R(2,2)+p3sq*EijR(3,2)+p4sq*EijR(4,3)-2*(Dij2345R(3,1)-2*(Eij
     -   R(11,2)-EijR(22,3)-EijR(23,3))-4*EijR(24,3)+EijR(39,4)+EijR(
     -   43,4)-EijR(44,4)-EijR(45,4))+(-EijR(9,3)+EijR(10,3))*P(8)-Ei
     -   jR(5,2)*P(18)+EijR(7,2)*P(36)+EijR(6,2)*P(45)+Is45*((Dij1234
     -   R(2,1)+Dij1234R(2,2)+Dij1234R(8,3))*P(9)-(Dij1234R(3,1)+Dij1
     -   234R(3,2)+Dij1234R(9,3))*P(78))+EijR(12,3)*P(93)+EijR(13,3)*
     -   P(139)-EijR(9,2)*P(161)+(-EijR(14,3)-EijR(17,3)+EijR(18,3)+E
     -   ijR(19,3))*P(316)-EijR(15,3)*P(317)+EijR(16,3)*P(318)-EijR(4
     -   ,2)*P(322)-EijR(8,2)*P(389)-EijR(10,2)*P(390)+EijR(20,3)*P(3
     -   91))
       FI(49) = 8*(-Dij1234I(6,2)+Dij2345I(1,1)+Dij2345I(2,1)+p2sq*Eij
     -   I(2,2)+p3sq*EijI(3,2)+p4sq*EijI(4,3)-2*(Dij2345I(3,1)-2*(Eij
     -   I(11,2)-EijI(22,3)-EijI(23,3))-4*EijI(24,3)+EijI(39,4)+EijI(
     -   43,4)-EijI(44,4)-EijI(45,4))+(-EijI(9,3)+EijI(10,3))*P(8)-Ei
     -   jI(5,2)*P(18)+EijI(7,2)*P(36)+EijI(6,2)*P(45)+Is45*((Dij1234
     -   I(2,1)+Dij1234I(2,2)+Dij1234I(8,3))*P(9)-(Dij1234I(3,1)+Dij1
     -   234I(3,2)+Dij1234I(9,3))*P(78))+EijI(12,3)*P(93)+EijI(13,3)*
     -   P(139)-EijI(9,2)*P(161)+(-EijI(14,3)-EijI(17,3)+EijI(18,3)+E
     -   ijI(19,3))*P(316)-EijI(15,3)*P(317)+EijI(16,3)*P(318)-EijI(4
     -   ,2)*P(322)-EijI(8,2)*P(389)-EijI(10,2)*P(390)+EijI(20,3)*P(3
     -   91))
       F(49)=DCMPLX(FR(49),FI(49))
       P(392) = p3sq-p4sq+s12+s15-s23-s34+s45
       P(393) = p2sq-p3sq-3*p4sq-s23+s34+s45
       FR(50) = 8*(D02345R-p2sq*EijR(2,2)+p4sq*(-2*EijR(4,2)+EijR(4,3)
     -   )-2*(4*EijR(11,2)+EijR(22,3)+2*EijR(23,3)-3*EijR(24,3)+EijR(
     -   39,4)+EijR(43,4)-EijR(44,4)-EijR(45,4))+(-EijR(9,3)+EijR(10,
     -   3))*P(8)+(-EijR(3,1)+EijR(4,1))*P(15)+EijR(5,2)*P(18)+Is45*(
     -   (Dij1234R(6,2)+Dij1234R(8,3))*P(9)-(Dij1234R(3,2)+Dij1234R(9
     -   ,3))*P(78))+EijR(12,3)*P(93)-EijR(6,2)*P(103)+EijR(13,3)*P(1
     -   39)-EijR(3,2)*P(204)+EijR(7,2)*P(313)+EijR(9,2)*P(314)+(-Eij
     -   R(14,3)-EijR(17,3)+EijR(18,3)+EijR(19,3))*P(316)-EijR(15,3)*
     -   P(317)+EijR(16,3)*P(318)+EijR(20,3)*P(391)+EijR(8,2)*P(392)-
     -   EijR(10,2)*P(393))
       FI(50) = 8*(D02345I-p2sq*EijI(2,2)+p4sq*(-2*EijI(4,2)+EijI(4,3)
     -   )-2*(4*EijI(11,2)+EijI(22,3)+2*EijI(23,3)-3*EijI(24,3)+EijI(
     -   39,4)+EijI(43,4)-EijI(44,4)-EijI(45,4))+(-EijI(9,3)+EijI(10,
     -   3))*P(8)+(-EijI(3,1)+EijI(4,1))*P(15)+EijI(5,2)*P(18)+Is45*(
     -   (Dij1234I(6,2)+Dij1234I(8,3))*P(9)-(Dij1234I(3,2)+Dij1234I(9
     -   ,3))*P(78))+EijI(12,3)*P(93)-EijI(6,2)*P(103)+EijI(13,3)*P(1
     -   39)-EijI(3,2)*P(204)+EijI(7,2)*P(313)+EijI(9,2)*P(314)+(-Eij
     -   I(14,3)-EijI(17,3)+EijI(18,3)+EijI(19,3))*P(316)-EijI(15,3)*
     -   P(317)+EijI(16,3)*P(318)+EijI(20,3)*P(391)+EijI(8,2)*P(392)-
     -   EijI(10,2)*P(393))
       F(50)=DCMPLX(FR(50),FI(50))
       P(394) = p2sq-s12-s23+s34+s45
       P(395) = -3*p4sq+2*P(394)
       P(396) = p2sq-3*p4sq-s12-s23+s34+s45
       FR(51) = 8*(Dij2345R(2,1)-Dij2345R(3,1)+p4sq*EijR(4,3)-2*(3*(Ei
     -   jR(23,3)-EijR(24,3))+EijR(38,4)+EijR(39,4)-2*EijR(45,4))-(Ei
     -   jR(12,3)+EijR(15,3)-2*EijR(20,3))*P(8)-(EijR(3,2)+EijR(4,2)-
     -   2*EijR(10,2))*P(15)+Is45*((Dij1234R(6,2)+Dij1234R(9,3))*P(9)
     -   -(Dij1234R(3,2)+Dij1234R(3,3))*P(78))+EijR(3,3)*P(93)-(EijR(
     -   11,3)+EijR(14,3)-2*EijR(19,3))*P(316)-EijR(13,3)*P(395)+EijR
     -   (16,3)*P(396))
       FI(51) = 8*(Dij2345I(2,1)-Dij2345I(3,1)+p4sq*EijI(4,3)-2*(3*(Ei
     -   jI(23,3)-EijI(24,3))+EijI(38,4)+EijI(39,4)-2*EijI(45,4))-(Ei
     -   jI(12,3)+EijI(15,3)-2*EijI(20,3))*P(8)-(EijI(3,2)+EijI(4,2)-
     -   2*EijI(10,2))*P(15)+Is45*((Dij1234I(6,2)+Dij1234I(9,3))*P(9)
     -   -(Dij1234I(3,2)+Dij1234I(3,3))*P(78))+EijI(3,3)*P(93)-(EijI(
     -   11,3)+EijI(14,3)-2*EijI(19,3))*P(316)-EijI(13,3)*P(395)+EijI
     -   (16,3)*P(396))
       F(51)=DCMPLX(FR(51),FI(51))
       FR(52) = -16*(EijR(2,1)-EijR(4,1)+EijR(4,4)+EijR(5,4)-EijR(7,4)
     -   -EijR(15,4)+6*(EijR(14,3)-EijR(18,3))+3*(EijR(4,2)-EijR(4,3)
     -   +EijR(5,2)+EijR(5,3)-EijR(7,2)-EijR(7,3)-EijR(9,2)-EijR(14,4
     -   )+EijR(15,3)+EijR(19,4)-EijR(24,4)+EijR(32,4)))
       FI(52) = -16*(EijI(2,1)-EijI(4,1)+EijI(4,4)+EijI(5,4)-EijI(7,4)
     -   -EijI(15,4)+6*(EijI(14,3)-EijI(18,3))+3*(EijI(4,2)-EijI(4,3)
     -   +EijI(5,2)+EijI(5,3)-EijI(7,2)-EijI(7,3)-EijI(9,2)-EijI(14,4
     -   )+EijI(15,3)+EijI(19,4)-EijI(24,4)+EijI(32,4)))
       F(52)=DCMPLX(FR(52),FI(52))
       FR(53) = -16*(EijR(2,2)+EijR(4,2)+EijR(4,4)+EijR(17,4)+4*(EijR(
     -   15,3)-EijR(18,3))+EijR(19,4)+EijR(21,4)-2*(EijR(4,3)-EijR(8,
     -   3)+EijR(9,2)+EijR(10,3)-EijR(14,3)+EijR(14,4)+EijR(15,4)+Eij
     -   R(24,4)+EijR(27,4)-2*EijR(32,4)))
       FI(53) = -16*(EijI(2,2)+EijI(4,2)+EijI(4,4)+EijI(17,4)+4*(EijI(
     -   15,3)-EijI(18,3))+EijI(19,4)+EijI(21,4)-2*(EijI(4,3)-EijI(8,
     -   3)+EijI(9,2)+EijI(10,3)-EijI(14,3)+EijI(14,4)+EijI(15,4)+Eij
     -   I(24,4)+EijI(27,4)-2*EijI(32,4)))
       F(53)=DCMPLX(FR(53),FI(53))
       FR(54) = 16*(EijR(4,3)-EijR(4,4)+EijR(6,2)+EijR(6,3)-EijR(7,2)-
     -   EijR(7,3)-EijR(8,2)+EijR(9,2)+EijR(15,4)-EijR(16,3)+EijR(16,
     -   4)-EijR(19,4)-EijR(23,4)+EijR(24,4)+EijR(25,4)-EijR(34,4)+2*
     -   (EijR(14,4)-EijR(15,3)-EijR(17,3)+EijR(18,3)+EijR(20,3)-EijR
     -   (32,4)-EijR(33,4)+EijR(35,4)))
       FI(54) = 16*(EijI(4,3)-EijI(4,4)+EijI(6,2)+EijI(6,3)-EijI(7,2)-
     -   EijI(7,3)-EijI(8,2)+EijI(9,2)+EijI(15,4)-EijI(16,3)+EijI(16,
     -   4)-EijI(19,4)-EijI(23,4)+EijI(24,4)+EijI(25,4)-EijI(34,4)+2*
     -   (EijI(14,4)-EijI(15,3)-EijI(17,3)+EijI(18,3)+EijI(20,3)-EijI
     -   (32,4)-EijI(33,4)+EijI(35,4)))
       F(54)=DCMPLX(FR(54),FI(54))
       FR(55) = 16*(EijR(4,3)-EijR(4,4)-EijR(7,2)-EijR(7,3)+EijR(9,2)+
     -   EijR(15,4)-EijR(19,4)+EijR(24,4)+2*(EijR(14,4)-EijR(15,3)+Ei
     -   jR(18,3)-EijR(32,4)))
       FI(55) = 16*(EijI(4,3)-EijI(4,4)-EijI(7,2)-EijI(7,3)+EijI(9,2)+
     -   EijI(15,4)-EijI(19,4)+EijI(24,4)+2*(EijI(14,4)-EijI(15,3)+Ei
     -   jI(18,3)-EijI(32,4)))
       F(55)=DCMPLX(FR(55),FI(55))
       FR(56) = -16*(EijR(2,1)+EijR(2,2)-EijR(4,1)+3*(EijR(4,2)-EijR(4
     -   ,3))+EijR(4,4)+EijR(5,3)-EijR(7,3)-4*(EijR(9,2)-EijR(14,3))+
     -   5*EijR(15,3)+EijR(17,4)-6*EijR(18,3)+EijR(19,4)+EijR(21,4)+2
     -   *(EijR(5,2)-EijR(7,2)+EijR(8,3)-EijR(10,3)-EijR(14,4)-EijR(1
     -   5,4)-EijR(24,4)-EijR(27,4)+2*EijR(32,4)))
       FI(56) = -16*(EijI(2,1)+EijI(2,2)-EijI(4,1)+3*(EijI(4,2)-EijI(4
     -   ,3))+EijI(4,4)+EijI(5,3)-EijI(7,3)-4*(EijI(9,2)-EijI(14,3))+
     -   5*EijI(15,3)+EijI(17,4)-6*EijI(18,3)+EijI(19,4)+EijI(21,4)+2
     -   *(EijI(5,2)-EijI(7,2)+EijI(8,3)-EijI(10,3)-EijI(14,4)-EijI(1
     -   5,4)-EijI(24,4)-EijI(27,4)+2*EijI(32,4)))
       F(56)=DCMPLX(FR(56),FI(56))
       FR(57) = -16*(EijR(4,4)+EijR(5,2)+EijR(5,3)-EijR(7,2)-EijR(7,3)
     -   +EijR(8,2)-EijR(10,2)-3*(EijR(4,3)-EijR(15,3))-EijR(15,4)-Ei
     -   jR(16,4)+4*(EijR(14,3)-EijR(18,3))+EijR(19,4)+EijR(23,4)-Eij
     -   R(24,4)-EijR(25,4)+EijR(34,4)+2*(EijR(4,2)-EijR(9,2)-EijR(14
     -   ,4)+EijR(16,3)+EijR(17,3)-EijR(19,3)-EijR(20,3)+EijR(32,4)+E
     -   ijR(33,4)-EijR(35,4)))
       FI(57) = -16*(EijI(4,4)+EijI(5,2)+EijI(5,3)-EijI(7,2)-EijI(7,3)
     -   +EijI(8,2)-EijI(10,2)-3*(EijI(4,3)-EijI(15,3))-EijI(15,4)-Ei
     -   jI(16,4)+4*(EijI(14,3)-EijI(18,3))+EijI(19,4)+EijI(23,4)-Eij
     -   I(24,4)-EijI(25,4)+EijI(34,4)+2*(EijI(4,2)-EijI(9,2)-EijI(14
     -   ,4)+EijI(16,3)+EijI(17,3)-EijI(19,3)-EijI(20,3)+EijI(32,4)+E
     -   ijI(33,4)-EijI(35,4)))
       F(57)=DCMPLX(FR(57),FI(57))
       FR(58) = -16*(EijR(4,2)+EijR(4,4)-EijR(9,2)-EijR(15,4)+EijR(19,
     -   4)-EijR(24,4)-2*(EijR(4,3)-EijR(14,3)+EijR(14,4)-EijR(15,3)+
     -   EijR(18,3)-EijR(32,4)))
       FI(58) = -16*(EijI(4,2)+EijI(4,4)-EijI(9,2)-EijI(15,4)+EijI(19,
     -   4)-EijI(24,4)-2*(EijI(4,3)-EijI(14,3)+EijI(14,4)-EijI(15,3)+
     -   EijI(18,3)-EijI(32,4)))
       F(58)=DCMPLX(FR(58),FI(58))
       FR(59) = -16*(EijR(2,1)+EijR(2,3)-EijR(4,1)+EijR(4,4)+EijR(5,2)
     -   -EijR(7,2)+EijR(8,4)-5*(EijR(9,2)+EijR(10,3))-EijR(10,4)+2*(
     -   EijR(2,2)+EijR(8,3)+EijR(14,3))-EijR(14,4)+7*EijR(15,3)-4*Ei
     -   jR(18,3)+3*(EijR(4,2)-EijR(4,3)-EijR(15,4)+EijR(21,4)-EijR(2
     -   7,4)+EijR(32,4)))
       FI(59) = -16*(EijI(2,1)+EijI(2,3)-EijI(4,1)+EijI(4,4)+EijI(5,2)
     -   -EijI(7,2)+EijI(8,4)-5*(EijI(9,2)+EijI(10,3))-EijI(10,4)+2*(
     -   EijI(2,2)+EijI(8,3)+EijI(14,3))-EijI(14,4)+7*EijI(15,3)-4*Ei
     -   jI(18,3)+3*(EijI(4,2)-EijI(4,3)-EijI(15,4)+EijI(21,4)-EijI(2
     -   7,4)+EijI(32,4)))
       F(59)=DCMPLX(FR(59),FI(59))
       FR(60) = -16*(EijR(2,2)+EijR(2,4)+EijR(4,2)+EijR(4,4)+2*(EijR(2
     -   ,3)-EijR(4,3)-EijR(9,2))-4*(EijR(10,4)+EijR(15,4))-6*(EijR(1
     -   0,3)-EijR(15,3)-EijR(21,4)))
       FI(60) = -16*(EijI(2,2)+EijI(2,4)+EijI(4,2)+EijI(4,4)+2*(EijI(2
     -   ,3)-EijI(4,3)-EijI(9,2))-4*(EijI(10,4)+EijI(15,4))-6*(EijI(1
     -   0,3)-EijI(15,3)-EijI(21,4)))
       F(60)=DCMPLX(FR(60),FI(60))
       FR(61) = -16*(EijR(3,1)-EijR(4,1)-EijR(4,3)+EijR(4,4)+EijR(6,2)
     -   +EijR(7,2)-EijR(9,2)-EijR(10,2)-EijR(14,4)+EijR(15,3)-EijR(1
     -   5,4)+EijR(16,3)-EijR(16,4)-EijR(20,3)+EijR(32,4)+EijR(33,4)+
     -   EijR(34,4)-EijR(35,4))
       FI(61) = -16*(EijI(3,1)-EijI(4,1)-EijI(4,3)+EijI(4,4)+EijI(6,2)
     -   +EijI(7,2)-EijI(9,2)-EijI(10,2)-EijI(14,4)+EijI(15,3)-EijI(1
     -   5,4)+EijI(16,3)-EijI(16,4)-EijI(20,3)+EijI(32,4)+EijI(33,4)+
     -   EijI(34,4)-EijI(35,4))
       F(61)=DCMPLX(FR(61),FI(61))
       FR(62) = 16*(EijR(2,1)-EijR(4,1)-EijR(4,2)-EijR(4,4)+EijR(5,2)+
     -   EijR(7,2)-EijR(9,2)-EijR(14,3)+EijR(14,4)-EijR(15,3)+EijR(15
     -   ,4)+2*(EijR(4,3)-EijR(16,3))+EijR(16,4)+EijR(19,3)+EijR(20,3
     -   )-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4))
       FI(62) = 16*(EijI(2,1)-EijI(4,1)-EijI(4,2)-EijI(4,4)+EijI(5,2)+
     -   EijI(7,2)-EijI(9,2)-EijI(14,3)+EijI(14,4)-EijI(15,3)+EijI(15
     -   ,4)+2*(EijI(4,3)-EijI(16,3))+EijI(16,4)+EijI(19,3)+EijI(20,3
     -   )-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4))
       F(62)=DCMPLX(FR(62),FI(62))
       FR(63) = 16*(EijR(2,1)-EijR(4,1)-EijR(4,4)+EijR(5,3)-EijR(7,2)-
     -   EijR(9,2)-EijR(15,3)+EijR(15,4)-EijR(19,4)+EijR(24,4)+2*(Eij
     -   R(4,3)+EijR(5,2)-EijR(14,3)+EijR(14,4)-EijR(32,4)))
       FI(63) = 16*(EijI(2,1)-EijI(4,1)-EijI(4,4)+EijI(5,3)-EijI(7,2)-
     -   EijI(9,2)-EijI(15,3)+EijI(15,4)-EijI(19,4)+EijI(24,4)+2*(Eij
     -   I(4,3)+EijI(5,2)-EijI(14,3)+EijI(14,4)-EijI(32,4)))
       F(63)=DCMPLX(FR(63),FI(63))
       FR(64) = 16*(EijR(2,1)+EijR(2,2)-EijR(4,1)-EijR(4,4)+EijR(5,2)+
     -   EijR(8,3)-EijR(14,3)+EijR(14,4)-EijR(21,4)+EijR(27,4)+2*(Eij
     -   R(4,3)-EijR(9,2)-EijR(15,3)+EijR(15,4)-EijR(32,4)))
       FI(64) = 16*(EijI(2,1)+EijI(2,2)-EijI(4,1)-EijI(4,4)+EijI(5,2)+
     -   EijI(8,3)-EijI(14,3)+EijI(14,4)-EijI(21,4)+EijI(27,4)+2*(Eij
     -   I(4,3)-EijI(9,2)-EijI(15,3)+EijI(15,4)-EijI(32,4)))
       F(64)=DCMPLX(FR(64),FI(64))
       FR(65) = 16*(EijR(2,1)+2*EijR(2,2)+EijR(2,3)-EijR(4,1)+EijR(4,2
     -   )+EijR(4,3)-EijR(4,4)-EijR(10,3)+EijR(10,4)-EijR(15,3)-3*(Ei
     -   jR(9,2)-EijR(15,4)+EijR(21,4)))
       FI(65) = 16*(EijI(2,1)+2*EijI(2,2)+EijI(2,3)-EijI(4,1)+EijI(4,2
     -   )+EijI(4,3)-EijI(4,4)-EijI(10,3)+EijI(10,4)-EijI(15,3)-3*(Ei
     -   jI(9,2)-EijI(15,4)+EijI(21,4)))
       F(65)=DCMPLX(FR(65),FI(65))
       FR(66) = 16*(EijR(3,1)-EijR(4,1)+EijR(4,2)-EijR(4,4)-EijR(9,2)+
     -   EijR(9,3)-EijR(10,3)+EijR(15,3)+EijR(16,4)-EijR(20,3)-EijR(2
     -   1,4)+EijR(28,4)+2*(EijR(8,2)-EijR(10,2)+EijR(15,4)-EijR(34,4
     -   )))
       FI(66) = 16*(EijI(3,1)-EijI(4,1)+EijI(4,2)-EijI(4,4)-EijI(9,2)+
     -   EijI(9,3)-EijI(10,3)+EijI(15,3)+EijI(16,4)-EijI(20,3)-EijI(2
     -   1,4)+EijI(28,4)+2*(EijI(8,2)-EijI(10,2)+EijI(15,4)-EijI(34,4
     -   )))
       F(66)=DCMPLX(FR(66),FI(66))
       P(397) = p2sq-p4sq-3*s12+s45+2*P(6)
       FR(67) = 8*(-Dij2345R(1,1)+Dij2345R(3,1)-s12*EijR(2,1)+p3sq*Eij
     -   R(3,1)-4*EijR(11,2)-2*(EijR(39,4)-EijR(44,4))+EijR(5,2)*P(18
     -   )+EijR(8,2)*P(25)-EijR(4,1)*P(42)+EijR(10,2)*P(43)-EijR(2,2)
     -   *P(68)+EijR(7,2)*P(72)-EijR(4,2)*P(174)-EijR(9,2)*P(397))
       FI(67) = 8*(-Dij2345I(1,1)+Dij2345I(3,1)-s12*EijI(2,1)+p3sq*Eij
     -   I(3,1)-4*EijI(11,2)-2*(EijI(39,4)-EijI(44,4))+EijI(5,2)*P(18
     -   )+EijI(8,2)*P(25)-EijI(4,1)*P(42)+EijI(10,2)*P(43)-EijI(2,2)
     -   *P(68)+EijI(7,2)*P(72)-EijI(4,2)*P(174)-EijI(9,2)*P(397))
       F(67)=DCMPLX(FR(67),FI(67))
       FR(68) = 8*(-2*(EijR(22,3)-EijR(24,3)+EijR(37,4)+EijR(39,4)-2*E
     -   ijR(44,4))-(EijR(5,2)-EijR(7,2))*P(45)-(-EijR(2,2)+EijR(9,2)
     -   )*P(53)+Is45*(Cij123R(2,1)+Cij123R(2,2)+p3sq*(Dij1234R(2,1)-
     -   Dij1234R(3,1)-Dij1234R(6,2))+2*(Dij1234R(7,2)+Dij1234R(12,3)
     -   )-Dij1234R(4,2)*P(45)+Dij1234R(2,2)*P(78))-(EijR(2,1)-EijR(4
     -   ,1))*P(235))
       FI(68) = 8*(-2*(EijI(22,3)-EijI(24,3)+EijI(37,4)+EijI(39,4)-2*E
     -   ijI(44,4))-(EijI(5,2)-EijI(7,2))*P(45)-(-EijI(2,2)+EijI(9,2)
     -   )*P(53)+Is45*(Cij123I(2,1)+Cij123I(2,2)+p3sq*(Dij1234I(2,1)-
     -   Dij1234I(3,1)-Dij1234I(6,2))+2*(Dij1234I(7,2)+Dij1234I(12,3)
     -   )-Dij1234I(4,2)*P(45)+Dij1234I(2,2)*P(78))-(EijI(2,1)-EijI(4
     -   ,1))*P(235))
       F(68)=DCMPLX(FR(68),FI(68))
       P(398) = s23-2*s45
       P(399) = p2sq-p3sq-s23
       P(400) = p3sq+p4sq+s12-s34-2*s45
       P(401) = p4sq+s12-s34-2*s45
       P(402) = p2sq+p3sq+p4sq+s12-s34-2*s45
       P(403) = -p2sq+p3sq+p4sq+s12-s34
       P(404) = p3sq-s34
       P(405) = p2sq-p4sq-s12+3*s45-2*P(404)
       FR(69) = 8*(-D01234R+Dij2345R(1,1)-Dij2345R(3,1)+p3sq*(EijR(3,1
     -   )+2*EijR(6,2))+2*(EijR(11,2)-EijR(21,3)-EijR(22,3)+2*EijR(24
     -   ,3)-EijR(39,4)-EijR(40,4)+EijR(42,4)+EijR(44,4))-EijR(8,2)*P
     -   (25)-EijR(10,2)*P(43)+EijR(2,2)*P(68)+EijR(4,2)*P(174)+Dij12
     -   34R(1,2)*P(274)+Is45*(C0123R-C0134R+Cij123R(1,1)+Cij123R(2,1
     -   )+Cij123R(3,2)+6*Dij1234R(7,2)+2*(p2sq*Dij1234R(2,2)+Dij1234
     -   R(11,3))-Dij1234R(6,2)*P(1)+Dij1234R(5,2)*P(10)-Dij1234R(2,1
     -   )*P(238)-Dij1234R(4,2)*P(239)+Dij1234R(1,1)*P(398)+Dij1234R(
     -   3,1)*P(399))-EijR(2,1)*P(400)+EijR(4,1)*P(401)-EijR(5,2)*P(4
     -   02)-EijR(7,2)*P(403)-EijR(9,2)*P(405))
       FI(69) = 8*(-D01234I+Dij2345I(1,1)-Dij2345I(3,1)+p3sq*(EijI(3,1
     -   )+2*EijI(6,2))+2*(EijI(11,2)-EijI(21,3)-EijI(22,3)+2*EijI(24
     -   ,3)-EijI(39,4)-EijI(40,4)+EijI(42,4)+EijI(44,4))-EijI(8,2)*P
     -   (25)-EijI(10,2)*P(43)+EijI(2,2)*P(68)+EijI(4,2)*P(174)+Dij12
     -   34I(1,2)*P(274)+Is45*(C0123I-C0134I+Cij123I(1,1)+Cij123I(2,1
     -   )+Cij123I(3,2)+6*Dij1234I(7,2)+2*(p2sq*Dij1234I(2,2)+Dij1234
     -   I(11,3))-Dij1234I(6,2)*P(1)+Dij1234I(5,2)*P(10)-Dij1234I(2,1
     -   )*P(238)-Dij1234I(4,2)*P(239)+Dij1234I(1,1)*P(398)+Dij1234I(
     -   3,1)*P(399))-EijI(2,1)*P(400)+EijI(4,1)*P(401)-EijI(5,2)*P(4
     -   02)-EijI(7,2)*P(403)-EijI(9,2)*P(405))
       F(69)=DCMPLX(FR(69),FI(69))
       P(406) = -p4sq-s12+s34+s45+2*P(1)
       P(407) = 2*p3sq+p4sq-s12-s34+s45
       FR(70) = 8*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)+EijR(5,2)*P(18)
     -   -2*(p3sq*EijR(3,2)+4*EijR(11,2)+EijR(39,4)+EijR(43,4)-EijR(4
     -   4,4)-EijR(45,4)-EijR(6,2)*P(46))-EijR(2,2)*P(68)+EijR(7,2)*P
     -   (72)+Is45*(C0134R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,2)-2*(D
     -   ij1234R(7,2)-Dij1234R(13,3))-Dij1234R(3,1)*P(10)-Dij1234R(5,
     -   2)*P(45)+Dij1234R(6,2)*P(78))-EijR(4,2)*P(174)+(-EijR(3,1)+E
     -   ijR(4,1))*P(235)-EijR(9,2)*P(397)+EijR(8,2)*P(406)+EijR(10,2
     -   )*P(407))
       FI(70) = 8*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)+EijI(5,2)*P(18)
     -   -2*(p3sq*EijI(3,2)+4*EijI(11,2)+EijI(39,4)+EijI(43,4)-EijI(4
     -   4,4)-EijI(45,4)-EijI(6,2)*P(46))-EijI(2,2)*P(68)+EijI(7,2)*P
     -   (72)+Is45*(C0134I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,2)-2*(D
     -   ij1234I(7,2)-Dij1234I(13,3))-Dij1234I(3,1)*P(10)-Dij1234I(5,
     -   2)*P(45)+Dij1234I(6,2)*P(78))-EijI(4,2)*P(174)+(-EijI(3,1)+E
     -   ijI(4,1))*P(235)-EijI(9,2)*P(397)+EijI(8,2)*P(406)+EijI(10,2
     -   )*P(407))
       F(70)=DCMPLX(FR(70),FI(70))
       P(408) = p2sq-s23+s34+s45-2*P(137)
       P(409) = p2sq-p4sq+s34+s45
       P(410) = p3sq+3*s12+s23-2*P(409)
       P(411) = -2*p3sq+s12+s23+s34+s45
       FR(71) = 8*(Dij2345R(1,1)-Dij2345R(3,1)+p3sq*(EijR(3,1)+2*EijR(
     -   6,2))+2*(2*EijR(11,2)-EijR(22,3)+EijR(24,3)-EijR(39,4)-EijR(
     -   40,4)+EijR(42,4)+EijR(44,4))-EijR(2,1)*P(9)-EijR(8,2)*P(25)-
     -   EijR(10,2)*P(43)+EijR(4,1)*P(66)+EijR(2,2)*P(68)-EijR(5,2)*P
     -   (124)+Is12*(Cij145R(1,1)+Cij145R(1,2)-Cij145R(2,1)+Cij145R(2
     -   ,2)+p3sq*(Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R(4,2)-Dij1345R
     -   (6,2))-2*(Cij145R(3,2)-Dij1345R(7,2)-Dij1345R(11,3)+Dij1345R
     -   (13,3))+(Dij1345R(3,2)-Dij1345R(5,2))*P(203))-EijR(4,2)*P(40
     -   8)-EijR(7,2)*P(410)-EijR(9,2)*P(411))
       FI(71) = 8*(Dij2345I(1,1)-Dij2345I(3,1)+p3sq*(EijI(3,1)+2*EijI(
     -   6,2))+2*(2*EijI(11,2)-EijI(22,3)+EijI(24,3)-EijI(39,4)-EijI(
     -   40,4)+EijI(42,4)+EijI(44,4))-EijI(2,1)*P(9)-EijI(8,2)*P(25)-
     -   EijI(10,2)*P(43)+EijI(4,1)*P(66)+EijI(2,2)*P(68)-EijI(5,2)*P
     -   (124)+Is12*(Cij145I(1,1)+Cij145I(1,2)-Cij145I(2,1)+Cij145I(2
     -   ,2)+p3sq*(Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I(4,2)-Dij1345I
     -   (6,2))-2*(Cij145I(3,2)-Dij1345I(7,2)-Dij1345I(11,3)+Dij1345I
     -   (13,3))+(Dij1345I(3,2)-Dij1345I(5,2))*P(203))-EijI(4,2)*P(40
     -   8)-EijI(7,2)*P(410)-EijI(9,2)*P(411))
       F(71)=DCMPLX(FR(71),FI(71))
       P(412) = -p2sq+p4sq+s23-s34+2*P(52)
       FR(72) = 8*(-2*(EijR(22,3)-EijR(24,3)+EijR(37,4)+EijR(39,4)-2*E
     -   ijR(44,4))+(EijR(4,1)+EijR(7,2))*P(45)+(EijR(2,1)+EijR(5,2))
     -   *P(46)+EijR(9,2)*P(53)+Is12*(Cij145R(1,1)+Cij145R(1,2)-Cij14
     -   5R(2,1)+Cij145R(2,2)+p3sq*(Dij1345R(2,1)-Dij1345R(3,1)+Dij13
     -   45R(4,2)-Dij1345R(6,2))-2*(Cij145R(3,2)-Dij1345R(7,2)-Dij134
     -   5R(11,3)+Dij1345R(13,3))+(Dij1345R(3,2)-Dij1345R(5,2))*P(203
     -   ))+EijR(4,2)*P(412))
       FI(72) = 8*(-2*(EijI(22,3)-EijI(24,3)+EijI(37,4)+EijI(39,4)-2*E
     -   ijI(44,4))+(EijI(4,1)+EijI(7,2))*P(45)+(EijI(2,1)+EijI(5,2))
     -   *P(46)+EijI(9,2)*P(53)+Is12*(Cij145I(1,1)+Cij145I(1,2)-Cij14
     -   5I(2,1)+Cij145I(2,2)+p3sq*(Dij1345I(2,1)-Dij1345I(3,1)+Dij13
     -   45I(4,2)-Dij1345I(6,2))-2*(Cij145I(3,2)-Dij1345I(7,2)-Dij134
     -   5I(11,3)+Dij1345I(13,3))+(Dij1345I(3,2)-Dij1345I(5,2))*P(203
     -   ))+EijI(4,2)*P(412))
       F(72)=DCMPLX(FR(72),FI(72))
       P(413) = p2sq+2*p3sq+s12-s23-s45
       P(414) = p2sq+2*p3sq-s23-3*P(52)
       FR(73) = 4*Is12s45*(-3*B013R-B014R+4*(Cij123R(4,2)+Cij134R(4,2)
     -   +Cij145R(4,2))-2*(p2sq*Cij123R(2,1)-(C0145R+Cij145R(2,1))*P(
     -   15)-(C0123R+Cij123R(1,1))*P(18)-Cij134R(2,1)*P(30)-C0134R*P(
     -   42)))-8*(-D02345R+Is45*(Cij134R(1,1)+2*(Dij1234R(7,2)+Dij123
     -   4R(12,3)-Dij1234R(13,3)))+Dij2345R(1,1)-Dij2345R(3,1)-s15*(E
     -   E0R+EijR(1,1))-EijR(5,2)*P(18)-(-EijR(3,1)+EijR(4,1))*P(45)+
     -   2*(p3sq*EijR(3,2)+5*EijR(11,2)+EijR(23,3)-EijR(24,3)+EijR(39
     -   ,4)+EijR(43,4)-EijR(44,4)-EijR(45,4)-EijR(6,2)*P(46))+EijR(2
     -   ,2)*P(68)-EijR(7,2)*P(72)-(D01345R+Dij1345R(1,1))*P(92)-Is12
     -   *(Cij145R(1,2)+Cij145R(2,2)-2*(Cij145R(3,2)+Dij1345R(7,2)-Di
     -   j1345R(12,3)+Dij1345R(13,3))+(Dij1345R(2,1)-Dij1345R(3,1))*P
     -   (30)+(Dij1345R(4,2)-Dij1345R(5,2))*P(95)-Dij1345R(2,2)*P(128
     -   )+Dij1345R(3,2)*P(260)-Dij1345R(6,2)*P(264))+EijR(4,2)*P(272
     -   )-(D01234R+Dij1234R(1,1))*P(274)-EijR(9,2)*P(411)-EijR(8,2)*
     -   P(413)-EijR(10,2)*P(414))
       FI(73) = 4*Is12s45*(-3*B013I-B014I+4*(Cij123I(4,2)+Cij134I(4,2)
     -   +Cij145I(4,2))-2*(p2sq*Cij123I(2,1)-(C0145I+Cij145I(2,1))*P(
     -   15)-(C0123I+Cij123I(1,1))*P(18)-Cij134I(2,1)*P(30)-C0134I*P(
     -   42)))-8*(-D02345I+Is45*(Cij134I(1,1)+2*(Dij1234I(7,2)+Dij123
     -   4I(12,3)-Dij1234I(13,3)))+Dij2345I(1,1)-Dij2345I(3,1)-s15*(E
     -   E0I+EijI(1,1))-EijI(5,2)*P(18)-(-EijI(3,1)+EijI(4,1))*P(45)+
     -   2*(p3sq*EijI(3,2)+5*EijI(11,2)+EijI(23,3)-EijI(24,3)+EijI(39
     -   ,4)+EijI(43,4)-EijI(44,4)-EijI(45,4)-EijI(6,2)*P(46))+EijI(2
     -   ,2)*P(68)-EijI(7,2)*P(72)-(D01345I+Dij1345I(1,1))*P(92)-Is12
     -   *(Cij145I(1,2)+Cij145I(2,2)-2*(Cij145I(3,2)+Dij1345I(7,2)-Di
     -   j1345I(12,3)+Dij1345I(13,3))+(Dij1345I(2,1)-Dij1345I(3,1))*P
     -   (30)+(Dij1345I(4,2)-Dij1345I(5,2))*P(95)-Dij1345I(2,2)*P(128
     -   )+Dij1345I(3,2)*P(260)-Dij1345I(6,2)*P(264))+EijI(4,2)*P(272
     -   )-(D01234I+Dij1234I(1,1))*P(274)-EijI(9,2)*P(411)-EijI(8,2)*
     -   P(413)-EijI(10,2)*P(414))
       F(73)=DCMPLX(FR(73),FI(73))
       FR(74) = 8*(D01345R-Dij2345R(1,1)+Dij2345R(3,1)-s12*EijR(2,1)+p
     -   3sq*EijR(3,1)-2*(EijR(11,2)-EijR(22,3)+EijR(39,4)-EijR(44,4)
     -   )+EijR(5,2)*P(18)+EijR(8,2)*P(25)-EijR(4,1)*P(42)+EijR(10,2)
     -   *P(43)-EijR(2,2)*P(68)+EijR(7,2)*P(72)-EijR(4,2)*P(272)+Is12
     -   *(C0134R-C0145R-Cij145R(1,1)+Cij145R(2,2)-Cij145R(3,2)-4*Dij
     -   1345R(7,2)-2*Dij1345R(13,3)-Dij1345R(5,2)*P(95)+(Dij1345R(2,
     -   1)+Dij1345R(6,2))*P(128)+Dij1345R(1,1)*P(136)+Dij1345R(3,2)*
     -   P(260)-Dij1345R(3,1)*P(277))+EijR(9,2)*P(411))
       FI(74) = 8*(D01345I-Dij2345I(1,1)+Dij2345I(3,1)-s12*EijI(2,1)+p
     -   3sq*EijI(3,1)-2*(EijI(11,2)-EijI(22,3)+EijI(39,4)-EijI(44,4)
     -   )+EijI(5,2)*P(18)+EijI(8,2)*P(25)-EijI(4,1)*P(42)+EijI(10,2)
     -   *P(43)-EijI(2,2)*P(68)+EijI(7,2)*P(72)-EijI(4,2)*P(272)+Is12
     -   *(C0134I-C0145I-Cij145I(1,1)+Cij145I(2,2)-Cij145I(3,2)-4*Dij
     -   1345I(7,2)-2*Dij1345I(13,3)-Dij1345I(5,2)*P(95)+(Dij1345I(2,
     -   1)+Dij1345I(6,2))*P(128)+Dij1345I(1,1)*P(136)+Dij1345I(3,2)*
     -   P(260)-Dij1345I(3,1)*P(277))+EijI(9,2)*P(411))
       F(74)=DCMPLX(FR(74),FI(74))
       P(415) = -1+Is45*s12
       P(416) = p2sq-s12+s15-s23+s45
       FR(75) = 4*Is12s45*(-3*B013R-B014R+4*(Cij123R(4,2)+Cij134R(4,2)
     -   +Cij145R(4,2))-2*(p2sq*Cij123R(2,1)-p4sq*Cij145R(2,1)-C0145R
     -   *P(15)-(C0123R+Cij123R(1,1))*P(18)-C0134R*P(42)-Cij134R(2,1)
     -   *P(43)))+8*(-Dij1345R(1,2)+s15*(EE0R+EijR(1,1))-2*(EijR(11,2
     -   )+EijR(37,4)+EijR(39,4)+2*(EijR(22,3)-EijR(24,3)-EijR(44,4))
     -   )+Is45*(Cij134R(1,2)-s12*Dij1234R(2,1)-p2sq*Dij1234R(2,2)+2*
     -   Dij1234R(12,3)+Dij1234R(4,2)*P(18))+Cij134R(3,2)*P(44)-(EijR
     -   (3,1)+EijR(8,2)-EijR(10,2))*P(53)+D01345R*P(92)-EijR(4,1)*P(
     -   99)-Is12*(Cij145R(1,1)-s45*Dij1345R(3,1)+p4sq*Dij1345R(3,2)-
     -   2*(Dij1345R(11,3)-Dij1345R(13,3))-Dij1345R(5,2)*P(137))+(D01
     -   234R+Dij1234R(1,1))*P(274)+(Dij1345R(2,1)+Dij1345R(4,2))*P(3
     -   40)+Dij1345R(1,1)*P(341)+Dij1345R(6,2)*P(342)+(Dij1234R(3,1)
     -   +Dij1234R(6,2))*P(415)+(EijR(2,1)+EijR(5,2)-EijR(7,2))*P(416
     -   ))
       FI(75) = 4*Is12s45*(-3*B013I-B014I+4*(Cij123I(4,2)+Cij134I(4,2)
     -   +Cij145I(4,2))-2*(p2sq*Cij123I(2,1)-p4sq*Cij145I(2,1)-C0145I
     -   *P(15)-(C0123I+Cij123I(1,1))*P(18)-C0134I*P(42)-Cij134I(2,1)
     -   *P(43)))+8*(-Dij1345I(1,2)+s15*(EE0I+EijI(1,1))-2*(EijI(11,2
     -   )+EijI(37,4)+EijI(39,4)+2*(EijI(22,3)-EijI(24,3)-EijI(44,4))
     -   )+Is45*(Cij134I(1,2)-s12*Dij1234I(2,1)-p2sq*Dij1234I(2,2)+2*
     -   Dij1234I(12,3)+Dij1234I(4,2)*P(18))+Cij134I(3,2)*P(44)-(EijI
     -   (3,1)+EijI(8,2)-EijI(10,2))*P(53)+D01345I*P(92)-EijI(4,1)*P(
     -   99)-Is12*(Cij145I(1,1)-s45*Dij1345I(3,1)+p4sq*Dij1345I(3,2)-
     -   2*(Dij1345I(11,3)-Dij1345I(13,3))-Dij1345I(5,2)*P(137))+(D01
     -   234I+Dij1234I(1,1))*P(274)+(Dij1345I(2,1)+Dij1345I(4,2))*P(3
     -   40)+Dij1345I(1,1)*P(341)+Dij1345I(6,2)*P(342)+(Dij1234I(3,1)
     -   +Dij1234I(6,2))*P(415)+(EijI(2,1)+EijI(5,2)-EijI(7,2))*P(416
     -   ))
       F(75)=DCMPLX(FR(75),FI(75))
       P(417) = p4sq+s12+s15-s34-2*s45
       P(418) = -p2sq+p4sq+s23-s34-2*P(43)
       P(419) = p2sq-p3sq-p4sq-s12+s15+s34
       FR(76) = 4*Is12s45*(-3*B013R-B014R+4*(Cij123R(4,2)+Cij134R(4,2)
     -   +Cij145R(4,2))+2*(C0134R*p3sq-p2sq*Cij123R(2,1)+p4sq*Cij145R
     -   (2,1)+C0145R*P(15)+(C0123R+Cij123R(1,1))*P(18)+Cij134R(2,1)*
     -   P(43)))+8*(-Dij1345R(1,2)-Dij2345R(1,1)+Dij2345R(3,1)+s15*(E
     -   E0R+2*EijR(1,1)+EijR(1,2))-2*(3*EijR(11,2)+EijR(21,3)+EijR(2
     -   2,3)-2*EijR(24,3)+EijR(39,4)+EijR(40,4)-EijR(42,4)-EijR(44,4
     -   ))+EijR(2,1)*P(9)+EijR(8,2)*P(25)+Cij134R(3,2)*P(44)-EijR(2,
     -   2)*P(68)-EijR(3,1)*P(83)+Is45*(Cij134R(1,2)+s23*(D01234R+Dij
     -   1234R(1,1))-2*(p2sq*Dij1234R(2,2)+Dij1234R(7,2)-Dij1234R(11,
     -   3))+Dij1234R(6,2)*P(1)+Dij1234R(2,1)*P(18)-Dij1234R(3,1)*P(4
     -   5)-Dij1234R(5,2)*P(78)+Dij1234R(4,2)*P(86))+EijR(10,2)*P(89)
     -   +D01345R*P(92)+EijR(5,2)*P(124)-Is12*(Cij145R(1,1)-s45*Dij13
     -   45R(3,1)+p4sq*Dij1345R(3,2)-2*(Dij1345R(11,3)-Dij1345R(13,3)
     -   )-Dij1345R(5,2)*P(137))-EijR(4,2)*P(174)+(Dij1345R(2,1)+Dij1
     -   345R(4,2))*P(340)+Dij1345R(1,1)*P(341)+Dij1345R(6,2)*P(342)+
     -   EijR(9,2)*P(411)-EijR(4,1)*P(417)+EijR(6,2)*P(418)-EijR(7,2)
     -   *P(419))
       FI(76) = 4*Is12s45*(-3*B013I-B014I+4*(Cij123I(4,2)+Cij134I(4,2)
     -   +Cij145I(4,2))+2*(C0134I*p3sq-p2sq*Cij123I(2,1)+p4sq*Cij145I
     -   (2,1)+C0145I*P(15)+(C0123I+Cij123I(1,1))*P(18)+Cij134I(2,1)*
     -   P(43)))+8*(-Dij1345I(1,2)-Dij2345I(1,1)+Dij2345I(3,1)+s15*(E
     -   E0I+2*EijI(1,1)+EijI(1,2))-2*(3*EijI(11,2)+EijI(21,3)+EijI(2
     -   2,3)-2*EijI(24,3)+EijI(39,4)+EijI(40,4)-EijI(42,4)-EijI(44,4
     -   ))+EijI(2,1)*P(9)+EijI(8,2)*P(25)+Cij134I(3,2)*P(44)-EijI(2,
     -   2)*P(68)-EijI(3,1)*P(83)+Is45*(Cij134I(1,2)+s23*(D01234I+Dij
     -   1234I(1,1))-2*(p2sq*Dij1234I(2,2)+Dij1234I(7,2)-Dij1234I(11,
     -   3))+Dij1234I(6,2)*P(1)+Dij1234I(2,1)*P(18)-Dij1234I(3,1)*P(4
     -   5)-Dij1234I(5,2)*P(78)+Dij1234I(4,2)*P(86))+EijI(10,2)*P(89)
     -   +D01345I*P(92)+EijI(5,2)*P(124)-Is12*(Cij145I(1,1)-s45*Dij13
     -   45I(3,1)+p4sq*Dij1345I(3,2)-2*(Dij1345I(11,3)-Dij1345I(13,3)
     -   )-Dij1345I(5,2)*P(137))-EijI(4,2)*P(174)+(Dij1345I(2,1)+Dij1
     -   345I(4,2))*P(340)+Dij1345I(1,1)*P(341)+Dij1345I(6,2)*P(342)+
     -   EijI(9,2)*P(411)-EijI(4,1)*P(417)+EijI(6,2)*P(418)-EijI(7,2)
     -   *P(419))
       F(76)=DCMPLX(FR(76),FI(76))
       P(420) = -1+Is45*p2sq
       P(421) = p2sq-p4sq-s23+s34-2*P(25)
       P(422) = p2sq-s12-s23+s45-2*P(287)
       FR(77) = 8*(-D02345R+Dij2345R(1,1)-Dij2345R(3,1)+8*EijR(11,2)+4
     -   *EijR(24,3)-2*(EijR(22,3)+EijR(23,3)+EijR(39,4)+EijR(43,4)-E
     -   ijR(44,4)-EijR(45,4))-EijR(5,2)*P(18)+Is45*(Cij134R(1,1)+Cij
     -   134R(3,2)-p2sq*(Dij1234R(2,1)+Dij1234R(6,2))+2*(Dij1234R(7,2
     -   )+Dij1234R(13,3))+Dij1234R(5,2)*P(18))+(Cij134R(2,1)+Cij134R
     -   (2,2))*P(44)+EijR(2,2)*P(68)-(-EijR(3,1)+EijR(4,1))*P(99)+Ei
     -   jR(6,2)*P(102)-EijR(7,2)*P(103)+EijR(4,2)*P(174)-Is12*(C0134
     -   R+s34*Dij1345R(1,1)-s34*Dij1345R(2,1)-p4sq*Dij1345R(3,2)-4*D
     -   ij1345R(7,2)-2*(Dij1345R(12,3)-Dij1345R(13,3))+(Dij1345R(4,2
     -   )-Dij1345R(5,2))*P(112)-Dij1345R(2,2)*P(203)+Dij1345R(6,2)*P
     -   (378))-EijR(9,2)*P(411)-EijR(8,2)*P(413)+Dij1234R(3,2)*P(415
     -   )+Dij1234R(3,1)*P(420)-EijR(3,2)*P(421)+EijR(10,2)*P(422))
       FI(77) = 8*(-D02345I+Dij2345I(1,1)-Dij2345I(3,1)+8*EijI(11,2)+4
     -   *EijI(24,3)-2*(EijI(22,3)+EijI(23,3)+EijI(39,4)+EijI(43,4)-E
     -   ijI(44,4)-EijI(45,4))-EijI(5,2)*P(18)+Is45*(Cij134I(1,1)+Cij
     -   134I(3,2)-p2sq*(Dij1234I(2,1)+Dij1234I(6,2))+2*(Dij1234I(7,2
     -   )+Dij1234I(13,3))+Dij1234I(5,2)*P(18))+(Cij134I(2,1)+Cij134I
     -   (2,2))*P(44)+EijI(2,2)*P(68)-(-EijI(3,1)+EijI(4,1))*P(99)+Ei
     -   jI(6,2)*P(102)-EijI(7,2)*P(103)+EijI(4,2)*P(174)-Is12*(C0134
     -   I+s34*Dij1345I(1,1)-s34*Dij1345I(2,1)-p4sq*Dij1345I(3,2)-4*D
     -   ij1345I(7,2)-2*(Dij1345I(12,3)-Dij1345I(13,3))+(Dij1345I(4,2
     -   )-Dij1345I(5,2))*P(112)-Dij1345I(2,2)*P(203)+Dij1345I(6,2)*P
     -   (378))-EijI(9,2)*P(411)-EijI(8,2)*P(413)+Dij1234I(3,2)*P(415
     -   )+Dij1234I(3,1)*P(420)-EijI(3,2)*P(421)+EijI(10,2)*P(422))
       F(77)=DCMPLX(FR(77),FI(77))
       P(423) = p3sq-s12-s15
       FR(78) = -8*(D01345R-Dij2345R(1,1)+Dij2345R(3,1)-s12*EijR(2,1)+
     -   p3sq*EijR(3,1)-2*(2*EijR(11,2)+EijR(24,3)-EijR(39,4)+EijR(44
     -   ,4))+EijR(5,2)*P(18)+EijR(8,2)*P(25)-EijR(2,2)*P(68)+EijR(7,
     -   2)*P(103)-EijR(10,2)*P(143)-EijR(4,2)*P(174)+Is12*(C0134R-p4
     -   sq*(Dij1345R(3,1)+Dij1345R(3,2))-2*(Dij1345R(7,2)-Dij1345R(1
     -   3,3))-Dij1345R(5,2)*P(112)+Dij1345R(2,1)*P(128)+Dij1345R(1,1
     -   )*P(136)+Dij1345R(6,2)*P(203))+EijR(9,2)*P(411)-EijR(4,1)*P(
     -   423))
       FI(78) = -8*(D01345I-Dij2345I(1,1)+Dij2345I(3,1)-s12*EijI(2,1)+
     -   p3sq*EijI(3,1)-2*(2*EijI(11,2)+EijI(24,3)-EijI(39,4)+EijI(44
     -   ,4))+EijI(5,2)*P(18)+EijI(8,2)*P(25)-EijI(2,2)*P(68)+EijI(7,
     -   2)*P(103)-EijI(10,2)*P(143)-EijI(4,2)*P(174)+Is12*(C0134I-p4
     -   sq*(Dij1345I(3,1)+Dij1345I(3,2))-2*(Dij1345I(7,2)-Dij1345I(1
     -   3,3))-Dij1345I(5,2)*P(112)+Dij1345I(2,1)*P(128)+Dij1345I(1,1
     -   )*P(136)+Dij1345I(6,2)*P(203))+EijI(9,2)*P(411)-EijI(4,1)*P(
     -   423))
       F(78)=DCMPLX(FR(78),FI(78))
       FR(79) = -16*Is45*(2*(Dij1234R(2,1)-Dij1234R(3,1))+Dij1234R(4,3
     -   )+3*(Dij1234R(4,2)-Dij1234R(5,2))-Dij1234R(5,3)+s45*(EijR(2,
     -   1)-EijR(4,1)+EijR(4,4)+EijR(5,2)+EijR(6,2)+EijR(6,3)-EijR(7,
     -   3)+EijR(8,2)-EijR(15,4)+3*(EijR(4,2)-EijR(4,3)+EijR(16,3))-E
     -   ijR(16,4)+4*(EijR(14,3)-EijR(19,3))+EijR(19,4)+EijR(23,4)-Ei
     -   jR(24,4)-EijR(25,4)+EijR(34,4)-2*(EijR(7,2)+EijR(9,2)+EijR(1
     -   0,2)+EijR(14,4)-EijR(15,3)-EijR(17,3)+EijR(18,3)+EijR(20,3)-
     -   EijR(32,4)-EijR(33,4)+EijR(35,4))))
       FI(79) = -16*Is45*(2*(Dij1234I(2,1)-Dij1234I(3,1))+Dij1234I(4,3
     -   )+3*(Dij1234I(4,2)-Dij1234I(5,2))-Dij1234I(5,3)+s45*(EijI(2,
     -   1)-EijI(4,1)+EijI(4,4)+EijI(5,2)+EijI(6,2)+EijI(6,3)-EijI(7,
     -   3)+EijI(8,2)-EijI(15,4)+3*(EijI(4,2)-EijI(4,3)+EijI(16,3))-E
     -   ijI(16,4)+4*(EijI(14,3)-EijI(19,3))+EijI(19,4)+EijI(23,4)-Ei
     -   jI(24,4)-EijI(25,4)+EijI(34,4)-2*(EijI(7,2)+EijI(9,2)+EijI(1
     -   0,2)+EijI(14,4)-EijI(15,3)-EijI(17,3)+EijI(18,3)+EijI(20,3)-
     -   EijI(32,4)-EijI(33,4)+EijI(35,4))))
       F(79)=DCMPLX(FR(79),FI(79))
       FR(80) = -16*Is45*(Dij1234R(2,2)-Dij1234R(5,2)+2*(Dij1234R(2,1)
     -   -Dij1234R(3,1)+Dij1234R(4,2)-Dij1234R(6,2))+Dij1234R(6,3)-Di
     -   j1234R(10,3)+s45*(EijR(2,1)-EijR(4,1)+EijR(4,4)+EijR(5,2)+Ei
     -   jR(6,2)+EijR(8,2)+EijR(9,3)-EijR(10,3)-EijR(14,4)+3*(EijR(4,
     -   2)-EijR(4,3)+EijR(16,3))-EijR(16,4)+4*(EijR(15,3)-EijR(20,3)
     -   )+EijR(21,4)+EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(33,4)-2*(
     -   EijR(7,2)+EijR(9,2)+EijR(10,2)-EijR(14,3)+EijR(15,4)-EijR(17
     -   ,3)+EijR(18,3)+EijR(19,3)-EijR(32,4)-EijR(34,4)+EijR(35,4)))
     -   )
       FI(80) = -16*Is45*(Dij1234I(2,2)-Dij1234I(5,2)+2*(Dij1234I(2,1)
     -   -Dij1234I(3,1)+Dij1234I(4,2)-Dij1234I(6,2))+Dij1234I(6,3)-Di
     -   j1234I(10,3)+s45*(EijI(2,1)-EijI(4,1)+EijI(4,4)+EijI(5,2)+Ei
     -   jI(6,2)+EijI(8,2)+EijI(9,3)-EijI(10,3)-EijI(14,4)+3*(EijI(4,
     -   2)-EijI(4,3)+EijI(16,3))-EijI(16,4)+4*(EijI(15,3)-EijI(20,3)
     -   )+EijI(21,4)+EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(33,4)-2*(
     -   EijI(7,2)+EijI(9,2)+EijI(10,2)-EijI(14,3)+EijI(15,4)-EijI(17
     -   ,3)+EijI(18,3)+EijI(19,3)-EijI(32,4)-EijI(34,4)+EijI(35,4)))
     -   )
       F(80)=DCMPLX(FR(80),FI(80))
       FR(81) = -16*Is45*(Dij1234R(2,1)+2*Dij1234R(2,2)-Dij1234R(3,1)-
     -   Dij1234R(5,2)-Dij1234R(6,2)+Dij1234R(6,3)-Dij1234R(10,3)+s45
     -   *(EijR(2,2)+EijR(4,2)+EijR(4,4)-EijR(5,2)+EijR(7,2)+EijR(8,2
     -   )+EijR(9,3)-EijR(10,2)-EijR(10,3)+EijR(14,3)-EijR(14,4)-EijR
     -   (16,4)+EijR(17,3)-EijR(18,3)-EijR(19,3)+3*(EijR(15,3)-EijR(2
     -   0,3))+EijR(21,4)+EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(33,4)
     -   -2*(EijR(4,3)+EijR(9,2)+EijR(15,4)-EijR(16,3)-EijR(32,4)-Eij
     -   R(34,4)+EijR(35,4))))
       FI(81) = -16*Is45*(Dij1234I(2,1)+2*Dij1234I(2,2)-Dij1234I(3,1)-
     -   Dij1234I(5,2)-Dij1234I(6,2)+Dij1234I(6,3)-Dij1234I(10,3)+s45
     -   *(EijI(2,2)+EijI(4,2)+EijI(4,4)-EijI(5,2)+EijI(7,2)+EijI(8,2
     -   )+EijI(9,3)-EijI(10,2)-EijI(10,3)+EijI(14,3)-EijI(14,4)-EijI
     -   (16,4)+EijI(17,3)-EijI(18,3)-EijI(19,3)+3*(EijI(15,3)-EijI(2
     -   0,3))+EijI(21,4)+EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(33,4)
     -   -2*(EijI(4,3)+EijI(9,2)+EijI(15,4)-EijI(16,3)-EijI(32,4)-Eij
     -   I(34,4)+EijI(35,4))))
       F(81)=DCMPLX(FR(81),FI(81))
       FR(82) = -16*Is45*(Dij1234R(2,1)+Dij1234R(2,3)-Dij1234R(3,1)+2*
     -   (Dij1234R(2,2)-Dij1234R(6,2))-Dij1234R(8,3)+s45*(EijR(4,2)+E
     -   ijR(4,4)+EijR(8,2)-EijR(9,2)+EijR(9,4)-EijR(10,2)-EijR(10,4)
     -   -2*(EijR(4,3)-EijR(9,3)+EijR(10,3)-EijR(16,3))-EijR(16,4)+4*
     -   (EijR(15,3)-EijR(20,3))-3*(EijR(15,4)-EijR(21,4)+EijR(28,4)-
     -   EijR(34,4))))
       FI(82) = -16*Is45*(Dij1234I(2,1)+Dij1234I(2,3)-Dij1234I(3,1)+2*
     -   (Dij1234I(2,2)-Dij1234I(6,2))-Dij1234I(8,3)+s45*(EijI(4,2)+E
     -   ijI(4,4)+EijI(8,2)-EijI(9,2)+EijI(9,4)-EijI(10,2)-EijI(10,4)
     -   -2*(EijI(4,3)-EijI(9,3)+EijI(10,3)-EijI(16,3))-EijI(16,4)+4*
     -   (EijI(15,3)-EijI(20,3))-3*(EijI(15,4)-EijI(21,4)+EijI(28,4)-
     -   EijI(34,4))))
       F(82)=DCMPLX(FR(82),FI(82))
       FR(83) = -16*Is45*(Dij1234R(2,1)-Dij1234R(3,1)-2*Dij1234R(3,2)+
     -   Dij1234R(4,2)+Dij1234R(6,2)-Dij1234R(7,3)+Dij1234R(10,3)+s45
     -   *(EijR(4,4)+EijR(6,2)-EijR(7,2)+EijR(8,2)-EijR(9,2)+EijR(11,
     -   3)+EijR(12,3)-EijR(14,4)-EijR(15,4)+5*EijR(16,3)+EijR(17,3)-
     -   EijR(18,3)-3*(EijR(4,3)+EijR(19,3)+EijR(20,3))+EijR(22,4)+Ei
     -   jR(29,4)-EijR(30,4)-EijR(31,4)+EijR(32,4)+2*(EijR(4,2)-EijR(
     -   10,2)-EijR(13,3)+EijR(14,3)+EijR(15,3)-EijR(16,4)+EijR(33,4)
     -   +EijR(34,4)-EijR(35,4))))
       FI(83) = -16*Is45*(Dij1234I(2,1)-Dij1234I(3,1)-2*Dij1234I(3,2)+
     -   Dij1234I(4,2)+Dij1234I(6,2)-Dij1234I(7,3)+Dij1234I(10,3)+s45
     -   *(EijI(4,4)+EijI(6,2)-EijI(7,2)+EijI(8,2)-EijI(9,2)+EijI(11,
     -   3)+EijI(12,3)-EijI(14,4)-EijI(15,4)+5*EijI(16,3)+EijI(17,3)-
     -   EijI(18,3)-3*(EijI(4,3)+EijI(19,3)+EijI(20,3))+EijI(22,4)+Ei
     -   jI(29,4)-EijI(30,4)-EijI(31,4)+EijI(32,4)+2*(EijI(4,2)-EijI(
     -   10,2)-EijI(13,3)+EijI(14,3)+EijI(15,3)-EijI(16,4)+EijI(33,4)
     -   +EijI(34,4)-EijI(35,4))))
       F(83)=DCMPLX(FR(83),FI(83))
       FR(84) = 16*Is45*(Dij1234R(3,2)+Dij1234R(5,2)-2*Dij1234R(6,2)+D
     -   ij1234R(7,3)-Dij1234R(10,3)+s45*(EijR(4,3)-EijR(4,4)+EijR(6,
     -   2)-EijR(7,2)-EijR(8,2)+EijR(9,2)-EijR(12,3)+EijR(13,3)+EijR(
     -   14,4)-EijR(15,3)+EijR(15,4)-EijR(22,4)-EijR(29,4)+EijR(30,4)
     -   +EijR(31,4)-EijR(32,4)-2*(EijR(16,3)-EijR(16,4)-EijR(20,3)+E
     -   ijR(33,4)+EijR(34,4)-EijR(35,4))))
       FI(84) = 16*Is45*(Dij1234I(3,2)+Dij1234I(5,2)-2*Dij1234I(6,2)+D
     -   ij1234I(7,3)-Dij1234I(10,3)+s45*(EijI(4,3)-EijI(4,4)+EijI(6,
     -   2)-EijI(7,2)-EijI(8,2)+EijI(9,2)-EijI(12,3)+EijI(13,3)+EijI(
     -   14,4)-EijI(15,3)+EijI(15,4)-EijI(22,4)-EijI(29,4)+EijI(30,4)
     -   +EijI(31,4)-EijI(32,4)-2*(EijI(16,3)-EijI(16,4)-EijI(20,3)+E
     -   ijI(33,4)+EijI(34,4)-EijI(35,4))))
       F(84)=DCMPLX(FR(84),FI(84))
       FR(85) = EijR(2,2)+EijR(2,3)+EijR(4,2)+EijR(4,4)+EijR(8,3)+EijR
     -   (8,4)-4*EijR(10,3)-EijR(10,4)+EijR(14,3)-EijR(14,4)+5*EijR(1
     -   5,3)-2*(EijR(4,3)+EijR(9,2)+EijR(18,3))-3*(EijR(15,4)-EijR(2
     -   1,4)+EijR(27,4)-EijR(32,4))
       FI(85) = EijI(2,2)+EijI(2,3)+EijI(4,2)+EijI(4,4)+EijI(8,3)+EijI
     -   (8,4)-4*EijI(10,3)-EijI(10,4)+EijI(14,3)-EijI(14,4)+5*EijI(1
     -   5,3)-2*(EijI(4,3)+EijI(9,2)+EijI(18,3))-3*(EijI(15,4)-EijI(2
     -   1,4)+EijI(27,4)-EijI(32,4))
       F(85)=DCMPLX(FR(85),FI(85))
       FR(86) = EijR(4,3)-EijR(4,4)+EijR(6,2)-EijR(7,2)-EijR(8,2)+EijR
     -   (9,2)-EijR(9,3)+EijR(10,3)+EijR(14,4)-EijR(16,3)+EijR(16,4)-
     -   EijR(21,4)-EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(33,4)-2*(Ei
     -   jR(15,3)-EijR(15,4)-EijR(20,3)+EijR(32,4)+EijR(34,4)-EijR(35
     -   ,4))
       FI(86) = EijI(4,3)-EijI(4,4)+EijI(6,2)-EijI(7,2)-EijI(8,2)+EijI
     -   (9,2)-EijI(9,3)+EijI(10,3)+EijI(14,4)-EijI(16,3)+EijI(16,4)-
     -   EijI(21,4)-EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(33,4)-2*(Ei
     -   jI(15,3)-EijI(15,4)-EijI(20,3)+EijI(32,4)+EijI(34,4)-EijI(35
     -   ,4))
       F(86)=DCMPLX(FR(86),FI(86))
       FR(87) = EijR(4,3)-EijR(4,4)-EijR(7,2)+EijR(9,2)+EijR(10,3)+Eij
     -   R(14,4)-EijR(21,4)+EijR(27,4)-2*(EijR(15,3)-EijR(15,4)+EijR(
     -   32,4))
       FI(87) = EijI(4,3)-EijI(4,4)-EijI(7,2)+EijI(9,2)+EijI(10,3)+Eij
     -   I(14,4)-EijI(21,4)+EijI(27,4)-2*(EijI(15,3)-EijI(15,4)+EijI(
     -   32,4))
       F(87)=DCMPLX(FR(87),FI(87))
       FR(88) = EijR(4,2)+EijR(4,4)+EijR(5,2)-EijR(7,2)+EijR(8,3)-EijR
     -   (9,2)+EijR(9,3)+EijR(14,3)-EijR(14,4)+4*EijR(15,3)+EijR(16,3
     -   )-EijR(16,4)+EijR(21,4)+EijR(26,4)-EijR(27,4)-EijR(28,4)+Eij
     -   R(33,4)-2*(EijR(4,3)+EijR(10,3)+EijR(15,4)+EijR(18,3)+EijR(2
     -   0,3)-EijR(32,4)-EijR(34,4)+EijR(35,4))
       FI(88) = EijI(4,2)+EijI(4,4)+EijI(5,2)-EijI(7,2)+EijI(8,3)-EijI
     -   (9,2)+EijI(9,3)+EijI(14,3)-EijI(14,4)+4*EijI(15,3)+EijI(16,3
     -   )-EijI(16,4)+EijI(21,4)+EijI(26,4)-EijI(27,4)-EijI(28,4)+Eij
     -   I(33,4)-2*(EijI(4,3)+EijI(10,3)+EijI(15,4)+EijI(18,3)+EijI(2
     -   0,3)-EijI(32,4)-EijI(34,4)+EijI(35,4))
       F(88)=DCMPLX(FR(88),FI(88))
       FR(89) = -EijR(4,3)+EijR(4,4)-EijR(11,3)+EijR(12,3)-EijR(14,4)-
     -   EijR(15,4)+EijR(16,3)+EijR(17,3)-EijR(18,3)+EijR(19,3)-3*Eij
     -   R(20,3)+EijR(22,4)+EijR(29,4)-EijR(30,4)-EijR(31,4)+EijR(32,
     -   4)+2*(EijR(15,3)-EijR(16,4)+EijR(33,4)+EijR(34,4)-EijR(35,4)
     -   )
       FI(89) = -EijI(4,3)+EijI(4,4)-EijI(11,3)+EijI(12,3)-EijI(14,4)-
     -   EijI(15,4)+EijI(16,3)+EijI(17,3)-EijI(18,3)+EijI(19,3)-3*Eij
     -   I(20,3)+EijI(22,4)+EijI(29,4)-EijI(30,4)-EijI(31,4)+EijI(32,
     -   4)+2*(EijI(15,3)-EijI(16,4)+EijI(33,4)+EijI(34,4)-EijI(35,4)
     -   )
       F(89)=DCMPLX(FR(89),FI(89))
       FR(90) = EijR(4,2)-EijR(4,3)+EijR(4,4)-EijR(7,2)-EijR(14,4)+2*E
     -   ijR(15,3)-EijR(15,4)-EijR(16,4)-EijR(18,3)+EijR(19,3)-EijR(2
     -   0,3)+EijR(32,4)+EijR(33,4)+EijR(34,4)-EijR(35,4)
       FI(90) = EijI(4,2)-EijI(4,3)+EijI(4,4)-EijI(7,2)-EijI(14,4)+2*E
     -   ijI(15,3)-EijI(15,4)-EijI(16,4)-EijI(18,3)+EijI(19,3)-EijI(2
     -   0,3)+EijI(32,4)+EijI(33,4)+EijI(34,4)-EijI(35,4)
       F(90)=DCMPLX(FR(90),FI(90))
       FR(91) = -EijR(4,3)+EijR(4,4)-EijR(10,3)-EijR(14,4)+EijR(21,4)-
     -   EijR(27,4)+2*(EijR(15,3)-EijR(15,4)+EijR(32,4))
       FI(91) = -EijI(4,3)+EijI(4,4)-EijI(10,3)-EijI(14,4)+EijI(21,4)-
     -   EijI(27,4)+2*(EijI(15,3)-EijI(15,4)+EijI(32,4))
       F(91)=DCMPLX(FR(91),FI(91))
       FR(92) = EijR(4,2)-EijR(4,4)-EijR(7,2)+EijR(14,3)+EijR(14,4)-Ei
     -   jR(15,3)+EijR(15,4)+EijR(16,4)-EijR(19,3)+EijR(20,3)-EijR(32
     -   ,4)-EijR(33,4)-EijR(34,4)+EijR(35,4)
       FI(92) = EijI(4,2)-EijI(4,4)-EijI(7,2)+EijI(14,3)+EijI(14,4)-Ei
     -   jI(15,3)+EijI(15,4)+EijI(16,4)-EijI(19,3)+EijI(20,3)-EijI(32
     -   ,4)-EijI(33,4)-EijI(34,4)+EijI(35,4)
       F(92)=DCMPLX(FR(92),FI(92))
       FR(93) = -EijR(4,4)+EijR(14,3)+EijR(14,4)-EijR(15,3)+EijR(15,4)
     -   -EijR(32,4)
       FI(93) = -EijI(4,4)+EijI(14,3)+EijI(14,4)-EijI(15,3)+EijI(15,4)
     -   -EijI(32,4)
       F(93)=DCMPLX(FR(93),FI(93))
       FR(94) = -EijR(4,3)+EijR(4,4)+EijR(9,3)-EijR(10,3)-EijR(14,4)+E
     -   ijR(16,3)-EijR(16,4)+EijR(21,4)+EijR(26,4)-EijR(27,4)-EijR(2
     -   8,4)+EijR(33,4)+2*(EijR(15,3)-EijR(15,4)-EijR(20,3)+EijR(32,
     -   4)+EijR(34,4)-EijR(35,4))
       FI(94) = -EijI(4,3)+EijI(4,4)+EijI(9,3)-EijI(10,3)-EijI(14,4)+E
     -   ijI(16,3)-EijI(16,4)+EijI(21,4)+EijI(26,4)-EijI(27,4)-EijI(2
     -   8,4)+EijI(33,4)+2*(EijI(15,3)-EijI(15,4)-EijI(20,3)+EijI(32,
     -   4)+EijI(34,4)-EijI(35,4))
       F(94)=DCMPLX(FR(94),FI(94))
       FR(95) = -EijR(4,3)+EijR(4,4)+EijR(9,3)+EijR(9,4)-EijR(10,3)-Ei
     -   jR(10,4)+EijR(16,3)-EijR(16,4)+2*(EijR(15,3)-EijR(20,3))-3*(
     -   EijR(15,4)-EijR(21,4)+EijR(28,4)-EijR(34,4))
       FI(95) = -EijI(4,3)+EijI(4,4)+EijI(9,3)+EijI(9,4)-EijI(10,3)-Ei
     -   jI(10,4)+EijI(16,3)-EijI(16,4)+2*(EijI(15,3)-EijI(20,3))-3*(
     -   EijI(15,4)-EijI(21,4)+EijI(28,4)-EijI(34,4))
       F(95)=DCMPLX(FR(95),FI(95))
       FR(96) = -EijR(4,3)+EijR(4,4)-EijR(10,3)-EijR(10,4)+2*EijR(15,3
     -   )-3*(EijR(15,4)-EijR(21,4))
       FI(96) = -EijI(4,3)+EijI(4,4)-EijI(10,3)-EijI(10,4)+2*EijI(15,3
     -   )-3*(EijI(15,4)-EijI(21,4))
       F(96)=DCMPLX(FR(96),FI(96))
       FR(97) = EijR(2,2)+EijR(4,4)+EijR(8,2)+EijR(8,3)+EijR(9,3)-EijR
     -   (10,2)-EijR(14,4)+5*EijR(15,3)-EijR(16,4)+EijR(17,3)-EijR(19
     -   ,3)-3*(EijR(4,3)+EijR(9,2)+EijR(18,3)+EijR(20,3))+EijR(21,4)
     -   +EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(33,4)+2*(EijR(4,2)-Ei
     -   jR(10,3)+EijR(14,3)-EijR(15,4)+EijR(16,3)+EijR(32,4)+EijR(34
     -   ,4)-EijR(35,4))
       FI(97) = EijI(2,2)+EijI(4,4)+EijI(8,2)+EijI(8,3)+EijI(9,3)-EijI
     -   (10,2)-EijI(14,4)+5*EijI(15,3)-EijI(16,4)+EijI(17,3)-EijI(19
     -   ,3)-3*(EijI(4,3)+EijI(9,2)+EijI(18,3)+EijI(20,3))+EijI(21,4)
     -   +EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(33,4)+2*(EijI(4,2)-Ei
     -   jI(10,3)+EijI(14,3)-EijI(15,4)+EijI(16,3)+EijI(32,4)+EijI(34
     -   ,4)-EijI(35,4))
       F(97)=DCMPLX(FR(97),FI(97))
       FR(98) = EijR(2,2)+EijR(2,3)+EijR(4,2)+EijR(4,4)+EijR(9,3)+EijR
     -   (9,4)-4*EijR(10,3)-EijR(10,4)+5*EijR(15,3)+EijR(16,3)-EijR(1
     -   6,4)-2*(EijR(4,3)+EijR(9,2)+EijR(20,3))-3*(EijR(15,4)-EijR(2
     -   1,4)+EijR(28,4)-EijR(34,4))
       FI(98) = EijI(2,2)+EijI(2,3)+EijI(4,2)+EijI(4,4)+EijI(9,3)+EijI
     -   (9,4)-4*EijI(10,3)-EijI(10,4)+5*EijI(15,3)+EijI(16,3)-EijI(1
     -   6,4)-2*(EijI(4,3)+EijI(9,2)+EijI(20,3))-3*(EijI(15,4)-EijI(2
     -   1,4)+EijI(28,4)-EijI(34,4))
       F(98)=DCMPLX(FR(98),FI(98))
       FR(99) = -EijR(4,3)+EijR(4,4)+EijR(9,3)-EijR(10,3)+EijR(16,3)+E
     -   ijR(20,4)+EijR(21,4)+EijR(22,4)+2*(EijR(15,3)-EijR(15,4)-Eij
     -   R(16,4)-EijR(20,3)-EijR(28,4)-EijR(31,4)+2*EijR(34,4))
       FI(99) = -EijI(4,3)+EijI(4,4)+EijI(9,3)-EijI(10,3)+EijI(16,3)+E
     -   ijI(20,4)+EijI(21,4)+EijI(22,4)+2*(EijI(15,3)-EijI(15,4)-Eij
     -   I(16,4)-EijI(20,3)-EijI(28,4)-EijI(31,4)+2*EijI(34,4))
       F(99)=DCMPLX(FR(99),FI(99))
       FR(100) = EijR(4,2)-EijR(4,3)+EijR(4,4)-EijR(9,2)-EijR(10,3)-Ei
     -   jR(16,4)+EijR(21,4)-EijR(28,4)+2*(EijR(15,3)-EijR(15,4)+EijR
     -   (34,4))
       FI(100) = EijI(4,2)-EijI(4,3)+EijI(4,4)-EijI(9,2)-EijI(10,3)-Ei
     -   jI(16,4)+EijI(21,4)-EijI(28,4)+2*(EijI(15,3)-EijI(15,4)+EijI
     -   (34,4))
       F(100)=DCMPLX(FR(100),FI(100))
       FR(101) = EijR(4,2)+EijR(4,4)-EijR(9,2)-EijR(10,3)+EijR(14,3)-E
     -   ijR(14,4)+3*EijR(15,3)-EijR(18,3)+EijR(21,4)-EijR(27,4)-2*(E
     -   ijR(4,3)+EijR(15,4)-EijR(32,4))
       FI(101) = EijI(4,2)+EijI(4,4)-EijI(9,2)-EijI(10,3)+EijI(14,3)-E
     -   ijI(14,4)+3*EijI(15,3)-EijI(18,3)+EijI(21,4)-EijI(27,4)-2*(E
     -   ijI(4,3)+EijI(15,4)-EijI(32,4))
       F(101)=DCMPLX(FR(101),FI(101))
       FR(102) = EijR(4,2)-EijR(4,4)-EijR(9,2)+EijR(16,4)-EijR(21,4)+E
     -   ijR(28,4)+2*(EijR(15,4)-EijR(34,4))
       FI(102) = EijI(4,2)-EijI(4,4)-EijI(9,2)+EijI(16,4)-EijI(21,4)+E
     -   ijI(28,4)+2*(EijI(15,4)-EijI(34,4))
       F(102)=DCMPLX(FR(102),FI(102))
       FR(103) = EijR(3,1)-EijR(4,1)-EijR(4,3)+EijR(4,4)+EijR(8,2)-Eij
     -   R(10,2)+EijR(15,3)+EijR(16,3)-EijR(16,4)-EijR(20,3)+EijR(21,
     -   4)-EijR(28,4)-2*(EijR(15,4)-EijR(34,4))
       FI(103) = EijI(3,1)-EijI(4,1)-EijI(4,3)+EijI(4,4)+EijI(8,2)-Eij
     -   I(10,2)+EijI(15,3)+EijI(16,3)-EijI(16,4)-EijI(20,3)+EijI(21,
     -   4)-EijI(28,4)-2*(EijI(15,4)-EijI(34,4))
       F(103)=DCMPLX(FR(103),FI(103))
       FR(104) = EijR(3,1)+EijR(3,2)-EijR(4,1)+EijR(4,2)-EijR(4,3)+Eij
     -   R(4,4)+EijR(15,3)-EijR(15,4)+EijR(16,3)-EijR(20,3)+EijR(22,4
     -   )-EijR(31,4)-2*(EijR(10,2)+EijR(16,4)-EijR(34,4))
       FI(104) = EijI(3,1)+EijI(3,2)-EijI(4,1)+EijI(4,2)-EijI(4,3)+Eij
     -   I(4,4)+EijI(15,3)-EijI(15,4)+EijI(16,3)-EijI(20,3)+EijI(22,4
     -   )-EijI(31,4)-2*(EijI(10,2)+EijI(16,4)-EijI(34,4))
       F(104)=DCMPLX(FR(104),FI(104))
       FR(105) = EijR(2,1)+EijR(2,2)-EijR(4,1)+EijR(4,3)-EijR(4,4)-Eij
     -   R(9,2)-EijR(15,3)-EijR(16,3)+EijR(16,4)+EijR(20,3)-EijR(21,4
     -   )+EijR(28,4)+2*(EijR(15,4)-EijR(34,4))
       FI(105) = EijI(2,1)+EijI(2,2)-EijI(4,1)+EijI(4,3)-EijI(4,4)-Eij
     -   I(9,2)-EijI(15,3)-EijI(16,3)+EijI(16,4)+EijI(20,3)-EijI(21,4
     -   )+EijI(28,4)+2*(EijI(15,4)-EijI(34,4))
       F(105)=DCMPLX(FR(105),FI(105))
       FR(106) = EijR(3,1)-EijR(4,1)+EijR(4,2)-EijR(4,4)+EijR(8,2)-Eij
     -   R(9,2)-EijR(10,2)+EijR(15,4)-EijR(22,4)+EijR(31,4)+2*(EijR(1
     -   6,4)-EijR(34,4))
       FI(106) = EijI(3,1)-EijI(4,1)+EijI(4,2)-EijI(4,4)+EijI(8,2)-Eij
     -   I(9,2)-EijI(10,2)+EijI(15,4)-EijI(22,4)+EijI(31,4)+2*(EijI(1
     -   6,4)-EijI(34,4))
       F(106)=DCMPLX(FR(106),FI(106))
       FR(107) = -EijR(4,4)-EijR(9,2)+EijR(10,2)+EijR(15,4)+EijR(16,4)
     -   -EijR(34,4)
       FI(107) = -EijI(4,4)-EijI(9,2)+EijI(10,2)+EijI(15,4)+EijI(16,4)
     -   -EijI(34,4)
       F(107)=DCMPLX(FR(107),FI(107))
       FR(108) = EijR(2,1)+EijR(2,2)-EijR(4,1)+EijR(4,2)+EijR(4,3)-Eij
     -   R(4,4)+EijR(5,2)-EijR(7,2)+EijR(8,3)+EijR(14,4)-EijR(15,3)-E
     -   ijR(18,3)-EijR(21,4)+EijR(27,4)-2*(EijR(9,2)-EijR(15,4)+EijR
     -   (32,4))
       FI(108) = EijI(2,1)+EijI(2,2)-EijI(4,1)+EijI(4,2)+EijI(4,3)-Eij
     -   I(4,4)+EijI(5,2)-EijI(7,2)+EijI(8,3)+EijI(14,4)-EijI(15,3)-E
     -   ijI(18,3)-EijI(21,4)+EijI(27,4)-2*(EijI(9,2)-EijI(15,4)+EijI
     -   (32,4))
       F(108)=DCMPLX(FR(108),FI(108))
       FR(109) = EijR(3,1)-EijR(4,1)+EijR(4,2)-EijR(4,4)+EijR(6,2)+Eij
     -   R(8,2)-EijR(9,2)-2*EijR(10,2)+EijR(14,3)+EijR(14,4)+EijR(15,
     -   4)+EijR(16,4)+EijR(17,3)-EijR(18,3)-EijR(19,3)-EijR(32,4)-Ei
     -   jR(33,4)-EijR(34,4)+EijR(35,4)
       FI(109) = EijI(3,1)-EijI(4,1)+EijI(4,2)-EijI(4,4)+EijI(6,2)+Eij
     -   I(8,2)-EijI(9,2)-2*EijI(10,2)+EijI(14,3)+EijI(14,4)+EijI(15,
     -   4)+EijI(16,4)+EijI(17,3)-EijI(18,3)-EijI(19,3)-EijI(32,4)-Ei
     -   jI(33,4)-EijI(34,4)+EijI(35,4)
       F(109)=DCMPLX(FR(109),FI(109))
       FR(110) = EijR(4,2)-EijR(4,4)-EijR(9,2)+EijR(14,3)+EijR(14,4)+E
     -   ijR(15,4)-EijR(18,3)-EijR(32,4)
       FI(110) = EijI(4,2)-EijI(4,4)-EijI(9,2)+EijI(14,3)+EijI(14,4)+E
     -   ijI(15,4)-EijI(18,3)-EijI(32,4)
       F(110)=DCMPLX(FR(110),FI(110))
       FR(111) = EijR(4,2)-EijR(4,4)-EijR(9,2)-EijR(10,3)+EijR(15,3)+2
     -   *EijR(15,4)-EijR(21,4)
       FI(111) = EijI(4,2)-EijI(4,4)-EijI(9,2)-EijI(10,3)+EijI(15,3)+2
     -   *EijI(15,4)-EijI(21,4)
       F(111)=DCMPLX(FR(111),FI(111))
       FR(112) = EijR(4,2)+2*EijR(4,3)-EijR(4,4)-EijR(7,2)+EijR(8,2)-E
     -   ijR(10,2)-EijR(14,3)+EijR(14,4)-EijR(15,3)+EijR(15,4)-EijR(1
     -   6,3)+EijR(16,4)+EijR(17,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+
     -   EijR(35,4)
       FI(112) = EijI(4,2)+2*EijI(4,3)-EijI(4,4)-EijI(7,2)+EijI(8,2)-E
     -   ijI(10,2)-EijI(14,3)+EijI(14,4)-EijI(15,3)+EijI(15,4)-EijI(1
     -   6,3)+EijI(16,4)+EijI(17,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+
     -   EijI(35,4)
       F(112)=DCMPLX(FR(112),FI(112))
       FR(113) = EijR(4,2)+EijR(4,3)-EijR(4,4)+EijR(8,2)-EijR(9,2)+Eij
     -   R(9,3)-EijR(10,2)-EijR(15,3)+EijR(16,4)-EijR(20,3)-EijR(21,4
     -   )+EijR(28,4)+2*(EijR(15,4)-EijR(34,4))
       FI(113) = EijI(4,2)+EijI(4,3)-EijI(4,4)+EijI(8,2)-EijI(9,2)+Eij
     -   I(9,3)-EijI(10,2)-EijI(15,3)+EijI(16,4)-EijI(20,3)-EijI(21,4
     -   )+EijI(28,4)+2*(EijI(15,4)-EijI(34,4))
       F(113)=DCMPLX(FR(113),FI(113))
       FR(114) = EijR(3,2)+EijR(4,2)-EijR(4,4)+EijR(12,3)-EijR(13,3)+E
     -   ijR(15,4)+EijR(16,3)-EijR(20,3)-EijR(22,4)+EijR(31,4)-2*(Eij
     -   R(10,2)-EijR(16,4)+EijR(34,4))
       FI(114) = EijI(3,2)+EijI(4,2)-EijI(4,4)+EijI(12,3)-EijI(13,3)+E
     -   ijI(15,4)+EijI(16,3)-EijI(20,3)-EijI(22,4)+EijI(31,4)-2*(Eij
     -   I(10,2)-EijI(16,4)+EijI(34,4))
       F(114)=DCMPLX(FR(114),FI(114))
       FR(115) = EijR(4,4)-EijR(15,4)-EijR(16,3)-EijR(16,4)+EijR(20,3)
     -   +EijR(34,4)
       FI(115) = EijI(4,4)-EijI(15,4)-EijI(16,3)-EijI(16,4)+EijI(20,3)
     -   +EijI(34,4)
       F(115)=DCMPLX(FR(115),FI(115))
       FR(116) = EijR(4,2)+EijR(4,3)-EijR(4,4)-EijR(9,2)+EijR(14,4)+Ei
     -   jR(15,4)-EijR(18,3)-EijR(32,4)
       FI(116) = EijI(4,2)+EijI(4,3)-EijI(4,4)-EijI(9,2)+EijI(14,4)+Ei
     -   jI(15,4)-EijI(18,3)-EijI(32,4)
       F(116)=DCMPLX(FR(116),FI(116))
       FR(117) = EijR(4,2)-EijR(4,3)-EijR(4,4)-EijR(10,2)+EijR(15,3)+E
     -   ijR(15,4)+EijR(16,3)+EijR(16,4)-EijR(20,3)-EijR(34,4)
       FI(117) = EijI(4,2)-EijI(4,3)-EijI(4,4)-EijI(10,2)+EijI(15,3)+E
     -   ijI(15,4)+EijI(16,3)+EijI(16,4)-EijI(20,3)-EijI(34,4)
       F(117)=DCMPLX(FR(117),FI(117))
       FR(118) = -EijR(4,3)-EijR(4,4)+EijR(15,3)+EijR(15,4)
       FI(118) = -EijI(4,3)-EijI(4,4)+EijI(15,3)+EijI(15,4)
       F(118)=DCMPLX(FR(118),FI(118))
       FR(119) = Is45*(Dij1234R(2,1)+Dij1234R(2,2)-Dij1234R(3,1)-Dij12
     -   34R(3,2)+Dij1234R(8,3)-Dij1234R(9,3)+s45*(EijR(4,2)+EijR(4,4
     -   )+EijR(8,2)-EijR(9,2)+EijR(9,3)-EijR(10,2)-EijR(10,3)+EijR(1
     -   2,3)-EijR(13,3)+3*(EijR(15,3)+EijR(16,3))-4*EijR(20,3)+EijR(
     -   20,4)+EijR(21,4)+EijR(22,4)-2*(EijR(4,3)+EijR(15,4)+EijR(16,
     -   4)+EijR(28,4)+EijR(31,4)-2*EijR(34,4))))
       FI(119) = Is45*(Dij1234I(2,1)+Dij1234I(2,2)-Dij1234I(3,1)-Dij12
     -   34I(3,2)+Dij1234I(8,3)-Dij1234I(9,3)+s45*(EijI(4,2)+EijI(4,4
     -   )+EijI(8,2)-EijI(9,2)+EijI(9,3)-EijI(10,2)-EijI(10,3)+EijI(1
     -   2,3)-EijI(13,3)+3*(EijI(15,3)+EijI(16,3))-4*EijI(20,3)+EijI(
     -   20,4)+EijI(21,4)+EijI(22,4)-2*(EijI(4,3)+EijI(15,4)+EijI(16,
     -   4)+EijI(28,4)+EijI(31,4)-2*EijI(34,4))))
       F(119)=DCMPLX(FR(119),FI(119))
       FR(120) = Is45*(-Dij1234R(3,2)+Dij1234R(6,2)+Dij1234R(8,3)-Dij1
     -   234R(9,3)+s45*(-EijR(4,3)+EijR(4,4)+EijR(12,3)-EijR(13,3)+Ei
     -   jR(15,3)+EijR(20,4)+EijR(21,4)+EijR(22,4)-2*(EijR(15,4)-EijR
     -   (16,3)+EijR(16,4)+EijR(20,3)+EijR(28,4)+EijR(31,4)-2*EijR(34
     -   ,4))))
       FI(120) = Is45*(-Dij1234I(3,2)+Dij1234I(6,2)+Dij1234I(8,3)-Dij1
     -   234I(9,3)+s45*(-EijI(4,3)+EijI(4,4)+EijI(12,3)-EijI(13,3)+Ei
     -   jI(15,3)+EijI(20,4)+EijI(21,4)+EijI(22,4)-2*(EijI(15,4)-EijI
     -   (16,3)+EijI(16,4)+EijI(20,3)+EijI(28,4)+EijI(31,4)-2*EijI(34
     -   ,4))))
       F(120)=DCMPLX(FR(120),FI(120))
       FR(121) = Is45*(-Dij1234R(3,2)-Dij1234R(3,3)+Dij1234R(6,2)+Dij1
     -   234R(9,3)+s45*(-EijR(4,3)+EijR(4,4)+EijR(12,3)+EijR(12,4)-Ei
     -   jR(13,3)-EijR(13,4)+EijR(15,3)-EijR(15,4)+2*(EijR(16,3)-EijR
     -   (20,3))-3*(EijR(16,4)-EijR(22,4)+EijR(31,4)-EijR(34,4))))
       FI(121) = Is45*(-Dij1234I(3,2)-Dij1234I(3,3)+Dij1234I(6,2)+Dij1
     -   234I(9,3)+s45*(-EijI(4,3)+EijI(4,4)+EijI(12,3)+EijI(12,4)-Ei
     -   jI(13,3)-EijI(13,4)+EijI(15,3)-EijI(15,4)+2*(EijI(16,3)-EijI
     -   (20,3))-3*(EijI(16,4)-EijI(22,4)+EijI(31,4)-EijI(34,4))))
       F(121)=DCMPLX(FR(121),FI(121))
       FR(122) = -16*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,3)+D
     -   ij1345R(4,3)-Dij1345R(5,3)+Dij1345R(9,3)+2*(Dij1345R(3,2)+Di
     -   j1345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij134
     -   5R(10,3))+s12*(EijR(3,1)-EijR(4,1)+EijR(4,4)+EijR(6,4)-EijR(
     -   7,4)-EijR(16,4)+6*(EijR(14,3)-EijR(19,3))+3*(EijR(4,2)-EijR(
     -   4,3)+EijR(6,2)+EijR(6,3)-EijR(7,2)-EijR(7,3)-EijR(10,2)-EijR
     -   (14,4)+EijR(16,3)+EijR(19,4)-EijR(25,4)+EijR(33,4))))
       FI(122) = -16*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,3)+D
     -   ij1345I(4,3)-Dij1345I(5,3)+Dij1345I(9,3)+2*(Dij1345I(3,2)+Di
     -   j1345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij134
     -   5I(10,3))+s12*(EijI(3,1)-EijI(4,1)+EijI(4,4)+EijI(6,4)-EijI(
     -   7,4)-EijI(16,4)+6*(EijI(14,3)-EijI(19,3))+3*(EijI(4,2)-EijI(
     -   4,3)+EijI(6,2)+EijI(6,3)-EijI(7,2)-EijI(7,3)-EijI(10,2)-EijI
     -   (14,4)+EijI(16,3)+EijI(19,4)-EijI(25,4)+EijI(33,4))))
       F(122)=DCMPLX(FR(122),FI(122))
       FR(123) = -16*(EijR(4,4)+EijR(6,2)+EijR(6,3)-EijR(7,2)-EijR(7,3
     -   )+EijR(8,2)-EijR(9,2)-EijR(15,4)-3*(EijR(4,3)-EijR(16,3))-Ei
     -   jR(16,4)+4*(EijR(14,3)-EijR(19,3))+EijR(19,4)+EijR(23,4)-Eij
     -   R(24,4)-EijR(25,4)+EijR(34,4)+2*(EijR(4,2)-EijR(10,2)-EijR(1
     -   4,4)+EijR(15,3)+EijR(17,3)-EijR(18,3)-EijR(20,3)+EijR(32,4)+
     -   EijR(33,4)-EijR(35,4)))
       FI(123) = -16*(EijI(4,4)+EijI(6,2)+EijI(6,3)-EijI(7,2)-EijI(7,3
     -   )+EijI(8,2)-EijI(9,2)-EijI(15,4)-3*(EijI(4,3)-EijI(16,3))-Ei
     -   jI(16,4)+4*(EijI(14,3)-EijI(19,3))+EijI(19,4)+EijI(23,4)-Eij
     -   I(24,4)-EijI(25,4)+EijI(34,4)+2*(EijI(4,2)-EijI(10,2)-EijI(1
     -   4,4)+EijI(15,3)+EijI(17,3)-EijI(18,3)-EijI(20,3)+EijI(32,4)+
     -   EijI(33,4)-EijI(35,4)))
       F(123)=DCMPLX(FR(123),FI(123))
       FR(124) = -16*(EijR(3,2)+EijR(4,2)+EijR(4,4)+EijR(18,4)+4*(EijR
     -   (16,3)-EijR(19,3))+EijR(19,4)+EijR(22,4)-2*(EijR(4,3)+EijR(1
     -   0,2)-EijR(11,3)+EijR(13,3)-EijR(14,3)+EijR(14,4)+EijR(16,4)+
     -   EijR(25,4)+EijR(30,4)-2*EijR(33,4)))
       FI(124) = -16*(EijI(3,2)+EijI(4,2)+EijI(4,4)+EijI(18,4)+4*(EijI
     -   (16,3)-EijI(19,3))+EijI(19,4)+EijI(22,4)-2*(EijI(4,3)+EijI(1
     -   0,2)-EijI(11,3)+EijI(13,3)-EijI(14,3)+EijI(14,4)+EijI(16,4)+
     -   EijI(25,4)+EijI(30,4)-2*EijI(33,4)))
       F(124)=DCMPLX(FR(124),FI(124))
       FR(125) = -16*(EijR(4,2)+EijR(4,4)-EijR(10,2)-EijR(16,4)+EijR(1
     -   9,4)-EijR(25,4)-2*(EijR(4,3)-EijR(14,3)+EijR(14,4)-EijR(16,3
     -   )+EijR(19,3)-EijR(33,4)))
       FI(125) = -16*(EijI(4,2)+EijI(4,4)-EijI(10,2)-EijI(16,4)+EijI(1
     -   9,4)-EijI(25,4)-2*(EijI(4,3)-EijI(14,3)+EijI(14,4)-EijI(16,3
     -   )+EijI(19,3)-EijI(33,4)))
       F(125)=DCMPLX(FR(125),FI(125))
       FR(126) = -16*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,3)+D
     -   ij1345R(4,3)-Dij1345R(5,3)+Dij1345R(9,3)+2*(Dij1345R(3,2)+Di
     -   j1345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij134
     -   5R(10,3))+s12*(EijR(3,1)-EijR(4,1)+EijR(4,4)-EijR(5,2)-EijR(
     -   5,3)+EijR(6,3)-EijR(7,2)+EijR(8,2)+EijR(15,3)-EijR(15,4)-3*(
     -   EijR(10,2)-EijR(16,3))-EijR(16,4)-4*EijR(19,3)+EijR(19,4)+Ei
     -   jR(23,4)-EijR(24,4)-EijR(25,4)+EijR(34,4)+2*(EijR(4,2)-EijR(
     -   4,3)+EijR(6,2)+EijR(14,3)-EijR(14,4)+EijR(17,3)-EijR(20,3)+E
     -   ijR(32,4)+EijR(33,4)-EijR(35,4))))
       FI(126) = -16*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,3)+D
     -   ij1345I(4,3)-Dij1345I(5,3)+Dij1345I(9,3)+2*(Dij1345I(3,2)+Di
     -   j1345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij134
     -   5I(10,3))+s12*(EijI(3,1)-EijI(4,1)+EijI(4,4)-EijI(5,2)-EijI(
     -   5,3)+EijI(6,3)-EijI(7,2)+EijI(8,2)+EijI(15,3)-EijI(15,4)-3*(
     -   EijI(10,2)-EijI(16,3))-EijI(16,4)-4*EijI(19,3)+EijI(19,4)+Ei
     -   jI(23,4)-EijI(24,4)-EijI(25,4)+EijI(34,4)+2*(EijI(4,2)-EijI(
     -   4,3)+EijI(6,2)+EijI(14,3)-EijI(14,4)+EijI(17,3)-EijI(20,3)+E
     -   ijI(32,4)+EijI(33,4)-EijI(35,4))))
       F(126)=DCMPLX(FR(126),FI(126))
       FR(127) = -16*(EijR(4,2)+EijR(4,4)-EijR(5,2)+EijR(6,2)+EijR(8,2
     -   )-EijR(8,3)+EijR(9,3)+EijR(14,3)-EijR(14,4)+3*EijR(16,3)-Eij
     -   R(16,4)-4*EijR(20,3)+EijR(21,4)+EijR(26,4)-EijR(27,4)-EijR(2
     -   8,4)+EijR(33,4)-2*(EijR(4,3)+EijR(10,2)-EijR(15,3)+EijR(15,4
     -   )-EijR(17,3)+EijR(19,3)-EijR(32,4)-EijR(34,4)+EijR(35,4)))
       FI(127) = -16*(EijI(4,2)+EijI(4,4)-EijI(5,2)+EijI(6,2)+EijI(8,2
     -   )-EijI(8,3)+EijI(9,3)+EijI(14,3)-EijI(14,4)+3*EijI(16,3)-Eij
     -   I(16,4)-4*EijI(20,3)+EijI(21,4)+EijI(26,4)-EijI(27,4)-EijI(2
     -   8,4)+EijI(33,4)-2*(EijI(4,3)+EijI(10,2)-EijI(15,3)+EijI(15,4
     -   )-EijI(17,3)+EijI(19,3)-EijI(32,4)-EijI(34,4)+EijI(35,4)))
       F(127)=DCMPLX(FR(127),FI(127))
       FR(128) = 16*(EijR(4,2)+EijR(4,3)-EijR(4,4)+EijR(10,2)+EijR(14,
     -   4)+EijR(15,4)-2*(EijR(7,2)+EijR(16,3))+EijR(16,4)-EijR(18,3)
     -   +EijR(19,3)+EijR(20,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR
     -   (35,4))
       FI(128) = 16*(EijI(4,2)+EijI(4,3)-EijI(4,4)+EijI(10,2)+EijI(14,
     -   4)+EijI(15,4)-2*(EijI(7,2)+EijI(16,3))+EijI(16,4)-EijI(18,3)
     -   +EijI(19,3)+EijI(20,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI
     -   (35,4))
       F(128)=DCMPLX(FR(128),FI(128))
       FR(129) = -16*Is12*(Dij1345R(2,2)+Dij1345R(3,2)-Dij1345R(3,3)+D
     -   ij1345R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)-2*(Dij1345R(6,2)-Di
     -   j1345R(9,3)+Dij1345R(10,3))+s12*(EijR(3,2)+EijR(4,2)+EijR(4,
     -   4)+EijR(18,4)+4*(EijR(16,3)-EijR(19,3))+EijR(19,4)+EijR(22,4
     -   )-2*(EijR(4,3)+EijR(10,2)-EijR(11,3)+EijR(13,3)-EijR(14,3)+E
     -   ijR(14,4)+EijR(16,4)+EijR(25,4)+EijR(30,4)-2*EijR(33,4))))
       FI(129) = -16*Is12*(Dij1345I(2,2)+Dij1345I(3,2)-Dij1345I(3,3)+D
     -   ij1345I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)-2*(Dij1345I(6,2)-Di
     -   j1345I(9,3)+Dij1345I(10,3))+s12*(EijI(3,2)+EijI(4,2)+EijI(4,
     -   4)+EijI(18,4)+4*(EijI(16,3)-EijI(19,3))+EijI(19,4)+EijI(22,4
     -   )-2*(EijI(4,3)+EijI(10,2)-EijI(11,3)+EijI(13,3)-EijI(14,3)+E
     -   ijI(14,4)+EijI(16,4)+EijI(25,4)+EijI(30,4)-2*EijI(33,4))))
       F(129)=DCMPLX(FR(129),FI(129))
       FR(130) = 16*Is12*(Dij1345R(3,3)-Dij1345R(5,2)+Dij1345R(6,2)-Di
     -   j1345R(7,3)-Dij1345R(9,3)+Dij1345R(10,3)+s12*(EijR(4,3)-EijR
     -   (4,4)-EijR(7,2)-EijR(7,3)+EijR(10,2)+EijR(16,4)-EijR(19,4)+E
     -   ijR(25,4)+2*(EijR(14,4)-EijR(16,3)+EijR(19,3)-EijR(33,4))))
       FI(130) = 16*Is12*(Dij1345I(3,3)-Dij1345I(5,2)+Dij1345I(6,2)-Di
     -   j1345I(7,3)-Dij1345I(9,3)+Dij1345I(10,3)+s12*(EijI(4,3)-EijI
     -   (4,4)-EijI(7,2)-EijI(7,3)+EijI(10,2)+EijI(16,4)-EijI(19,4)+E
     -   ijI(25,4)+2*(EijI(14,4)-EijI(16,3)+EijI(19,3)-EijI(33,4))))
       F(130)=DCMPLX(FR(130),FI(130))
       FR(131) = 16*(EijR(3,1)-EijR(4,1)-EijR(4,2)-EijR(4,4)+EijR(6,2)
     -   +EijR(7,2)-EijR(10,2)-EijR(14,3)+EijR(14,4)+2*(EijR(4,3)-Eij
     -   R(15,3))+EijR(15,4)-EijR(16,3)+EijR(16,4)+EijR(18,3)+EijR(20
     -   ,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4))
       FI(131) = 16*(EijI(3,1)-EijI(4,1)-EijI(4,2)-EijI(4,4)+EijI(6,2)
     -   +EijI(7,2)-EijI(10,2)-EijI(14,3)+EijI(14,4)+2*(EijI(4,3)-Eij
     -   I(15,3))+EijI(15,4)-EijI(16,3)+EijI(16,4)+EijI(18,3)+EijI(20
     -   ,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4))
       F(131)=DCMPLX(FR(131),FI(131))
       FR(132) = 16*Is12*(Dij1345R(3,3)-Dij1345R(5,2)+Dij1345R(6,2)-Di
     -   j1345R(7,3)-Dij1345R(9,3)+Dij1345R(10,3)+s12*(-EijR(2,1)+Eij
     -   R(4,1)+EijR(4,3)-EijR(4,4)-EijR(5,2)-EijR(7,2)+EijR(9,2)+Eij
     -   R(10,2)+EijR(14,4)-EijR(15,3)+EijR(15,4)-EijR(16,3)+EijR(16,
     -   4)+EijR(20,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4)))
       FI(132) = 16*Is12*(Dij1345I(3,3)-Dij1345I(5,2)+Dij1345I(6,2)-Di
     -   j1345I(7,3)-Dij1345I(9,3)+Dij1345I(10,3)+s12*(-EijI(2,1)+Eij
     -   I(4,1)+EijI(4,3)-EijI(4,4)-EijI(5,2)-EijI(7,2)+EijI(9,2)+Eij
     -   I(10,2)+EijI(14,4)-EijI(15,3)+EijI(15,4)-EijI(16,3)+EijI(16,
     -   4)+EijI(20,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4)))
       F(132)=DCMPLX(FR(132),FI(132))
       FR(133) = -16*Is12*(Dij1345R(2,1)+Dij1345R(2,2)-Dij1345R(3,1)-D
     -   ij1345R(3,3)+Dij1345R(4,2)-Dij1345R(5,2)-3*Dij1345R(6,2)+Dij
     -   1345R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)+2*(Dij1345R(3,2)+Dij1
     -   345R(9,3)-Dij1345R(10,3))+s12*(EijR(3,1)+EijR(3,2)-EijR(4,1)
     -   +3*(EijR(4,2)-EijR(4,3))+EijR(4,4)+EijR(6,3)-EijR(7,3)-4*(Ei
     -   jR(10,2)-EijR(14,3))+5*EijR(16,3)+EijR(18,4)-6*EijR(19,3)+Ei
     -   jR(19,4)+EijR(22,4)+2*(EijR(6,2)-EijR(7,2)+EijR(11,3)-EijR(1
     -   3,3)-EijR(14,4)-EijR(16,4)-EijR(25,4)-EijR(30,4)+2*EijR(33,4
     -   ))))
       FI(133) = -16*Is12*(Dij1345I(2,1)+Dij1345I(2,2)-Dij1345I(3,1)-D
     -   ij1345I(3,3)+Dij1345I(4,2)-Dij1345I(5,2)-3*Dij1345I(6,2)+Dij
     -   1345I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)+2*(Dij1345I(3,2)+Dij1
     -   345I(9,3)-Dij1345I(10,3))+s12*(EijI(3,1)+EijI(3,2)-EijI(4,1)
     -   +3*(EijI(4,2)-EijI(4,3))+EijI(4,4)+EijI(6,3)-EijI(7,3)-4*(Ei
     -   jI(10,2)-EijI(14,3))+5*EijI(16,3)+EijI(18,4)-6*EijI(19,3)+Ei
     -   jI(19,4)+EijI(22,4)+2*(EijI(6,2)-EijI(7,2)+EijI(11,3)-EijI(1
     -   3,3)-EijI(14,4)-EijI(16,4)-EijI(25,4)-EijI(30,4)+2*EijI(33,4
     -   ))))
       F(133)=DCMPLX(FR(133),FI(133))
       FR(134) = -16*(EijR(3,2)+EijR(4,2)+EijR(4,4)+4*(EijR(16,3)-EijR
     -   (20,3))+EijR(20,4)+EijR(21,4)+EijR(22,4)-2*(EijR(4,3)+EijR(1
     -   0,2)-EijR(12,3)+EijR(13,3)-EijR(15,3)+EijR(15,4)+EijR(16,4)+
     -   EijR(28,4)+EijR(31,4)-2*EijR(34,4)))
       FI(134) = -16*(EijI(3,2)+EijI(4,2)+EijI(4,4)+4*(EijI(16,3)-EijI
     -   (20,3))+EijI(20,4)+EijI(21,4)+EijI(22,4)-2*(EijI(4,3)+EijI(1
     -   0,2)-EijI(12,3)+EijI(13,3)-EijI(15,3)+EijI(15,4)+EijI(16,4)+
     -   EijI(28,4)+EijI(31,4)-2*EijI(34,4)))
       F(134)=DCMPLX(FR(134),FI(134))
       FR(135) = 16*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,2)+Di
     -   j1345R(3,3)+Dij1345R(4,2)-Dij1345R(7,3)-Dij1345R(9,3)+Dij134
     -   5R(10,3)-s12*(-EijR(3,1)+EijR(4,1)+EijR(4,4)-EijR(6,3)+EijR(
     -   7,2)+EijR(10,2)+EijR(16,3)-EijR(16,4)+EijR(19,4)-EijR(25,4)-
     -   2*(EijR(4,3)+EijR(6,2)-EijR(14,3)+EijR(14,4)-EijR(33,4))))
       FI(135) = 16*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,2)+Di
     -   j1345I(3,3)+Dij1345I(4,2)-Dij1345I(7,3)-Dij1345I(9,3)+Dij134
     -   5I(10,3)-s12*(-EijI(3,1)+EijI(4,1)+EijI(4,4)-EijI(6,3)+EijI(
     -   7,2)+EijI(10,2)+EijI(16,3)-EijI(16,4)+EijI(19,4)-EijI(25,4)-
     -   2*(EijI(4,3)+EijI(6,2)-EijI(14,3)+EijI(14,4)-EijI(33,4))))
       F(135)=DCMPLX(FR(135),FI(135))
       FR(136) = 16*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,2)+Di
     -   j1345R(3,3)+Dij1345R(4,2)-Dij1345R(7,3)-Dij1345R(9,3)+Dij134
     -   5R(10,3)+s12*(EijR(2,1)+EijR(3,1)+EijR(4,3)-EijR(4,4)+EijR(5
     -   ,2)+EijR(6,2)+EijR(8,2)-2*(EijR(4,1)+EijR(9,2))-EijR(10,2)+E
     -   ijR(14,4)+EijR(15,4)-EijR(16,3)+EijR(16,4)+EijR(17,3)-EijR(1
     -   8,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4)))
       FI(136) = 16*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,2)+Di
     -   j1345I(3,3)+Dij1345I(4,2)-Dij1345I(7,3)-Dij1345I(9,3)+Dij134
     -   5I(10,3)+s12*(EijI(2,1)+EijI(3,1)+EijI(4,3)-EijI(4,4)+EijI(5
     -   ,2)+EijI(6,2)+EijI(8,2)-2*(EijI(4,1)+EijI(9,2))-EijI(10,2)+E
     -   ijI(14,4)+EijI(15,4)-EijI(16,3)+EijI(16,4)+EijI(17,3)-EijI(1
     -   8,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4)))
       F(136)=DCMPLX(FR(136),FI(136))
       FR(137) = 16*(EijR(2,1)+EijR(2,2)-EijR(4,1)+EijR(4,2)+EijR(4,3)
     -   -EijR(4,4)+EijR(8,2)-3*EijR(9,2)+EijR(9,3)-EijR(10,3)-EijR(1
     -   6,3)+EijR(16,4)-EijR(21,4)+EijR(28,4)+2*(EijR(15,4)-EijR(34,
     -   4)))
       FI(137) = 16*(EijI(2,1)+EijI(2,2)-EijI(4,1)+EijI(4,2)+EijI(4,3)
     -   -EijI(4,4)+EijI(8,2)-3*EijI(9,2)+EijI(9,3)-EijI(10,3)-EijI(1
     -   6,3)+EijI(16,4)-EijI(21,4)+EijI(28,4)+2*(EijI(15,4)-EijI(34,
     -   4)))
       F(137)=DCMPLX(FR(137),FI(137))
       P(424) = p2sq-s15+s34
       FR(138) = -8*(D02345R+Dij2345R(2,1)+s15*EijR(7,2)-p2sq*(EijR(6,
     -   2)+EijR(7,2)-EijR(8,2))+s12*(EijR(2,1)-EijR(3,2)-EijR(4,1)+E
     -   ijR(6,2)+EijR(7,2)+EijR(8,2)-2*EijR(9,2))-2*(EijR(11,2)+EijR
     -   (24,3)-EijR(39,4)+EijR(45,4))-(EijR(3,1)+EijR(3,2)-EijR(4,1)
     -   )*P(30)+EijR(10,2)*P(95)+EijR(9,2)*P(424))
       FI(138) = -8*(D02345I+Dij2345I(2,1)+s15*EijI(7,2)-p2sq*(EijI(6,
     -   2)+EijI(7,2)-EijI(8,2))+s12*(EijI(2,1)-EijI(3,2)-EijI(4,1)+E
     -   ijI(6,2)+EijI(7,2)+EijI(8,2)-2*EijI(9,2))-2*(EijI(11,2)+EijI
     -   (24,3)-EijI(39,4)+EijI(45,4))-(EijI(3,1)+EijI(3,2)-EijI(4,1)
     -   )*P(30)+EijI(10,2)*P(95)+EijI(9,2)*P(424))
       F(138)=DCMPLX(FR(138),FI(138))
       P(425) = -1+Is45*p3sq
       P(426) = p2sq-s12-s15+s23+s45-2*P(404)
       FR(139) = 8*(-D02345R+Dij2345R(1,1)-Dij2345R(3,1)+2*(p3sq*EijR(
     -   3,2)+4*EijR(11,2)-EijR(22,3)-EijR(23,3)+2*EijR(24,3)-EijR(39
     -   ,4)-EijR(43,4)+EijR(44,4)+EijR(45,4))+EijR(6,2)*P(45)+Is45*(
     -   Cij123R(2,1)+Cij123R(2,2)-p3sq*(Dij1234R(3,1)+Dij1234R(6,2))
     -   +2*(Dij1234R(7,2)+Dij1234R(12,3))+Dij1234R(2,2)*P(9)-Dij1234
     -   R(4,2)*P(45))+EijR(2,2)*P(90)+EijR(4,2)*P(174)-(EijR(2,1)-Ei
     -   jR(4,1))*P(230)-EijR(5,2)*P(231)+EijR(7,2)*P(232)-EijR(10,2)
     -   *P(407)-EijR(8,2)*P(413)+Dij1234R(2,1)*P(425)-EijR(9,2)*P(42
     -   6))
       FI(139) = 8*(-D02345I+Dij2345I(1,1)-Dij2345I(3,1)+2*(p3sq*EijI(
     -   3,2)+4*EijI(11,2)-EijI(22,3)-EijI(23,3)+2*EijI(24,3)-EijI(39
     -   ,4)-EijI(43,4)+EijI(44,4)+EijI(45,4))+EijI(6,2)*P(45)+Is45*(
     -   Cij123I(2,1)+Cij123I(2,2)-p3sq*(Dij1234I(3,1)+Dij1234I(6,2))
     -   +2*(Dij1234I(7,2)+Dij1234I(12,3))+Dij1234I(2,2)*P(9)-Dij1234
     -   I(4,2)*P(45))+EijI(2,2)*P(90)+EijI(4,2)*P(174)-(EijI(2,1)-Ei
     -   jI(4,1))*P(230)-EijI(5,2)*P(231)+EijI(7,2)*P(232)-EijI(10,2)
     -   *P(407)-EijI(8,2)*P(413)+Dij1234I(2,1)*P(425)-EijI(9,2)*P(42
     -   6))
       F(139)=DCMPLX(FR(139),FI(139))
       FR(140) = 8*(-2*(EijR(23,3)-EijR(24,3)+EijR(38,4)+EijR(39,4)-2*
     -   EijR(45,4))+Is45*(C0134R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,
     -   2)-2*(Dij1234R(7,2)-Dij1234R(13,3))+Dij1234R(6,2)*P(9)-Dij12
     -   34R(5,2)*P(45)-Dij1234R(3,1)*P(214))-(EijR(6,2)-EijR(7,2))*P
     -   (221)-(-EijR(8,2)+EijR(9,2))*P(222)-(EijR(3,1)-EijR(4,1))*P(
     -   230))
       FI(140) = 8*(-2*(EijI(23,3)-EijI(24,3)+EijI(38,4)+EijI(39,4)-2*
     -   EijI(45,4))+Is45*(C0134I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,
     -   2)-2*(Dij1234I(7,2)-Dij1234I(13,3))+Dij1234I(6,2)*P(9)-Dij12
     -   34I(5,2)*P(45)-Dij1234I(3,1)*P(214))-(EijI(6,2)-EijI(7,2))*P
     -   (221)-(-EijI(8,2)+EijI(9,2))*P(222)-(EijI(3,1)-EijI(4,1))*P(
     -   230))
       F(140)=DCMPLX(FR(140),FI(140))
       P(427) = p2sq-p3sq-2*s12+s45
       P(428) = p2sq-p3sq-s12-s15+s34+s45
       FR(141) = 4*Is12s45*(-3*B013R-B014R+2*(p2sq*(C0123R+Cij123R(1,1
     -   ))+p4sq*Cij145R(2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,
     -   2))+C0145R*P(15)-Cij123R(2,1)*P(18)+C0134R*P(42)+Cij134R(2,1
     -   )*P(147)))-8*(-D02345R+Is12*(Cij145R(1,1)+2*(Dij1345R(7,2)+D
     -   ij1345R(12,3)-Dij1345R(13,3)))-Dij2345R(2,1)-s15*(EE0R+2*Eij
     -   R(1,1)+EijR(1,2))-s12*EijR(3,1)+2*(3*EijR(11,2)+EijR(21,3)+E
     -   ijR(23,3)-2*EijR(24,3)+EijR(39,4)+EijR(41,4)-EijR(42,4)-EijR
     -   (45,4))+EijR(3,2)*P(25)-EijR(8,2)*P(68)-(D01345R+Dij1345R(1,
     -   1))*P(92)+EijR(9,2)*P(222)-EijR(4,1)*P(230)+EijR(10,2)*P(261
     -   )-Is45*(Cij123R(3,2)-Dij1234R(3,2)*P(1)-2*(Cij134R(1,1)-p2sq
     -   *Dij1234R(6,2)-Dij1234R(11,3)-Dij1234R(1,1)*P(36))-Dij1234R(
     -   3,1)*P(42)+Dij1234R(2,1)*P(210)-Dij1234R(4,2)*P(262)+Dij1234
     -   R(5,2)*P(263))+EijR(2,1)*P(265)+EijR(5,2)*P(268)-(D01234R+Di
     -   j1234R(1,2))*P(274)+EijR(6,2)*P(427)-EijR(7,2)*P(428))
       FI(141) = 4*Is12s45*(-3*B013I-B014I+2*(p2sq*(C0123I+Cij123I(1,1
     -   ))+p4sq*Cij145I(2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,
     -   2))+C0145I*P(15)-Cij123I(2,1)*P(18)+C0134I*P(42)+Cij134I(2,1
     -   )*P(147)))-8*(-D02345I+Is12*(Cij145I(1,1)+2*(Dij1345I(7,2)+D
     -   ij1345I(12,3)-Dij1345I(13,3)))-Dij2345I(2,1)-s15*(EE0I+2*Eij
     -   I(1,1)+EijI(1,2))-s12*EijI(3,1)+2*(3*EijI(11,2)+EijI(21,3)+E
     -   ijI(23,3)-2*EijI(24,3)+EijI(39,4)+EijI(41,4)-EijI(42,4)-EijI
     -   (45,4))+EijI(3,2)*P(25)-EijI(8,2)*P(68)-(D01345I+Dij1345I(1,
     -   1))*P(92)+EijI(9,2)*P(222)-EijI(4,1)*P(230)+EijI(10,2)*P(261
     -   )-Is45*(Cij123I(3,2)-Dij1234I(3,2)*P(1)-2*(Cij134I(1,1)-p2sq
     -   *Dij1234I(6,2)-Dij1234I(11,3)-Dij1234I(1,1)*P(36))-Dij1234I(
     -   3,1)*P(42)+Dij1234I(2,1)*P(210)-Dij1234I(4,2)*P(262)+Dij1234
     -   I(5,2)*P(263))+EijI(2,1)*P(265)+EijI(5,2)*P(268)-(D01234I+Di
     -   j1234I(1,2))*P(274)+EijI(6,2)*P(427)-EijI(7,2)*P(428))
       F(141)=DCMPLX(FR(141),FI(141))
       FR(142) = 8*(Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R(3,2)+Dij1345R
     -   (4,2)-Dij1345R(5,2)-Dij1345R(6,2)-2*(EijR(38,4)+EijR(39,4)-2
     -   *EijR(45,4))+(-EijR(3,1)+EijR(4,1)-EijR(6,2)+EijR(7,2))*P(18
     -   )+(-EijR(4,2)+EijR(10,2))*P(222))
       FI(142) = 8*(Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I(3,2)+Dij1345I
     -   (4,2)-Dij1345I(5,2)-Dij1345I(6,2)-2*(EijI(38,4)+EijI(39,4)-2
     -   *EijI(45,4))+(-EijI(3,1)+EijI(4,1)-EijI(6,2)+EijI(7,2))*P(18
     -   )+(-EijI(4,2)+EijI(10,2))*P(222))
       F(142)=DCMPLX(FR(142),FI(142))
       P(429) = p2sq+s34
       P(430) = -p3sq-3*s12-s15+s45+2*P(429)
       FR(143) = 8*(D02345R+Dij1345R(1,1)+Dij1345R(1,2)-Dij1345R(3,1)+
     -   Dij1345R(3,2)+Dij2345R(2,1)-p2sq*EijR(2,1)+s12*EijR(3,1)-2*(
     -   Dij1345R(5,2)+EijR(11,2)-EijR(21,3)+EijR(23,3)+EijR(39,4)+Ei
     -   jR(41,4)-EijR(42,4)-EijR(45,4))+EijR(4,1)*P(18)-EijR(3,2)*P(
     -   25)+(-EijR(5,2)+EijR(8,2))*P(68)-EijR(4,2)*P(222)-EijR(10,2)
     -   *P(261)-EijR(6,2)*P(427)+EijR(7,2)*P(430))
       FI(143) = 8*(D02345I+Dij1345I(1,1)+Dij1345I(1,2)-Dij1345I(3,1)+
     -   Dij1345I(3,2)+Dij2345I(2,1)-p2sq*EijI(2,1)+s12*EijI(3,1)-2*(
     -   Dij1345I(5,2)+EijI(11,2)-EijI(21,3)+EijI(23,3)+EijI(39,4)+Ei
     -   jI(41,4)-EijI(42,4)-EijI(45,4))+EijI(4,1)*P(18)-EijI(3,2)*P(
     -   25)+(-EijI(5,2)+EijI(8,2))*P(68)-EijI(4,2)*P(222)-EijI(10,2)
     -   *P(261)-EijI(6,2)*P(427)+EijI(7,2)*P(430))
       F(143)=DCMPLX(FR(143),FI(143))
       P(431) = -p4sq+s12+s45
       P(432) = 2*p4sq+s12
       P(433) = -p2sq+p4sq+2*s12+s15-s34+s45
       P(434) = p2sq-s12+s23-s45
       P(435) = p2sq+2*p3sq-3*s12-s15-s23-s45
       FR(144) = 4*Is12s45*(3*B013R+B014R+2*(p2sq*Cij123R(2,1)-p4sq*Ci
     -   j145R(2,1)-2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2))-C0145R
     -   *P(15)-(C0123R+Cij123R(1,1))*P(18)-Cij134R(2,1)*P(30)-C0134R
     -   *P(43)))+8*(-D02345R+Is45*(Cij134R(1,1)+2*(Dij1234R(7,2)+Dij
     -   1234R(12,3)-Dij1234R(13,3)))+Dij1345R(1,2)+Dij2345R(1,1)-Dij
     -   2345R(3,1)-s15*(EE0R+EijR(1,1))+2*(p3sq*EijR(3,2)+5*EijR(11,
     -   2)-EijR(39,4)-EijR(43,4)+EijR(44,4)+EijR(45,4))+(-EijR(2,1)+
     -   EijR(4,1)-2*EijR(5,2))*P(18)+D01345R*P(39)+(D01234R+Dij1234R
     -   (1,1))*P(40)+EijR(6,2)*P(45)+EijR(2,2)*P(68)-EijR(10,2)*P(40
     -   7)-EijR(8,2)*P(413)+Is12*(Cij145R(1,1)+6*Dij1345R(7,2)-Dij13
     -   45R(5,2)*P(96)+2*Dij1345R(1,1)*P(230)+Dij1345R(4,2)*P(261)+D
     -   ij1345R(2,1)*P(284)+Dij1345R(2,2)*P(287)-Dij1345R(6,2)*P(293
     -   )-Dij1345R(3,1)*P(431)+Dij1345R(3,2)*P(432))+EijR(4,2)*P(433
     -   )+EijR(7,2)*P(434)+EijR(9,2)*P(435))
       FI(144) = 4*Is12s45*(3*B013I+B014I+2*(p2sq*Cij123I(2,1)-p4sq*Ci
     -   j145I(2,1)-2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2))-C0145I
     -   *P(15)-(C0123I+Cij123I(1,1))*P(18)-Cij134I(2,1)*P(30)-C0134I
     -   *P(43)))+8*(-D02345I+Is45*(Cij134I(1,1)+2*(Dij1234I(7,2)+Dij
     -   1234I(12,3)-Dij1234I(13,3)))+Dij1345I(1,2)+Dij2345I(1,1)-Dij
     -   2345I(3,1)-s15*(EE0I+EijI(1,1))+2*(p3sq*EijI(3,2)+5*EijI(11,
     -   2)-EijI(39,4)-EijI(43,4)+EijI(44,4)+EijI(45,4))+(-EijI(2,1)+
     -   EijI(4,1)-2*EijI(5,2))*P(18)+D01345I*P(39)+(D01234I+Dij1234I
     -   (1,1))*P(40)+EijI(6,2)*P(45)+EijI(2,2)*P(68)-EijI(10,2)*P(40
     -   7)-EijI(8,2)*P(413)+Is12*(Cij145I(1,1)+6*Dij1345I(7,2)-Dij13
     -   45I(5,2)*P(96)+2*Dij1345I(1,1)*P(230)+Dij1345I(4,2)*P(261)+D
     -   ij1345I(2,1)*P(284)+Dij1345I(2,2)*P(287)-Dij1345I(6,2)*P(293
     -   )-Dij1345I(3,1)*P(431)+Dij1345I(3,2)*P(432))+EijI(4,2)*P(433
     -   )+EijI(7,2)*P(434)+EijI(9,2)*P(435))
       F(144)=DCMPLX(FR(144),FI(144))
       FR(145) = 8*(D01345R-D02345R+Dij1345R(1,1)+Dij1345R(3,2)-Dij134
     -   5R(5,2)-Dij2345R(2,1)-s12*EijR(2,1)+2*(EijR(11,2)+EijR(23,3)
     -   -EijR(24,3)-EijR(39,4)+EijR(45,4))+(EijR(6,2)+EijR(7,2))*P(1
     -   8)+EijR(3,2)*P(25)+EijR(3,1)*P(30)-EijR(8,2)*P(68)-EijR(10,2
     -   )*P(95)+EijR(4,1)*P(144)-EijR(4,2)*P(222))
       FI(145) = 8*(D01345I-D02345I+Dij1345I(1,1)+Dij1345I(3,2)-Dij134
     -   5I(5,2)-Dij2345I(2,1)-s12*EijI(2,1)+2*(EijI(11,2)+EijI(23,3)
     -   -EijI(24,3)-EijI(39,4)+EijI(45,4))+(EijI(6,2)+EijI(7,2))*P(1
     -   8)+EijI(3,2)*P(25)+EijI(3,1)*P(30)-EijI(8,2)*P(68)-EijI(10,2
     -   )*P(95)+EijI(4,1)*P(144)-EijI(4,2)*P(222))
       F(145)=DCMPLX(FR(145),FI(145))
       P(436) = -p3sq+s12+s34+s45
       P(437) = -p3sq-3*p4sq+s12+s34
       P(438) = -2*p3sq-3*s12-s15+s23+s34+s45
       P(439) = p2sq+2*p3sq+p4sq-3*s12-s15+s45
       FR(146) = 4*Is12s45*(-3*B013R-B014R+4*(Cij123R(4,2)+Cij134R(4,2
     -   )+Cij145R(4,2))-2*(p2sq*Cij123R(2,1)-p4sq*Cij145R(2,1)-C0145
     -   R*P(15)-(C0123R+Cij123R(1,1))*P(18)-C0134R*P(43)+Cij134R(2,1
     -   )*P(61)))+8*(D02345R-Dij1345R(1,2)-Dij2345R(1,1)+Dij2345R(3,
     -   1)+s15*(EE0R+EijR(1,1))-2*(p3sq*EijR(3,2)+5*EijR(11,2)+EijR(
     -   23,3)-EijR(24,3)+EijR(39,4)+EijR(43,4)-EijR(44,4)-EijR(45,4)
     -   )+(EijR(2,1)+2*EijR(5,2))*P(18)+Is45*(Cij134R(1,2)-Cij134R(3
     -   ,2)-p2sq*Dij1234R(2,2)+s12*(-Dij1234R(2,1)+Dij1234R(3,1)+Dij
     -   1234R(6,2))+2*Dij1234R(12,3)+Dij1234R(4,2)*P(18))-EijR(4,1)*
     -   P(28)-EijR(6,2)*P(45)-EijR(2,2)*P(68)+(D01345R+2*Dij1345R(1,
     -   1))*P(92)-EijR(4,2)*P(174)-EijR(3,1)*P(222)+(D01234R+Dij1234
     -   R(1,1))*P(274)+EijR(9,2)*P(411)-EijR(7,2)*P(434)-Is12*(Cij14
     -   5R(1,1)+2*p4sq*Dij1345R(3,2)+6*Dij1345R(7,2)+Dij1345R(3,1)*P
     -   (15)-Dij1345R(4,2)*P(112)-Dij1345R(2,1)*P(235)+Dij1345R(2,2)
     -   *P(287)-Dij1345R(5,2)*P(436)+Dij1345R(6,2)*P(437))-EijR(8,2)
     -   *P(438)+EijR(10,2)*P(439))
       FI(146) = 4*Is12s45*(-3*B013I-B014I+4*(Cij123I(4,2)+Cij134I(4,2
     -   )+Cij145I(4,2))-2*(p2sq*Cij123I(2,1)-p4sq*Cij145I(2,1)-C0145
     -   I*P(15)-(C0123I+Cij123I(1,1))*P(18)-C0134I*P(43)+Cij134I(2,1
     -   )*P(61)))+8*(D02345I-Dij1345I(1,2)-Dij2345I(1,1)+Dij2345I(3,
     -   1)+s15*(EE0I+EijI(1,1))-2*(p3sq*EijI(3,2)+5*EijI(11,2)+EijI(
     -   23,3)-EijI(24,3)+EijI(39,4)+EijI(43,4)-EijI(44,4)-EijI(45,4)
     -   )+(EijI(2,1)+2*EijI(5,2))*P(18)+Is45*(Cij134I(1,2)-Cij134I(3
     -   ,2)-p2sq*Dij1234I(2,2)+s12*(-Dij1234I(2,1)+Dij1234I(3,1)+Dij
     -   1234I(6,2))+2*Dij1234I(12,3)+Dij1234I(4,2)*P(18))-EijI(4,1)*
     -   P(28)-EijI(6,2)*P(45)-EijI(2,2)*P(68)+(D01345I+2*Dij1345I(1,
     -   1))*P(92)-EijI(4,2)*P(174)-EijI(3,1)*P(222)+(D01234I+Dij1234
     -   I(1,1))*P(274)+EijI(9,2)*P(411)-EijI(7,2)*P(434)-Is12*(Cij14
     -   5I(1,1)+2*p4sq*Dij1345I(3,2)+6*Dij1345I(7,2)+Dij1345I(3,1)*P
     -   (15)-Dij1345I(4,2)*P(112)-Dij1345I(2,1)*P(235)+Dij1345I(2,2)
     -   *P(287)-Dij1345I(5,2)*P(436)+Dij1345I(6,2)*P(437))-EijI(8,2)
     -   *P(438)+EijI(10,2)*P(439))
       F(146)=DCMPLX(FR(146),FI(146))
       FR(147) = 8*(Dij1345R(2,2)-Dij1345R(4,2)+Dij1345R(5,2)-Dij1345R
     -   (6,2)-2*(EijR(23,3)-EijR(24,3)+EijR(38,4)+EijR(39,4)-2*EijR(
     -   45,4))-EijR(7,2)*P(18)+Is45*(Cij134R(1,1)-Cij134R(2,1)-Cij13
     -   4R(2,2)+Cij134R(3,2)+s12*(Dij1234R(3,2)-Dij1234R(5,2))+p2sq*
     -   (-Dij1234R(2,1)+Dij1234R(3,1)+Dij1234R(5,2)-Dij1234R(6,2))+2
     -   *(Dij1234R(7,2)+Dij1234R(13,3))+s45*EijR(6,2)*P(18))+(EijR(3
     -   ,1)-EijR(4,1))*P(28)+(-EijR(3,2)+EijR(10,2))*P(222))
       FI(147) = 8*(Dij1345I(2,2)-Dij1345I(4,2)+Dij1345I(5,2)-Dij1345I
     -   (6,2)-2*(EijI(23,3)-EijI(24,3)+EijI(38,4)+EijI(39,4)-2*EijI(
     -   45,4))-EijI(7,2)*P(18)+Is45*(Cij134I(1,1)-Cij134I(2,1)-Cij13
     -   4I(2,2)+Cij134I(3,2)+s12*(Dij1234I(3,2)-Dij1234I(5,2))+p2sq*
     -   (-Dij1234I(2,1)+Dij1234I(3,1)+Dij1234I(5,2)-Dij1234I(6,2))+2
     -   *(Dij1234I(7,2)+Dij1234I(13,3))+s45*EijI(6,2)*P(18))+(EijI(3
     -   ,1)-EijI(4,1))*P(28)+(-EijI(3,2)+EijI(10,2))*P(222))
       F(147)=DCMPLX(FR(147),FI(147))
       P(440) = p3sq-s15+s34-s45
       P(441) = p2sq-p3sq-s12+s34+s45
       P(442) = p2sq-p3sq-s15+s45-2*P(230)
       FR(148) = 8*(-D02345R-Dij1345R(1,1)-Dij1345R(1,2)+Dij1345R(2,1)
     -   +Dij1345R(4,2)+Dij1345R(5,2)-Dij1345R(6,2)-Dij2345R(2,1)+p2s
     -   q*EijR(2,1)+2*(EijR(11,2)-EijR(21,3)-EijR(23,3)+2*EijR(24,3)
     -   -EijR(39,4)-EijR(41,4)+EijR(42,4)+EijR(45,4))+Is45*(Cij134R(
     -   1,2)-Cij134R(3,2)+p2sq*Dij1234R(4,2)+4*Dij1234R(7,2)+2*(Cij1
     -   34R(1,1)-Cij134R(2,1)-p2sq*Dij1234R(6,2)+Dij1234R(11,3))+Dij
     -   1234R(3,2)*P(1)-Dij1234R(5,2)*P(6))-EijR(3,1)*P(8)+EijR(3,2)
     -   *P(25)-EijR(4,1)*P(28)-(-EijR(5,2)+EijR(8,2))*P(68)-EijR(6,2
     -   )*P(440)-EijR(7,2)*P(441)+EijR(10,2)*P(442))
       FI(148) = 8*(-D02345I-Dij1345I(1,1)-Dij1345I(1,2)+Dij1345I(2,1)
     -   +Dij1345I(4,2)+Dij1345I(5,2)-Dij1345I(6,2)-Dij2345I(2,1)+p2s
     -   q*EijI(2,1)+2*(EijI(11,2)-EijI(21,3)-EijI(23,3)+2*EijI(24,3)
     -   -EijI(39,4)-EijI(41,4)+EijI(42,4)+EijI(45,4))+Is45*(Cij134I(
     -   1,2)-Cij134I(3,2)+p2sq*Dij1234I(4,2)+4*Dij1234I(7,2)+2*(Cij1
     -   34I(1,1)-Cij134I(2,1)-p2sq*Dij1234I(6,2)+Dij1234I(11,3))+Dij
     -   1234I(3,2)*P(1)-Dij1234I(5,2)*P(6))-EijI(3,1)*P(8)+EijI(3,2)
     -   *P(25)-EijI(4,1)*P(28)-(-EijI(5,2)+EijI(8,2))*P(68)-EijI(6,2
     -   )*P(440)-EijI(7,2)*P(441)+EijI(10,2)*P(442))
       F(148)=DCMPLX(FR(148),FI(148))
       FR(149) = -8*(D01345R-D02345R+Dij1345R(1,1)-Dij1345R(5,2)+Dij13
     -   45R(6,2)-Dij2345R(2,1)-s12*EijR(2,1)+2*(EijR(11,2)+EijR(39,4
     -   )-EijR(45,4))+(EijR(6,2)+EijR(7,2))*P(18)-EijR(10,2)*P(24)+E
     -   ijR(3,2)*P(25)+EijR(3,1)*P(30)-EijR(8,2)*P(68)+EijR(4,1)*P(1
     -   44))
       FI(149) = -8*(D01345I-D02345I+Dij1345I(1,1)-Dij1345I(5,2)+Dij13
     -   45I(6,2)-Dij2345I(2,1)-s12*EijI(2,1)+2*(EijI(11,2)+EijI(39,4
     -   )-EijI(45,4))+(EijI(6,2)+EijI(7,2))*P(18)-EijI(10,2)*P(24)+E
     -   ijI(3,2)*P(25)+EijI(3,1)*P(30)-EijI(8,2)*P(68)+EijI(4,1)*P(1
     -   44))
       F(149)=DCMPLX(FR(149),FI(149))
       FR(150) = -16*Is12s45*(s45*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345
     -   R(3,3)+Dij1345R(4,3)-Dij1345R(5,3)+Dij1345R(9,3)+2*(Dij1345R
     -   (3,2)+Dij1345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3
     -   )-Dij1345R(10,3)))+s12*(-2*(Dij1234R(2,1)-Dij1234R(3,1))-Dij
     -   1234R(4,3)-3*(Dij1234R(4,2)-Dij1234R(5,2))+Dij1234R(5,3)+s45
     -   *(EijR(3,1)-EijR(4,1)+EijR(4,4)+EijR(5,2)+EijR(5,3)+EijR(6,2
     -   )-EijR(7,3)+EijR(8,2)+3*(EijR(4,2)-EijR(4,3)+EijR(15,3))-Eij
     -   R(15,4)-EijR(16,4)+4*(EijR(14,3)-EijR(18,3))+EijR(19,4)+EijR
     -   (23,4)-EijR(24,4)-EijR(25,4)+EijR(34,4)-2*(EijR(7,2)+EijR(9,
     -   2)+EijR(10,2)+EijR(14,4)-EijR(16,3)-EijR(17,3)+EijR(19,3)+Ei
     -   jR(20,3)-EijR(32,4)-EijR(33,4)+EijR(35,4)))))
       FI(150) = -16*Is12s45*(s45*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345
     -   I(3,3)+Dij1345I(4,3)-Dij1345I(5,3)+Dij1345I(9,3)+2*(Dij1345I
     -   (3,2)+Dij1345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3
     -   )-Dij1345I(10,3)))+s12*(-2*(Dij1234I(2,1)-Dij1234I(3,1))-Dij
     -   1234I(4,3)-3*(Dij1234I(4,2)-Dij1234I(5,2))+Dij1234I(5,3)+s45
     -   *(EijI(3,1)-EijI(4,1)+EijI(4,4)+EijI(5,2)+EijI(5,3)+EijI(6,2
     -   )-EijI(7,3)+EijI(8,2)+3*(EijI(4,2)-EijI(4,3)+EijI(15,3))-Eij
     -   I(15,4)-EijI(16,4)+4*(EijI(14,3)-EijI(18,3))+EijI(19,4)+EijI
     -   (23,4)-EijI(24,4)-EijI(25,4)+EijI(34,4)-2*(EijI(7,2)+EijI(9,
     -   2)+EijI(10,2)+EijI(14,4)-EijI(16,3)-EijI(17,3)+EijI(19,3)+Ei
     -   jI(20,3)-EijI(32,4)-EijI(33,4)+EijI(35,4)))))
       F(150)=DCMPLX(FR(150),FI(150))
       FR(151) = -16*Is12s45*(s45*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345
     -   R(3,3)+Dij1345R(4,3)-Dij1345R(5,3)+Dij1345R(9,3)+2*(Dij1345R
     -   (3,2)+Dij1345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3
     -   )-Dij1345R(10,3)))+s12*(-Dij1234R(2,2)+Dij1234R(5,2)-2*(Dij1
     -   234R(2,1)-Dij1234R(3,1)+Dij1234R(4,2)-Dij1234R(6,2))-Dij1234
     -   R(6,3)+Dij1234R(10,3)+s45*(EijR(3,1)-EijR(4,1)+EijR(4,4)+Eij
     -   R(9,3)-EijR(10,3)+EijR(14,3)-EijR(14,4)-EijR(16,4)+EijR(17,3
     -   )-EijR(18,3)-EijR(19,3)+3*(EijR(15,3)-EijR(20,3))+EijR(21,4)
     -   +EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(33,4)+2*(EijR(4,2)-Ei
     -   jR(4,3)+EijR(8,2)-EijR(9,2)-EijR(10,2)-EijR(15,4)+EijR(16,3)
     -   +EijR(32,4)+EijR(34,4)-EijR(35,4)))))
       FI(151) = -16*Is12s45*(s45*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345
     -   I(3,3)+Dij1345I(4,3)-Dij1345I(5,3)+Dij1345I(9,3)+2*(Dij1345I
     -   (3,2)+Dij1345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3
     -   )-Dij1345I(10,3)))+s12*(-Dij1234I(2,2)+Dij1234I(5,2)-2*(Dij1
     -   234I(2,1)-Dij1234I(3,1)+Dij1234I(4,2)-Dij1234I(6,2))-Dij1234
     -   I(6,3)+Dij1234I(10,3)+s45*(EijI(3,1)-EijI(4,1)+EijI(4,4)+Eij
     -   I(9,3)-EijI(10,3)+EijI(14,3)-EijI(14,4)-EijI(16,4)+EijI(17,3
     -   )-EijI(18,3)-EijI(19,3)+3*(EijI(15,3)-EijI(20,3))+EijI(21,4)
     -   +EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(33,4)+2*(EijI(4,2)-Ei
     -   jI(4,3)+EijI(8,2)-EijI(9,2)-EijI(10,2)-EijI(15,4)+EijI(16,3)
     -   +EijI(32,4)+EijI(34,4)-EijI(35,4)))))
       F(151)=DCMPLX(FR(151),FI(151))
       FR(152) = -16*Is45*(-Dij1234R(2,1)-2*Dij1234R(2,2)+Dij1234R(3,1
     -   )+Dij1234R(5,2)+Dij1234R(6,2)-Dij1234R(6,3)+Dij1234R(10,3)+s
     -   45*(EijR(4,4)+EijR(5,2)-EijR(7,2)+EijR(8,2)+EijR(8,3)+EijR(9
     -   ,3)-EijR(10,2)-EijR(14,4)+5*EijR(15,3)-EijR(16,4)+EijR(17,3)
     -   -EijR(19,3)-3*(EijR(4,3)+EijR(18,3)+EijR(20,3))+EijR(21,4)+E
     -   ijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(33,4)+2*(EijR(4,2)-EijR
     -   (9,2)-EijR(10,3)+EijR(14,3)-EijR(15,4)+EijR(16,3)+EijR(32,4)
     -   +EijR(34,4)-EijR(35,4))))
       FI(152) = -16*Is45*(-Dij1234I(2,1)-2*Dij1234I(2,2)+Dij1234I(3,1
     -   )+Dij1234I(5,2)+Dij1234I(6,2)-Dij1234I(6,3)+Dij1234I(10,3)+s
     -   45*(EijI(4,4)+EijI(5,2)-EijI(7,2)+EijI(8,2)+EijI(8,3)+EijI(9
     -   ,3)-EijI(10,2)-EijI(14,4)+5*EijI(15,3)-EijI(16,4)+EijI(17,3)
     -   -EijI(19,3)-3*(EijI(4,3)+EijI(18,3)+EijI(20,3))+EijI(21,4)+E
     -   ijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(33,4)+2*(EijI(4,2)-EijI
     -   (9,2)-EijI(10,3)+EijI(14,3)-EijI(15,4)+EijI(16,3)+EijI(32,4)
     -   +EijI(34,4)-EijI(35,4))))
       F(152)=DCMPLX(FR(152),FI(152))
       FR(153) = 16*Is45*(Dij1234R(2,1)+Dij1234R(2,3)-Dij1234R(3,1)+2*
     -   (Dij1234R(2,2)-Dij1234R(6,2))-Dij1234R(8,3)+s45*(-EijR(4,2)-
     -   EijR(4,4)-EijR(8,2)+EijR(9,2)-EijR(9,4)+EijR(10,2)+EijR(10,4
     -   )+2*(EijR(4,3)-EijR(9,3)+EijR(10,3)-EijR(16,3))+EijR(16,4)-4
     -   *(EijR(15,3)-EijR(20,3))+3*(EijR(15,4)-EijR(21,4)+EijR(28,4)
     -   -EijR(34,4))))
       FI(153) = 16*Is45*(Dij1234I(2,1)+Dij1234I(2,3)-Dij1234I(3,1)+2*
     -   (Dij1234I(2,2)-Dij1234I(6,2))-Dij1234I(8,3)+s45*(-EijI(4,2)-
     -   EijI(4,4)-EijI(8,2)+EijI(9,2)-EijI(9,4)+EijI(10,2)+EijI(10,4
     -   )+2*(EijI(4,3)-EijI(9,3)+EijI(10,3)-EijI(16,3))+EijI(16,4)-4
     -   *(EijI(15,3)-EijI(20,3))+3*(EijI(15,4)-EijI(21,4)+EijI(28,4)
     -   -EijI(34,4))))
       F(153)=DCMPLX(FR(153),FI(153))
       FR(154) = -16*Is12s45*(s45*(Dij1345R(2,2)+Dij1345R(3,2)-Dij1345
     -   R(3,3)+Dij1345R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)-2*(Dij1345R
     -   (6,2)-Dij1345R(9,3)+Dij1345R(10,3)))-s12*(Dij1234R(2,1)-Dij1
     -   234R(3,1)-2*Dij1234R(3,2)+Dij1234R(4,2)+Dij1234R(6,2)-Dij123
     -   4R(7,3)+Dij1234R(10,3)-s45*(EijR(3,2)+EijR(4,2)+EijR(4,4)-Ei
     -   jR(6,2)+EijR(7,2)+EijR(8,2)-EijR(9,2)+EijR(12,3)-EijR(13,3)+
     -   EijR(14,3)-EijR(14,4)-EijR(15,4)+EijR(17,3)-EijR(18,3)-EijR(
     -   19,3)+3*(EijR(16,3)-EijR(20,3))+EijR(22,4)+EijR(29,4)-EijR(3
     -   0,4)-EijR(31,4)+EijR(32,4)-2*(EijR(4,3)+EijR(10,2)-EijR(15,3
     -   )+EijR(16,4)-EijR(33,4)-EijR(34,4)+EijR(35,4)))))
       FI(154) = -16*Is12s45*(s45*(Dij1345I(2,2)+Dij1345I(3,2)-Dij1345
     -   I(3,3)+Dij1345I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)-2*(Dij1345I
     -   (6,2)-Dij1345I(9,3)+Dij1345I(10,3)))-s12*(Dij1234I(2,1)-Dij1
     -   234I(3,1)-2*Dij1234I(3,2)+Dij1234I(4,2)+Dij1234I(6,2)-Dij123
     -   4I(7,3)+Dij1234I(10,3)-s45*(EijI(3,2)+EijI(4,2)+EijI(4,4)-Ei
     -   jI(6,2)+EijI(7,2)+EijI(8,2)-EijI(9,2)+EijI(12,3)-EijI(13,3)+
     -   EijI(14,3)-EijI(14,4)-EijI(15,4)+EijI(17,3)-EijI(18,3)-EijI(
     -   19,3)+3*(EijI(16,3)-EijI(20,3))+EijI(22,4)+EijI(29,4)-EijI(3
     -   0,4)-EijI(31,4)+EijI(32,4)-2*(EijI(4,3)+EijI(10,2)-EijI(15,3
     -   )+EijI(16,4)-EijI(33,4)-EijI(34,4)+EijI(35,4)))))
       F(154)=DCMPLX(FR(154),FI(154))
       FR(155) = -16*Is45*(Dij1234R(3,2)+Dij1234R(5,2)-2*Dij1234R(6,2)
     -   +Dij1234R(7,3)-Dij1234R(10,3)+s45*(EijR(4,2)+EijR(4,4)+EijR(
     -   6,2)-EijR(7,2)-EijR(10,2)+EijR(12,3)-EijR(13,3)+EijR(14,3)-E
     -   ijR(14,4)-EijR(15,4)+EijR(17,3)-EijR(18,3)-EijR(19,3)+3*(Eij
     -   R(16,3)-EijR(20,3))+EijR(22,4)+EijR(29,4)-EijR(30,4)-EijR(31
     -   ,4)+EijR(32,4)-2*(EijR(4,3)-EijR(15,3)+EijR(16,4)-EijR(33,4)
     -   -EijR(34,4)+EijR(35,4))))
       FI(155) = -16*Is45*(Dij1234I(3,2)+Dij1234I(5,2)-2*Dij1234I(6,2)
     -   +Dij1234I(7,3)-Dij1234I(10,3)+s45*(EijI(4,2)+EijI(4,4)+EijI(
     -   6,2)-EijI(7,2)-EijI(10,2)+EijI(12,3)-EijI(13,3)+EijI(14,3)-E
     -   ijI(14,4)-EijI(15,4)+EijI(17,3)-EijI(18,3)-EijI(19,3)+3*(Eij
     -   I(16,3)-EijI(20,3))+EijI(22,4)+EijI(29,4)-EijI(30,4)-EijI(31
     -   ,4)+EijI(32,4)-2*(EijI(4,3)-EijI(15,3)+EijI(16,4)-EijI(33,4)
     -   -EijI(34,4)+EijI(35,4))))
       F(155)=DCMPLX(FR(155),FI(155))
       FR(156) = EijR(3,2)-EijR(4,3)+EijR(4,4)-EijR(6,2)+EijR(7,2)-Eij
     -   R(10,2)+EijR(11,3)+EijR(12,3)-EijR(14,4)-EijR(15,4)+3*EijR(1
     -   6,3)-EijR(17,3)+EijR(18,3)-EijR(19,3)-EijR(20,3)+EijR(22,4)+
     -   EijR(29,4)-EijR(30,4)-EijR(31,4)+EijR(32,4)-2*(EijR(13,3)+Ei
     -   jR(16,4)-EijR(33,4)-EijR(34,4)+EijR(35,4))
       FI(156) = EijI(3,2)-EijI(4,3)+EijI(4,4)-EijI(6,2)+EijI(7,2)-Eij
     -   I(10,2)+EijI(11,3)+EijI(12,3)-EijI(14,4)-EijI(15,4)+3*EijI(1
     -   6,3)-EijI(17,3)+EijI(18,3)-EijI(19,3)-EijI(20,3)+EijI(22,4)+
     -   EijI(29,4)-EijI(30,4)-EijI(31,4)+EijI(32,4)-2*(EijI(13,3)+Ei
     -   jI(16,4)-EijI(33,4)-EijI(34,4)+EijI(35,4))
       F(156)=DCMPLX(FR(156),FI(156))
       FR(157) = EijR(4,2)+EijR(4,4)+EijR(6,2)-EijR(7,2)-EijR(10,2)+Ei
     -   jR(11,3)+EijR(12,3)+EijR(14,3)-EijR(14,4)+EijR(15,3)-EijR(15
     -   ,4)+4*EijR(16,3)+EijR(22,4)+EijR(29,4)-EijR(30,4)-EijR(31,4)
     -   +EijR(32,4)-2*(EijR(4,3)+EijR(13,3)+EijR(16,4)+EijR(19,3)+Ei
     -   jR(20,3)-EijR(33,4)-EijR(34,4)+EijR(35,4))
       FI(157) = EijI(4,2)+EijI(4,4)+EijI(6,2)-EijI(7,2)-EijI(10,2)+Ei
     -   jI(11,3)+EijI(12,3)+EijI(14,3)-EijI(14,4)+EijI(15,3)-EijI(15
     -   ,4)+4*EijI(16,3)+EijI(22,4)+EijI(29,4)-EijI(30,4)-EijI(31,4)
     -   +EijI(32,4)-2*(EijI(4,3)+EijI(13,3)+EijI(16,4)+EijI(19,3)+Ei
     -   jI(20,3)-EijI(33,4)-EijI(34,4)+EijI(35,4))
       F(157)=DCMPLX(FR(157),FI(157))
       FR(158) = EijR(3,3)-EijR(4,3)+EijR(4,4)+EijR(11,4)-EijR(13,4)-E
     -   ijR(14,4)-3*(EijR(13,3)-EijR(16,3)+EijR(16,4)-EijR(22,4)+Eij
     -   R(30,4)-EijR(33,4))
       FI(158) = EijI(3,3)-EijI(4,3)+EijI(4,4)+EijI(11,4)-EijI(13,4)-E
     -   ijI(14,4)-3*(EijI(13,3)-EijI(16,3)+EijI(16,4)-EijI(22,4)+Eij
     -   I(30,4)-EijI(33,4))
       F(158)=DCMPLX(FR(158),FI(158))
       FR(159) = -EijR(4,3)+EijR(4,4)-EijR(13,3)-EijR(14,4)+EijR(22,4)
     -   -EijR(30,4)+2*(EijR(16,3)-EijR(16,4)+EijR(33,4))
       FI(159) = -EijI(4,3)+EijI(4,4)-EijI(13,3)-EijI(14,4)+EijI(22,4)
     -   -EijI(30,4)+2*(EijI(16,3)-EijI(16,4)+EijI(33,4))
       F(159)=DCMPLX(FR(159),FI(159))
       FR(160) = EijR(4,2)-EijR(4,3)+EijR(4,4)-EijR(7,2)-EijR(14,4)-Ei
     -   jR(15,4)+2*EijR(16,3)-EijR(16,4)+EijR(18,3)-EijR(19,3)-EijR(
     -   20,3)+EijR(32,4)+EijR(33,4)+EijR(34,4)-EijR(35,4)
       FI(160) = EijI(4,2)-EijI(4,3)+EijI(4,4)-EijI(7,2)-EijI(14,4)-Ei
     -   jI(15,4)+2*EijI(16,3)-EijI(16,4)+EijI(18,3)-EijI(19,3)-EijI(
     -   20,3)+EijI(32,4)+EijI(33,4)+EijI(34,4)-EijI(35,4)
       F(160)=DCMPLX(FR(160),FI(160))
       FR(161) = -EijR(4,4)+EijR(13,3)+EijR(14,3)+EijR(14,4)-EijR(16,3
     -   )-EijR(19,3)-EijR(22,4)+EijR(30,4)+2*(EijR(16,4)-EijR(33,4))
       FI(161) = -EijI(4,4)+EijI(13,3)+EijI(14,3)+EijI(14,4)-EijI(16,3
     -   )-EijI(19,3)-EijI(22,4)+EijI(30,4)+2*(EijI(16,4)-EijI(33,4))
       F(161)=DCMPLX(FR(161),FI(161))
       FR(162) = -EijR(4,4)+EijR(14,3)+EijR(14,4)-EijR(16,3)+EijR(16,4
     -   )-EijR(33,4)
       FI(162) = -EijI(4,4)+EijI(14,3)+EijI(14,4)-EijI(16,3)+EijI(16,4
     -   )-EijI(33,4)
       F(162)=DCMPLX(FR(162),FI(162))
       FR(163) = EijR(3,1)-EijR(4,1)+EijR(4,2)+EijR(4,3)-EijR(4,4)+Eij
     -   R(8,2)-EijR(9,2)-EijR(10,2)-EijR(15,3)-EijR(16,3)+EijR(16,4)
     -   +EijR(20,3)-EijR(21,4)+EijR(28,4)+2*(EijR(15,4)-EijR(34,4))
       FI(163) = EijI(3,1)-EijI(4,1)+EijI(4,2)+EijI(4,3)-EijI(4,4)+Eij
     -   I(8,2)-EijI(9,2)-EijI(10,2)-EijI(15,3)-EijI(16,3)+EijI(16,4)
     -   +EijI(20,3)-EijI(21,4)+EijI(28,4)+2*(EijI(15,4)-EijI(34,4))
       F(163)=DCMPLX(FR(163),FI(163))
       FR(164) = EijR(3,1)+EijR(3,2)-EijR(4,1)+EijR(4,3)-EijR(4,4)-Eij
     -   R(10,2)-EijR(15,3)+EijR(15,4)-EijR(16,3)+EijR(20,3)-EijR(22,
     -   4)+EijR(31,4)+2*(EijR(16,4)-EijR(34,4))
       FI(164) = EijI(3,1)+EijI(3,2)-EijI(4,1)+EijI(4,3)-EijI(4,4)-Eij
     -   I(10,2)-EijI(15,3)+EijI(15,4)-EijI(16,3)+EijI(20,3)-EijI(22,
     -   4)+EijI(31,4)+2*(EijI(16,4)-EijI(34,4))
       F(164)=DCMPLX(FR(164),FI(164))
       FR(165) = EijR(2,1)+EijR(2,2)-EijR(4,1)+EijR(4,2)-EijR(4,3)+Eij
     -   R(4,4)+EijR(15,3)+EijR(16,3)-EijR(16,4)-EijR(20,3)+EijR(21,4
     -   )-EijR(28,4)-2*(EijR(9,2)+EijR(15,4)-EijR(34,4))
       FI(165) = EijI(2,1)+EijI(2,2)-EijI(4,1)+EijI(4,2)-EijI(4,3)+Eij
     -   I(4,4)+EijI(15,3)+EijI(16,3)-EijI(16,4)-EijI(20,3)+EijI(21,4
     -   )-EijI(28,4)-2*(EijI(9,2)+EijI(15,4)-EijI(34,4))
       F(165)=DCMPLX(FR(165),FI(165))
       FR(166) = EijR(3,1)-EijR(4,1)+EijR(4,4)+EijR(8,2)-EijR(9,2)-Eij
     -   R(15,4)+EijR(22,4)-EijR(31,4)-2*(EijR(16,4)-EijR(34,4))
       FI(166) = EijI(3,1)-EijI(4,1)+EijI(4,4)+EijI(8,2)-EijI(9,2)-Eij
     -   I(15,4)+EijI(22,4)-EijI(31,4)-2*(EijI(16,4)-EijI(34,4))
       F(166)=DCMPLX(FR(166),FI(166))
       FR(167) = EijR(4,4)-EijR(9,2)+EijR(10,2)-EijR(15,4)-EijR(16,4)+
     -   EijR(34,4)
       FI(167) = EijI(4,4)-EijI(9,2)+EijI(10,2)-EijI(15,4)-EijI(16,4)+
     -   EijI(34,4)
       F(167)=DCMPLX(FR(167),FI(167))
       FR(168) = EijR(3,2)+EijR(4,4)+EijR(8,2)-EijR(9,2)+EijR(11,3)+Ei
     -   jR(12,3)-EijR(14,4)-EijR(15,4)+5*EijR(16,3)+EijR(17,3)-EijR(
     -   18,3)-3*(EijR(4,3)+EijR(10,2)+EijR(19,3)+EijR(20,3))+EijR(22
     -   ,4)+EijR(29,4)-EijR(30,4)-EijR(31,4)+EijR(32,4)+2*(EijR(4,2)
     -   -EijR(13,3)+EijR(14,3)+EijR(15,3)-EijR(16,4)+EijR(33,4)+EijR
     -   (34,4)-EijR(35,4))
       FI(168) = EijI(3,2)+EijI(4,4)+EijI(8,2)-EijI(9,2)+EijI(11,3)+Ei
     -   jI(12,3)-EijI(14,4)-EijI(15,4)+5*EijI(16,3)+EijI(17,3)-EijI(
     -   18,3)-3*(EijI(4,3)+EijI(10,2)+EijI(19,3)+EijI(20,3))+EijI(22
     -   ,4)+EijI(29,4)-EijI(30,4)-EijI(31,4)+EijI(32,4)+2*(EijI(4,2)
     -   -EijI(13,3)+EijI(14,3)+EijI(15,3)-EijI(16,4)+EijI(33,4)+EijI
     -   (34,4)-EijI(35,4))
       F(168)=DCMPLX(FR(168),FI(168))
       FR(169) = EijR(3,2)+EijR(3,3)+EijR(4,2)+EijR(4,4)+EijR(11,3)+Ei
     -   jR(11,4)-4*EijR(13,3)-EijR(13,4)+EijR(14,3)-EijR(14,4)+5*Eij
     -   R(16,3)-2*(EijR(4,3)+EijR(10,2)+EijR(19,3))-3*(EijR(16,4)-Ei
     -   jR(22,4)+EijR(30,4)-EijR(33,4))
       FI(169) = EijI(3,2)+EijI(3,3)+EijI(4,2)+EijI(4,4)+EijI(11,3)+Ei
     -   jI(11,4)-4*EijI(13,3)-EijI(13,4)+EijI(14,3)-EijI(14,4)+5*Eij
     -   I(16,3)-2*(EijI(4,3)+EijI(10,2)+EijI(19,3))-3*(EijI(16,4)-Ei
     -   jI(22,4)+EijI(30,4)-EijI(33,4))
       F(169)=DCMPLX(FR(169),FI(169))
       FR(170) = EijR(4,2)+EijR(4,4)-EijR(10,2)-EijR(13,3)+EijR(14,3)-
     -   EijR(14,4)+3*EijR(16,3)-EijR(19,3)+EijR(22,4)-EijR(30,4)-2*(
     -   EijR(4,3)+EijR(16,4)-EijR(33,4))
       FI(170) = EijI(4,2)+EijI(4,4)-EijI(10,2)-EijI(13,3)+EijI(14,3)-
     -   EijI(14,4)+3*EijI(16,3)-EijI(19,3)+EijI(22,4)-EijI(30,4)-2*(
     -   EijI(4,3)+EijI(16,4)-EijI(33,4))
       F(170)=DCMPLX(FR(170),FI(170))
       FR(171) = Is12*(Dij1345R(2,1)+Dij1345R(2,2)-Dij1345R(3,1)-Dij13
     -   45R(3,3)+Dij1345R(4,2)-Dij1345R(5,2)-3*Dij1345R(6,2)+Dij1345
     -   R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)+2*(Dij1345R(3,2)+Dij1345R
     -   (9,3)-Dij1345R(10,3))+s12*(EijR(3,1)+EijR(3,2)-EijR(4,1)+Eij
     -   R(4,4)+EijR(6,2)-EijR(7,2)-3*EijR(10,2)+EijR(11,3)+EijR(12,3
     -   )+EijR(14,3)-EijR(14,4)+EijR(15,3)-EijR(15,4)+4*EijR(16,3)+E
     -   ijR(22,4)+EijR(29,4)-EijR(30,4)-EijR(31,4)+EijR(32,4)+2*(Eij
     -   R(4,2)-EijR(4,3)-EijR(13,3)-EijR(16,4)-EijR(19,3)-EijR(20,3)
     -   +EijR(33,4)+EijR(34,4)-EijR(35,4))))
       FI(171) = Is12*(Dij1345I(2,1)+Dij1345I(2,2)-Dij1345I(3,1)-Dij13
     -   45I(3,3)+Dij1345I(4,2)-Dij1345I(5,2)-3*Dij1345I(6,2)+Dij1345
     -   I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)+2*(Dij1345I(3,2)+Dij1345I
     -   (9,3)-Dij1345I(10,3))+s12*(EijI(3,1)+EijI(3,2)-EijI(4,1)+Eij
     -   I(4,4)+EijI(6,2)-EijI(7,2)-3*EijI(10,2)+EijI(11,3)+EijI(12,3
     -   )+EijI(14,3)-EijI(14,4)+EijI(15,3)-EijI(15,4)+4*EijI(16,3)+E
     -   ijI(22,4)+EijI(29,4)-EijI(30,4)-EijI(31,4)+EijI(32,4)+2*(Eij
     -   I(4,2)-EijI(4,3)-EijI(13,3)-EijI(16,4)-EijI(19,3)-EijI(20,3)
     -   +EijI(33,4)+EijI(34,4)-EijI(35,4))))
       F(171)=DCMPLX(FR(171),FI(171))
       FR(172) = EijR(3,3)-EijR(4,3)+EijR(4,4)+EijR(12,4)-EijR(13,4)-E
     -   ijR(15,4)-3*(EijR(13,3)-EijR(16,3)+EijR(16,4)-EijR(22,4)+Eij
     -   R(31,4)-EijR(34,4))
       FI(172) = EijI(3,3)-EijI(4,3)+EijI(4,4)+EijI(12,4)-EijI(13,4)-E
     -   ijI(15,4)-3*(EijI(13,3)-EijI(16,3)+EijI(16,4)-EijI(22,4)+Eij
     -   I(31,4)-EijI(34,4))
       F(172)=DCMPLX(FR(172),FI(172))
       FR(173) = EijR(4,2)+EijR(4,3)-EijR(4,4)-EijR(10,2)+EijR(13,3)+E
     -   ijR(15,4)-EijR(22,4)+EijR(31,4)-2*(EijR(16,3)-EijR(16,4)+Eij
     -   R(34,4))
       FI(173) = EijI(4,2)+EijI(4,3)-EijI(4,4)-EijI(10,2)+EijI(13,3)+E
     -   ijI(15,4)-EijI(22,4)+EijI(31,4)-2*(EijI(16,3)-EijI(16,4)+Eij
     -   I(34,4))
       F(173)=DCMPLX(FR(173),FI(173))
       FR(174) = Is12*(Dij1345R(2,2)+Dij1345R(2,3)+Dij1345R(3,2)-Dij13
     -   45R(3,3)-2*Dij1345R(6,2)-3*(Dij1345R(8,3)-Dij1345R(9,3))+s12
     -   *(EijR(3,2)+EijR(3,3)+EijR(4,2)+EijR(4,4)+EijR(11,3)+EijR(11
     -   ,4)-4*EijR(13,3)-EijR(13,4)+EijR(14,3)-EijR(14,4)+5*EijR(16,
     -   3)-2*(EijR(4,3)+EijR(10,2)+EijR(19,3))-3*(EijR(16,4)-EijR(22
     -   ,4)+EijR(30,4)-EijR(33,4))))
       FI(174) = Is12*(Dij1345I(2,2)+Dij1345I(2,3)+Dij1345I(3,2)-Dij13
     -   45I(3,3)-2*Dij1345I(6,2)-3*(Dij1345I(8,3)-Dij1345I(9,3))+s12
     -   *(EijI(3,2)+EijI(3,3)+EijI(4,2)+EijI(4,4)+EijI(11,3)+EijI(11
     -   ,4)-4*EijI(13,3)-EijI(13,4)+EijI(14,3)-EijI(14,4)+5*EijI(16,
     -   3)-2*(EijI(4,3)+EijI(10,2)+EijI(19,3))-3*(EijI(16,4)-EijI(22
     -   ,4)+EijI(30,4)-EijI(33,4))))
       F(174)=DCMPLX(FR(174),FI(174))
       FR(175) = EijR(3,2)+EijR(3,3)+EijR(4,2)+EijR(4,4)+EijR(12,3)+Ei
     -   jR(12,4)-4*EijR(13,3)-EijR(13,4)+EijR(15,3)-EijR(15,4)+5*Eij
     -   R(16,3)-2*(EijR(4,3)+EijR(10,2)+EijR(20,3))-3*(EijR(16,4)-Ei
     -   jR(22,4)+EijR(31,4)-EijR(34,4))
       FI(175) = EijI(3,2)+EijI(3,3)+EijI(4,2)+EijI(4,4)+EijI(12,3)+Ei
     -   jI(12,4)-4*EijI(13,3)-EijI(13,4)+EijI(15,3)-EijI(15,4)+5*Eij
     -   I(16,3)-2*(EijI(4,3)+EijI(10,2)+EijI(20,3))-3*(EijI(16,4)-Ei
     -   jI(22,4)+EijI(31,4)-EijI(34,4))
       F(175)=DCMPLX(FR(175),FI(175))
       FR(176) = EijR(3,3)+EijR(3,4)-EijR(4,3)+EijR(4,4)-3*(EijR(13,3)
     -   -EijR(16,3))-4*(EijR(13,4)+EijR(16,4))+6*EijR(22,4)
       FI(176) = EijI(3,3)+EijI(3,4)-EijI(4,3)+EijI(4,4)-3*(EijI(13,3)
     -   -EijI(16,3))-4*(EijI(13,4)+EijI(16,4))+6*EijI(22,4)
       F(176)=DCMPLX(FR(176),FI(176))
       FR(177) = -EijR(4,3)+EijR(4,4)-EijR(13,3)-EijR(13,4)+2*EijR(16,
     -   3)-3*(EijR(16,4)-EijR(22,4))
       FI(177) = -EijI(4,3)+EijI(4,4)-EijI(13,3)-EijI(13,4)+2*EijI(16,
     -   3)-3*(EijI(16,4)-EijI(22,4))
       F(177)=DCMPLX(FR(177),FI(177))
       FR(178) = Is12*(-Dij1345R(3,3)-Dij1345R(8,3)+2*Dij1345R(9,3)-s1
     -   2*(EijR(4,3)-EijR(4,4)+EijR(13,3)+EijR(14,4)-EijR(22,4)+EijR
     -   (30,4)-2*(EijR(16,3)-EijR(16,4)+EijR(33,4))))
       FI(178) = Is12*(-Dij1345I(3,3)-Dij1345I(8,3)+2*Dij1345I(9,3)-s1
     -   2*(EijI(4,3)-EijI(4,4)+EijI(13,3)+EijI(14,4)-EijI(22,4)+EijI
     -   (30,4)-2*(EijI(16,3)-EijI(16,4)+EijI(33,4))))
       F(178)=DCMPLX(FR(178),FI(178))
       FR(179) = EijR(4,2)-EijR(4,3)+EijR(4,4)-EijR(10,2)-EijR(13,3)-E
     -   ijR(15,4)+EijR(22,4)-EijR(31,4)+2*(EijR(16,3)-EijR(16,4)+Eij
     -   R(34,4))
       FI(179) = EijI(4,2)-EijI(4,3)+EijI(4,4)-EijI(10,2)-EijI(13,3)-E
     -   ijI(15,4)+EijI(22,4)-EijI(31,4)+2*(EijI(16,3)-EijI(16,4)+Eij
     -   I(34,4))
       F(179)=DCMPLX(FR(179),FI(179))
       FR(180) = EijR(4,2)+2*EijR(4,3)-EijR(4,4)-EijR(7,2)+EijR(8,2)-E
     -   ijR(9,2)-EijR(14,3)+EijR(14,4)-EijR(15,3)+EijR(15,4)-EijR(16
     -   ,3)+EijR(16,4)+EijR(17,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+E
     -   ijR(35,4)
       FI(180) = EijI(4,2)+2*EijI(4,3)-EijI(4,4)-EijI(7,2)+EijI(8,2)-E
     -   ijI(9,2)-EijI(14,3)+EijI(14,4)-EijI(15,3)+EijI(15,4)-EijI(16
     -   ,3)+EijI(16,4)+EijI(17,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+E
     -   ijI(35,4)
       F(180)=DCMPLX(FR(180),FI(180))
       FR(181) = EijR(3,2)+EijR(4,2)+EijR(4,3)-EijR(4,4)+EijR(11,3)+Ei
     -   jR(14,4)-EijR(16,3)-EijR(19,3)-EijR(22,4)+EijR(30,4)-2*(EijR
     -   (10,2)-EijR(16,4)+EijR(33,4))
       FI(181) = EijI(3,2)+EijI(4,2)+EijI(4,3)-EijI(4,4)+EijI(11,3)+Ei
     -   jI(14,4)-EijI(16,3)-EijI(19,3)-EijI(22,4)+EijI(30,4)-2*(EijI
     -   (10,2)-EijI(16,4)+EijI(33,4))
       F(181)=DCMPLX(FR(181),FI(181))
       FR(182) = EijR(4,2)+EijR(4,3)-EijR(4,4)-EijR(10,2)+EijR(14,4)+E
     -   ijR(16,4)-EijR(19,3)-EijR(33,4)
       FI(182) = EijI(4,2)+EijI(4,3)-EijI(4,4)-EijI(10,2)+EijI(14,4)+E
     -   ijI(16,4)-EijI(19,3)-EijI(33,4)
       F(182)=DCMPLX(FR(182),FI(182))
       FR(183) = EijR(3,1)+EijR(3,2)-EijR(4,1)+EijR(4,2)-EijR(4,4)+Eij
     -   R(8,2)-EijR(9,2)+EijR(12,3)+EijR(15,3)+EijR(15,4)-EijR(22,4)
     -   +EijR(31,4)-2*(EijR(10,2)-EijR(16,4)+EijR(20,3)+EijR(34,4))
       FI(183) = EijI(3,1)+EijI(3,2)-EijI(4,1)+EijI(4,2)-EijI(4,4)+Eij
     -   I(8,2)-EijI(9,2)+EijI(12,3)+EijI(15,3)+EijI(15,4)-EijI(22,4)
     -   +EijI(31,4)-2*(EijI(10,2)-EijI(16,4)+EijI(20,3)+EijI(34,4))
       F(183)=DCMPLX(FR(183),FI(183))
       FR(184) = 2*EijR(4,2)-EijR(4,4)-EijR(9,2)-EijR(10,2)+EijR(15,3)
     -   +EijR(15,4)+EijR(16,4)-EijR(20,3)-EijR(34,4)
       FI(184) = 2*EijI(4,2)-EijI(4,4)-EijI(9,2)-EijI(10,2)+EijI(15,3)
     -   +EijI(15,4)+EijI(16,4)-EijI(20,3)-EijI(34,4)
       F(184)=DCMPLX(FR(184),FI(184))
       FR(185) = Is12*(Dij1345R(2,1)+Dij1345R(2,2)-Dij1345R(3,1)+Dij13
     -   45R(3,3)-Dij1345R(6,2)+Dij1345R(8,3)-2*Dij1345R(9,3)+s12*(Ei
     -   jR(3,1)+EijR(3,2)-EijR(4,1)+EijR(4,2)+EijR(4,3)-EijR(4,4)+Ei
     -   jR(6,2)-EijR(7,2)+EijR(11,3)+EijR(14,4)-EijR(16,3)-EijR(19,3
     -   )-EijR(22,4)+EijR(30,4)-2*(EijR(10,2)-EijR(16,4)+EijR(33,4))
     -   ))
       FI(185) = Is12*(Dij1345I(2,1)+Dij1345I(2,2)-Dij1345I(3,1)+Dij13
     -   45I(3,3)-Dij1345I(6,2)+Dij1345I(8,3)-2*Dij1345I(9,3)+s12*(Ei
     -   jI(3,1)+EijI(3,2)-EijI(4,1)+EijI(4,2)+EijI(4,3)-EijI(4,4)+Ei
     -   jI(6,2)-EijI(7,2)+EijI(11,3)+EijI(14,4)-EijI(16,3)-EijI(19,3
     -   )-EijI(22,4)+EijI(30,4)-2*(EijI(10,2)-EijI(16,4)+EijI(33,4))
     -   ))
       F(185)=DCMPLX(FR(185),FI(185))
       FR(186) = EijR(4,2)+EijR(4,3)-EijR(4,4)+EijR(8,2)-EijR(9,2)-Eij
     -   R(10,2)+EijR(12,3)+EijR(15,4)-EijR(16,3)-EijR(20,3)-EijR(22,
     -   4)+EijR(31,4)+2*(EijR(16,4)-EijR(34,4))
       FI(186) = EijI(4,2)+EijI(4,3)-EijI(4,4)+EijI(8,2)-EijI(9,2)-Eij
     -   I(10,2)+EijI(12,3)+EijI(15,4)-EijI(16,3)-EijI(20,3)-EijI(22,
     -   4)+EijI(31,4)+2*(EijI(16,4)-EijI(34,4))
       F(186)=DCMPLX(FR(186),FI(186))
       FR(187) = EijR(3,2)+EijR(3,3)+EijR(4,2)-EijR(4,4)-2*(EijR(10,2)
     -   +EijR(13,3))+EijR(13,4)+EijR(16,3)+3*(EijR(16,4)-EijR(22,4))
       FI(187) = EijI(3,2)+EijI(3,3)+EijI(4,2)-EijI(4,4)-2*(EijI(10,2)
     -   +EijI(13,3))+EijI(13,4)+EijI(16,3)+3*(EijI(16,4)-EijI(22,4))
       F(187)=DCMPLX(FR(187),FI(187))
       FR(188) = EijR(4,2)-EijR(4,4)-EijR(10,2)-EijR(13,3)+EijR(16,3)+
     -   2*EijR(16,4)-EijR(22,4)
       FI(188) = EijI(4,2)-EijI(4,4)-EijI(10,2)-EijI(13,3)+EijI(16,3)+
     -   2*EijI(16,4)-EijI(22,4)
       F(188)=DCMPLX(FR(188),FI(188))
       FR(189) = Is12*(Dij1345R(3,2)+Dij1345R(3,3)-Dij1345R(6,2)-Dij13
     -   45R(9,3)+s12*(EijR(4,2)-EijR(4,4)-EijR(10,2)+EijR(14,3)+EijR
     -   (14,4)+EijR(16,4)-EijR(19,3)-EijR(33,4)))
       FI(189) = Is12*(Dij1345I(3,2)+Dij1345I(3,3)-Dij1345I(6,2)-Dij13
     -   45I(9,3)+s12*(EijI(4,2)-EijI(4,4)-EijI(10,2)+EijI(14,3)+EijI
     -   (14,4)+EijI(16,4)-EijI(19,3)-EijI(33,4)))
       F(189)=DCMPLX(FR(189),FI(189))
       FR(190) = EijR(4,4)-EijR(15,3)-EijR(15,4)-EijR(16,4)+EijR(20,3)
     -   +EijR(34,4)
       FI(190) = EijI(4,4)-EijI(15,3)-EijI(15,4)-EijI(16,4)+EijI(20,3)
     -   +EijI(34,4)
       F(190)=DCMPLX(FR(190),FI(190))
       FR(191) = -EijR(4,3)-EijR(4,4)-EijR(13,3)+2*(EijR(16,3)+EijR(16
     -   ,4))-EijR(22,4)
       FI(191) = -EijI(4,3)-EijI(4,4)-EijI(13,3)+2*(EijI(16,3)+EijI(16
     -   ,4))-EijI(22,4)
       F(191)=DCMPLX(FR(191),FI(191))
       FR(192) = -EijR(4,3)-EijR(4,4)+EijR(16,3)+EijR(16,4)
       FI(192) = -EijI(4,3)-EijI(4,4)+EijI(16,3)+EijI(16,4)
       F(192)=DCMPLX(FR(192),FI(192))
       FR(193) = Is45*(-Dij1234R(2,1)-Dij1234R(2,2)+Dij1234R(3,1)+Dij1
     -   234R(3,2)-Dij1234R(8,3)+Dij1234R(9,3)+s45*(EijR(4,2)+EijR(4,
     -   4)+EijR(8,2)-EijR(9,2)+EijR(9,3)-EijR(10,2)-EijR(10,3)+EijR(
     -   12,3)-EijR(13,3)+3*(EijR(15,3)+EijR(16,3))-4*EijR(20,3)+EijR
     -   (20,4)+EijR(21,4)+EijR(22,4)-2*(EijR(4,3)+EijR(15,4)+EijR(16
     -   ,4)+EijR(28,4)+EijR(31,4)-2*EijR(34,4))))
       FI(193) = Is45*(-Dij1234I(2,1)-Dij1234I(2,2)+Dij1234I(3,1)+Dij1
     -   234I(3,2)-Dij1234I(8,3)+Dij1234I(9,3)+s45*(EijI(4,2)+EijI(4,
     -   4)+EijI(8,2)-EijI(9,2)+EijI(9,3)-EijI(10,2)-EijI(10,3)+EijI(
     -   12,3)-EijI(13,3)+3*(EijI(15,3)+EijI(16,3))-4*EijI(20,3)+EijI
     -   (20,4)+EijI(21,4)+EijI(22,4)-2*(EijI(4,3)+EijI(15,4)+EijI(16
     -   ,4)+EijI(28,4)+EijI(31,4)-2*EijI(34,4))))
       F(193)=DCMPLX(FR(193),FI(193))
       FR(194) = Is45*(Dij1234R(3,2)-Dij1234R(6,2)-Dij1234R(8,3)+Dij12
     -   34R(9,3)+s45*(-EijR(4,3)+EijR(4,4)+EijR(12,3)-EijR(13,3)+Eij
     -   R(15,3)+EijR(20,4)+EijR(21,4)+EijR(22,4)-2*(EijR(15,4)-EijR(
     -   16,3)+EijR(16,4)+EijR(20,3)+EijR(28,4)+EijR(31,4)-2*EijR(34,
     -   4))))
       FI(194) = Is45*(Dij1234I(3,2)-Dij1234I(6,2)-Dij1234I(8,3)+Dij12
     -   34I(9,3)+s45*(-EijI(4,3)+EijI(4,4)+EijI(12,3)-EijI(13,3)+Eij
     -   I(15,3)+EijI(20,4)+EijI(21,4)+EijI(22,4)-2*(EijI(15,4)-EijI(
     -   16,3)+EijI(16,4)+EijI(20,3)+EijI(28,4)+EijI(31,4)-2*EijI(34,
     -   4))))
       F(194)=DCMPLX(FR(194),FI(194))
       FR(195) = Is45*(Dij1234R(3,2)+Dij1234R(3,3)-Dij1234R(6,2)-Dij12
     -   34R(9,3)+s45*(-EijR(4,3)+EijR(4,4)+EijR(12,3)+EijR(12,4)-Eij
     -   R(13,3)-EijR(13,4)+EijR(15,3)-EijR(15,4)+2*(EijR(16,3)-EijR(
     -   20,3))-3*(EijR(16,4)-EijR(22,4)+EijR(31,4)-EijR(34,4))))
       FI(195) = Is45*(Dij1234I(3,2)+Dij1234I(3,3)-Dij1234I(6,2)-Dij12
     -   34I(9,3)+s45*(-EijI(4,3)+EijI(4,4)+EijI(12,3)+EijI(12,4)-Eij
     -   I(13,3)-EijI(13,4)+EijI(15,3)-EijI(15,4)+2*(EijI(16,3)-EijI(
     -   20,3))-3*(EijI(16,4)-EijI(22,4)+EijI(31,4)-EijI(34,4))))
       F(195)=DCMPLX(FR(195),FI(195))
       P(443) = 2*p3sq-p4sq+s15-4*s23+s34+3*P(377)
       P(444) = 3*p2sq+p3sq-s23-2*P(28)
       P(445) = p2sq-s23-2*P(235)
       P(446) = s23-s34
       P(447) = p3sq+p4sq-s15+s45+3*P(18)-2*P(446)
       P(448) = s15-s34-2*P(6)
       P(449) = -p3sq+2*s12+s15+s23-s34
       P(450) = -3*p2sq**2+p2sq*P(448)+s12*P(449)
       P(451) = p4sq+5*s23-s34+3*P(52)
       P(452) = -5*p3sq-s12+s45
       P(453) = -3*p3sq**2+P(52)*P(168)+p3sq*P(451)+p2sq*P(452)
       P(454) = -3*s12-s15+s34+2*P(6)
       P(455) = p3sq-3*s12-s15-2*P(446)
       P(456) = s45**2+3*p2sq*P(174)+p4sq*P(454)+s45*P(455)
       P(457) = -3*p3sq+4*s12+5*s23-s34-2*s45
       P(458) = s15+s23
       P(459) = 3*p3sq-7*s12+s34-s45-2*P(458)
       P(460) = 3*p2sq**2+s12*P(457)+p2sq*P(459)
       P(461) = 3*s12+s15
       P(462) = -3*p3sq+p4sq+7*s12+5*s23-s34-6*s45
       P(463) = s12**2-s23*s45
       P(464) = s23**2+s45**2
       P(465) = -3*p2sq**2-p4sq*s12+s12*s15-5*s12*s23+s12*s34+6*s12*s4
     -   5-s15*s45+p3sq*P(461)+p2sq*P(462)-4*P(463)-2*P(464)
       P(466) = p2sq-s12+s34
       P(467) = s23*P(466)
       P(468) = 3*s34+4*P(66)
       P(469) = p2sq*p4sq-p4sq*s12+s23**2+s45**2
       P(470) = -(p3sq*s15)+s15*s23+s15*s34-3*P(467)+s45*P(468)+2*P(46
     -   9)
       P(471) = p2sq**2+p3sq**2
       P(472) = 7*p3sq-p4sq+s12-5*s23+s34
       P(473) = s12*s23-s23**2-s12*s34
       P(474) = 2*s12+s15+5*s23-s34+s45
       P(475) = -(p4sq*s12)-4*s12**2-s12*s15+s45*P(189)+3*P(471)+p2sq*
     -   P(472)-2*P(473)-p3sq*P(474)
       P(476) = p3sq**2+s12**2
       P(477) = s12+s15+s23-s34
       P(478) = 7*s12+2*s15+5*s23-s34+s45
       P(479) = 6*p3sq+p4sq-3*P(290)
       P(480) = s12*s15+5*s12*s23+s15*s23+2*s23**2+s12*s34+s15*s34-s34
     -   **2+3*P(476)+s45*P(477)-p3sq*P(478)+p2sq*P(479)
       P(481) = p3sq**2+s12**2+s23*s34
       P(482) = s15+s23+2*P(136)
       P(483) = -2*s23+s34+s45
       P(484) = p4sq-7*s12-s15-5*s23+s34+4*s45
       P(485) = 6*p3sq-3*P(136)+2*P(174)
       P(486) = s12*s15+2*s12*s23+s12*s34-s34**2-s45**2+3*P(481)-s45*P
     -   (482)+p4sq*P(483)+p3sq*P(484)+p2sq*P(485)
       P(487) = 14*p2sq+9*p3sq-p4sq-11*s12-s15-12*s23+3*s34+4*s45
       P(488) = p2sq-p3sq+2*s12+s23
       P(489) = p3sq-p4sq-s34-2*s45
       FR(196) = Is12s45*((3*B013R+B014R)*P(5)+2*(P(9)*(-(p2sq*Cij123R
     -   (2,1))+p4sq*Cij145R(2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij145
     -   R(4,2))+C0145R*P(15)+(C0123R+Cij123R(1,1))*P(18)+Cij134R(2,1
     -   )*P(25))+C0134R*P(280)))+2*(-C0234R+EE0R*s15*P(9)+s15*EijR(1
     -   ,2)*P(18)+4*(Dij2345R(7,2)-EijR(23,3)*P(52))+s15*EijR(1,1)*P
     -   (159)+Is45*(P(9)*(-2*p2sq*Dij1234R(2,2)+Dij1234R(4,2)*P(18)-
     -   Dij1234R(3,1)*P(45))-(Dij1234R(3,2)*P(1)+Dij1234R(5,2)*P(18)
     -   )*P(78)+(D01234R+Dij1234R(1,1))*P(151)+Dij1234R(2,1)*P(156)+
     -   Dij1234R(6,2)*P(158)+2*(Cij134R(1,1)*P(5)-Dij1234R(7,2)*P(19
     -   9)))+Dij1345R(1,2)*P(261)+D01345R*P(281)+Is12*((Cij145R(1,1)
     -   +2*p4sq*Dij1345R(3,2))*P(5)+Dij1345R(1,1)*P(283)+Dij1345R(2,
     -   1)*P(285)-Dij1345R(3,1)*P(286)-Dij1345R(2,2)*P(1)*P(287)+Dij
     -   1345R(4,2)*P(289)-Dij1345R(5,2)*P(291)+Dij1345R(6,2)*P(295))
     -   +EijR(2,1)*P(297)+EijR(3,1)*P(299)+EijR(4,1)*P(305)+D02345R*
     -   P(443)-Dij2345R(1,1)*P(444)+Dij2345R(2,1)*P(445)+Dij2345R(3,
     -   1)*P(447)+EijR(2,2)*P(450)+EijR(3,2)*P(453)-EijR(4,2)*P(456)
     -   +EijR(5,2)*P(460)+EijR(6,2)*P(465)+EijR(7,2)*P(470)+EijR(8,2
     -   )*P(475)-EijR(9,2)*P(480)+EijR(10,2)*P(486)+2*(Dij1345R(7,2)
     -   *P(109)-EijR(21,3)*P(434)-EijR(11,2)*P(487)+EijR(22,3)*P(488
     -   )+EijR(24,3)*P(489)))
       FI(196) = Is12s45*((3*B013I+B014I)*P(5)+2*(P(9)*(-(p2sq*Cij123I
     -   (2,1))+p4sq*Cij145I(2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij145
     -   I(4,2))+C0145I*P(15)+(C0123I+Cij123I(1,1))*P(18)+Cij134I(2,1
     -   )*P(25))+C0134I*P(280)))+2*(-C0234I+EE0I*s15*P(9)+s15*EijI(1
     -   ,2)*P(18)+4*(Dij2345I(7,2)-EijI(23,3)*P(52))+s15*EijI(1,1)*P
     -   (159)+Is45*(P(9)*(-2*p2sq*Dij1234I(2,2)+Dij1234I(4,2)*P(18)-
     -   Dij1234I(3,1)*P(45))-(Dij1234I(3,2)*P(1)+Dij1234I(5,2)*P(18)
     -   )*P(78)+(D01234I+Dij1234I(1,1))*P(151)+Dij1234I(2,1)*P(156)+
     -   Dij1234I(6,2)*P(158)+2*(Cij134I(1,1)*P(5)-Dij1234I(7,2)*P(19
     -   9)))+Dij1345I(1,2)*P(261)+D01345I*P(281)+Is12*((Cij145I(1,1)
     -   +2*p4sq*Dij1345I(3,2))*P(5)+Dij1345I(1,1)*P(283)+Dij1345I(2,
     -   1)*P(285)-Dij1345I(3,1)*P(286)-Dij1345I(2,2)*P(1)*P(287)+Dij
     -   1345I(4,2)*P(289)-Dij1345I(5,2)*P(291)+Dij1345I(6,2)*P(295))
     -   +EijI(2,1)*P(297)+EijI(3,1)*P(299)+EijI(4,1)*P(305)+D02345I*
     -   P(443)-Dij2345I(1,1)*P(444)+Dij2345I(2,1)*P(445)+Dij2345I(3,
     -   1)*P(447)+EijI(2,2)*P(450)+EijI(3,2)*P(453)-EijI(4,2)*P(456)
     -   +EijI(5,2)*P(460)+EijI(6,2)*P(465)+EijI(7,2)*P(470)+EijI(8,2
     -   )*P(475)-EijI(9,2)*P(480)+EijI(10,2)*P(486)+2*(Dij1345I(7,2)
     -   *P(109)-EijI(21,3)*P(434)-EijI(11,2)*P(487)+EijI(22,3)*P(488
     -   )+EijI(24,3)*P(489)))
       F(196)=DCMPLX(FR(196),FI(196))
       P(490) = -p2sq+p3sq+3*p4sq+s23+2*P(219)
       P(491) = 2*p2sq+p3sq-s23
       P(492) = p2sq-p3sq
       P(493) = s23+3*P(235)-2*P(492)
       P(494) = 2*p3sq-s23-s34-s45
       P(495) = 2*p3sq-s12+s45
       P(496) = p2sq+p3sq-s12-s34-s45
       P(497) = p2sq-p3sq+s12+s34+s45
       FR(197) = -4*(-Dij2345R(3,2)+Dij2345R(5,2)+p3sq*EijR(3,1)+3*p3s
     -   q*EijR(6,2)+4*(EijR(11,2)-EijR(21,3))+2*(Dij2345R(1,1)-Dij23
     -   45R(3,1)-EijR(22,3))+6*EijR(24,3)-EijR(2,1)*P(9)+(EijR(5,3)-
     -   EijR(7,3))*P(18)-(EijR(8,2)-EijR(16,3)-EijR(17,3)+EijR(19,3)
     -   +EijR(20,3))*P(25)+EijR(4,1)*P(66)+(EijR(2,2)-EijR(8,3)+EijR
     -   (10,3))*P(68)+EijR(18,3)*P(96)+Is12*(Cij145R(1,1)+Cij145R(1,
     -   2)-Cij145R(2,1)+Cij145R(2,2)+p3sq*(Dij1345R(2,1)-Dij1345R(3,
     -   1)+Dij1345R(4,2)-Dij1345R(6,2))-2*(Cij145R(3,2)-Dij1345R(7,2
     -   )-Dij1345R(11,3)+Dij1345R(13,3))+(Dij1345R(3,2)-Dij1345R(5,2
     -   ))*P(203))+EijR(4,3)*P(261)+EijR(4,2)*P(490)-EijR(5,2)*P(491
     -   )-EijR(7,2)*P(493)+EijR(9,2)*P(494)-EijR(10,2)*P(495)+EijR(1
     -   4,3)*P(496)-EijR(15,3)*P(497))
       FI(197) = -4*(-Dij2345I(3,2)+Dij2345I(5,2)+p3sq*EijI(3,1)+3*p3s
     -   q*EijI(6,2)+4*(EijI(11,2)-EijI(21,3))+2*(Dij2345I(1,1)-Dij23
     -   45I(3,1)-EijI(22,3))+6*EijI(24,3)-EijI(2,1)*P(9)+(EijI(5,3)-
     -   EijI(7,3))*P(18)-(EijI(8,2)-EijI(16,3)-EijI(17,3)+EijI(19,3)
     -   +EijI(20,3))*P(25)+EijI(4,1)*P(66)+(EijI(2,2)-EijI(8,3)+EijI
     -   (10,3))*P(68)+EijI(18,3)*P(96)+Is12*(Cij145I(1,1)+Cij145I(1,
     -   2)-Cij145I(2,1)+Cij145I(2,2)+p3sq*(Dij1345I(2,1)-Dij1345I(3,
     -   1)+Dij1345I(4,2)-Dij1345I(6,2))-2*(Cij145I(3,2)-Dij1345I(7,2
     -   )-Dij1345I(11,3)+Dij1345I(13,3))+(Dij1345I(3,2)-Dij1345I(5,2
     -   ))*P(203))+EijI(4,3)*P(261)+EijI(4,2)*P(490)-EijI(5,2)*P(491
     -   )-EijI(7,2)*P(493)+EijI(9,2)*P(494)-EijI(10,2)*P(495)+EijI(1
     -   4,3)*P(496)-EijI(15,3)*P(497))
       F(197)=DCMPLX(FR(197),FI(197))
       P(498) = p4sq+s12-s34
       P(499) = -p2sq+p3sq+s23-3*s45+2*P(498)
       P(500) = -p2sq+p3sq+s12+s23-3*s45+2*P(128)
       P(501) = -p3sq+s34+s45+2*P(68)
       P(502) = p2sq+s12-2*P(95)
       FR(198) = -4*(-Dij2345R(1,2)-Dij2345R(3,2)+2*Dij2345R(5,2)-s12*
     -   EijR(2,2)+p3sq*EijR(8,2)-p3sq*EijR(10,2)-6*(EijR(22,3)-EijR(
     -   24,3))+(EijR(8,3)+EijR(14,3)-2*EijR(18,3))*P(18)+(EijR(9,3)+
     -   EijR(16,3)-2*EijR(20,3))*P(25)-(EijR(2,1)+EijR(5,2))*P(45)-(
     -   EijR(4,1)+EijR(7,2))*P(46)-EijR(2,3)*P(68)+Is12*(Cij145R(1,1
     -   )+Cij145R(1,2)-Cij145R(2,1)+Cij145R(2,2)+p3sq*(Dij1345R(2,1)
     -   -Dij1345R(3,1)+Dij1345R(4,2)-Dij1345R(6,2))-2*(Cij145R(3,2)-
     -   Dij1345R(7,2)-Dij1345R(11,3)+Dij1345R(13,3))+(Dij1345R(3,2)-
     -   Dij1345R(5,2))*P(203))+EijR(4,3)*P(261)+EijR(4,2)*P(499)-Eij
     -   R(9,2)*P(500)+EijR(10,3)*P(501)-EijR(15,3)*P(502))
       FI(198) = -4*(-Dij2345I(1,2)-Dij2345I(3,2)+2*Dij2345I(5,2)-s12*
     -   EijI(2,2)+p3sq*EijI(8,2)-p3sq*EijI(10,2)-6*(EijI(22,3)-EijI(
     -   24,3))+(EijI(8,3)+EijI(14,3)-2*EijI(18,3))*P(18)+(EijI(9,3)+
     -   EijI(16,3)-2*EijI(20,3))*P(25)-(EijI(2,1)+EijI(5,2))*P(45)-(
     -   EijI(4,1)+EijI(7,2))*P(46)-EijI(2,3)*P(68)+Is12*(Cij145I(1,1
     -   )+Cij145I(1,2)-Cij145I(2,1)+Cij145I(2,2)+p3sq*(Dij1345I(2,1)
     -   -Dij1345I(3,1)+Dij1345I(4,2)-Dij1345I(6,2))-2*(Cij145I(3,2)-
     -   Dij1345I(7,2)-Dij1345I(11,3)+Dij1345I(13,3))+(Dij1345I(3,2)-
     -   Dij1345I(5,2))*P(203))+EijI(4,3)*P(261)+EijI(4,2)*P(499)-Eij
     -   I(9,2)*P(500)+EijI(10,3)*P(501)-EijI(15,3)*P(502))
       F(198)=DCMPLX(FR(198),FI(198))
       P(503) = -p2sq+p3sq+s23-5*s45+2*P(230)
       P(504) = p2sq-s12-3*P(36)
       P(505) = 4*p3sq+s12+2*P(226)
       P(506) = s23+s34+s45
       P(507) = -4*p3sq+3*s12+2*P(506)
       P(508) = p2sq+2*p3sq-4*s12-s23+5*s45
       P(509) = s12-s34+2*P(30)
       P(510) = p2sq+s34-2*P(30)
       FR(199) = 2*Is12s45*(3*B013R+B014R-2*(C0134R*p3sq-p2sq*Cij123R(
     -   2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2))+(C0145R+Cij1
     -   45R(2,1))*P(15)+(C0123R+Cij123R(1,1))*P(18)+Cij134R(2,1)*P(2
     -   5)))+4*(-2*(D02345R-Dij2345R(1,1)+Dij2345R(3,1))+Dij2345R(3,
     -   2)+Dij2345R(4,2)-Dij2345R(5,2)-Dij2345R(6,2)-s15*(EE0R+EijR(
     -   1,1))+3*p3sq*EijR(3,2)+18*EijR(11,2)+4*(EijR(23,3)-EijR(24,3
     -   ))+(-2*EijR(5,2)-EijR(14,3)-EijR(17,3)+EijR(18,3)+EijR(19,3)
     -   )*P(18)+(-EijR(12,3)+EijR(13,3))*P(25)+(D01345R+Dij1345R(1,1
     -   ))*P(39)+(D01234R+Dij1234R(1,1))*P(40)+(EijR(3,1)+3*EijR(6,2
     -   ))*P(45)+EijR(4,1)*P(46)+(2*EijR(2,2)+EijR(9,3)-EijR(10,3))*
     -   P(68)+Is45*(2*(Cij134R(1,1)+p2sq*Dij1234R(2,2))+6*Dij1234R(7
     -   ,2)+Dij1234R(3,2)*P(1)-(Dij1234R(2,1)+Dij1234R(4,2)-Dij1234R
     -   (5,2))*P(18)+Dij1234R(3,1)*P(66)-Dij1234R(6,2)*P(157))-EijR(
     -   4,3)*P(261)-Is12*(Cij145R(1,2)+Cij145R(2,2)-2*(Cij145R(3,2)+
     -   Dij1345R(7,2)-Dij1345R(12,3)+Dij1345R(13,3))+(Dij1345R(2,1)-
     -   Dij1345R(3,1))*P(30)-Dij1345R(2,2)*P(128)+Dij1345R(3,2)*P(26
     -   0)+(-Dij1345R(4,2)+Dij1345R(5,2))*P(261)-Dij1345R(6,2)*P(264
     -   ))+EijR(15,3)*P(497)-EijR(4,2)*P(503)-EijR(7,2)*P(504)-EijR(
     -   8,2)*P(505)-EijR(9,2)*P(507)-EijR(10,2)*P(508)-EijR(16,3)*P(
     -   509)-EijR(20,3)*P(510))
       FI(199) = 2*Is12s45*(3*B013I+B014I-2*(C0134I*p3sq-p2sq*Cij123I(
     -   2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2))+(C0145I+Cij1
     -   45I(2,1))*P(15)+(C0123I+Cij123I(1,1))*P(18)+Cij134I(2,1)*P(2
     -   5)))+4*(-2*(D02345I-Dij2345I(1,1)+Dij2345I(3,1))+Dij2345I(3,
     -   2)+Dij2345I(4,2)-Dij2345I(5,2)-Dij2345I(6,2)-s15*(EE0I+EijI(
     -   1,1))+3*p3sq*EijI(3,2)+18*EijI(11,2)+4*(EijI(23,3)-EijI(24,3
     -   ))+(-2*EijI(5,2)-EijI(14,3)-EijI(17,3)+EijI(18,3)+EijI(19,3)
     -   )*P(18)+(-EijI(12,3)+EijI(13,3))*P(25)+(D01345I+Dij1345I(1,1
     -   ))*P(39)+(D01234I+Dij1234I(1,1))*P(40)+(EijI(3,1)+3*EijI(6,2
     -   ))*P(45)+EijI(4,1)*P(46)+(2*EijI(2,2)+EijI(9,3)-EijI(10,3))*
     -   P(68)+Is45*(2*(Cij134I(1,1)+p2sq*Dij1234I(2,2))+6*Dij1234I(7
     -   ,2)+Dij1234I(3,2)*P(1)-(Dij1234I(2,1)+Dij1234I(4,2)-Dij1234I
     -   (5,2))*P(18)+Dij1234I(3,1)*P(66)-Dij1234I(6,2)*P(157))-EijI(
     -   4,3)*P(261)-Is12*(Cij145I(1,2)+Cij145I(2,2)-2*(Cij145I(3,2)+
     -   Dij1345I(7,2)-Dij1345I(12,3)+Dij1345I(13,3))+(Dij1345I(2,1)-
     -   Dij1345I(3,1))*P(30)-Dij1345I(2,2)*P(128)+Dij1345I(3,2)*P(26
     -   0)+(-Dij1345I(4,2)+Dij1345I(5,2))*P(261)-Dij1345I(6,2)*P(264
     -   ))+EijI(15,3)*P(497)-EijI(4,2)*P(503)-EijI(7,2)*P(504)-EijI(
     -   8,2)*P(505)-EijI(9,2)*P(507)-EijI(10,2)*P(508)-EijI(16,3)*P(
     -   509)-EijI(20,3)*P(510))
       F(199)=DCMPLX(FR(199),FI(199))
       P(511) = p2sq-2*s12-s23+s34+4*s45
       P(512) = 2*P(18)-3*P(36)
       P(513) = p2sq+3*p3sq-s34-s45-2*P(177)
       FR(200) = -4*(D01345R-Dij2345R(1,1)+Dij2345R(3,1)-Dij2345R(3,2)
     -   +Dij2345R(5,2)-s12*EijR(2,1)+p3sq*EijR(3,1)-4*(EijR(11,2)-Ei
     -   jR(24,3))+(EijR(5,2)+EijR(14,3)-EijR(18,3))*P(18)+(EijR(8,2)
     -   +EijR(16,3)-EijR(20,3))*P(25)-EijR(4,1)*P(42)+(-EijR(2,2)+Ei
     -   jR(10,3))*P(68)+EijR(4,3)*P(261)+Is12*(C0134R-C0145R-Cij145R
     -   (1,1)+Cij145R(2,2)-Cij145R(3,2)-4*Dij1345R(7,2)-2*Dij1345R(1
     -   3,3)+(Dij1345R(2,1)+Dij1345R(6,2))*P(128)+Dij1345R(1,1)*P(13
     -   6)+Dij1345R(3,2)*P(260)+Dij1345R(5,2)*P(261)-Dij1345R(3,1)*P
     -   (277))+EijR(10,2)*P(495)-EijR(15,3)*P(497)-EijR(4,2)*P(511)+
     -   EijR(7,2)*P(512)-EijR(9,2)*P(513))
       FI(200) = -4*(D01345I-Dij2345I(1,1)+Dij2345I(3,1)-Dij2345I(3,2)
     -   +Dij2345I(5,2)-s12*EijI(2,1)+p3sq*EijI(3,1)-4*(EijI(11,2)-Ei
     -   jI(24,3))+(EijI(5,2)+EijI(14,3)-EijI(18,3))*P(18)+(EijI(8,2)
     -   +EijI(16,3)-EijI(20,3))*P(25)-EijI(4,1)*P(42)+(-EijI(2,2)+Ei
     -   jI(10,3))*P(68)+EijI(4,3)*P(261)+Is12*(C0134I-C0145I-Cij145I
     -   (1,1)+Cij145I(2,2)-Cij145I(3,2)-4*Dij1345I(7,2)-2*Dij1345I(1
     -   3,3)+(Dij1345I(2,1)+Dij1345I(6,2))*P(128)+Dij1345I(1,1)*P(13
     -   6)+Dij1345I(3,2)*P(260)+Dij1345I(5,2)*P(261)-Dij1345I(3,1)*P
     -   (277))+EijI(10,2)*P(495)-EijI(15,3)*P(497)-EijI(4,2)*P(511)+
     -   EijI(7,2)*P(512)-EijI(9,2)*P(513))
       F(200)=DCMPLX(FR(200),FI(200))
       P(514) = 2-Is12*s34
       P(515) = s15-s23+s45+2*P(18)
       P(516) = p2sq-p4sq-s23+s34-3*P(52)
       P(517) = p2sq+p4sq
       FR(201) = 2*Is12s45*(3*B013R+B014R+2*(p2sq*Cij123R(2,1)-p4sq*Ci
     -   j145R(2,1)-2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2))-C0145R
     -   *P(15)-(C0123R+Cij123R(1,1))*P(18)-C0134R*P(42)-Cij134R(2,1)
     -   *P(43)))-4*(D02345R-Dij1345R(1,2)+s15*(EE0R+EijR(1,1))-p2sq*
     -   EijR(2,2)-p4sq*EijR(4,2)-6*EijR(11,2)-4*(EijR(22,3)-EijR(24,
     -   3))+Is45*(Cij134R(1,2)-s12*Dij1234R(2,1)-p2sq*Dij1234R(2,2)+
     -   2*Dij1234R(12,3)+Dij1234R(4,2)*P(18))-D01345R*P(39)-(D01234R
     -   +Dij1234R(1,1))*P(40)-(Dij1234R(3,1)+Dij1234R(6,2))*P(41)+Ci
     -   j134R(3,2)*P(44)-EijR(3,1)*P(53)-EijR(4,1)*P(99)-Is12*(Cij14
     -   5R(1,1)-s45*Dij1345R(3,1)+p4sq*Dij1345R(3,2)-2*(Dij1345R(11,
     -   3)-Dij1345R(13,3))-Dij1345R(5,2)*P(137))-Dij1345R(6,2)*P(340
     -   )-(Dij1345R(2,1)+Dij1345R(4,2))*P(342)+EijR(2,1)*P(416)-Dij1
     -   345R(1,1)*P(514)-(-EijR(5,2)+EijR(7,2))*P(515)-(EijR(8,2)-Ei
     -   jR(10,2))*P(516)+EijR(9,2)*P(517))
       FI(201) = 2*Is12s45*(3*B013I+B014I+2*(p2sq*Cij123I(2,1)-p4sq*Ci
     -   j145I(2,1)-2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2))-C0145I
     -   *P(15)-(C0123I+Cij123I(1,1))*P(18)-C0134I*P(42)-Cij134I(2,1)
     -   *P(43)))-4*(D02345I-Dij1345I(1,2)+s15*(EE0I+EijI(1,1))-p2sq*
     -   EijI(2,2)-p4sq*EijI(4,2)-6*EijI(11,2)-4*(EijI(22,3)-EijI(24,
     -   3))+Is45*(Cij134I(1,2)-s12*Dij1234I(2,1)-p2sq*Dij1234I(2,2)+
     -   2*Dij1234I(12,3)+Dij1234I(4,2)*P(18))-D01345I*P(39)-(D01234I
     -   +Dij1234I(1,1))*P(40)-(Dij1234I(3,1)+Dij1234I(6,2))*P(41)+Ci
     -   j134I(3,2)*P(44)-EijI(3,1)*P(53)-EijI(4,1)*P(99)-Is12*(Cij14
     -   5I(1,1)-s45*Dij1345I(3,1)+p4sq*Dij1345I(3,2)-2*(Dij1345I(11,
     -   3)-Dij1345I(13,3))-Dij1345I(5,2)*P(137))-Dij1345I(6,2)*P(340
     -   )-(Dij1345I(2,1)+Dij1345I(4,2))*P(342)+EijI(2,1)*P(416)-Dij1
     -   345I(1,1)*P(514)-(-EijI(5,2)+EijI(7,2))*P(515)-(EijI(8,2)-Ei
     -   jI(10,2))*P(516)+EijI(9,2)*P(517))
       F(201)=DCMPLX(FR(201),FI(201))
       P(518) = 2*p4sq+s45
       P(519) = 3*p2sq+p3sq-s12-s23
       P(520) = p2sq+3*p3sq-p4sq-s23+s34-2*P(52)
       P(521) = p2sq-p3sq-p4sq-s12+s34
       P(522) = s15+s45+2*P(521)
       P(523) = p2sq-3*p3sq+s12+s23+2*P(200)
       P(524) = p2sq+p3sq-p4sq-s23+s34-4*P(52)
       FR(202) = 2*Is12s45*(3*B013R+B014R-2*(C0134R*p3sq-p2sq*Cij123R(
     -   2,1)+p4sq*Cij145R(2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(
     -   4,2))+C0145R*P(15)+(C0123R+Cij123R(1,1))*P(18)+Cij134R(2,1)*
     -   P(43)))-4*(-Dij1345R(1,2)+Dij2345R(3,1)+s15*(EE0R+2*EijR(1,1
     -   )+EijR(1,2))-10*EijR(11,2)-4*(EijR(21,3)-EijR(24,3))+EijR(2,
     -   1)*P(9)-D01345R*P(39)+Cij134R(3,2)*P(44)-2*(Dij2345R(1,1)-Ei
     -   jR(8,2)*P(25)+EijR(2,2)*P(68))-EijR(3,1)*P(83)+Is45*(Cij134R
     -   (1,2)+s23*(D01234R+Dij1234R(1,1))-2*(p2sq*Dij1234R(2,2)+Dij1
     -   234R(7,2)-Dij1234R(11,3))+Dij1234R(6,2)*P(1)+Dij1234R(2,1)*P
     -   (18)-Dij1234R(3,1)*P(45)-Dij1234R(5,2)*P(78)+Dij1234R(4,2)*P
     -   (86))-Is12*(Cij145R(1,1)-s45*Dij1345R(3,1)+p4sq*Dij1345R(3,2
     -   )-2*(Dij1345R(11,3)-Dij1345R(13,3))-Dij1345R(5,2)*P(137))-Di
     -   j1345R(6,2)*P(340)-(Dij1345R(2,1)+Dij1345R(4,2))*P(342)-EijR
     -   (4,1)*P(417)-Dij1345R(1,1)*P(514)-EijR(4,2)*P(518)+EijR(5,2)
     -   *P(519)-EijR(6,2)*P(520)-EijR(7,2)*P(522)+EijR(9,2)*P(523)+E
     -   ijR(10,2)*P(524))
       FI(202) = 2*Is12s45*(3*B013I+B014I-2*(C0134I*p3sq-p2sq*Cij123I(
     -   2,1)+p4sq*Cij145I(2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(
     -   4,2))+C0145I*P(15)+(C0123I+Cij123I(1,1))*P(18)+Cij134I(2,1)*
     -   P(43)))-4*(-Dij1345I(1,2)+Dij2345I(3,1)+s15*(EE0I+2*EijI(1,1
     -   )+EijI(1,2))-10*EijI(11,2)-4*(EijI(21,3)-EijI(24,3))+EijI(2,
     -   1)*P(9)-D01345I*P(39)+Cij134I(3,2)*P(44)-2*(Dij2345I(1,1)-Ei
     -   jI(8,2)*P(25)+EijI(2,2)*P(68))-EijI(3,1)*P(83)+Is45*(Cij134I
     -   (1,2)+s23*(D01234I+Dij1234I(1,1))-2*(p2sq*Dij1234I(2,2)+Dij1
     -   234I(7,2)-Dij1234I(11,3))+Dij1234I(6,2)*P(1)+Dij1234I(2,1)*P
     -   (18)-Dij1234I(3,1)*P(45)-Dij1234I(5,2)*P(78)+Dij1234I(4,2)*P
     -   (86))-Is12*(Cij145I(1,1)-s45*Dij1345I(3,1)+p4sq*Dij1345I(3,2
     -   )-2*(Dij1345I(11,3)-Dij1345I(13,3))-Dij1345I(5,2)*P(137))-Di
     -   j1345I(6,2)*P(340)-(Dij1345I(2,1)+Dij1345I(4,2))*P(342)-EijI
     -   (4,1)*P(417)-Dij1345I(1,1)*P(514)-EijI(4,2)*P(518)+EijI(5,2)
     -   *P(519)-EijI(6,2)*P(520)-EijI(7,2)*P(522)+EijI(9,2)*P(523)+E
     -   ijI(10,2)*P(524))
       F(202)=DCMPLX(FR(202),FI(202))
       P(525) = p2sq-3*p3sq-p4sq-s23+s34-2*P(52)
       P(526) = p4sq+2*s45
       P(527) = s15+3*P(45)
       P(528) = s15+2*P(18)-3*P(36)
       P(529) = 3*p3sq-s45+2*P(10)
       P(530) = p2sq-3*p3sq-s12-s23-2*P(128)
       FR(203) = -4*(-D02345R+Dij2345R(1,1)+Dij2345R(2,1)-2*Dij2345R(3
     -   ,1)+12*EijR(11,2)-4*(EijR(23,3)-EijR(24,3))-EijR(5,2)*P(18)+
     -   Is45*(Cij134R(1,1)+Cij134R(3,2)-p2sq*(Dij1234R(2,1)+Dij1234R
     -   (6,2))+2*(Dij1234R(7,2)+Dij1234R(13,3))+Dij1234R(5,2)*P(18))
     -   -Dij1234R(3,2)*P(41)+(Cij134R(2,1)+Cij134R(2,2))*P(44)+EijR(
     -   2,2)*P(68)-Dij1234R(3,1)*P(91)-(-EijR(3,1)+EijR(4,1))*P(99)-
     -   Is12*(C0134R+s34*Dij1345R(1,1)-s34*Dij1345R(2,1)-p4sq*Dij134
     -   5R(3,2)-4*Dij1345R(7,2)-2*(Dij1345R(12,3)-Dij1345R(13,3))+(D
     -   ij1345R(4,2)-Dij1345R(5,2))*P(112)-Dij1345R(2,2)*P(203)+Dij1
     -   345R(6,2)*P(378))+EijR(9,2)*P(513)-EijR(3,2)*P(525)+EijR(4,2
     -   )*P(526)+EijR(6,2)*P(527)-EijR(7,2)*P(528)-EijR(8,2)*P(529)+
     -   EijR(10,2)*P(530))
       FI(203) = -4*(-D02345I+Dij2345I(1,1)+Dij2345I(2,1)-2*Dij2345I(3
     -   ,1)+12*EijI(11,2)-4*(EijI(23,3)-EijI(24,3))-EijI(5,2)*P(18)+
     -   Is45*(Cij134I(1,1)+Cij134I(3,2)-p2sq*(Dij1234I(2,1)+Dij1234I
     -   (6,2))+2*(Dij1234I(7,2)+Dij1234I(13,3))+Dij1234I(5,2)*P(18))
     -   -Dij1234I(3,2)*P(41)+(Cij134I(2,1)+Cij134I(2,2))*P(44)+EijI(
     -   2,2)*P(68)-Dij1234I(3,1)*P(91)-(-EijI(3,1)+EijI(4,1))*P(99)-
     -   Is12*(C0134I+s34*Dij1345I(1,1)-s34*Dij1345I(2,1)-p4sq*Dij134
     -   5I(3,2)-4*Dij1345I(7,2)-2*(Dij1345I(12,3)-Dij1345I(13,3))+(D
     -   ij1345I(4,2)-Dij1345I(5,2))*P(112)-Dij1345I(2,2)*P(203)+Dij1
     -   345I(6,2)*P(378))+EijI(9,2)*P(513)-EijI(3,2)*P(525)+EijI(4,2
     -   )*P(526)+EijI(6,2)*P(527)-EijI(7,2)*P(528)-EijI(8,2)*P(529)+
     -   EijI(10,2)*P(530))
       F(203)=DCMPLX(FR(203),FI(203))
       P(531) = p2sq-2*p3sq-p4sq-s12-s23+s34+s45
       FR(204) = 4*(D01345R-Dij2345R(1,1)+2*Dij2345R(3,1)-s12*EijR(2,1
     -   )+p3sq*EijR(3,1)-4*(EijR(11,2)+EijR(24,3))+EijR(5,2)*P(18)+E
     -   ijR(8,2)*P(25)-EijR(2,2)*P(68)+EijR(4,1)*P(140)+Is12*(C0134R
     -   -p4sq*(Dij1345R(3,1)+Dij1345R(3,2))-2*(Dij1345R(7,2)-Dij1345
     -   R(13,3))-Dij1345R(5,2)*P(112)+Dij1345R(2,1)*P(128)+Dij1345R(
     -   1,1)*P(136)+Dij1345R(6,2)*P(203))-EijR(9,2)*P(513)-EijR(4,2)
     -   *P(526)+EijR(7,2)*P(528)-EijR(10,2)*P(531))
       FI(204) = 4*(D01345I-Dij2345I(1,1)+2*Dij2345I(3,1)-s12*EijI(2,1
     -   )+p3sq*EijI(3,1)-4*(EijI(11,2)+EijI(24,3))+EijI(5,2)*P(18)+E
     -   ijI(8,2)*P(25)-EijI(2,2)*P(68)+EijI(4,1)*P(140)+Is12*(C0134I
     -   -p4sq*(Dij1345I(3,1)+Dij1345I(3,2))-2*(Dij1345I(7,2)-Dij1345
     -   I(13,3))-Dij1345I(5,2)*P(112)+Dij1345I(2,1)*P(128)+Dij1345I(
     -   1,1)*P(136)+Dij1345I(6,2)*P(203))-EijI(9,2)*P(513)-EijI(4,2)
     -   *P(526)+EijI(7,2)*P(528)-EijI(10,2)*P(531))
       F(204)=DCMPLX(FR(204),FI(204))
       P(532) = 2*p2sq-p4sq-4*s12+s45+3*P(6)
       P(533) = p3sq-p4sq-s12-s23-s45
       FR(205) = -4*(-Dij2345R(1,1)+Dij2345R(3,1)-Dij2345R(3,2)+Dij234
     -   5R(5,2)-s12*EijR(2,1)+p3sq*EijR(3,1)-4*(EijR(11,2)-EijR(24,3
     -   ))+EijR(5,2)*P(18)+EijR(8,2)*P(25)+(EijR(14,3)-EijR(18,3))*P
     -   (36)-EijR(4,1)*P(42)+(-EijR(16,3)+EijR(20,3))*P(43)-EijR(2,2
     -   )*P(68)+EijR(4,3)*P(174)-EijR(10,3)*P(210)+EijR(10,2)*P(495)
     -   +EijR(7,2)*P(512)-EijR(4,2)*P(526)-EijR(9,2)*P(532)+EijR(15,
     -   3)*P(533))
       FI(205) = -4*(-Dij2345I(1,1)+Dij2345I(3,1)-Dij2345I(3,2)+Dij234
     -   5I(5,2)-s12*EijI(2,1)+p3sq*EijI(3,1)-4*(EijI(11,2)-EijI(24,3
     -   ))+EijI(5,2)*P(18)+EijI(8,2)*P(25)+(EijI(14,3)-EijI(18,3))*P
     -   (36)-EijI(4,1)*P(42)+(-EijI(16,3)+EijI(20,3))*P(43)-EijI(2,2
     -   )*P(68)+EijI(4,3)*P(174)-EijI(10,3)*P(210)+EijI(10,2)*P(495)
     -   +EijI(7,2)*P(512)-EijI(4,2)*P(526)-EijI(9,2)*P(532)+EijI(15,
     -   3)*P(533))
       F(205)=DCMPLX(FR(205),FI(205))
       P(534) = p3sq-p4sq-3*s12+s34+2*P(214)
       P(535) = -p4sq-s45+2*P(210)
       P(536) = p3sq-s12-s23-2*P(174)
       FR(206) = -4*(-Dij2345R(1,2)-Dij2345R(3,2)-s45*EijR(4,2)-p3sq*E
     -   ijR(8,2)+p3sq*EijR(10,2)-6*(EijR(22,3)-EijR(24,3))-EijR(9,2)
     -   *P(12)+(EijR(8,3)+EijR(14,3)-2*EijR(18,3))*P(36)-(EijR(9,3)+
     -   EijR(16,3)-2*EijR(20,3))*P(43)+2*(Dij2345R(5,2)-EijR(5,2)*P(
     -   45)-EijR(7,2)*P(46))+Is45*(Cij123R(2,1)+Cij123R(2,2)+p3sq*(D
     -   ij1234R(2,1)-Dij1234R(3,1)-Dij1234R(6,2))+2*(Dij1234R(7,2)+D
     -   ij1234R(12,3))-Dij1234R(4,2)*P(45)+Dij1234R(2,2)*P(78))+EijR
     -   (4,3)*P(174)+EijR(2,3)*P(210)-EijR(2,1)*P(235)-EijR(4,1)*P(2
     -   55)+EijR(2,2)*P(534)-EijR(10,3)*P(535)+EijR(15,3)*P(536))
       FI(206) = -4*(-Dij2345I(1,2)-Dij2345I(3,2)-s45*EijI(4,2)-p3sq*E
     -   ijI(8,2)+p3sq*EijI(10,2)-6*(EijI(22,3)-EijI(24,3))-EijI(9,2)
     -   *P(12)+(EijI(8,3)+EijI(14,3)-2*EijI(18,3))*P(36)-(EijI(9,3)+
     -   EijI(16,3)-2*EijI(20,3))*P(43)+2*(Dij2345I(5,2)-EijI(5,2)*P(
     -   45)-EijI(7,2)*P(46))+Is45*(Cij123I(2,1)+Cij123I(2,2)+p3sq*(D
     -   ij1234I(2,1)-Dij1234I(3,1)-Dij1234I(6,2))+2*(Dij1234I(7,2)+D
     -   ij1234I(12,3))-Dij1234I(4,2)*P(45)+Dij1234I(2,2)*P(78))+EijI
     -   (4,3)*P(174)+EijI(2,3)*P(210)-EijI(2,1)*P(235)-EijI(4,1)*P(2
     -   55)+EijI(2,2)*P(534)-EijI(10,3)*P(535)+EijI(15,3)*P(536))
       F(206)=DCMPLX(FR(206),FI(206))
       P(537) = 2-Is45*s23
       P(538) = 2*p2sq+p3sq+s12-s23
       P(539) = 2*p2sq+p3sq+p4sq+s12-s23-s34-s45
       P(540) = p2sq-p3sq-p4sq-s12+s34+s45
       P(541) = s23-2*P(540)
       P(542) = p2sq-p3sq-p4sq-s12
       P(543) = -s23+3*P(200)+2*P(542)
       P(544) = p4sq-s23+2*s45
       P(545) = p3sq-p4sq-s12-2*s45
       FR(207) = -4*(-D01234R-Dij2345R(3,2)+Dij2345R(5,2)+p3sq*EijR(3,
     -   1)+p4sq*EijR(4,2)+3*p3sq*EijR(6,2)+4*(EijR(11,2)-EijR(21,3))
     -   +2*(Dij2345R(1,1)-Dij2345R(3,1)-EijR(22,3))+6*EijR(24,3)-(-E
     -   ijR(5,3)+EijR(7,3))*P(36)-Dij1234R(1,2)*P(40)-(EijR(16,3)+Ei
     -   jR(17,3)-EijR(19,3)-EijR(20,3))*P(43)+EijR(10,2)*P(52)-EijR(
     -   8,2)*P(169)+EijR(4,3)*P(174)-(-EijR(8,3)+EijR(10,3))*P(210)+
     -   Is45*(C0123R-C0134R+Cij123R(1,1)+Cij123R(2,1)+Cij123R(3,2)+6
     -   *Dij1234R(7,2)+2*(p2sq*Dij1234R(2,2)+Dij1234R(11,3))-Dij1234
     -   R(6,2)*P(1)+Dij1234R(5,2)*P(10)-Dij1234R(3,1)*P(22)-Dij1234R
     -   (2,1)*P(238)-Dij1234R(4,2)*P(239))-EijR(2,1)*P(400)+EijR(4,1
     -   )*P(401)+EijR(15,3)*P(533)-Dij1234R(1,1)*P(537)+EijR(2,2)*P(
     -   538)-EijR(5,2)*P(539)-EijR(7,2)*P(541)-EijR(9,2)*P(543)-EijR
     -   (14,3)*P(544)-EijR(18,3)*P(545))
       FI(207) = -4*(-D01234I-Dij2345I(3,2)+Dij2345I(5,2)+p3sq*EijI(3,
     -   1)+p4sq*EijI(4,2)+3*p3sq*EijI(6,2)+4*(EijI(11,2)-EijI(21,3))
     -   +2*(Dij2345I(1,1)-Dij2345I(3,1)-EijI(22,3))+6*EijI(24,3)-(-E
     -   ijI(5,3)+EijI(7,3))*P(36)-Dij1234I(1,2)*P(40)-(EijI(16,3)+Ei
     -   jI(17,3)-EijI(19,3)-EijI(20,3))*P(43)+EijI(10,2)*P(52)-EijI(
     -   8,2)*P(169)+EijI(4,3)*P(174)-(-EijI(8,3)+EijI(10,3))*P(210)+
     -   Is45*(C0123I-C0134I+Cij123I(1,1)+Cij123I(2,1)+Cij123I(3,2)+6
     -   *Dij1234I(7,2)+2*(p2sq*Dij1234I(2,2)+Dij1234I(11,3))-Dij1234
     -   I(6,2)*P(1)+Dij1234I(5,2)*P(10)-Dij1234I(3,1)*P(22)-Dij1234I
     -   (2,1)*P(238)-Dij1234I(4,2)*P(239))-EijI(2,1)*P(400)+EijI(4,1
     -   )*P(401)+EijI(15,3)*P(533)-Dij1234I(1,1)*P(537)+EijI(2,2)*P(
     -   538)-EijI(5,2)*P(539)-EijI(7,2)*P(541)-EijI(9,2)*P(543)-EijI
     -   (14,3)*P(544)-EijI(18,3)*P(545))
       F(207)=DCMPLX(FR(207),FI(207))
       P(546) = -p4sq-2*s12+s34+s45+3*P(1)
       P(547) = 3*p3sq+p4sq-s12-s34+2*s45
       P(548) = p3sq+p4sq-s12+2*s45
       P(549) = p4sq+s23+2*s45
       FR(208) = -4*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)-Dij2345R(3,2)
     -   -Dij2345R(4,2)+Dij2345R(5,2)+Dij2345R(6,2)-2*(6*EijR(11,2)+E
     -   ijR(22,3)+2*EijR(23,3)-3*EijR(24,3))+EijR(5,2)*P(18)+(EijR(1
     -   4,3)+EijR(17,3)-EijR(18,3)-EijR(19,3))*P(36)+(-EijR(12,3)+Ei
     -   jR(13,3))*P(43)-3*(p3sq*EijR(3,2)+EijR(6,2)*P(45))-EijR(2,2)
     -   *P(68)+Is45*(C0134R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,2)-2*
     -   (Dij1234R(7,2)-Dij1234R(13,3))-Dij1234R(3,1)*P(10)-Dij1234R(
     -   5,2)*P(45)+Dij1234R(6,2)*P(78))+EijR(4,3)*P(174)+(EijR(9,3)-
     -   EijR(10,3))*P(210)-EijR(3,1)*P(235)-EijR(4,1)*P(255)+EijR(7,
     -   2)*P(512)-EijR(4,2)*P(526)-EijR(9,2)*P(532)+EijR(15,3)*P(533
     -   )+EijR(8,2)*P(546)+EijR(10,2)*P(547)-EijR(16,3)*P(548)+EijR(
     -   20,3)*P(549))
       FI(208) = -4*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)-Dij2345I(3,2)
     -   -Dij2345I(4,2)+Dij2345I(5,2)+Dij2345I(6,2)-2*(6*EijI(11,2)+E
     -   ijI(22,3)+2*EijI(23,3)-3*EijI(24,3))+EijI(5,2)*P(18)+(EijI(1
     -   4,3)+EijI(17,3)-EijI(18,3)-EijI(19,3))*P(36)+(-EijI(12,3)+Ei
     -   jI(13,3))*P(43)-3*(p3sq*EijI(3,2)+EijI(6,2)*P(45))-EijI(2,2)
     -   *P(68)+Is45*(C0134I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,2)-2*
     -   (Dij1234I(7,2)-Dij1234I(13,3))-Dij1234I(3,1)*P(10)-Dij1234I(
     -   5,2)*P(45)+Dij1234I(6,2)*P(78))+EijI(4,3)*P(174)+(EijI(9,3)-
     -   EijI(10,3))*P(210)-EijI(3,1)*P(235)-EijI(4,1)*P(255)+EijI(7,
     -   2)*P(512)-EijI(4,2)*P(526)-EijI(9,2)*P(532)+EijI(15,3)*P(533
     -   )+EijI(8,2)*P(546)+EijI(10,2)*P(547)-EijI(16,3)*P(548)+EijI(
     -   20,3)*P(549))
       F(208)=DCMPLX(FR(208),FI(208))
       FR(209) = 8*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,3)+Dij
     -   1345R(4,3)-Dij1345R(5,3)+Dij1345R(9,3)+2*(Dij1345R(3,2)+Dij1
     -   345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij1345R
     -   (10,3))+s12*(EijR(3,1)-EijR(4,1)+EijR(4,3)-EijR(6,3)+EijR(7,
     -   3)-EijR(16,3)-2*(EijR(14,3)-EijR(19,3))))
       FI(209) = 8*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,3)+Dij
     -   1345I(4,3)-Dij1345I(5,3)+Dij1345I(9,3)+2*(Dij1345I(3,2)+Dij1
     -   345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij1345I
     -   (10,3))+s12*(EijI(3,1)-EijI(4,1)+EijI(4,3)-EijI(6,3)+EijI(7,
     -   3)-EijI(16,3)-2*(EijI(14,3)-EijI(19,3))))
       F(209)=DCMPLX(FR(209),FI(209))
       FR(210) = -8*Is12*(-Dij1345R(2,1)+Dij1345R(3,1)+Dij1345R(3,3)-D
     -   ij1345R(4,3)+Dij1345R(5,3)-Dij1345R(9,3)-2*(Dij1345R(3,2)+Di
     -   j1345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij134
     -   5R(10,3))+s12*(EijR(2,1)-EijR(3,1)+EijR(4,2)-EijR(4,3)+EijR(
     -   5,2)+EijR(6,2)-2*EijR(7,2)-EijR(8,2)+EijR(14,3)+EijR(15,3)+E
     -   ijR(16,3)+EijR(17,3)-EijR(18,3)-EijR(19,3)-EijR(20,3)))
       FI(210) = -8*Is12*(-Dij1345I(2,1)+Dij1345I(3,1)+Dij1345I(3,3)-D
     -   ij1345I(4,3)+Dij1345I(5,3)-Dij1345I(9,3)-2*(Dij1345I(3,2)+Di
     -   j1345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij134
     -   5I(10,3))+s12*(EijI(2,1)-EijI(3,1)+EijI(4,2)-EijI(4,3)+EijI(
     -   5,2)+EijI(6,2)-2*EijI(7,2)-EijI(8,2)+EijI(14,3)+EijI(15,3)+E
     -   ijI(16,3)+EijI(17,3)-EijI(18,3)-EijI(19,3)-EijI(20,3)))
       F(210)=DCMPLX(FR(210),FI(210))
       FR(211) = 8*Is12*(Dij1345R(2,2)+Dij1345R(3,2)-Dij1345R(3,3)+Dij
     -   1345R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)-2*(Dij1345R(6,2)-Dij1
     -   345R(9,3)+Dij1345R(10,3))+s12*(EijR(3,2)-EijR(4,2)+EijR(4,3)
     -   -EijR(11,3)+EijR(13,3)-EijR(14,3)-2*(EijR(6,2)-EijR(7,2)+Eij
     -   R(16,3)-EijR(19,3))))
       FI(211) = 8*Is12*(Dij1345I(2,2)+Dij1345I(3,2)-Dij1345I(3,3)+Dij
     -   1345I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)-2*(Dij1345I(6,2)-Dij1
     -   345I(9,3)+Dij1345I(10,3))+s12*(EijI(3,2)-EijI(4,2)+EijI(4,3)
     -   -EijI(11,3)+EijI(13,3)-EijI(14,3)-2*(EijI(6,2)-EijI(7,2)+Eij
     -   I(16,3)-EijI(19,3))))
       F(211)=DCMPLX(FR(211),FI(211))
       FR(212) = -8*Is12*(Dij1345R(3,3)-Dij1345R(5,2)+Dij1345R(6,2)-Di
     -   j1345R(7,3)-Dij1345R(9,3)+Dij1345R(10,3)+s12*(2*EijR(4,2)-Ei
     -   jR(4,3)-3*EijR(7,2)+EijR(10,2)+EijR(14,3)+EijR(16,3)-EijR(19
     -   ,3)))
       FI(212) = -8*Is12*(Dij1345I(3,3)-Dij1345I(5,2)+Dij1345I(6,2)-Di
     -   j1345I(7,3)-Dij1345I(9,3)+Dij1345I(10,3)+s12*(2*EijI(4,2)-Ei
     -   jI(4,3)-3*EijI(7,2)+EijI(10,2)+EijI(14,3)+EijI(16,3)-EijI(19
     -   ,3)))
       F(212)=DCMPLX(FR(212),FI(212))
       FR(213) = 8*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,3)+Dij
     -   1345R(4,3)-Dij1345R(5,3)+Dij1345R(9,3)+2*(Dij1345R(3,2)+Dij1
     -   345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij1345R
     -   (10,3))+s12*(EijR(2,1)-EijR(4,1)+EijR(4,3)+EijR(5,2)-EijR(7,
     -   2)-EijR(8,2)+EijR(10,2)-EijR(14,3)-EijR(15,3)-EijR(16,3)-Eij
     -   R(17,3)+EijR(18,3)+EijR(19,3)+EijR(20,3)))
       FI(213) = 8*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,3)+Dij
     -   1345I(4,3)-Dij1345I(5,3)+Dij1345I(9,3)+2*(Dij1345I(3,2)+Dij1
     -   345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij1345I
     -   (10,3))+s12*(EijI(2,1)-EijI(4,1)+EijI(4,3)+EijI(5,2)-EijI(7,
     -   2)-EijI(8,2)+EijI(10,2)-EijI(14,3)-EijI(15,3)-EijI(16,3)-Eij
     -   I(17,3)+EijI(18,3)+EijI(19,3)+EijI(20,3)))
       F(213)=DCMPLX(FR(213),FI(213))
       FR(214) = 8*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,3)+Dij
     -   1345R(4,3)-Dij1345R(5,3)+Dij1345R(9,3)+2*(Dij1345R(3,2)+Dij1
     -   345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij1345R
     -   (10,3))+s12*(-EijR(4,2)+EijR(4,3)-EijR(8,2)+EijR(9,2)-EijR(9
     -   ,3)+EijR(10,2)+EijR(10,3)-EijR(16,3)-2*(EijR(15,3)-EijR(20,3
     -   ))))
       FI(214) = 8*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,3)+Dij
     -   1345I(4,3)-Dij1345I(5,3)+Dij1345I(9,3)+2*(Dij1345I(3,2)+Dij1
     -   345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij1345I
     -   (10,3))+s12*(-EijI(4,2)+EijI(4,3)-EijI(8,2)+EijI(9,2)-EijI(9
     -   ,3)+EijI(10,2)+EijI(10,3)-EijI(16,3)-2*(EijI(15,3)-EijI(20,3
     -   ))))
       F(214)=DCMPLX(FR(214),FI(214))
       FR(215) = -8*Is12*(Dij1345R(3,3)-Dij1345R(5,2)+Dij1345R(6,2)-Di
     -   j1345R(7,3)-Dij1345R(9,3)+Dij1345R(10,3)+s12*(-EijR(4,3)+2*(
     -   EijR(4,2)-EijR(9,2))+EijR(15,3)+EijR(16,3)-EijR(20,3)))
       FI(215) = -8*Is12*(Dij1345I(3,3)-Dij1345I(5,2)+Dij1345I(6,2)-Di
     -   j1345I(7,3)-Dij1345I(9,3)+Dij1345I(10,3)+s12*(-EijI(4,3)+2*(
     -   EijI(4,2)-EijI(9,2))+EijI(15,3)+EijI(16,3)-EijI(20,3)))
       F(215)=DCMPLX(FR(215),FI(215))
       FR(216) = -8*Is12*(Dij1345R(3,3)+Dij1345R(8,3)-2*Dij1345R(9,3)-
     -   s12*(EijR(4,3)+EijR(13,3)-2*(EijR(4,2)-EijR(10,2)+EijR(16,3)
     -   )))
       FI(216) = -8*Is12*(Dij1345I(3,3)+Dij1345I(8,3)-2*Dij1345I(9,3)-
     -   s12*(EijI(4,3)+EijI(13,3)-2*(EijI(4,2)-EijI(10,2)+EijI(16,3)
     -   )))
       F(216)=DCMPLX(FR(216),FI(216))
       FR(217) = -8*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,2)+Di
     -   j1345R(3,3)+Dij1345R(4,2)-Dij1345R(7,3)-Dij1345R(9,3)+Dij134
     -   5R(10,3)+s12*(EijR(3,1)-EijR(4,1)-EijR(4,2)-EijR(4,3)+EijR(6
     -   ,2)+2*(EijR(7,2)-EijR(10,2))+EijR(14,3)+EijR(16,3)-EijR(19,3
     -   )))
       FI(217) = -8*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,2)+Di
     -   j1345I(3,3)+Dij1345I(4,2)-Dij1345I(7,3)-Dij1345I(9,3)+Dij134
     -   5I(10,3)+s12*(EijI(3,1)-EijI(4,1)-EijI(4,2)-EijI(4,3)+EijI(6
     -   ,2)+2*(EijI(7,2)-EijI(10,2))+EijI(14,3)+EijI(16,3)-EijI(19,3
     -   )))
       F(217)=DCMPLX(FR(217),FI(217))
       FR(218) = -8*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,2)+Di
     -   j1345R(3,3)+Dij1345R(4,2)-Dij1345R(7,3)-Dij1345R(9,3)+Dij134
     -   5R(10,3)+s12*(EijR(3,1)-EijR(4,1)-EijR(4,3)+EijR(8,2)+EijR(9
     -   ,2)-2*EijR(10,2)+EijR(15,3)+EijR(16,3)-EijR(20,3)))
       FI(218) = -8*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,2)+Di
     -   j1345I(3,3)+Dij1345I(4,2)-Dij1345I(7,3)-Dij1345I(9,3)+Dij134
     -   5I(10,3)+s12*(EijI(3,1)-EijI(4,1)-EijI(4,3)+EijI(8,2)+EijI(9
     -   ,2)-2*EijI(10,2)+EijI(15,3)+EijI(16,3)-EijI(20,3)))
       F(218)=DCMPLX(FR(218),FI(218))
       FR(219) = Is12*(Dij1345R(2,2)+Dij1345R(3,2)-Dij1345R(3,3)+Dij13
     -   45R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)-2*(Dij1345R(6,2)-Dij134
     -   5R(9,3)+Dij1345R(10,3))+s12*(-EijR(4,2)+EijR(4,3)-EijR(8,2)+
     -   EijR(9,2)+EijR(10,2)-EijR(12,3)+EijR(13,3)-EijR(15,3)-2*(Eij
     -   R(16,3)-EijR(20,3))))
       FI(219) = Is12*(Dij1345I(2,2)+Dij1345I(3,2)-Dij1345I(3,3)+Dij13
     -   45I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)-2*(Dij1345I(6,2)-Dij134
     -   5I(9,3)+Dij1345I(10,3))+s12*(-EijI(4,2)+EijI(4,3)-EijI(8,2)+
     -   EijI(9,2)+EijI(10,2)-EijI(12,3)+EijI(13,3)-EijI(15,3)-2*(Eij
     -   I(16,3)-EijI(20,3))))
       F(219)=DCMPLX(FR(219),FI(219))
       FR(220) = Is12*(Dij1345R(2,1)+Dij1345R(2,2)-Dij1345R(3,1)-Dij13
     -   45R(3,3)+Dij1345R(4,2)-Dij1345R(5,2)-3*Dij1345R(6,2)+Dij1345
     -   R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)+2*(Dij1345R(3,2)+Dij1345R
     -   (9,3)-Dij1345R(10,3))+s12*(EijR(3,1)-EijR(3,2)-EijR(4,1)+Eij
     -   R(4,3)+EijR(6,2)-EijR(7,2)+EijR(10,2)-EijR(11,3)+EijR(13,3)-
     -   EijR(14,3)-2*(EijR(16,3)-EijR(19,3))))
       FI(220) = Is12*(Dij1345I(2,1)+Dij1345I(2,2)-Dij1345I(3,1)-Dij13
     -   45I(3,3)+Dij1345I(4,2)-Dij1345I(5,2)-3*Dij1345I(6,2)+Dij1345
     -   I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)+2*(Dij1345I(3,2)+Dij1345I
     -   (9,3)-Dij1345I(10,3))+s12*(EijI(3,1)-EijI(3,2)-EijI(4,1)+Eij
     -   I(4,3)+EijI(6,2)-EijI(7,2)+EijI(10,2)-EijI(11,3)+EijI(13,3)-
     -   EijI(14,3)-2*(EijI(16,3)-EijI(19,3))))
       F(220)=DCMPLX(FR(220),FI(220))
       FR(221) = Is12*(Dij1345R(2,1)+Dij1345R(2,2)-Dij1345R(3,1)-Dij13
     -   45R(3,3)+Dij1345R(4,2)-Dij1345R(5,2)-3*Dij1345R(6,2)+Dij1345
     -   R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)+2*(Dij1345R(3,2)+Dij1345R
     -   (9,3)-Dij1345R(10,3))-s12*(EijR(3,2)+EijR(4,2)-EijR(4,3)+Eij
     -   R(12,3)-EijR(13,3)+EijR(15,3)-2*(EijR(10,2)-EijR(16,3)+EijR(
     -   20,3))))
       FI(221) = Is12*(Dij1345I(2,1)+Dij1345I(2,2)-Dij1345I(3,1)-Dij13
     -   45I(3,3)+Dij1345I(4,2)-Dij1345I(5,2)-3*Dij1345I(6,2)+Dij1345
     -   I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)+2*(Dij1345I(3,2)+Dij1345I
     -   (9,3)-Dij1345I(10,3))-s12*(EijI(3,2)+EijI(4,2)-EijI(4,3)+Eij
     -   I(12,3)-EijI(13,3)+EijI(15,3)-2*(EijI(10,2)-EijI(16,3)+EijI(
     -   20,3))))
       F(221)=DCMPLX(FR(221),FI(221))
       FR(222) = Is12*(Dij1345R(2,2)+Dij1345R(2,3)+Dij1345R(3,2)-Dij13
     -   45R(3,3)-2*Dij1345R(6,2)-3*(Dij1345R(8,3)-Dij1345R(9,3))-s12
     -   *(EijR(3,2)+EijR(3,3)+EijR(4,2)-EijR(4,3)-2*EijR(10,2)-3*(Ei
     -   jR(13,3)-EijR(16,3))))
       FI(222) = Is12*(Dij1345I(2,2)+Dij1345I(2,3)+Dij1345I(3,2)-Dij13
     -   45I(3,3)-2*Dij1345I(6,2)-3*(Dij1345I(8,3)-Dij1345I(9,3))-s12
     -   *(EijI(3,2)+EijI(3,3)+EijI(4,2)-EijI(4,3)-2*EijI(10,2)-3*(Ei
     -   jI(13,3)-EijI(16,3))))
       F(222)=DCMPLX(FR(222),FI(222))
       FR(223) = Is12*(Dij1345R(2,1)+Dij1345R(2,2)-Dij1345R(3,1)+Dij13
     -   45R(3,3)-Dij1345R(6,2)+Dij1345R(8,3)-2*Dij1345R(9,3)-s12*(-E
     -   ijR(3,1)-EijR(3,2)+EijR(4,1)+EijR(4,3)+EijR(10,2)+EijR(13,3)
     -   -2*EijR(16,3)))
       FI(223) = Is12*(Dij1345I(2,1)+Dij1345I(2,2)-Dij1345I(3,1)+Dij13
     -   45I(3,3)-Dij1345I(6,2)+Dij1345I(8,3)-2*Dij1345I(9,3)-s12*(-E
     -   ijI(3,1)-EijI(3,2)+EijI(4,1)+EijI(4,3)+EijI(10,2)+EijI(13,3)
     -   -2*EijI(16,3)))
       F(223)=DCMPLX(FR(223),FI(223))
       FR(224) = Is12*(Dij1345R(3,2)+Dij1345R(3,3)-Dij1345R(6,2)-Dij13
     -   45R(9,3)+s12*EijR(4,2)-s12*EijR(4,3)-s12*EijR(10,2)+s12*EijR
     -   (16,3))
       FI(224) = Is12*(Dij1345I(3,2)+Dij1345I(3,3)-Dij1345I(6,2)-Dij13
     -   45I(9,3)+s12*EijI(4,2)-s12*EijI(4,3)-s12*EijI(10,2)+s12*EijI
     -   (16,3))
       F(224)=DCMPLX(FR(224),FI(224))
       FR(225) = 2*Is12s45*(3*B013R+B014R-2*(C0134R*p3sq-p2sq*Cij123R(
     -   2,1)+p4sq*Cij145R(2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(
     -   4,2))+C0145R*P(15)+(C0123R+Cij123R(1,1))*P(18)+Cij134R(2,1)*
     -   P(147)))+4*(Is12*(Cij145R(1,1)+2*(Dij1345R(7,2)+Dij1345R(12,
     -   3)-Dij1345R(13,3)))-s15*(EE0R+EijR(1,1))+(D01345R+Dij1345R(1
     -   ,1))*P(39)+(D01234R+Dij1234R(1,1))*P(40)+Is45*(2*(Cij134R(1,
     -   1)+p2sq*Dij1234R(2,2))+6*Dij1234R(7,2)+Dij1234R(3,2)*P(1)-(D
     -   ij1234R(2,1)+Dij1234R(4,2)-Dij1234R(5,2))*P(18)+Dij1234R(3,1
     -   )*P(66)-Dij1234R(6,2)*P(157))-2*(D02345R-Dij2345R(1,1)+Dij23
     -   45R(3,1)-2*p3sq*EijR(3,2)-9*EijR(11,2)+EijR(23,3)-EijR(24,3)
     -   +EijR(5,2)*P(18)-EijR(7,2)*P(36)-EijR(6,2)*P(45)-EijR(2,2)*P
     -   (68)-EijR(4,2)*P(174)+EijR(10,2)*P(407)+EijR(9,2)*P(411)+Eij
     -   R(8,2)*P(413)))
       FI(225) = 2*Is12s45*(3*B013I+B014I-2*(C0134I*p3sq-p2sq*Cij123I(
     -   2,1)+p4sq*Cij145I(2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(
     -   4,2))+C0145I*P(15)+(C0123I+Cij123I(1,1))*P(18)+Cij134I(2,1)*
     -   P(147)))+4*(Is12*(Cij145I(1,1)+2*(Dij1345I(7,2)+Dij1345I(12,
     -   3)-Dij1345I(13,3)))-s15*(EE0I+EijI(1,1))+(D01345I+Dij1345I(1
     -   ,1))*P(39)+(D01234I+Dij1234I(1,1))*P(40)+Is45*(2*(Cij134I(1,
     -   1)+p2sq*Dij1234I(2,2))+6*Dij1234I(7,2)+Dij1234I(3,2)*P(1)-(D
     -   ij1234I(2,1)+Dij1234I(4,2)-Dij1234I(5,2))*P(18)+Dij1234I(3,1
     -   )*P(66)-Dij1234I(6,2)*P(157))-2*(D02345I-Dij2345I(1,1)+Dij23
     -   45I(3,1)-2*p3sq*EijI(3,2)-9*EijI(11,2)+EijI(23,3)-EijI(24,3)
     -   +EijI(5,2)*P(18)-EijI(7,2)*P(36)-EijI(6,2)*P(45)-EijI(2,2)*P
     -   (68)-EijI(4,2)*P(174)+EijI(10,2)*P(407)+EijI(9,2)*P(411)+Eij
     -   I(8,2)*P(413)))
       F(225)=DCMPLX(FR(225),FI(225))
       FR(226) = -8*(EijR(2,2)-EijR(4,2)+EijR(4,3)-EijR(8,3)+EijR(10,3
     -   )-EijR(14,3)-2*(EijR(5,2)-EijR(7,2)+EijR(15,3)-EijR(18,3)))
       FI(226) = -8*(EijI(2,2)-EijI(4,2)+EijI(4,3)-EijI(8,3)+EijI(10,3
     -   )-EijI(14,3)-2*(EijI(5,2)-EijI(7,2)+EijI(15,3)-EijI(18,3)))
       F(226)=DCMPLX(FR(226),FI(226))
       FR(227) = 8*(EijR(4,2)-EijR(4,3)+2*(EijR(6,2)-EijR(7,2))-EijR(8
     -   ,2)+EijR(9,2)-EijR(10,2)+EijR(14,3)+EijR(15,3)+EijR(16,3)+Ei
     -   jR(17,3)-EijR(18,3)-EijR(19,3)-EijR(20,3))
       FI(227) = 8*(EijI(4,2)-EijI(4,3)+2*(EijI(6,2)-EijI(7,2))-EijI(8
     -   ,2)+EijI(9,2)-EijI(10,2)+EijI(14,3)+EijI(15,3)+EijI(16,3)+Ei
     -   jI(17,3)-EijI(18,3)-EijI(19,3)-EijI(20,3))
       F(227)=DCMPLX(FR(227),FI(227))
       FR(228) = 8*(2*EijR(4,2)-EijR(4,3)-3*EijR(7,2)+EijR(9,2)+EijR(1
     -   4,3)+EijR(15,3)-EijR(18,3))
       FI(228) = 8*(2*EijI(4,2)-EijI(4,3)-3*EijI(7,2)+EijI(9,2)+EijI(1
     -   4,3)+EijI(15,3)-EijI(18,3))
       F(228)=DCMPLX(FR(228),FI(228))
       FR(229) = 8*(-EijR(4,3)-EijR(10,3)+2*(EijR(4,2)-EijR(9,2)+EijR(
     -   15,3)))
       FI(229) = 8*(-EijI(4,3)-EijI(10,3)+2*(EijI(4,2)-EijI(9,2)+EijI(
     -   15,3)))
       F(229)=DCMPLX(FR(229),FI(229))
       FR(230) = 8*(-EijR(4,3)+2*(EijR(4,2)-EijR(10,2))+EijR(15,3)+Eij
     -   R(16,3)-EijR(20,3))
       FI(230) = 8*(-EijI(4,3)+2*(EijI(4,2)-EijI(10,2))+EijI(15,3)+Eij
     -   I(16,3)-EijI(20,3))
       F(230)=DCMPLX(FR(230),FI(230))
       FR(231) = 8*(EijR(2,1)-EijR(4,1)-EijR(4,2)-EijR(4,3)+EijR(5,2)+
     -   2*(EijR(7,2)-EijR(9,2))+EijR(14,3)+EijR(15,3)-EijR(18,3))
       FI(231) = 8*(EijI(2,1)-EijI(4,1)-EijI(4,2)-EijI(4,3)+EijI(5,2)+
     -   2*(EijI(7,2)-EijI(9,2))+EijI(14,3)+EijI(15,3)-EijI(18,3))
       F(231)=DCMPLX(FR(231),FI(231))
       FR(232) = EijR(2,1)-EijR(4,1)+EijR(4,3)-EijR(5,3)+EijR(7,3)-Eij
     -   R(15,3)-2*(EijR(14,3)-EijR(18,3))
       FI(232) = EijI(2,1)-EijI(4,1)+EijI(4,3)-EijI(5,3)+EijI(7,3)-Eij
     -   I(15,3)-2*(EijI(14,3)-EijI(18,3))
       F(232)=DCMPLX(FR(232),FI(232))
       FR(233) = EijR(2,1)-EijR(2,2)-EijR(4,1)+EijR(4,3)+EijR(5,2)-Eij
     -   R(7,2)-EijR(8,3)+EijR(9,2)+EijR(10,3)-EijR(14,3)-2*(EijR(15,
     -   3)-EijR(18,3))
       FI(233) = EijI(2,1)-EijI(2,2)-EijI(4,1)+EijI(4,3)+EijI(5,2)-Eij
     -   I(7,2)-EijI(8,3)+EijI(9,2)+EijI(10,3)-EijI(14,3)-2*(EijI(15,
     -   3)-EijI(18,3))
       F(233)=DCMPLX(FR(233),FI(233))
       FR(234) = EijR(2,2)+EijR(2,3)+EijR(4,2)-EijR(4,3)-2*EijR(9,2)-3
     -   *(EijR(10,3)-EijR(15,3))
       FI(234) = EijI(2,2)+EijI(2,3)+EijI(4,2)-EijI(4,3)-2*EijI(9,2)-3
     -   *(EijI(10,3)-EijI(15,3))
       F(234)=DCMPLX(FR(234),FI(234))
       FR(235) = EijR(4,2)-EijR(4,3)+EijR(8,2)-EijR(9,2)+EijR(9,3)-Eij
     -   R(10,2)-EijR(10,3)+EijR(16,3)+2*(EijR(15,3)-EijR(20,3))
       FI(235) = EijI(4,2)-EijI(4,3)+EijI(8,2)-EijI(9,2)+EijI(9,3)-Eij
     -   I(10,2)-EijI(10,3)+EijI(16,3)+2*(EijI(15,3)-EijI(20,3))
       F(235)=DCMPLX(FR(235),FI(235))
       FR(236) = EijR(3,1)-EijR(4,1)+EijR(4,3)+EijR(6,2)-EijR(7,2)-Eij
     -   R(8,2)+EijR(9,2)-EijR(14,3)-EijR(15,3)-EijR(16,3)-EijR(17,3)
     -   +EijR(18,3)+EijR(19,3)+EijR(20,3)
       FI(236) = EijI(3,1)-EijI(4,1)+EijI(4,3)+EijI(6,2)-EijI(7,2)-Eij
     -   I(8,2)+EijI(9,2)-EijI(14,3)-EijI(15,3)-EijI(16,3)-EijI(17,3)
     -   +EijI(18,3)+EijI(19,3)+EijI(20,3)
       F(236)=DCMPLX(FR(236),FI(236))
       FR(237) = EijR(3,2)+EijR(4,2)-EijR(4,3)+EijR(12,3)-EijR(13,3)+E
     -   ijR(15,3)-2*(EijR(10,2)-EijR(16,3)+EijR(20,3))
       FI(237) = EijI(3,2)+EijI(4,2)-EijI(4,3)+EijI(12,3)-EijI(13,3)+E
     -   ijI(15,3)-2*(EijI(10,2)-EijI(16,3)+EijI(20,3))
       F(237)=DCMPLX(FR(237),FI(237))
       FR(238) = EijR(2,1)+EijR(2,2)-EijR(4,1)-EijR(4,3)-EijR(9,2)-Eij
     -   R(10,3)+2*EijR(15,3)
       FI(238) = EijI(2,1)+EijI(2,2)-EijI(4,1)-EijI(4,3)-EijI(9,2)-Eij
     -   I(10,3)+2*EijI(15,3)
       F(238)=DCMPLX(FR(238),FI(238))
       FR(239) = EijR(3,1)-EijR(4,1)-EijR(4,3)+EijR(8,2)-EijR(9,2)+Eij
     -   R(15,3)+EijR(16,3)-EijR(20,3)
       FI(239) = EijI(3,1)-EijI(4,1)-EijI(4,3)+EijI(8,2)-EijI(9,2)+Eij
     -   I(15,3)+EijI(16,3)-EijI(20,3)
       F(239)=DCMPLX(FR(239),FI(239))
       FR(240) = EijR(4,2)-EijR(4,3)-EijR(9,2)+EijR(15,3)
       FI(240) = EijI(4,2)-EijI(4,3)-EijI(9,2)+EijI(15,3)
       F(240)=DCMPLX(FR(240),FI(240))
       P(550) = 1-2*Is12*s34
       P(551) = p3sq+s45
       FR(241) = -4*(-(s15*(EE0R+EijR(1,1)))-D01345R*P(92)+Is45*(2*(Ci
     -   j134R(1,1)+p2sq*Dij1234R(2,2))+6*Dij1234R(7,2)+Dij1234R(3,2)
     -   *P(1)-Dij1234R(3,1)*P(11)-(Dij1234R(2,1)+Dij1234R(4,2)-Dij12
     -   34R(5,2))*P(18)-Dij1234R(6,2)*P(157))-(D01234R+Dij1234R(1,1)
     -   )*P(274)+Is12*(Cij145R(1,1)+2*p4sq*Dij1345R(3,2)+6*Dij1345R(
     -   7,2)+Dij1345R(3,1)*P(15)+(-Dij1345R(4,2)+Dij1345R(5,2))*P(95
     -   )-Dij1345R(2,1)*P(106)+Dij1345R(2,2)*P(287)-Dij1345R(6,2)*P(
     -   293))-2*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)-2*p3sq*EijR(3,2
     -   )-9*EijR(11,2)+EijR(22,3)-EijR(24,3)+EijR(5,2)*P(18)-EijR(7,
     -   2)*P(36)+EijR(6,2)*P(46)-EijR(2,2)*P(68)-EijR(4,2)*P(174)+Ei
     -   jR(10,2)*P(407)+EijR(9,2)*P(411)+EijR(8,2)*P(413))+Dij1345R(
     -   1,1)*P(550))+Is12s45*(-6*B013R-2*B014R+8*(Cij123R(4,2)+Cij13
     -   4R(4,2)+Cij145R(4,2))-4*(p2sq*Cij123R(2,1)-p4sq*Cij145R(2,1)
     -   -C0145R*P(15)-(C0123R+Cij123R(1,1))*P(18)-Cij134R(2,1)*P(25)
     -   -C0134R*P(551)))
       FI(241) = -4*(-(s15*(EE0I+EijI(1,1)))-D01345I*P(92)+Is45*(2*(Ci
     -   j134I(1,1)+p2sq*Dij1234I(2,2))+6*Dij1234I(7,2)+Dij1234I(3,2)
     -   *P(1)-Dij1234I(3,1)*P(11)-(Dij1234I(2,1)+Dij1234I(4,2)-Dij12
     -   34I(5,2))*P(18)-Dij1234I(6,2)*P(157))-(D01234I+Dij1234I(1,1)
     -   )*P(274)+Is12*(Cij145I(1,1)+2*p4sq*Dij1345I(3,2)+6*Dij1345I(
     -   7,2)+Dij1345I(3,1)*P(15)+(-Dij1345I(4,2)+Dij1345I(5,2))*P(95
     -   )-Dij1345I(2,1)*P(106)+Dij1345I(2,2)*P(287)-Dij1345I(6,2)*P(
     -   293))-2*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)-2*p3sq*EijI(3,2
     -   )-9*EijI(11,2)+EijI(22,3)-EijI(24,3)+EijI(5,2)*P(18)-EijI(7,
     -   2)*P(36)+EijI(6,2)*P(46)-EijI(2,2)*P(68)-EijI(4,2)*P(174)+Ei
     -   jI(10,2)*P(407)+EijI(9,2)*P(411)+EijI(8,2)*P(413))+Dij1345I(
     -   1,1)*P(550))+Is12s45*(-6*B013I-2*B014I+8*(Cij123I(4,2)+Cij13
     -   4I(4,2)+Cij145I(4,2))-4*(p2sq*Cij123I(2,1)-p4sq*Cij145I(2,1)
     -   -C0145I*P(15)-(C0123I+Cij123I(1,1))*P(18)-Cij134I(2,1)*P(25)
     -   -C0134I*P(551)))
       F(241)=DCMPLX(FR(241),FI(241))
       FR(242) = 8*(EijR(3,1)-EijR(4,1)+EijR(6,2)-EijR(9,2)+2*(EijR(7,
     -   2)-EijR(10,2))-EijR(15,3)+EijR(16,3)+EijR(18,3)-EijR(19,3))
       FI(242) = 8*(EijI(3,1)-EijI(4,1)+EijI(6,2)-EijI(9,2)+2*(EijI(7,
     -   2)-EijI(10,2))-EijI(15,3)+EijI(16,3)+EijI(18,3)-EijI(19,3))
       F(242)=DCMPLX(FR(242),FI(242))
       FR(243) = 8*(EijR(3,1)-EijR(4,1)+EijR(8,2)+EijR(9,2)-2*EijR(10,
     -   2)+EijR(10,3)-EijR(15,3)+EijR(16,3)-EijR(20,3))
       FI(243) = 8*(EijI(3,1)-EijI(4,1)+EijI(8,2)+EijI(9,2)-2*EijI(10,
     -   2)+EijI(10,3)-EijI(15,3)+EijI(16,3)-EijI(20,3))
       F(243)=DCMPLX(FR(243),FI(243))
       FR(244) = -8*(EijR(2,1)-EijR(4,1)+EijR(5,2)+2*(EijR(7,2)-EijR(9
     -   ,2))-EijR(10,2)+EijR(15,3)-EijR(16,3)-EijR(18,3)+EijR(19,3))
       FI(244) = -8*(EijI(2,1)-EijI(4,1)+EijI(5,2)+2*(EijI(7,2)-EijI(9
     -   ,2))-EijI(10,2)+EijI(15,3)-EijI(16,3)-EijI(18,3)+EijI(19,3))
       F(244)=DCMPLX(FR(244),FI(244))
       FR(245) = 2*Is12s45*(3*B013R+B014R+2*(p2sq*Cij123R(2,1)-p4sq*Ci
     -   j145R(2,1)-2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2))-C0145R
     -   *P(15)-(C0123R+Cij123R(1,1))*P(18)-Cij134R(2,1)*P(30)-C0134R
     -   *P(43)))+4*(Is45*(Cij134R(1,1)+2*(Dij1234R(7,2)+Dij1234R(12,
     -   3)-Dij1234R(13,3)))-s15*(EE0R+EijR(1,1))+D01345R*P(39)+(D012
     -   34R+Dij1234R(1,1))*P(40)+Is12*(Cij145R(1,1)+2*p4sq*Dij1345R(
     -   3,2)+6*Dij1345R(7,2)+Dij1345R(3,1)*P(15)+(Dij1345R(4,2)-Dij1
     -   345R(5,2))*P(261)+Dij1345R(2,1)*P(284)+Dij1345R(2,2)*P(287)-
     -   Dij1345R(6,2)*P(293))-2*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)
     -   -2*p3sq*EijR(3,2)-9*EijR(11,2)+EijR(22,3)-EijR(23,3)+EijR(5,
     -   2)*P(18)-EijR(7,2)*P(36)-EijR(6,2)*P(45)-EijR(2,2)*P(68)-Eij
     -   R(4,2)*P(174)+EijR(10,2)*P(407)+EijR(9,2)*P(411)+EijR(8,2)*P
     -   (413))+Dij1345R(1,1)*P(550))
       FI(245) = 2*Is12s45*(3*B013I+B014I+2*(p2sq*Cij123I(2,1)-p4sq*Ci
     -   j145I(2,1)-2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2))-C0145I
     -   *P(15)-(C0123I+Cij123I(1,1))*P(18)-Cij134I(2,1)*P(30)-C0134I
     -   *P(43)))+4*(Is45*(Cij134I(1,1)+2*(Dij1234I(7,2)+Dij1234I(12,
     -   3)-Dij1234I(13,3)))-s15*(EE0I+EijI(1,1))+D01345I*P(39)+(D012
     -   34I+Dij1234I(1,1))*P(40)+Is12*(Cij145I(1,1)+2*p4sq*Dij1345I(
     -   3,2)+6*Dij1345I(7,2)+Dij1345I(3,1)*P(15)+(Dij1345I(4,2)-Dij1
     -   345I(5,2))*P(261)+Dij1345I(2,1)*P(284)+Dij1345I(2,2)*P(287)-
     -   Dij1345I(6,2)*P(293))-2*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)
     -   -2*p3sq*EijI(3,2)-9*EijI(11,2)+EijI(22,3)-EijI(23,3)+EijI(5,
     -   2)*P(18)-EijI(7,2)*P(36)-EijI(6,2)*P(45)-EijI(2,2)*P(68)-Eij
     -   I(4,2)*P(174)+EijI(10,2)*P(407)+EijI(9,2)*P(411)+EijI(8,2)*P
     -   (413))+Dij1345I(1,1)*P(550))
       F(245)=DCMPLX(FR(245),FI(245))
       FR(246) = 8*Is45*(2*(Dij1234R(2,1)-Dij1234R(3,1))+Dij1234R(4,3)
     -   +3*(Dij1234R(4,2)-Dij1234R(5,2))-Dij1234R(5,3)+s45*(EijR(2,1
     -   )-EijR(3,1)-EijR(5,3)+EijR(6,3)-EijR(15,3)+EijR(16,3)+2*(Eij
     -   R(18,3)-EijR(19,3))))
       FI(246) = 8*Is45*(2*(Dij1234I(2,1)-Dij1234I(3,1))+Dij1234I(4,3)
     -   +3*(Dij1234I(4,2)-Dij1234I(5,2))-Dij1234I(5,3)+s45*(EijI(2,1
     -   )-EijI(3,1)-EijI(5,3)+EijI(6,3)-EijI(15,3)+EijI(16,3)+2*(Eij
     -   I(18,3)-EijI(19,3))))
       F(246)=DCMPLX(FR(246),FI(246))
       FR(247) = 8*Is45*(Dij1234R(2,2)-Dij1234R(5,2)+2*(Dij1234R(2,1)-
     -   Dij1234R(3,1)+Dij1234R(4,2)-Dij1234R(6,2))+Dij1234R(6,3)-Dij
     -   1234R(10,3)+s45*(EijR(2,1)-EijR(2,2)-EijR(3,1)+EijR(5,2)+Eij
     -   R(6,2)-EijR(8,2)-EijR(8,3)-2*(EijR(7,2)-EijR(9,2))+EijR(10,3
     -   )-EijR(15,3)+EijR(16,3)+EijR(17,3)+EijR(18,3)-EijR(19,3)-Eij
     -   R(20,3)))
       FI(247) = 8*Is45*(Dij1234I(2,2)-Dij1234I(5,2)+2*(Dij1234I(2,1)-
     -   Dij1234I(3,1)+Dij1234I(4,2)-Dij1234I(6,2))+Dij1234I(6,3)-Dij
     -   1234I(10,3)+s45*(EijI(2,1)-EijI(2,2)-EijI(3,1)+EijI(5,2)+Eij
     -   I(6,2)-EijI(8,2)-EijI(8,3)-2*(EijI(7,2)-EijI(9,2))+EijI(10,3
     -   )-EijI(15,3)+EijI(16,3)+EijI(17,3)+EijI(18,3)-EijI(19,3)-Eij
     -   I(20,3)))
       F(247)=DCMPLX(FR(247),FI(247))
       FR(248) = 8*Is45*(Dij1234R(2,1)+2*Dij1234R(2,2)-Dij1234R(3,1)-D
     -   ij1234R(5,2)-Dij1234R(6,2)+Dij1234R(6,3)-Dij1234R(10,3)+s45*
     -   (EijR(2,2)-2*(EijR(5,2)-EijR(7,2))+EijR(8,2)-EijR(8,3)-EijR(
     -   9,2)-EijR(10,2)+EijR(10,3)-EijR(15,3)+EijR(16,3)+EijR(17,3)+
     -   EijR(18,3)-EijR(19,3)-EijR(20,3)))
       FI(248) = 8*Is45*(Dij1234I(2,1)+2*Dij1234I(2,2)-Dij1234I(3,1)-D
     -   ij1234I(5,2)-Dij1234I(6,2)+Dij1234I(6,3)-Dij1234I(10,3)+s45*
     -   (EijI(2,2)-2*(EijI(5,2)-EijI(7,2))+EijI(8,2)-EijI(8,3)-EijI(
     -   9,2)-EijI(10,2)+EijI(10,3)-EijI(15,3)+EijI(16,3)+EijI(17,3)+
     -   EijI(18,3)-EijI(19,3)-EijI(20,3)))
       F(248)=DCMPLX(FR(248),FI(248))
       FR(249) = 8*Is45*(Dij1234R(2,1)+Dij1234R(2,3)-Dij1234R(3,1)+2*(
     -   Dij1234R(2,2)-Dij1234R(6,2))-Dij1234R(8,3)+s45*(-EijR(2,2)-E
     -   ijR(2,3)+EijR(8,2)+EijR(9,2)+EijR(9,3)-EijR(10,2)-EijR(15,3)
     -   +EijR(16,3)+2*(EijR(10,3)-EijR(20,3))))
       FI(249) = 8*Is45*(Dij1234I(2,1)+Dij1234I(2,3)-Dij1234I(3,1)+2*(
     -   Dij1234I(2,2)-Dij1234I(6,2))-Dij1234I(8,3)+s45*(-EijI(2,2)-E
     -   ijI(2,3)+EijI(8,2)+EijI(9,2)+EijI(9,3)-EijI(10,2)-EijI(15,3)
     -   +EijI(16,3)+2*(EijI(10,3)-EijI(20,3))))
       F(249)=DCMPLX(FR(249),FI(249))
       FR(250) = -8*Is45*(Dij1234R(3,2)+Dij1234R(5,2)-2*Dij1234R(6,2)+
     -   Dij1234R(7,3)-Dij1234R(10,3)+s45*(-EijR(3,2)+2*(EijR(6,2)-Ei
     -   jR(7,2))-EijR(8,2)+EijR(9,2)+EijR(10,2)-EijR(11,3)+EijR(13,3
     -   )+EijR(15,3)-EijR(16,3)+EijR(17,3)-EijR(18,3)+EijR(19,3)-Eij
     -   R(20,3)))
       FI(250) = -8*Is45*(Dij1234I(3,2)+Dij1234I(5,2)-2*Dij1234I(6,2)+
     -   Dij1234I(7,3)-Dij1234I(10,3)+s45*(-EijI(3,2)+2*(EijI(6,2)-Ei
     -   jI(7,2))-EijI(8,2)+EijI(9,2)+EijI(10,2)-EijI(11,3)+EijI(13,3
     -   )+EijI(15,3)-EijI(16,3)+EijI(17,3)-EijI(18,3)+EijI(19,3)-Eij
     -   I(20,3)))
       F(250)=DCMPLX(FR(250),FI(250))
       FR(251) = 8*Is45*(Dij1234R(2,1)-Dij1234R(3,1)-2*Dij1234R(3,2)+D
     -   ij1234R(4,2)+Dij1234R(6,2)-Dij1234R(7,3)+Dij1234R(10,3)+s45*
     -   (-EijR(3,2)+2*(EijR(6,2)-EijR(7,2))-EijR(8,2)+EijR(9,2)+EijR
     -   (10,2)+EijR(11,3)-EijR(13,3)-EijR(15,3)+EijR(16,3)-EijR(17,3
     -   )+EijR(18,3)-EijR(19,3)+EijR(20,3)))
       FI(251) = 8*Is45*(Dij1234I(2,1)-Dij1234I(3,1)-2*Dij1234I(3,2)+D
     -   ij1234I(4,2)+Dij1234I(6,2)-Dij1234I(7,3)+Dij1234I(10,3)+s45*
     -   (-EijI(3,2)+2*(EijI(6,2)-EijI(7,2))-EijI(8,2)+EijI(9,2)+EijI
     -   (10,2)+EijI(11,3)-EijI(13,3)-EijI(15,3)+EijI(16,3)-EijI(17,3
     -   )+EijI(18,3)-EijI(19,3)+EijI(20,3)))
       F(251)=DCMPLX(FR(251),FI(251))
       FR(252) = EijR(3,1)+EijR(3,2)-EijR(4,1)-EijR(10,2)-EijR(13,3)-E
     -   ijR(15,3)+EijR(16,3)+EijR(20,3)
       FI(252) = EijI(3,1)+EijI(3,2)-EijI(4,1)-EijI(10,2)-EijI(13,3)-E
     -   ijI(15,3)+EijI(16,3)+EijI(20,3)
       F(252)=DCMPLX(FR(252),FI(252))
       FR(253) = EijR(2,1)+EijR(2,2)-EijR(4,1)-EijR(9,2)-EijR(10,3)+Ei
     -   jR(15,3)-EijR(16,3)+EijR(20,3)
       FI(253) = EijI(2,1)+EijI(2,2)-EijI(4,1)-EijI(9,2)-EijI(10,3)+Ei
     -   jI(15,3)-EijI(16,3)+EijI(20,3)
       F(253)=DCMPLX(FR(253),FI(253))
       FR(254) = EijR(3,1)-EijR(4,1)+EijR(8,2)-EijR(9,2)+EijR(13,3)+Ei
     -   jR(15,3)-EijR(16,3)-EijR(20,3)
       FI(254) = EijI(3,1)-EijI(4,1)+EijI(8,2)-EijI(9,2)+EijI(13,3)+Ei
     -   jI(15,3)-EijI(16,3)-EijI(20,3)
       F(254)=DCMPLX(FR(254),FI(254))
       FR(255) = -EijR(9,2)+EijR(10,2)+EijR(15,3)-EijR(16,3)
       FI(255) = -EijI(9,2)+EijI(10,2)+EijI(15,3)-EijI(16,3)
       F(255)=DCMPLX(FR(255),FI(255))
       FR(256) = Is45*(-Dij1234R(3,2)+Dij1234R(6,2)+Dij1234R(8,3)-Dij1
     -   234R(9,3)+s45*(EijR(3,2)-EijR(8,2)+EijR(9,2)-EijR(9,3)-EijR(
     -   10,2)+EijR(10,3)+EijR(12,3)-EijR(13,3)-EijR(15,3)+EijR(16,3)
     -   ))
       FI(256) = Is45*(-Dij1234I(3,2)+Dij1234I(6,2)+Dij1234I(8,3)-Dij1
     -   234I(9,3)+s45*(EijI(3,2)-EijI(8,2)+EijI(9,2)-EijI(9,3)-EijI(
     -   10,2)+EijI(10,3)+EijI(12,3)-EijI(13,3)-EijI(15,3)+EijI(16,3)
     -   ))
       F(256)=DCMPLX(FR(256),FI(256))
       FR(257) = Is45*(-Dij1234R(2,1)-Dij1234R(2,2)+Dij1234R(3,1)+Dij1
     -   234R(3,2)-Dij1234R(8,3)+Dij1234R(9,3)+s45*(EijR(9,3)-EijR(10
     -   ,3)-EijR(12,3)+EijR(13,3)+EijR(15,3)-EijR(16,3)))
       FI(257) = Is45*(-Dij1234I(2,1)-Dij1234I(2,2)+Dij1234I(3,1)+Dij1
     -   234I(3,2)-Dij1234I(8,3)+Dij1234I(9,3)+s45*(EijI(9,3)-EijI(10
     -   ,3)-EijI(12,3)+EijI(13,3)+EijI(15,3)-EijI(16,3)))
       F(257)=DCMPLX(FR(257),FI(257))
       FR(258) = Is45*(Dij1234R(3,2)+Dij1234R(3,3)-Dij1234R(6,2)-Dij12
     -   34R(9,3)+s45*(-EijR(3,3)+EijR(12,3)+EijR(15,3)-EijR(16,3)+2*
     -   (EijR(13,3)-EijR(20,3))))
       FI(258) = Is45*(Dij1234I(3,2)+Dij1234I(3,3)-Dij1234I(6,2)-Dij12
     -   34I(9,3)+s45*(-EijI(3,3)+EijI(12,3)+EijI(15,3)-EijI(16,3)+2*
     -   (EijI(13,3)-EijI(20,3))))
       F(258)=DCMPLX(FR(258),FI(258))
       P(552) = -p2sq+s15+2*P(112)
       P(553) = 2*s12-s34+3*P(30)
       P(554) = s12-2*s34+3*P(30)
       FR(259) = 4*(-Dij1345R(2,1)+Dij1345R(3,1)-Dij1345R(3,2)-Dij1345
     -   R(4,2)+Dij1345R(5,2)+Dij1345R(6,2)+Dij2345R(2,1)+Dij2345R(2,
     -   2)-Dij2345R(3,1)+Dij2345R(3,2)+s12*EijR(8,2)-s12*EijR(9,2)-2
     -   *(Dij2345R(6,2)+EijR(23,3)-EijR(24,3))+(EijR(3,1)+EijR(6,2)-
     -   EijR(11,3)-EijR(14,3)+2*EijR(19,3))*P(18)-EijR(3,3)*P(25)-Ei
     -   jR(3,2)*P(30)+(EijR(12,3)+EijR(15,3)-2*EijR(20,3))*P(68)+(Ei
     -   jR(4,1)+EijR(7,2))*P(160)-EijR(4,3)*P(261)+EijR(4,2)*P(442)+
     -   EijR(10,2)*P(552)+EijR(13,3)*P(553)-EijR(16,3)*P(554))
       FI(259) = 4*(-Dij1345I(2,1)+Dij1345I(3,1)-Dij1345I(3,2)-Dij1345
     -   I(4,2)+Dij1345I(5,2)+Dij1345I(6,2)+Dij2345I(2,1)+Dij2345I(2,
     -   2)-Dij2345I(3,1)+Dij2345I(3,2)+s12*EijI(8,2)-s12*EijI(9,2)-2
     -   *(Dij2345I(6,2)+EijI(23,3)-EijI(24,3))+(EijI(3,1)+EijI(6,2)-
     -   EijI(11,3)-EijI(14,3)+2*EijI(19,3))*P(18)-EijI(3,3)*P(25)-Ei
     -   jI(3,2)*P(30)+(EijI(12,3)+EijI(15,3)-2*EijI(20,3))*P(68)+(Ei
     -   jI(4,1)+EijI(7,2))*P(160)-EijI(4,3)*P(261)+EijI(4,2)*P(442)+
     -   EijI(10,2)*P(552)+EijI(13,3)*P(553)-EijI(16,3)*P(554))
       F(259)=DCMPLX(FR(259),FI(259))
       P(555) = p2sq+2*s12
       P(556) = p2sq-2*P(25)
       P(557) = p2sq-p3sq+s45
       P(558) = -s15-3*P(230)+2*P(557)
       P(559) = p2sq-s34+2*P(30)
       FR(260) = -4*(Dij1345R(1,1)+Dij1345R(1,2)-Dij1345R(3,1)+Dij1345
     -   R(3,2)-Dij2345R(3,2)+Dij2345R(6,2)-p2sq*EijR(2,1)+s12*(EijR(
     -   3,1)+EijR(9,2))-s34*EijR(10,2)-6*EijR(11,2)+2*(D02345R-Dij13
     -   45R(5,2)+Dij2345R(2,1)+EijR(21,3)-EijR(23,3))+(EijR(4,1)+Eij
     -   R(6,3)-EijR(7,3))*P(18)-(EijR(3,2)-EijR(11,3)+EijR(13,3))*P(
     -   25)+(EijR(8,2)-EijR(15,3)-EijR(17,3)+EijR(18,3)+EijR(20,3))*
     -   P(68)+EijR(4,3)*P(261)-EijR(4,2)*P(442)+EijR(14,3)*P(496)+Ei
     -   jR(16,3)*P(509)-EijR(5,2)*P(555)-EijR(6,2)*P(556)+EijR(7,2)*
     -   P(558)-EijR(19,3)*P(559))
       FI(260) = -4*(Dij1345I(1,1)+Dij1345I(1,2)-Dij1345I(3,1)+Dij1345
     -   I(3,2)-Dij2345I(3,2)+Dij2345I(6,2)-p2sq*EijI(2,1)+s12*(EijI(
     -   3,1)+EijI(9,2))-s34*EijI(10,2)-6*EijI(11,2)+2*(D02345I-Dij13
     -   45I(5,2)+Dij2345I(2,1)+EijI(21,3)-EijI(23,3))+(EijI(4,1)+Eij
     -   I(6,3)-EijI(7,3))*P(18)-(EijI(3,2)-EijI(11,3)+EijI(13,3))*P(
     -   25)+(EijI(8,2)-EijI(15,3)-EijI(17,3)+EijI(18,3)+EijI(20,3))*
     -   P(68)+EijI(4,3)*P(261)-EijI(4,2)*P(442)+EijI(14,3)*P(496)+Ei
     -   jI(16,3)*P(509)-EijI(5,2)*P(555)-EijI(6,2)*P(556)+EijI(7,2)*
     -   P(558)-EijI(19,3)*P(559))
       F(260)=DCMPLX(FR(260),FI(260))
       P(560) = p2sq+s12-s23
       P(561) = 3*p3sq-s45+2*P(560)
       P(562) = 2*p2sq+s12
       P(563) = -p2sq+p3sq+s15+s45+2*P(498)
       P(564) = p2sq-s15-2*s23-s45+3*P(42)
       P(565) = p4sq-s12-s34
       P(566) = 5*p3sq+s45+2*P(565)
       FR(261) = Is12s45*(-6*B013R-2*B014R+8*(Cij123R(4,2)+Cij134R(4,2
     -   )+Cij145R(4,2))-4*(p2sq*Cij123R(2,1)-p4sq*Cij145R(2,1)-C0145
     -   R*P(15)-(C0123R+Cij123R(1,1))*P(18)-Cij134R(2,1)*P(25)-C0134
     -   R*P(551)))-4*(Dij1345R(1,2)+Dij2345R(1,1)-Dij2345R(3,1)-Dij2
     -   345R(3,2)-Dij2345R(4,2)+Dij2345R(5,2)+Dij2345R(6,2)-s15*(EE0
     -   R+EijR(1,1))+4*p3sq*EijR(3,2)+18*EijR(11,2)-(EijR(2,1)+3*Eij
     -   R(5,2)-EijR(14,3)-EijR(17,3)+EijR(18,3)+EijR(19,3))*P(18)-(-
     -   EijR(12,3)+EijR(13,3))*P(25)-2*(D02345R-EijR(22,3)+EijR(23,3
     -   )+EijR(6,2)*P(46))-(EijR(9,3)-EijR(10,3))*P(68)+EijR(7,2)*P(
     -   73)-(D01345R+2*Dij1345R(1,1))*P(92)+Is45*(2*(Cij134R(1,1)+p2
     -   sq*Dij1234R(2,2))+6*Dij1234R(7,2)+Dij1234R(3,2)*P(1)-Dij1234
     -   R(3,1)*P(11)-(Dij1234R(2,1)+Dij1234R(4,2)-Dij1234R(5,2))*P(1
     -   8)-Dij1234R(6,2)*P(157))-EijR(4,1)*P(160)+EijR(4,3)*P(261)-(
     -   D01234R+Dij1234R(1,1))*P(274)+Is12*(Cij145R(1,1)+6*Dij1345R(
     -   7,2)-Dij1345R(4,2)*P(95)-Dij1345R(5,2)*P(96)-Dij1345R(2,1)*P
     -   (106)+Dij1345R(2,2)*P(287)-Dij1345R(6,2)*P(293)-Dij1345R(3,1
     -   )*P(431)+Dij1345R(3,2)*P(432))-EijR(15,3)*P(497)+EijR(16,3)*
     -   P(509)+EijR(20,3)*P(510)-EijR(8,2)*P(561)+EijR(2,2)*P(562)+E
     -   ijR(4,2)*P(563)+EijR(9,2)*P(564)-EijR(10,2)*P(566))
       FI(261) = Is12s45*(-6*B013I-2*B014I+8*(Cij123I(4,2)+Cij134I(4,2
     -   )+Cij145I(4,2))-4*(p2sq*Cij123I(2,1)-p4sq*Cij145I(2,1)-C0145
     -   I*P(15)-(C0123I+Cij123I(1,1))*P(18)-Cij134I(2,1)*P(25)-C0134
     -   I*P(551)))-4*(Dij1345I(1,2)+Dij2345I(1,1)-Dij2345I(3,1)-Dij2
     -   345I(3,2)-Dij2345I(4,2)+Dij2345I(5,2)+Dij2345I(6,2)-s15*(EE0
     -   I+EijI(1,1))+4*p3sq*EijI(3,2)+18*EijI(11,2)-(EijI(2,1)+3*Eij
     -   I(5,2)-EijI(14,3)-EijI(17,3)+EijI(18,3)+EijI(19,3))*P(18)-(-
     -   EijI(12,3)+EijI(13,3))*P(25)-2*(D02345I-EijI(22,3)+EijI(23,3
     -   )+EijI(6,2)*P(46))-(EijI(9,3)-EijI(10,3))*P(68)+EijI(7,2)*P(
     -   73)-(D01345I+2*Dij1345I(1,1))*P(92)+Is45*(2*(Cij134I(1,1)+p2
     -   sq*Dij1234I(2,2))+6*Dij1234I(7,2)+Dij1234I(3,2)*P(1)-Dij1234
     -   I(3,1)*P(11)-(Dij1234I(2,1)+Dij1234I(4,2)-Dij1234I(5,2))*P(1
     -   8)-Dij1234I(6,2)*P(157))-EijI(4,1)*P(160)+EijI(4,3)*P(261)-(
     -   D01234I+Dij1234I(1,1))*P(274)+Is12*(Cij145I(1,1)+6*Dij1345I(
     -   7,2)-Dij1345I(4,2)*P(95)-Dij1345I(5,2)*P(96)-Dij1345I(2,1)*P
     -   (106)+Dij1345I(2,2)*P(287)-Dij1345I(6,2)*P(293)-Dij1345I(3,1
     -   )*P(431)+Dij1345I(3,2)*P(432))-EijI(15,3)*P(497)+EijI(16,3)*
     -   P(509)+EijI(20,3)*P(510)-EijI(8,2)*P(561)+EijI(2,2)*P(562)+E
     -   ijI(4,2)*P(563)+EijI(9,2)*P(564)-EijI(10,2)*P(566))
       F(261)=DCMPLX(FR(261),FI(261))
       FR(262) = -4*(D01345R-D02345R+Dij1345R(1,1)+Dij1345R(3,2)-Dij13
     -   45R(5,2)-Dij2345R(2,1)-Dij2345R(3,2)+Dij2345R(6,2)-s12*EijR(
     -   2,1)-p2sq*EijR(9,2)+2*(EijR(11,2)-EijR(24,3))+(EijR(6,2)+2*E
     -   ijR(7,2)+EijR(14,3)-EijR(19,3))*P(18)+(EijR(3,2)-EijR(13,3))
     -   *P(25)+EijR(3,1)*P(30)-(EijR(8,2)+EijR(15,3)-EijR(20,3))*P(6
     -   8)+EijR(4,1)*P(144)-EijR(4,2)*P(222)+EijR(4,3)*P(261)+EijR(1
     -   0,2)*P(436)+EijR(16,3)*P(509))
       FI(262) = -4*(D01345I-D02345I+Dij1345I(1,1)+Dij1345I(3,2)-Dij13
     -   45I(5,2)-Dij2345I(2,1)-Dij2345I(3,2)+Dij2345I(6,2)-s12*EijI(
     -   2,1)-p2sq*EijI(9,2)+2*(EijI(11,2)-EijI(24,3))+(EijI(6,2)+2*E
     -   ijI(7,2)+EijI(14,3)-EijI(19,3))*P(18)+(EijI(3,2)-EijI(13,3))
     -   *P(25)+EijI(3,1)*P(30)-(EijI(8,2)+EijI(15,3)-EijI(20,3))*P(6
     -   8)+EijI(4,1)*P(144)-EijI(4,2)*P(222)+EijI(4,3)*P(261)+EijI(1
     -   0,2)*P(436)+EijI(16,3)*P(509))
       F(262)=DCMPLX(FR(262),FI(262))
       P(567) = 3*p3sq-s12+s45
       P(568) = s15-2*s23-s34-s45+3*P(147)
       P(569) = p2sq-3*p3sq+s34+s45+2*P(177)
       P(570) = p2sq+2*p4sq-5*s12-s15+3*P(551)
       FR(263) = 2*Is12s45*(3*B013R+B014R+2*(p2sq*Cij123R(2,1)-p4sq*Ci
     -   j145R(2,1)-2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2))-C0145R
     -   *P(15)-(C0123R+Cij123R(1,1))*P(18)-C0134R*P(43)+Cij134R(2,1)
     -   *P(61)))-4*(-Dij1345R(1,2)-Dij2345R(1,1)-Dij2345R(2,1)+s15*(
     -   EE0R+EijR(1,1))+EijR(2,1)*P(18)+Is45*(Cij134R(1,2)-Cij134R(3
     -   ,2)-p2sq*Dij1234R(2,2)+s12*(-Dij1234R(2,1)+Dij1234R(3,1)+Dij
     -   1234R(6,2))+2*Dij1234R(12,3)+Dij1234R(4,2)*P(18))-EijR(4,1)*
     -   P(28)-(D01345R+2*Dij1345R(1,1))*P(39)-(D01234R+Dij1234R(1,1)
     -   )*P(40)-EijR(6,2)*P(72)-3*EijR(5,2)*P(160)-EijR(3,1)*P(222)+
     -   2*(D02345R+Dij2345R(3,1)-9*EijR(11,2)-EijR(22,3)+EijR(24,3)-
     -   EijR(4,2)*P(174)-EijR(7,2)*P(434))-Is12*(Cij145R(1,1)+2*p4sq
     -   *Dij1345R(3,2)+6*Dij1345R(7,2)+Dij1345R(3,1)*P(15)-Dij1345R(
     -   4,2)*P(112)-Dij1345R(2,1)*P(235)+Dij1345R(2,2)*P(287)-Dij134
     -   5R(5,2)*P(436)+Dij1345R(6,2)*P(437))-EijR(2,2)*P(562)-EijR(3
     -   ,2)*P(567)+EijR(8,2)*P(568)+EijR(9,2)*P(569)+EijR(10,2)*P(57
     -   0))
       FI(263) = 2*Is12s45*(3*B013I+B014I+2*(p2sq*Cij123I(2,1)-p4sq*Ci
     -   j145I(2,1)-2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2))-C0145I
     -   *P(15)-(C0123I+Cij123I(1,1))*P(18)-C0134I*P(43)+Cij134I(2,1)
     -   *P(61)))-4*(-Dij1345I(1,2)-Dij2345I(1,1)-Dij2345I(2,1)+s15*(
     -   EE0I+EijI(1,1))+EijI(2,1)*P(18)+Is45*(Cij134I(1,2)-Cij134I(3
     -   ,2)-p2sq*Dij1234I(2,2)+s12*(-Dij1234I(2,1)+Dij1234I(3,1)+Dij
     -   1234I(6,2))+2*Dij1234I(12,3)+Dij1234I(4,2)*P(18))-EijI(4,1)*
     -   P(28)-(D01345I+2*Dij1345I(1,1))*P(39)-(D01234I+Dij1234I(1,1)
     -   )*P(40)-EijI(6,2)*P(72)-3*EijI(5,2)*P(160)-EijI(3,1)*P(222)+
     -   2*(D02345I+Dij2345I(3,1)-9*EijI(11,2)-EijI(22,3)+EijI(24,3)-
     -   EijI(4,2)*P(174)-EijI(7,2)*P(434))-Is12*(Cij145I(1,1)+2*p4sq
     -   *Dij1345I(3,2)+6*Dij1345I(7,2)+Dij1345I(3,1)*P(15)-Dij1345I(
     -   4,2)*P(112)-Dij1345I(2,1)*P(235)+Dij1345I(2,2)*P(287)-Dij134
     -   5I(5,2)*P(436)+Dij1345I(6,2)*P(437))-EijI(2,2)*P(562)-EijI(3
     -   ,2)*P(567)+EijI(8,2)*P(568)+EijI(9,2)*P(569)+EijI(10,2)*P(57
     -   0))
       F(263)=DCMPLX(FR(263),FI(263))
       P(571) = p2sq-3*s12-s15+s34
       FR(264) = -4*(Dij1345R(2,2)-Dij1345R(4,2)+Dij1345R(5,2)-Dij1345
     -   R(6,2)-p2sq*EijR(8,2)+p2sq*EijR(9,2)+Is45*(Cij134R(1,1)-Cij1
     -   34R(2,1)-Cij134R(2,2)+Cij134R(3,2)+s12*Dij1234R(3,2)-p2sq*(D
     -   ij1234R(2,1)-Dij1234R(3,1)+Dij1234R(6,2))+2*(Dij1234R(7,2)+D
     -   ij1234R(13,3))+Dij1234R(5,2)*P(18))-2*(EijR(23,3)-EijR(24,3)
     -   +(-EijR(6,2)+EijR(7,2))*P(18))-(-EijR(3,1)+EijR(4,1))*P(28)-
     -   (EijR(3,2)-EijR(10,2))*P(571))
       FI(264) = -4*(Dij1345I(2,2)-Dij1345I(4,2)+Dij1345I(5,2)-Dij1345
     -   I(6,2)-p2sq*EijI(8,2)+p2sq*EijI(9,2)+Is45*(Cij134I(1,1)-Cij1
     -   34I(2,1)-Cij134I(2,2)+Cij134I(3,2)+s12*Dij1234I(3,2)-p2sq*(D
     -   ij1234I(2,1)-Dij1234I(3,1)+Dij1234I(6,2))+2*(Dij1234I(7,2)+D
     -   ij1234I(13,3))+Dij1234I(5,2)*P(18))-2*(EijI(23,3)-EijI(24,3)
     -   +(-EijI(6,2)+EijI(7,2))*P(18))-(-EijI(3,1)+EijI(4,1))*P(28)-
     -   (EijI(3,2)-EijI(10,2))*P(571))
       F(264)=DCMPLX(FR(264),FI(264))
       P(572) = p2sq-s12+s15-s34-2*P(30)
       P(573) = p2sq-s15-2*P(30)-3*P(230)
       FR(265) = -4*(-Dij1345R(1,1)-Dij1345R(1,2)+Dij1345R(2,1)+Dij134
     -   5R(4,2)+Dij1345R(5,2)-Dij1345R(6,2)+p2sq*EijR(2,1)+p2sq*EijR
     -   (9,2)+6*EijR(11,2)+Is45*(Cij134R(1,2)-Cij134R(3,2)+p2sq*Dij1
     -   234R(4,2)+4*Dij1234R(7,2)+2*(Cij134R(1,1)-Cij134R(2,1)-p2sq*
     -   Dij1234R(6,2)+Dij1234R(11,3))+Dij1234R(3,2)*P(1)-Dij1234R(5,
     -   2)*P(6))-EijR(3,1)*P(8)-EijR(4,1)*P(28)-2*(D02345R+Dij2345R(
     -   2,1)+EijR(21,3)-EijR(24,3)-EijR(3,2)*P(25)+EijR(8,2)*P(68)+E
     -   ijR(7,2)*P(441))+EijR(5,2)*P(555)+EijR(6,2)*P(572)+EijR(10,2
     -   )*P(573))
       FI(265) = -4*(-Dij1345I(1,1)-Dij1345I(1,2)+Dij1345I(2,1)+Dij134
     -   5I(4,2)+Dij1345I(5,2)-Dij1345I(6,2)+p2sq*EijI(2,1)+p2sq*EijI
     -   (9,2)+6*EijI(11,2)+Is45*(Cij134I(1,2)-Cij134I(3,2)+p2sq*Dij1
     -   234I(4,2)+4*Dij1234I(7,2)+2*(Cij134I(1,1)-Cij134I(2,1)-p2sq*
     -   Dij1234I(6,2)+Dij1234I(11,3))+Dij1234I(3,2)*P(1)-Dij1234I(5,
     -   2)*P(6))-EijI(3,1)*P(8)-EijI(4,1)*P(28)-2*(D02345I+Dij2345I(
     -   2,1)+EijI(21,3)-EijI(24,3)-EijI(3,2)*P(25)+EijI(8,2)*P(68)+E
     -   ijI(7,2)*P(441))+EijI(5,2)*P(555)+EijI(6,2)*P(572)+EijI(10,2
     -   )*P(573))
       F(265)=DCMPLX(FR(265),FI(265))
       FR(266) = 4*(D01345R-D02345R+Dij1345R(1,1)-Dij1345R(5,2)+Dij134
     -   5R(6,2)-Dij2345R(2,1)-s12*EijR(2,1)-p2sq*EijR(9,2)+2*(EijR(1
     -   1,2)-EijR(24,3))+(EijR(6,2)+2*EijR(7,2))*P(18)+EijR(3,2)*P(2
     -   5)+EijR(3,1)*P(30)-EijR(10,2)*P(58)-EijR(8,2)*P(68)+EijR(4,1
     -   )*P(144))
       FI(266) = 4*(D01345I-D02345I+Dij1345I(1,1)-Dij1345I(5,2)+Dij134
     -   5I(6,2)-Dij2345I(2,1)-s12*EijI(2,1)-p2sq*EijI(9,2)+2*(EijI(1
     -   1,2)-EijI(24,3))+(EijI(6,2)+2*EijI(7,2))*P(18)+EijI(3,2)*P(2
     -   5)+EijI(3,1)*P(30)-EijI(10,2)*P(58)-EijI(8,2)*P(68)+EijI(4,1
     -   )*P(144))
       F(266)=DCMPLX(FR(266),FI(266))
       P(574) = -s15+s34+2*P(18)
       FR(267) = 4*(D02345R+Dij2345R(2,1)+Dij2345R(3,1)+Dij2345R(3,2)-
     -   Dij2345R(6,2)+s12*EijR(2,1)-2*EijR(11,2)-6*EijR(24,3)-EijR(6
     -   ,2)*P(18)-EijR(3,2)*P(25)-EijR(3,1)*P(30)+(-EijR(14,3)+EijR(
     -   19,3))*P(36)-EijR(13,3)*P(43)-EijR(4,1)*P(61)+EijR(8,2)*P(68
     -   )-EijR(4,3)*P(174)+(-EijR(15,3)+EijR(20,3))*P(210)-EijR(7,2)
     -   *P(231)-EijR(10,2)*P(436)+EijR(16,3)*P(548)+EijR(9,2)*P(574)
     -   )
       FI(267) = 4*(D02345I+Dij2345I(2,1)+Dij2345I(3,1)+Dij2345I(3,2)-
     -   Dij2345I(6,2)+s12*EijI(2,1)-2*EijI(11,2)-6*EijI(24,3)-EijI(6
     -   ,2)*P(18)-EijI(3,2)*P(25)-EijI(3,1)*P(30)+(-EijI(14,3)+EijI(
     -   19,3))*P(36)-EijI(13,3)*P(43)-EijI(4,1)*P(61)+EijI(8,2)*P(68
     -   )-EijI(4,3)*P(174)+(-EijI(15,3)+EijI(20,3))*P(210)-EijI(7,2)
     -   *P(231)-EijI(10,2)*P(436)+EijI(16,3)*P(548)+EijI(9,2)*P(574)
     -   )
       F(267)=DCMPLX(FR(267),FI(267))
       P(575) = 3*p2sq-s12-s15+s34
       P(576) = -s15+3*P(18)
       P(577) = -s15+s23-s45+2*P(18)
       P(578) = -p2sq+s23+s45-2*P(147)
       P(579) = p2sq-p3sq+s34
       P(580) = -s12-s15+s23+s45+2*P(579)
       P(581) = p4sq-s34+s45+2*P(42)
       FR(268) = -4*(-D02345R-Dij2345R(3,2)-Dij2345R(4,2)+Dij2345R(5,2
     -   )+Dij2345R(6,2)+12*EijR(11,2)-6*EijR(22,3)+2*(Dij2345R(1,1)-
     -   Dij2345R(3,1)+p3sq*EijR(3,2)-EijR(23,3))+8*EijR(24,3)-(-EijR
     -   (14,3)-EijR(17,3)+EijR(18,3)+EijR(19,3))*P(36)-(EijR(12,3)-E
     -   ijR(13,3))*P(43)+EijR(6,2)*P(45)+Is45*(Cij123R(2,1)+Cij123R(
     -   2,2)-p3sq*(Dij1234R(3,1)+Dij1234R(6,2))+2*(Dij1234R(7,2)+Dij
     -   1234R(12,3))+Dij1234R(2,2)*P(9)-Dij1234R(4,2)*P(45))-EijR(4,
     -   1)*P(153)+(EijR(4,2)+EijR(4,3))*P(174)-(-EijR(9,3)+EijR(10,3
     -   ))*P(210)-Dij1234R(2,1)*P(229)-EijR(2,1)*P(230)+EijR(15,3)*P
     -   (533)-EijR(16,3)*P(548)+EijR(20,3)*P(549)+EijR(2,2)*P(575)-E
     -   ijR(5,2)*P(576)+EijR(7,2)*P(577)+EijR(8,2)*P(578)-EijR(9,2)*
     -   P(580)-EijR(10,2)*P(581))
       FI(268) = -4*(-D02345I-Dij2345I(3,2)-Dij2345I(4,2)+Dij2345I(5,2
     -   )+Dij2345I(6,2)+12*EijI(11,2)-6*EijI(22,3)+2*(Dij2345I(1,1)-
     -   Dij2345I(3,1)+p3sq*EijI(3,2)-EijI(23,3))+8*EijI(24,3)-(-EijI
     -   (14,3)-EijI(17,3)+EijI(18,3)+EijI(19,3))*P(36)-(EijI(12,3)-E
     -   ijI(13,3))*P(43)+EijI(6,2)*P(45)+Is45*(Cij123I(2,1)+Cij123I(
     -   2,2)-p3sq*(Dij1234I(3,1)+Dij1234I(6,2))+2*(Dij1234I(7,2)+Dij
     -   1234I(12,3))+Dij1234I(2,2)*P(9)-Dij1234I(4,2)*P(45))-EijI(4,
     -   1)*P(153)+(EijI(4,2)+EijI(4,3))*P(174)-(-EijI(9,3)+EijI(10,3
     -   ))*P(210)-Dij1234I(2,1)*P(229)-EijI(2,1)*P(230)+EijI(15,3)*P
     -   (533)-EijI(16,3)*P(548)+EijI(20,3)*P(549)+EijI(2,2)*P(575)-E
     -   ijI(5,2)*P(576)+EijI(7,2)*P(577)+EijI(8,2)*P(578)-EijI(9,2)*
     -   P(580)-EijI(10,2)*P(581))
       F(268)=DCMPLX(FR(268),FI(268))
       P(582) = 2*p2sq-P(180)
       P(583) = s15-2*P(18)
       P(584) = p4sq+3*s45+2*P(42)
       P(585) = p3sq+2*p4sq-s12+3*s45
       FR(269) = -4*(Dij2345R(2,1)-Dij2345R(2,2)-Dij2345R(3,1)-Dij2345
     -   R(3,2)+2*Dij2345R(6,2)-s12*EijR(3,2)+s12*EijR(10,2)-8*(EijR(
     -   23,3)-EijR(24,3))+(EijR(11,3)+EijR(14,3)-2*EijR(19,3))*P(36)
     -   -EijR(3,3)*P(43)-EijR(4,1)*P(153)+EijR(4,3)*P(174)+(EijR(12,
     -   3)+EijR(15,3)-2*EijR(20,3))*P(210)+Is45*(C0134R+p2sq*Dij1234
     -   R(2,1)-p3sq*Dij1234R(3,2)-2*(Dij1234R(7,2)-Dij1234R(13,3))+D
     -   ij1234R(6,2)*P(9)-Dij1234R(5,2)*P(45)-Dij1234R(3,1)*P(214))-
     -   EijR(3,1)*P(230)-(-EijR(8,2)+EijR(9,2))*P(574)-EijR(6,2)*P(5
     -   82)-EijR(7,2)*P(583)+EijR(13,3)*P(584)-EijR(16,3)*P(585))
       FI(269) = -4*(Dij2345I(2,1)-Dij2345I(2,2)-Dij2345I(3,1)-Dij2345
     -   I(3,2)+2*Dij2345I(6,2)-s12*EijI(3,2)+s12*EijI(10,2)-8*(EijI(
     -   23,3)-EijI(24,3))+(EijI(11,3)+EijI(14,3)-2*EijI(19,3))*P(36)
     -   -EijI(3,3)*P(43)-EijI(4,1)*P(153)+EijI(4,3)*P(174)+(EijI(12,
     -   3)+EijI(15,3)-2*EijI(20,3))*P(210)+Is45*(C0134I+p2sq*Dij1234
     -   I(2,1)-p3sq*Dij1234I(3,2)-2*(Dij1234I(7,2)-Dij1234I(13,3))+D
     -   ij1234I(6,2)*P(9)-Dij1234I(5,2)*P(45)-Dij1234I(3,1)*P(214))-
     -   EijI(3,1)*P(230)-(-EijI(8,2)+EijI(9,2))*P(574)-EijI(6,2)*P(5
     -   82)-EijI(7,2)*P(583)+EijI(13,3)*P(584)-EijI(16,3)*P(585))
       F(269)=DCMPLX(FR(269),FI(269))
       P(586) = 4*s12+s15-s34
       P(587) = -s23-3*P(52)+2*P(492)
       P(588) = p2sq-p3sq-s12+s34
       P(589) = s15+s23-3*s45-2*P(588)
       P(590) = 2*p2sq+p3sq-3*s12-s15-s23+s34
       P(591) = 3*p3sq+p4sq-2*s34
       P(592) = p3sq+p4sq-s12-s23+3*s45
       FR(270) = 2*Is12s45*(3*B013R+B014R-2*(p2sq*(C0123R+Cij123R(1,1)
     -   )+p4sq*Cij145R(2,1)+2*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2
     -   ))+C0145R*P(15)-Cij123R(2,1)*P(18)+Cij134R(2,1)*P(25)+C0134R
     -   *P(43)))+4*(Dij2345R(3,1)+Dij2345R(3,2)-Dij2345R(6,2)-s15*(E
     -   E0R+2*EijR(1,1)+EijR(1,2))-s12*EijR(3,1)-2*(D02345R+Dij2345R
     -   (2,1)-6*EijR(11,2)-3*EijR(21,3)-EijR(23,3)+4*EijR(24,3))+(-E
     -   ijR(6,3)+EijR(7,3))*P(36)+D01345R*P(39)+(D01234R+2*Dij1234R(
     -   1,1)+Dij1234R(1,2))*P(40)+(EijR(11,3)-EijR(13,3))*P(43)+EijR
     -   (3,2)*P(169)+(EijR(4,2)-EijR(4,3))*P(174)+(-EijR(15,3)-EijR(
     -   17,3)+EijR(18,3)+EijR(20,3))*P(210)-EijR(4,1)*P(230)+Is45*(-
     -   Cij123R(3,2)+2*(Cij134R(1,1)-p2sq*Dij1234R(6,2)-Dij1234R(11,
     -   3))+Dij1234R(3,2)*P(1)+Dij1234R(3,1)*P(42)+Dij1234R(2,1)*P(2
     -   08)+Dij1234R(4,2)*P(262)-Dij1234R(5,2)*P(263))+EijR(2,1)*P(2
     -   65)+Is12*(Cij145R(1,1)+2*p4sq*Dij1345R(3,2)+6*Dij1345R(7,2)+
     -   Dij1345R(3,1)*P(15)+(Dij1345R(4,2)-Dij1345R(5,2))*P(261)+Dij
     -   1345R(2,1)*P(284)+Dij1345R(2,2)*P(287)-Dij1345R(6,2)*P(293))
     -   -EijR(8,2)*P(538)+EijR(14,3)*P(544)+EijR(16,3)*P(548)+Dij134
     -   5R(1,1)*P(550)+EijR(5,2)*P(586)+EijR(6,2)*P(587)+EijR(7,2)*P
     -   (589)+EijR(9,2)*P(590)-EijR(10,2)*P(591)-EijR(19,3)*P(592))
       FI(270) = 2*Is12s45*(3*B013I+B014I-2*(p2sq*(C0123I+Cij123I(1,1)
     -   )+p4sq*Cij145I(2,1)+2*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2
     -   ))+C0145I*P(15)-Cij123I(2,1)*P(18)+Cij134I(2,1)*P(25)+C0134I
     -   *P(43)))+4*(Dij2345I(3,1)+Dij2345I(3,2)-Dij2345I(6,2)-s15*(E
     -   E0I+2*EijI(1,1)+EijI(1,2))-s12*EijI(3,1)-2*(D02345I+Dij2345I
     -   (2,1)-6*EijI(11,2)-3*EijI(21,3)-EijI(23,3)+4*EijI(24,3))+(-E
     -   ijI(6,3)+EijI(7,3))*P(36)+D01345I*P(39)+(D01234I+2*Dij1234I(
     -   1,1)+Dij1234I(1,2))*P(40)+(EijI(11,3)-EijI(13,3))*P(43)+EijI
     -   (3,2)*P(169)+(EijI(4,2)-EijI(4,3))*P(174)+(-EijI(15,3)-EijI(
     -   17,3)+EijI(18,3)+EijI(20,3))*P(210)-EijI(4,1)*P(230)+Is45*(-
     -   Cij123I(3,2)+2*(Cij134I(1,1)-p2sq*Dij1234I(6,2)-Dij1234I(11,
     -   3))+Dij1234I(3,2)*P(1)+Dij1234I(3,1)*P(42)+Dij1234I(2,1)*P(2
     -   08)+Dij1234I(4,2)*P(262)-Dij1234I(5,2)*P(263))+EijI(2,1)*P(2
     -   65)+Is12*(Cij145I(1,1)+2*p4sq*Dij1345I(3,2)+6*Dij1345I(7,2)+
     -   Dij1345I(3,1)*P(15)+(Dij1345I(4,2)-Dij1345I(5,2))*P(261)+Dij
     -   1345I(2,1)*P(284)+Dij1345I(2,2)*P(287)-Dij1345I(6,2)*P(293))
     -   -EijI(8,2)*P(538)+EijI(14,3)*P(544)+EijI(16,3)*P(548)+Dij134
     -   5I(1,1)*P(550)+EijI(5,2)*P(586)+EijI(6,2)*P(587)+EijI(7,2)*P
     -   (589)+EijI(9,2)*P(590)-EijI(10,2)*P(591)-EijI(19,3)*P(592))
       F(270)=DCMPLX(FR(270),FI(270))
       FR(271) = -2*(-(s15*(EE0R+EijR(1,1)))-D01345R*P(92)+Is45*(2*(Ci
     -   j134R(1,1)+p2sq*Dij1234R(2,2))+6*Dij1234R(7,2)+Dij1234R(3,2)
     -   *P(1)-Dij1234R(3,1)*P(11)-(Dij1234R(2,1)+Dij1234R(4,2)-Dij12
     -   34R(5,2))*P(18)-Dij1234R(6,2)*P(157))-(D01234R+Dij1234R(1,1)
     -   )*P(274)+Is12*(Cij145R(1,1)+2*p4sq*Dij1345R(3,2)+6*Dij1345R(
     -   7,2)+Dij1345R(3,1)*P(15)+(-Dij1345R(4,2)+Dij1345R(5,2))*P(95
     -   )-Dij1345R(2,1)*P(106)+Dij1345R(2,2)*P(287)-Dij1345R(6,2)*P(
     -   293))-3*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)-2*p3sq*EijR(3,2
     -   )-10*EijR(11,2)+EijR(5,2)*P(18)-EijR(7,2)*P(36)-EijR(6,2)*P(
     -   45)-EijR(2,2)*P(68)-EijR(4,2)*P(174)+EijR(10,2)*P(407)+EijR(
     -   9,2)*P(411)+EijR(8,2)*P(413))+Dij1345R(1,1)*P(550))+Is12s45*
     -   (-3*B013R-B014R+4*(Cij123R(4,2)+Cij134R(4,2)+Cij145R(4,2))-2
     -   *(p2sq*Cij123R(2,1)-p4sq*Cij145R(2,1)-C0145R*P(15)-(C0123R+C
     -   ij123R(1,1))*P(18)-Cij134R(2,1)*P(25)-C0134R*P(551)))
       FI(271) = -2*(-(s15*(EE0I+EijI(1,1)))-D01345I*P(92)+Is45*(2*(Ci
     -   j134I(1,1)+p2sq*Dij1234I(2,2))+6*Dij1234I(7,2)+Dij1234I(3,2)
     -   *P(1)-Dij1234I(3,1)*P(11)-(Dij1234I(2,1)+Dij1234I(4,2)-Dij12
     -   34I(5,2))*P(18)-Dij1234I(6,2)*P(157))-(D01234I+Dij1234I(1,1)
     -   )*P(274)+Is12*(Cij145I(1,1)+2*p4sq*Dij1345I(3,2)+6*Dij1345I(
     -   7,2)+Dij1345I(3,1)*P(15)+(-Dij1345I(4,2)+Dij1345I(5,2))*P(95
     -   )-Dij1345I(2,1)*P(106)+Dij1345I(2,2)*P(287)-Dij1345I(6,2)*P(
     -   293))-3*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)-2*p3sq*EijI(3,2
     -   )-10*EijI(11,2)+EijI(5,2)*P(18)-EijI(7,2)*P(36)-EijI(6,2)*P(
     -   45)-EijI(2,2)*P(68)-EijI(4,2)*P(174)+EijI(10,2)*P(407)+EijI(
     -   9,2)*P(411)+EijI(8,2)*P(413))+Dij1345I(1,1)*P(550))+Is12s45*
     -   (-3*B013I-B014I+4*(Cij123I(4,2)+Cij134I(4,2)+Cij145I(4,2))-2
     -   *(p2sq*Cij123I(2,1)-p4sq*Cij145I(2,1)-C0145I*P(15)-(C0123I+C
     -   ij123I(1,1))*P(18)-Cij134I(2,1)*P(25)-C0134I*P(551)))
       F(271)=DCMPLX(FR(271),FI(271))
       KI(1) = EijI(4,4)-2*EijI(15,4)+EijI(21,4)
       KR(1) = EijR(4,4)-2*EijR(15,4)+EijR(21,4)
       K(1)=DCMPLX(KR(1),KI(1))
       KI(2) = EijI(4,4)-EijI(13,4)-3*(EijI(16,4)-EijI(22,4))
       KR(2) = EijR(4,4)-EijR(13,4)-3*(EijR(16,4)-EijR(22,4))
       K(2)=DCMPLX(KR(2),KI(2))
       KI(3) = EijI(4,4)-2*EijI(16,4)+EijI(22,4)
       KR(3) = EijR(4,4)-2*EijR(16,4)+EijR(22,4)
       K(3)=DCMPLX(KR(3),KI(3))
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
       p1mup4 = dotrc(p1,mup4)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup4 = dotrc(p2,mup4)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup4 = dotrc(p3,mup4)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup4 = dotrc(p4,mup4)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrc(p5,mup3)
       p5mup4 = dotrc(p5,mup4)
       mup2mup3=dotcc(mup2,mup3)
       mup2mup4=dotcc(mup2,mup4)
       mup3mup4=dotcc(mup3,mup4)
c      Print*," p1mup2 ",  p1mup2 
c      Print*," p1mup3 ",  p1mup3 
c      Print*," p1mup4 ",  p1mup4 
c      Print*," p2mup2 ",  p2mup2 
c      Print*," p2mup3 ",  p2mup3 
c      Print*," p2mup4 ",  p2mup4 
c      Print*," p3mup2 ",  p3mup2 
c      Print*," p3mup3 ",  p3mup3 
c      Print*," p3mup4 ",  p3mup4 
c      Print*," p4mup2 ",  p4mup2 
c      Print*," p4mup3 ",  p4mup3 
c      Print*," p4mup4 ",  p4mup4 
c      Print*," p5mup2 ",  p5mup2 
c      Print*," p5mup3 ",  p5mup3 
c      Print*," p5mup4 ",  p5mup4 
c       Print*," mup2mup3", mup2mup3  
c       Print*," mup2mup4",mup2mup4   
c       Print*,"mup3mup4", mup3mup4  
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call FaFunctionmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2mup
     -   4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,p5m
     -   up3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
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
       SMB(1) = SC1c(barpsi_p5,mup2,psi_p1,alpha)
       SMB(2) = SC1c(barpsi_p5,mup3,psi_p1,alpha)
       SMB(3) = SC1c(barpsi_p5,mup4,psi_p1,alpha)
       SMB(4) = SC1r(barpsi_p5,p2,psi_p1,alpha)
       SMB(5) = SC1r(barpsi_p5,p3,psi_p1,alpha)
       SMB(6) = SC3ccc(barpsi_p5,mup2,mup3,mup4,psi_p1,alpha)
       SMB(7) = SC3rcc(barpsi_p5,p2,mup2,mup3,psi_p1,alpha)
       SMB(8) = SC3rcc(barpsi_p5,p2,mup2,mup4,psi_p1,alpha)
       SMB(9) = SC3rcc(barpsi_p5,p2,mup3,mup4,psi_p1,alpha)
       SMB(10) = SC3rrc(barpsi_p5,p2,p3,mup2,psi_p1,alpha)
       SMB(11) = SC3rrc(barpsi_p5,p2,p3,mup3,psi_p1,alpha)
       SMB(12) = SC3rrc(barpsi_p5,p2,p3,mup4,psi_p1,alpha)
       SMB(13) = SC3rcc(barpsi_p5,p3,mup2,mup3,psi_p1,alpha)
       SMB(14) = SC3rcc(barpsi_p5,p3,mup2,mup4,psi_p1,alpha)
       SMB(15) = SC3rcc(barpsi_p5,p3,mup3,mup4,psi_p1,alpha)
       SMB(16) = SC5rrccc(barpsi_p5,p2,p3,mup2,mup3,mup4,psi_p1,alpha)
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
     -   (5)*SMB(5)+F(196)*SMB(6)+Fa(6)*SMB(7)+Fa(7)*SMB(8)+Fa(8)*SMB
     -   (9)+Fa(9)*SMB(10)+Fa(10)*SMB(11)+Fa(11)*SMB(12)+Fa(12)*SMB(1
     -   3)+Fa(13)*SMB(14)+Fa(14)*SMB(15)+F(271)*SMB(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       result =-result
       resultb = (2*(mup3mup4*SMB(1)-mup2mup4*SMB(2)+mup2mup3*SMB(3))-
     -   SMB(6))/s45+(4*((mup2mup4*(p1mup3+p2mup3)+mup2mup3*p3mup4)*S
     -   MB(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2mup3*p2mup4)*SMB(
     -   5))+(p2sq+p3sq-s23)*SMB(6)-2*((2*(p1mup3*p2mup4+p2mup3*(p2mu
     -   p4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMB(1)-(-2*p1mup2*p3mu
     -   p4+mup2mup4*(p2sq+p3sq-s23))*SMB(2)+(2*p1mup2*(p1mup3+p2mup3
     -   )+mup2mup3*(p2sq+p3sq-s23))*SMB(3)+p3mup4*SMB(7)+(p1mup3+p2m
     -   up3)*SMB(8)+mup3mup4*SMB(10)-mup2mup4*SMB(11)+mup2mup3*SMB(1
     -   2)-p2mup4*SMB(13)+p2mup3*SMB(14)+p1mup2*SMB(15))+SMB(16))/(s
     -   12*s45)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c                         CHECK OF GAUGE INVARIANCE             
      If(gaugetest.eq.0) then 
cc      The mu_p2,mu_p3,mu_p4 must be replaced for the incoming moment 
c      mup2->p2,mup3->p2+p3+p4,mup4->p4. SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2 = dotrr(p1,p2)
       p1mup3 = dotrr(p1,p2)+dotrr(p1,p3) +dotrr(p1,p4)
       p1mup4 = dotrr(p1,p4)
       p2mup2 = dotrr(p2,p2)
       p2mup3 = dotrr(p2,p2)+dotrr(p2,p3) +dotrr(p2,p4)
       p2mup4 = dotrr(p2,p4)
       p3mup2 = dotrr(p3,p2)
       p3mup3 = dotrr(p3,p2)+dotrr(p3,p3) +dotrr(p3,p4)
       p3mup4 = dotrr(p3,p4)
       p4mup2 = dotrr(p4,p2)
       p4mup3 = dotrr(p4,p2)+dotrr(p4,p3) +dotrr(p4,p4)
       p4mup4 = dotrr(p4,p4)
       p5mup2 = dotrr(p5,p2)
       p5mup3 = dotrr(p5,p2)+dotrr(p5,p3) +dotrr(p5,p4)
       p5mup4 = dotrr(p5,p4)
       mup2mup3=dotrr(p2,p2)+dotrr(p2,p3)+dotrr(p2,p4)
       mup2mup4=dotrr(p2,p4)
       mup3mup4=dotrr(p2,p4)+dotrr(p3,p4)+dotrr(p4,p4)
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
c************** Calling the Fa functions**********************************************************************
       call FaFunctionmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2mup
     -   4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,p5m
     -   up3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2,4) and mu_p3=p2+p3+p4
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(4)
       SMG(2) = 0
       SMG(3) = -SMB(4)-SMB(5)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = (-s12-s23+s34+s45)*SMB(4)-(s12+s15-s34)*SMB(5)
       SMG(7) = 0
       SMG(8) = -(p2sq*(SMB(4)+SMB(5)))
       SMG(9) = (-s12-s23+s34+s45)*SMB(4)-(s12+s15-s34)*SMB(5)
       SMG(10) = -((p2sq+p3sq-s23)*SMB(4))-p2sq*SMB(5)
       SMG(11) = (-p4sq-s12+s34+s45)*SMB(4)-(s12+s15-s34)*SMB(5)
       SMG(12) = (p2sq-p4sq-s12-s23+s34+s45)*SMB(4)+(p2sq-s12-s15+s34)
     -   *SMB(5)
       SMG(13) = (p4sq+s12-s34-s45)*SMB(4)+(s12+s15-s34)*SMB(5)
       SMG(14) = (p3sq+p4sq+s12-s34-s45)*SMB(4)+(p3sq+s12+s15-s23-s34)
     -   *SMB(5)
       SMG(15) = (p4sq+s12-s34-s45)*SMB(4)+(p4sq+s12+s15-s23-s34)*SMB(
     -   5)
       SMG(16) = (p2sq*(-p4sq+s23)+(p3sq-s23)*(s12+s23-s34-s45))*SMB(4
     -   )+(p2sq*(-p4sq+s23)+(p3sq-s23)*(s12+s15-s34))*SMB(5)
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
c    This should be zero:
c************************************************************************************
c************************************************************************************
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*S
     -   MG(4)+Fa(5)*SMG(5)+F(196)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(271)*SMG(16)
       resultgauge(1) =resultgauge(1)/(Sqrt(abs(p2sq)*abs(s15)*abs(p4s
     -   q)))
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
       resultgaugeb(1) = (2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2mup3*
     -   SMG(3))-SMG(6))/s45+(4*((mup2mup4*(p1mup3+p2mup3)+mup2mup3*p
     -   3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2mup3*p2mu
     -   p4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2mup4+p2mu
     -   p3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1)-(-2*p1m
     -   up2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mup2*(p1mup
     -   3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*SMG(7)+(p1
     -   mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(11)+mup2mu
     -   p3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG(15))+SMG
     -   (16))/(s12*s45)
       resultgaugeb(1) =resultgaugeb(1)/(Sqrt(abs(p2sq)*abs(s15)*abs(p
     -   4sq)))
           resultgauge(2)=(0d0,0d0)      
           resultgauge(3)=(0d0,0d0)      
           resultgaugeb(2)=(0d0,0d0)      
           resultgaugeb(3)=(0d0,0d0)     
           return      
           endif      
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
       mup2mup4=dotrc(p2,mup4)
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
       call FaFunctionmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2mup
     -   4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,p5m
     -   up3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(4)
       SMG(2) = SMB(2)
       SMG(3) = SMB(3)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = SMB(9)
       SMG(7) = p2sq*SMB(2)
       SMG(8) = p2sq*SMB(3)
       SMG(9) = SMB(9)
       SMG(10) = -((p2sq+p3sq-s23)*SMB(4))-p2sq*SMB(5)
       SMG(11) = SMB(11)
       SMG(12) = SMB(12)
       SMG(13) = -((p2sq+p3sq-s23)*SMB(2))-SMB(11)
       SMG(14) = -((p2sq+p3sq-s23)*SMB(3))-SMB(12)
       SMG(15) = SMB(15)
       SMG(16) = -((p2sq+p3sq-s23)*SMB(9))-p2sq*SMB(15)
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
     -   MG(4)+Fa(5)*SMG(5)+F(196)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(271)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
       resultgaugeb(1) = (2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2mup3*
     -   SMG(3))-SMG(6))/s45+(4*((mup2mup4*(p1mup3+p2mup3)+mup2mup3*p
     -   3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2mup3*p2mu
     -   p4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2mup4+p2mu
     -   p3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1)-(-2*p1m
     -   up2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mup2*(p1mup
     -   3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*SMG(7)+(p1
     -   mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(11)+mup2mu
     -   p3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG(15))+SMG
     -   (16))/(s12*s45)
           resultgauge(2)=(0d0,0d0)      
           resultgauge(3)=(0d0,0d0)      
           resultgaugeb(2)=(0d0,0d0)      
           resultgaugeb(3)=(0d0,0d0)      
           return      
           endif      
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      If(gaugetest.eq.2) then 
cc      The mu_p3 is replaced for the incoming moment 
c      mup3->p3 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup3=dotrr(p1,p3)
       p2mup3=dotrr(p2,p3)
       p3mup3=dotrr(p3,p3)
       p4mup3=dotrr(p4,p3)
       p5mup3=dotrr(p5,p3)
       mup2mup3=dotrc(p3,mup2)
       mup3mup4=dotrc(p3,mup4)
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
       call FaFunctionmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2mup
     -   4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,p5m
     -   up3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:3)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = SMB(5)
       SMG(3) = SMB(3)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = 2*p3mup2*SMB(3)-SMB(14)
       SMG(7) = 2*p3mup2*SMB(4)-SMB(10)
       SMG(8) = SMB(8)
       SMG(9) = SMB(12)
       SMG(10) = SMB(10)
       SMG(11) = p3sq*SMB(4)
       SMG(12) = SMB(12)
       SMG(13) = -(p3sq*SMB(1))+2*p3mup2*SMB(5)
       SMG(14) = SMB(14)
       SMG(15) = p3sq*SMB(3)
       SMG(16) = -(p3sq*SMB(8))+2*p3mup2*SMB(12)
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
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*S
     -   MG(4)+Fa(5)*SMG(5)+F(196)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(271)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
       resultgaugeb(1) = (2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2mup3*
     -   SMG(3))-SMG(6))/s45+(4*((mup2mup4*(p1mup3+p2mup3)+mup2mup3*p
     -   3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2mup3*p2mu
     -   p4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2mup4+p2mu
     -   p3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1)-(-2*p1m
     -   up2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mup2*(p1mup
     -   3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*SMG(7)+(p1
     -   mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(11)+mup2mu
     -   p3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG(15))+SMG
     -   (16))/(s12*s45)
           resultgauge(2)=(0d0,0d0)      
           resultgauge(3)=(0d0,0d0)      
           resultgaugeb(2)=(0d0,0d0)      
           resultgaugeb(3)=(0d0,0d0)      
          return      
           endif      
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      If(gaugetest.eq.3) then 
cc      The mu_p4 is replaced for the incoming moment 
c      mup4->p4 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup4=dotrr(p1,p4)
       p2mup4=dotrr(p2,p4)
       p3mup4=dotrr(p3,p4)
       p4mup4=dotrr(p4,p4)
       p5mup4=dotrr(p5,p4)
       mup2mup4=dotrc(p4,mup2)
       mup3mup4=dotrc(p4,mup3)
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
       call FaFunctionmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2mup
     -   4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,p5m
     -   up3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:3)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = SMB(2)
       SMG(3) = -SMB(4)-SMB(5)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = 2*(p1mup3+p4mup3)*SMB(1)-2*(p1mup2+p4mup2)*SMB(2)-SMB(
     -   7)-SMB(13)
       SMG(7) = SMB(7)
       SMG(8) = (p2sq-s12-s15+s34)*SMB(1)+2*(p1mup2+p4mup2)*SMB(4)+SMB
     -   (10)
       SMG(9) = (p2sq-s12-s15+s34)*SMB(2)+2*(p1mup3+p4mup3)*SMB(4)+SMB
     -   (11)
       SMG(10) = SMB(10)
       SMG(11) = SMB(11)
       SMG(12) = (p2sq-p4sq-s12-s23+s34+s45)*SMB(4)+(p2sq-s12-s15+s34)
     -   *SMB(5)
       SMG(13) = SMB(13)
       SMG(14) = (-p2sq+p4sq+s12+s23-s34-s45)*SMB(1)+2*(p1mup2+p4mup2)
     -   *SMB(5)-SMB(10)
       SMG(15) = (-p2sq+p4sq+s12+s23-s34-s45)*SMB(2)+2*(p1mup3+p4mup3)
     -   *SMB(5)-SMB(11)
       SMG(16) = (p2sq-p4sq-s12-s23+s34+s45)*SMB(7)+2*((p1mup3+p4mup3)
     -   *SMB(10)-(p1mup2+p4mup2)*SMB(11))+(p2sq-s12-s15+s34)*SMB(13)
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
c    This should be equal to the subsctraction of boxline with momenta -(p1,p2,p3+p4,p5):
c************************************************************************************
c************************************************************************************
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*S
     -   MG(4)+Fa(5)*SMG(5)+F(196)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(271)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
       resultgaugeb(1) = (2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2mup3*
     -   SMG(3))-SMG(6))/s45+(4*((mup2mup4*(p1mup3+p2mup3)+mup2mup3*p
     -   3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2mup3*p2mu
     -   p4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2mup4+p2mu
     -   p3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1)-(-2*p1m
     -   up2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mup2*(p1mup
     -   3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*SMG(7)+(p1
     -   mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(11)+mup2mu
     -   p3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG(15))+SMG
     -   (16))/(s12*s45)
           resultgauge(2)=(0d0,0d0)      
           resultgauge(3)=(0d0,0d0)      
           resultgaugeb(2)=(0d0,0d0)      
           resultgaugeb(3)=(0d0,0d0)      
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
       mup2mup4=dotrc(p2,mup4)
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
       call FaFunctionmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2mup
     -   4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,p5m
     -   up3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(4)
       SMG(2) = SMB(2)
       SMG(3) = SMB(3)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = SMB(9)
       SMG(7) = p2sq*SMB(2)
       SMG(8) = p2sq*SMB(3)
       SMG(9) = SMB(9)
       SMG(10) = -((p2sq+p3sq-s23)*SMB(4))-p2sq*SMB(5)
       SMG(11) = SMB(11)
       SMG(12) = SMB(12)
       SMG(13) = -((p2sq+p3sq-s23)*SMB(2))-SMB(11)
       SMG(14) = -((p2sq+p3sq-s23)*SMB(3))-SMB(12)
       SMG(15) = SMB(15)
       SMG(16) = -((p2sq+p3sq-s23)*SMB(9))-p2sq*SMB(15)
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
     -   MG(4)+Fa(5)*SMG(5)+F(196)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(271)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
       resultgaugeb(1) = (2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2mup3*
     -   SMG(3))-SMG(6))/s45+(4*((mup2mup4*(p1mup3+p2mup3)+mup2mup3*p
     -   3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2mup3*p2mu
     -   p4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2mup4+p2mu
     -   p3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1)-(-2*p1m
     -   up2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mup2*(p1mup
     -   3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*SMG(7)+(p1
     -   mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(11)+mup2mu
     -   p3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG(15))+SMG
     -   (16))/(s12*s45)
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
       p1mup4 = dotrc(p1,mup4)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrr(p2,p3)
       p2mup4 = dotrc(p2,mup4)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrr(p3,p3)
       p3mup4 = dotrc(p3,mup4)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrr(p4,p3)
       p4mup4 = dotrc(p4,mup4)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrr(p5,p3)
       p5mup4 = dotrc(p5,mup4)
       mup2mup3=dotrc(p3,mup2)
       mup3mup4=dotrc(p3,mup4)
       mup2mup4=dotcc(mup2,mup4)
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
       call FaFunctionmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2mup
     -   4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,p5m
     -   up3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:3)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = SMB(5)
       SMG(3) = SMB(3)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = 2*p3mup2*SMB(3)-SMB(14)
       SMG(7) = 2*p3mup2*SMB(4)-SMB(10)
       SMG(8) = SMB(8)
       SMG(9) = SMB(12)
       SMG(10) = SMB(10)
       SMG(11) = p3sq*SMB(4)
       SMG(12) = SMB(12)
       SMG(13) = -(p3sq*SMB(1))+2*p3mup2*SMB(5)
       SMG(14) = SMB(14)
       SMG(15) = p3sq*SMB(3)
       SMG(16) = -(p3sq*SMB(8))+2*p3mup2*SMB(12)
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
     -   MG(4)+Fa(5)*SMG(5)+F(196)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(271)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(2) =-resultgauge(2)
       resultgaugeb(2) = (2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2mup3*
     -   SMG(3))-SMG(6))/s45+(4*((mup2mup4*(p1mup3+p2mup3)+mup2mup3*p
     -   3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2mup3*p2mu
     -   p4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2mup4+p2mu
     -   p3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1)-(-2*p1m
     -   up2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mup2*(p1mup
     -   3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*SMG(7)+(p1
     -   mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(11)+mup2mu
     -   p3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG(15))+SMG
     -   (16))/(s12*s45)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c      The mu_p4 is replaced for the incoming moment 
c      mup4->p4 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrc(p1,mup3)
       p1mup4 = dotrr(p1,p4)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p2mup4 = dotrr(p2,p4)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p3mup4 = dotrr(p3,p4)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       p4mup4 = dotrr(p4,p4)
       p5mup2 = dotrc(p5,mup2)
       p5mup3 = dotrc(p5,mup3)
       p5mup4 = dotrr(p5,p4)
       mup2mup3=dotcc(mup2,mup3)
       mup3mup4=dotrc(p4,mup3)
       mup2mup4=dotrc(p4,mup2)
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
       call FaFunctionmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2mup
     -   4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,p5m
     -   up3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:3)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = SMB(2)
       SMG(3) = -SMB(4)-SMB(5)
       SMG(4) = SMB(4)
       SMG(5) = SMB(5)
       SMG(6) = 2*(p1mup3+p4mup3)*SMB(1)-2*(p1mup2+p4mup2)*SMB(2)-SMB(
     -   7)-SMB(13)
       SMG(7) = SMB(7)
       SMG(8) = (p2sq-s12-s15+s34)*SMB(1)+2*(p1mup2+p4mup2)*SMB(4)+SMB
     -   (10)
       SMG(9) = (p2sq-s12-s15+s34)*SMB(2)+2*(p1mup3+p4mup3)*SMB(4)+SMB
     -   (11)
       SMG(10) = SMB(10)
       SMG(11) = SMB(11)
       SMG(12) = (p2sq-p4sq-s12-s23+s34+s45)*SMB(4)+(p2sq-s12-s15+s34)
     -   *SMB(5)
       SMG(13) = SMB(13)
       SMG(14) = (-p2sq+p4sq+s12+s23-s34-s45)*SMB(1)+2*(p1mup2+p4mup2)
     -   *SMB(5)-SMB(10)
       SMG(15) = (-p2sq+p4sq+s12+s23-s34-s45)*SMB(2)+2*(p1mup3+p4mup3)
     -   *SMB(5)-SMB(11)
       SMG(16) = (p2sq-p4sq-s12-s23+s34+s45)*SMB(7)+2*((p1mup3+p4mup3)
     -   *SMB(10)-(p1mup2+p4mup2)*SMB(11))+(p2sq-s12-s15+s34)*SMB(13)
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
     -   MG(4)+Fa(5)*SMG(5)+F(196)*SMG(6)+Fa(6)*SMG(7)+Fa(7)*SMG(8)+F
     -   a(8)*SMG(9)+Fa(9)*SMG(10)+Fa(10)*SMG(11)+Fa(11)*SMG(12)+Fa(1
     -   2)*SMG(13)+Fa(13)*SMG(14)+Fa(14)*SMG(15)+F(271)*SMG(16)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(3) =-resultgauge(3)
       resultgaugeb(3) = (2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2mup3*
     -   SMG(3))-SMG(6))/s45+(4*((mup2mup4*(p1mup3+p2mup3)+mup2mup3*p
     -   3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2mup3*p2mu
     -   p4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2mup4+p2mu
     -   p3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1)-(-2*p1m
     -   up2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mup2*(p1mup
     -   3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*SMG(7)+(p1
     -   mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(11)+mup2mu
     -   p3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG(15))+SMG
     -   (16))/(s12*s45)
           endif      
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
c************** Calling the Fa functions**********************************************************************
       subroutine FaFunctionmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3
     -   ,p2mup4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mu
     -   p2,p5mup3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p2mup2, p2mup3, p2mup4, 
     -          p3mup2, p3mup3, p3mup4, p4mup2, p4mup3, p4mup4, 
     -          p5mup2, p5mup3, p5mup4
       Complex*16 mup2mup3,mup2mup4,mup3mup4
        common/Ffunctionsmmm/F,P
       common/Kfunctionsmmm/K
      Complex*16 Fa(14),F(271),K(3)
       Real*8 P(592) 
        Fa(1) = p3mup3*(p5mup4*F(3)+p1mup4*F(4)+p2mup4*F(8)+p3mup4*F(9)
     -   )+p1mup4*p5mup3*F(11)+p2mup3*(p1mup4*F(2)+p2mup4*F(5)+p3mup4
     -   *F(6)+p5mup4*F(12))+p1mup3*(p1mup4*F(1)+p2mup4*F(7)+p3mup4*F
     -   (10)+p5mup4*F(13))+p2mup4*p5mup3*F(14)+p3mup4*p5mup3*F(15)+m
     -   up3mup4*F(16)+p5mup3*p5mup4*F(17)
       Fa(2) = p2mup2*(p1mup4*F(19)+p5mup4*F(21)+p3mup4*F(26)+p2mup4*F
     -   (27))+p3mup2*(p5mup4*F(22)+p1mup4*F(25)+p3mup4*F(29)+p2mup4*
     -   F(30))+p5mup2*(p2mup4*F(20)+p5mup4*F(23)+p3mup4*F(24)+p1mup4
     -   *F(32))+p1mup2*(p1mup4*F(18)+p2mup4*F(28)+p3mup4*F(31)+p5mup
     -   4*F(33))+mup2mup4*F(34)
       Fa(3) = p5mup2*(p2mup3*F(35)+p5mup3*F(37)+p1mup3*F(40))+mup2mup
     -   3*F(42)+p1mup2*(p5mup3*F(41)+p1mup3*F(43)+p2mup3*F(44))+p2mu
     -   p2*(p5mup3*F(38)+p1mup3*F(45)+p2mup3*F(46)+p3mup3*F(49))+p3m
     -   up2*(p5mup3*F(36)+p1mup3*F(48)+p2mup3*F(50))+p3mup3*(p5mup2*
     -   F(39)+p1mup2*F(47)+p3mup2*F(51))
       Fa(4) = mup3mup4*(p5mup2*F(67)+p1mup2*F(69)+p3mup2*F(70))+mup2m
     -   up3*(p1mup4*F(71)+p2mup4*F(72)+p3mup4*F(73)+p5mup4*F(74))+mu
     -   p2mup4*(p2mup3*F(75)+p1mup3*F(76)+p5mup3*F(78))+p1mup2*(p1mu
     -   p4*p5mup3*F(58)+p3mup4*p5mup3*F(62)+p1mup3*(p1mup4*F(52)+p2m
     -   up4*F(56)+p5mup4*F(63)+p3mup4*F(79))+p2mup3*(p1mup4*F(56)+p2
     -   mup4*F(59)+p5mup4*F(64)+p3mup4*F(80))-16*(p2mup4*p5mup3*F(10
     -   1)-p5mup3*p5mup4*F(116)))+p2mup2*(mup3mup4*F(68)+p2mup3*(p2m
     -   up4*F(60)+p5mup4*F(65)+p3mup4*F(82)-16*p1mup4*F(85))+p1mup3*
     -   (p1mup4*F(53)+p3mup4*F(81)-16*(p2mup4*F(85)-p5mup4*F(108)))-
     -   16*(p1mup4*p5mup3*F(91)+p2mup4*p5mup3*F(96)-p3mup4*p5mup3*F(
     -   105)-p5mup3*p5mup4*F(111)+p3mup3*(p1mup4*F(88)+p2mup4*F(98)-
     -   p5mup4*F(113)+p3mup4*F(119))))+p3mup2*(p1mup3*(p1mup4*F(54)+
     -   p3mup4*F(84)-16*(p2mup4*F(94)-p5mup4*F(109)))+16*(p1mup4*p5m
     -   up3*F(92)+p2mup4*p5mup3*F(102)+p3mup4*p5mup3*F(106)+p5mup3*p
     -   5mup4*F(117))+p2mup3*(p5mup4*F(66)+16*(p1mup4*F(86)-p2mup4*F
     -   (95)-p3mup4*F(120))))+p3mup3*(mup2mup4*F(77)+p1mup2*(p1mup4*
     -   F(57)+p3mup4*F(83)-16*(p2mup4*F(97)-p5mup4*F(112)))-16*(p5mu
     -   p2*(p1mup4*F(90)+p2mup4*F(100)+p3mup4*F(104)+p5mup4*F(115))+
     -   p3mup2*(p1mup4*F(89)+p2mup4*F(99)-p5mup4*F(114)+p3mup4*F(121
     -   ))))+p5mup2*(p1mup3*(p1mup4*F(55)+p3mup4*F(61)-16*(p2mup4*F(
     -   91)-p5mup4*F(110)))+16*(p1mup4*p5mup3*F(93)+p3mup4*p5mup3*F(
     -   107)+p2mup3*(p1mup4*F(87)-p2mup4*F(96)-p3mup4*F(103)+p5mup4*
     -   F(111))+p5mup3*p5mup4*F(118)-p2mup4*p5mup3*K(1)))
       Fa(5) = mup3mup4*(p5mup2*F(138)+p3mup2*F(140)+p1mup2*F(141))+mu
     -   p2mup3*(p3mup4*F(142)+p1mup4*F(143)+p2mup4*F(144)+p5mup4*F(1
     -   45))+mup2mup4*(p2mup3*F(146)+p1mup3*F(148)+p5mup3*F(149))+p1
     -   mup2*(p1mup4*p5mup3*F(130)+p2mup4*p5mup3*F(132)+p1mup3*(p1mu
     -   p4*F(122)+p3mup4*F(133)+p5mup4*F(135)+p2mup4*F(150))+p2mup3*
     -   (p1mup4*F(126)+p5mup4*F(136)+p2mup4*F(151)-16*p3mup4*F(171))
     -   -16*(p3mup4*p5mup3*F(178)-p5mup3*p5mup4*F(189)))+p2mup2*(mup
     -   3mup4*F(139)+p2mup3*(p1mup4*F(127)+p3mup4*F(134)+p5mup4*F(13
     -   7)+p2mup4*F(153))+p1mup3*(p1mup4*F(123)+p2mup4*F(152)-16*(p3
     -   mup4*F(168)-p5mup4*F(180)))-16*(p1mup4*p5mup3*F(160)+p2mup4*
     -   p5mup3*F(165)+p3mup4*p5mup3*F(179)+p5mup3*p5mup4*F(190)+p3mu
     -   p3*(p1mup4*F(157)+p3mup4*F(175)-p5mup4*F(186)+p2mup4*F(193))
     -   ))+p3mup3*(mup2mup4*F(147)+p1mup2*(p1mup4*F(129)+p2mup4*F(15
     -   4)-16*(p3mup4*F(174)-p5mup4*F(185)))-16*(p5mup2*(p1mup4*F(15
     -   9)-p2mup4*F(164)+p3mup4*F(177)-p5mup4*F(188))+p3mup2*(p1mup4
     -   *F(158)+p3mup4*F(176)-p5mup4*F(187)+p2mup4*F(195))))+p3mup2*
     -   (p1mup3*(p1mup4*F(124)+p2mup4*F(155)-16*(p3mup4*F(169)-p5mup
     -   4*F(181)))-16*(-(p1mup4*p5mup3*F(161))+p2mup4*p5mup3*F(166)-
     -   p5mup3*p5mup4*F(191)+p2mup3*(p1mup4*F(156)+p3mup4*F(172)-p5m
     -   up4*F(183)+p2mup4*F(194))+p3mup4*p5mup3*K(2)))+p5mup2*(p1mup
     -   3*(p1mup4*F(125)+p2mup4*F(131)-16*(p3mup4*F(170)-p5mup4*F(18
     -   2)))+p2mup3*(p1mup4*F(128)+16*(p2mup4*F(163)+p3mup4*F(173)+p
     -   5mup4*F(184)))+16*(p1mup4*p5mup3*F(162)-p2mup4*p5mup3*F(167)
     -   +p5mup3*p5mup4*F(192)-p3mup4*p5mup3*K(3)))
       Fa(6) = p1mup4*F(197)+p2mup4*F(198)+p3mup4*F(199)+p5mup4*F(200)
       Fa(7) = p2mup3*F(201)+p1mup3*F(202)+p3mup3*F(203)+p5mup3*F(204)
       Fa(8) = p5mup2*F(205)+p2mup2*F(206)+p1mup2*F(207)+p3mup2*F(208)
       Fa(9) = p1mup4*p5mup3*F(212)+p2mup4*p5mup3*F(215)+p3mup4*p5mup3
     -   *F(216)+p1mup3*(p1mup4*F(209)+p2mup4*F(213)+p5mup4*F(217)+8*
     -   p3mup4*F(220))+p2mup3*(p1mup4*F(210)+p2mup4*F(214)+p5mup4*F(
     -   218)+8*p3mup4*F(221))+p3mup3*(p1mup4*F(211)+8*(p2mup4*F(219)
     -   +p3mup4*F(222)-p5mup4*F(223)))-8*p5mup3*p5mup4*F(224)+mup3mu
     -   p4*F(225)
       Fa(10) = p1mup2*(p5mup4*F(231)-8*(p1mup4*F(232)+p2mup4*F(233)+p
     -   3mup4*F(236)))+p2mup2*(p1mup4*F(226)+8*(p2mup4*F(234)+p3mup4
     -   *F(235)+p5mup4*F(238)))+p3mup2*(p1mup4*F(227)+8*(p2mup4*F(23
     -   5)+p3mup4*F(237)+p5mup4*F(239)))+p5mup2*(p1mup4*F(228)+p2mup
     -   4*F(229)+p3mup4*F(230)+8*p5mup4*F(240))+mup2mup4*F(241)
       Fa(11) = mup2mup3*F(245)+p1mup2*(p5mup3*F(244)+p1mup3*F(246)+p2
     -   mup3*F(247))+p5mup2*(p1mup3*F(242)+p2mup3*F(243)-8*p5mup3*F(
     -   255))+p3mup2*(p1mup3*F(250)-8*(p5mup3*F(254)-p2mup3*F(256)))
     -   +p2mup2*(p1mup3*F(248)+p2mup3*F(249)-8*(p5mup3*F(253)+p3mup3
     -   *F(257)))+p3mup3*(p1mup2*F(251)+8*(p5mup2*F(252)-p3mup2*F(25
     -   8)))
       Fa(12) = p3mup4*F(259)+p1mup4*F(260)+p2mup4*F(261)+p5mup4*F(262
     -   )
       Fa(13) = p2mup3*F(263)+p3mup3*F(264)+p1mup3*F(265)+p5mup3*F(266
     -   )
       Fa(14) = p5mup2*F(267)+p2mup2*F(268)+p3mup2*F(269)+p1mup2*F(270
     -   )
       Return
       End
