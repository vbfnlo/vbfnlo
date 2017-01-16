       subroutine penlineABEmmm(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,m
     -   up3,mup4,alpha,musqIn,ngluon,posgluon,gaugetest,comp,resultgau
     -   ge,result,resultgaugeb,resultb)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
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
       Complex*16 B0fin,C0fin,D02m_fin,D02m_fine,D01m_fin,E0fin,D0fin
       EXTERNAL dotrr,B0fin,C0fin,D01m_fin,D02m_fine,D0fin,D02m_fin,E0fin
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
       Complex*16 SMB(16),SMG(16) ,Fa(62),F(448),K(4)
       Real*8 FI(448),FR(448),PaR(23),PaI(23),KI(4),KR(4)
       Complex*16 barpsi_p5(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3) 
       Complex*16 SC1c,SC3ccc,SC3rcc,SC3rrc,SC5rrccc,SC1r,dotrc,dotcc,result(8),resultg(8)
       Complex*16 resultgauge(3),resultgaugeb(3),resultb 
       Real*8 musq,Is12,Is45,Is12s45,P(554),musqIn 
       EXTERNAL   dotrc,dotcc,SC1c,SC1r,SC3ccc,SC3rcc,SC3rrc,SC5rrccc
       Integer alpha,comp,gaugetest,ngluon,posgluon
       common/invarianetsmmm/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,Is12
       common/Ffunctionsmmm/F
       common/Pfunctionsmmm/P
       common/Kfunctionsmmm/K
       SAVE/Ffunctionsmmm/
       SAVE/Pfunctionsmmm/
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

cfc       print*, "p1", p1
cfc       print*, "p2", p2
cfc       print*, "p3", p3
cfc       print*, "p4", p4
cfc       print*, "p5", p5

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
       Is12=1d0/s12
       Is45=1d0/s45
       Is12s45=1d0/(s12*s45)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
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
c$$$       D02345=D02m_fin(s23,s34,p4sq,s15,musq)
c$$$       D01345=D02m_fine(s45,s34,s12,p4sq,musq)
c$$$       D01245=D02m_fin(s15,s45,s23,p4sq,musq)
c$$$       D01235=D01m_fin(s12,s15,s34,musq)
c$$$       D01234=D01m_fin(s12,s23,s45,musq)
       call tens_red4_new_Re_Com(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C0345R,
     -   C0245R,C0235R,C0234R,Cij345R,Cij245R,Cij235R,Cij234R,C0345I,
     -   C0245I,C0235I,C0234I,Cij345I,Cij245I,Cij235I,Cij234I,D02345,
     -   D02345R,D02345I,Dij2345R,Dij2345I)
        call tens_red4_new_Re_Com(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4,p3p
     -   4,C0345R,C0145R,C0135R,C0134R,Cij345R,Cij145R,Cij135R,Cij134
     -   R,C0345I,C0145I,C0135I,C0134I,Cij345I,Cij145I,Cij135I,Cij134
     -   I,D01345,D01345R,D01345I,Dij1345R,Dij1345I)
        call tens_red4_new_Re_Com(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4+p3p
     -   4,C0245R,C0145R,C0125R,C0124R,Cij245R,Cij145R,Cij125R,Cij124
     -   R,C0245I,C0145I,C0125I,C0124I,Cij245I,Cij145I,Cij125I,Cij124
     -   I,D01245,D01245R,D01245I,Dij1245R,Dij1245I)
        call tens_red4_new_Re_Com(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3+p2p
     -   4,C0235R,C0135R,C0125R,C0123R,Cij235R,Cij135R,Cij125R,Cij123
     -   R,C0235I,C0135I,C0125I,C0123I,Cij235I,Cij135I,Cij125I,Cij123
     -   I,D01235,D01235R,D01235I,Dij1235R,Dij1235I)
        call tens_red4_new_Re_Com(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C0234R,
     -   C0134R,C0124R,C0123R,Cij234R,Cij134R,Cij124R,Cij123R,C0234I,
     -   C0134I,C0124I,C0123I,Cij234I,Cij134I,Cij124I,Cij123I,D01234,
     -   D01234R,D01234I,Dij1234R,Dij1234I)


c$$$       write(*,*)"B13"
c$$$       write(*,*) "B13=", B013
c$$$
c$$$       write(*,*)"B14"
c$$$       write(*,*) "B14=", B014
c$$$
c$$$       write(*,*)"B15"
c$$$       write(*,*) "B15=", B015
c$$$
c$$$       write(*,*)"B34"
c$$$       write(*,*) "B34=", B034
c$$$
c$$$
c$$$       write(*,*)"B35"
c$$$       write(*,*) "B35=", B035
c$$$
c$$$
c$$$       write(*,*)"B45"
c$$$       write(*,*) "B45=", B045
c$$$
c$$$
c$$$
c$$$
c$$$
c$$$       write(*,*)"C134"
c$$$       call print_CijR(C0134,Cij134R,Cij134I)
c$$$ 
c$$$
c$$$       write(*,*)"C135"
c$$$       call print_CijR(C0135,Cij135R,Cij135I)
c$$$
c$$$
c$$$       write(*,*)"C145"
c$$$       call print_CijR(C0145,Cij145R,Cij145I)
c$$$
c$$$
c$$$       write(*,*)"C345"
c$$$       call print_CijR(C0345,Cij345R,Cij345I)
c$$$
c$$$       write(*,*)"D1345"
c$$$       call print_DijR(D01345,Dij1345R,Dij1345I)
c************************************************************************************
c************************************************************************************
       EE0=E0fin(p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
c ,p1p2,p1p3,p1p4,p2p3,p2p4,p3p4
     -   ,D02345
     -   ,D01345,D01245,D01235,D01234)
       EE0R=Dble(EE0) 
       EE0I=DIMAG(EE0) 
c************************************************************************************
c************************************************************************************
       call tens_red5_new_Re_Com(p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
     -   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,Dij2345
     -   R,Dij1345R,Dij1245R,Dij1235R,Dij1234R,D02345I,D01345I,D01245
     -   I,D01235I,D01234I,Dij2345I,Dij1345I,Dij1245I,Dij1235I,Dij123
     -   4I,EijR,EijI)
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
       P(1) = -p2sq-p3sq+s12+s23
       PaR(1) = B013R+B014R
       PaI(1) = B013I+B014I
       FR(1) = Is12s45*P(1)*PaR(1)
       FI(1) = Is12s45*P(1)*PaI(1)
       F(1)=DCMPLX(FR(1),FI(1))
       FR(2) = Is12s45*PaR(1)
       FI(2) = Is12s45*PaI(1)
       F(2)=DCMPLX(FR(2),FI(2))
       P(2) = p2sq+p3sq-s12-s23
       FR(3) = Is12s45*P(2)*PaR(1)
       FI(3) = Is12s45*P(2)*PaI(1)
       F(3)=DCMPLX(FR(3),FI(3))
       P(3) = p3sq-s12
       P(4) = -(s12*s45)+P(2)*P(3)
       P(5) = p2sq+p3sq-s12-s23+s45
       PaR(2) = Is12*(2*Cij134R(4,2)*P(1)+(B013R-p3sq*Cij134R(2,1))*P(
     -   2)-C0134R*P(4))+Cij134R(1,1)*P(5)
       PaI(2) = Is12*(2*Cij134I(4,2)*P(1)+(B013I-p3sq*Cij134I(2,1))*P(
     -   2)-C0134I*P(4))+Cij134I(1,1)*P(5)
       FR(4) = 4*Is45*PaR(2)
       FI(4) = 4*Is45*PaI(2)
       F(4)=DCMPLX(FR(4),FI(4))
       P(6) = p3sq+s45
       PaR(3) = Cij134R(1,2)-Cij134R(3,2)
       PaI(3) = Cij134I(1,2)-Cij134I(3,2)
       PaR(4) = Cij134R(1,1)+Cij134R(2,1)+Cij134R(3,2)
       PaI(4) = Cij134I(1,1)+Cij134I(2,1)+Cij134I(3,2)
       PaR(5) = -B013R+2*Cij134R(4,2)+Cij134R(2,1)*P(3)+C0134R*P(6)+s1
     -   2*(Cij134R(1,1)+PaR(3))+s45*PaR(4)
       PaI(5) = -B013I+2*Cij134I(4,2)+Cij134I(2,1)*P(3)+C0134I*P(6)+s1
     -   2*(Cij134I(1,1)+PaI(3))+s45*PaI(4)
       FR(5) = -8*Is12s45*PaR(5)
       FI(5) = -8*Is12s45*PaI(5)
       F(5)=DCMPLX(FR(5),FI(5))
       PaR(6) = C0134R+2*Cij134R(1,1)+Cij134R(1,2)-Cij134R(2,1)-Cij134
     -   R(3,2)
       PaI(6) = C0134I+2*Cij134I(1,1)+Cij134I(1,2)-Cij134I(2,1)-Cij134
     -   I(3,2)
       FR(6) = -8*Is45*PaR(6)
       FI(6) = -8*Is45*PaI(6)
       F(6)=DCMPLX(FR(6),FI(6))
       PaR(7) = Cij134R(1,1)+Cij134R(1,2)-Cij134R(2,1)-Cij134R(3,2)
       PaI(7) = Cij134I(1,1)+Cij134I(1,2)-Cij134I(2,1)-Cij134I(3,2)
       PaR(8) = -B013R+p3sq*(C0134R+Cij134R(2,1))+2*Cij134R(4,2)+s12*P
     -   aR(7)
       PaI(8) = -B013I+p3sq*(C0134I+Cij134I(2,1))+2*Cij134I(4,2)+s12*P
     -   aI(7)
       FR(7) = -8*Is12s45*PaR(8)
       FI(7) = -8*Is12s45*PaI(8)
       F(7)=DCMPLX(FR(7),FI(7))
       P(7) = s12-s45
       PaR(9) = C0134R+Cij134R(1,1)+Cij134R(3,2)
       PaI(9) = C0134I+Cij134I(1,1)+Cij134I(3,2)
       PaR(10) = s45*Cij134R(2,1)-Cij134R(2,2)*P(7)+s12*PaR(9)
       PaI(10) = s45*Cij134I(2,1)-Cij134I(2,2)*P(7)+s12*PaI(9)
       FR(8) = Is12s45*PaR(10)
       FI(8) = Is12s45*PaI(10)
       F(8)=DCMPLX(FR(8),FI(8))
       PaR(11) = C0134R+Cij134R(1,1)-Cij134R(2,2)+Cij134R(3,2)
       PaI(11) = C0134I+Cij134I(1,1)-Cij134I(2,2)+Cij134I(3,2)
       FR(9) = Is45*PaR(11)
       FI(9) = Is45*PaI(11)
       F(9)=DCMPLX(FR(9),FI(9))
       P(8) = -p3sq+s12
       PaR(12) = B013R+s12*Cij134R(1,1)-p3sq*Cij134R(2,1)-2*Cij134R(4,
     -   2)+C0134R*P(8)
       PaI(12) = B013I+s12*Cij134I(1,1)-p3sq*Cij134I(2,1)-2*Cij134I(4,
     -   2)+C0134I*P(8)
       FR(10) = 8*Is12s45*PaR(12)
       FI(10) = 8*Is12s45*PaI(12)
       F(10)=DCMPLX(FR(10),FI(10))
       PaR(13) = Is12*(B013R*P(1)+(p3sq*Cij134R(2,1)+2*Cij134R(4,2))*P
     -   (2)+C0134R*P(4))-Cij134R(1,1)*P(5)
       PaI(13) = Is12*(B013I*P(1)+(p3sq*Cij134I(2,1)+2*Cij134I(4,2))*P
     -   (2)+C0134I*P(4))-Cij134I(1,1)*P(5)
       FR(11) = 4*Is45*PaR(13)
       FI(11) = 4*Is45*PaI(13)
       F(11)=DCMPLX(FR(11),FI(11))
       FR(12) = -8*Is12s45*PaR(12)
       FI(12) = -8*Is12s45*PaI(12)
       F(12)=DCMPLX(FR(12),FI(12))
       FR(13) = 8*Is12s45*PaR(5)
       FI(13) = 8*Is12s45*PaI(5)
       F(13)=DCMPLX(FR(13),FI(13))
       FR(14) = 8*Is45*PaR(6)
       FI(14) = 8*Is45*PaI(6)
       F(14)=DCMPLX(FR(14),FI(14))
       FR(15) = 8*Is12s45*PaR(8)
       FI(15) = 8*Is12s45*PaI(8)
       F(15)=DCMPLX(FR(15),FI(15))
       FR(16) = 2*Is45*PaR(13)
       FI(16) = 2*Is45*PaI(13)
       F(16)=DCMPLX(FR(16),FI(16))
       FR(17) = 4*Is12s45*PaR(12)
       FI(17) = 4*Is12s45*PaI(12)
       F(17)=DCMPLX(FR(17),FI(17))
       FR(18) = -4*Is12s45*PaR(5)
       FI(18) = -4*Is12s45*PaI(5)
       F(18)=DCMPLX(FR(18),FI(18))
       FR(19) = -4*Is12s45*PaR(12)
       FI(19) = -4*Is12s45*PaI(12)
       F(19)=DCMPLX(FR(19),FI(19))
       FR(20) = -4*Is45*PaR(6)
       FI(20) = -4*Is45*PaI(6)
       F(20)=DCMPLX(FR(20),FI(20))
       FR(21) = -4*Is12s45*PaR(8)
       FI(21) = -4*Is12s45*PaI(8)
       F(21)=DCMPLX(FR(21),FI(21))
       FR(22) = -2*Is12s45*PaR(12)
       FI(22) = -2*Is12s45*PaI(12)
       F(22)=DCMPLX(FR(22),FI(22))
       P(9) = p2sq-s12
       PaR(14) = C0123R+Cij123R(1,1)
       PaI(14) = C0123I+Cij123I(1,1)
       PaR(15) = B013R+p2sq*Cij123R(2,1)-2*Cij123R(4,2)-P(9)*PaR(14)
       PaI(15) = B013I+p2sq*Cij123I(2,1)-2*Cij123I(4,2)-P(9)*PaI(14)
       FR(23) = 8*Is12s45*PaR(15)
       FI(23) = 8*Is12s45*PaI(15)
       F(23)=DCMPLX(FR(23),FI(23))
       FR(24) = -4*Is12s45*P(1)*PaR(15)
       FI(24) = -4*Is12s45*P(1)*PaI(15)
       F(24)=DCMPLX(FR(24),FI(24))
       FR(25) = 4*Is12s45*P(1)*PaR(15)
       FI(25) = 4*Is12s45*P(1)*PaI(15)
       F(25)=DCMPLX(FR(25),FI(25))
       PaR(16) = -B013R+C0123R*p2sq+p2sq*(Cij123R(1,1)-Cij123R(2,1))+s
     -   12*(Cij123R(2,1)+Cij123R(3,2))+2*Cij123R(4,2)
       PaI(16) = -B013I+C0123I*p2sq+p2sq*(Cij123I(1,1)-Cij123I(2,1))+s
     -   12*(Cij123I(2,1)+Cij123I(3,2))+2*Cij123I(4,2)
       FR(26) = -8*Is12s45*PaR(16)
       FI(26) = -8*Is12s45*PaI(16)
       F(26)=DCMPLX(FR(26),FI(26))
       PaR(17) = C0123R+Cij123R(1,1)+Cij123R(2,1)+Cij123R(3,2)
       PaI(17) = C0123I+Cij123I(1,1)+Cij123I(2,1)+Cij123I(3,2)
       FR(27) = Is45*PaR(17)
       FI(27) = Is45*PaI(17)
       F(27)=DCMPLX(FR(27),FI(27))
       PaR(18) = Cij123R(2,1)+Cij123R(2,2)
       PaI(18) = Cij123I(2,1)+Cij123I(2,2)
       FR(28) = Is45*PaR(18)
       FI(28) = Is45*PaI(18)
       F(28)=DCMPLX(FR(28),FI(28))
       FR(29) = -8*Is12s45*PaR(15)
       FI(29) = -8*Is12s45*PaI(15)
       F(29)=DCMPLX(FR(29),FI(29))
       FR(30) = 8*Is12s45*PaR(16)
       FI(30) = 8*Is12s45*PaI(16)
       F(30)=DCMPLX(FR(30),FI(30))
       FR(31) = 2*Is12s45*P(1)*PaR(15)
       FI(31) = 2*Is12s45*P(1)*PaI(15)
       F(31)=DCMPLX(FR(31),FI(31))
       FR(32) = 4*Is12s45*PaR(15)
       FI(32) = 4*Is12s45*PaI(15)
       F(32)=DCMPLX(FR(32),FI(32))
       FR(33) = -4*Is12s45*PaR(15)
       FI(33) = -4*Is12s45*PaI(15)
       F(33)=DCMPLX(FR(33),FI(33))
       FR(34) = -4*Is12s45*PaR(16)
       FI(34) = -4*Is12s45*PaI(16)
       F(34)=DCMPLX(FR(34),FI(34))
       FR(35) = -2*Is12s45*PaR(15)
       FI(35) = -2*Is12s45*PaI(15)
       F(35)=DCMPLX(FR(35),FI(35))
       PaR(19) = Cij145R(1,1)+Cij145R(1,2)-Cij145R(2,1)+Cij145R(2,2)-2
     -   *Cij145R(3,2)
       PaI(19) = Cij145I(1,1)+Cij145I(1,2)-Cij145I(2,1)+Cij145I(2,2)-2
     -   *Cij145I(3,2)
       FR(36) = Is12*PaR(19)
       FI(36) = Is12*PaI(19)
       F(36)=DCMPLX(FR(36),FI(36))
       PaR(20) = C0145R+Cij145R(1,1)-Cij145R(2,2)+Cij145R(3,2)
       PaI(20) = C0145I+Cij145I(1,1)-Cij145I(2,2)+Cij145I(3,2)
       FR(37) = Is12*PaR(20)
       FI(37) = Is12*PaI(20)
       F(37)=DCMPLX(FR(37),FI(37))
       P(10) = -p4sq+s45
       PaR(21) = B014R+s45*Cij145R(1,1)-p4sq*Cij145R(2,1)-2*Cij145R(4,
     -   2)+C0145R*P(10)
       PaI(21) = B014I+s45*Cij145I(1,1)-p4sq*Cij145I(2,1)-2*Cij145I(4,
     -   2)+C0145I*P(10)
       FR(38) = 8*Is12s45*PaR(21)
       FI(38) = 8*Is12s45*PaI(21)
       F(38)=DCMPLX(FR(38),FI(38))
       P(11) = p4sq-s45
       PaR(22) = Cij145R(1,2)+Cij145R(2,2)-2*Cij145R(3,2)
       PaI(22) = Cij145I(1,2)+Cij145I(2,2)-2*Cij145I(3,2)
       PaR(23) = -B014R+2*Cij145R(4,2)+(C0145R+Cij145R(2,1))*P(11)+s45
     -   *PaR(22)
       PaI(23) = -B014I+2*Cij145I(4,2)+(C0145I+Cij145I(2,1))*P(11)+s45
     -   *PaI(22)
       FR(39) = -8*Is12s45*PaR(23)
       FI(39) = -8*Is12s45*PaI(23)
       F(39)=DCMPLX(FR(39),FI(39))
       FR(40) = 4*Is12s45*P(2)*PaR(21)
       FI(40) = 4*Is12s45*P(2)*PaI(21)
       F(40)=DCMPLX(FR(40),FI(40))
       FR(41) = 4*Is12s45*P(1)*PaR(21)
       FI(41) = 4*Is12s45*P(1)*PaI(21)
       F(41)=DCMPLX(FR(41),FI(41))
       FR(42) = -8*Is12s45*PaR(21)
       FI(42) = -8*Is12s45*PaI(21)
       F(42)=DCMPLX(FR(42),FI(42))
       FR(43) = 8*Is12s45*PaR(23)
       FI(43) = 8*Is12s45*PaI(23)
       F(43)=DCMPLX(FR(43),FI(43))
       FR(44) = 2*Is12s45*P(1)*PaR(21)
       FI(44) = 2*Is12s45*P(1)*PaI(21)
       F(44)=DCMPLX(FR(44),FI(44))
       FR(45) = -4*Is12s45*PaR(23)
       FI(45) = -4*Is12s45*PaI(23)
       F(45)=DCMPLX(FR(45),FI(45))
       FR(46) = 4*Is12s45*PaR(21)
       FI(46) = 4*Is12s45*PaI(21)
       F(46)=DCMPLX(FR(46),FI(46))
       FR(47) = -4*Is12s45*PaR(21)
       FI(47) = -4*Is12s45*PaI(21)
       F(47)=DCMPLX(FR(47),FI(47))
       FR(48) = -2*Is12s45*PaR(21)
       FI(48) = -2*Is12s45*PaI(21)
       F(48)=DCMPLX(FR(48),FI(48))
       P(12) = -s23+s45
       FR(49) = 8*Is45*(C0134R+Cij134R(1,1)+p2sq*(Dij1234R(2,2)-Dij123
     -   4R(4,2))+s12*(Dij1234R(2,1)-Dij1234R(3,1)+Dij1234R(4,2)-Dij1
     -   234R(6,2))+s45*(Dij1234R(3,1)+Dij1234R(6,2))-2*Dij1234R(12,3
     -   )+(D01234R+Dij1234R(1,1))*P(12))
       FI(49) = 8*Is45*(C0134I+Cij134I(1,1)+p2sq*(Dij1234I(2,2)-Dij123
     -   4I(4,2))+s12*(Dij1234I(2,1)-Dij1234I(3,1)+Dij1234I(4,2)-Dij1
     -   234I(6,2))+s45*(Dij1234I(3,1)+Dij1234I(6,2))-2*Dij1234I(12,3
     -   )+(D01234I+Dij1234I(1,1))*P(12))
       F(49)=DCMPLX(FR(49),FI(49))
       P(13) = s23-s45
       FR(50) = 8*Is45*(C0134R+Cij134R(1,1)+p2sq*Dij1234R(2,2)-p2sq*Di
     -   j1234R(4,2)+s12*(Dij1234R(2,1)-Dij1234R(3,1)+Dij1234R(4,2)-D
     -   ij1234R(6,2))-2*Dij1234R(12,3)-(D01234R+Dij1234R(1,1))*P(13)
     -   )
       FI(50) = 8*Is45*(C0134I+Cij134I(1,1)+p2sq*Dij1234I(2,2)-p2sq*Di
     -   j1234I(4,2)+s12*(Dij1234I(2,1)-Dij1234I(3,1)+Dij1234I(4,2)-D
     -   ij1234I(6,2))-2*Dij1234I(12,3)-(D01234I+Dij1234I(1,1))*P(13)
     -   )
       F(50)=DCMPLX(FR(50),FI(50))
       P(14) = p2sq-s12-s23+s45
       P(15) = -2*p2sq+s12
       P(16) = p2sq+p3sq-s23
       FR(51) = 8*Is45*(Cij134R(1,1)-s23*(D01234R+Dij1234R(1,1))+2*(p2
     -   sq*Dij1234R(2,2)+Dij1234R(7,2)-Dij1234R(11,3))+Dij1234R(5,2)
     -   *P(5)-Dij1234R(2,1)*P(9)+Dij1234R(3,1)*P(14)+Dij1234R(4,2)*P
     -   (15)-Dij1234R(6,2)*P(16))
       FI(51) = 8*Is45*(Cij134I(1,1)-s23*(D01234I+Dij1234I(1,1))+2*(p2
     -   sq*Dij1234I(2,2)+Dij1234I(7,2)-Dij1234I(11,3))+Dij1234I(5,2)
     -   *P(5)-Dij1234I(2,1)*P(9)+Dij1234I(3,1)*P(14)+Dij1234I(4,2)*P
     -   (15)-Dij1234I(6,2)*P(16))
       F(51)=DCMPLX(FR(51),FI(51))
       P(17) = -p2sq+s45
       FR(52) = 8*Is45*(C0134R+Cij134R(2,1)+(-s12+s45)*Dij1234R(3,2)+s
     -   12*Dij1234R(5,2)+p2sq*(Dij1234R(2,1)-Dij1234R(5,2)+Dij1234R(
     -   6,2))-2*(Dij1234R(7,2)+Dij1234R(13,3))+Dij1234R(3,1)*P(17))
       FI(52) = 8*Is45*(C0134I+Cij134I(2,1)+(-s12+s45)*Dij1234I(3,2)+s
     -   12*Dij1234I(5,2)+p2sq*(Dij1234I(2,1)-Dij1234I(5,2)+Dij1234I(
     -   6,2))-2*(Dij1234I(7,2)+Dij1234I(13,3))+Dij1234I(3,1)*P(17))
       F(52)=DCMPLX(FR(52),FI(52))
       FR(53) = 8*Is45*(C0134R+Cij134R(2,1)+s12*(-Dij1234R(3,2)+Dij123
     -   4R(5,2))+p2sq*(Dij1234R(2,1)-Dij1234R(3,1)-Dij1234R(5,2)+Dij
     -   1234R(6,2))-2*(Dij1234R(7,2)+Dij1234R(13,3)))
       FI(53) = 8*Is45*(C0134I+Cij134I(2,1)+s12*(-Dij1234I(3,2)+Dij123
     -   4I(5,2))+p2sq*(Dij1234I(2,1)-Dij1234I(3,1)-Dij1234I(5,2)+Dij
     -   1234I(6,2))-2*(Dij1234I(7,2)+Dij1234I(13,3)))
       F(53)=DCMPLX(FR(53),FI(53))
       P(18) = p3sq-s23
       FR(54) = 8*Is45*(C0134R+Cij134R(2,1)-p2sq*Dij1234R(4,2)-4*Dij12
     -   34R(7,2)+2*(p2sq*Dij1234R(6,2)-Dij1234R(11,3))-Dij1234R(3,2)
     -   *P(16)+Dij1234R(5,2)*P(18))
       FI(54) = 8*Is45*(C0134I+Cij134I(2,1)-p2sq*Dij1234I(4,2)-4*Dij12
     -   34I(7,2)+2*(p2sq*Dij1234I(6,2)-Dij1234I(11,3))-Dij1234I(3,2)
     -   *P(16)+Dij1234I(5,2)*P(18))
       F(54)=DCMPLX(FR(54),FI(54))
       P(19) = -p3sq+s23
       P(20) = s23*P(2)+s45*P(19)
       P(21) = p2sq*s45+P(2)*P(9)
       P(22) = 3*p2sq+p3sq-s23
       P(23) = 2*p2sq*s45+P(2)*P(22)
       P(24) = s45+3*P(2)
       FR(55) = -4*Is45*(-(P(2)*(Cij134R(1,1)+2*p2sq*Dij1234R(2,2)-Dij
     -   1234R(4,2)*P(9)+Dij1234R(3,1)*P(14)))+P(5)*(C0134R+Cij134R(2
     -   ,1)-Dij1234R(5,2)*P(9)-Dij1234R(3,2)*P(16))+(D01234R+Dij1234
     -   R(1,1))*P(20)+Dij1234R(2,1)*P(21)+Dij1234R(6,2)*P(23)-2*Dij1
     -   234R(7,2)*P(24))
       FI(55) = -4*Is45*(-(P(2)*(Cij134I(1,1)+2*p2sq*Dij1234I(2,2)-Dij
     -   1234I(4,2)*P(9)+Dij1234I(3,1)*P(14)))+P(5)*(C0134I+Cij134I(2
     -   ,1)-Dij1234I(5,2)*P(9)-Dij1234I(3,2)*P(16))+(D01234I+Dij1234
     -   I(1,1))*P(20)+Dij1234I(2,1)*P(21)+Dij1234I(6,2)*P(23)-2*Dij1
     -   234I(7,2)*P(24))
       F(55)=DCMPLX(FR(55),FI(55))
       FR(56) = -8*Is45*(p3sq*(Dij1234R(2,1)-Dij1234R(3,1)-Dij1234R(6,
     -   2))+2*(Dij1234R(7,2)+Dij1234R(12,3))+Dij1234R(2,2)*P(5)-Dij1
     -   234R(4,2)*P(14))
       FI(56) = -8*Is45*(p3sq*(Dij1234I(2,1)-Dij1234I(3,1)-Dij1234I(6,
     -   2))+2*(Dij1234I(7,2)+Dij1234I(12,3))+Dij1234I(2,2)*P(5)-Dij1
     -   234I(4,2)*P(14))
       F(56)=DCMPLX(FR(56),FI(56))
       P(25) = -p3sq+s45
       P(26) = p3sq-s12-s23
       FR(57) = -8*Is45*(p2sq*Dij1234R(2,2)-p3sq*(Dij1234R(3,1)+Dij123
     -   4R(6,2))+2*(Dij1234R(7,2)+Dij1234R(12,3))-Dij1234R(4,2)*P(14
     -   )-Dij1234R(2,1)*P(25)+Dij1234R(2,2)*P(26))
       FI(57) = -8*Is45*(p2sq*Dij1234I(2,2)-p3sq*(Dij1234I(3,1)+Dij123
     -   4I(6,2))+2*(Dij1234I(7,2)+Dij1234I(12,3))-Dij1234I(4,2)*P(14
     -   )-Dij1234I(2,1)*P(25)+Dij1234I(2,2)*P(26))
       F(57)=DCMPLX(FR(57),FI(57))
       P(27) = p2sq-p3sq
       P(28) = -p2sq+p3sq
       P(29) = p2sq+p3sq
       FR(58) = 8*Is45*(C0134R+s45*(D01234R+2*Dij1234R(1,1)+Dij1234R(1
     -   ,2)-Dij1234R(2,1)-Dij1234R(4,2))-p3sq*Dij1234R(4,2)-p2sq*Dij
     -   1234R(5,2)+s23*(-Dij1234R(1,1)-Dij1234R(1,2)+Dij1234R(2,1)+D
     -   ij1234R(3,1)+Dij1234R(4,2)+Dij1234R(5,2)-Dij1234R(6,2))-6*Di
     -   j1234R(7,2)-2*(p2sq*Dij1234R(2,2)-p2sq*Dij1234R(4,2)+Dij1234
     -   R(11,3))+Dij1234R(2,1)*P(27)+Dij1234R(3,1)*P(28)+Dij1234R(6,
     -   2)*P(29))
       FI(58) = 8*Is45*(C0134I+s45*(D01234I+2*Dij1234I(1,1)+Dij1234I(1
     -   ,2)-Dij1234I(2,1)-Dij1234I(4,2))-p3sq*Dij1234I(4,2)-p2sq*Dij
     -   1234I(5,2)+s23*(-Dij1234I(1,1)-Dij1234I(1,2)+Dij1234I(2,1)+D
     -   ij1234I(3,1)+Dij1234I(4,2)+Dij1234I(5,2)-Dij1234I(6,2))-6*Di
     -   j1234I(7,2)-2*(p2sq*Dij1234I(2,2)-p2sq*Dij1234I(4,2)+Dij1234
     -   I(11,3))+Dij1234I(2,1)*P(27)+Dij1234I(3,1)*P(28)+Dij1234I(6,
     -   2)*P(29))
       F(58)=DCMPLX(FR(58),FI(58))
       P(30) = p2sq-s23+s45
       FR(59) = -8*Is45*(C0134R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,2)-
     -   2*(Dij1234R(7,2)-Dij1234R(13,3))+Dij1234R(6,2)*P(2)-Dij1234R
     -   (5,2)*P(14)-Dij1234R(3,1)*P(30))
       FI(59) = -8*Is45*(C0134I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,2)-
     -   2*(Dij1234I(7,2)-Dij1234I(13,3))+Dij1234I(6,2)*P(2)-Dij1234I
     -   (5,2)*P(14)-Dij1234I(3,1)*P(30))
       F(59)=DCMPLX(FR(59),FI(59))
       P(31) = -p2sq+s23
       FR(60) = -8*Is45*(C0134R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,2)-
     -   2*(Dij1234R(7,2)-Dij1234R(13,3))+Dij1234R(6,2)*P(5)-Dij1234R
     -   (5,2)*P(14)+Dij1234R(3,1)*P(31))
       FI(60) = -8*Is45*(C0134I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,2)-
     -   2*(Dij1234I(7,2)-Dij1234I(13,3))+Dij1234I(6,2)*P(5)-Dij1234I
     -   (5,2)*P(14)+Dij1234I(3,1)*P(31))
       F(60)=DCMPLX(FR(60),FI(60))
       FR(61) = 8*Is45*(Cij134R(1,1)-Cij134R(2,1)+s45*Dij1234R(1,2)+p3
     -   sq*Dij1234R(3,2)-p3sq*Dij1234R(4,2)+s12*Dij1234R(4,2)-s12*Di
     -   j1234R(5,2)+s23*(-Dij1234R(1,2)+Dij1234R(2,1)-Dij1234R(3,2)+
     -   Dij1234R(4,2)+Dij1234R(5,2))+p2sq*(Dij1234R(3,2)+Dij1234R(4,
     -   2)-2*Dij1234R(6,2))+Dij1234R(3,1)*P(3)+Dij1234R(2,1)*P(8)+D0
     -   1234R*P(12)-2*(Dij1234R(11,3)+Dij1234R(1,1)*P(13)))
       FI(61) = 8*Is45*(Cij134I(1,1)-Cij134I(2,1)+s45*Dij1234I(1,2)+p3
     -   sq*Dij1234I(3,2)-p3sq*Dij1234I(4,2)+s12*Dij1234I(4,2)-s12*Di
     -   j1234I(5,2)+s23*(-Dij1234I(1,2)+Dij1234I(2,1)-Dij1234I(3,2)+
     -   Dij1234I(4,2)+Dij1234I(5,2))+p2sq*(Dij1234I(3,2)+Dij1234I(4,
     -   2)-2*Dij1234I(6,2))+Dij1234I(3,1)*P(3)+Dij1234I(2,1)*P(8)+D0
     -   1234I*P(12)-2*(Dij1234I(11,3)+Dij1234I(1,1)*P(13)))
       F(61)=DCMPLX(FR(61),FI(61))
       FR(62) = 4*Is45*(Cij134R(1,1)*P(1)+P(2)*(-2*p2sq*Dij1234R(2,2)+
     -   Dij1234R(4,2)*P(9)-Dij1234R(3,1)*P(14))+P(5)*(C0134R+Cij134R
     -   (2,1)-Dij1234R(5,2)*P(9)-Dij1234R(3,2)*P(16))+(D01234R+Dij12
     -   34R(1,1))*P(20)+Dij1234R(2,1)*P(21)+Dij1234R(6,2)*P(23)-2*Di
     -   j1234R(7,2)*P(24))
       FI(62) = 4*Is45*(Cij134I(1,1)*P(1)+P(2)*(-2*p2sq*Dij1234I(2,2)+
     -   Dij1234I(4,2)*P(9)-Dij1234I(3,1)*P(14))+P(5)*(C0134I+Cij134I
     -   (2,1)-Dij1234I(5,2)*P(9)-Dij1234I(3,2)*P(16))+(D01234I+Dij12
     -   34I(1,1))*P(20)+Dij1234I(2,1)*P(21)+Dij1234I(6,2)*P(23)-2*Di
     -   j1234I(7,2)*P(24))
       F(62)=DCMPLX(FR(62),FI(62))
       P(32) = s12+s23-s45
       P(33) = -(p2sq*s23)+p3sq*P(12)+s23*P(32)
       FR(63) = 4*Is45*(2*((Dij1234R(7,2)+Dij1234R(12,3))*P(2)-Dij1234
     -   R(13,3)*P(5))+(D01234R+Dij1234R(1,1))*P(33))
       FI(63) = 4*Is45*(2*((Dij1234I(7,2)+Dij1234I(12,3))*P(2)-Dij1234
     -   I(13,3)*P(5))+(D01234I+Dij1234I(1,1))*P(33))
       F(63)=DCMPLX(FR(63),FI(63))
       P(34) = 2*p2sq+p3sq-s12-s23
       P(35) = p2sq-s23
       P(36) = p3sq-s12+s45+2*P(35)
       P(37) = 3*p2sq+2*P(26)
       P(38) = p3sq-s12+s45
       P(39) = 3*P(35)+2*P(38)
       FR(64) = 8*Is45*(C0134R+Cij134R(1,1)-s23*(D01234R+2*Dij1234R(1,
     -   1)+Dij1234R(1,2))-4*(Dij1234R(7,2)+Dij1234R(11,3))+Dij1234R(
     -   4,3)*P(2)-Dij1234R(5,3)*P(5)+Dij1234R(2,1)*P(34)-Dij1234R(3,
     -   1)*P(36)+Dij1234R(4,2)*P(37)-Dij1234R(5,2)*P(39))
       FI(64) = 8*Is45*(C0134I+Cij134I(1,1)-s23*(D01234I+2*Dij1234I(1,
     -   1)+Dij1234I(1,2))-4*(Dij1234I(7,2)+Dij1234I(11,3))+Dij1234I(
     -   4,3)*P(2)-Dij1234I(5,3)*P(5)+Dij1234I(2,1)*P(34)-Dij1234I(3,
     -   1)*P(36)+Dij1234I(4,2)*P(37)-Dij1234I(5,2)*P(39))
       F(64)=DCMPLX(FR(64),FI(64))
       P(40) = s12-2*P(16)
       P(41) = s12-s45-2*P(16)
       P(42) = p2sq-s12+2*P(18)
       FR(65) = 8*Is45*(C0134R+Cij134R(1,1)-s23*(D01234R+Dij1234R(1,1)
     -   )-2*Dij1234R(12,3)+Dij1234R(6,3)*P(2)-(Dij1234R(5,2)+Dij1234
     -   R(10,3))*P(5)+Dij1234R(2,2)*P(34)-Dij1234R(2,1)*P(40)+(Dij12
     -   34R(3,1)+Dij1234R(6,2))*P(41)+Dij1234R(4,2)*P(42))
       FI(65) = 8*Is45*(C0134I+Cij134I(1,1)-s23*(D01234I+Dij1234I(1,1)
     -   )-2*Dij1234I(12,3)+Dij1234I(6,3)*P(2)-(Dij1234I(5,2)+Dij1234
     -   I(10,3))*P(5)+Dij1234I(2,2)*P(34)-Dij1234I(2,1)*P(40)+(Dij12
     -   34I(3,1)+Dij1234I(6,2))*P(41)+Dij1234I(4,2)*P(42))
       F(65)=DCMPLX(FR(65),FI(65))
       P(43) = p2sq-s12-s23
       FR(66) = 8*Is45*(-2*(Dij1234R(7,2)+Dij1234R(12,3))+(Dij1234R(2,
     -   2)+Dij1234R(6,3))*P(2)-(Dij1234R(5,2)+Dij1234R(10,3))*P(5)-(
     -   Dij1234R(3,1)+Dij1234R(6,2))*P(14)+(Dij1234R(2,1)+Dij1234R(4
     -   ,2))*P(43))
       FI(66) = 8*Is45*(-2*(Dij1234I(7,2)+Dij1234I(12,3))+(Dij1234I(2,
     -   2)+Dij1234I(6,3))*P(2)-(Dij1234I(5,2)+Dij1234I(10,3))*P(5)-(
     -   Dij1234I(3,1)+Dij1234I(6,2))*P(14)+(Dij1234I(2,1)+Dij1234I(4
     -   ,2))*P(43))
       F(66)=DCMPLX(FR(66),FI(66))
       FR(67) = 8*Is45*((Dij1234R(2,1)+2*Dij1234R(2,2)+Dij1234R(2,3))*
     -   P(2)-(Dij1234R(3,1)+2*Dij1234R(6,2)+Dij1234R(8,3))*P(5))
       FI(67) = 8*Is45*((Dij1234I(2,1)+2*Dij1234I(2,2)+Dij1234I(2,3))*
     -   P(2)-(Dij1234I(3,1)+2*Dij1234I(6,2)+Dij1234I(8,3))*P(5))
       F(67)=DCMPLX(FR(67),FI(67))
       P(44) = 2*p2sq+p3sq-s12-s23+s45
       P(45) = p2sq-s12+s45
       FR(68) = 8*Is45*(C0134R+Cij134R(2,1)-2*(Dij1234R(7,2)+Dij1234R(
     -   13,3))+(Dij1234R(4,2)+Dij1234R(10,3))*P(2)-Dij1234R(7,3)*P(5
     -   )+(Dij1234R(2,1)+Dij1234R(6,2))*P(34)+Dij1234R(3,2)*P(41)-Di
     -   j1234R(3,1)*P(44)-Dij1234R(5,2)*P(45))
       FI(68) = 8*Is45*(C0134I+Cij134I(2,1)-2*(Dij1234I(7,2)+Dij1234I(
     -   13,3))+(Dij1234I(4,2)+Dij1234I(10,3))*P(2)-Dij1234I(7,3)*P(5
     -   )+(Dij1234I(2,1)+Dij1234I(6,2))*P(34)+Dij1234I(3,2)*P(41)-Di
     -   j1234I(3,1)*P(44)-Dij1234I(5,2)*P(45))
       F(68)=DCMPLX(FR(68),FI(68))
       FR(69) = -8*Is45*(C0134R+p2sq*Dij1234R(2,1)+p3sq*Dij1234R(5,2)-
     -   2*(Dij1234R(7,2)-Dij1234R(13,3))-(Dij1234R(6,2)+Dij1234R(10,
     -   3))*P(2)+Dij1234R(7,3)*P(5)+Dij1234R(3,2)*P(14)+Dij1234R(3,1
     -   )*P(31))
       FI(69) = -8*Is45*(C0134I+p2sq*Dij1234I(2,1)+p3sq*Dij1234I(5,2)-
     -   2*(Dij1234I(7,2)-Dij1234I(13,3))-(Dij1234I(6,2)+Dij1234I(10,
     -   3))*P(2)+Dij1234I(7,3)*P(5)+Dij1234I(3,2)*P(14)+Dij1234I(3,1
     -   )*P(31))
       F(69)=DCMPLX(FR(69),FI(69))
       FR(70) = -Dij1234R(3,1)-Dij1234R(3,2)-Dij1234R(6,2)-Dij1234R(9,
     -   3)+Is45*(Dij1234R(2,1)+Dij1234R(2,2)-Dij1234R(3,1)-Dij1234R(
     -   3,2)+Dij1234R(8,3)-Dij1234R(9,3))*P(2)
       FI(70) = -Dij1234I(3,1)-Dij1234I(3,2)-Dij1234I(6,2)-Dij1234I(9,
     -   3)+Is45*(Dij1234I(2,1)+Dij1234I(2,2)-Dij1234I(3,1)-Dij1234I(
     -   3,2)+Dij1234I(8,3)-Dij1234I(9,3))*P(2)
       F(70)=DCMPLX(FR(70),FI(70))
       FR(71) = Is45*((Dij1234R(6,2)+Dij1234R(8,3))*P(2)-(Dij1234R(3,2
     -   )+Dij1234R(9,3))*P(5))
       FI(71) = Is45*((Dij1234I(6,2)+Dij1234I(8,3))*P(2)-(Dij1234I(3,2
     -   )+Dij1234I(9,3))*P(5))
       F(71)=DCMPLX(FR(71),FI(71))
       FR(72) = Is45*((Dij1234R(6,2)+Dij1234R(9,3))*P(2)-(Dij1234R(3,2
     -   )+Dij1234R(3,3))*P(5))
       FI(72) = Is45*((Dij1234I(6,2)+Dij1234I(9,3))*P(2)-(Dij1234I(3,2
     -   )+Dij1234I(3,3))*P(5))
       F(72)=DCMPLX(FR(72),FI(72))
       FR(73) = 8*Is45*(-2*(Dij1234R(7,2)+Dij1234R(12,3)-Dij1234R(13,3
     -   ))+(D01234R+Dij1234R(1,1))*P(13))
       FI(73) = 8*Is45*(-2*(Dij1234I(7,2)+Dij1234I(12,3)-Dij1234I(13,3
     -   ))+(D01234I+Dij1234I(1,1))*P(13))
       F(73)=DCMPLX(FR(73),FI(73))
       FR(74) = -8*Is45*(C0134R+Cij134R(1,1)+p2sq*(Dij1234R(2,2)-Dij12
     -   34R(4,2))+s12*(Dij1234R(2,1)-Dij1234R(3,1)+Dij1234R(4,2)-Dij
     -   1234R(6,2))+s45*(Dij1234R(3,1)+Dij1234R(6,2))-2*Dij1234R(12,
     -   3)+(D01234R+Dij1234R(1,1))*P(12))
       FI(74) = -8*Is45*(C0134I+Cij134I(1,1)+p2sq*(Dij1234I(2,2)-Dij12
     -   34I(4,2))+s12*(Dij1234I(2,1)-Dij1234I(3,1)+Dij1234I(4,2)-Dij
     -   1234I(6,2))+s45*(Dij1234I(3,1)+Dij1234I(6,2))-2*Dij1234I(12,
     -   3)+(D01234I+Dij1234I(1,1))*P(12))
       F(74)=DCMPLX(FR(74),FI(74))
       P(46) = 2*p2sq-s12
       FR(75) = 8*Is45*(-Cij134R(1,1)+s23*(D01234R+Dij1234R(1,1))-2*(p
     -   2sq*Dij1234R(2,2)+Dij1234R(7,2)-Dij1234R(11,3))-Dij1234R(5,2
     -   )*P(5)+Dij1234R(2,1)*P(9)-Dij1234R(3,1)*P(14)+Dij1234R(6,2)*
     -   P(16)+Dij1234R(4,2)*P(46))
       FI(75) = 8*Is45*(-Cij134I(1,1)+s23*(D01234I+Dij1234I(1,1))-2*(p
     -   2sq*Dij1234I(2,2)+Dij1234I(7,2)-Dij1234I(11,3))-Dij1234I(5,2
     -   )*P(5)+Dij1234I(2,1)*P(9)-Dij1234I(3,1)*P(14)+Dij1234I(6,2)*
     -   P(16)+Dij1234I(4,2)*P(46))
       F(75)=DCMPLX(FR(75),FI(75))
       FR(76) = -8*Is45*(C0134R+Cij134R(2,1)+(-s12+s45)*Dij1234R(3,2)+
     -   s12*Dij1234R(5,2)+p2sq*(Dij1234R(2,1)-Dij1234R(5,2)+Dij1234R
     -   (6,2))-2*(Dij1234R(7,2)+Dij1234R(13,3))+Dij1234R(3,1)*P(17))
       FI(76) = -8*Is45*(C0134I+Cij134I(2,1)+(-s12+s45)*Dij1234I(3,2)+
     -   s12*Dij1234I(5,2)+p2sq*(Dij1234I(2,1)-Dij1234I(5,2)+Dij1234I
     -   (6,2))-2*(Dij1234I(7,2)+Dij1234I(13,3))+Dij1234I(3,1)*P(17))
       F(76)=DCMPLX(FR(76),FI(76))
       FR(77) = 8*Is45*(p3sq*(Dij1234R(2,1)-Dij1234R(3,1)-Dij1234R(6,2
     -   ))+2*(Dij1234R(7,2)+Dij1234R(12,3))+Dij1234R(2,2)*P(5)-Dij12
     -   34R(4,2)*P(14))
       FI(77) = 8*Is45*(p3sq*(Dij1234I(2,1)-Dij1234I(3,1)-Dij1234I(6,2
     -   ))+2*(Dij1234I(7,2)+Dij1234I(12,3))+Dij1234I(2,2)*P(5)-Dij12
     -   34I(4,2)*P(14))
       F(77)=DCMPLX(FR(77),FI(77))
       FR(78) = -8*Is45*(C0134R+s45*(D01234R+2*Dij1234R(1,1)+Dij1234R(
     -   1,2)-Dij1234R(2,1)-Dij1234R(4,2))-p3sq*Dij1234R(4,2)-p2sq*Di
     -   j1234R(5,2)+s23*(-Dij1234R(1,1)-Dij1234R(1,2)+Dij1234R(2,1)+
     -   Dij1234R(3,1)+Dij1234R(4,2)+Dij1234R(5,2)-Dij1234R(6,2))-6*D
     -   ij1234R(7,2)-2*(p2sq*Dij1234R(2,2)-p2sq*Dij1234R(4,2)+Dij123
     -   4R(11,3))+Dij1234R(2,1)*P(27)+Dij1234R(3,1)*P(28)+Dij1234R(6
     -   ,2)*P(29))
       FI(78) = -8*Is45*(C0134I+s45*(D01234I+2*Dij1234I(1,1)+Dij1234I(
     -   1,2)-Dij1234I(2,1)-Dij1234I(4,2))-p3sq*Dij1234I(4,2)-p2sq*Di
     -   j1234I(5,2)+s23*(-Dij1234I(1,1)-Dij1234I(1,2)+Dij1234I(2,1)+
     -   Dij1234I(3,1)+Dij1234I(4,2)+Dij1234I(5,2)-Dij1234I(6,2))-6*D
     -   ij1234I(7,2)-2*(p2sq*Dij1234I(2,2)-p2sq*Dij1234I(4,2)+Dij123
     -   4I(11,3))+Dij1234I(2,1)*P(27)+Dij1234I(3,1)*P(28)+Dij1234I(6
     -   ,2)*P(29))
       F(78)=DCMPLX(FR(78),FI(78))
       FR(79) = 8*Is45*(C0134R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,2)-2
     -   *(Dij1234R(7,2)-Dij1234R(13,3))+Dij1234R(6,2)*P(5)-Dij1234R(
     -   5,2)*P(14)+Dij1234R(3,1)*P(31))
       FI(79) = 8*Is45*(C0134I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,2)-2
     -   *(Dij1234I(7,2)-Dij1234I(13,3))+Dij1234I(6,2)*P(5)-Dij1234I(
     -   5,2)*P(14)+Dij1234I(3,1)*P(31))
       F(79)=DCMPLX(FR(79),FI(79))
       FR(80) = -16*Is45*(2*(Dij1234R(2,1)-Dij1234R(3,1))+Dij1234R(4,3
     -   )+3*(Dij1234R(4,2)-Dij1234R(5,2))-Dij1234R(5,3))
       FI(80) = -16*Is45*(2*(Dij1234I(2,1)-Dij1234I(3,1))+Dij1234I(4,3
     -   )+3*(Dij1234I(4,2)-Dij1234I(5,2))-Dij1234I(5,3))
       F(80)=DCMPLX(FR(80),FI(80))
       FR(81) = -16*Is45*(Dij1234R(2,2)-Dij1234R(5,2)+2*(Dij1234R(2,1)
     -   -Dij1234R(3,1)+Dij1234R(4,2)-Dij1234R(6,2))+Dij1234R(6,3)-Di
     -   j1234R(10,3))
       FI(81) = -16*Is45*(Dij1234I(2,2)-Dij1234I(5,2)+2*(Dij1234I(2,1)
     -   -Dij1234I(3,1)+Dij1234I(4,2)-Dij1234I(6,2))+Dij1234I(6,3)-Di
     -   j1234I(10,3))
       F(81)=DCMPLX(FR(81),FI(81))
       FR(82) = -16*Is45*(Dij1234R(2,1)+2*Dij1234R(2,2)-Dij1234R(3,1)-
     -   Dij1234R(5,2)-Dij1234R(6,2)+Dij1234R(6,3)-Dij1234R(10,3))
       FI(82) = -16*Is45*(Dij1234I(2,1)+2*Dij1234I(2,2)-Dij1234I(3,1)-
     -   Dij1234I(5,2)-Dij1234I(6,2)+Dij1234I(6,3)-Dij1234I(10,3))
       F(82)=DCMPLX(FR(82),FI(82))
       FR(83) = -16*Is45*(Dij1234R(2,1)+Dij1234R(2,3)-Dij1234R(3,1)+2*
     -   (Dij1234R(2,2)-Dij1234R(6,2))-Dij1234R(8,3))
       FI(83) = -16*Is45*(Dij1234I(2,1)+Dij1234I(2,3)-Dij1234I(3,1)+2*
     -   (Dij1234I(2,2)-Dij1234I(6,2))-Dij1234I(8,3))
       F(83)=DCMPLX(FR(83),FI(83))
       FR(84) = -16*Is45*(Dij1234R(2,1)-Dij1234R(3,1)-2*Dij1234R(3,2)+
     -   Dij1234R(4,2)+Dij1234R(6,2)-Dij1234R(7,3)+Dij1234R(10,3))
       FI(84) = -16*Is45*(Dij1234I(2,1)-Dij1234I(3,1)-2*Dij1234I(3,2)+
     -   Dij1234I(4,2)+Dij1234I(6,2)-Dij1234I(7,3)+Dij1234I(10,3))
       F(84)=DCMPLX(FR(84),FI(84))
       FR(85) = 16*Is45*(Dij1234R(3,2)+Dij1234R(5,2)-2*Dij1234R(6,2)+D
     -   ij1234R(7,3)-Dij1234R(10,3))
       FI(85) = 16*Is45*(Dij1234I(3,2)+Dij1234I(5,2)-2*Dij1234I(6,2)+D
     -   ij1234I(7,3)-Dij1234I(10,3))
       F(85)=DCMPLX(FR(85),FI(85))
       FR(86) = Is45*(Dij1234R(2,1)+Dij1234R(2,2)-Dij1234R(3,1)-Dij123
     -   4R(3,2)+Dij1234R(8,3)-Dij1234R(9,3))
       FI(86) = Is45*(Dij1234I(2,1)+Dij1234I(2,2)-Dij1234I(3,1)-Dij123
     -   4I(3,2)+Dij1234I(8,3)-Dij1234I(9,3))
       F(86)=DCMPLX(FR(86),FI(86))
       FR(87) = Is45*(-Dij1234R(3,2)+Dij1234R(6,2)+Dij1234R(8,3)-Dij12
     -   34R(9,3))
       FI(87) = Is45*(-Dij1234I(3,2)+Dij1234I(6,2)+Dij1234I(8,3)-Dij12
     -   34I(9,3))
       F(87)=DCMPLX(FR(87),FI(87))
       FR(88) = Is45*(-Dij1234R(3,2)-Dij1234R(3,3)+Dij1234R(6,2)+Dij12
     -   34R(9,3))
       FI(88) = Is45*(-Dij1234I(3,2)-Dij1234I(3,3)+Dij1234I(6,2)+Dij12
     -   34I(9,3))
       F(88)=DCMPLX(FR(88),FI(88))
       FR(89) = Is45*(16*(Dij1234R(7,2)+Dij1234R(12,3)-Dij1234R(13,3))
     -   -8*(D01234R+Dij1234R(1,1))*P(13))
       FI(89) = Is45*(16*(Dij1234I(7,2)+Dij1234I(12,3)-Dij1234I(13,3))
     -   -8*(D01234I+Dij1234I(1,1))*P(13))
       F(89)=DCMPLX(FR(89),FI(89))
       FR(90) = -8*Is45*(C0134R+Cij134R(1,1)+p2sq*Dij1234R(2,2)-p2sq*D
     -   ij1234R(4,2)+s12*(Dij1234R(2,1)-Dij1234R(3,1)+Dij1234R(4,2)-
     -   Dij1234R(6,2))-2*Dij1234R(12,3)-(D01234R+Dij1234R(1,1))*P(13
     -   ))
       FI(90) = -8*Is45*(C0134I+Cij134I(1,1)+p2sq*Dij1234I(2,2)-p2sq*D
     -   ij1234I(4,2)+s12*(Dij1234I(2,1)-Dij1234I(3,1)+Dij1234I(4,2)-
     -   Dij1234I(6,2))-2*Dij1234I(12,3)-(D01234I+Dij1234I(1,1))*P(13
     -   ))
       F(90)=DCMPLX(FR(90),FI(90))
       FR(91) = -8*Is45*(C0134R+Cij134R(2,1)+s12*(-Dij1234R(3,2)+Dij12
     -   34R(5,2))+p2sq*(Dij1234R(2,1)-Dij1234R(3,1)-Dij1234R(5,2)+Di
     -   j1234R(6,2))-2*(Dij1234R(7,2)+Dij1234R(13,3)))
       FI(91) = -8*Is45*(C0134I+Cij134I(2,1)+s12*(-Dij1234I(3,2)+Dij12
     -   34I(5,2))+p2sq*(Dij1234I(2,1)-Dij1234I(3,1)-Dij1234I(5,2)+Di
     -   j1234I(6,2))-2*(Dij1234I(7,2)+Dij1234I(13,3)))
       F(91)=DCMPLX(FR(91),FI(91))
       FR(92) = -8*Is45*(C0134R+Cij134R(2,1)-p2sq*Dij1234R(4,2)-4*Dij1
     -   234R(7,2)+2*(p2sq*Dij1234R(6,2)-Dij1234R(11,3))-Dij1234R(3,2
     -   )*P(16)-Dij1234R(5,2)*P(19))
       FI(92) = -8*Is45*(C0134I+Cij134I(2,1)-p2sq*Dij1234I(4,2)-4*Dij1
     -   234I(7,2)+2*(p2sq*Dij1234I(6,2)-Dij1234I(11,3))-Dij1234I(3,2
     -   )*P(16)-Dij1234I(5,2)*P(19))
       F(92)=DCMPLX(FR(92),FI(92))
       FR(93) = 8*Is45*(p2sq*Dij1234R(2,2)-p3sq*(Dij1234R(3,1)+Dij1234
     -   R(6,2))+2*(Dij1234R(7,2)+Dij1234R(12,3))-Dij1234R(4,2)*P(14)
     -   -Dij1234R(2,1)*P(25)+Dij1234R(2,2)*P(26))
       FI(93) = 8*Is45*(p2sq*Dij1234I(2,2)-p3sq*(Dij1234I(3,1)+Dij1234
     -   I(6,2))+2*(Dij1234I(7,2)+Dij1234I(12,3))-Dij1234I(4,2)*P(14)
     -   -Dij1234I(2,1)*P(25)+Dij1234I(2,2)*P(26))
       F(93)=DCMPLX(FR(93),FI(93))
       FR(94) = 8*Is45*(C0134R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,2)-2
     -   *(Dij1234R(7,2)-Dij1234R(13,3))+Dij1234R(6,2)*P(2)-Dij1234R(
     -   5,2)*P(14)-Dij1234R(3,1)*P(30))
       FI(94) = 8*Is45*(C0134I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,2)-2
     -   *(Dij1234I(7,2)-Dij1234I(13,3))+Dij1234I(6,2)*P(2)-Dij1234I(
     -   5,2)*P(14)-Dij1234I(3,1)*P(30))
       F(94)=DCMPLX(FR(94),FI(94))
       FR(95) = -8*Is45*(Cij134R(1,1)-Cij134R(2,1)+s45*Dij1234R(1,2)+p
     -   3sq*Dij1234R(3,2)-p3sq*Dij1234R(4,2)+s12*Dij1234R(4,2)-s12*D
     -   ij1234R(5,2)+s23*(-Dij1234R(1,2)+Dij1234R(2,1)-Dij1234R(3,2)
     -   +Dij1234R(4,2)+Dij1234R(5,2))+p2sq*(Dij1234R(3,2)+Dij1234R(4
     -   ,2)-2*Dij1234R(6,2))+Dij1234R(3,1)*P(3)+Dij1234R(2,1)*P(8)+D
     -   01234R*P(12)-2*(Dij1234R(11,3)+Dij1234R(1,1)*P(13)))
       FI(95) = -8*Is45*(Cij134I(1,1)-Cij134I(2,1)+s45*Dij1234I(1,2)+p
     -   3sq*Dij1234I(3,2)-p3sq*Dij1234I(4,2)+s12*Dij1234I(4,2)-s12*D
     -   ij1234I(5,2)+s23*(-Dij1234I(1,2)+Dij1234I(2,1)-Dij1234I(3,2)
     -   +Dij1234I(4,2)+Dij1234I(5,2))+p2sq*(Dij1234I(3,2)+Dij1234I(4
     -   ,2)-2*Dij1234I(6,2))+Dij1234I(3,1)*P(3)+Dij1234I(2,1)*P(8)+D
     -   01234I*P(12)-2*(Dij1234I(11,3)+Dij1234I(1,1)*P(13)))
       F(95)=DCMPLX(FR(95),FI(95))
       FR(96) = 16*Is45*(2*(Dij1234R(2,1)-Dij1234R(3,1))+Dij1234R(4,3)
     -   +3*(Dij1234R(4,2)-Dij1234R(5,2))-Dij1234R(5,3))
       FI(96) = 16*Is45*(2*(Dij1234I(2,1)-Dij1234I(3,1))+Dij1234I(4,3)
     -   +3*(Dij1234I(4,2)-Dij1234I(5,2))-Dij1234I(5,3))
       F(96)=DCMPLX(FR(96),FI(96))
       FR(97) = 16*Is45*(Dij1234R(2,2)-Dij1234R(5,2)+2*(Dij1234R(2,1)-
     -   Dij1234R(3,1)+Dij1234R(4,2)-Dij1234R(6,2))+Dij1234R(6,3)-Dij
     -   1234R(10,3))
       FI(97) = 16*Is45*(Dij1234I(2,2)-Dij1234I(5,2)+2*(Dij1234I(2,1)-
     -   Dij1234I(3,1)+Dij1234I(4,2)-Dij1234I(6,2))+Dij1234I(6,3)-Dij
     -   1234I(10,3))
       F(97)=DCMPLX(FR(97),FI(97))
       FR(98) = 16*Is45*(Dij1234R(2,1)+2*Dij1234R(2,2)-Dij1234R(3,1)-D
     -   ij1234R(5,2)-Dij1234R(6,2)+Dij1234R(6,3)-Dij1234R(10,3))
       FI(98) = 16*Is45*(Dij1234I(2,1)+2*Dij1234I(2,2)-Dij1234I(3,1)-D
     -   ij1234I(5,2)-Dij1234I(6,2)+Dij1234I(6,3)-Dij1234I(10,3))
       F(98)=DCMPLX(FR(98),FI(98))
       FR(99) = 16*Is45*(Dij1234R(2,1)+Dij1234R(2,3)-Dij1234R(3,1)+2*(
     -   Dij1234R(2,2)-Dij1234R(6,2))-Dij1234R(8,3))
       FI(99) = 16*Is45*(Dij1234I(2,1)+Dij1234I(2,3)-Dij1234I(3,1)+2*(
     -   Dij1234I(2,2)-Dij1234I(6,2))-Dij1234I(8,3))
       F(99)=DCMPLX(FR(99),FI(99))
       FR(100) = 16*Is45*(Dij1234R(2,1)-Dij1234R(3,1)-2*Dij1234R(3,2)+
     -   Dij1234R(4,2)+Dij1234R(6,2)-Dij1234R(7,3)+Dij1234R(10,3))
       FI(100) = 16*Is45*(Dij1234I(2,1)-Dij1234I(3,1)-2*Dij1234I(3,2)+
     -   Dij1234I(4,2)+Dij1234I(6,2)-Dij1234I(7,3)+Dij1234I(10,3))
       F(100)=DCMPLX(FR(100),FI(100))
       FR(101) = -16*Is45*(Dij1234R(3,2)+Dij1234R(5,2)-2*Dij1234R(6,2)
     -   +Dij1234R(7,3)-Dij1234R(10,3))
       FI(101) = -16*Is45*(Dij1234I(3,2)+Dij1234I(5,2)-2*Dij1234I(6,2)
     -   +Dij1234I(7,3)-Dij1234I(10,3))
       F(101)=DCMPLX(FR(101),FI(101))
       FR(102) = 2*Is45*(Cij134R(1,1)*P(1)+P(2)*(-2*p2sq*Dij1234R(2,2)
     -   +Dij1234R(4,2)*P(9)-Dij1234R(3,1)*P(14))+P(5)*(C0134R+Cij134
     -   R(2,1)-Dij1234R(5,2)*P(9)-Dij1234R(3,2)*P(16))+(D01234R+Dij1
     -   234R(1,1))*P(20)+Dij1234R(2,1)*P(21)+Dij1234R(6,2)*P(23)-2*D
     -   ij1234R(7,2)*P(24))
       FI(102) = 2*Is45*(Cij134I(1,1)*P(1)+P(2)*(-2*p2sq*Dij1234I(2,2)
     -   +Dij1234I(4,2)*P(9)-Dij1234I(3,1)*P(14))+P(5)*(C0134I+Cij134
     -   I(2,1)-Dij1234I(5,2)*P(9)-Dij1234I(3,2)*P(16))+(D01234I+Dij1
     -   234I(1,1))*P(20)+Dij1234I(2,1)*P(21)+Dij1234I(6,2)*P(23)-2*D
     -   ij1234I(7,2)*P(24))
       F(102)=DCMPLX(FR(102),FI(102))
       FR(103) = -4*Is45*(C0134R-Cij134R(1,1)+Cij134R(2,1)-2*p2sq*Dij1
     -   234R(2,2)+s23*Dij1234R(3,1)-6*Dij1234R(7,2)+(Dij1234R(2,1)-D
     -   ij1234R(3,1)+Dij1234R(4,2)-Dij1234R(5,2))*P(9)+(D01234R+Dij1
     -   234R(1,1))*P(13)-Dij1234R(3,2)*P(16)+Dij1234R(6,2)*P(22))
       FI(103) = -4*Is45*(C0134I-Cij134I(1,1)+Cij134I(2,1)-2*p2sq*Dij1
     -   234I(2,2)+s23*Dij1234I(3,1)-6*Dij1234I(7,2)+(Dij1234I(2,1)-D
     -   ij1234I(3,1)+Dij1234I(4,2)-Dij1234I(5,2))*P(9)+(D01234I+Dij1
     -   234I(1,1))*P(13)-Dij1234I(3,2)*P(16)+Dij1234I(6,2)*P(22))
       F(103)=DCMPLX(FR(103),FI(103))
       FR(104) = 4*Is45*(C0134R+Cij134R(1,1)+p2sq*(Dij1234R(2,2)-Dij12
     -   34R(4,2))+s12*(Dij1234R(2,1)-Dij1234R(3,1)+Dij1234R(4,2)-Dij
     -   1234R(6,2))+s45*(Dij1234R(3,1)+Dij1234R(6,2))-2*Dij1234R(12,
     -   3)+(D01234R+Dij1234R(1,1))*P(12))
       FI(104) = 4*Is45*(C0134I+Cij134I(1,1)+p2sq*(Dij1234I(2,2)-Dij12
     -   34I(4,2))+s12*(Dij1234I(2,1)-Dij1234I(3,1)+Dij1234I(4,2)-Dij
     -   1234I(6,2))+s45*(Dij1234I(3,1)+Dij1234I(6,2))-2*Dij1234I(12,
     -   3)+(D01234I+Dij1234I(1,1))*P(12))
       F(104)=DCMPLX(FR(104),FI(104))
       FR(105) = -4*Is45*(-Cij134R(1,1)+s23*(D01234R+Dij1234R(1,1))-2*
     -   (p2sq*Dij1234R(2,2)+Dij1234R(7,2)-Dij1234R(11,3))-Dij1234R(5
     -   ,2)*P(5)+Dij1234R(2,1)*P(9)-Dij1234R(3,1)*P(14)+Dij1234R(6,2
     -   )*P(16)+Dij1234R(4,2)*P(46))
       FI(105) = -4*Is45*(-Cij134I(1,1)+s23*(D01234I+Dij1234I(1,1))-2*
     -   (p2sq*Dij1234I(2,2)+Dij1234I(7,2)-Dij1234I(11,3))-Dij1234I(5
     -   ,2)*P(5)+Dij1234I(2,1)*P(9)-Dij1234I(3,1)*P(14)+Dij1234I(6,2
     -   )*P(16)+Dij1234I(4,2)*P(46))
       F(105)=DCMPLX(FR(105),FI(105))
       FR(106) = 4*Is45*(C0134R+Cij134R(2,1)+(-s12+s45)*Dij1234R(3,2)+
     -   s12*Dij1234R(5,2)+p2sq*(Dij1234R(2,1)-Dij1234R(5,2)+Dij1234R
     -   (6,2))-2*(Dij1234R(7,2)+Dij1234R(13,3))+Dij1234R(3,1)*P(17))
       FI(106) = 4*Is45*(C0134I+Cij134I(2,1)+(-s12+s45)*Dij1234I(3,2)+
     -   s12*Dij1234I(5,2)+p2sq*(Dij1234I(2,1)-Dij1234I(5,2)+Dij1234I
     -   (6,2))-2*(Dij1234I(7,2)+Dij1234I(13,3))+Dij1234I(3,1)*P(17))
       F(106)=DCMPLX(FR(106),FI(106))
       FR(107) = 4*Is45*(p3sq*(-Dij1234R(2,1)+Dij1234R(3,1)+Dij1234R(6
     -   ,2))-2*(Dij1234R(7,2)+Dij1234R(12,3))-Dij1234R(2,2)*P(5)+Dij
     -   1234R(4,2)*P(14))
       FI(107) = -4*Is45*(p3sq*(Dij1234I(2,1)-Dij1234I(3,1)-Dij1234I(6
     -   ,2))+2*(Dij1234I(7,2)+Dij1234I(12,3))+Dij1234I(2,2)*P(5)-Dij
     -   1234I(4,2)*P(14))
       F(107)=DCMPLX(FR(107),FI(107))
       FR(108) = 4*Is45*(C0134R+s45*(D01234R+2*Dij1234R(1,1)+Dij1234R(
     -   1,2)-Dij1234R(2,1)-Dij1234R(4,2))-p3sq*Dij1234R(4,2)-p2sq*Di
     -   j1234R(5,2)+s23*(-Dij1234R(1,1)-Dij1234R(1,2)+Dij1234R(2,1)+
     -   Dij1234R(3,1)+Dij1234R(4,2)+Dij1234R(5,2)-Dij1234R(6,2))-6*D
     -   ij1234R(7,2)-2*(p2sq*Dij1234R(2,2)-p2sq*Dij1234R(4,2)+Dij123
     -   4R(11,3))+Dij1234R(2,1)*P(27)+Dij1234R(3,1)*P(28)+Dij1234R(6
     -   ,2)*P(29))
       FI(108) = 4*Is45*(C0134I+s45*(D01234I+2*Dij1234I(1,1)+Dij1234I(
     -   1,2)-Dij1234I(2,1)-Dij1234I(4,2))-p3sq*Dij1234I(4,2)-p2sq*Di
     -   j1234I(5,2)+s23*(-Dij1234I(1,1)-Dij1234I(1,2)+Dij1234I(2,1)+
     -   Dij1234I(3,1)+Dij1234I(4,2)+Dij1234I(5,2)-Dij1234I(6,2))-6*D
     -   ij1234I(7,2)-2*(p2sq*Dij1234I(2,2)-p2sq*Dij1234I(4,2)+Dij123
     -   4I(11,3))+Dij1234I(2,1)*P(27)+Dij1234I(3,1)*P(28)+Dij1234I(6
     -   ,2)*P(29))
       F(108)=DCMPLX(FR(108),FI(108))
       FR(109) = -4*Is45*(C0134R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,2)
     -   -2*(Dij1234R(7,2)-Dij1234R(13,3))+Dij1234R(6,2)*P(5)-Dij1234
     -   R(5,2)*P(14)+Dij1234R(3,1)*P(31))
       FI(109) = -4*Is45*(C0134I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,2)
     -   -2*(Dij1234I(7,2)-Dij1234I(13,3))+Dij1234I(6,2)*P(5)-Dij1234
     -   I(5,2)*P(14)+Dij1234I(3,1)*P(31))
       F(109)=DCMPLX(FR(109),FI(109))
       P(47) = -p2sq+s12
       FR(110) = 4*Is45*(C0134R-Cij134R(1,1)+Cij134R(2,1)-2*p2sq*Dij12
     -   34R(2,2)+s23*(Dij1234R(3,1)+Dij1234R(3,2))+p2sq*Dij1234R(4,2
     -   )-s12*Dij1234R(4,2)-p2sq*Dij1234R(5,2)+s12*Dij1234R(5,2)+(3*
     -   p2sq+p3sq-s23)*Dij1234R(6,2)-6*Dij1234R(7,2)+Dij1234R(2,1)*P
     -   (9)+(D01234R+Dij1234R(1,1))*P(13)-Dij1234R(3,2)*P(29)+Dij123
     -   4R(3,1)*P(47))
       FI(110) = 4*Is45*(C0134I-Cij134I(1,1)+Cij134I(2,1)-2*p2sq*Dij12
     -   34I(2,2)+s23*(Dij1234I(3,1)+Dij1234I(3,2))+p2sq*Dij1234I(4,2
     -   )-s12*Dij1234I(4,2)-p2sq*Dij1234I(5,2)+s12*Dij1234I(5,2)+(3*
     -   p2sq+p3sq-s23)*Dij1234I(6,2)-6*Dij1234I(7,2)+Dij1234I(2,1)*P
     -   (9)+(D01234I+Dij1234I(1,1))*P(13)-Dij1234I(3,2)*P(29)+Dij123
     -   4I(3,1)*P(47))
       F(110)=DCMPLX(FR(110),FI(110))
       FR(111) = Is45*(8*(Dij1234R(7,2)+Dij1234R(12,3)-Dij1234R(13,3))
     -   -4*(D01234R+Dij1234R(1,1))*P(13))
       FI(111) = Is45*(8*(Dij1234I(7,2)+Dij1234I(12,3)-Dij1234I(13,3))
     -   -4*(D01234I+Dij1234I(1,1))*P(13))
       F(111)=DCMPLX(FR(111),FI(111))
       FR(112) = 8*Is45*(2*(Dij1234R(2,1)-Dij1234R(3,1))+Dij1234R(4,3)
     -   +3*(Dij1234R(4,2)-Dij1234R(5,2))-Dij1234R(5,3))
       FI(112) = 8*Is45*(2*(Dij1234I(2,1)-Dij1234I(3,1))+Dij1234I(4,3)
     -   +3*(Dij1234I(4,2)-Dij1234I(5,2))-Dij1234I(5,3))
       F(112)=DCMPLX(FR(112),FI(112))
       FR(113) = 8*Is45*(Dij1234R(2,2)-Dij1234R(5,2)+2*(Dij1234R(2,1)-
     -   Dij1234R(3,1)+Dij1234R(4,2)-Dij1234R(6,2))+Dij1234R(6,3)-Dij
     -   1234R(10,3))
       FI(113) = 8*Is45*(Dij1234I(2,2)-Dij1234I(5,2)+2*(Dij1234I(2,1)-
     -   Dij1234I(3,1)+Dij1234I(4,2)-Dij1234I(6,2))+Dij1234I(6,3)-Dij
     -   1234I(10,3))
       F(113)=DCMPLX(FR(113),FI(113))
       FR(114) = 8*Is45*(Dij1234R(2,1)+2*Dij1234R(2,2)-Dij1234R(3,1)-D
     -   ij1234R(5,2)-Dij1234R(6,2)+Dij1234R(6,3)-Dij1234R(10,3))
       FI(114) = 8*Is45*(Dij1234I(2,1)+2*Dij1234I(2,2)-Dij1234I(3,1)-D
     -   ij1234I(5,2)-Dij1234I(6,2)+Dij1234I(6,3)-Dij1234I(10,3))
       F(114)=DCMPLX(FR(114),FI(114))
       FR(115) = 8*Is45*(Dij1234R(2,1)+Dij1234R(2,3)-Dij1234R(3,1)+2*(
     -   Dij1234R(2,2)-Dij1234R(6,2))-Dij1234R(8,3))
       FI(115) = 8*Is45*(Dij1234I(2,1)+Dij1234I(2,3)-Dij1234I(3,1)+2*(
     -   Dij1234I(2,2)-Dij1234I(6,2))-Dij1234I(8,3))
       F(115)=DCMPLX(FR(115),FI(115))
       FR(116) = 8*Is45*(Dij1234R(2,1)-Dij1234R(3,1)-2*Dij1234R(3,2)+D
     -   ij1234R(4,2)+Dij1234R(6,2)-Dij1234R(7,3)+Dij1234R(10,3))
       FI(116) = 8*Is45*(Dij1234I(2,1)-Dij1234I(3,1)-2*Dij1234I(3,2)+D
     -   ij1234I(4,2)+Dij1234I(6,2)-Dij1234I(7,3)+Dij1234I(10,3))
       F(116)=DCMPLX(FR(116),FI(116))
       FR(117) = -8*Is45*(Dij1234R(3,2)+Dij1234R(5,2)-2*Dij1234R(6,2)+
     -   Dij1234R(7,3)-Dij1234R(10,3))
       FI(117) = -8*Is45*(Dij1234I(3,2)+Dij1234I(5,2)-2*Dij1234I(6,2)+
     -   Dij1234I(7,3)-Dij1234I(10,3))
       F(117)=DCMPLX(FR(117),FI(117))
       FR(118) = 4*Is45*(C0134R+Cij134R(1,1)+p2sq*Dij1234R(2,2)-p2sq*D
     -   ij1234R(4,2)+s12*(Dij1234R(2,1)-Dij1234R(3,1)+Dij1234R(4,2)-
     -   Dij1234R(6,2))-2*Dij1234R(12,3)-(D01234R+Dij1234R(1,1))*P(13
     -   ))
       FI(118) = 4*Is45*(C0134I+Cij134I(1,1)+p2sq*Dij1234I(2,2)-p2sq*D
     -   ij1234I(4,2)+s12*(Dij1234I(2,1)-Dij1234I(3,1)+Dij1234I(4,2)-
     -   Dij1234I(6,2))-2*Dij1234I(12,3)-(D01234I+Dij1234I(1,1))*P(13
     -   ))
       F(118)=DCMPLX(FR(118),FI(118))
       FR(119) = 4*Is45*(C0134R+Cij134R(2,1)+s12*(-Dij1234R(3,2)+Dij12
     -   34R(5,2))+p2sq*(Dij1234R(2,1)-Dij1234R(3,1)-Dij1234R(5,2)+Di
     -   j1234R(6,2))-2*(Dij1234R(7,2)+Dij1234R(13,3)))
       FI(119) = 4*Is45*(C0134I+Cij134I(2,1)+s12*(-Dij1234I(3,2)+Dij12
     -   34I(5,2))+p2sq*(Dij1234I(2,1)-Dij1234I(3,1)-Dij1234I(5,2)+Di
     -   j1234I(6,2))-2*(Dij1234I(7,2)+Dij1234I(13,3)))
       F(119)=DCMPLX(FR(119),FI(119))
       FR(120) = 4*Is45*(C0134R+Cij134R(2,1)-p2sq*Dij1234R(4,2)-4*Dij1
     -   234R(7,2)+2*(p2sq*Dij1234R(6,2)-Dij1234R(11,3))-Dij1234R(3,2
     -   )*P(16)+Dij1234R(5,2)*P(18))
       FI(120) = 4*Is45*(C0134I+Cij134I(2,1)-p2sq*Dij1234I(4,2)-4*Dij1
     -   234I(7,2)+2*(p2sq*Dij1234I(6,2)-Dij1234I(11,3))-Dij1234I(3,2
     -   )*P(16)+Dij1234I(5,2)*P(18))
       F(120)=DCMPLX(FR(120),FI(120))
       FR(121) = 4*Is45*(-(p2sq*Dij1234R(2,2))+p3sq*(Dij1234R(3,1)+Dij
     -   1234R(6,2))-2*(Dij1234R(7,2)+Dij1234R(12,3))+Dij1234R(4,2)*P
     -   (14)+Dij1234R(2,1)*P(25)-Dij1234R(2,2)*P(26))
       FI(121) = 4*Is45*(-(p2sq*Dij1234I(2,2))+p3sq*(Dij1234I(3,1)+Dij
     -   1234I(6,2))-2*(Dij1234I(7,2)+Dij1234I(12,3))+Dij1234I(4,2)*P
     -   (14)+Dij1234I(2,1)*P(25)-Dij1234I(2,2)*P(26))
       F(121)=DCMPLX(FR(121),FI(121))
       FR(122) = -4*Is45*(C0134R+p2sq*Dij1234R(2,1)-p3sq*Dij1234R(3,2)
     -   -2*(Dij1234R(7,2)-Dij1234R(13,3))+Dij1234R(6,2)*P(2)-Dij1234
     -   R(5,2)*P(14)-Dij1234R(3,1)*P(30))
       FI(122) = -4*Is45*(C0134I+p2sq*Dij1234I(2,1)-p3sq*Dij1234I(3,2)
     -   -2*(Dij1234I(7,2)-Dij1234I(13,3))+Dij1234I(6,2)*P(2)-Dij1234
     -   I(5,2)*P(14)-Dij1234I(3,1)*P(30))
       F(122)=DCMPLX(FR(122),FI(122))
       FR(123) = 4*Is45*(Cij134R(1,1)-Cij134R(2,1)+s45*Dij1234R(1,2)+p
     -   3sq*Dij1234R(3,2)-p3sq*Dij1234R(4,2)+s12*Dij1234R(4,2)-s12*D
     -   ij1234R(5,2)+s23*(-Dij1234R(1,2)+Dij1234R(2,1)-Dij1234R(3,2)
     -   +Dij1234R(4,2)+Dij1234R(5,2))+p2sq*(Dij1234R(3,2)+Dij1234R(4
     -   ,2)-2*Dij1234R(6,2))+Dij1234R(3,1)*P(3)+Dij1234R(2,1)*P(8)+D
     -   01234R*P(12)-2*(Dij1234R(11,3)+Dij1234R(1,1)*P(13)))
       FI(123) = 4*Is45*(Cij134I(1,1)-Cij134I(2,1)+s45*Dij1234I(1,2)+p
     -   3sq*Dij1234I(3,2)-p3sq*Dij1234I(4,2)+s12*Dij1234I(4,2)-s12*D
     -   ij1234I(5,2)+s23*(-Dij1234I(1,2)+Dij1234I(2,1)-Dij1234I(3,2)
     -   +Dij1234I(4,2)+Dij1234I(5,2))+p2sq*(Dij1234I(3,2)+Dij1234I(4
     -   ,2)-2*Dij1234I(6,2))+Dij1234I(3,1)*P(3)+Dij1234I(2,1)*P(8)+D
     -   01234I*P(12)-2*(Dij1234I(11,3)+Dij1234I(1,1)*P(13)))
       F(123)=DCMPLX(FR(123),FI(123))
       FR(124) = 2*Is45*(C0134R-Cij134R(1,1)+Cij134R(2,1)-2*p2sq*Dij12
     -   34R(2,2)+s23*(Dij1234R(3,1)+Dij1234R(3,2))+p2sq*Dij1234R(4,2
     -   )-s12*Dij1234R(4,2)-p2sq*Dij1234R(5,2)+s12*Dij1234R(5,2)+(3*
     -   p2sq+p3sq-s23)*Dij1234R(6,2)-6*Dij1234R(7,2)+Dij1234R(2,1)*P
     -   (9)+(D01234R+Dij1234R(1,1))*P(13)-Dij1234R(3,2)*P(29)+Dij123
     -   4R(3,1)*P(47))
       FI(124) = 2*Is45*(C0134I-Cij134I(1,1)+Cij134I(2,1)-2*p2sq*Dij12
     -   34I(2,2)+s23*(Dij1234I(3,1)+Dij1234I(3,2))+p2sq*Dij1234I(4,2
     -   )-s12*Dij1234I(4,2)-p2sq*Dij1234I(5,2)+s12*Dij1234I(5,2)+(3*
     -   p2sq+p3sq-s23)*Dij1234I(6,2)-6*Dij1234I(7,2)+Dij1234I(2,1)*P
     -   (9)+(D01234I+Dij1234I(1,1))*P(13)-Dij1234I(3,2)*P(29)+Dij123
     -   4I(3,1)*P(47))
       F(124)=DCMPLX(FR(124),FI(124))
       P(48) = s12+s23
       FR(125) = -8*Is12*(s12*(Dij1345R(1,1)+Dij1345R(1,3)+2*(Dij1345R
     -   (1,2)+Dij1345R(3,2))-Dij1345R(3,3)-4*Dij1345R(5,2)-3*(Dij134
     -   5R(5,3)-Dij1345R(7,3)))+s23*(-Dij1345R(3,3)+Dij1345R(4,3)-Di
     -   j1345R(5,3)+Dij1345R(9,3)+2*(Dij1345R(3,2)+Dij1345R(4,2)-Dij
     -   1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij1345R(10,3)))-Dij1
     -   345R(2,1)*P(16)+(Dij1345R(3,1)+Dij1345R(3,3)-Dij1345R(4,3)+D
     -   ij1345R(5,3)-Dij1345R(9,3)-2*(Dij1345R(3,2)+Dij1345R(4,2)-Di
     -   j1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij1345R(10,3)))*P(2
     -   9)-Dij1345R(3,1)*P(48))
       FI(125) = -8*Is12*(s12*(Dij1345I(1,1)+Dij1345I(1,3)+2*(Dij1345I
     -   (1,2)+Dij1345I(3,2))-Dij1345I(3,3)-4*Dij1345I(5,2)-3*(Dij134
     -   5I(5,3)-Dij1345I(7,3)))+s23*(-Dij1345I(3,3)+Dij1345I(4,3)-Di
     -   j1345I(5,3)+Dij1345I(9,3)+2*(Dij1345I(3,2)+Dij1345I(4,2)-Dij
     -   1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij1345I(10,3)))-Dij1
     -   345I(2,1)*P(16)+(Dij1345I(3,1)+Dij1345I(3,3)-Dij1345I(4,3)+D
     -   ij1345I(5,3)-Dij1345I(9,3)-2*(Dij1345I(3,2)+Dij1345I(4,2)-Di
     -   j1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij1345I(10,3)))*P(2
     -   9)-Dij1345I(3,1)*P(48))
       F(125)=DCMPLX(FR(125),FI(125))
       P(49) = 3-Is12*P(16)
       P(50) = 3-2*Is12*P(16)
       P(51) = -p2sq+s12+s23
       P(52) = p3sq-p4sq-3*s12+s34+s45+2*P(35)
       P(53) = p3sq+2*P(35)
       P(54) = p3sq-p4sq-5*s12+s34+s45+2*P(35)
       FR(126) = -8*(Dij1345R(1,1)+2*Dij1345R(1,2)+Dij1345R(1,3)-Dij13
     -   45R(5,3)*P(49)+Dij1345R(7,3)*P(50)-Is12*(-2*(Dij1345R(7,2)+D
     -   ij1345R(11,3)-Dij1345R(13,3))+Dij1345R(3,3)*P(1)+(Dij1345R(4
     -   ,3)+Dij1345R(9,3)-2*Dij1345R(10,3))*P(16)+Dij1345R(2,1)*P(35
     -   )+Dij1345R(3,1)*P(51)+Dij1345R(3,2)*P(52)+(Dij1345R(4,2)-Dij
     -   1345R(6,2))*P(53)-Dij1345R(5,2)*P(54)))
       FI(126) = -8*(Dij1345I(1,1)+2*Dij1345I(1,2)+Dij1345I(1,3)-Dij13
     -   45I(5,3)*P(49)+Dij1345I(7,3)*P(50)-Is12*(-2*(Dij1345I(7,2)+D
     -   ij1345I(11,3)-Dij1345I(13,3))+Dij1345I(3,3)*P(1)+(Dij1345I(4
     -   ,3)+Dij1345I(9,3)-2*Dij1345I(10,3))*P(16)+Dij1345I(2,1)*P(35
     -   )+Dij1345I(3,1)*P(51)+Dij1345I(3,2)*P(52)+(Dij1345I(4,2)-Dij
     -   1345I(6,2))*P(53)-Dij1345I(5,2)*P(54)))
       F(126)=DCMPLX(FR(126),FI(126))
       P(55) = s12-s34
       P(56) = -p2sq+p4sq
       P(57) = p2sq+p3sq-p4sq
       P(58) = s34+s45
       FR(127) = 8*Is12*(C0134R-s34*Dij1345R(2,1)+s12*(-Dij1345R(1,2)-
     -   Dij1345R(2,2)+Dij1345R(3,1)+Dij1345R(3,3)-Dij1345R(4,3)+Dij1
     -   345R(5,3)-2*(Dij1345R(3,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij13
     -   45R(7,3)-Dij1345R(10,3)))-p3sq*(Dij1345R(3,3)-Dij1345R(4,2)+
     -   Dij1345R(5,2)+Dij1345R(6,2)-Dij1345R(6,3)-Dij1345R(7,3)+Dij1
     -   345R(8,3)+2*Dij1345R(10,3))-s23*(Dij1345R(2,2)+Dij1345R(3,2)
     -   -Dij1345R(3,3)+Dij1345R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)-2*(
     -   Dij1345R(6,2)+Dij1345R(10,3)))-p2sq*(Dij1345R(3,3)-Dij1345R(
     -   6,3)-Dij1345R(7,3)+Dij1345R(8,3)+2*(Dij1345R(6,2)+Dij1345R(1
     -   0,3)))+2*(p4sq*Dij1345R(6,2)-2*Dij1345R(7,2)-Dij1345R(12,3)+
     -   Dij1345R(13,3))-Dij1345R(9,3)*P(40)-Dij1345R(1,1)*P(55)-Dij1
     -   345R(2,2)*P(56)+Dij1345R(3,2)*P(57)+(Dij1345R(2,2)-Dij1345R(
     -   4,2)+Dij1345R(5,2)-Dij1345R(6,2))*P(58))
       FI(127) = 8*Is12*(C0134I-s34*Dij1345I(2,1)+s12*(-Dij1345I(1,2)-
     -   Dij1345I(2,2)+Dij1345I(3,1)+Dij1345I(3,3)-Dij1345I(4,3)+Dij1
     -   345I(5,3)-2*(Dij1345I(3,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij13
     -   45I(7,3)-Dij1345I(10,3)))-p3sq*(Dij1345I(3,3)-Dij1345I(4,2)+
     -   Dij1345I(5,2)+Dij1345I(6,2)-Dij1345I(6,3)-Dij1345I(7,3)+Dij1
     -   345I(8,3)+2*Dij1345I(10,3))-s23*(Dij1345I(2,2)+Dij1345I(3,2)
     -   -Dij1345I(3,3)+Dij1345I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)-2*(
     -   Dij1345I(6,2)+Dij1345I(10,3)))-p2sq*(Dij1345I(3,3)-Dij1345I(
     -   6,3)-Dij1345I(7,3)+Dij1345I(8,3)+2*(Dij1345I(6,2)+Dij1345I(1
     -   0,3)))+2*(p4sq*Dij1345I(6,2)-2*Dij1345I(7,2)-Dij1345I(12,3)+
     -   Dij1345I(13,3))-Dij1345I(9,3)*P(40)-Dij1345I(1,1)*P(55)-Dij1
     -   345I(2,2)*P(56)+Dij1345I(3,2)*P(57)+(Dij1345I(2,2)-Dij1345I(
     -   4,2)+Dij1345I(5,2)-Dij1345I(6,2))*P(58))
       F(127)=DCMPLX(FR(127),FI(127))
       FR(128) = 8*Is12*(D01345R*s12+s12*(Dij1345R(1,2)-Dij1345R(3,1)-
     -   Dij1345R(3,2)+Dij1345R(3,3)+Dij1345R(5,3)+2*(Dij1345R(1,1)-D
     -   ij1345R(7,3)))-s23*(Dij1345R(3,1)+Dij1345R(3,2)-Dij1345R(3,3
     -   )-Dij1345R(4,2)+Dij1345R(7,3)+Dij1345R(9,3)-Dij1345R(10,3))-
     -   Dij1345R(2,1)*P(16)+(Dij1345R(3,1)+Dij1345R(3,2)-Dij1345R(3,
     -   3)-Dij1345R(4,2)+Dij1345R(7,3)+Dij1345R(9,3)-Dij1345R(10,3))
     -   *P(29))
       FI(128) = 8*Is12*(D01345I*s12+s12*(Dij1345I(1,2)-Dij1345I(3,1)-
     -   Dij1345I(3,2)+Dij1345I(3,3)+Dij1345I(5,3)+2*(Dij1345I(1,1)-D
     -   ij1345I(7,3)))-s23*(Dij1345I(3,1)+Dij1345I(3,2)-Dij1345I(3,3
     -   )-Dij1345I(4,2)+Dij1345I(7,3)+Dij1345I(9,3)-Dij1345I(10,3))-
     -   Dij1345I(2,1)*P(16)+(Dij1345I(3,1)+Dij1345I(3,2)-Dij1345I(3,
     -   3)-Dij1345I(4,2)+Dij1345I(7,3)+Dij1345I(9,3)-Dij1345I(10,3))
     -   *P(29))
       F(128)=DCMPLX(FR(128),FI(128))
       P(59) = 1-Is12*s34
       P(60) = -p2sq+s12+s23-s45
       P(61) = p3sq-3*s12+s34+s45+2*P(35)
       P(62) = p3sq-4*s12+s34+s45+2*P(35)
       FR(129) = 8*(-Dij1345R(1,2)-Dij1345R(1,3)+Dij1345R(5,3)*P(49)-D
     -   ij1345R(7,3)*P(50)+(D01345R+Dij1345R(1,1))*P(59)+Is12*(C0134
     -   R+Cij134R(1,1)-2*Dij1345R(7,2)-4*(Dij1345R(11,3)-Dij1345R(13
     -   ,3))+Dij1345R(3,3)*P(1)+Dij1345R(2,1)*P(14)+(Dij1345R(4,3)+D
     -   ij1345R(9,3)-2*Dij1345R(10,3))*P(16)+(Dij1345R(4,2)-Dij1345R
     -   (6,2))*P(36)+Dij1345R(3,1)*P(60)+Dij1345R(3,2)*P(61)-Dij1345
     -   R(5,2)*P(62)))
       FI(129) = 8*(-Dij1345I(1,2)-Dij1345I(1,3)+Dij1345I(5,3)*P(49)-D
     -   ij1345I(7,3)*P(50)+(D01345I+Dij1345I(1,1))*P(59)+Is12*(C0134
     -   I+Cij134I(1,1)-2*Dij1345I(7,2)-4*(Dij1345I(11,3)-Dij1345I(13
     -   ,3))+Dij1345I(3,3)*P(1)+Dij1345I(2,1)*P(14)+(Dij1345I(4,3)+D
     -   ij1345I(9,3)-2*Dij1345I(10,3))*P(16)+(Dij1345I(4,2)-Dij1345I
     -   (6,2))*P(36)+Dij1345I(3,1)*P(60)+Dij1345I(3,2)*P(61)-Dij1345
     -   I(5,2)*P(62)))
       F(129)=DCMPLX(FR(129),FI(129))
       P(63) = p2sq+p3sq-s23+s45
       FR(130) = 8*(D01345R+Dij1345R(1,1)-Dij1345R(1,2)-Dij1345R(1,3)-
     -   Dij1345R(2,1)+Dij1345R(3,1)-2*Dij1345R(3,2)+Dij1345R(3,3)-Di
     -   j1345R(4,2)+Dij1345R(6,2)+3*(Dij1345R(5,2)+Dij1345R(5,3)-Dij
     -   1345R(7,3))+Is12*(C0134R+Cij134R(1,1)-s34*(D01345R+Dij1345R(
     -   1,1))+p4sq*Dij1345R(3,2)+s23*Dij1345R(3,3)+s45*Dij1345R(4,2)
     -   -s23*Dij1345R(4,3)-p4sq*Dij1345R(5,2)+s23*Dij1345R(5,3)-s45*
     -   Dij1345R(6,2)+p2sq*(-Dij1345R(3,3)+Dij1345R(4,3)-Dij1345R(5,
     -   3)+2*(Dij1345R(3,2)+Dij1345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2
     -   )+Dij1345R(7,3)-Dij1345R(10,3)))+p3sq*(-Dij1345R(3,3)+Dij134
     -   5R(4,3)-Dij1345R(5,3)+2*(Dij1345R(3,2)+Dij1345R(4,2)-Dij1345
     -   R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij1345R(10,3)))-2*(s23*D
     -   ij1345R(3,2)+s23*Dij1345R(4,2)-s23*Dij1345R(5,2)-s23*Dij1345
     -   R(6,2)+s23*Dij1345R(7,3)-s23*Dij1345R(10,3)+Dij1345R(11,3)-D
     -   ij1345R(13,3))+Dij1345R(9,3)*P(16)+(Dij1345R(2,1)-Dij1345R(3
     -   ,1))*P(63)))
       FI(130) = 8*(D01345I+Dij1345I(1,1)-Dij1345I(1,2)-Dij1345I(1,3)-
     -   Dij1345I(2,1)+Dij1345I(3,1)-2*Dij1345I(3,2)+Dij1345I(3,3)-Di
     -   j1345I(4,2)+Dij1345I(6,2)+3*(Dij1345I(5,2)+Dij1345I(5,3)-Dij
     -   1345I(7,3))+Is12*(C0134I+Cij134I(1,1)-s34*(D01345I+Dij1345I(
     -   1,1))+p4sq*Dij1345I(3,2)+s23*Dij1345I(3,3)+s45*Dij1345I(4,2)
     -   -s23*Dij1345I(4,3)-p4sq*Dij1345I(5,2)+s23*Dij1345I(5,3)-s45*
     -   Dij1345I(6,2)+p2sq*(-Dij1345I(3,3)+Dij1345I(4,3)-Dij1345I(5,
     -   3)+2*(Dij1345I(3,2)+Dij1345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2
     -   )+Dij1345I(7,3)-Dij1345I(10,3)))+p3sq*(-Dij1345I(3,3)+Dij134
     -   5I(4,3)-Dij1345I(5,3)+2*(Dij1345I(3,2)+Dij1345I(4,2)-Dij1345
     -   I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij1345I(10,3)))-2*(s23*D
     -   ij1345I(3,2)+s23*Dij1345I(4,2)-s23*Dij1345I(5,2)-s23*Dij1345
     -   I(6,2)+s23*Dij1345I(7,3)-s23*Dij1345I(10,3)+Dij1345I(11,3)-D
     -   ij1345I(13,3))+Dij1345I(9,3)*P(16)+(Dij1345I(2,1)-Dij1345I(3
     -   ,1))*P(63)))
       F(130)=DCMPLX(FR(130),FI(130))
       P(64) = p2sq+p3sq+p4sq-s23-s34
       P(65) = p3sq+p4sq-3*s12+s34+s45+2*P(35)
       P(66) = p2sq-2*s12-s23+s34+s45
       P(67) = p3sq+p4sq
       P(68) = s45+3*P(43)+2*P(67)
       P(69) = p2sq+p3sq-2*s12-s23
       P(70) = 1-2*Is12*P(16)
       FR(131) = 8*(-Dij1345R(4,3)+Dij1345R(5,3)+(D01345R+Dij1345R(1,1
     -   ))*P(59)+Is12*(Cij134R(2,1)+2*(Dij1345R(7,2)-Dij1345R(12,3)+
     -   Dij1345R(13,3))+(Dij1345R(3,3)+2*Dij1345R(10,3))*P(1)+Dij134
     -   5R(2,1)*P(14)+(Dij1345R(6,3)-Dij1345R(8,3))*P(16)+Dij1345R(3
     -   ,1)*P(60)+Dij1345R(2,2)*P(64)+Dij1345R(3,2)*P(65)+(Dij1345R(
     -   4,2)-Dij1345R(5,2))*P(66)-Dij1345R(6,2)*P(68)+Dij1345R(7,3)*
     -   P(69))-Dij1345R(9,3)*P(70))
       FI(131) = 8*(-Dij1345I(4,3)+Dij1345I(5,3)+(D01345I+Dij1345I(1,1
     -   ))*P(59)+Is12*(Cij134I(2,1)+2*(Dij1345I(7,2)-Dij1345I(12,3)+
     -   Dij1345I(13,3))+(Dij1345I(3,3)+2*Dij1345I(10,3))*P(1)+Dij134
     -   5I(2,1)*P(14)+(Dij1345I(6,3)-Dij1345I(8,3))*P(16)+Dij1345I(3
     -   ,1)*P(60)+Dij1345I(2,2)*P(64)+Dij1345I(3,2)*P(65)+(Dij1345I(
     -   4,2)-Dij1345I(5,2))*P(66)-Dij1345I(6,2)*P(68)+Dij1345I(7,3)*
     -   P(69))-Dij1345I(9,3)*P(70))
       F(131)=DCMPLX(FR(131),FI(131))
       P(71) = s12+s34
       FR(132) = 8*Is12*(C0134R+D01345R*s12-s34*Dij1345R(2,1)+(s34+s45
     -   )*Dij1345R(5,2)-p2sq*Dij1345R(6,2)+p4sq*(Dij1345R(2,1)-Dij13
     -   45R(3,1)-Dij1345R(3,2)+Dij1345R(6,2))+s12*(-Dij1345R(3,2)+Di
     -   j1345R(3,3)+Dij1345R(5,3)+Dij1345R(6,2)-2*Dij1345R(7,3))+s23
     -   *(Dij1345R(3,3)+Dij1345R(6,2)-Dij1345R(7,3)-Dij1345R(9,3)+Di
     -   j1345R(10,3))-2*(Dij1345R(7,2)-Dij1345R(13,3))+(-Dij1345R(3,
     -   3)+Dij1345R(7,3)+Dij1345R(9,3)-Dij1345R(10,3))*P(29)+Dij1345
     -   R(5,2)*P(35)-Dij1345R(6,2)*P(58)+Dij1345R(1,1)*P(71))
       FI(132) = 8*Is12*(C0134I+D01345I*s12-s34*Dij1345I(2,1)+(s34+s45
     -   )*Dij1345I(5,2)-p2sq*Dij1345I(6,2)+p4sq*(Dij1345I(2,1)-Dij13
     -   45I(3,1)-Dij1345I(3,2)+Dij1345I(6,2))+s12*(-Dij1345I(3,2)+Di
     -   j1345I(3,3)+Dij1345I(5,3)+Dij1345I(6,2)-2*Dij1345I(7,3))+s23
     -   *(Dij1345I(3,3)+Dij1345I(6,2)-Dij1345I(7,3)-Dij1345I(9,3)+Di
     -   j1345I(10,3))-2*(Dij1345I(7,2)-Dij1345I(13,3))+(-Dij1345I(3,
     -   3)+Dij1345I(7,3)+Dij1345I(9,3)-Dij1345I(10,3))*P(29)+Dij1345
     -   I(5,2)*P(35)-Dij1345I(6,2)*P(58)+Dij1345I(1,1)*P(71))
       F(132)=DCMPLX(FR(132),FI(132))
       P(72) = p3sq*s12-p4sq*s12-s12**2-s12*s23-p3sq*s34+s12*s34+s23*s
     -   34+s12*s45+p2sq*P(55)
       P(73) = p2sq*s12+p3sq*s12-p4sq*s12-s12**2-s12*s23-p2sq*s34-p3sq
     -   *s34+s12*s34+s23*s34+s12*s45
       FR(133) = -8*(Dij1345R(7,2)+Dij1345R(11,3)-Dij1345R(13,3))+4*Is
     -   12*(2*(Dij1345R(7,2)+Dij1345R(12,3)-Dij1345R(13,3))*P(16)+D0
     -   1345R*P(72)+Dij1345R(1,1)*P(73))
       FI(133) = -8*(Dij1345I(7,2)+Dij1345I(11,3)-Dij1345I(13,3))+4*Is
     -   12*(2*(Dij1345I(7,2)+Dij1345I(12,3)-Dij1345I(13,3))*P(16)+D0
     -   1345I*P(72)+Dij1345I(1,1)*P(73))
       F(133)=DCMPLX(FR(133),FI(133))
       P(74) = s12+2*s23
       FR(134) = Is12*(s12*(Dij1345R(1,1)+Dij1345R(1,2)-Dij1345R(3,1)-
     -   Dij1345R(3,3)+Dij1345R(4,2)+Dij1345R(4,3)-3*Dij1345R(5,2)-Di
     -   j1345R(5,3)-Dij1345R(6,2)+2*(Dij1345R(3,2)+Dij1345R(7,3)-Dij
     -   1345R(10,3)))+s23*(Dij1345R(3,2)-Dij1345R(3,3)+Dij1345R(6,3)
     -   +Dij1345R(7,3)-Dij1345R(8,3)-2*(Dij1345R(6,2)+Dij1345R(10,3)
     -   ))-Dij1345R(2,2)*P(16)-(Dij1345R(3,2)-Dij1345R(3,3)+Dij1345R
     -   (6,3)+Dij1345R(7,3)-Dij1345R(8,3)-2*(Dij1345R(6,2)-Dij1345R(
     -   9,3)+Dij1345R(10,3)))*P(29)+Dij1345R(9,3)*P(74))
       FI(134) = Is12*(s12*(Dij1345I(1,1)+Dij1345I(1,2)-Dij1345I(3,1)-
     -   Dij1345I(3,3)+Dij1345I(4,2)+Dij1345I(4,3)-3*Dij1345I(5,2)-Di
     -   j1345I(5,3)-Dij1345I(6,2)+2*(Dij1345I(3,2)+Dij1345I(7,3)-Dij
     -   1345I(10,3)))+s23*(Dij1345I(3,2)-Dij1345I(3,3)+Dij1345I(6,3)
     -   +Dij1345I(7,3)-Dij1345I(8,3)-2*(Dij1345I(6,2)+Dij1345I(10,3)
     -   ))-Dij1345I(2,2)*P(16)-(Dij1345I(3,2)-Dij1345I(3,3)+Dij1345I
     -   (6,3)+Dij1345I(7,3)-Dij1345I(8,3)-2*(Dij1345I(6,2)-Dij1345I(
     -   9,3)+Dij1345I(10,3)))*P(29)+Dij1345I(9,3)*P(74))
       F(134)=DCMPLX(FR(134),FI(134))
       P(75) = p2sq+p3sq+s12-s23
       FR(135) = Dij1345R(3,2)-Dij1345R(3,3)-Dij1345R(5,3)+2*Dij1345R(
     -   7,3)+Is12*((Dij1345R(3,3)+Dij1345R(6,2)-Dij1345R(7,3)-Dij134
     -   5R(9,3)+Dij1345R(10,3))*P(16)-Dij1345R(5,2)*P(75))
       FI(135) = Dij1345I(3,2)-Dij1345I(3,3)-Dij1345I(5,3)+2*Dij1345I(
     -   7,3)+Is12*((Dij1345I(3,3)+Dij1345I(6,2)-Dij1345I(7,3)-Dij134
     -   5I(9,3)+Dij1345I(10,3))*P(16)-Dij1345I(5,2)*P(75))
       F(135)=DCMPLX(FR(135),FI(135))
       FR(136) = Dij1345R(3,1)+Dij1345R(3,3)-Dij1345R(4,3)+Dij1345R(5,
     -   3)-Dij1345R(9,3)-2*(Dij1345R(3,2)+Dij1345R(4,2)-Dij1345R(5,2
     -   )-Dij1345R(6,2)+Dij1345R(7,3)-Dij1345R(10,3))+Is12*(Dij1345R
     -   (2,1)*P(2)+(Dij1345R(2,2)-Dij1345R(3,1)-Dij1345R(3,3)+Dij134
     -   5R(4,2)-Dij1345R(5,2)-3*Dij1345R(6,2)+Dij1345R(6,3)+Dij1345R
     -   (7,3)-Dij1345R(8,3)+2*(Dij1345R(3,2)+Dij1345R(9,3)-Dij1345R(
     -   10,3)))*P(16))
       FI(136) = Dij1345I(3,1)+Dij1345I(3,3)-Dij1345I(4,3)+Dij1345I(5,
     -   3)-Dij1345I(9,3)-2*(Dij1345I(3,2)+Dij1345I(4,2)-Dij1345I(5,2
     -   )-Dij1345I(6,2)+Dij1345I(7,3)-Dij1345I(10,3))+Is12*(Dij1345I
     -   (2,1)*P(2)+(Dij1345I(2,2)-Dij1345I(3,1)-Dij1345I(3,3)+Dij134
     -   5I(4,2)-Dij1345I(5,2)-3*Dij1345I(6,2)+Dij1345I(6,3)+Dij1345I
     -   (7,3)-Dij1345I(8,3)+2*(Dij1345I(3,2)+Dij1345I(9,3)-Dij1345I(
     -   10,3)))*P(16))
       F(136)=DCMPLX(FR(136),FI(136))
       P(76) = 2*s12+3*s23
       FR(137) = Is12*(s23*(Dij1345R(2,2)+Dij1345R(2,3)+Dij1345R(3,2)-
     -   Dij1345R(3,3)-2*Dij1345R(6,2)-3*Dij1345R(8,3))+s12*(Dij1345R
     -   (2,1)+Dij1345R(2,2)-Dij1345R(3,1)-Dij1345R(3,3)+Dij1345R(4,2
     -   )-Dij1345R(5,2)-3*Dij1345R(6,2)+Dij1345R(6,3)+Dij1345R(7,3)-
     -   Dij1345R(8,3)+2*(Dij1345R(3,2)-Dij1345R(10,3)))-(Dij1345R(2,
     -   2)+Dij1345R(2,3)+Dij1345R(3,2)-Dij1345R(3,3)-2*Dij1345R(6,2)
     -   -3*(Dij1345R(8,3)-Dij1345R(9,3)))*P(29)+Dij1345R(9,3)*P(76))
       FI(137) = Is12*(s23*(Dij1345I(2,2)+Dij1345I(2,3)+Dij1345I(3,2)-
     -   Dij1345I(3,3)-2*Dij1345I(6,2)-3*Dij1345I(8,3))+s12*(Dij1345I
     -   (2,1)+Dij1345I(2,2)-Dij1345I(3,1)-Dij1345I(3,3)+Dij1345I(4,2
     -   )-Dij1345I(5,2)-3*Dij1345I(6,2)+Dij1345I(6,3)+Dij1345I(7,3)-
     -   Dij1345I(8,3)+2*(Dij1345I(3,2)-Dij1345I(10,3)))-(Dij1345I(2,
     -   2)+Dij1345I(2,3)+Dij1345I(3,2)-Dij1345I(3,3)-2*Dij1345I(6,2)
     -   -3*(Dij1345I(8,3)-Dij1345I(9,3)))*P(29)+Dij1345I(9,3)*P(76))
       F(137)=DCMPLX(FR(137),FI(137))
       FR(138) = -Dij1345R(3,1)-Dij1345R(3,2)+Dij1345R(3,3)+Dij1345R(4
     -   ,2)-Dij1345R(7,3)-Dij1345R(9,3)+Dij1345R(10,3)+Is12*(Dij1345
     -   R(2,1)*P(1)-(Dij1345R(2,2)-Dij1345R(3,1)+Dij1345R(3,3)-Dij13
     -   45R(6,2)+Dij1345R(8,3)-2*Dij1345R(9,3))*P(16))
       FI(138) = -Dij1345I(3,1)-Dij1345I(3,2)+Dij1345I(3,3)+Dij1345I(4
     -   ,2)-Dij1345I(7,3)-Dij1345I(9,3)+Dij1345I(10,3)+Is12*(Dij1345
     -   I(2,1)*P(1)-(Dij1345I(2,2)-Dij1345I(3,1)+Dij1345I(3,3)-Dij13
     -   45I(6,2)+Dij1345I(8,3)-2*Dij1345I(9,3))*P(16))
       F(138)=DCMPLX(FR(138),FI(138))
       FR(139) = D01345R+Dij1345R(1,1)-Dij1345R(3,2)+Dij1345R(6,2)-Dij
     -   1345R(7,3)-Dij1345R(9,3)+Dij1345R(10,3)+Is12*(Dij1345R(3,3)*
     -   P(1)-(Dij1345R(8,3)-2*Dij1345R(9,3))*P(16))
       FI(139) = D01345I+Dij1345I(1,1)-Dij1345I(3,2)+Dij1345I(6,2)-Dij
     -   1345I(7,3)-Dij1345I(9,3)+Dij1345I(10,3)+Is12*(Dij1345I(3,3)*
     -   P(1)-(Dij1345I(8,3)-2*Dij1345I(9,3))*P(16))
       F(139)=DCMPLX(FR(139),FI(139))
       FR(140) = Is12*(-(s12*(Dij1345R(3,1)+Dij1345R(5,2)+Dij1345R(7,3
     -   )))-Dij1345R(3,3)*P(2)+(-Dij1345R(3,2)+Dij1345R(6,2)+Dij1345
     -   R(9,3))*P(16))
       FI(140) = Is12*(-(s12*(Dij1345I(3,1)+Dij1345I(5,2)+Dij1345I(7,3
     -   )))-Dij1345I(3,3)*P(2)+(-Dij1345I(3,2)+Dij1345I(6,2)+Dij1345
     -   I(9,3))*P(16))
       F(140)=DCMPLX(FR(140),FI(140))
       P(77) = p2sq+p3sq+p4sq-s23
       FR(141) = 8*(Dij1345R(1,1)+Dij1345R(1,2)+Dij1345R(3,3)+Dij1345R
     -   (5,3)-2*(Dij1345R(3,1)+Dij1345R(3,2)+Dij1345R(7,3))-Is12*(C0
     -   134R+s45*(-Dij1345R(3,2)+Dij1345R(5,2))+s34*(Dij1345R(1,1)-D
     -   ij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,2)+Dij1345R(5,2)-Dij13
     -   45R(6,2))+p4sq*(-Dij1345R(3,2)+Dij1345R(6,2))+s23*(Dij1345R(
     -   3,2)-Dij1345R(3,3)-Dij1345R(4,2)+Dij1345R(7,3)-Dij1345R(10,3
     -   ))+p2sq*(-Dij1345R(3,2)+Dij1345R(3,3)+Dij1345R(4,2)-Dij1345R
     -   (7,3)+Dij1345R(10,3))+p3sq*(Dij1345R(3,3)+Dij1345R(4,2)-Dij1
     -   345R(5,2)-Dij1345R(7,3)+Dij1345R(10,3))-2*(2*Dij1345R(7,2)+D
     -   ij1345R(13,3))-Dij1345R(9,3)*P(16)+(Dij1345R(2,1)-Dij1345R(3
     -   ,1))*P(77)))
       FI(141) = 8*(Dij1345I(1,1)+Dij1345I(1,2)+Dij1345I(3,3)+Dij1345I
     -   (5,3)-2*(Dij1345I(3,1)+Dij1345I(3,2)+Dij1345I(7,3))-Is12*(C0
     -   134I+s45*(-Dij1345I(3,2)+Dij1345I(5,2))+s34*(Dij1345I(1,1)-D
     -   ij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,2)+Dij1345I(5,2)-Dij13
     -   45I(6,2))+p4sq*(-Dij1345I(3,2)+Dij1345I(6,2))+s23*(Dij1345I(
     -   3,2)-Dij1345I(3,3)-Dij1345I(4,2)+Dij1345I(7,3)-Dij1345I(10,3
     -   ))+p2sq*(-Dij1345I(3,2)+Dij1345I(3,3)+Dij1345I(4,2)-Dij1345I
     -   (7,3)+Dij1345I(10,3))+p3sq*(Dij1345I(3,3)+Dij1345I(4,2)-Dij1
     -   345I(5,2)-Dij1345I(7,3)+Dij1345I(10,3))-2*(2*Dij1345I(7,2)+D
     -   ij1345I(13,3))-Dij1345I(9,3)*P(16)+(Dij1345I(2,1)-Dij1345I(3
     -   ,1))*P(77)))
       F(141)=DCMPLX(FR(141),FI(141))
       P(78) = p3sq+p4sq+s12-s34-s45
       FR(142) = -8*Is12*(p3sq*(Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R(4
     -   ,2)-Dij1345R(6,2))+2*(Dij1345R(7,2)+Dij1345R(11,3)-Dij1345R(
     -   13,3))+(Dij1345R(3,2)-Dij1345R(5,2))*P(78))
       FI(142) = -8*Is12*(p3sq*(Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I(4
     -   ,2)-Dij1345I(6,2))+2*(Dij1345I(7,2)+Dij1345I(11,3)-Dij1345I(
     -   13,3))+(Dij1345I(3,2)-Dij1345I(5,2))*P(78))
       F(142)=DCMPLX(FR(142),FI(142))
       P(79) = p4sq-s34
       P(80) = p3sq-2*p4sq+s12-s45
       FR(143) = 8*Is12*(Cij134R(2,1)+p4sq*Dij1345R(3,2)-s12*Dij1345R(
     -   3,2)+s34*Dij1345R(3,2)+s34*Dij1345R(4,2)-p3sq*(Dij1345R(2,1)
     -   -Dij1345R(3,1)+Dij1345R(3,2)+Dij1345R(4,2)-Dij1345R(5,2))+s4
     -   5*(Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R(3,2)+Dij1345R(4,2)-D
     -   ij1345R(5,2))-s34*Dij1345R(5,2)+2*(Dij1345R(7,2)-Dij1345R(12
     -   ,3)+Dij1345R(13,3))+(D01345R+Dij1345R(1,1))*P(55)+Dij1345R(2
     -   ,2)*P(79)+Dij1345R(6,2)*P(80))
       FI(143) = 8*Is12*(Cij134I(2,1)+p4sq*Dij1345I(3,2)-s12*Dij1345I(
     -   3,2)+s34*Dij1345I(3,2)+s34*Dij1345I(4,2)-p3sq*(Dij1345I(2,1)
     -   -Dij1345I(3,1)+Dij1345I(3,2)+Dij1345I(4,2)-Dij1345I(5,2))+s4
     -   5*(Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I(3,2)+Dij1345I(4,2)-D
     -   ij1345I(5,2))-s34*Dij1345I(5,2)+2*(Dij1345I(7,2)-Dij1345I(12
     -   ,3)+Dij1345I(13,3))+(D01345I+Dij1345I(1,1))*P(55)+Dij1345I(2
     -   ,2)*P(79)+Dij1345I(6,2)*P(80))
       F(143)=DCMPLX(FR(143),FI(143))
       P(81) = -p2sq-p3sq+p4sq+s12+s23-s34-s45+Is12*s34*P(16)
       P(82) = -p2sq-2*p3sq+p4sq+s12+s23
       P(83) = -s12+2*P(16)
       P(84) = s12*P(82)+s34*P(83)
       P(85) = p4sq-s34-s45
       P(86) = p4sq+s12-s34-s45
       P(87) = -p4sq+s34+s45
       P(88) = p2sq*P(85)+p3sq*P(86)+s23*P(87)
       P(89) = -(s45*P(2))+p4sq*P(16)
       P(90) = -p3sq+s34+s45
       P(91) = p3sq+p4sq-s34
       P(92) = p3sq-s34-s45
       P(93) = s12-s23-s34-s45
       P(94) = p3sq**2+p4sq*s12-s12*s34+s23*s34+s23*s45+p2sq*P(92)+p3s
     -   q*P(93)
       P(95) = s12+s23+s34+s45
       P(96) = p3sq**2+2*p4sq*s12+s12*s34+s23*s34+s45*P(48)+p2sq*P(92)
     -   -p3sq*P(95)
       P(97) = s12+3*s23
       P(98) = p3sq+3*p4sq-s34
       P(99) = 3*p4sq-s12-s23-s34
       P(100) = p3sq**2+s34*P(48)-p4sq*P(97)+p2sq*P(98)+p3sq*P(99)
       P(101) = 2-3*Is12*P(16)
       FR(144) = 4*(Cij134R(1,1)+D01345R*P(81)+Dij1345R(1,2)*P(90)+Is1
     -   2*(2*p4sq*Dij1345R(3,2)*P(1)+Dij1345R(1,1)*P(84)+Dij1345R(2,
     -   1)*P(88)-Dij1345R(3,1)*P(89)+P(16)*(C0134R-Cij134R(2,1)-Dij1
     -   345R(2,2)*P(91))+Dij1345R(4,2)*P(94)-Dij1345R(5,2)*P(96)+Dij
     -   1345R(6,2)*P(100))+2*Dij1345R(7,2)*P(101))
       FI(144) = 4*(Cij134I(1,1)+D01345I*P(81)+Dij1345I(1,2)*P(90)+Is1
     -   2*(2*p4sq*Dij1345I(3,2)*P(1)+Dij1345I(1,1)*P(84)+Dij1345I(2,
     -   1)*P(88)-Dij1345I(3,1)*P(89)+P(16)*(C0134I-Cij134I(2,1)-Dij1
     -   345I(2,2)*P(91))+Dij1345I(4,2)*P(94)-Dij1345I(5,2)*P(96)+Dij
     -   1345I(6,2)*P(100))+2*Dij1345I(7,2)*P(101))
       F(144)=DCMPLX(FR(144),FI(144))
       P(102) = -p4sq+s34
       P(103) = p4sq-s12+s34
       P(104) = p3sq-p4sq+s12-s34-s45
       FR(145) = -8*Is12*(C0134R+D01345R*s12-4*Dij1345R(7,2)-2*Dij1345
     -   R(13,3)+Dij1345R(1,1)*P(71)+Dij1345R(5,2)*P(90)-(Dij1345R(2,
     -   1)+Dij1345R(6,2))*P(102)-Dij1345R(3,1)*P(103)+Dij1345R(3,2)*
     -   P(104))
       FI(145) = -8*Is12*(C0134I+D01345I*s12-4*Dij1345I(7,2)-2*Dij1345
     -   I(13,3)+Dij1345I(1,1)*P(71)+Dij1345I(5,2)*P(90)-(Dij1345I(2,
     -   1)+Dij1345I(6,2))*P(102)-Dij1345I(3,1)*P(103)+Dij1345I(3,2)*
     -   P(104))
       F(145)=DCMPLX(FR(145),FI(145))
       P(105) = p3sq+p4sq+s12
       P(106) = p3sq+s12
       P(107) = p3sq+2*p4sq+s12
       FR(146) = 8*Is12*(C0134R+s34*Dij1345R(1,1)-s34*Dij1345R(2,1)-p4
     -   sq*Dij1345R(3,2)-4*Dij1345R(7,2)-2*(Dij1345R(12,3)-Dij1345R(
     -   13,3))+(Dij1345R(2,2)-Dij1345R(4,2)+Dij1345R(5,2)-Dij1345R(6
     -   ,2))*P(58)-Dij1345R(2,2)*P(105)+(Dij1345R(4,2)-Dij1345R(5,2)
     -   )*P(106)+Dij1345R(6,2)*P(107))
       FI(146) = 8*Is12*(C0134I+s34*Dij1345I(1,1)-s34*Dij1345I(2,1)-p4
     -   sq*Dij1345I(3,2)-4*Dij1345I(7,2)-2*(Dij1345I(12,3)-Dij1345I(
     -   13,3))+(Dij1345I(2,2)-Dij1345I(4,2)+Dij1345I(5,2)-Dij1345I(6
     -   ,2))*P(58)-Dij1345I(2,2)*P(105)+(Dij1345I(4,2)-Dij1345I(5,2)
     -   )*P(106)+Dij1345I(6,2)*P(107))
       F(146)=DCMPLX(FR(146),FI(146))
       P(108) = 2*s12-s34
       FR(147) = 8*Is12*(C0134R+Cij134R(1,1)+p4sq*Dij1345R(3,2)+s45*(D
     -   ij1345R(2,1)-Dij1345R(3,1)+Dij1345R(4,2))+s12*(Dij1345R(1,2)
     -   -Dij1345R(2,1)-Dij1345R(4,2)-Dij1345R(5,2))-p4sq*Dij1345R(5,
     -   2)-2*(Dij1345R(11,3)-Dij1345R(13,3))+Dij1345R(6,2)*P(7)+D013
     -   45R*P(55)+Dij1345R(1,1)*P(108))
       FI(147) = 8*Is12*(C0134I+Cij134I(1,1)+p4sq*Dij1345I(3,2)+s45*(D
     -   ij1345I(2,1)-Dij1345I(3,1)+Dij1345I(4,2))+s12*(Dij1345I(1,2)
     -   -Dij1345I(2,1)-Dij1345I(4,2)-Dij1345I(5,2))-p4sq*Dij1345I(5,
     -   2)-2*(Dij1345I(11,3)-Dij1345I(13,3))+Dij1345I(6,2)*P(7)+D013
     -   45I*P(55)+Dij1345I(1,1)*P(108))
       F(147)=DCMPLX(FR(147),FI(147))
       FR(148) = 8*Is12*(C0134R+D01345R*s12-s34*Dij1345R(2,1)+p4sq*(Di
     -   j1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,2))+(s34+s45)*Dij1345R(
     -   5,2)+p3sq*Dij1345R(6,2)-2*(Dij1345R(7,2)-Dij1345R(13,3))+Dij
     -   1345R(1,1)*P(71)+Dij1345R(6,2)*P(86)-Dij1345R(5,2)*P(106))
       FI(148) = 8*Is12*(C0134I+D01345I*s12-s34*Dij1345I(2,1)+p4sq*(Di
     -   j1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,2))+(s34+s45)*Dij1345I(
     -   5,2)+p3sq*Dij1345I(6,2)-2*(Dij1345I(7,2)-Dij1345I(13,3))+Dij
     -   1345I(1,1)*P(71)+Dij1345I(6,2)*P(86)-Dij1345I(5,2)*P(106))
       F(148)=DCMPLX(FR(148),FI(148))
       P(109) = p2sq+p3sq-p4sq-s12-s23+s34+s45-Is12*s34*P(16)
       P(110) = -p3sq-3*p4sq+s34
       P(111) = -3*p4sq+s12+s23+s34
       P(112) = -p3sq**2-s34*P(48)+p4sq*P(97)+p2sq*P(110)+p3sq*P(111)
       P(113) = -2*s12+3*P(16)
       FR(149) = -4*(Cij134R(1,1)+Dij1345R(1,2)*P(90)-D01345R*P(109)-I
     -   s12*(-(Dij1345R(1,1)*P(84))-Dij1345R(2,1)*P(88)+Dij1345R(3,1
     -   )*P(89)+P(16)*(-C0134R+Cij134R(2,1)+Dij1345R(2,2)*P(91))-Dij
     -   1345R(4,2)*P(94)+Dij1345R(5,2)*P(96)+Dij1345R(6,2)*P(112)+2*
     -   (p4sq*Dij1345R(3,2)*P(2)+Dij1345R(7,2)*P(113))))
       FI(149) = -4*(Cij134I(1,1)+Dij1345I(1,2)*P(90)-D01345I*P(109)-I
     -   s12*(-(Dij1345I(1,1)*P(84))-Dij1345I(2,1)*P(88)+Dij1345I(3,1
     -   )*P(89)+P(16)*(-C0134I+Cij134I(2,1)+Dij1345I(2,2)*P(91))-Dij
     -   1345I(4,2)*P(94)+Dij1345I(5,2)*P(96)+Dij1345I(6,2)*P(112)+2*
     -   (p4sq*Dij1345I(3,2)*P(2)+Dij1345I(7,2)*P(113))))
       F(149)=DCMPLX(FR(149),FI(149))
       FR(150) = -8*Is12*(C0134R+s34*Dij1345R(1,1)-s34*Dij1345R(2,1)-p
     -   4sq*Dij1345R(3,2)-4*Dij1345R(7,2)-2*(Dij1345R(12,3)-Dij1345R
     -   (13,3))+(Dij1345R(2,2)-Dij1345R(4,2)+Dij1345R(5,2)-Dij1345R(
     -   6,2))*P(58)-Dij1345R(2,2)*P(105)+(Dij1345R(4,2)-Dij1345R(5,2
     -   ))*P(106)+Dij1345R(6,2)*P(107))
       FI(150) = -8*Is12*(C0134I+s34*Dij1345I(1,1)-s34*Dij1345I(2,1)-p
     -   4sq*Dij1345I(3,2)-4*Dij1345I(7,2)-2*(Dij1345I(12,3)-Dij1345I
     -   (13,3))+(Dij1345I(2,2)-Dij1345I(4,2)+Dij1345I(5,2)-Dij1345I(
     -   6,2))*P(58)-Dij1345I(2,2)*P(105)+(Dij1345I(4,2)-Dij1345I(5,2
     -   ))*P(106)+Dij1345I(6,2)*P(107))
       F(150)=DCMPLX(FR(150),FI(150))
       FR(151) = -8*Is12*(C0134R+Cij134R(1,1)+p4sq*Dij1345R(3,2)+s45*(
     -   Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R(4,2))+s12*(Dij1345R(1,2
     -   )-Dij1345R(2,1)-Dij1345R(4,2)-Dij1345R(5,2))-p4sq*Dij1345R(5
     -   ,2)-2*(Dij1345R(11,3)-Dij1345R(13,3))+Dij1345R(6,2)*P(7)+D01
     -   345R*P(55)+Dij1345R(1,1)*P(108))
       FI(151) = -8*Is12*(C0134I+Cij134I(1,1)+p4sq*Dij1345I(3,2)+s45*(
     -   Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I(4,2))+s12*(Dij1345I(1,2
     -   )-Dij1345I(2,1)-Dij1345I(4,2)-Dij1345I(5,2))-p4sq*Dij1345I(5
     -   ,2)-2*(Dij1345I(11,3)-Dij1345I(13,3))+Dij1345I(6,2)*P(7)+D01
     -   345I*P(55)+Dij1345I(1,1)*P(108))
       F(151)=DCMPLX(FR(151),FI(151))
       FR(152) = -8*Is12*(C0134R+D01345R*s12-s34*Dij1345R(2,1)+p4sq*(D
     -   ij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,2))+(s34+s45)*Dij1345R
     -   (5,2)+p3sq*Dij1345R(6,2)-2*(Dij1345R(7,2)-Dij1345R(13,3))+Di
     -   j1345R(1,1)*P(71)+Dij1345R(6,2)*P(86)-Dij1345R(5,2)*P(106))
       FI(152) = -8*Is12*(C0134I+D01345I*s12-s34*Dij1345I(2,1)+p4sq*(D
     -   ij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,2))+(s34+s45)*Dij1345I
     -   (5,2)+p3sq*Dij1345I(6,2)-2*(Dij1345I(7,2)-Dij1345I(13,3))+Di
     -   j1345I(1,1)*P(71)+Dij1345I(6,2)*P(86)-Dij1345I(5,2)*P(106))
       F(152)=DCMPLX(FR(152),FI(152))
       FR(153) = 8*Is12*(p3sq*(Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R(4,
     -   2)-Dij1345R(6,2))+2*(Dij1345R(7,2)+Dij1345R(11,3)-Dij1345R(1
     -   3,3))+(Dij1345R(3,2)-Dij1345R(5,2))*P(78))
       FI(153) = 8*Is12*(p3sq*(Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I(4,
     -   2)-Dij1345I(6,2))+2*(Dij1345I(7,2)+Dij1345I(11,3)-Dij1345I(1
     -   3,3))+(Dij1345I(3,2)-Dij1345I(5,2))*P(78))
       F(153)=DCMPLX(FR(153),FI(153))
       P(114) = p3sq-s45
       FR(154) = -8*Is12*(Cij134R(2,1)+2*(Dij1345R(7,2)-Dij1345R(12,3)
     -   +Dij1345R(13,3))+(D01345R+Dij1345R(1,1))*P(55)+Dij1345R(6,2)
     -   *P(80)+(Dij1345R(4,2)-Dij1345R(5,2))*P(90)-Dij1345R(2,2)*P(1
     -   02)-Dij1345R(3,2)*P(104)+(-Dij1345R(2,1)+Dij1345R(3,1))*P(11
     -   4))
       FI(154) = -8*Is12*(Cij134I(2,1)+2*(Dij1345I(7,2)-Dij1345I(12,3)
     -   +Dij1345I(13,3))+(D01345I+Dij1345I(1,1))*P(55)+Dij1345I(6,2)
     -   *P(80)+(Dij1345I(4,2)-Dij1345I(5,2))*P(90)-Dij1345I(2,2)*P(1
     -   02)-Dij1345I(3,2)*P(104)+(-Dij1345I(2,1)+Dij1345I(3,1))*P(11
     -   4))
       F(154)=DCMPLX(FR(154),FI(154))
       FR(155) = 8*Is12*(C0134R+D01345R*s12-s34*Dij1345R(2,1)-s34*Dij1
     -   345R(3,1)+p3sq*Dij1345R(3,2)-s34*Dij1345R(3,2)-s45*Dij1345R(
     -   3,2)+s12*(Dij1345R(3,1)+Dij1345R(3,2))-p3sq*Dij1345R(5,2)+s3
     -   4*Dij1345R(5,2)+s45*Dij1345R(5,2)-s34*Dij1345R(6,2)+p4sq*(Di
     -   j1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,2)+Dij1345R(6,2))-4*Dij
     -   1345R(7,2)-2*Dij1345R(13,3)+Dij1345R(1,1)*P(71))
       FI(155) = 8*Is12*(C0134I+D01345I*s12-s34*Dij1345I(2,1)-s34*Dij1
     -   345I(3,1)+p3sq*Dij1345I(3,2)-s34*Dij1345I(3,2)-s45*Dij1345I(
     -   3,2)+s12*(Dij1345I(3,1)+Dij1345I(3,2))-p3sq*Dij1345I(5,2)+s3
     -   4*Dij1345I(5,2)+s45*Dij1345I(5,2)-s34*Dij1345I(6,2)+p4sq*(Di
     -   j1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,2)+Dij1345I(6,2))-4*Dij
     -   1345I(7,2)-2*Dij1345I(13,3)+Dij1345I(1,1)*P(71))
       F(155)=DCMPLX(FR(155),FI(155))
       FR(156) = -16*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,3)+D
     -   ij1345R(4,3)-Dij1345R(5,3)+Dij1345R(9,3)+2*(Dij1345R(3,2)+Di
     -   j1345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij134
     -   5R(10,3)))
       FI(156) = -16*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,3)+D
     -   ij1345I(4,3)-Dij1345I(5,3)+Dij1345I(9,3)+2*(Dij1345I(3,2)+Di
     -   j1345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij134
     -   5I(10,3)))
       F(156)=DCMPLX(FR(156),FI(156))
       FR(157) = Is12*(Dij1345R(2,2)+Dij1345R(3,2)-Dij1345R(3,3)+Dij13
     -   45R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)-2*(Dij1345R(6,2)-Dij134
     -   5R(9,3)+Dij1345R(10,3)))
       FI(157) = Is12*(Dij1345I(2,2)+Dij1345I(3,2)-Dij1345I(3,3)+Dij13
     -   45I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)-2*(Dij1345I(6,2)-Dij134
     -   5I(9,3)+Dij1345I(10,3)))
       F(157)=DCMPLX(FR(157),FI(157))
       FR(158) = Is12*(Dij1345R(3,3)-Dij1345R(5,2)+Dij1345R(6,2)-Dij13
     -   45R(7,3)-Dij1345R(9,3)+Dij1345R(10,3))
       FI(158) = Is12*(Dij1345I(3,3)-Dij1345I(5,2)+Dij1345I(6,2)-Dij13
     -   45I(7,3)-Dij1345I(9,3)+Dij1345I(10,3))
       F(158)=DCMPLX(FR(158),FI(158))
       FR(159) = Is12*(Dij1345R(2,1)+Dij1345R(2,2)-Dij1345R(3,1)-Dij13
     -   45R(3,3)+Dij1345R(4,2)-Dij1345R(5,2)-3*Dij1345R(6,2)+Dij1345
     -   R(6,3)+Dij1345R(7,3)-Dij1345R(8,3)+2*(Dij1345R(3,2)+Dij1345R
     -   (9,3)-Dij1345R(10,3)))
       FI(159) = Is12*(Dij1345I(2,1)+Dij1345I(2,2)-Dij1345I(3,1)-Dij13
     -   45I(3,3)+Dij1345I(4,2)-Dij1345I(5,2)-3*Dij1345I(6,2)+Dij1345
     -   I(6,3)+Dij1345I(7,3)-Dij1345I(8,3)+2*(Dij1345I(3,2)+Dij1345I
     -   (9,3)-Dij1345I(10,3)))
       F(159)=DCMPLX(FR(159),FI(159))
       FR(160) = Is12*(Dij1345R(2,2)+Dij1345R(2,3)+Dij1345R(3,2)-Dij13
     -   45R(3,3)-2*Dij1345R(6,2)-3*(Dij1345R(8,3)-Dij1345R(9,3)))
       FI(160) = Is12*(Dij1345I(2,2)+Dij1345I(2,3)+Dij1345I(3,2)-Dij13
     -   45I(3,3)-2*Dij1345I(6,2)-3*(Dij1345I(8,3)-Dij1345I(9,3)))
       F(160)=DCMPLX(FR(160),FI(160))
       FR(161) = Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,2)+Dij13
     -   45R(3,3)+Dij1345R(4,2)-Dij1345R(7,3)-Dij1345R(9,3)+Dij1345R(
     -   10,3))
       FI(161) = Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,2)+Dij13
     -   45I(3,3)+Dij1345I(4,2)-Dij1345I(7,3)-Dij1345I(9,3)+Dij1345I(
     -   10,3))
       F(161)=DCMPLX(FR(161),FI(161))
       FR(162) = Is12*(Dij1345R(2,1)+Dij1345R(2,2)-Dij1345R(3,1)+Dij13
     -   45R(3,3)-Dij1345R(6,2)+Dij1345R(8,3)-2*Dij1345R(9,3))
       FI(162) = Is12*(Dij1345I(2,1)+Dij1345I(2,2)-Dij1345I(3,1)+Dij13
     -   45I(3,3)-Dij1345I(6,2)+Dij1345I(8,3)-2*Dij1345I(9,3))
       F(162)=DCMPLX(FR(162),FI(162))
       FR(163) = Is12*(Dij1345R(3,2)+Dij1345R(3,3)-Dij1345R(6,2)-Dij13
     -   45R(9,3))
       FI(163) = Is12*(Dij1345I(3,2)+Dij1345I(3,3)-Dij1345I(6,2)-Dij13
     -   45I(9,3))
       F(163)=DCMPLX(FR(163),FI(163))
       FR(164) = -8*Is12*(2*(Dij1345R(7,2)+Dij1345R(12,3)-Dij1345R(13,
     -   3))+(D01345R+Dij1345R(1,1))*P(55))
       FI(164) = -8*Is12*(2*(Dij1345I(7,2)+Dij1345I(12,3)-Dij1345I(13,
     -   3))+(D01345I+Dij1345I(1,1))*P(55))
       F(164)=DCMPLX(FR(164),FI(164))
       FR(165) = -8*Is12*(C0134R-Cij134R(2,1)-s45*Dij1345R(2,1)-s34*(D
     -   ij1345R(2,1)-Dij1345R(2,2)+Dij1345R(4,2))-s12*(Dij1345R(1,2)
     -   -Dij1345R(3,1)+Dij1345R(3,2)-2*Dij1345R(5,2))+s45*(Dij1345R(
     -   3,1)-Dij1345R(4,2)+Dij1345R(5,2))+s34*(Dij1345R(5,2)-Dij1345
     -   R(6,2))+p3sq*(-Dij1345R(2,2)+Dij1345R(4,2)-Dij1345R(5,2)+Dij
     -   1345R(6,2))+p4sq*(Dij1345R(2,1)-Dij1345R(2,2)-Dij1345R(3,1)-
     -   2*Dij1345R(3,2)+3*Dij1345R(6,2))-6*Dij1345R(7,2)-(D01345R+2*
     -   Dij1345R(1,1))*P(55))
       FI(165) = -8*Is12*(C0134I-Cij134I(2,1)-s45*Dij1345I(2,1)-s34*(D
     -   ij1345I(2,1)-Dij1345I(2,2)+Dij1345I(4,2))-s12*(Dij1345I(1,2)
     -   -Dij1345I(3,1)+Dij1345I(3,2)-2*Dij1345I(5,2))+s45*(Dij1345I(
     -   3,1)-Dij1345I(4,2)+Dij1345I(5,2))+s34*(Dij1345I(5,2)-Dij1345
     -   I(6,2))+p3sq*(-Dij1345I(2,2)+Dij1345I(4,2)-Dij1345I(5,2)+Dij
     -   1345I(6,2))+p4sq*(Dij1345I(2,1)-Dij1345I(2,2)-Dij1345I(3,1)-
     -   2*Dij1345I(3,2)+3*Dij1345I(6,2))-6*Dij1345I(7,2)-(D01345I+2*
     -   Dij1345I(1,1))*P(55))
       F(165)=DCMPLX(FR(165),FI(165))
       FR(166) = Dij1345R(1,1)+Dij1345R(1,2)-Dij1345R(3,1)+Dij1345R(3,
     -   2)-2*Dij1345R(5,2)
       FI(166) = Dij1345I(1,1)+Dij1345I(1,2)-Dij1345I(3,1)+Dij1345I(3,
     -   2)-2*Dij1345I(5,2)
       F(166)=DCMPLX(FR(166),FI(166))
       FR(167) = Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R(3,2)+Dij1345R(4,
     -   2)-Dij1345R(5,2)-Dij1345R(6,2)
       FI(167) = Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I(3,2)+Dij1345I(4,
     -   2)-Dij1345I(5,2)-Dij1345I(6,2)
       F(167)=DCMPLX(FR(167),FI(167))
       FR(168) = D01345R+Dij1345R(1,1)+Dij1345R(3,2)-Dij1345R(5,2)
       FI(168) = D01345I+Dij1345I(1,1)+Dij1345I(3,2)-Dij1345I(5,2)
       F(168)=DCMPLX(FR(168),FI(168))
       FR(169) = 8*Is12*(C0134R-Cij134R(2,1)+s34*(-Dij1345R(2,1)+Dij13
     -   45R(2,2)-Dij1345R(4,2)+Dij1345R(5,2))+s45*(-Dij1345R(2,1)+Di
     -   j1345R(3,1)-Dij1345R(4,2)+Dij1345R(5,2))+s12*(-Dij1345R(1,2)
     -   +Dij1345R(2,1)+Dij1345R(4,2)+Dij1345R(5,2))+p3sq*(-Dij1345R(
     -   2,2)+Dij1345R(4,2)-Dij1345R(5,2)+Dij1345R(6,2))+p4sq*(Dij134
     -   5R(2,1)-Dij1345R(2,2)-Dij1345R(3,1)-2*Dij1345R(3,2)+3*Dij134
     -   5R(6,2))-6*Dij1345R(7,2)-(D01345R+2*Dij1345R(1,1))*P(55)-Dij
     -   1345R(6,2)*P(71))
       FI(169) = 8*Is12*(C0134I-Cij134I(2,1)+s34*(-Dij1345I(2,1)+Dij13
     -   45I(2,2)-Dij1345I(4,2)+Dij1345I(5,2))+s45*(-Dij1345I(2,1)+Di
     -   j1345I(3,1)-Dij1345I(4,2)+Dij1345I(5,2))+s12*(-Dij1345I(1,2)
     -   +Dij1345I(2,1)+Dij1345I(4,2)+Dij1345I(5,2))+p3sq*(-Dij1345I(
     -   2,2)+Dij1345I(4,2)-Dij1345I(5,2)+Dij1345I(6,2))+p4sq*(Dij134
     -   5I(2,1)-Dij1345I(2,2)-Dij1345I(3,1)-2*Dij1345I(3,2)+3*Dij134
     -   5I(6,2))-6*Dij1345I(7,2)-(D01345I+2*Dij1345I(1,1))*P(55)-Dij
     -   1345I(6,2)*P(71))
       F(169)=DCMPLX(FR(169),FI(169))
       FR(170) = Dij1345R(1,1)+Dij1345R(1,2)-Dij1345R(2,1)-Dij1345R(4,
     -   2)-Dij1345R(5,2)+Dij1345R(6,2)
       FI(170) = Dij1345I(1,1)+Dij1345I(1,2)-Dij1345I(2,1)-Dij1345I(4,
     -   2)-Dij1345I(5,2)+Dij1345I(6,2)
       F(170)=DCMPLX(FR(170),FI(170))
       FR(171) = Dij1345R(2,2)-Dij1345R(4,2)+Dij1345R(5,2)-Dij1345R(6,
     -   2)
       FI(171) = Dij1345I(2,2)-Dij1345I(4,2)+Dij1345I(5,2)-Dij1345I(6,
     -   2)
       F(171)=DCMPLX(FR(171),FI(171))
       FR(172) = D01345R+Dij1345R(1,1)-Dij1345R(5,2)+Dij1345R(6,2)
       FI(172) = D01345I+Dij1345I(1,1)-Dij1345I(5,2)+Dij1345I(6,2)
       F(172)=DCMPLX(FR(172),FI(172))
       FR(173) = 2*(Cij134R(1,1)+D01345R*P(81)+Dij1345R(1,2)*P(90)+Is1
     -   2*(2*p4sq*Dij1345R(3,2)*P(1)+Dij1345R(1,1)*P(84)+Dij1345R(2,
     -   1)*P(88)-Dij1345R(3,1)*P(89)+P(16)*(C0134R-Cij134R(2,1)-Dij1
     -   345R(2,2)*P(91))+Dij1345R(4,2)*P(94)-Dij1345R(5,2)*P(96)+Dij
     -   1345R(6,2)*P(100))+2*Dij1345R(7,2)*P(101))
       FI(173) = 2*(Cij134I(1,1)+D01345I*P(81)+Dij1345I(1,2)*P(90)+Is1
     -   2*(2*p4sq*Dij1345I(3,2)*P(1)+Dij1345I(1,1)*P(84)+Dij1345I(2,
     -   1)*P(88)-Dij1345I(3,1)*P(89)+P(16)*(C0134I-Cij134I(2,1)-Dij1
     -   345I(2,2)*P(91))+Dij1345I(4,2)*P(94)-Dij1345I(5,2)*P(96)+Dij
     -   1345I(6,2)*P(100))+2*Dij1345I(7,2)*P(101))
       F(173)=DCMPLX(FR(173),FI(173))
       FR(174) = -4*Is12*(p3sq*(Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R(4
     -   ,2)-Dij1345R(6,2))+2*(Dij1345R(7,2)+Dij1345R(11,3)-Dij1345R(
     -   13,3))+(Dij1345R(3,2)-Dij1345R(5,2))*P(78))
       FI(174) = -4*Is12*(p3sq*(Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I(4
     -   ,2)-Dij1345I(6,2))+2*(Dij1345I(7,2)+Dij1345I(11,3)-Dij1345I(
     -   13,3))+(Dij1345I(3,2)-Dij1345I(5,2))*P(78))
       F(174)=DCMPLX(FR(174),FI(174))
       FR(175) = 4*Is12*(Cij134R(2,1)+p4sq*Dij1345R(3,2)-s12*Dij1345R(
     -   3,2)+s34*Dij1345R(3,2)+s34*Dij1345R(4,2)-p3sq*(Dij1345R(2,1)
     -   -Dij1345R(3,1)+Dij1345R(3,2)+Dij1345R(4,2)-Dij1345R(5,2))+s4
     -   5*(Dij1345R(2,1)-Dij1345R(3,1)+Dij1345R(3,2)+Dij1345R(4,2)-D
     -   ij1345R(5,2))-s34*Dij1345R(5,2)+2*(Dij1345R(7,2)-Dij1345R(12
     -   ,3)+Dij1345R(13,3))+(D01345R+Dij1345R(1,1))*P(55)+Dij1345R(2
     -   ,2)*P(79)+Dij1345R(6,2)*P(80))
       FI(175) = 4*Is12*(Cij134I(2,1)+p4sq*Dij1345I(3,2)-s12*Dij1345I(
     -   3,2)+s34*Dij1345I(3,2)+s34*Dij1345I(4,2)-p3sq*(Dij1345I(2,1)
     -   -Dij1345I(3,1)+Dij1345I(3,2)+Dij1345I(4,2)-Dij1345I(5,2))+s4
     -   5*(Dij1345I(2,1)-Dij1345I(3,1)+Dij1345I(3,2)+Dij1345I(4,2)-D
     -   ij1345I(5,2))-s34*Dij1345I(5,2)+2*(Dij1345I(7,2)-Dij1345I(12
     -   ,3)+Dij1345I(13,3))+(D01345I+Dij1345I(1,1))*P(55)+Dij1345I(2
     -   ,2)*P(79)+Dij1345I(6,2)*P(80))
       F(175)=DCMPLX(FR(175),FI(175))
       FR(176) = -4*Is12*(C0134R+D01345R*s12-4*Dij1345R(7,2)-2*Dij1345
     -   R(13,3)+Dij1345R(1,1)*P(71)+Dij1345R(5,2)*P(90)-(Dij1345R(2,
     -   1)+Dij1345R(6,2))*P(102)-Dij1345R(3,1)*P(103)+Dij1345R(3,2)*
     -   P(104))
       FI(176) = -4*Is12*(C0134I+D01345I*s12-4*Dij1345I(7,2)-2*Dij1345
     -   I(13,3)+Dij1345I(1,1)*P(71)+Dij1345I(5,2)*P(90)-(Dij1345I(2,
     -   1)+Dij1345I(6,2))*P(102)-Dij1345I(3,1)*P(103)+Dij1345I(3,2)*
     -   P(104))
       F(176)=DCMPLX(FR(176),FI(176))
       FR(177) = 4*Is12*(C0134R+s34*Dij1345R(1,1)-s34*Dij1345R(2,1)-p4
     -   sq*Dij1345R(3,2)-4*Dij1345R(7,2)-2*(Dij1345R(12,3)-Dij1345R(
     -   13,3))+(Dij1345R(2,2)-Dij1345R(4,2)+Dij1345R(5,2)-Dij1345R(6
     -   ,2))*P(58)-Dij1345R(2,2)*P(105)+(Dij1345R(4,2)-Dij1345R(5,2)
     -   )*P(106)+Dij1345R(6,2)*P(107))
       FI(177) = 4*Is12*(C0134I+s34*Dij1345I(1,1)-s34*Dij1345I(2,1)-p4
     -   sq*Dij1345I(3,2)-4*Dij1345I(7,2)-2*(Dij1345I(12,3)-Dij1345I(
     -   13,3))+(Dij1345I(2,2)-Dij1345I(4,2)+Dij1345I(5,2)-Dij1345I(6
     -   ,2))*P(58)-Dij1345I(2,2)*P(105)+(Dij1345I(4,2)-Dij1345I(5,2)
     -   )*P(106)+Dij1345I(6,2)*P(107))
       F(177)=DCMPLX(FR(177),FI(177))
       FR(178) = 4*Is12*(C0134R+Cij134R(1,1)+p4sq*Dij1345R(3,2)+s45*(D
     -   ij1345R(2,1)-Dij1345R(3,1)+Dij1345R(4,2))+s12*(Dij1345R(1,2)
     -   -Dij1345R(2,1)-Dij1345R(4,2)-Dij1345R(5,2))-p4sq*Dij1345R(5,
     -   2)-2*(Dij1345R(11,3)-Dij1345R(13,3))+Dij1345R(6,2)*P(7)+D013
     -   45R*P(55)+Dij1345R(1,1)*P(108))
       FI(178) = 4*Is12*(C0134I+Cij134I(1,1)+p4sq*Dij1345I(3,2)+s45*(D
     -   ij1345I(2,1)-Dij1345I(3,1)+Dij1345I(4,2))+s12*(Dij1345I(1,2)
     -   -Dij1345I(2,1)-Dij1345I(4,2)-Dij1345I(5,2))-p4sq*Dij1345I(5,
     -   2)-2*(Dij1345I(11,3)-Dij1345I(13,3))+Dij1345I(6,2)*P(7)+D013
     -   45I*P(55)+Dij1345I(1,1)*P(108))
       F(178)=DCMPLX(FR(178),FI(178))
       FR(179) = 4*Is12*(C0134R+D01345R*s12-s34*Dij1345R(2,1)+p4sq*(Di
     -   j1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,2))+(s34+s45)*Dij1345R(
     -   5,2)+p3sq*Dij1345R(6,2)-2*(Dij1345R(7,2)-Dij1345R(13,3))+Dij
     -   1345R(1,1)*P(71)+Dij1345R(6,2)*P(86)-Dij1345R(5,2)*P(106))
       FI(179) = 4*Is12*(C0134I+D01345I*s12-s34*Dij1345I(2,1)+p4sq*(Di
     -   j1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,2))+(s34+s45)*Dij1345I(
     -   5,2)+p3sq*Dij1345I(6,2)-2*(Dij1345I(7,2)-Dij1345I(13,3))+Dij
     -   1345I(1,1)*P(71)+Dij1345I(6,2)*P(86)-Dij1345I(5,2)*P(106))
       F(179)=DCMPLX(FR(179),FI(179))
       FR(180) = 8*Is12*(Dij1345R(2,1)-Dij1345R(3,1)-Dij1345R(3,3)+Dij
     -   1345R(4,3)-Dij1345R(5,3)+Dij1345R(9,3)+2*(Dij1345R(3,2)+Dij1
     -   345R(4,2)-Dij1345R(5,2)-Dij1345R(6,2)+Dij1345R(7,3)-Dij1345R
     -   (10,3)))
       FI(180) = 8*Is12*(Dij1345I(2,1)-Dij1345I(3,1)-Dij1345I(3,3)+Dij
     -   1345I(4,3)-Dij1345I(5,3)+Dij1345I(9,3)+2*(Dij1345I(3,2)+Dij1
     -   345I(4,2)-Dij1345I(5,2)-Dij1345I(6,2)+Dij1345I(7,3)-Dij1345I
     -   (10,3)))
       F(180)=DCMPLX(FR(180),FI(180))
       FR(181) = Is12*(8*(Dij1345R(7,2)+Dij1345R(12,3)-Dij1345R(13,3))
     -   +4*(D01345R+Dij1345R(1,1))*P(55))
       FI(181) = Is12*(8*(Dij1345I(7,2)+Dij1345I(12,3)-Dij1345I(13,3))
     -   +4*(D01345I+Dij1345I(1,1))*P(55))
       F(181)=DCMPLX(FR(181),FI(181))
       P(115) = s12-2*s34
       FR(182) = 4*Is12*(C0134R-Cij134R(2,1)+s45*(-Dij1345R(2,1)+Dij13
     -   45R(3,1)-Dij1345R(4,2)+Dij1345R(5,2))+p3sq*(-Dij1345R(2,2)+D
     -   ij1345R(4,2)-Dij1345R(5,2)+Dij1345R(6,2))-s34*(Dij1345R(2,1)
     -   -Dij1345R(2,2)+Dij1345R(4,2)-Dij1345R(5,2)+Dij1345R(6,2))+p4
     -   sq*(Dij1345R(2,1)-Dij1345R(2,2)-Dij1345R(3,1)-2*Dij1345R(3,2
     -   )+3*Dij1345R(6,2))-6*Dij1345R(7,2)-D01345R*P(55)-Dij1345R(1,
     -   1)*P(115))
       FI(182) = 4*Is12*(C0134I-Cij134I(2,1)+s45*(-Dij1345I(2,1)+Dij13
     -   45I(3,1)-Dij1345I(4,2)+Dij1345I(5,2))+p3sq*(-Dij1345I(2,2)+D
     -   ij1345I(4,2)-Dij1345I(5,2)+Dij1345I(6,2))-s34*(Dij1345I(2,1)
     -   -Dij1345I(2,2)+Dij1345I(4,2)-Dij1345I(5,2)+Dij1345I(6,2))+p4
     -   sq*(Dij1345I(2,1)-Dij1345I(2,2)-Dij1345I(3,1)-2*Dij1345I(3,2
     -   )+3*Dij1345I(6,2))-6*Dij1345I(7,2)-D01345I*P(55)-Dij1345I(1,
     -   1)*P(115))
       F(182)=DCMPLX(FR(182),FI(182))
       FR(183) = -4*Is12*(C0134R-Cij134R(2,1)+s45*(-Dij1345R(2,1)+Dij1
     -   345R(3,1)-Dij1345R(4,2)+Dij1345R(5,2))+p3sq*(-Dij1345R(2,2)+
     -   Dij1345R(4,2)-Dij1345R(5,2)+Dij1345R(6,2))-s34*(Dij1345R(2,1
     -   )-Dij1345R(2,2)+Dij1345R(4,2)-Dij1345R(5,2)+Dij1345R(6,2))+p
     -   4sq*(Dij1345R(2,1)-Dij1345R(2,2)-Dij1345R(3,1)-2*Dij1345R(3,
     -   2)+3*Dij1345R(6,2))-6*Dij1345R(7,2)-D01345R*P(55)-Dij1345R(1
     -   ,1)*P(115))
       FI(183) = -4*Is12*(C0134I-Cij134I(2,1)+s45*(-Dij1345I(2,1)+Dij1
     -   345I(3,1)-Dij1345I(4,2)+Dij1345I(5,2))+p3sq*(-Dij1345I(2,2)+
     -   Dij1345I(4,2)-Dij1345I(5,2)+Dij1345I(6,2))-s34*(Dij1345I(2,1
     -   )-Dij1345I(2,2)+Dij1345I(4,2)-Dij1345I(5,2)+Dij1345I(6,2))+p
     -   4sq*(Dij1345I(2,1)-Dij1345I(2,2)-Dij1345I(3,1)-2*Dij1345I(3,
     -   2)+3*Dij1345I(6,2))-6*Dij1345I(7,2)-D01345I*P(55)-Dij1345I(1
     -   ,1)*P(115))
       F(183)=DCMPLX(FR(183),FI(183))
       P(116) = 2*p4sq+s12
       P(117) = -p3sq+2*s12+s34+s45
       FR(184) = 4*Is12*(C0134R-Cij134R(2,1)-s12*Dij1345R(1,2)+p4sq*Di
     -   j1345R(2,1)-p4sq*Dij1345R(3,1)+(s12+s45)*Dij1345R(3,1)-6*Dij
     -   1345R(7,2)-(D01345R+2*Dij1345R(1,1))*P(55)-Dij1345R(2,1)*P(5
     -   8)-Dij1345R(4,2)*P(90)-Dij1345R(2,2)*P(91)+Dij1345R(6,2)*P(9
     -   8)-Dij1345R(3,2)*P(116)+Dij1345R(5,2)*P(117))
       FI(184) = 4*Is12*(C0134I-Cij134I(2,1)-s12*Dij1345I(1,2)+p4sq*Di
     -   j1345I(2,1)-p4sq*Dij1345I(3,1)+(s12+s45)*Dij1345I(3,1)-6*Dij
     -   1345I(7,2)-(D01345I+2*Dij1345I(1,1))*P(55)-Dij1345I(2,1)*P(5
     -   8)-Dij1345I(4,2)*P(90)-Dij1345I(2,2)*P(91)+Dij1345I(6,2)*P(9
     -   8)-Dij1345I(3,2)*P(116)+Dij1345I(5,2)*P(117))
       F(184)=DCMPLX(FR(184),FI(184))
       FR(185) = -4*Is12*(C0134R-Cij134R(2,1)+s34*(-Dij1345R(2,1)+Dij1
     -   345R(2,2)-Dij1345R(4,2)+Dij1345R(5,2))+s45*(-Dij1345R(2,1)+D
     -   ij1345R(3,1)-Dij1345R(4,2)+Dij1345R(5,2))+s12*(-Dij1345R(1,2
     -   )+Dij1345R(2,1)+Dij1345R(4,2)+Dij1345R(5,2))+p3sq*(-Dij1345R
     -   (2,2)+Dij1345R(4,2)-Dij1345R(5,2)+Dij1345R(6,2))+p4sq*(Dij13
     -   45R(2,1)-Dij1345R(2,2)-Dij1345R(3,1)-2*Dij1345R(3,2)+3*Dij13
     -   45R(6,2))-6*Dij1345R(7,2)-(D01345R+2*Dij1345R(1,1))*P(55)-Di
     -   j1345R(6,2)*P(71))
       FI(185) = -4*Is12*(C0134I-Cij134I(2,1)+s34*(-Dij1345I(2,1)+Dij1
     -   345I(2,2)-Dij1345I(4,2)+Dij1345I(5,2))+s45*(-Dij1345I(2,1)+D
     -   ij1345I(3,1)-Dij1345I(4,2)+Dij1345I(5,2))+s12*(-Dij1345I(1,2
     -   )+Dij1345I(2,1)+Dij1345I(4,2)+Dij1345I(5,2))+p3sq*(-Dij1345I
     -   (2,2)+Dij1345I(4,2)-Dij1345I(5,2)+Dij1345I(6,2))+p4sq*(Dij13
     -   45I(2,1)-Dij1345I(2,2)-Dij1345I(3,1)-2*Dij1345I(3,2)+3*Dij13
     -   45I(6,2))-6*Dij1345I(7,2)-(D01345I+2*Dij1345I(1,1))*P(55)-Di
     -   j1345I(6,2)*P(71))
       F(185)=DCMPLX(FR(185),FI(185))
       FR(186) = 2*Is12*(C0134R-Cij134R(2,1)+s45*(-Dij1345R(2,1)+Dij13
     -   45R(3,1)-Dij1345R(4,2)+Dij1345R(5,2))+p3sq*(-Dij1345R(2,2)+D
     -   ij1345R(4,2)-Dij1345R(5,2)+Dij1345R(6,2))-s34*(Dij1345R(2,1)
     -   -Dij1345R(2,2)+Dij1345R(4,2)-Dij1345R(5,2)+Dij1345R(6,2))+p4
     -   sq*(Dij1345R(2,1)-Dij1345R(2,2)-Dij1345R(3,1)-2*Dij1345R(3,2
     -   )+3*Dij1345R(6,2))-6*Dij1345R(7,2)-D01345R*P(55)-Dij1345R(1,
     -   1)*P(115))
       FI(186) = 2*Is12*(C0134I-Cij134I(2,1)+s45*(-Dij1345I(2,1)+Dij13
     -   45I(3,1)-Dij1345I(4,2)+Dij1345I(5,2))+p3sq*(-Dij1345I(2,2)+D
     -   ij1345I(4,2)-Dij1345I(5,2)+Dij1345I(6,2))-s34*(Dij1345I(2,1)
     -   -Dij1345I(2,2)+Dij1345I(4,2)-Dij1345I(5,2)+Dij1345I(6,2))+p4
     -   sq*(Dij1345I(2,1)-Dij1345I(2,2)-Dij1345I(3,1)-2*Dij1345I(3,2
     -   )+3*Dij1345I(6,2))-6*Dij1345I(7,2)-D01345I*P(55)-Dij1345I(1,
     -   1)*P(115))
       F(186)=DCMPLX(FR(186),FI(186))
       P(118) = p2sq-p3sq-s12-s15+s23+s34
       P(119) = p2sq-s12-s15+s34
       FR(187) = 8*(-(p2sq*(EijR(2,1)+EijR(5,3)+EijR(15,3)+2*(EijR(5,2
     -   )-EijR(9,2)-EijR(18,3))))-2*(EijR(11,2)+EijR(36,4)+EijR(39,4
     -   )+2*(EijR(21,3)-EijR(24,3)-EijR(42,4)))+(EijR(3,1)+EijR(6,2)
     -   -EijR(10,2))*P(16)-EijR(4,1)*P(18)-(EijR(4,2)-EijR(7,2))*P(1
     -   18)+(EijR(4,3)+EijR(7,3)-2*EijR(14,3))*P(119))
       FI(187) = 8*(-(p2sq*(EijI(2,1)+EijI(5,3)+EijI(15,3)+2*(EijI(5,2
     -   )-EijI(9,2)-EijI(18,3))))-2*(EijI(11,2)+EijI(36,4)+EijI(39,4
     -   )+2*(EijI(21,3)-EijI(24,3)-EijI(42,4)))+(EijI(3,1)+EijI(6,2)
     -   -EijI(10,2))*P(16)-EijI(4,1)*P(18)-(EijI(4,2)-EijI(7,2))*P(1
     -   18)+(EijI(4,3)+EijI(7,3)-2*EijI(14,3))*P(119))
       F(187)=DCMPLX(FR(187),FI(187))
       P(120) = p2sq+p3sq-s15+s23-s34+2*P(11)
       P(121) = p2sq+p3sq-s12-s15-s45
       P(122) = p3sq+p4sq-s45
       P(123) = p2sq-s15-s34+2*P(122)
       P(124) = 2*p2sq-p3sq-s12+s23
       P(125) = -p2sq+p3sq+s23
       P(126) = p3sq-s15-s45+2*P(9)
       P(127) = p2sq+p3sq-2*s12-s15-s45
       P(128) = p3sq+s12-s45
       P(129) = 2*p2sq+p3sq-3*s12-s15-s45
       FR(188) = -8*(Dij2345R(1,1)-Dij2345R(3,1)-Dij2345R(3,2)+Dij2345
     -   R(5,2)+p2sq*(EijR(2,1)+EijR(2,2))-s23*EijR(4,1)+s12*(-EijR(8
     -   ,3)+EijR(10,3))+2*(p3sq*EijR(6,2)+2*EijR(11,2)+EijR(39,4)+Ei
     -   jR(40,4)-EijR(42,4)-EijR(44,4))+(EijR(5,2)+EijR(5,3)-EijR(7,
     -   3))*P(9)-EijR(8,2)*P(16)-EijR(3,1)*P(35)+EijR(4,2)*P(120)-Ei
     -   jR(4,3)*P(121)-EijR(7,2)*P(123)-EijR(9,2)*P(124)-EijR(10,2)*
     -   P(125)+EijR(14,3)*P(126)+EijR(15,3)*P(127)+(EijR(16,3)+EijR(
     -   17,3)-EijR(19,3)-EijR(20,3))*P(128)-EijR(18,3)*P(129))
       FI(188) = -8*(Dij2345I(1,1)-Dij2345I(3,1)-Dij2345I(3,2)+Dij2345
     -   I(5,2)+p2sq*(EijI(2,1)+EijI(2,2))-s23*EijI(4,1)+s12*(-EijI(8
     -   ,3)+EijI(10,3))+2*(p3sq*EijI(6,2)+2*EijI(11,2)+EijI(39,4)+Ei
     -   jI(40,4)-EijI(42,4)-EijI(44,4))+(EijI(5,2)+EijI(5,3)-EijI(7,
     -   3))*P(9)-EijI(8,2)*P(16)-EijI(3,1)*P(35)+EijI(4,2)*P(120)-Ei
     -   jI(4,3)*P(121)-EijI(7,2)*P(123)-EijI(9,2)*P(124)-EijI(10,2)*
     -   P(125)+EijI(14,3)*P(126)+EijI(15,3)*P(127)+(EijI(16,3)+EijI(
     -   17,3)-EijI(19,3)-EijI(20,3))*P(128)-EijI(18,3)*P(129))
       F(188)=DCMPLX(FR(188),FI(188))
       P(130) = p4sq+s12+s15-s34-s45
       P(131) = p2sq-p4sq-s23+s34-2*P(128)
       P(132) = p2sq+p3sq+p4sq-s34
       P(133) = -p2sq-p3sq+s12+s15+s45
       P(134) = s15+2*P(14)
       P(135) = p2sq-s12+s15-2*P(13)
       P(136) = p3sq+s12-s23
       P(137) = p2sq-p3sq+s23
       P(138) = 3*p3sq+s12+s23+2*P(85)
       P(139) = p2sq-s15+2*P(114)
       P(140) = p2sq-s12-s15+2*P(114)
       FR(189) = 8*(D02345R+Dij2345R(3,2)+Dij2345R(4,2)-Dij2345R(5,2)-
     -   Dij2345R(6,2)-p2sq*EijR(2,2)+s12*EijR(9,3)-s12*EijR(10,3)-8*
     -   EijR(11,2)-2*(EijR(22,3)-EijR(23,3)+EijR(39,4)+EijR(43,4)-Ei
     -   jR(44,4)-EijR(45,4))+(EijR(5,2)-EijR(14,3)-EijR(17,3)+EijR(1
     -   8,3)+EijR(19,3))*P(9)-EijR(15,3)*P(127)+(-EijR(12,3)+EijR(13
     -   ,3))*P(128)+(-EijR(3,1)+EijR(4,1))*P(130)+EijR(3,2)*P(131)-E
     -   ijR(4,2)*P(132)-EijR(4,3)*P(133)-EijR(6,2)*P(134)+EijR(7,2)*
     -   P(135)+EijR(8,2)*P(136)+EijR(9,2)*P(137)+EijR(10,2)*P(138)-E
     -   ijR(16,3)*P(139)+EijR(20,3)*P(140))
       FI(189) = 8*(D02345I+Dij2345I(3,2)+Dij2345I(4,2)-Dij2345I(5,2)-
     -   Dij2345I(6,2)-p2sq*EijI(2,2)+s12*EijI(9,3)-s12*EijI(10,3)-8*
     -   EijI(11,2)-2*(EijI(22,3)-EijI(23,3)+EijI(39,4)+EijI(43,4)-Ei
     -   jI(44,4)-EijI(45,4))+(EijI(5,2)-EijI(14,3)-EijI(17,3)+EijI(1
     -   8,3)+EijI(19,3))*P(9)-EijI(15,3)*P(127)+(-EijI(12,3)+EijI(13
     -   ,3))*P(128)+(-EijI(3,1)+EijI(4,1))*P(130)+EijI(3,2)*P(131)-E
     -   ijI(4,2)*P(132)-EijI(4,3)*P(133)-EijI(6,2)*P(134)+EijI(7,2)*
     -   P(135)+EijI(8,2)*P(136)+EijI(9,2)*P(137)+EijI(10,2)*P(138)-E
     -   ijI(16,3)*P(139)+EijI(20,3)*P(140))
       F(189)=DCMPLX(FR(189),FI(189))
       P(141) = s12+s15-s34
       P(142) = p3sq+2*s12+s15-s34-s45
       P(143) = p2sq+p3sq-s34-s45
       P(144) = p3sq+s12-s34-s45
       P(145) = p2sq+s15+2*P(144)
       P(146) = p2sq+s12-s15+3*P(114)
       P(147) = p2sq-s15
       P(148) = -s12+3*P(114)+2*P(147)
       FR(190) = 8*(Dij2345R(2,1)+Dij2345R(2,2)-Dij2345R(3,1)+Dij2345R
     -   (3,2)+s12*(EijR(8,2)-EijR(9,2)+EijR(12,3)+EijR(15,3)-2*EijR(
     -   20,3))-2*(Dij2345R(6,2)+EijR(23,3)-EijR(24,3)+EijR(38,4)+Eij
     -   R(39,4)-2*EijR(45,4))-(EijR(6,2)-EijR(7,2)+EijR(11,3)+EijR(1
     -   4,3)-2*EijR(19,3))*P(9)-EijR(3,3)*P(128)-EijR(4,3)*P(133)+(-
     -   EijR(3,1)+EijR(4,1))*P(141)-EijR(3,2)*P(142)-EijR(4,2)*P(143
     -   )+EijR(10,2)*P(145)+EijR(13,3)*P(146)-EijR(16,3)*P(148))
       FI(190) = 8*(Dij2345I(2,1)+Dij2345I(2,2)-Dij2345I(3,1)+Dij2345I
     -   (3,2)+s12*(EijI(8,2)-EijI(9,2)+EijI(12,3)+EijI(15,3)-2*EijI(
     -   20,3))-2*(Dij2345I(6,2)+EijI(23,3)-EijI(24,3)+EijI(38,4)+Eij
     -   I(39,4)-2*EijI(45,4))-(EijI(6,2)-EijI(7,2)+EijI(11,3)+EijI(1
     -   4,3)-2*EijI(19,3))*P(9)-EijI(3,3)*P(128)-EijI(4,3)*P(133)+(-
     -   EijI(3,1)+EijI(4,1))*P(141)-EijI(3,2)*P(142)-EijI(4,2)*P(143
     -   )+EijI(10,2)*P(145)+EijI(13,3)*P(146)-EijI(16,3)*P(148))
       F(190)=DCMPLX(FR(190),FI(190))
       P(149) = p2sq+p3sq+s15-s23-s34
       FR(191) = 8*(Dij2345R(3,2)-Dij2345R(6,2)+p2sq*EijR(8,2)+s12*(Ei
     -   jR(15,3)-EijR(20,3))+2*(EijR(23,3)-EijR(39,4)+EijR(45,4))+(-
     -   EijR(3,1)+EijR(4,1))*P(2)-(EijR(6,2)+EijR(14,3)-EijR(19,3))*
     -   P(9)-EijR(3,2)*P(16)+EijR(4,3)*P(121)+EijR(13,3)*P(128)-EijR
     -   (16,3)*P(139)-EijR(4,2)*P(141)+EijR(10,2)*P(149))
       FI(191) = 8*(Dij2345I(3,2)-Dij2345I(6,2)+p2sq*EijI(8,2)+s12*(Ei
     -   jI(15,3)-EijI(20,3))+2*(EijI(23,3)-EijI(39,4)+EijI(45,4))+(-
     -   EijI(3,1)+EijI(4,1))*P(2)-(EijI(6,2)+EijI(14,3)-EijI(19,3))*
     -   P(9)-EijI(3,2)*P(16)+EijI(4,3)*P(121)+EijI(13,3)*P(128)-EijI
     -   (16,3)*P(139)-EijI(4,2)*P(141)+EijI(10,2)*P(149))
       F(191)=DCMPLX(FR(191),FI(191))
       P(150) = p2sq+p3sq-s12-s15-s23+s34
       P(151) = p3sq-s12-s15-s23+s34
       P(152) = p2sq+s12
       P(153) = 2*p3sq-s15-s23+s34-s45
       P(154) = p2sq-2*p3sq-s12+s23+s34+s45
       P(155) = s12+s23-s34-s45
       FR(192) = 8*(D02345R+Dij2345R(2,1)+s12*EijR(8,2)+p2sq*(-EijR(2,
     -   1)+EijR(9,2)-EijR(15,3)-EijR(17,3)+EijR(18,3)+EijR(20,3))-2*
     -   (2*EijR(11,2)+EijR(23,3)-EijR(24,3)+EijR(39,4)+EijR(41,4)-Ei
     -   jR(42,4)-EijR(45,4))-EijR(4,2)*P(118)+(EijR(4,3)-EijR(14,3)-
     -   EijR(16,3)+EijR(19,3))*P(119)-EijR(3,2)*P(128)+EijR(3,1)*P(1
     -   50)-EijR(4,1)*P(151)-EijR(5,2)*P(152)+EijR(6,2)*P(153)+EijR(
     -   7,2)*P(154)+EijR(10,2)*P(155))
       FI(192) = 8*(D02345I+Dij2345I(2,1)+s12*EijI(8,2)+p2sq*(-EijI(2,
     -   1)+EijI(9,2)-EijI(15,3)-EijI(17,3)+EijI(18,3)+EijI(20,3))-2*
     -   (2*EijI(11,2)+EijI(23,3)-EijI(24,3)+EijI(39,4)+EijI(41,4)-Ei
     -   jI(42,4)-EijI(45,4))-EijI(4,2)*P(118)+(EijI(4,3)-EijI(14,3)-
     -   EijI(16,3)+EijI(19,3))*P(119)-EijI(3,2)*P(128)+EijI(3,1)*P(1
     -   50)-EijI(4,1)*P(151)-EijI(5,2)*P(152)+EijI(6,2)*P(153)+EijI(
     -   7,2)*P(154)+EijI(10,2)*P(155))
       F(192)=DCMPLX(FR(192),FI(192))
       P(156) = p2sq+2*p3sq-s23-s45
       P(157) = -p2sq-2*p3sq+s23+s34+s45
       P(158) = p2sq+p3sq-s45
       P(159) = -s12-s15+2*P(158)
       FR(193) = -8*(D02345R+Dij2345R(2,1)-Dij2345R(3,2)+Dij2345R(6,2)
     -   +p2sq*EijR(8,2)-s12*(EijR(5,2)+EijR(15,3)+EijR(17,3)-EijR(18
     -   ,3)-EijR(20,3))-2*(EijR(11,2)-EijR(21,3)+EijR(24,3)-EijR(39,
     -   4)-EijR(41,4)+EijR(42,4)+EijR(45,4))+(EijR(6,3)-EijR(7,3)-Ei
     -   jR(9,2))*P(9)-EijR(10,2)*P(13)-EijR(3,2)*P(16)+EijR(14,3)*P(
     -   126)+(EijR(11,3)-EijR(13,3))*P(128)+EijR(4,3)*P(133)+EijR(16
     -   ,3)*P(139)+EijR(4,2)*P(143)+EijR(6,2)*P(156)+EijR(7,2)*P(157
     -   )-EijR(19,3)*P(159))
       FI(193) = -8*(D02345I+Dij2345I(2,1)-Dij2345I(3,2)+Dij2345I(6,2)
     -   +p2sq*EijI(8,2)-s12*(EijI(5,2)+EijI(15,3)+EijI(17,3)-EijI(18
     -   ,3)-EijI(20,3))-2*(EijI(11,2)-EijI(21,3)+EijI(24,3)-EijI(39,
     -   4)-EijI(41,4)+EijI(42,4)+EijI(45,4))+(EijI(6,3)-EijI(7,3)-Ei
     -   jI(9,2))*P(9)-EijI(10,2)*P(13)-EijI(3,2)*P(16)+EijI(14,3)*P(
     -   126)+(EijI(11,3)-EijI(13,3))*P(128)+EijI(4,3)*P(133)+EijI(16
     -   ,3)*P(139)+EijI(4,2)*P(143)+EijI(6,2)*P(156)+EijI(7,2)*P(157
     -   )-EijI(19,3)*P(159))
       F(193)=DCMPLX(FR(193),FI(193))
       P(160) = p2sq+p3sq-p4sq-s23+s34-2*P(7)
       P(161) = -p3sq+p4sq+s12+s15+s23-s34-2*s45
       P(162) = p2sq-p3sq-p4sq-s12-s15+s23+s34-s45
       P(163) = p2sq-p4sq-s23+s34+2*P(38)
       P(164) = p2sq-2*p3sq-p4sq-s12+s15+s23+s34
       P(165) = p2sq+p3sq-p4sq-s23+s34-3*P(7)
       P(166) = 2*p2sq-s12-s15+s34
       FR(194) = 8*(Dij2345R(1,1)-Dij2345R(3,1)-s15*(EE0R+2*EijR(1,1)+
     -   EijR(1,2))+s12*EijR(2,2)+p2sq*(-EijR(8,3)+EijR(10,3))+2*(2*E
     -   ijR(11,2)-EijR(22,3)+EijR(24,3)-EijR(39,4)-EijR(40,4)+EijR(4
     -   2,4)+EijR(44,4))-EijR(2,1)*P(9)+EijR(5,2)*P(15)+(EijR(4,3)-E
     -   ijR(14,3))*P(119)+EijR(9,2)*P(126)-EijR(8,2)*P(128)+EijR(3,1
     -   )*P(160)+EijR(4,1)*P(161)-EijR(4,2)*P(162)+EijR(6,2)*P(163)+
     -   EijR(7,2)*P(164)-EijR(10,2)*P(165)+(-EijR(15,3)+EijR(18,3))*
     -   P(166))
       FI(194) = 8*(Dij2345I(1,1)-Dij2345I(3,1)-s15*(EE0I+2*EijI(1,1)+
     -   EijI(1,2))+s12*EijI(2,2)+p2sq*(-EijI(8,3)+EijI(10,3))+2*(2*E
     -   ijI(11,2)-EijI(22,3)+EijI(24,3)-EijI(39,4)-EijI(40,4)+EijI(4
     -   2,4)+EijI(44,4))-EijI(2,1)*P(9)+EijI(5,2)*P(15)+(EijI(4,3)-E
     -   ijI(14,3))*P(119)+EijI(9,2)*P(126)-EijI(8,2)*P(128)+EijI(3,1
     -   )*P(160)+EijI(4,1)*P(161)-EijI(4,2)*P(162)+EijI(6,2)*P(163)+
     -   EijI(7,2)*P(164)-EijI(10,2)*P(165)+(-EijI(15,3)+EijI(18,3))*
     -   P(166))
       F(194)=DCMPLX(FR(194),FI(194))
       P(167) = p2sq-s12+s15
       P(168) = p2sq-2*s12
       P(169) = p2sq-p4sq-s23+s34-2*P(7)
       P(170) = p4sq+s12+s15+s23-s34-2*s45
       P(171) = p2sq+p3sq+p4sq-s15+s23-s34-3*s45
       P(172) = -p2sq+p3sq+p4sq+s23-s34+3*P(7)
       P(173) = p3sq+p4sq-s15+s23-s34-3*s45+2*P(9)
       P(174) = p2sq+p3sq-3*s12-s15-s45
       P(175) = p2sq+p3sq-s15-s45
       P(176) = -3*s12+2*P(175)
       FR(195) = -8*(EE0R*s15-Dij2345R(1,1)-Dij2345R(1,2)+Dij2345R(3,1
     -   )-Dij2345R(3,2)+s15*EijR(1,1)-s12*EijR(2,3)+2*(Dij2345R(5,2)
     -   -EijR(11,2)-EijR(22,3)+EijR(24,3)+EijR(37,4)+EijR(39,4))-4*E
     -   ijR(44,4)+(EijR(8,3)+EijR(14,3)-2*EijR(18,3))*P(9)+(EijR(9,3
     -   )+EijR(16,3)-2*EijR(20,3))*P(128)+EijR(4,3)*P(133)+(EijR(2,1
     -   )+EijR(5,2)-EijR(7,2))*P(167)+EijR(2,2)*P(168)-EijR(3,1)*P(1
     -   69)-EijR(4,1)*P(170)+EijR(4,2)*P(171)+(EijR(8,2)-EijR(10,2))
     -   *P(172)-EijR(9,2)*P(173)-EijR(10,3)*P(174)+EijR(15,3)*P(176)
     -   )
       FI(195) = -8*(EE0I*s15-Dij2345I(1,1)-Dij2345I(1,2)+Dij2345I(3,1
     -   )-Dij2345I(3,2)+s15*EijI(1,1)-s12*EijI(2,3)+2*(Dij2345I(5,2)
     -   -EijI(11,2)-EijI(22,3)+EijI(24,3)+EijI(37,4)+EijI(39,4))-4*E
     -   ijI(44,4)+(EijI(8,3)+EijI(14,3)-2*EijI(18,3))*P(9)+(EijI(9,3
     -   )+EijI(16,3)-2*EijI(20,3))*P(128)+EijI(4,3)*P(133)+(EijI(2,1
     -   )+EijI(5,2)-EijI(7,2))*P(167)+EijI(2,2)*P(168)-EijI(3,1)*P(1
     -   69)-EijI(4,1)*P(170)+EijI(4,2)*P(171)+(EijI(8,2)-EijI(10,2))
     -   *P(172)-EijI(9,2)*P(173)-EijI(10,3)*P(174)+EijI(15,3)*P(176)
     -   )
       F(195)=DCMPLX(FR(195),FI(195))
       P(177) = p2sq-2*s12-s15-s23+s34+s45
       P(178) = s12+s15+s23-s34-s45
       P(179) = p2sq+p3sq-s15+s23-s34-4*s45
       P(180) = p2sq-s12-2*P(13)
       P(181) = p2sq-s12+2*P(13)
       P(182) = p2sq+2*p3sq+s12+s15-s23-s34-s45
       P(183) = p2sq+2*p3sq-3*s12-s23-s34-s45
       P(184) = 4*s12+s23-s34-5*s45
       FR(196) = -8*(D02345R-Dij2345R(1,1)-Dij2345R(2,1)-Dij2345R(3,2)
     -   -Dij2345R(4,2)+Dij2345R(5,2)+Dij2345R(6,2)+s15*(EE0R+EijR(1,
     -   1))-s12*EijR(9,3)+s12*EijR(10,3)-10*EijR(11,2)+2*(Dij2345R(3
     -   ,1)+EijR(39,4)+EijR(43,4)-EijR(44,4)-EijR(45,4))+(EijR(2,1)+
     -   2*EijR(5,2)+EijR(14,3)+EijR(17,3)-EijR(18,3)-EijR(19,3))*P(9
     -   )-EijR(3,2)*P(38)+EijR(15,3)*P(127)+(EijR(12,3)-EijR(13,3))*
     -   P(128)+EijR(4,3)*P(133)+EijR(16,3)*P(139)-EijR(20,3)*P(140)-
     -   EijR(2,2)*P(152)-EijR(3,1)*P(177)-EijR(4,1)*P(178)+EijR(4,2)
     -   *P(179)-EijR(6,2)*P(180)-EijR(7,2)*P(181)+EijR(8,2)*P(182)-E
     -   ijR(9,2)*P(183)-EijR(10,2)*P(184))
       FI(196) = -8*(D02345I-Dij2345I(1,1)-Dij2345I(2,1)-Dij2345I(3,2)
     -   -Dij2345I(4,2)+Dij2345I(5,2)+Dij2345I(6,2)+s15*(EE0I+EijI(1,
     -   1))-s12*EijI(9,3)+s12*EijI(10,3)-10*EijI(11,2)+2*(Dij2345I(3
     -   ,1)+EijI(39,4)+EijI(43,4)-EijI(44,4)-EijI(45,4))+(EijI(2,1)+
     -   2*EijI(5,2)+EijI(14,3)+EijI(17,3)-EijI(18,3)-EijI(19,3))*P(9
     -   )-EijI(3,2)*P(38)+EijI(15,3)*P(127)+(EijI(12,3)-EijI(13,3))*
     -   P(128)+EijI(4,3)*P(133)+EijI(16,3)*P(139)-EijI(20,3)*P(140)-
     -   EijI(2,2)*P(152)-EijI(3,1)*P(177)-EijI(4,1)*P(178)+EijI(4,2)
     -   *P(179)-EijI(6,2)*P(180)-EijI(7,2)*P(181)+EijI(8,2)*P(182)-E
     -   ijI(9,2)*P(183)-EijI(10,2)*P(184))
       F(196)=DCMPLX(FR(196),FI(196))
       FR(197) = 8*(Dij2345R(3,1)+Dij2345R(3,2)+p2sq*EijR(9,2)+s12*Eij
     -   R(15,3)-s12*EijR(18,3)-2*(EijR(39,4)-EijR(42,4))+EijR(7,3)*P
     -   (9)-EijR(10,2)*P(16)-EijR(4,2)*P(34)+EijR(7,2)*P(83)+EijR(4,
     -   3)*P(121)-EijR(14,3)*P(126)+(-EijR(16,3)+EijR(19,3))*P(128))
       FI(197) = 8*(Dij2345I(3,1)+Dij2345I(3,2)+p2sq*EijI(9,2)+s12*Eij
     -   I(15,3)-s12*EijI(18,3)-2*(EijI(39,4)-EijI(42,4))+EijI(7,3)*P
     -   (9)-EijI(10,2)*P(16)-EijI(4,2)*P(34)+EijI(7,2)*P(83)+EijI(4,
     -   3)*P(121)-EijI(14,3)*P(126)+(-EijI(16,3)+EijI(19,3))*P(128))
       F(197)=DCMPLX(FR(197),FI(197))
       P(185) = p2sq+s15-s23+s45
       P(186) = -2*p3sq+s15+s23+s45
       FR(198) = -8*(s12*(EijR(2,1)+EijR(5,2))+p2sq*(EijR(15,3)-EijR(1
     -   8,3))-2*(EijR(21,3)-EijR(39,4)+EijR(42,4))+(EijR(3,1)+EijR(6
     -   ,2))*P(14)-EijR(10,2)*P(16)-EijR(4,2)*P(18)+(-EijR(4,3)+EijR
     -   (14,3))*P(119)-EijR(4,1)*P(185)-EijR(7,2)*P(186))
       FI(198) = -8*(s12*(EijI(2,1)+EijI(5,2))+p2sq*(EijI(15,3)-EijI(1
     -   8,3))-2*(EijI(21,3)-EijI(39,4)+EijI(42,4))+(EijI(3,1)+EijI(6
     -   ,2))*P(14)-EijI(10,2)*P(16)-EijI(4,2)*P(18)+(-EijI(4,3)+EijI
     -   (14,3))*P(119)-EijI(4,1)*P(185)-EijI(7,2)*P(186))
       F(198)=DCMPLX(FR(198),FI(198))
       P(187) = -p3sq+s12+s15
       P(188) = 2*p2sq+p3sq+p4sq-s12-s23+s45
       P(189) = 2*p2sq-p3sq+s34+s45
       P(190) = p2sq-p3sq-p4sq-s12-s23+s34+s45
       FR(199) = 8*(-Dij2345R(1,1)+Dij2345R(3,1)+Dij2345R(3,2)-Dij2345
     -   R(5,2)+p3sq*EijR(3,1)-s12*(EijR(2,1)+EijR(10,3))-2*(2*EijR(1
     -   1,2)+EijR(24,3)+EijR(39,4)-EijR(44,4))+(EijR(5,2)-EijR(14,3)
     -   +EijR(18,3))*P(9)-EijR(15,3)*P(127)+(EijR(8,2)-EijR(16,3)+Ei
     -   jR(20,3))*P(128)-EijR(4,3)*P(133)+EijR(7,2)*P(135)-EijR(2,2)
     -   *P(152)+EijR(4,1)*P(187)-EijR(4,2)*P(188)+EijR(9,2)*P(189)-E
     -   ijR(10,2)*P(190))
       FI(199) = 8*(-Dij2345I(1,1)+Dij2345I(3,1)+Dij2345I(3,2)-Dij2345
     -   I(5,2)+p3sq*EijI(3,1)-s12*(EijI(2,1)+EijI(10,3))-2*(2*EijI(1
     -   1,2)+EijI(24,3)+EijI(39,4)-EijI(44,4))+(EijI(5,2)-EijI(14,3)
     -   +EijI(18,3))*P(9)-EijI(15,3)*P(127)+(EijI(8,2)-EijI(16,3)+Ei
     -   jI(20,3))*P(128)-EijI(4,3)*P(133)+EijI(7,2)*P(135)-EijI(2,2)
     -   *P(152)+EijI(4,1)*P(187)-EijI(4,2)*P(188)+EijI(9,2)*P(189)-E
     -   ijI(10,2)*P(190))
       F(199)=DCMPLX(FR(199),FI(199))
       P(191) = p2sq+p3sq-s12+s15-s23+s45
       P(192) = s12+s23-s34-3*s45
       P(193) = -p3sq+s12-s15+s23+s34
       P(194) = p2sq+s12-s23-s45
       FR(200) = -8*(D02345R+Dij2345R(3,1)-Dij2345R(3,2)+Dij2345R(5,2)
     -   -p2sq*EijR(2,2)+s12*EijR(10,3)+EijR(3,1)*P(5)+(EijR(5,2)+Eij
     -   R(14,3)-EijR(18,3))*P(9)-2*(2*EijR(11,2)+EijR(22,3)-EijR(24,
     -   3)-EijR(39,4)+EijR(44,4)+EijR(7,2)*P(13))+EijR(8,2)*P(16)-Ei
     -   jR(4,3)*P(121)+EijR(15,3)*P(127)+(EijR(16,3)-EijR(20,3))*P(1
     -   28)-EijR(4,1)*P(191)+EijR(4,2)*P(192)+EijR(9,2)*P(193)-EijR(
     -   10,2)*P(194))
       FI(200) = -8*(D02345I+Dij2345I(3,1)-Dij2345I(3,2)+Dij2345I(5,2)
     -   -p2sq*EijI(2,2)+s12*EijI(10,3)+EijI(3,1)*P(5)+(EijI(5,2)+Eij
     -   I(14,3)-EijI(18,3))*P(9)-2*(2*EijI(11,2)+EijI(22,3)-EijI(24,
     -   3)-EijI(39,4)+EijI(44,4)+EijI(7,2)*P(13))+EijI(8,2)*P(16)-Ei
     -   jI(4,3)*P(121)+EijI(15,3)*P(127)+(EijI(16,3)-EijI(20,3))*P(1
     -   28)-EijI(4,1)*P(191)+EijI(4,2)*P(192)+EijI(9,2)*P(193)-EijI(
     -   10,2)*P(194))
       F(200)=DCMPLX(FR(200),FI(200))
       P(195) = -p3sq+s12+s15+s45
       P(196) = p2sq+s12+s15-s23+s45
       FR(201) = -8*(D02345R+Dij2345R(2,1)-Dij2345R(3,2)+Dij2345R(6,2)
     -   +s12*(EijR(2,1)-EijR(15,3)+EijR(20,3))-2*(EijR(11,2)-EijR(39
     -   ,4)+EijR(45,4))-(EijR(6,2)+EijR(7,2)-EijR(14,3)+EijR(19,3))*
     -   P(9)+EijR(4,2)*P(34)-EijR(3,1)*P(114)-EijR(4,3)*P(121)-(EijR
     -   (3,2)+EijR(13,3))*P(128)+EijR(16,3)*P(139)+EijR(8,2)*P(152)-
     -   EijR(4,1)*P(195)-EijR(10,2)*P(196))
       FI(201) = -8*(D02345I+Dij2345I(2,1)-Dij2345I(3,2)+Dij2345I(6,2)
     -   +s12*(EijI(2,1)-EijI(15,3)+EijI(20,3))-2*(EijI(11,2)-EijI(39
     -   ,4)+EijI(45,4))-(EijI(6,2)+EijI(7,2)-EijI(14,3)+EijI(19,3))*
     -   P(9)+EijI(4,2)*P(34)-EijI(3,1)*P(114)-EijI(4,3)*P(121)-(EijI
     -   (3,2)+EijI(13,3))*P(128)+EijI(16,3)*P(139)+EijI(8,2)*P(152)-
     -   EijI(4,1)*P(195)-EijI(10,2)*P(196))
       F(201)=DCMPLX(FR(201),FI(201))
       P(197) = p3sq-p4sq+s34+2*P(14)
       P(198) = p3sq-s15-s23+s34+2*P(9)
       P(199) = p3sq-s23+2*P(9)
       P(200) = -p3sq-p4sq+s12+s23+s45
       P(201) = p2sq*P(195)-s12*P(200)
       P(202) = s12-s15+s23
       P(203) = p4sq+s12-s45
       P(204) = p4sq+s15-s34
       P(205) = -s45**2+p2sq*P(3)+s12*P(202)-p3sq*P(203)+s45*P(204)
       P(206) = -p4sq-3*s23+s34-2*P(7)
       P(207) = p4sq+s23-s34+2*P(7)
       P(208) = 3*p3sq+s12-s45
       P(209) = 2*p3sq**2+p3sq*P(206)-P(7)*P(207)+p2sq*P(208)
       P(210) = s12+s15+s45
       P(211) = -p2sq+s12+s34+s45
       P(212) = p3sq*P(11)-p4sq*P(210)+s45*P(211)
       P(213) = p4sq+s45
       P(214) = p3sq-s15
       P(215) = 5*s12+s23-s34+s45-2*P(214)
       P(216) = 2*p3sq+s34+s45-3*P(48)
       P(217) = -2*p2sq**2+p2sq*P(215)+s12*P(216)
       P(218) = 2*s12+s15
       P(219) = 2*p3sq-p4sq-5*s12-3*s23+s34+4*s45
       P(220) = p2sq**2-s23*s45
       P(221) = p4sq*s12-s12*s15+s23**2-s12*s34-4*s12*s45+s15*s45+s45*
     -   *2+3*s12*P(48)-p3sq*P(218)+p2sq*P(219)+2*P(220)
       P(222) = p4sq-2*s23+3*s45
       P(223) = s12*s23-s23*s34-s23*s45+s34*s45
       P(224) = -(p4sq*s12)-p3sq*s15+s15*s23+s23**2+s15*s34-3*s12*s45+
     -   s45**2+p2sq*P(222)+2*P(223)
       P(225) = 4*p3sq-p4sq+s12-3*s23+s34
       P(226) = p2sq**2+p3sq**2-s12*s23+s12*s34
       P(227) = 4*s12+s15+s23-s34
       P(228) = s12+s15+3*s23-s34+s45
       P(229) = p4sq*s12+4*s12**2+s12*s15-s23**2-p2sq*P(225)-2*P(226)-
     -   s45*P(227)+p3sq*P(228)
       P(230) = 2*p3sq-s12-s23-s34-s45
       P(231) = -5*s12-s15+s34-3*P(13)
       P(232) = s15+2*s34
       P(233) = -s23+s34+s45
       P(234) = 4*p3sq+p4sq+s45-2*P(71)
       P(235) = p3sq**2+s12**2+s23*s34-s45**2
       P(236) = s12*s15+s12*s23+s12*s34-s34**2+p3sq*P(231)-s45*P(232)+
     -   p4sq*P(233)+p2sq*P(234)+2*P(235)
       P(237) = 8*p2sq+5*p3sq-p4sq-s15-7*P(48)+3*P(58)
       FR(202) = 4*(-4*EijR(46,4)+EE0R*s15*P(1)+s15*EijR(1,2)*P(47)-Di
     -   j2345R(2,1)*P(169)-D02345R*P(197)-s15*EijR(1,1)*P(199)+EijR(
     -   2,1)*P(201)+EijR(3,1)*P(205)+EijR(3,2)*P(209)+EijR(4,1)*P(21
     -   2)+EijR(5,2)*P(217)+EijR(6,2)*P(221)-EijR(7,2)*P(224)+EijR(8
     -   ,2)*P(229)+P(198)*(Dij2345R(1,1)-Dij2345R(3,1)+EijR(2,2)*P(1
     -   52)+EijR(4,2)*P(213)+EijR(9,2)*P(230))-EijR(10,2)*P(236)-2*(
     -   p2sq*EijR(22,3)-EijR(24,3)*P(119)-EijR(11,2)*P(237)))
       FI(202) = 4*(-4*EijI(46,4)+EE0I*s15*P(1)+s15*EijI(1,2)*P(47)-Di
     -   j2345I(2,1)*P(169)-D02345I*P(197)-s15*EijI(1,1)*P(199)+EijI(
     -   2,1)*P(201)+EijI(3,1)*P(205)+EijI(3,2)*P(209)+EijI(4,1)*P(21
     -   2)+EijI(5,2)*P(217)+EijI(6,2)*P(221)-EijI(7,2)*P(224)+EijI(8
     -   ,2)*P(229)+P(198)*(Dij2345I(1,1)-Dij2345I(3,1)+EijI(2,2)*P(1
     -   52)+EijI(4,2)*P(213)+EijI(9,2)*P(230))-EijI(10,2)*P(236)-2*(
     -   p2sq*EijI(22,3)-EijI(24,3)*P(119)-EijI(11,2)*P(237)))
       F(202)=DCMPLX(FR(202),FI(202))
       FR(203) = 8*(Dij2345R(3,1)+Dij2345R(3,2)+s12*(EijR(9,2)+EijR(15
     -   ,3))-2*(EijR(24,3)+EijR(39,4))-EijR(14,3)*P(9)+EijR(10,2)*P(
     -   14)+EijR(4,3)*P(121)-EijR(16,3)*P(128)-EijR(4,2)*P(185))
       FI(203) = 8*(Dij2345I(3,1)+Dij2345I(3,2)+s12*(EijI(9,2)+EijI(15
     -   ,3))-2*(EijI(24,3)+EijI(39,4))-EijI(14,3)*P(9)+EijI(10,2)*P(
     -   14)+EijI(4,3)*P(121)-EijI(16,3)*P(128)-EijI(4,2)*P(185))
       F(203)=DCMPLX(FR(203),FI(203))
       P(238) = -p2sq+p3sq+p4sq+s12+s23-s34-s45
       P(239) = p2sq-p4sq-s12-s23+s34+s45
       FR(204) = -8*(p3sq*(EijR(3,1)-EijR(4,1)+EijR(6,3)+EijR(16,3)+2*
     -   (EijR(6,2)-EijR(10,2)-EijR(19,3)))+2*(EijR(11,2)+EijR(36,4)+
     -   EijR(39,4)+2*(EijR(21,3)-EijR(24,3)-EijR(42,4)))-(EijR(5,2)+
     -   EijR(5,3)-EijR(9,2)+EijR(15,3)-2*EijR(18,3))*P(16)+(EijR(4,2
     -   )-EijR(7,2))*P(238)+(EijR(4,3)+EijR(7,3)-2*EijR(14,3))*P(239
     -   ))
       FI(204) = -8*(p3sq*(EijI(3,1)-EijI(4,1)+EijI(6,3)+EijI(16,3)+2*
     -   (EijI(6,2)-EijI(10,2)-EijI(19,3)))+2*(EijI(11,2)+EijI(36,4)+
     -   EijI(39,4)+2*(EijI(21,3)-EijI(24,3)-EijI(42,4)))-(EijI(5,2)+
     -   EijI(5,3)-EijI(9,2)+EijI(15,3)-2*EijI(18,3))*P(16)+(EijI(4,2
     -   )-EijI(7,2))*P(238)+(EijI(4,3)+EijI(7,3)-2*EijI(14,3))*P(239
     -   ))
       F(204)=DCMPLX(FR(204),FI(204))
       P(240) = p2sq+p3sq+p4sq-s23-s34-s45
       P(241) = p2sq+p3sq-p4sq-s12-s23
       P(242) = p2sq+2*p3sq+p4sq-s23-s34-s45
       P(243) = -p3sq+s12+s23
       P(244) = p3sq-p4sq-s23+2*P(9)
       P(245) = p2sq-p4sq+2*P(26)
       P(246) = -p4sq-3*s12+2*P(16)
       P(247) = -s12+s45
       FR(205) = 8*(Dij2345R(1,1)-Dij2345R(3,1)-Dij2345R(3,2)+Dij2345R
     -   (5,2)+p3sq*EijR(6,2)-p3sq*EijR(8,2)-4*(EijR(21,3)+EijR(22,3)
     -   )+8*EijR(24,3)-2*(EijR(39,4)+EijR(40,4)-EijR(42,4)-EijR(44,4
     -   ))+EijR(5,2)*P(2)+(EijR(16,3)+EijR(17,3))*P(7)+EijR(5,3)*P(9
     -   )-EijR(9,2)*P(43)+EijR(7,3)*P(47)+EijR(4,2)*P(240)-EijR(4,3)
     -   *P(241)-EijR(7,2)*P(242)+(-EijR(8,3)+EijR(10,3))*P(243)+EijR
     -   (14,3)*P(244)+EijR(15,3)*P(245)-EijR(18,3)*P(246)+(EijR(19,3
     -   )+EijR(20,3))*P(247))
       FI(205) = 8*(Dij2345I(1,1)-Dij2345I(3,1)-Dij2345I(3,2)+Dij2345I
     -   (5,2)+p3sq*EijI(6,2)-p3sq*EijI(8,2)-4*(EijI(21,3)+EijI(22,3)
     -   )+8*EijI(24,3)-2*(EijI(39,4)+EijI(40,4)-EijI(42,4)-EijI(44,4
     -   ))+EijI(5,2)*P(2)+(EijI(16,3)+EijI(17,3))*P(7)+EijI(5,3)*P(9
     -   )-EijI(9,2)*P(43)+EijI(7,3)*P(47)+EijI(4,2)*P(240)-EijI(4,3)
     -   *P(241)-EijI(7,2)*P(242)+(-EijI(8,3)+EijI(10,3))*P(243)+EijI
     -   (14,3)*P(244)+EijI(15,3)*P(245)-EijI(18,3)*P(246)+(EijI(19,3
     -   )+EijI(20,3))*P(247))
       F(205)=DCMPLX(FR(205),FI(205))
       P(248) = p3sq-p4sq-s12+s34+s45
       P(249) = -p2sq+p3sq+p4sq+s12+s23-s34-2*s45
       P(250) = -p2sq+s12+s23-s34-2*s45
       P(251) = p2sq+p4sq+s12-s23-s34
       P(252) = -p2sq+p4sq+s12+s23-s34+2*P(114)
       P(253) = p2sq+p3sq+p4sq+s12-s23-s34
       P(254) = p2sq-s12-2*s23+s34+3*s45
       P(255) = p2sq+s45
       P(256) = -s23+s34+2*P(255)
       P(257) = s34+3*s45+2*P(35)
       FR(206) = -8*(Dij2345R(1,1)-Dij2345R(3,1)-Dij2345R(3,2)+Dij2345
     -   R(5,2)+p3sq*(EijR(3,1)+2*EijR(6,2)-EijR(8,2)-EijR(10,2))+2*(
     -   EijR(11,2)-EijR(21,3)+EijR(24,3)+EijR(39,4)+EijR(40,4)-EijR(
     -   42,4)-EijR(44,4))+(EijR(16,3)+EijR(17,3)-EijR(19,3)-EijR(20,
     -   3))*P(7)+(EijR(5,3)-EijR(7,3))*P(13)-EijR(2,1)*P(86)+(-EijR(
     -   8,3)+EijR(10,3))*P(152)-EijR(4,1)*P(248)+EijR(4,2)*P(249)-Ei
     -   jR(4,3)*P(250)-EijR(5,2)*P(251)-EijR(7,2)*P(252)+EijR(9,2)*P
     -   (253)-EijR(14,3)*P(254)-EijR(15,3)*P(256)+EijR(18,3)*P(257))
       FI(206) = -8*(Dij2345I(1,1)-Dij2345I(3,1)-Dij2345I(3,2)+Dij2345
     -   I(5,2)+p3sq*(EijI(3,1)+2*EijI(6,2)-EijI(8,2)-EijI(10,2))+2*(
     -   EijI(11,2)-EijI(21,3)+EijI(24,3)+EijI(39,4)+EijI(40,4)-EijI(
     -   42,4)-EijI(44,4))+(EijI(16,3)+EijI(17,3)-EijI(19,3)-EijI(20,
     -   3))*P(7)+(EijI(5,3)-EijI(7,3))*P(13)-EijI(2,1)*P(86)+(-EijI(
     -   8,3)+EijI(10,3))*P(152)-EijI(4,1)*P(248)+EijI(4,2)*P(249)-Ei
     -   jI(4,3)*P(250)-EijI(5,2)*P(251)-EijI(7,2)*P(252)+EijI(9,2)*P
     -   (253)-EijI(14,3)*P(254)-EijI(15,3)*P(256)+EijI(18,3)*P(257))
       F(206)=DCMPLX(FR(206),FI(206))
       P(258) = p2sq+p3sq+s12-s23-s45+2*P(79)
       FR(207) = 8*(p3sq*(EijR(8,2)+EijR(9,3)-EijR(10,2)+EijR(16,3)-2*
     -   EijR(20,3))-2*(EijR(22,3)-EijR(24,3)+EijR(37,4)+EijR(39,4)-2
     -   *EijR(44,4))-(EijR(4,3)+EijR(10,3)-2*EijR(15,3))*P(5)+(EijR(
     -   5,2)-EijR(7,2)+EijR(8,3)+EijR(14,3)-2*EijR(18,3))*P(14)+EijR
     -   (4,2)*P(64)+(EijR(2,1)+EijR(2,2)-EijR(4,1))*P(86)-EijR(9,2)*
     -   P(258))
       FI(207) = 8*(p3sq*(EijI(8,2)+EijI(9,3)-EijI(10,2)+EijI(16,3)-2*
     -   EijI(20,3))-2*(EijI(22,3)-EijI(24,3)+EijI(37,4)+EijI(39,4)-2
     -   *EijI(44,4))-(EijI(4,3)+EijI(10,3)-2*EijI(15,3))*P(5)+(EijI(
     -   5,2)-EijI(7,2)+EijI(8,3)+EijI(14,3)-2*EijI(18,3))*P(14)+EijI
     -   (4,2)*P(64)+(EijI(2,1)+EijI(2,2)-EijI(4,1))*P(86)-EijI(9,2)*
     -   P(258))
       F(207)=DCMPLX(FR(207),FI(207))
       P(259) = -p4sq-s12+s34+s45
       P(260) = p2sq+p3sq+p4sq-s23-s34+s45
       P(261) = p2sq+p3sq-p4sq-s12-s23+s34+s45
       P(262) = p3sq-p4sq-2*s12-s23+s34+s45
       P(263) = p2sq+3*p3sq+p4sq-s12-s23-s34+s45
       P(264) = p2sq+2*p3sq-s12-s23+s45
       FR(208) = -8*(D02345R-p2sq*EijR(2,2)+p3sq*(-2*EijR(3,2)-EijR(12
     -   ,3)+EijR(13,3))-2*(4*EijR(11,2)+EijR(22,3)-EijR(24,3)-EijR(3
     -   9,4)-EijR(43,4)+EijR(44,4)+EijR(45,4))+(EijR(4,3)-EijR(15,3)
     -   )*P(5)+EijR(5,2)*P(9)-(2*EijR(6,2)+EijR(14,3)+EijR(17,3))*P(
     -   14)-(EijR(18,3)+EijR(19,3))*P(60)-EijR(3,1)*P(86)+EijR(7,2)*
     -   P(180)-EijR(4,1)*P(259)-EijR(4,2)*P(260)+EijR(8,2)*P(261)-Ei
     -   jR(9,2)*P(262)+EijR(10,2)*P(263)+(-EijR(16,3)+EijR(20,3))*P(
     -   264))
       FI(208) = -8*(D02345I-p2sq*EijI(2,2)+p3sq*(-2*EijI(3,2)-EijI(12
     -   ,3)+EijI(13,3))-2*(4*EijI(11,2)+EijI(22,3)-EijI(24,3)-EijI(3
     -   9,4)-EijI(43,4)+EijI(44,4)+EijI(45,4))+(EijI(4,3)-EijI(15,3)
     -   )*P(5)+EijI(5,2)*P(9)-(2*EijI(6,2)+EijI(14,3)+EijI(17,3))*P(
     -   14)-(EijI(18,3)+EijI(19,3))*P(60)-EijI(3,1)*P(86)+EijI(7,2)*
     -   P(180)-EijI(4,1)*P(259)-EijI(4,2)*P(260)+EijI(8,2)*P(261)-Ei
     -   jI(9,2)*P(262)+EijI(10,2)*P(263)+(-EijI(16,3)+EijI(20,3))*P(
     -   264))
       F(208)=DCMPLX(FR(208),FI(208))
       P(265) = p3sq+p4sq-s12+2*P(30)
       P(266) = p2sq-p3sq+p4sq+2*s12
       FR(209) = 8*(Dij2345R(1,1)-Dij2345R(3,1)+s12*EijR(2,1)-p3sq*(Ei
     -   jR(3,1)-EijR(16,3)+EijR(20,3))+2*(2*EijR(11,2)-EijR(39,4)+Ei
     -   jR(44,4))+EijR(4,1)*P(3)+(-EijR(4,3)+EijR(15,3))*P(5)-EijR(5
     -   ,2)*P(9)+EijR(14,3)*P(14)-EijR(10,2)*P(38)+EijR(18,3)*P(60)-
     -   EijR(8,2)*P(128)+EijR(2,2)*P(152)-EijR(7,2)*P(180)+EijR(4,2)
     -   *P(265)-EijR(9,2)*P(266))
       FI(209) = 8*(Dij2345I(1,1)-Dij2345I(3,1)+s12*EijI(2,1)-p3sq*(Ei
     -   jI(3,1)-EijI(16,3)+EijI(20,3))+2*(2*EijI(11,2)-EijI(39,4)+Ei
     -   jI(44,4))+EijI(4,1)*P(3)+(-EijI(4,3)+EijI(15,3))*P(5)-EijI(5
     -   ,2)*P(9)+EijI(14,3)*P(14)-EijI(10,2)*P(38)+EijI(18,3)*P(60)-
     -   EijI(8,2)*P(128)+EijI(2,2)*P(152)-EijI(7,2)*P(180)+EijI(4,2)
     -   *P(265)-EijI(9,2)*P(266))
       F(209)=DCMPLX(FR(209),FI(209))
       P(267) = -s12+s34
       P(268) = p2sq+p3sq-s23-s34-s45
       P(269) = -s15+2*P(9)
       P(270) = p2sq-s12-s15+s23-s45
       P(271) = p3sq-s23-s34
       P(272) = s15+2*P(271)
       FR(210) = 8*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)-p3sq*(EijR(3,2
     -   )-EijR(12,3)+EijR(13,3))+(-EijR(4,3)+EijR(15,3))*P(5)+(-EijR
     -   (6,2)+EijR(14,3)+EijR(17,3))*P(14)-2*(4*EijR(11,2)-EijR(22,3
     -   )+EijR(23,3)+EijR(39,4)+EijR(43,4)-EijR(44,4)-EijR(45,4)-Eij
     -   R(8,2)*P(16))+EijR(2,1)*P(55)+(EijR(18,3)+EijR(19,3))*P(60)-
     -   EijR(2,2)*P(166)-EijR(10,2)*P(194)+(EijR(16,3)-EijR(20,3))*P
     -   (264)+EijR(4,1)*P(267)+EijR(4,2)*P(268)+EijR(5,2)*P(269)-Eij
     -   R(7,2)*P(270)-EijR(9,2)*P(272))
       FI(210) = 8*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)-p3sq*(EijI(3,2
     -   )-EijI(12,3)+EijI(13,3))+(-EijI(4,3)+EijI(15,3))*P(5)+(-EijI
     -   (6,2)+EijI(14,3)+EijI(17,3))*P(14)-2*(4*EijI(11,2)-EijI(22,3
     -   )+EijI(23,3)+EijI(39,4)+EijI(43,4)-EijI(44,4)-EijI(45,4)-Eij
     -   I(8,2)*P(16))+EijI(2,1)*P(55)+(EijI(18,3)+EijI(19,3))*P(60)-
     -   EijI(2,2)*P(166)-EijI(10,2)*P(194)+(EijI(16,3)-EijI(20,3))*P
     -   (264)+EijI(4,1)*P(267)+EijI(4,2)*P(268)+EijI(5,2)*P(269)-Eij
     -   I(7,2)*P(270)-EijI(9,2)*P(272))
       F(210)=DCMPLX(FR(210),FI(210))
       P(273) = p2sq+p3sq-s23-s34
       P(274) = -p2sq+s12+s15
       P(275) = -s34+2*P(16)
       P(276) = p2sq+3*p3sq-s12-s23+s45
       P(277) = 3*p3sq+2*P(14)
       FR(211) = 8*(-Dij2345R(2,1)+Dij2345R(3,1)+p3sq*EijR(3,3)+4*EijR
     -   (23,3)-EijR(4,3)*P(5)+(EijR(11,3)+EijR(14,3))*P(14)+EijR(3,2
     -   )*P(16)+EijR(3,1)*P(55)-2*(EijR(38,4)+EijR(39,4)+2*(EijR(24,
     -   3)-EijR(45,4))-EijR(19,3)*P(60))+(-EijR(8,2)+EijR(9,2))*P(11
     -   9)+EijR(4,1)*P(267)+EijR(4,2)*P(273)+(-EijR(6,2)+EijR(7,2))*
     -   P(274)-EijR(10,2)*P(275)-EijR(13,3)*P(276)+EijR(16,3)*P(277)
     -   )
       FI(211) = 8*(-Dij2345I(2,1)+Dij2345I(3,1)+p3sq*EijI(3,3)+4*EijI
     -   (23,3)-EijI(4,3)*P(5)+(EijI(11,3)+EijI(14,3))*P(14)+EijI(3,2
     -   )*P(16)+EijI(3,1)*P(55)-2*(EijI(38,4)+EijI(39,4)+2*(EijI(24,
     -   3)-EijI(45,4))-EijI(19,3)*P(60))+(-EijI(8,2)+EijI(9,2))*P(11
     -   9)+EijI(4,1)*P(267)+EijI(4,2)*P(273)+(-EijI(6,2)+EijI(7,2))*
     -   P(274)-EijI(10,2)*P(275)-EijI(13,3)*P(276)+EijI(16,3)*P(277)
     -   )
       F(211)=DCMPLX(FR(211),FI(211))
       P(278) = s12-s34-2*s45
       P(279) = p2sq-s12-s34-s45+2*P(18)
       FR(212) = 8*(-Dij2345R(1,1)+Dij2345R(3,1)+p3sq*(EijR(16,3)-EijR
     -   (20,3))-2*(EijR(11,2)-EijR(22,3)+EijR(39,4)-EijR(44,4))+EijR
     -   (4,1)*P(1)+EijR(2,1)*P(2)+(-EijR(4,3)+EijR(15,3))*P(5)+EijR(
     -   8,2)*P(7)+EijR(5,2)*P(9)-EijR(7,2)*P(13)+EijR(14,3)*P(14)+Ei
     -   jR(10,2)*P(38)+EijR(18,3)*P(60)-EijR(2,2)*P(243)+EijR(4,2)*P
     -   (278)-EijR(9,2)*P(279))
       FI(212) = 8*(-Dij2345I(1,1)+Dij2345I(3,1)+p3sq*(EijI(16,3)-EijI
     -   (20,3))-2*(EijI(11,2)-EijI(22,3)+EijI(39,4)-EijI(44,4))+EijI
     -   (4,1)*P(1)+EijI(2,1)*P(2)+(-EijI(4,3)+EijI(15,3))*P(5)+EijI(
     -   8,2)*P(7)+EijI(5,2)*P(9)-EijI(7,2)*P(13)+EijI(14,3)*P(14)+Ei
     -   jI(10,2)*P(38)+EijI(18,3)*P(60)-EijI(2,2)*P(243)+EijI(4,2)*P
     -   (278)-EijI(9,2)*P(279))
       F(212)=DCMPLX(FR(212),FI(212))
       P(280) = s12-s34-s45
       FR(213) = -8*(Dij2345R(2,1)-Dij2345R(3,1)+p3sq*EijR(13,3)+4*Eij
     -   R(24,3)-2*(EijR(23,3)-EijR(39,4)+EijR(45,4))-EijR(4,1)*P(1)+
     -   (-EijR(3,1)+EijR(9,2))*P(2)+EijR(4,3)*P(5)-EijR(3,2)*P(7)-Ei
     -   jR(6,2)*P(9)-EijR(14,3)*P(14)-EijR(19,3)*P(60)+EijR(8,2)*P(2
     -   43)-EijR(16,3)*P(264)+EijR(10,2)*P(278)-EijR(4,2)*P(280))
       FI(213) = -8*(Dij2345I(2,1)-Dij2345I(3,1)+p3sq*EijI(13,3)+4*Eij
     -   I(24,3)-2*(EijI(23,3)-EijI(39,4)+EijI(45,4))-EijI(4,1)*P(1)+
     -   (-EijI(3,1)+EijI(9,2))*P(2)+EijI(4,3)*P(5)-EijI(3,2)*P(7)-Ei
     -   jI(6,2)*P(9)-EijI(14,3)*P(14)-EijI(19,3)*P(60)+EijI(8,2)*P(2
     -   43)-EijI(16,3)*P(264)+EijI(10,2)*P(278)-EijI(4,2)*P(280))
       F(213)=DCMPLX(FR(213),FI(213))
       FR(214) = 8*(p3sq*(EijR(10,2)+EijR(16,3))-2*(EijR(24,3)+EijR(39
     -   ,4))-EijR(4,3)*P(5)+EijR(14,3)*P(14)-EijR(9,2)*P(16)+EijR(4,
     -   2)*P(35))
       FI(214) = 8*(p3sq*(EijI(10,2)+EijI(16,3))-2*(EijI(24,3)+EijI(39
     -   ,4))-EijI(4,3)*P(5)+EijI(14,3)*P(14)-EijI(9,2)*P(16)+EijI(4,
     -   2)*P(35))
       F(214)=DCMPLX(FR(214),FI(214))
       P(281) = -p3sq+s12+s45
       P(282) = p2sq-s12-s15
       P(283) = p2sq-2*s12-s15+s34
       P(284) = -s12+s34+2*P(30)
       FR(215) = 8*(D02345R+Dij2345R(2,1)+s12*EijR(2,1)-p3sq*EijR(13,3
     -   )-2*(EijR(11,2)+EijR(24,3)+EijR(39,4)-EijR(45,4))-EijR(4,3)*
     -   P(5)-EijR(6,2)*P(9)+EijR(14,3)*P(14)+EijR(4,2)*P(36)+EijR(19
     -   ,3)*P(60)-EijR(3,1)*P(114)-EijR(3,2)*P(128)+EijR(8,2)*P(152)
     -   +EijR(16,3)*P(264)-EijR(4,1)*P(281)-EijR(7,2)*P(282)+EijR(9,
     -   2)*P(283)-EijR(10,2)*P(284))
       FI(215) = 8*(D02345I+Dij2345I(2,1)+s12*EijI(2,1)-p3sq*EijI(13,3
     -   )-2*(EijI(11,2)+EijI(24,3)+EijI(39,4)-EijI(45,4))-EijI(4,3)*
     -   P(5)-EijI(6,2)*P(9)+EijI(14,3)*P(14)+EijI(4,2)*P(36)+EijI(19
     -   ,3)*P(60)-EijI(3,1)*P(114)-EijI(3,2)*P(128)+EijI(8,2)*P(152)
     -   +EijI(16,3)*P(264)-EijI(4,1)*P(281)-EijI(7,2)*P(282)+EijI(9,
     -   2)*P(283)-EijI(10,2)*P(284))
       F(215)=DCMPLX(FR(215),FI(215))
       P(285) = -p2sq-p3sq+p4sq+s12+s23
       P(286) = p2sq-s23-s45
       P(287) = p2sq+p3sq-p4sq-s23-s45
       P(288) = 2*p2sq+p3sq-p4sq-s12-s23-s45
       FR(216) = 8*(D02345R+Dij2345R(2,1)-Dij2345R(3,2)+Dij2345R(6,2)-
     -   p3sq*EijR(3,2)+s12*(-EijR(5,2)+EijR(9,2))-2*(2*(EijR(11,2)+E
     -   ijR(23,3)-EijR(24,3))+EijR(39,4)+EijR(41,4)-EijR(42,4)-EijR(
     -   45,4))+EijR(11,3)*P(7)+EijR(6,3)*P(9)+EijR(7,3)*P(47)+EijR(6
     -   ,2)*P(156)+EijR(7,2)*P(157)+(-EijR(15,3)-EijR(17,3)+EijR(18,
     -   3)+EijR(20,3))*P(243)+EijR(14,3)*P(244)+EijR(13,3)*P(247)+Ei
     -   jR(4,2)*P(268)+EijR(4,3)*P(285)-EijR(10,2)*P(286)+EijR(16,3)
     -   *P(287)-EijR(19,3)*P(288))
       FI(216) = 8*(D02345I+Dij2345I(2,1)-Dij2345I(3,2)+Dij2345I(6,2)-
     -   p3sq*EijI(3,2)+s12*(-EijI(5,2)+EijI(9,2))-2*(2*(EijI(11,2)+E
     -   ijI(23,3)-EijI(24,3))+EijI(39,4)+EijI(41,4)-EijI(42,4)-EijI(
     -   45,4))+EijI(11,3)*P(7)+EijI(6,3)*P(9)+EijI(7,3)*P(47)+EijI(6
     -   ,2)*P(156)+EijI(7,2)*P(157)+(-EijI(15,3)-EijI(17,3)+EijI(18,
     -   3)+EijI(20,3))*P(243)+EijI(14,3)*P(244)+EijI(13,3)*P(247)+Ei
     -   jI(4,2)*P(268)+EijI(4,3)*P(285)-EijI(10,2)*P(286)+EijI(16,3)
     -   *P(287)-EijI(19,3)*P(288))
       F(216)=DCMPLX(FR(216),FI(216))
       P(289) = 2*s12+s15-s34
       P(290) = p3sq-s12+s34-s45
       P(291) = p2sq-p3sq-s12-s23+s34+3*s45
       P(292) = 3*s12+s15-s34
       P(293) = p2sq-s23+3*s45-2*P(106)
       P(294) = -p2sq+2*p3sq+s12+s15+s23-s34-3*s45
       P(295) = p2sq-2*s12-s23+3*s45
       P(296) = p2sq-2*s12-s23+s34+3*s45
       P(297) = p2sq+s34+4*s45-2*P(48)
       FR(217) = -8*(D02345R+Dij2345R(2,1)-Dij2345R(3,2)+Dij2345R(6,2)
     -   +s15*(EE0R+2*EijR(1,1)+EijR(1,2))-p3sq*EijR(3,2)-2*(3*EijR(1
     -   1,2)+EijR(21,3)-EijR(24,3)-EijR(39,4)-EijR(41,4)+EijR(42,4)+
     -   EijR(45,4))+(EijR(11,3)-EijR(13,3))*P(7)+(EijR(6,3)-EijR(7,3
     -   ))*P(13)+EijR(3,1)*P(128)+(-EijR(15,3)-EijR(17,3)+EijR(18,3)
     -   +EijR(20,3))*P(152)-EijR(4,3)*P(250)-EijR(14,3)*P(254)-EijR(
     -   2,1)*P(289)-EijR(4,1)*P(290)-EijR(4,2)*P(291)+(-EijR(5,2)+Ei
     -   jR(9,2))*P(292)-EijR(6,2)*P(293)-EijR(7,2)*P(294)+EijR(10,2)
     -   *P(295)-EijR(16,3)*P(296)+EijR(19,3)*P(297))
       FI(217) = -8*(D02345I+Dij2345I(2,1)-Dij2345I(3,2)+Dij2345I(6,2)
     -   +s15*(EE0I+2*EijI(1,1)+EijI(1,2))-p3sq*EijI(3,2)-2*(3*EijI(1
     -   1,2)+EijI(21,3)-EijI(24,3)-EijI(39,4)-EijI(41,4)+EijI(42,4)+
     -   EijI(45,4))+(EijI(11,3)-EijI(13,3))*P(7)+(EijI(6,3)-EijI(7,3
     -   ))*P(13)+EijI(3,1)*P(128)+(-EijI(15,3)-EijI(17,3)+EijI(18,3)
     -   +EijI(20,3))*P(152)-EijI(4,3)*P(250)-EijI(14,3)*P(254)-EijI(
     -   2,1)*P(289)-EijI(4,1)*P(290)-EijI(4,2)*P(291)+(-EijI(5,2)+Ei
     -   jI(9,2))*P(292)-EijI(6,2)*P(293)-EijI(7,2)*P(294)+EijI(10,2)
     -   *P(295)-EijI(16,3)*P(296)+EijI(19,3)*P(297))
       F(217)=DCMPLX(FR(217),FI(217))
       P(298) = p2sq-s12-s23+s34+2*s45
       P(299) = 2*p3sq-s45
       FR(218) = 8*(Dij2345R(3,1)+Dij2345R(3,2)-p3sq*(EijR(3,1)+EijR(6
     -   ,2))+2*(EijR(11,2)+EijR(21,3)-EijR(24,3)-EijR(39,4)+EijR(42,
     -   4))+(-EijR(16,3)+EijR(19,3))*P(7)+EijR(7,3)*P(13)+(EijR(2,1)
     -   +EijR(5,2)-EijR(9,2))*P(16)-EijR(4,1)*P(35)-EijR(4,2)*P(114)
     -   +(EijR(15,3)-EijR(18,3))*P(152)+EijR(14,3)*P(254)-EijR(4,3)*
     -   P(298)+EijR(7,2)*P(299))
       FI(218) = 8*(Dij2345I(3,1)+Dij2345I(3,2)-p3sq*(EijI(3,1)+EijI(6
     -   ,2))+2*(EijI(11,2)+EijI(21,3)-EijI(24,3)-EijI(39,4)+EijI(42,
     -   4))+(-EijI(16,3)+EijI(19,3))*P(7)+EijI(7,3)*P(13)+(EijI(2,1)
     -   +EijI(5,2)-EijI(9,2))*P(16)-EijI(4,1)*P(35)-EijI(4,2)*P(114)
     -   +(EijI(15,3)-EijI(18,3))*P(152)+EijI(14,3)*P(254)-EijI(4,3)*
     -   P(298)+EijI(7,2)*P(299))
       F(218)=DCMPLX(FR(218),FI(218))
       P(300) = p3sq-s12+2*P(35)
       FR(219) = 8*(-Dij2345R(3,1)-Dij2345R(3,2)+p3sq*EijR(10,2)+4*Eij
     -   R(24,3)-2*(EijR(39,4)-EijR(42,4))+EijR(16,3)*P(7)+EijR(7,2)*
     -   P(40)+EijR(7,3)*P(47)-EijR(4,3)*P(241)+(-EijR(15,3)+EijR(18,
     -   3))*P(243)+EijR(14,3)*P(244)+EijR(19,3)*P(247)+EijR(4,2)*P(3
     -   00))
       FI(219) = 8*(-Dij2345I(3,1)-Dij2345I(3,2)+p3sq*EijI(10,2)+4*Eij
     -   I(24,3)-2*(EijI(39,4)-EijI(42,4))+EijI(16,3)*P(7)+EijI(7,2)*
     -   P(40)+EijI(7,3)*P(47)-EijI(4,3)*P(241)+(-EijI(15,3)+EijI(18,
     -   3))*P(243)+EijI(14,3)*P(244)+EijI(19,3)*P(247)+EijI(4,2)*P(3
     -   00))
       F(219)=DCMPLX(FR(219),FI(219))
       P(301) = s12+s34+s45
       P(302) = p2sq*s45+p3sq*P(10)+p4sq*P(210)-s45*P(301)
       P(303) = 8*p2sq+5*p3sq-p4sq-7*s12-s15-6*s23+3*s34+2*s45
       FR(220) = 4*(-4*EijR(46,4)+Dij2345R(2,1)*P(169)+D02345R*P(197)+
     -   P(152)*(2*EijR(22,3)-EijR(2,2)*P(198))+s15*(EE0R*P(2)+EijR(1
     -   ,2)*P(9)+EijR(1,1)*P(199))-EijR(2,1)*P(201)-EijR(3,1)*P(205)
     -   -EijR(3,2)*P(209)-EijR(5,2)*P(217)-EijR(6,2)*P(221)+EijR(7,2
     -   )*P(224)-EijR(8,2)*P(229)-P(198)*(Dij2345R(1,1)-Dij2345R(3,1
     -   )+EijR(4,2)*P(213)+EijR(9,2)*P(230))+EijR(10,2)*P(236)+EijR(
     -   4,1)*P(302)+2*(Dij2345R(7,2)-EijR(23,3)*P(7)-EijR(21,3)*P(13
     -   )-EijR(24,3)*P(298)-EijR(11,2)*P(303)))
       FI(220) = 4*(-4*EijI(46,4)+Dij2345I(2,1)*P(169)+D02345I*P(197)+
     -   P(152)*(2*EijI(22,3)-EijI(2,2)*P(198))+s15*(EE0I*P(2)+EijI(1
     -   ,2)*P(9)+EijI(1,1)*P(199))-EijI(2,1)*P(201)-EijI(3,1)*P(205)
     -   -EijI(3,2)*P(209)-EijI(5,2)*P(217)-EijI(6,2)*P(221)+EijI(7,2
     -   )*P(224)-EijI(8,2)*P(229)-P(198)*(Dij2345I(1,1)-Dij2345I(3,1
     -   )+EijI(4,2)*P(213)+EijI(9,2)*P(230))+EijI(10,2)*P(236)+EijI(
     -   4,1)*P(302)+2*(Dij2345I(7,2)-EijI(23,3)*P(7)-EijI(21,3)*P(13
     -   )-EijI(24,3)*P(298)-EijI(11,2)*P(303)))
       F(220)=DCMPLX(FR(220),FI(220))
       P(304) = p2sq+p4sq-s12-s15+s34-s45
       P(305) = -2*p4sq+s45
       P(306) = s15-s23+s45
       P(307) = p2sq+3*p4sq-s12-s15+s34-2*s45
       P(308) = -p2sq+p4sq+s12+s23-s34-s45
       P(309) = p4sq+2*P(119)
       P(310) = p2sq+2*p4sq-s12-s15+s34
       FR(221) = 8*(D02345R+Dij2345R(1,1)-Dij2345R(3,1)+p4sq*EijR(4,3)
     -   -2*(4*(EijR(22,3)-EijR(24,3))+EijR(37,4)+EijR(39,4)+2*(EijR(
     -   11,2)-EijR(44,4)))+(-EijR(2,1)+EijR(4,1))*P(11)-EijR(2,3)*P(
     -   119)+(EijR(8,2)+EijR(9,3)+EijR(16,3))*P(239)-EijR(2,2)*P(304
     -   )+EijR(4,2)*P(305)-(EijR(5,2)-EijR(7,2)+EijR(8,3)+EijR(14,3)
     -   -2*EijR(18,3))*P(306)+EijR(9,2)*P(307)+(EijR(10,2)+2*EijR(20
     -   ,3))*P(308)+EijR(10,3)*P(309)-EijR(15,3)*P(310))
       FI(221) = 8*(D02345I+Dij2345I(1,1)-Dij2345I(3,1)+p4sq*EijI(4,3)
     -   -2*(4*(EijI(22,3)-EijI(24,3))+EijI(37,4)+EijI(39,4)+2*(EijI(
     -   11,2)-EijI(44,4)))+(-EijI(2,1)+EijI(4,1))*P(11)-EijI(2,3)*P(
     -   119)+(EijI(8,2)+EijI(9,3)+EijI(16,3))*P(239)-EijI(2,2)*P(304
     -   )+EijI(4,2)*P(305)-(EijI(5,2)-EijI(7,2)+EijI(8,3)+EijI(14,3)
     -   -2*EijI(18,3))*P(306)+EijI(9,2)*P(307)+(EijI(10,2)+2*EijI(20
     -   ,3))*P(308)+EijI(10,3)*P(309)-EijI(15,3)*P(310))
       F(221)=DCMPLX(FR(221),FI(221))
       P(311) = s15-2*P(13)
       P(312) = p3sq-p4sq+s12+s15-s23-s34+s45
       P(313) = p2sq-p3sq+p4sq-s15+s23+s34-s45
       P(314) = p2sq-p3sq-3*p4sq-s23+s34+s45
       P(315) = p2sq+p4sq-s12-s15+s34
       P(316) = p2sq-2*p4sq-s12-s23+s34+s45
       P(317) = 2*p4sq-s15+s23-s45
       FR(222) = 8*(D02345R-p2sq*EijR(2,2)+p4sq*(-2*EijR(4,2)+EijR(4,3
     -   ))-2*(4*EijR(11,2)+EijR(22,3)+2*EijR(23,3)-3*EijR(24,3)+EijR
     -   (39,4)+EijR(43,4)-EijR(44,4)-EijR(45,4))+EijR(5,2)*P(9)+(-Ei
     -   jR(3,1)+EijR(4,1))*P(11)+(-EijR(9,3)+EijR(10,3))*P(119)-EijR
     -   (6,2)*P(135)-EijR(3,2)*P(238)+EijR(12,3)*P(239)+(-EijR(14,3)
     -   -EijR(17,3)+EijR(18,3)+EijR(19,3))*P(306)+EijR(13,3)*P(308)+
     -   EijR(7,2)*P(311)+EijR(8,2)*P(312)+EijR(9,2)*P(313)-EijR(10,2
     -   )*P(314)-EijR(15,3)*P(315)+EijR(16,3)*P(316)+EijR(20,3)*P(31
     -   7))
       FI(222) = 8*(D02345I-p2sq*EijI(2,2)+p4sq*(-2*EijI(4,2)+EijI(4,3
     -   ))-2*(4*EijI(11,2)+EijI(22,3)+2*EijI(23,3)-3*EijI(24,3)+EijI
     -   (39,4)+EijI(43,4)-EijI(44,4)-EijI(45,4))+EijI(5,2)*P(9)+(-Ei
     -   jI(3,1)+EijI(4,1))*P(11)+(-EijI(9,3)+EijI(10,3))*P(119)-EijI
     -   (6,2)*P(135)-EijI(3,2)*P(238)+EijI(12,3)*P(239)+(-EijI(14,3)
     -   -EijI(17,3)+EijI(18,3)+EijI(19,3))*P(306)+EijI(13,3)*P(308)+
     -   EijI(7,2)*P(311)+EijI(8,2)*P(312)+EijI(9,2)*P(313)-EijI(10,2
     -   )*P(314)-EijI(15,3)*P(315)+EijI(16,3)*P(316)+EijI(20,3)*P(31
     -   7))
       F(222)=DCMPLX(FR(222),FI(222))
       P(318) = -p4sq-s12+s34+s45+2*P(35)
       FR(223) = 8*(D02345R-p2sq*EijR(2,2)+p4sq*(-2*EijR(4,2)+EijR(4,3
     -   ))-2*(2*EijR(11,2)+EijR(22,3)-3*EijR(24,3)+EijR(39,4)-EijR(4
     -   4,4))+(EijR(3,1)-EijR(4,1))*P(5)+EijR(5,2)*P(9)+EijR(8,2)*P(
     -   16)+EijR(10,3)*P(119)+(EijR(16,3)-EijR(20,3))*P(239)+(-EijR(
     -   14,3)+EijR(18,3))*P(306)+EijR(7,2)*P(311)+EijR(9,2)*P(313)-E
     -   ijR(15,3)*P(315)-EijR(10,2)*P(318))
       FI(223) = 8*(D02345I-p2sq*EijI(2,2)+p4sq*(-2*EijI(4,2)+EijI(4,3
     -   ))-2*(2*EijI(11,2)+EijI(22,3)-3*EijI(24,3)+EijI(39,4)-EijI(4
     -   4,4))+(EijI(3,1)-EijI(4,1))*P(5)+EijI(5,2)*P(9)+EijI(8,2)*P(
     -   16)+EijI(10,3)*P(119)+(EijI(16,3)-EijI(20,3))*P(239)+(-EijI(
     -   14,3)+EijI(18,3))*P(306)+EijI(7,2)*P(311)+EijI(9,2)*P(313)-E
     -   ijI(15,3)*P(315)-EijI(10,2)*P(318))
       F(223)=DCMPLX(FR(223),FI(223))
       P(319) = -p4sq+s34-2*P(128)
       P(320) = p4sq-s34+s45-2*P(214)
       P(321) = p4sq+s45-2*P(214)
       P(322) = p3sq+s12-s34
       P(323) = 3*P(11)+2*P(322)
       P(324) = p4sq+s15-s23+s45
       P(325) = -p2sq+s12+s15-s34
       P(326) = p2sq-2*p4sq-s12-s15+s34
       FR(224) = 8*(D02345R+Dij2345R(2,1)-Dij2345R(3,1)+p4sq*EijR(4,3)
     -   -s12*EijR(5,2)-2*(2*EijR(21,3)+EijR(23,3)+3*(EijR(11,2)-EijR
     -   (24,3))+EijR(39,4)+EijR(41,4)-EijR(42,4)-EijR(45,4))+(EijR(1
     -   8,3)+EijR(20,3))*P(119)+(-EijR(3,1)+EijR(4,1))*P(130)-EijR(9
     -   ,2)*P(141)+(EijR(11,3)-EijR(13,3))*P(239)+EijR(8,2)*P(289)+E
     -   ijR(4,2)*P(305)+(-EijR(6,3)+EijR(7,3))*P(306)+EijR(16,3)*P(3
     -   16)+EijR(3,2)*P(319)-EijR(6,2)*P(320)+EijR(7,2)*P(321)+EijR(
     -   10,2)*P(323)-EijR(14,3)*P(324)+(EijR(15,3)+EijR(17,3))*P(325
     -   )-EijR(19,3)*P(326))
       FI(224) = 8*(D02345I+Dij2345I(2,1)-Dij2345I(3,1)+p4sq*EijI(4,3)
     -   -s12*EijI(5,2)-2*(2*EijI(21,3)+EijI(23,3)+3*(EijI(11,2)-EijI
     -   (24,3))+EijI(39,4)+EijI(41,4)-EijI(42,4)-EijI(45,4))+(EijI(1
     -   8,3)+EijI(20,3))*P(119)+(-EijI(3,1)+EijI(4,1))*P(130)-EijI(9
     -   ,2)*P(141)+(EijI(11,3)-EijI(13,3))*P(239)+EijI(8,2)*P(289)+E
     -   ijI(4,2)*P(305)+(-EijI(6,3)+EijI(7,3))*P(306)+EijI(16,3)*P(3
     -   16)+EijI(3,2)*P(319)-EijI(6,2)*P(320)+EijI(7,2)*P(321)+EijI(
     -   10,2)*P(323)-EijI(14,3)*P(324)+(EijI(15,3)+EijI(17,3))*P(325
     -   )-EijI(19,3)*P(326))
       F(224)=DCMPLX(FR(224),FI(224))
       P(327) = p4sq-2*s45
       P(328) = p2sq+p3sq+p4sq-s23-s45
       P(329) = p3sq-p4sq-s12+2*s45
       FR(225) = 8*(Dij2345R(1,1)+Dij2345R(2,1)+p2sq*EijR(2,2)+p3sq*Ei
     -   jR(3,2)+p4sq*EijR(4,3)-2*(Dij2345R(3,1)-2*(EijR(11,2)-EijR(2
     -   2,3)-EijR(23,3))-4*EijR(24,3)+EijR(39,4)+EijR(43,4)-EijR(44,
     -   4)-EijR(45,4))-EijR(5,2)*P(9)+EijR(7,2)*P(13)+EijR(6,2)*P(14
     -   )+(-EijR(9,3)+EijR(10,3))*P(119)-EijR(9,2)*P(200)+EijR(12,3)
     -   *P(239)+(-EijR(14,3)-EijR(17,3)+EijR(18,3)+EijR(19,3))*P(306
     -   )+EijR(13,3)*P(308)-EijR(15,3)*P(315)+EijR(16,3)*P(316)+EijR
     -   (20,3)*P(317)-EijR(4,2)*P(327)-EijR(8,2)*P(328)-EijR(10,2)*P
     -   (329))
       FI(225) = 8*(Dij2345I(1,1)+Dij2345I(2,1)+p2sq*EijI(2,2)+p3sq*Ei
     -   jI(3,2)+p4sq*EijI(4,3)-2*(Dij2345I(3,1)-2*(EijI(11,2)-EijI(2
     -   2,3)-EijI(23,3))-4*EijI(24,3)+EijI(39,4)+EijI(43,4)-EijI(44,
     -   4)-EijI(45,4))-EijI(5,2)*P(9)+EijI(7,2)*P(13)+EijI(6,2)*P(14
     -   )+(-EijI(9,3)+EijI(10,3))*P(119)-EijI(9,2)*P(200)+EijI(12,3)
     -   *P(239)+(-EijI(14,3)-EijI(17,3)+EijI(18,3)+EijI(19,3))*P(306
     -   )+EijI(13,3)*P(308)-EijI(15,3)*P(315)+EijI(16,3)*P(316)+EijI
     -   (20,3)*P(317)-EijI(4,2)*P(327)-EijI(8,2)*P(328)-EijI(10,2)*P
     -   (329))
       F(225)=DCMPLX(FR(225),FI(225))
       P(330) = p2sq-s12-s23+s34+s45
       P(331) = -3*p4sq+2*P(330)
       P(332) = p2sq-3*p4sq-s12-s23+s34+s45
       FR(226) = 8*(Dij2345R(2,1)-Dij2345R(3,1)+p4sq*EijR(4,3)-2*(3*(E
     -   ijR(23,3)-EijR(24,3))+EijR(38,4)+EijR(39,4)-2*EijR(45,4))-(E
     -   ijR(3,2)+EijR(4,2)-2*EijR(10,2))*P(11)-(EijR(12,3)+EijR(15,3
     -   )-2*EijR(20,3))*P(119)+EijR(3,3)*P(239)-(EijR(11,3)+EijR(14,
     -   3)-2*EijR(19,3))*P(306)-EijR(13,3)*P(331)+EijR(16,3)*P(332))
       FI(226) = 8*(Dij2345I(2,1)-Dij2345I(3,1)+p4sq*EijI(4,3)-2*(3*(E
     -   ijI(23,3)-EijI(24,3))+EijI(38,4)+EijI(39,4)-2*EijI(45,4))-(E
     -   ijI(3,2)+EijI(4,2)-2*EijI(10,2))*P(11)-(EijI(12,3)+EijI(15,3
     -   )-2*EijI(20,3))*P(119)+EijI(3,3)*P(239)-(EijI(11,3)+EijI(14,
     -   3)-2*EijI(19,3))*P(306)-EijI(13,3)*P(331)+EijI(16,3)*P(332))
       F(226)=DCMPLX(FR(226),FI(226))
       FR(227) = 8*(Dij2345R(2,1)-Dij2345R(3,1)+s45*EijR(4,2)+p4sq*Eij
     -   R(4,3)+4*EijR(24,3)-2*(s45*EijR(10,2)+EijR(23,3)+EijR(39,4)-
     -   EijR(45,4))+(-EijR(3,1)+EijR(4,1)+EijR(9,2))*P(2)-EijR(3,2)*
     -   P(7)-EijR(6,2)*P(9)-EijR(8,2)*P(26)+(-EijR(15,3)+EijR(20,3))
     -   *P(119)+(-EijR(14,3)+EijR(19,3))*P(306)+EijR(13,3)*P(308)+Ei
     -   jR(16,3)*P(316))
       FI(227) = 8*(Dij2345I(2,1)-Dij2345I(3,1)+s45*EijI(4,2)+p4sq*Eij
     -   I(4,3)+4*EijI(24,3)-2*(s45*EijI(10,2)+EijI(23,3)+EijI(39,4)-
     -   EijI(45,4))+(-EijI(3,1)+EijI(4,1)+EijI(9,2))*P(2)-EijI(3,2)*
     -   P(7)-EijI(6,2)*P(9)-EijI(8,2)*P(26)+(-EijI(15,3)+EijI(20,3))
     -   *P(119)+(-EijI(14,3)+EijI(19,3))*P(306)+EijI(13,3)*P(308)+Ei
     -   jI(16,3)*P(316))
       F(227)=DCMPLX(FR(227),FI(227))
       FR(228) = 8*(-Dij2345R(3,1)+s45*EijR(4,2)+p4sq*EijR(4,3)+4*EijR
     -   (24,3)-2*EijR(39,4)+EijR(9,2)*P(2)-EijR(10,2)*P(5)+EijR(16,3
     -   )*P(239)-EijR(14,3)*P(306)+EijR(15,3)*P(325))
       FI(228) = 8*(-Dij2345I(3,1)+s45*EijI(4,2)+p4sq*EijI(4,3)+4*EijI
     -   (24,3)-2*EijI(39,4)+EijI(9,2)*P(2)-EijI(10,2)*P(5)+EijI(16,3
     -   )*P(239)-EijI(14,3)*P(306)+EijI(15,3)*P(325))
       F(228)=DCMPLX(FR(228),FI(228))
       P(333) = p2sq-s45+2*P(26)
       FR(229) = 8*(Dij2345R(1,1)+p4sq*EijR(4,3)-2*(Dij2345R(3,1)-s45*
     -   EijR(4,2)-EijR(11,2)+EijR(22,3)-3*EijR(24,3)+EijR(39,4)-EijR
     -   (44,4))+(-EijR(2,1)+EijR(4,1))*P(2)-EijR(8,2)*P(7)-EijR(5,2)
     -   *P(9)+EijR(7,2)*P(13)-EijR(2,2)*P(26)-EijR(10,2)*P(38)+EijR(
     -   10,3)*P(119)+EijR(16,3)*P(239)+(-EijR(14,3)+EijR(18,3))*P(30
     -   6)+EijR(20,3)*P(308)-EijR(15,3)*P(315)+EijR(9,2)*P(333))
       FI(229) = 8*(Dij2345I(1,1)+p4sq*EijI(4,3)-2*(Dij2345I(3,1)-s45*
     -   EijI(4,2)-EijI(11,2)+EijI(22,3)-3*EijI(24,3)+EijI(39,4)-EijI
     -   (44,4))+(-EijI(2,1)+EijI(4,1))*P(2)-EijI(8,2)*P(7)-EijI(5,2)
     -   *P(9)+EijI(7,2)*P(13)-EijI(2,2)*P(26)-EijI(10,2)*P(38)+EijI(
     -   10,3)*P(119)+EijI(16,3)*P(239)+(-EijI(14,3)+EijI(18,3))*P(30
     -   6)+EijI(20,3)*P(308)-EijI(15,3)*P(315)+EijI(9,2)*P(333))
       F(229)=DCMPLX(FR(229),FI(229))
       P(334) = p3sq-s12+s15-s23+s34
       P(335) = 2*p3sq+p4sq+s15-s23-s45
       P(336) = p2sq+2*p3sq-p4sq-s12-s23+s34+s45
       P(337) = p2sq+p3sq-2*s12-s15-s23+s34
       P(338) = p2sq-p3sq+p4sq-s15-s45-2*P(55)
       P(339) = s15-s45
       P(340) = -p2sq-p3sq+2*s12+s15+s23-s34+s45
       P(341) = p2sq+p3sq-2*s12-s23+s34-s45
       FR(230) = -8*(Dij2345R(1,1)-Dij2345R(3,2)+Dij2345R(5,2)+s45*Eij
     -   R(4,3)+p3sq*EijR(6,2)+s15*(EijR(5,3)-EijR(7,3))+2*(2*EijR(11
     -   ,2)+EijR(21,3)-EijR(24,3)+EijR(39,4)+EijR(40,4)-EijR(42,4)-E
     -   ijR(44,4)+EijR(4,2)*P(11))-EijR(2,1)*P(55)+EijR(2,2)*P(119)+
     -   (-EijR(16,3)-EijR(17,3)+EijR(19,3)+EijR(20,3))*P(160)+EijR(1
     -   0,2)*P(261)-EijR(4,1)*P(267)+EijR(5,2)*P(334)-EijR(7,2)*P(33
     -   5)-EijR(8,2)*P(336)+(EijR(8,3)-EijR(10,3))*P(337)-EijR(9,2)*
     -   P(338)+EijR(14,3)*P(339)-EijR(15,3)*P(340)-EijR(18,3)*P(341)
     -   )
       FI(230) = -8*(Dij2345I(1,1)-Dij2345I(3,2)+Dij2345I(5,2)+s45*Eij
     -   I(4,3)+p3sq*EijI(6,2)+s15*(EijI(5,3)-EijI(7,3))+2*(2*EijI(11
     -   ,2)+EijI(21,3)-EijI(24,3)+EijI(39,4)+EijI(40,4)-EijI(42,4)-E
     -   ijI(44,4)+EijI(4,2)*P(11))-EijI(2,1)*P(55)+EijI(2,2)*P(119)+
     -   (-EijI(16,3)-EijI(17,3)+EijI(19,3)+EijI(20,3))*P(160)+EijI(1
     -   0,2)*P(261)-EijI(4,1)*P(267)+EijI(5,2)*P(334)-EijI(7,2)*P(33
     -   5)-EijI(8,2)*P(336)+(EijI(8,3)-EijI(10,3))*P(337)-EijI(9,2)*
     -   P(338)+EijI(14,3)*P(339)-EijI(15,3)*P(340)-EijI(18,3)*P(341)
     -   )
       F(230)=DCMPLX(FR(230),FI(230))
       P(342) = p2sq+p3sq-p4sq-s12-s23+s45
       FR(231) = -8*(Dij2345R(3,1)-p4sq*EijR(4,3)+p2sq*EijR(8,2)+2*(Ei
     -   jR(23,3)-3*EijR(24,3)+EijR(39,4)-EijR(45,4))-EijR(4,1)*P(1)-
     -   EijR(3,1)*P(2)-EijR(6,2)*P(9)+EijR(4,2)*P(11)-EijR(3,2)*P(16
     -   )+(EijR(15,3)-EijR(20,3))*P(119)+(EijR(14,3)-EijR(19,3))*P(3
     -   06)-EijR(13,3)*P(308)-EijR(16,3)*P(316)+EijR(10,2)*P(342))
       FI(231) = -8*(Dij2345I(3,1)-p4sq*EijI(4,3)+p2sq*EijI(8,2)+2*(Ei
     -   jI(23,3)-3*EijI(24,3)+EijI(39,4)-EijI(45,4))-EijI(4,1)*P(1)-
     -   EijI(3,1)*P(2)-EijI(6,2)*P(9)+EijI(4,2)*P(11)-EijI(3,2)*P(16
     -   )+(EijI(15,3)-EijI(20,3))*P(119)+(EijI(14,3)-EijI(19,3))*P(3
     -   06)-EijI(13,3)*P(308)-EijI(16,3)*P(316)+EijI(10,2)*P(342))
       F(231)=DCMPLX(FR(231),FI(231))
       P(343) = 2*p3sq+s15-s23+s34-s45
       P(344) = 2*p3sq+s15-s23-s45
       P(345) = p3sq-p4sq
       P(346) = p2sq-s12-s23+3*s45+2*P(345)
       P(347) = p2sq+p3sq-p4sq-2*s12-s23+s34+3*s45
       P(348) = -p2sq-p3sq+p4sq+2*s12+s15+s23-s34-3*s45
       FR(232) = -8*(D02345R+Dij2345R(2,1)-Dij2345R(3,2)+Dij2345R(6,2)
     -   +s45*EijR(4,3)-s12*EijR(5,2)+s15*(EijR(6,3)-EijR(7,3))-2*(2*
     -   EijR(11,2)-EijR(39,4)-EijR(41,4)+EijR(42,4)+EijR(45,4))-EijR
     -   (3,1)*P(55)+EijR(8,2)*P(119)+(-EijR(11,3)+EijR(13,3))*P(160)
     -   -EijR(4,1)*P(267)-EijR(9,2)*P(283)+EijR(4,2)*P(327)-EijR(3,2
     -   )*P(336)+(EijR(15,3)+EijR(17,3)-EijR(18,3)-EijR(20,3))*P(337
     -   )+EijR(14,3)*P(339)+EijR(6,2)*P(343)-EijR(7,2)*P(344)+EijR(1
     -   0,2)*P(346)-EijR(16,3)*P(347)-EijR(19,3)*P(348))
       FI(232) = -8*(D02345I+Dij2345I(2,1)-Dij2345I(3,2)+Dij2345I(6,2)
     -   +s45*EijI(4,3)-s12*EijI(5,2)+s15*(EijI(6,3)-EijI(7,3))-2*(2*
     -   EijI(11,2)-EijI(39,4)-EijI(41,4)+EijI(42,4)+EijI(45,4))-EijI
     -   (3,1)*P(55)+EijI(8,2)*P(119)+(-EijI(11,3)+EijI(13,3))*P(160)
     -   -EijI(4,1)*P(267)-EijI(9,2)*P(283)+EijI(4,2)*P(327)-EijI(3,2
     -   )*P(336)+(EijI(15,3)+EijI(17,3)-EijI(18,3)-EijI(20,3))*P(337
     -   )+EijI(14,3)*P(339)+EijI(6,2)*P(343)-EijI(7,2)*P(344)+EijI(1
     -   0,2)*P(346)-EijI(16,3)*P(347)-EijI(19,3)*P(348))
       F(232)=DCMPLX(FR(232),FI(232))
       P(349) = -p4sq+s12+s45
       P(350) = p4sq-s34+2*P(7)
       P(351) = 2*p4sq+s12+s15-s34-3*s45
       P(352) = p2sq+p4sq-s12+s15-s23
       P(353) = -2*p3sq+3*s15
       P(354) = p3sq+p4sq-s15-s45
       P(355) = p2sq-3*s12-s23+s34+2*P(354)
       P(356) = -p2sq+s23+4*P(7)+2*P(79)
       P(357) = p2sq+p4sq-s12-s23+s34+s45
       FR(233) = -8*(-Dij2345R(1,1)+s15*(EE0R+2*EijR(1,1)+EijR(1,2))-p
     -   4sq*EijR(4,3)+2*(Dij2345R(3,1)-EijR(11,2)+2*EijR(21,3)+EijR(
     -   22,3)-3*EijR(24,3)+EijR(39,4)+EijR(40,4)-EijR(42,4)-EijR(44,
     -   4))-EijR(10,3)*P(119)-EijR(6,2)*P(163)+(-EijR(16,3)-EijR(17,
     -   3)+EijR(19,3)+EijR(20,3))*P(239)-EijR(2,2)*P(289)+(EijR(5,3)
     -   -EijR(7,3))*P(306)+EijR(15,3)*P(315)-EijR(8,2)*P(319)+EijR(1
     -   4,3)*P(324)-EijR(8,3)*P(325)+EijR(4,2)*P(327)-EijR(2,1)*P(34
     -   9)+EijR(3,1)*P(350)-EijR(4,1)*P(351)+EijR(5,2)*P(352)-EijR(7
     -   ,2)*P(353)-EijR(9,2)*P(355)-EijR(10,2)*P(356)-EijR(18,3)*P(3
     -   57))
       FI(233) = -8*(-Dij2345I(1,1)+s15*(EE0I+2*EijI(1,1)+EijI(1,2))-p
     -   4sq*EijI(4,3)+2*(Dij2345I(3,1)-EijI(11,2)+2*EijI(21,3)+EijI(
     -   22,3)-3*EijI(24,3)+EijI(39,4)+EijI(40,4)-EijI(42,4)-EijI(44,
     -   4))-EijI(10,3)*P(119)-EijI(6,2)*P(163)+(-EijI(16,3)-EijI(17,
     -   3)+EijI(19,3)+EijI(20,3))*P(239)-EijI(2,2)*P(289)+(EijI(5,3)
     -   -EijI(7,3))*P(306)+EijI(15,3)*P(315)-EijI(8,2)*P(319)+EijI(1
     -   4,3)*P(324)-EijI(8,3)*P(325)+EijI(4,2)*P(327)-EijI(2,1)*P(34
     -   9)+EijI(3,1)*P(350)-EijI(4,1)*P(351)+EijI(5,2)*P(352)-EijI(7
     -   ,2)*P(353)-EijI(9,2)*P(355)-EijI(10,2)*P(356)-EijI(18,3)*P(3
     -   57))
       F(233)=DCMPLX(FR(233),FI(233))
       P(358) = -p4sq+s34+2*P(5)
       FR(234) = 8*(Dij2345R(3,2)-s45*EijR(4,3)+s12*(EijR(2,1)+EijR(5,
     -   2))+s15*EijR(7,3)-2*(EijR(21,3)-EijR(24,3)+EijR(39,4)-EijR(4
     -   2,4))+(EijR(3,1)+EijR(6,2))*P(14)-EijR(4,1)*P(30)+(EijR(16,3
     -   )-EijR(19,3))*P(160)+EijR(9,2)*P(283)-EijR(4,2)*P(327)+(-Eij
     -   R(15,3)+EijR(18,3))*P(337)-EijR(14,3)*P(339)+EijR(7,2)*P(344
     -   )-EijR(10,2)*P(358))
       FI(234) = 8*(Dij2345I(3,2)-s45*EijI(4,3)+s12*(EijI(2,1)+EijI(5,
     -   2))+s15*EijI(7,3)-2*(EijI(21,3)-EijI(24,3)+EijI(39,4)-EijI(4
     -   2,4))+(EijI(3,1)+EijI(6,2))*P(14)-EijI(4,1)*P(30)+(EijI(16,3
     -   )-EijI(19,3))*P(160)+EijI(9,2)*P(283)-EijI(4,2)*P(327)+(-Eij
     -   I(15,3)+EijI(18,3))*P(337)-EijI(14,3)*P(339)+EijI(7,2)*P(344
     -   )-EijI(10,2)*P(358))
       F(234)=DCMPLX(FR(234),FI(234))
       P(359) = p4sq+s15-2*s45
       P(360) = p4sq-3*s45
       P(361) = s15-s34
       P(362) = p2sq+p3sq-4*s12-s23-2*P(361)
       P(363) = p2sq+p3sq-s23-4*P(7)-2*P(79)
       P(364) = p4sq+3*P(339)
       P(365) = 2*s15-s45
       P(366) = s15-2*s45
       FR(235) = -8*(-Dij2345R(3,2)+s15*(EE0R+3*(EijR(1,1)+EijR(1,2))+
     -   EijR(1,3))+s45*EijR(4,3)+2*(-EijR(11,2)+EijR(36,4)+EijR(39,4
     -   )-2*EijR(42,4))-(EijR(6,3)+EijR(16,3)-2*EijR(19,3))*P(160)-E
     -   ijR(2,1)*P(289)+(EijR(5,3)+EijR(15,3)-2*EijR(18,3))*P(337)+E
     -   ijR(3,1)*P(350)-EijR(4,1)*P(359)+EijR(4,2)*P(360)+(EijR(5,2)
     -   -EijR(9,2))*P(362)+(-EijR(6,2)+EijR(10,2))*P(363)-EijR(7,2)*
     -   P(364)-EijR(7,3)*P(365)+EijR(14,3)*P(366))
       FI(235) = -8*(-Dij2345I(3,2)+s15*(EE0I+3*(EijI(1,1)+EijI(1,2))+
     -   EijI(1,3))+s45*EijI(4,3)+2*(-EijI(11,2)+EijI(36,4)+EijI(39,4
     -   )-2*EijI(42,4))-(EijI(6,3)+EijI(16,3)-2*EijI(19,3))*P(160)-E
     -   ijI(2,1)*P(289)+(EijI(5,3)+EijI(15,3)-2*EijI(18,3))*P(337)+E
     -   ijI(3,1)*P(350)-EijI(4,1)*P(359)+EijI(4,2)*P(360)+(EijI(5,2)
     -   -EijI(9,2))*P(362)+(-EijI(6,2)+EijI(10,2))*P(363)-EijI(7,2)*
     -   P(364)-EijI(7,3)*P(365)+EijI(14,3)*P(366))
       F(235)=DCMPLX(FR(235),FI(235))
       P(367) = p2sq+s15-s23
       P(368) = -s45+2*P(214)
       P(369) = p3sq+p4sq-s34+2*P(7)
       FR(236) = 8*(-Dij2345R(3,1)+p4sq*EijR(4,3)+p3sq*(EijR(3,1)+EijR
     -   (6,2))-2*(EijR(11,2)+EijR(21,3)-2*EijR(24,3)+EijR(39,4)-EijR
     -   (42,4))-EijR(4,2)*P(11)-(EijR(2,1)+EijR(5,2))*P(16)+(-EijR(1
     -   5,3)+EijR(18,3))*P(119)+EijR(16,3)*P(239)+EijR(7,3)*P(306)+E
     -   ijR(19,3)*P(308)-EijR(14,3)*P(324)+EijR(9,2)*P(337)+EijR(4,1
     -   )*P(367)-EijR(7,2)*P(368)+EijR(10,2)*P(369))
       FI(236) = 8*(-Dij2345I(3,1)+p4sq*EijI(4,3)+p3sq*(EijI(3,1)+EijI
     -   (6,2))-2*(EijI(11,2)+EijI(21,3)-2*EijI(24,3)+EijI(39,4)-EijI
     -   (42,4))-EijI(4,2)*P(11)-(EijI(2,1)+EijI(5,2))*P(16)+(-EijI(1
     -   5,3)+EijI(18,3))*P(119)+EijI(16,3)*P(239)+EijI(7,3)*P(306)+E
     -   ijI(19,3)*P(308)-EijI(14,3)*P(324)+EijI(9,2)*P(337)+EijI(4,1
     -   )*P(367)-EijI(7,2)*P(368)+EijI(10,2)*P(369))
       F(236)=DCMPLX(FR(236),FI(236))
       P(370) = 8*p2sq+5*p3sq-p4sq-s15-6*P(48)+2*P(58)
       FR(237) = 4*(EE0R*s15*P(1)+s15*EijR(1,2)*P(47)-Dij2345R(2,1)*P(
     -   169)-D02345R*P(197)-s15*EijR(1,1)*P(199)+EijR(2,1)*P(201)+Ei
     -   jR(3,1)*P(205)+EijR(3,2)*P(209)+EijR(4,1)*P(212)+EijR(5,2)*P
     -   (217)+EijR(6,2)*P(221)-EijR(7,2)*P(224)+EijR(8,2)*P(229)+P(1
     -   98)*(Dij2345R(1,1)-Dij2345R(3,1)+EijR(2,2)*P(152)+EijR(4,2)*
     -   P(213)+EijR(9,2)*P(230))-EijR(10,2)*P(236)+2*(p4sq*EijR(24,3
     -   )-2*EijR(46,4)+EijR(23,3)*P(239)-EijR(21,3)*P(306)+EijR(22,3
     -   )*P(325)+EijR(11,2)*P(370)))
       FI(237) = 4*(EE0I*s15*P(1)+s15*EijI(1,2)*P(47)-Dij2345I(2,1)*P(
     -   169)-D02345I*P(197)-s15*EijI(1,1)*P(199)+EijI(2,1)*P(201)+Ei
     -   jI(3,1)*P(205)+EijI(3,2)*P(209)+EijI(4,1)*P(212)+EijI(5,2)*P
     -   (217)+EijI(6,2)*P(221)-EijI(7,2)*P(224)+EijI(8,2)*P(229)+P(1
     -   98)*(Dij2345I(1,1)-Dij2345I(3,1)+EijI(2,2)*P(152)+EijI(4,2)*
     -   P(213)+EijI(9,2)*P(230))-EijI(10,2)*P(236)+2*(p4sq*EijI(24,3
     -   )-2*EijI(46,4)+EijI(23,3)*P(239)-EijI(21,3)*P(306)+EijI(22,3
     -   )*P(325)+EijI(11,2)*P(370)))
       F(237)=DCMPLX(FR(237),FI(237))
       P(371) = p3sq+p4sq+s12-s34-2*s45
       P(372) = p4sq+s12-s34-2*s45
       P(373) = p2sq+p3sq+p4sq+s12-s34-2*s45
       P(374) = -p2sq+p3sq+p4sq+s12-s34
       P(375) = p3sq-s34
       P(376) = p2sq-p4sq-s12+3*s45-2*P(375)
       FR(238) = 8*(Dij2345R(1,1)-Dij2345R(3,1)+p3sq*(EijR(3,1)+2*EijR
     -   (6,2))+2*(EijR(11,2)-EijR(21,3)-EijR(22,3)+2*EijR(24,3)-EijR
     -   (39,4)-EijR(40,4)+EijR(42,4)+EijR(44,4))-EijR(10,2)*P(38)-Ei
     -   jR(8,2)*P(128)+EijR(2,2)*P(152)+EijR(4,2)*P(213)-EijR(2,1)*P
     -   (371)+EijR(4,1)*P(372)-EijR(5,2)*P(373)-EijR(7,2)*P(374)-Eij
     -   R(9,2)*P(376))
       FI(238) = 8*(Dij2345I(1,1)-Dij2345I(3,1)+p3sq*(EijI(3,1)+2*EijI
     -   (6,2))+2*(EijI(11,2)-EijI(21,3)-EijI(22,3)+2*EijI(24,3)-EijI
     -   (39,4)-EijI(40,4)+EijI(42,4)+EijI(44,4))-EijI(10,2)*P(38)-Ei
     -   jI(8,2)*P(128)+EijI(2,2)*P(152)+EijI(4,2)*P(213)-EijI(2,1)*P
     -   (371)+EijI(4,1)*P(372)-EijI(5,2)*P(373)-EijI(7,2)*P(374)-Eij
     -   I(9,2)*P(376))
       F(238)=DCMPLX(FR(238),FI(238))
       FR(239) = -8*(2*(EijR(22,3)-EijR(24,3)+EijR(37,4)+EijR(39,4)-2*
     -   EijR(44,4))+(EijR(5,2)-EijR(7,2))*P(14)+(EijR(2,1)-EijR(4,1)
     -   )*P(86)+(-EijR(2,2)+EijR(9,2))*P(169))
       FI(239) = -8*(2*(EijI(22,3)-EijI(24,3)+EijI(37,4)+EijI(39,4)-2*
     -   EijI(44,4))+(EijI(5,2)-EijI(7,2))*P(14)+(EijI(2,1)-EijI(4,1)
     -   )*P(86)+(-EijI(2,2)+EijI(9,2))*P(169))
       F(239)=DCMPLX(FR(239),FI(239))
       P(377) = -p4sq-s12+s34+s45+2*P(16)
       P(378) = p2sq-p4sq-3*s12+s45+2*P(18)
       P(379) = 2*p3sq+p4sq-s12-s34+s45
       FR(240) = 8*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)+EijR(5,2)*P(9)
     -   -2*(p3sq*EijR(3,2)+4*EijR(11,2)+EijR(39,4)+EijR(43,4)-EijR(4
     -   4,4)-EijR(45,4)+EijR(6,2)*P(14))+(-EijR(3,1)+EijR(4,1))*P(86
     -   )-EijR(2,2)*P(152)+EijR(7,2)*P(180)-EijR(4,2)*P(213)+EijR(8,
     -   2)*P(377)-EijR(9,2)*P(378)+EijR(10,2)*P(379))
       FI(240) = 8*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)+EijI(5,2)*P(9)
     -   -2*(p3sq*EijI(3,2)+4*EijI(11,2)+EijI(39,4)+EijI(43,4)-EijI(4
     -   4,4)-EijI(45,4)+EijI(6,2)*P(14))+(-EijI(3,1)+EijI(4,1))*P(86
     -   )-EijI(2,2)*P(152)+EijI(7,2)*P(180)-EijI(4,2)*P(213)+EijI(8,
     -   2)*P(377)-EijI(9,2)*P(378)+EijI(10,2)*P(379))
       F(240)=DCMPLX(FR(240),FI(240))
       FR(241) = 8*(-Dij2345R(1,1)+Dij2345R(3,1)-s12*EijR(2,1)+p3sq*Ei
     -   jR(3,1)-4*EijR(11,2)-2*(EijR(39,4)-EijR(44,4))-EijR(4,1)*P(3
     -   )+EijR(5,2)*P(9)+EijR(10,2)*P(38)+EijR(8,2)*P(128)-EijR(2,2)
     -   *P(152)+EijR(7,2)*P(180)-EijR(4,2)*P(213)-EijR(9,2)*P(378))
       FI(241) = 8*(-Dij2345I(1,1)+Dij2345I(3,1)-s12*EijI(2,1)+p3sq*Ei
     -   jI(3,1)-4*EijI(11,2)-2*(EijI(39,4)-EijI(44,4))-EijI(4,1)*P(3
     -   )+EijI(5,2)*P(9)+EijI(10,2)*P(38)+EijI(8,2)*P(128)-EijI(2,2)
     -   *P(152)+EijI(7,2)*P(180)-EijI(4,2)*P(213)-EijI(9,2)*P(378))
       F(241)=DCMPLX(FR(241),FI(241))
       FR(242) = -16*(EijR(2,1)-EijR(4,1)+EijR(4,4)+EijR(5,4)-EijR(7,4
     -   )-EijR(15,4)+6*(EijR(14,3)-EijR(18,3))+3*(EijR(4,2)-EijR(4,3
     -   )+EijR(5,2)+EijR(5,3)-EijR(7,2)-EijR(7,3)-EijR(9,2)-EijR(14,
     -   4)+EijR(15,3)+EijR(19,4)-EijR(24,4)+EijR(32,4)))
       FI(242) = -16*(EijI(2,1)-EijI(4,1)+EijI(4,4)+EijI(5,4)-EijI(7,4
     -   )-EijI(15,4)+6*(EijI(14,3)-EijI(18,3))+3*(EijI(4,2)-EijI(4,3
     -   )+EijI(5,2)+EijI(5,3)-EijI(7,2)-EijI(7,3)-EijI(9,2)-EijI(14,
     -   4)+EijI(15,3)+EijI(19,4)-EijI(24,4)+EijI(32,4)))
       F(242)=DCMPLX(FR(242),FI(242))
       FR(243) = -16*(EijR(2,2)+EijR(4,2)+EijR(4,4)+EijR(17,4)+4*(EijR
     -   (15,3)-EijR(18,3))+EijR(19,4)+EijR(21,4)-2*(EijR(4,3)-EijR(8
     -   ,3)+EijR(9,2)+EijR(10,3)-EijR(14,3)+EijR(14,4)+EijR(15,4)+Ei
     -   jR(24,4)+EijR(27,4)-2*EijR(32,4)))
       FI(243) = -16*(EijI(2,2)+EijI(4,2)+EijI(4,4)+EijI(17,4)+4*(EijI
     -   (15,3)-EijI(18,3))+EijI(19,4)+EijI(21,4)-2*(EijI(4,3)-EijI(8
     -   ,3)+EijI(9,2)+EijI(10,3)-EijI(14,3)+EijI(14,4)+EijI(15,4)+Ei
     -   jI(24,4)+EijI(27,4)-2*EijI(32,4)))
       F(243)=DCMPLX(FR(243),FI(243))
       FR(244) = 16*(EijR(4,3)-EijR(4,4)+EijR(6,2)+EijR(6,3)-EijR(7,2)
     -   -EijR(7,3)-EijR(8,2)+EijR(9,2)+EijR(15,4)-EijR(16,3)+EijR(16
     -   ,4)-EijR(19,4)-EijR(23,4)+EijR(24,4)+EijR(25,4)-EijR(34,4)+2
     -   *(EijR(14,4)-EijR(15,3)-EijR(17,3)+EijR(18,3)+EijR(20,3)-Eij
     -   R(32,4)-EijR(33,4)+EijR(35,4)))
       FI(244) = 16*(EijI(4,3)-EijI(4,4)+EijI(6,2)+EijI(6,3)-EijI(7,2)
     -   -EijI(7,3)-EijI(8,2)+EijI(9,2)+EijI(15,4)-EijI(16,3)+EijI(16
     -   ,4)-EijI(19,4)-EijI(23,4)+EijI(24,4)+EijI(25,4)-EijI(34,4)+2
     -   *(EijI(14,4)-EijI(15,3)-EijI(17,3)+EijI(18,3)+EijI(20,3)-Eij
     -   I(32,4)-EijI(33,4)+EijI(35,4)))
       F(244)=DCMPLX(FR(244),FI(244))
       FR(245) = 16*(EijR(4,3)-EijR(4,4)-EijR(7,2)-EijR(7,3)+EijR(9,2)
     -   +EijR(15,4)-EijR(19,4)+EijR(24,4)+2*(EijR(14,4)-EijR(15,3)+E
     -   ijR(18,3)-EijR(32,4)))
       FI(245) = 16*(EijI(4,3)-EijI(4,4)-EijI(7,2)-EijI(7,3)+EijI(9,2)
     -   +EijI(15,4)-EijI(19,4)+EijI(24,4)+2*(EijI(14,4)-EijI(15,3)+E
     -   ijI(18,3)-EijI(32,4)))
       F(245)=DCMPLX(FR(245),FI(245))
       FR(246) = -16*(EijR(2,1)+EijR(2,2)-EijR(4,1)+3*(EijR(4,2)-EijR(
     -   4,3))+EijR(4,4)+EijR(5,3)-EijR(7,3)-4*(EijR(9,2)-EijR(14,3))
     -   +5*EijR(15,3)+EijR(17,4)-6*EijR(18,3)+EijR(19,4)+EijR(21,4)+
     -   2*(EijR(5,2)-EijR(7,2)+EijR(8,3)-EijR(10,3)-EijR(14,4)-EijR(
     -   15,4)-EijR(24,4)-EijR(27,4)+2*EijR(32,4)))
       FI(246) = -16*(EijI(2,1)+EijI(2,2)-EijI(4,1)+3*(EijI(4,2)-EijI(
     -   4,3))+EijI(4,4)+EijI(5,3)-EijI(7,3)-4*(EijI(9,2)-EijI(14,3))
     -   +5*EijI(15,3)+EijI(17,4)-6*EijI(18,3)+EijI(19,4)+EijI(21,4)+
     -   2*(EijI(5,2)-EijI(7,2)+EijI(8,3)-EijI(10,3)-EijI(14,4)-EijI(
     -   15,4)-EijI(24,4)-EijI(27,4)+2*EijI(32,4)))
       F(246)=DCMPLX(FR(246),FI(246))
       FR(247) = -16*(EijR(4,4)+EijR(5,2)+EijR(5,3)-EijR(7,2)-EijR(7,3
     -   )+EijR(8,2)-EijR(10,2)-3*(EijR(4,3)-EijR(15,3))-EijR(15,4)-E
     -   ijR(16,4)+4*(EijR(14,3)-EijR(18,3))+EijR(19,4)+EijR(23,4)-Ei
     -   jR(24,4)-EijR(25,4)+EijR(34,4)+2*(EijR(4,2)-EijR(9,2)-EijR(1
     -   4,4)+EijR(16,3)+EijR(17,3)-EijR(19,3)-EijR(20,3)+EijR(32,4)+
     -   EijR(33,4)-EijR(35,4)))
       FI(247) = -16*(EijI(4,4)+EijI(5,2)+EijI(5,3)-EijI(7,2)-EijI(7,3
     -   )+EijI(8,2)-EijI(10,2)-3*(EijI(4,3)-EijI(15,3))-EijI(15,4)-E
     -   ijI(16,4)+4*(EijI(14,3)-EijI(18,3))+EijI(19,4)+EijI(23,4)-Ei
     -   jI(24,4)-EijI(25,4)+EijI(34,4)+2*(EijI(4,2)-EijI(9,2)-EijI(1
     -   4,4)+EijI(16,3)+EijI(17,3)-EijI(19,3)-EijI(20,3)+EijI(32,4)+
     -   EijI(33,4)-EijI(35,4)))
       F(247)=DCMPLX(FR(247),FI(247))
       FR(248) = -16*(EijR(4,2)+EijR(4,4)-EijR(9,2)-EijR(15,4)+EijR(19
     -   ,4)-EijR(24,4)-2*(EijR(4,3)-EijR(14,3)+EijR(14,4)-EijR(15,3)
     -   +EijR(18,3)-EijR(32,4)))
       FI(248) = -16*(EijI(4,2)+EijI(4,4)-EijI(9,2)-EijI(15,4)+EijI(19
     -   ,4)-EijI(24,4)-2*(EijI(4,3)-EijI(14,3)+EijI(14,4)-EijI(15,3)
     -   +EijI(18,3)-EijI(32,4)))
       F(248)=DCMPLX(FR(248),FI(248))
       FR(249) = -16*(EijR(2,1)+EijR(2,3)-EijR(4,1)+EijR(4,4)+EijR(5,2
     -   )-EijR(7,2)+EijR(8,4)-5*(EijR(9,2)+EijR(10,3))-EijR(10,4)+2*
     -   (EijR(2,2)+EijR(8,3)+EijR(14,3))-EijR(14,4)+7*EijR(15,3)-4*E
     -   ijR(18,3)+3*(EijR(4,2)-EijR(4,3)-EijR(15,4)+EijR(21,4)-EijR(
     -   27,4)+EijR(32,4)))
       FI(249) = -16*(EijI(2,1)+EijI(2,3)-EijI(4,1)+EijI(4,4)+EijI(5,2
     -   )-EijI(7,2)+EijI(8,4)-5*(EijI(9,2)+EijI(10,3))-EijI(10,4)+2*
     -   (EijI(2,2)+EijI(8,3)+EijI(14,3))-EijI(14,4)+7*EijI(15,3)-4*E
     -   ijI(18,3)+3*(EijI(4,2)-EijI(4,3)-EijI(15,4)+EijI(21,4)-EijI(
     -   27,4)+EijI(32,4)))
       F(249)=DCMPLX(FR(249),FI(249))
       FR(250) = -16*(EijR(2,2)+EijR(2,4)+EijR(4,2)+EijR(4,4)+2*(EijR(
     -   2,3)-EijR(4,3)-EijR(9,2))-4*(EijR(10,4)+EijR(15,4))-6*(EijR(
     -   10,3)-EijR(15,3)-EijR(21,4)))
       FI(250) = -16*(EijI(2,2)+EijI(2,4)+EijI(4,2)+EijI(4,4)+2*(EijI(
     -   2,3)-EijI(4,3)-EijI(9,2))-4*(EijI(10,4)+EijI(15,4))-6*(EijI(
     -   10,3)-EijI(15,3)-EijI(21,4)))
       F(250)=DCMPLX(FR(250),FI(250))
       FR(251) = -16*(EijR(2,1)-EijR(4,1)+EijR(4,4)+EijR(5,2)+EijR(6,2
     -   )+EijR(6,3)-EijR(7,3)+EijR(8,2)-EijR(15,4)+3*(EijR(4,2)-EijR
     -   (4,3)+EijR(16,3))-EijR(16,4)+4*(EijR(14,3)-EijR(19,3))+EijR(
     -   19,4)+EijR(23,4)-EijR(24,4)-EijR(25,4)+EijR(34,4)-2*(EijR(7,
     -   2)+EijR(9,2)+EijR(10,2)+EijR(14,4)-EijR(15,3)-EijR(17,3)+Eij
     -   R(18,3)+EijR(20,3)-EijR(32,4)-EijR(33,4)+EijR(35,4)))
       FI(251) = -16*(EijI(2,1)-EijI(4,1)+EijI(4,4)+EijI(5,2)+EijI(6,2
     -   )+EijI(6,3)-EijI(7,3)+EijI(8,2)-EijI(15,4)+3*(EijI(4,2)-EijI
     -   (4,3)+EijI(16,3))-EijI(16,4)+4*(EijI(14,3)-EijI(19,3))+EijI(
     -   19,4)+EijI(23,4)-EijI(24,4)-EijI(25,4)+EijI(34,4)-2*(EijI(7,
     -   2)+EijI(9,2)+EijI(10,2)+EijI(14,4)-EijI(15,3)-EijI(17,3)+Eij
     -   I(18,3)+EijI(20,3)-EijI(32,4)-EijI(33,4)+EijI(35,4)))
       F(251)=DCMPLX(FR(251),FI(251))
       FR(252) = -16*(EijR(3,1)-EijR(4,1)-EijR(4,3)+EijR(4,4)+EijR(6,2
     -   )+EijR(7,2)-EijR(9,2)-EijR(10,2)-EijR(14,4)+EijR(15,3)-EijR(
     -   15,4)+EijR(16,3)-EijR(16,4)-EijR(20,3)+EijR(32,4)+EijR(33,4)
     -   +EijR(34,4)-EijR(35,4))
       FI(252) = -16*(EijI(3,1)-EijI(4,1)-EijI(4,3)+EijI(4,4)+EijI(6,2
     -   )+EijI(7,2)-EijI(9,2)-EijI(10,2)-EijI(14,4)+EijI(15,3)-EijI(
     -   15,4)+EijI(16,3)-EijI(16,4)-EijI(20,3)+EijI(32,4)+EijI(33,4)
     -   +EijI(34,4)-EijI(35,4))
       F(252)=DCMPLX(FR(252),FI(252))
       FR(253) = -16*(EijR(2,1)-EijR(4,1)+EijR(4,4)+EijR(5,2)+EijR(6,2
     -   )+EijR(8,2)+EijR(9,3)-EijR(10,3)-EijR(14,4)+3*(EijR(4,2)-Eij
     -   R(4,3)+EijR(16,3))-EijR(16,4)+4*(EijR(15,3)-EijR(20,3))+EijR
     -   (21,4)+EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(33,4)-2*(EijR(7
     -   ,2)+EijR(9,2)+EijR(10,2)-EijR(14,3)+EijR(15,4)-EijR(17,3)+Ei
     -   jR(18,3)+EijR(19,3)-EijR(32,4)-EijR(34,4)+EijR(35,4)))
       FI(253) = -16*(EijI(2,1)-EijI(4,1)+EijI(4,4)+EijI(5,2)+EijI(6,2
     -   )+EijI(8,2)+EijI(9,3)-EijI(10,3)-EijI(14,4)+3*(EijI(4,2)-Eij
     -   I(4,3)+EijI(16,3))-EijI(16,4)+4*(EijI(15,3)-EijI(20,3))+EijI
     -   (21,4)+EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(33,4)-2*(EijI(7
     -   ,2)+EijI(9,2)+EijI(10,2)-EijI(14,3)+EijI(15,4)-EijI(17,3)+Ei
     -   jI(18,3)+EijI(19,3)-EijI(32,4)-EijI(34,4)+EijI(35,4)))
       F(253)=DCMPLX(FR(253),FI(253))
       FR(254) = -16*(EijR(4,2)+EijR(4,4)+EijR(8,2)-EijR(9,2)+EijR(9,4
     -   )-EijR(10,2)-EijR(10,4)-2*(EijR(4,3)-EijR(9,3)+EijR(10,3)-Ei
     -   jR(16,3))-EijR(16,4)+4*(EijR(15,3)-EijR(20,3))-3*(EijR(15,4)
     -   -EijR(21,4)+EijR(28,4)-EijR(34,4)))
       FI(254) = -16*(EijI(4,2)+EijI(4,4)+EijI(8,2)-EijI(9,2)+EijI(9,4
     -   )-EijI(10,2)-EijI(10,4)-2*(EijI(4,3)-EijI(9,3)+EijI(10,3)-Ei
     -   jI(16,3))-EijI(16,4)+4*(EijI(15,3)-EijI(20,3))-3*(EijI(15,4)
     -   -EijI(21,4)+EijI(28,4)-EijI(34,4)))
       F(254)=DCMPLX(FR(254),FI(254))
       FR(255) = 16*(EijR(2,1)-EijR(4,1)-EijR(4,2)-EijR(4,4)+EijR(5,2)
     -   +EijR(7,2)-EijR(9,2)-EijR(14,3)+EijR(14,4)-EijR(15,3)+EijR(1
     -   5,4)+2*(EijR(4,3)-EijR(16,3))+EijR(16,4)+EijR(19,3)+EijR(20,
     -   3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4))
       FI(255) = 16*(EijI(2,1)-EijI(4,1)-EijI(4,2)-EijI(4,4)+EijI(5,2)
     -   +EijI(7,2)-EijI(9,2)-EijI(14,3)+EijI(14,4)-EijI(15,3)+EijI(1
     -   5,4)+2*(EijI(4,3)-EijI(16,3))+EijI(16,4)+EijI(19,3)+EijI(20,
     -   3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4))
       F(255)=DCMPLX(FR(255),FI(255))
       FR(256) = 16*(EijR(2,1)-EijR(4,1)-EijR(4,4)+EijR(5,3)-EijR(7,2)
     -   -EijR(9,2)-EijR(15,3)+EijR(15,4)-EijR(19,4)+EijR(24,4)+2*(Ei
     -   jR(4,3)+EijR(5,2)-EijR(14,3)+EijR(14,4)-EijR(32,4)))
       FI(256) = 16*(EijI(2,1)-EijI(4,1)-EijI(4,4)+EijI(5,3)-EijI(7,2)
     -   -EijI(9,2)-EijI(15,3)+EijI(15,4)-EijI(19,4)+EijI(24,4)+2*(Ei
     -   jI(4,3)+EijI(5,2)-EijI(14,3)+EijI(14,4)-EijI(32,4)))
       F(256)=DCMPLX(FR(256),FI(256))
       FR(257) = 16*(EijR(2,1)+EijR(2,2)-EijR(4,1)-EijR(4,4)+EijR(5,2)
     -   +EijR(8,3)-EijR(14,3)+EijR(14,4)-EijR(21,4)+EijR(27,4)+2*(Ei
     -   jR(4,3)-EijR(9,2)-EijR(15,3)+EijR(15,4)-EijR(32,4)))
       FI(257) = 16*(EijI(2,1)+EijI(2,2)-EijI(4,1)-EijI(4,4)+EijI(5,2)
     -   +EijI(8,3)-EijI(14,3)+EijI(14,4)-EijI(21,4)+EijI(27,4)+2*(Ei
     -   jI(4,3)-EijI(9,2)-EijI(15,3)+EijI(15,4)-EijI(32,4)))
       F(257)=DCMPLX(FR(257),FI(257))
       FR(258) = 16*(EijR(2,1)+2*EijR(2,2)+EijR(2,3)-EijR(4,1)+EijR(4,
     -   2)+EijR(4,3)-EijR(4,4)-EijR(10,3)+EijR(10,4)-EijR(15,3)-3*(E
     -   ijR(9,2)-EijR(15,4)+EijR(21,4)))
       FI(258) = 16*(EijI(2,1)+2*EijI(2,2)+EijI(2,3)-EijI(4,1)+EijI(4,
     -   2)+EijI(4,3)-EijI(4,4)-EijI(10,3)+EijI(10,4)-EijI(15,3)-3*(E
     -   ijI(9,2)-EijI(15,4)+EijI(21,4)))
       F(258)=DCMPLX(FR(258),FI(258))
       FR(259) = 16*(EijR(3,1)-EijR(4,1)+EijR(4,2)-EijR(4,4)-EijR(9,2)
     -   +EijR(9,3)-EijR(10,3)+EijR(15,3)+EijR(16,4)-EijR(20,3)-EijR(
     -   21,4)+EijR(28,4)+2*(EijR(8,2)-EijR(10,2)+EijR(15,4)-EijR(34,
     -   4)))
       FI(259) = 16*(EijI(3,1)-EijI(4,1)+EijI(4,2)-EijI(4,4)-EijI(9,2)
     -   +EijI(9,3)-EijI(10,3)+EijI(15,3)+EijI(16,4)-EijI(20,3)-EijI(
     -   21,4)+EijI(28,4)+2*(EijI(8,2)-EijI(10,2)+EijI(15,4)-EijI(34,
     -   4)))
       F(259)=DCMPLX(FR(259),FI(259))
       FR(260) = EijR(2,2)+EijR(2,3)+EijR(4,2)+EijR(4,4)+EijR(8,3)+Eij
     -   R(8,4)-4*EijR(10,3)-EijR(10,4)+EijR(14,3)-EijR(14,4)+5*EijR(
     -   15,3)-2*(EijR(4,3)+EijR(9,2)+EijR(18,3))-3*(EijR(15,4)-EijR(
     -   21,4)+EijR(27,4)-EijR(32,4))
       FI(260) = EijI(2,2)+EijI(2,3)+EijI(4,2)+EijI(4,4)+EijI(8,3)+Eij
     -   I(8,4)-4*EijI(10,3)-EijI(10,4)+EijI(14,3)-EijI(14,4)+5*EijI(
     -   15,3)-2*(EijI(4,3)+EijI(9,2)+EijI(18,3))-3*(EijI(15,4)-EijI(
     -   21,4)+EijI(27,4)-EijI(32,4))
       F(260)=DCMPLX(FR(260),FI(260))
       FR(261) = EijR(4,3)-EijR(4,4)+EijR(6,2)-EijR(7,2)-EijR(8,2)+Eij
     -   R(9,2)-EijR(9,3)+EijR(10,3)+EijR(14,4)-EijR(16,3)+EijR(16,4)
     -   -EijR(21,4)-EijR(26,4)+EijR(27,4)+EijR(28,4)-EijR(33,4)-2*(E
     -   ijR(15,3)-EijR(15,4)-EijR(20,3)+EijR(32,4)+EijR(34,4)-EijR(3
     -   5,4))
       FI(261) = EijI(4,3)-EijI(4,4)+EijI(6,2)-EijI(7,2)-EijI(8,2)+Eij
     -   I(9,2)-EijI(9,3)+EijI(10,3)+EijI(14,4)-EijI(16,3)+EijI(16,4)
     -   -EijI(21,4)-EijI(26,4)+EijI(27,4)+EijI(28,4)-EijI(33,4)-2*(E
     -   ijI(15,3)-EijI(15,4)-EijI(20,3)+EijI(32,4)+EijI(34,4)-EijI(3
     -   5,4))
       F(261)=DCMPLX(FR(261),FI(261))
       FR(262) = EijR(4,3)-EijR(4,4)-EijR(7,2)+EijR(9,2)+EijR(10,3)+Ei
     -   jR(14,4)-EijR(21,4)+EijR(27,4)-2*(EijR(15,3)-EijR(15,4)+EijR
     -   (32,4))
       FI(262) = EijI(4,3)-EijI(4,4)-EijI(7,2)+EijI(9,2)+EijI(10,3)+Ei
     -   jI(14,4)-EijI(21,4)+EijI(27,4)-2*(EijI(15,3)-EijI(15,4)+EijI
     -   (32,4))
       F(262)=DCMPLX(FR(262),FI(262))
       FR(263) = EijR(4,2)+EijR(4,4)+EijR(5,2)-EijR(7,2)+EijR(8,3)-Eij
     -   R(9,2)+EijR(9,3)+EijR(14,3)-EijR(14,4)+4*EijR(15,3)+EijR(16,
     -   3)-EijR(16,4)+EijR(21,4)+EijR(26,4)-EijR(27,4)-EijR(28,4)+Ei
     -   jR(33,4)-2*(EijR(4,3)+EijR(10,3)+EijR(15,4)+EijR(18,3)+EijR(
     -   20,3)-EijR(32,4)-EijR(34,4)+EijR(35,4))
       FI(263) = EijI(4,2)+EijI(4,4)+EijI(5,2)-EijI(7,2)+EijI(8,3)-Eij
     -   I(9,2)+EijI(9,3)+EijI(14,3)-EijI(14,4)+4*EijI(15,3)+EijI(16,
     -   3)-EijI(16,4)+EijI(21,4)+EijI(26,4)-EijI(27,4)-EijI(28,4)+Ei
     -   jI(33,4)-2*(EijI(4,3)+EijI(10,3)+EijI(15,4)+EijI(18,3)+EijI(
     -   20,3)-EijI(32,4)-EijI(34,4)+EijI(35,4))
       F(263)=DCMPLX(FR(263),FI(263))
       FR(264) = -EijR(4,3)+EijR(4,4)-EijR(11,3)+EijR(12,3)-EijR(14,4)
     -   -EijR(15,4)+EijR(16,3)+EijR(17,3)-EijR(18,3)+EijR(19,3)-3*Ei
     -   jR(20,3)+EijR(22,4)+EijR(29,4)-EijR(30,4)-EijR(31,4)+EijR(32
     -   ,4)+2*(EijR(15,3)-EijR(16,4)+EijR(33,4)+EijR(34,4)-EijR(35,4
     -   ))
       FI(264) = -EijI(4,3)+EijI(4,4)-EijI(11,3)+EijI(12,3)-EijI(14,4)
     -   -EijI(15,4)+EijI(16,3)+EijI(17,3)-EijI(18,3)+EijI(19,3)-3*Ei
     -   jI(20,3)+EijI(22,4)+EijI(29,4)-EijI(30,4)-EijI(31,4)+EijI(32
     -   ,4)+2*(EijI(15,3)-EijI(16,4)+EijI(33,4)+EijI(34,4)-EijI(35,4
     -   ))
       F(264)=DCMPLX(FR(264),FI(264))
       FR(265) = EijR(4,2)-EijR(4,3)+EijR(4,4)-EijR(7,2)-EijR(14,4)+2*
     -   EijR(15,3)-EijR(15,4)-EijR(16,4)-EijR(18,3)+EijR(19,3)-EijR(
     -   20,3)+EijR(32,4)+EijR(33,4)+EijR(34,4)-EijR(35,4)
       FI(265) = EijI(4,2)-EijI(4,3)+EijI(4,4)-EijI(7,2)-EijI(14,4)+2*
     -   EijI(15,3)-EijI(15,4)-EijI(16,4)-EijI(18,3)+EijI(19,3)-EijI(
     -   20,3)+EijI(32,4)+EijI(33,4)+EijI(34,4)-EijI(35,4)
       F(265)=DCMPLX(FR(265),FI(265))
       FR(266) = -EijR(4,3)+EijR(4,4)-EijR(10,3)-EijR(14,4)+EijR(21,4)
     -   -EijR(27,4)+2*(EijR(15,3)-EijR(15,4)+EijR(32,4))
       FI(266) = -EijI(4,3)+EijI(4,4)-EijI(10,3)-EijI(14,4)+EijI(21,4)
     -   -EijI(27,4)+2*(EijI(15,3)-EijI(15,4)+EijI(32,4))
       F(266)=DCMPLX(FR(266),FI(266))
       FR(267) = EijR(4,2)-EijR(4,4)-EijR(7,2)+EijR(14,3)+EijR(14,4)-E
     -   ijR(15,3)+EijR(15,4)+EijR(16,4)-EijR(19,3)+EijR(20,3)-EijR(3
     -   2,4)-EijR(33,4)-EijR(34,4)+EijR(35,4)
       FI(267) = EijI(4,2)-EijI(4,4)-EijI(7,2)+EijI(14,3)+EijI(14,4)-E
     -   ijI(15,3)+EijI(15,4)+EijI(16,4)-EijI(19,3)+EijI(20,3)-EijI(3
     -   2,4)-EijI(33,4)-EijI(34,4)+EijI(35,4)
       F(267)=DCMPLX(FR(267),FI(267))
       FR(268) = -EijR(4,4)+EijR(14,3)+EijR(14,4)-EijR(15,3)+EijR(15,4
     -   )-EijR(32,4)
       FI(268) = -EijI(4,4)+EijI(14,3)+EijI(14,4)-EijI(15,3)+EijI(15,4
     -   )-EijI(32,4)
       F(268)=DCMPLX(FR(268),FI(268))
       FR(269) = -EijR(4,3)+EijR(4,4)+EijR(9,3)-EijR(10,3)-EijR(14,4)+
     -   EijR(16,3)-EijR(16,4)+EijR(21,4)+EijR(26,4)-EijR(27,4)-EijR(
     -   28,4)+EijR(33,4)+2*(EijR(15,3)-EijR(15,4)-EijR(20,3)+EijR(32
     -   ,4)+EijR(34,4)-EijR(35,4))
       FI(269) = -EijI(4,3)+EijI(4,4)+EijI(9,3)-EijI(10,3)-EijI(14,4)+
     -   EijI(16,3)-EijI(16,4)+EijI(21,4)+EijI(26,4)-EijI(27,4)-EijI(
     -   28,4)+EijI(33,4)+2*(EijI(15,3)-EijI(15,4)-EijI(20,3)+EijI(32
     -   ,4)+EijI(34,4)-EijI(35,4))
       F(269)=DCMPLX(FR(269),FI(269))
       FR(270) = -EijR(4,3)+EijR(4,4)+EijR(9,3)+EijR(9,4)-EijR(10,3)-E
     -   ijR(10,4)+EijR(16,3)-EijR(16,4)+2*(EijR(15,3)-EijR(20,3))-3*
     -   (EijR(15,4)-EijR(21,4)+EijR(28,4)-EijR(34,4))
       FI(270) = -EijI(4,3)+EijI(4,4)+EijI(9,3)+EijI(9,4)-EijI(10,3)-E
     -   ijI(10,4)+EijI(16,3)-EijI(16,4)+2*(EijI(15,3)-EijI(20,3))-3*
     -   (EijI(15,4)-EijI(21,4)+EijI(28,4)-EijI(34,4))
       F(270)=DCMPLX(FR(270),FI(270))
       FR(271) = -EijR(4,3)+EijR(4,4)-EijR(10,3)-EijR(10,4)+2*EijR(15,
     -   3)-3*(EijR(15,4)-EijR(21,4))
       FI(271) = -EijI(4,3)+EijI(4,4)-EijI(10,3)-EijI(10,4)+2*EijI(15,
     -   3)-3*(EijI(15,4)-EijI(21,4))
       F(271)=DCMPLX(FR(271),FI(271))
       FR(272) = EijR(2,2)+EijR(4,4)+EijR(8,2)+EijR(8,3)+EijR(9,3)-Eij
     -   R(10,2)-EijR(14,4)+5*EijR(15,3)-EijR(16,4)+EijR(17,3)-EijR(1
     -   9,3)-3*(EijR(4,3)+EijR(9,2)+EijR(18,3)+EijR(20,3))+EijR(21,4
     -   )+EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(33,4)+2*(EijR(4,2)-E
     -   ijR(10,3)+EijR(14,3)-EijR(15,4)+EijR(16,3)+EijR(32,4)+EijR(3
     -   4,4)-EijR(35,4))
       FI(272) = EijI(2,2)+EijI(4,4)+EijI(8,2)+EijI(8,3)+EijI(9,3)-Eij
     -   I(10,2)-EijI(14,4)+5*EijI(15,3)-EijI(16,4)+EijI(17,3)-EijI(1
     -   9,3)-3*(EijI(4,3)+EijI(9,2)+EijI(18,3)+EijI(20,3))+EijI(21,4
     -   )+EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(33,4)+2*(EijI(4,2)-E
     -   ijI(10,3)+EijI(14,3)-EijI(15,4)+EijI(16,3)+EijI(32,4)+EijI(3
     -   4,4)-EijI(35,4))
       F(272)=DCMPLX(FR(272),FI(272))
       FR(273) = EijR(2,2)+EijR(2,3)+EijR(4,2)+EijR(4,4)+EijR(9,3)+Eij
     -   R(9,4)-4*EijR(10,3)-EijR(10,4)+5*EijR(15,3)+EijR(16,3)-EijR(
     -   16,4)-2*(EijR(4,3)+EijR(9,2)+EijR(20,3))-3*(EijR(15,4)-EijR(
     -   21,4)+EijR(28,4)-EijR(34,4))
       FI(273) = EijI(2,2)+EijI(2,3)+EijI(4,2)+EijI(4,4)+EijI(9,3)+Eij
     -   I(9,4)-4*EijI(10,3)-EijI(10,4)+5*EijI(15,3)+EijI(16,3)-EijI(
     -   16,4)-2*(EijI(4,3)+EijI(9,2)+EijI(20,3))-3*(EijI(15,4)-EijI(
     -   21,4)+EijI(28,4)-EijI(34,4))
       F(273)=DCMPLX(FR(273),FI(273))
       FR(274) = -EijR(4,3)+EijR(4,4)+EijR(9,3)-EijR(10,3)+EijR(16,3)+
     -   EijR(20,4)+EijR(21,4)+EijR(22,4)+2*(EijR(15,3)-EijR(15,4)-Ei
     -   jR(16,4)-EijR(20,3)-EijR(28,4)-EijR(31,4)+2*EijR(34,4))
       FI(274) = -EijI(4,3)+EijI(4,4)+EijI(9,3)-EijI(10,3)+EijI(16,3)+
     -   EijI(20,4)+EijI(21,4)+EijI(22,4)+2*(EijI(15,3)-EijI(15,4)-Ei
     -   jI(16,4)-EijI(20,3)-EijI(28,4)-EijI(31,4)+2*EijI(34,4))
       F(274)=DCMPLX(FR(274),FI(274))
       FR(275) = EijR(4,2)-EijR(4,3)+EijR(4,4)-EijR(9,2)-EijR(10,3)-Ei
     -   jR(16,4)+EijR(21,4)-EijR(28,4)+2*(EijR(15,3)-EijR(15,4)+EijR
     -   (34,4))
       FI(275) = EijI(4,2)-EijI(4,3)+EijI(4,4)-EijI(9,2)-EijI(10,3)-Ei
     -   jI(16,4)+EijI(21,4)-EijI(28,4)+2*(EijI(15,3)-EijI(15,4)+EijI
     -   (34,4))
       F(275)=DCMPLX(FR(275),FI(275))
       FR(276) = EijR(4,2)+EijR(4,4)-EijR(9,2)-EijR(10,3)+EijR(14,3)-E
     -   ijR(14,4)+3*EijR(15,3)-EijR(18,3)+EijR(21,4)-EijR(27,4)-2*(E
     -   ijR(4,3)+EijR(15,4)-EijR(32,4))
       FI(276) = EijI(4,2)+EijI(4,4)-EijI(9,2)-EijI(10,3)+EijI(14,3)-E
     -   ijI(14,4)+3*EijI(15,3)-EijI(18,3)+EijI(21,4)-EijI(27,4)-2*(E
     -   ijI(4,3)+EijI(15,4)-EijI(32,4))
       F(276)=DCMPLX(FR(276),FI(276))
       FR(277) = EijR(4,2)-EijR(4,4)-EijR(9,2)+EijR(16,4)-EijR(21,4)+E
     -   ijR(28,4)+2*(EijR(15,4)-EijR(34,4))
       FI(277) = EijI(4,2)-EijI(4,4)-EijI(9,2)+EijI(16,4)-EijI(21,4)+E
     -   ijI(28,4)+2*(EijI(15,4)-EijI(34,4))
       F(277)=DCMPLX(FR(277),FI(277))
       FR(278) = EijR(2,2)+EijR(4,2)+EijR(4,4)-EijR(5,2)+EijR(7,2)+Eij
     -   R(8,2)+EijR(9,3)-EijR(10,2)-EijR(10,3)+EijR(14,3)-EijR(14,4)
     -   -EijR(16,4)+EijR(17,3)-EijR(18,3)-EijR(19,3)+3*(EijR(15,3)-E
     -   ijR(20,3))+EijR(21,4)+EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(
     -   33,4)-2*(EijR(4,3)+EijR(9,2)+EijR(15,4)-EijR(16,3)-EijR(32,4
     -   )-EijR(34,4)+EijR(35,4))
       FI(278) = EijI(2,2)+EijI(4,2)+EijI(4,4)-EijI(5,2)+EijI(7,2)+Eij
     -   I(8,2)+EijI(9,3)-EijI(10,2)-EijI(10,3)+EijI(14,3)-EijI(14,4)
     -   -EijI(16,4)+EijI(17,3)-EijI(18,3)-EijI(19,3)+3*(EijI(15,3)-E
     -   ijI(20,3))+EijI(21,4)+EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(
     -   33,4)-2*(EijI(4,3)+EijI(9,2)+EijI(15,4)-EijI(16,3)-EijI(32,4
     -   )-EijI(34,4)+EijI(35,4))
       F(278)=DCMPLX(FR(278),FI(278))
       FR(279) = EijR(4,3)-EijR(4,4)+EijR(6,2)-EijR(7,2)-EijR(8,2)+Eij
     -   R(9,2)-EijR(12,3)+EijR(13,3)+EijR(14,4)-EijR(15,3)+EijR(15,4
     -   )-EijR(22,4)-EijR(29,4)+EijR(30,4)+EijR(31,4)-EijR(32,4)-2*(
     -   EijR(16,3)-EijR(16,4)-EijR(20,3)+EijR(33,4)+EijR(34,4)-EijR(
     -   35,4))
       FI(279) = EijI(4,3)-EijI(4,4)+EijI(6,2)-EijI(7,2)-EijI(8,2)+Eij
     -   I(9,2)-EijI(12,3)+EijI(13,3)+EijI(14,4)-EijI(15,3)+EijI(15,4
     -   )-EijI(22,4)-EijI(29,4)+EijI(30,4)+EijI(31,4)-EijI(32,4)-2*(
     -   EijI(16,3)-EijI(16,4)-EijI(20,3)+EijI(33,4)+EijI(34,4)-EijI(
     -   35,4))
       F(279)=DCMPLX(FR(279),FI(279))
       FR(280) = -EijR(4,3)+EijR(4,4)+EijR(12,3)-EijR(13,3)+EijR(15,3)
     -   +EijR(20,4)+EijR(21,4)+EijR(22,4)-2*(EijR(15,4)-EijR(16,3)+E
     -   ijR(16,4)+EijR(20,3)+EijR(28,4)+EijR(31,4)-2*EijR(34,4))
       FI(280) = -EijI(4,3)+EijI(4,4)+EijI(12,3)-EijI(13,3)+EijI(15,3)
     -   +EijI(20,4)+EijI(21,4)+EijI(22,4)-2*(EijI(15,4)-EijI(16,3)+E
     -   ijI(16,4)+EijI(20,3)+EijI(28,4)+EijI(31,4)-2*EijI(34,4))
       F(280)=DCMPLX(FR(280),FI(280))
       FR(281) = EijR(3,1)-EijR(4,1)-EijR(4,3)+EijR(4,4)+EijR(8,2)-Eij
     -   R(10,2)+EijR(15,3)+EijR(16,3)-EijR(16,4)-EijR(20,3)+EijR(21,
     -   4)-EijR(28,4)-2*(EijR(15,4)-EijR(34,4))
       FI(281) = EijI(3,1)-EijI(4,1)-EijI(4,3)+EijI(4,4)+EijI(8,2)-Eij
     -   I(10,2)+EijI(15,3)+EijI(16,3)-EijI(16,4)-EijI(20,3)+EijI(21,
     -   4)-EijI(28,4)-2*(EijI(15,4)-EijI(34,4))
       F(281)=DCMPLX(FR(281),FI(281))
       FR(282) = EijR(4,4)+EijR(6,2)-EijR(7,2)+EijR(8,2)-EijR(9,2)+Eij
     -   R(11,3)+EijR(12,3)-EijR(14,4)-EijR(15,4)+5*EijR(16,3)+EijR(1
     -   7,3)-EijR(18,3)-3*(EijR(4,3)+EijR(19,3)+EijR(20,3))+EijR(22,
     -   4)+EijR(29,4)-EijR(30,4)-EijR(31,4)+EijR(32,4)+2*(EijR(4,2)-
     -   EijR(10,2)-EijR(13,3)+EijR(14,3)+EijR(15,3)-EijR(16,4)+EijR(
     -   33,4)+EijR(34,4)-EijR(35,4))
       FI(282) = EijI(4,4)+EijI(6,2)-EijI(7,2)+EijI(8,2)-EijI(9,2)+Eij
     -   I(11,3)+EijI(12,3)-EijI(14,4)-EijI(15,4)+5*EijI(16,3)+EijI(1
     -   7,3)-EijI(18,3)-3*(EijI(4,3)+EijI(19,3)+EijI(20,3))+EijI(22,
     -   4)+EijI(29,4)-EijI(30,4)-EijI(31,4)+EijI(32,4)+2*(EijI(4,2)-
     -   EijI(10,2)-EijI(13,3)+EijI(14,3)+EijI(15,3)-EijI(16,4)+EijI(
     -   33,4)+EijI(34,4)-EijI(35,4))
       F(282)=DCMPLX(FR(282),FI(282))
       FR(283) = EijR(4,2)+EijR(4,4)+EijR(8,2)-EijR(9,2)+EijR(9,3)-Eij
     -   R(10,2)-EijR(10,3)+EijR(12,3)-EijR(13,3)+3*(EijR(15,3)+EijR(
     -   16,3))-4*EijR(20,3)+EijR(20,4)+EijR(21,4)+EijR(22,4)-2*(EijR
     -   (4,3)+EijR(15,4)+EijR(16,4)+EijR(28,4)+EijR(31,4)-2*EijR(34,
     -   4))
       FI(283) = EijI(4,2)+EijI(4,4)+EijI(8,2)-EijI(9,2)+EijI(9,3)-Eij
     -   I(10,2)-EijI(10,3)+EijI(12,3)-EijI(13,3)+3*(EijI(15,3)+EijI(
     -   16,3))-4*EijI(20,3)+EijI(20,4)+EijI(21,4)+EijI(22,4)-2*(EijI
     -   (4,3)+EijI(15,4)+EijI(16,4)+EijI(28,4)+EijI(31,4)-2*EijI(34,
     -   4))
       F(283)=DCMPLX(FR(283),FI(283))
       FR(284) = -EijR(4,3)+EijR(4,4)+EijR(12,3)+EijR(12,4)-EijR(13,3)
     -   -EijR(13,4)+EijR(15,3)-EijR(15,4)+2*(EijR(16,3)-EijR(20,3))-
     -   3*(EijR(16,4)-EijR(22,4)+EijR(31,4)-EijR(34,4))
       FI(284) = -EijI(4,3)+EijI(4,4)+EijI(12,3)+EijI(12,4)-EijI(13,3)
     -   -EijI(13,4)+EijI(15,3)-EijI(15,4)+2*(EijI(16,3)-EijI(20,3))-
     -   3*(EijI(16,4)-EijI(22,4)+EijI(31,4)-EijI(34,4))
       F(284)=DCMPLX(FR(284),FI(284))
       FR(285) = EijR(3,1)+EijR(3,2)-EijR(4,1)+EijR(4,2)-EijR(4,3)+Eij
     -   R(4,4)+EijR(15,3)-EijR(15,4)+EijR(16,3)-EijR(20,3)+EijR(22,4
     -   )-EijR(31,4)-2*(EijR(10,2)+EijR(16,4)-EijR(34,4))
       FI(285) = EijI(3,1)+EijI(3,2)-EijI(4,1)+EijI(4,2)-EijI(4,3)+Eij
     -   I(4,4)+EijI(15,3)-EijI(15,4)+EijI(16,3)-EijI(20,3)+EijI(22,4
     -   )-EijI(31,4)-2*(EijI(10,2)+EijI(16,4)-EijI(34,4))
       F(285)=DCMPLX(FR(285),FI(285))
       FR(286) = EijR(2,1)+EijR(2,2)-EijR(4,1)+EijR(4,3)-EijR(4,4)-Eij
     -   R(9,2)-EijR(15,3)-EijR(16,3)+EijR(16,4)+EijR(20,3)-EijR(21,4
     -   )+EijR(28,4)+2*(EijR(15,4)-EijR(34,4))
       FI(286) = EijI(2,1)+EijI(2,2)-EijI(4,1)+EijI(4,3)-EijI(4,4)-Eij
     -   I(9,2)-EijI(15,3)-EijI(16,3)+EijI(16,4)+EijI(20,3)-EijI(21,4
     -   )+EijI(28,4)+2*(EijI(15,4)-EijI(34,4))
       F(286)=DCMPLX(FR(286),FI(286))
       FR(287) = EijR(3,1)-EijR(4,1)+EijR(4,2)-EijR(4,4)+EijR(8,2)-Eij
     -   R(9,2)-EijR(10,2)+EijR(15,4)-EijR(22,4)+EijR(31,4)+2*(EijR(1
     -   6,4)-EijR(34,4))
       FI(287) = EijI(3,1)-EijI(4,1)+EijI(4,2)-EijI(4,4)+EijI(8,2)-Eij
     -   I(9,2)-EijI(10,2)+EijI(15,4)-EijI(22,4)+EijI(31,4)+2*(EijI(1
     -   6,4)-EijI(34,4))
       F(287)=DCMPLX(FR(287),FI(287))
       FR(288) = -EijR(4,4)-EijR(9,2)+EijR(10,2)+EijR(15,4)+EijR(16,4)
     -   -EijR(34,4)
       FI(288) = -EijI(4,4)-EijI(9,2)+EijI(10,2)+EijI(15,4)+EijI(16,4)
     -   -EijI(34,4)
       F(288)=DCMPLX(FR(288),FI(288))
       FR(289) = EijR(2,1)+EijR(2,2)-EijR(4,1)+EijR(4,2)+EijR(4,3)-Eij
     -   R(4,4)+EijR(5,2)-EijR(7,2)+EijR(8,3)+EijR(14,4)-EijR(15,3)-E
     -   ijR(18,3)-EijR(21,4)+EijR(27,4)-2*(EijR(9,2)-EijR(15,4)+EijR
     -   (32,4))
       FI(289) = EijI(2,1)+EijI(2,2)-EijI(4,1)+EijI(4,2)+EijI(4,3)-Eij
     -   I(4,4)+EijI(5,2)-EijI(7,2)+EijI(8,3)+EijI(14,4)-EijI(15,3)-E
     -   ijI(18,3)-EijI(21,4)+EijI(27,4)-2*(EijI(9,2)-EijI(15,4)+EijI
     -   (32,4))
       F(289)=DCMPLX(FR(289),FI(289))
       FR(290) = EijR(3,1)-EijR(4,1)+EijR(4,2)-EijR(4,4)+EijR(6,2)+Eij
     -   R(8,2)-EijR(9,2)-2*EijR(10,2)+EijR(14,3)+EijR(14,4)+EijR(15,
     -   4)+EijR(16,4)+EijR(17,3)-EijR(18,3)-EijR(19,3)-EijR(32,4)-Ei
     -   jR(33,4)-EijR(34,4)+EijR(35,4)
       FI(290) = EijI(3,1)-EijI(4,1)+EijI(4,2)-EijI(4,4)+EijI(6,2)+Eij
     -   I(8,2)-EijI(9,2)-2*EijI(10,2)+EijI(14,3)+EijI(14,4)+EijI(15,
     -   4)+EijI(16,4)+EijI(17,3)-EijI(18,3)-EijI(19,3)-EijI(32,4)-Ei
     -   jI(33,4)-EijI(34,4)+EijI(35,4)
       F(290)=DCMPLX(FR(290),FI(290))
       FR(291) = EijR(4,2)-EijR(4,4)-EijR(9,2)+EijR(14,3)+EijR(14,4)+E
     -   ijR(15,4)-EijR(18,3)-EijR(32,4)
       FI(291) = EijI(4,2)-EijI(4,4)-EijI(9,2)+EijI(14,3)+EijI(14,4)+E
     -   ijI(15,4)-EijI(18,3)-EijI(32,4)
       F(291)=DCMPLX(FR(291),FI(291))
       FR(292) = EijR(4,2)-EijR(4,4)-EijR(9,2)-EijR(10,3)+EijR(15,3)+2
     -   *EijR(15,4)-EijR(21,4)
       FI(292) = EijI(4,2)-EijI(4,4)-EijI(9,2)-EijI(10,3)+EijI(15,3)+2
     -   *EijI(15,4)-EijI(21,4)
       F(292)=DCMPLX(FR(292),FI(292))
       FR(293) = EijR(4,2)+2*EijR(4,3)-EijR(4,4)-EijR(7,2)+EijR(8,2)-E
     -   ijR(10,2)-EijR(14,3)+EijR(14,4)-EijR(15,3)+EijR(15,4)-EijR(1
     -   6,3)+EijR(16,4)+EijR(17,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+
     -   EijR(35,4)
       FI(293) = EijI(4,2)+2*EijI(4,3)-EijI(4,4)-EijI(7,2)+EijI(8,2)-E
     -   ijI(10,2)-EijI(14,3)+EijI(14,4)-EijI(15,3)+EijI(15,4)-EijI(1
     -   6,3)+EijI(16,4)+EijI(17,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+
     -   EijI(35,4)
       F(293)=DCMPLX(FR(293),FI(293))
       FR(294) = EijR(4,2)+EijR(4,3)-EijR(4,4)+EijR(8,2)-EijR(9,2)+Eij
     -   R(9,3)-EijR(10,2)-EijR(15,3)+EijR(16,4)-EijR(20,3)-EijR(21,4
     -   )+EijR(28,4)+2*(EijR(15,4)-EijR(34,4))
       FI(294) = EijI(4,2)+EijI(4,3)-EijI(4,4)+EijI(8,2)-EijI(9,2)+Eij
     -   I(9,3)-EijI(10,2)-EijI(15,3)+EijI(16,4)-EijI(20,3)-EijI(21,4
     -   )+EijI(28,4)+2*(EijI(15,4)-EijI(34,4))
       F(294)=DCMPLX(FR(294),FI(294))
       FR(295) = EijR(3,2)+EijR(4,2)-EijR(4,4)+EijR(12,3)-EijR(13,3)+E
     -   ijR(15,4)+EijR(16,3)-EijR(20,3)-EijR(22,4)+EijR(31,4)-2*(Eij
     -   R(10,2)-EijR(16,4)+EijR(34,4))
       FI(295) = EijI(3,2)+EijI(4,2)-EijI(4,4)+EijI(12,3)-EijI(13,3)+E
     -   ijI(15,4)+EijI(16,3)-EijI(20,3)-EijI(22,4)+EijI(31,4)-2*(Eij
     -   I(10,2)-EijI(16,4)+EijI(34,4))
       F(295)=DCMPLX(FR(295),FI(295))
       FR(296) = EijR(4,4)-EijR(15,4)-EijR(16,3)-EijR(16,4)+EijR(20,3)
     -   +EijR(34,4)
       FI(296) = EijI(4,4)-EijI(15,4)-EijI(16,3)-EijI(16,4)+EijI(20,3)
     -   +EijI(34,4)
       F(296)=DCMPLX(FR(296),FI(296))
       FR(297) = EijR(4,2)+EijR(4,3)-EijR(4,4)-EijR(9,2)+EijR(14,4)+Ei
     -   jR(15,4)-EijR(18,3)-EijR(32,4)
       FI(297) = EijI(4,2)+EijI(4,3)-EijI(4,4)-EijI(9,2)+EijI(14,4)+Ei
     -   jI(15,4)-EijI(18,3)-EijI(32,4)
       F(297)=DCMPLX(FR(297),FI(297))
       FR(298) = EijR(4,2)-EijR(4,3)-EijR(4,4)-EijR(10,2)+EijR(15,3)+E
     -   ijR(15,4)+EijR(16,3)+EijR(16,4)-EijR(20,3)-EijR(34,4)
       FI(298) = EijI(4,2)-EijI(4,3)-EijI(4,4)-EijI(10,2)+EijI(15,3)+E
     -   ijI(15,4)+EijI(16,3)+EijI(16,4)-EijI(20,3)-EijI(34,4)
       F(298)=DCMPLX(FR(298),FI(298))
       FR(299) = -EijR(4,3)-EijR(4,4)+EijR(15,3)+EijR(15,4)
       FI(299) = -EijI(4,3)-EijI(4,4)+EijI(15,3)+EijI(15,4)
       F(299)=DCMPLX(FR(299),FI(299))
       P(380) = p4sq+s12
       P(381) = p2sq-s23+s34+s45-2*P(380)
       P(382) = p2sq-p4sq+s34+s45
       P(383) = p3sq+3*s12+s23-2*P(382)
       P(384) = -2*p3sq+s12+s23+s34+s45
       FR(300) = 8*(Dij2345R(1,1)-Dij2345R(3,1)+p3sq*(EijR(3,1)+2*EijR
     -   (6,2))+2*(2*EijR(11,2)-EijR(22,3)+EijR(24,3)-EijR(39,4)-EijR
     -   (40,4)+EijR(42,4)+EijR(44,4))-EijR(2,1)*P(2)-EijR(5,2)*P(34)
     -   -EijR(10,2)*P(38)+EijR(4,1)*P(43)-EijR(8,2)*P(128)+EijR(2,2)
     -   *P(152)-EijR(4,2)*P(381)-EijR(7,2)*P(383)-EijR(9,2)*P(384))
       FI(300) = 8*(Dij2345I(1,1)-Dij2345I(3,1)+p3sq*(EijI(3,1)+2*EijI
     -   (6,2))+2*(2*EijI(11,2)-EijI(22,3)+EijI(24,3)-EijI(39,4)-EijI
     -   (40,4)+EijI(42,4)+EijI(44,4))-EijI(2,1)*P(2)-EijI(5,2)*P(34)
     -   -EijI(10,2)*P(38)+EijI(4,1)*P(43)-EijI(8,2)*P(128)+EijI(2,2)
     -   *P(152)-EijI(4,2)*P(381)-EijI(7,2)*P(383)-EijI(9,2)*P(384))
       F(300)=DCMPLX(FR(300),FI(300))
       FR(301) = -8*(2*(EijR(22,3)-EijR(24,3)+EijR(37,4)+EijR(39,4)-2*
     -   EijR(44,4))+(EijR(2,1)-EijR(4,1)+EijR(5,2)-EijR(7,2))*P(14)+
     -   (EijR(4,2)-EijR(9,2))*P(169))
       FI(301) = -8*(2*(EijI(22,3)-EijI(24,3)+EijI(37,4)+EijI(39,4)-2*
     -   EijI(44,4))+(EijI(2,1)-EijI(4,1)+EijI(5,2)-EijI(7,2))*P(14)+
     -   (EijI(4,2)-EijI(9,2))*P(169))
       F(301)=DCMPLX(FR(301),FI(301))
       FR(302) = 8*(-Dij2345R(1,1)+Dij2345R(3,1)-s12*EijR(2,1)+p3sq*Ei
     -   jR(3,1)-2*(EijR(11,2)-EijR(22,3)+EijR(39,4)-EijR(44,4))-EijR
     -   (4,1)*P(3)+EijR(5,2)*P(9)+EijR(10,2)*P(38)+EijR(8,2)*P(128)-
     -   EijR(2,2)*P(152)+EijR(7,2)*P(180)-EijR(4,2)*P(296)+EijR(9,2)
     -   *P(384))
       FI(302) = 8*(-Dij2345I(1,1)+Dij2345I(3,1)-s12*EijI(2,1)+p3sq*Ei
     -   jI(3,1)-2*(EijI(11,2)-EijI(22,3)+EijI(39,4)-EijI(44,4))-EijI
     -   (4,1)*P(3)+EijI(5,2)*P(9)+EijI(10,2)*P(38)+EijI(8,2)*P(128)-
     -   EijI(2,2)*P(152)+EijI(7,2)*P(180)-EijI(4,2)*P(296)+EijI(9,2)
     -   *P(384))
       F(302)=DCMPLX(FR(302),FI(302))
       P(385) = p2sq+2*p3sq+s12-s23-s45
       P(386) = p2sq+2*p3sq-s23-3*P(7)
       FR(303) = 8*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)+s15*(EE0R+EijR
     -   (1,1))-2*(p3sq*EijR(3,2)+5*EijR(11,2)+EijR(23,3)-EijR(24,3)+
     -   EijR(39,4)+EijR(43,4)-EijR(44,4)-EijR(45,4))+EijR(5,2)*P(9)+
     -   (-EijR(3,1)+EijR(4,1)-2*EijR(6,2))*P(14)-EijR(2,2)*P(152)+Ei
     -   jR(7,2)*P(180)-EijR(4,2)*P(296)+EijR(9,2)*P(384)+EijR(8,2)*P
     -   (385)+EijR(10,2)*P(386))
       FI(303) = 8*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)+s15*(EE0I+EijI
     -   (1,1))-2*(p3sq*EijI(3,2)+5*EijI(11,2)+EijI(23,3)-EijI(24,3)+
     -   EijI(39,4)+EijI(43,4)-EijI(44,4)-EijI(45,4))+EijI(5,2)*P(9)+
     -   (-EijI(3,1)+EijI(4,1)-2*EijI(6,2))*P(14)-EijI(2,2)*P(152)+Ei
     -   jI(7,2)*P(180)-EijI(4,2)*P(296)+EijI(9,2)*P(384)+EijI(8,2)*P
     -   (385)+EijI(10,2)*P(386))
       F(303)=DCMPLX(FR(303),FI(303))
       P(387) = p2sq-s12-s23+s45-2*P(91)
       FR(304) = -8*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)-8*EijR(11,2)-
     -   4*EijR(24,3)+2*(EijR(22,3)+EijR(23,3)+EijR(39,4)+EijR(43,4)-
     -   EijR(44,4)-EijR(45,4))+EijR(5,2)*P(9)+(-EijR(3,1)+EijR(4,1))
     -   *P(130)+EijR(3,2)*P(131)-EijR(6,2)*P(134)+EijR(7,2)*P(135)-E
     -   ijR(2,2)*P(152)-EijR(4,2)*P(213)+EijR(9,2)*P(384)+EijR(8,2)*
     -   P(385)-EijR(10,2)*P(387))
       FI(304) = -8*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)-8*EijI(11,2)-
     -   4*EijI(24,3)+2*(EijI(22,3)+EijI(23,3)+EijI(39,4)+EijI(43,4)-
     -   EijI(44,4)-EijI(45,4))+EijI(5,2)*P(9)+(-EijI(3,1)+EijI(4,1))
     -   *P(130)+EijI(3,2)*P(131)-EijI(6,2)*P(134)+EijI(7,2)*P(135)-E
     -   ijI(2,2)*P(152)-EijI(4,2)*P(213)+EijI(9,2)*P(384)+EijI(8,2)*
     -   P(385)-EijI(10,2)*P(387))
       F(304)=DCMPLX(FR(304),FI(304))
       P(388) = p4sq+s12+s15-s34-2*s45
       P(389) = p2sq-p3sq-p4sq-s12+s15+s34
       FR(305) = 8*(-Dij2345R(1,1)+Dij2345R(3,1)+s15*(EE0R+2*EijR(1,1)
     -   +EijR(1,2))-2*(3*EijR(11,2)+EijR(21,3)+EijR(22,3)-2*EijR(24,
     -   3)+EijR(39,4)+EijR(40,4)-EijR(42,4)-EijR(44,4))+EijR(2,1)*P(
     -   2)+EijR(5,2)*P(34)+EijR(8,2)*P(128)-EijR(2,2)*P(152)-EijR(3,
     -   1)*P(160)-EijR(6,2)*P(163)+EijR(10,2)*P(165)-EijR(4,2)*P(213
     -   )+EijR(9,2)*P(384)-EijR(4,1)*P(388)-EijR(7,2)*P(389))
       FI(305) = 8*(-Dij2345I(1,1)+Dij2345I(3,1)+s15*(EE0I+2*EijI(1,1)
     -   +EijI(1,2))-2*(3*EijI(11,2)+EijI(21,3)+EijI(22,3)-2*EijI(24,
     -   3)+EijI(39,4)+EijI(40,4)-EijI(42,4)-EijI(44,4))+EijI(2,1)*P(
     -   2)+EijI(5,2)*P(34)+EijI(8,2)*P(128)-EijI(2,2)*P(152)-EijI(3,
     -   1)*P(160)-EijI(6,2)*P(163)+EijI(10,2)*P(165)-EijI(4,2)*P(213
     -   )+EijI(9,2)*P(384)-EijI(4,1)*P(388)-EijI(7,2)*P(389))
       F(305)=DCMPLX(FR(305),FI(305))
       P(390) = p2sq-s12+s15-s23+s45
       FR(306) = 8*(s15*(EE0R+EijR(1,1))-2*(EijR(11,2)+EijR(37,4)+EijR
     -   (39,4)+2*(EijR(22,3)-EijR(24,3)-EijR(44,4)))-EijR(4,1)*P(130
     -   )-(EijR(3,1)+EijR(8,2)-EijR(10,2))*P(169)+(EijR(2,1)+EijR(5,
     -   2)-EijR(7,2))*P(390))
       FI(306) = 8*(s15*(EE0I+EijI(1,1))-2*(EijI(11,2)+EijI(37,4)+EijI
     -   (39,4)+2*(EijI(22,3)-EijI(24,3)-EijI(44,4)))-EijI(4,1)*P(130
     -   )-(EijI(3,1)+EijI(8,2)-EijI(10,2))*P(169)+(EijI(2,1)+EijI(5,
     -   2)-EijI(7,2))*P(390))
       F(306)=DCMPLX(FR(306),FI(306))
       P(391) = p3sq-s12-s15
       FR(307) = 8*(Dij2345R(1,1)-Dij2345R(3,1)+s12*EijR(2,1)-p3sq*Eij
     -   R(3,1)+2*(2*EijR(11,2)+EijR(24,3)-EijR(39,4)+EijR(44,4))-Eij
     -   R(5,2)*P(9)-EijR(8,2)*P(128)-EijR(7,2)*P(135)+EijR(2,2)*P(15
     -   2)+EijR(10,2)*P(190)+EijR(4,2)*P(213)-EijR(9,2)*P(384)+EijR(
     -   4,1)*P(391))
       FI(307) = 8*(Dij2345I(1,1)-Dij2345I(3,1)+s12*EijI(2,1)-p3sq*Eij
     -   I(3,1)+2*(2*EijI(11,2)+EijI(24,3)-EijI(39,4)+EijI(44,4))-Eij
     -   I(5,2)*P(9)-EijI(8,2)*P(128)-EijI(7,2)*P(135)+EijI(2,2)*P(15
     -   2)+EijI(10,2)*P(190)+EijI(4,2)*P(213)-EijI(9,2)*P(384)+EijI(
     -   4,1)*P(391))
       F(307)=DCMPLX(FR(307),FI(307))
       FR(308) = -16*(EijR(3,1)-EijR(4,1)+EijR(4,4)+EijR(6,4)-EijR(7,4
     -   )-EijR(16,4)+6*(EijR(14,3)-EijR(19,3))+3*(EijR(4,2)-EijR(4,3
     -   )+EijR(6,2)+EijR(6,3)-EijR(7,2)-EijR(7,3)-EijR(10,2)-EijR(14
     -   ,4)+EijR(16,3)+EijR(19,4)-EijR(25,4)+EijR(33,4)))
       FI(308) = -16*(EijI(3,1)-EijI(4,1)+EijI(4,4)+EijI(6,4)-EijI(7,4
     -   )-EijI(16,4)+6*(EijI(14,3)-EijI(19,3))+3*(EijI(4,2)-EijI(4,3
     -   )+EijI(6,2)+EijI(6,3)-EijI(7,2)-EijI(7,3)-EijI(10,2)-EijI(14
     -   ,4)+EijI(16,3)+EijI(19,4)-EijI(25,4)+EijI(33,4)))
       F(308)=DCMPLX(FR(308),FI(308))
       FR(309) = -16*(EijR(4,4)+EijR(6,2)+EijR(6,3)-EijR(7,2)-EijR(7,3
     -   )+EijR(8,2)-EijR(9,2)-EijR(15,4)-3*(EijR(4,3)-EijR(16,3))-Ei
     -   jR(16,4)+4*(EijR(14,3)-EijR(19,3))+EijR(19,4)+EijR(23,4)-Eij
     -   R(24,4)-EijR(25,4)+EijR(34,4)+2*(EijR(4,2)-EijR(10,2)-EijR(1
     -   4,4)+EijR(15,3)+EijR(17,3)-EijR(18,3)-EijR(20,3)+EijR(32,4)+
     -   EijR(33,4)-EijR(35,4)))
       FI(309) = -16*(EijI(4,4)+EijI(6,2)+EijI(6,3)-EijI(7,2)-EijI(7,3
     -   )+EijI(8,2)-EijI(9,2)-EijI(15,4)-3*(EijI(4,3)-EijI(16,3))-Ei
     -   jI(16,4)+4*(EijI(14,3)-EijI(19,3))+EijI(19,4)+EijI(23,4)-Eij
     -   I(24,4)-EijI(25,4)+EijI(34,4)+2*(EijI(4,2)-EijI(10,2)-EijI(1
     -   4,4)+EijI(15,3)+EijI(17,3)-EijI(18,3)-EijI(20,3)+EijI(32,4)+
     -   EijI(33,4)-EijI(35,4)))
       F(309)=DCMPLX(FR(309),FI(309))
       FR(310) = -16*(EijR(3,2)+EijR(4,2)+EijR(4,4)+EijR(18,4)+4*(EijR
     -   (16,3)-EijR(19,3))+EijR(19,4)+EijR(22,4)-2*(EijR(4,3)+EijR(1
     -   0,2)-EijR(11,3)+EijR(13,3)-EijR(14,3)+EijR(14,4)+EijR(16,4)+
     -   EijR(25,4)+EijR(30,4)-2*EijR(33,4)))
       FI(310) = -16*(EijI(3,2)+EijI(4,2)+EijI(4,4)+EijI(18,4)+4*(EijI
     -   (16,3)-EijI(19,3))+EijI(19,4)+EijI(22,4)-2*(EijI(4,3)+EijI(1
     -   0,2)-EijI(11,3)+EijI(13,3)-EijI(14,3)+EijI(14,4)+EijI(16,4)+
     -   EijI(25,4)+EijI(30,4)-2*EijI(33,4)))
       F(310)=DCMPLX(FR(310),FI(310))
       FR(311) = -16*(EijR(4,2)+EijR(4,4)-EijR(10,2)-EijR(16,4)+EijR(1
     -   9,4)-EijR(25,4)-2*(EijR(4,3)-EijR(14,3)+EijR(14,4)-EijR(16,3
     -   )+EijR(19,3)-EijR(33,4)))
       FI(311) = -16*(EijI(4,2)+EijI(4,4)-EijI(10,2)-EijI(16,4)+EijI(1
     -   9,4)-EijI(25,4)-2*(EijI(4,3)-EijI(14,3)+EijI(14,4)-EijI(16,3
     -   )+EijI(19,3)-EijI(33,4)))
       F(311)=DCMPLX(FR(311),FI(311))
       FR(312) = -16*(EijR(3,1)-EijR(4,1)+EijR(4,4)-EijR(5,2)-EijR(5,3
     -   )+EijR(6,3)-EijR(7,2)+EijR(8,2)+EijR(15,3)-EijR(15,4)-3*(Eij
     -   R(10,2)-EijR(16,3))-EijR(16,4)-4*EijR(19,3)+EijR(19,4)+EijR(
     -   23,4)-EijR(24,4)-EijR(25,4)+EijR(34,4)+2*(EijR(4,2)-EijR(4,3
     -   )+EijR(6,2)+EijR(14,3)-EijR(14,4)+EijR(17,3)-EijR(20,3)+EijR
     -   (32,4)+EijR(33,4)-EijR(35,4)))
       FI(312) = -16*(EijI(3,1)-EijI(4,1)+EijI(4,4)-EijI(5,2)-EijI(5,3
     -   )+EijI(6,3)-EijI(7,2)+EijI(8,2)+EijI(15,3)-EijI(15,4)-3*(Eij
     -   I(10,2)-EijI(16,3))-EijI(16,4)-4*EijI(19,3)+EijI(19,4)+EijI(
     -   23,4)-EijI(24,4)-EijI(25,4)+EijI(34,4)+2*(EijI(4,2)-EijI(4,3
     -   )+EijI(6,2)+EijI(14,3)-EijI(14,4)+EijI(17,3)-EijI(20,3)+EijI
     -   (32,4)+EijI(33,4)-EijI(35,4)))
       F(312)=DCMPLX(FR(312),FI(312))
       FR(313) = -16*(EijR(4,2)+EijR(4,4)-EijR(5,2)+EijR(6,2)+EijR(8,2
     -   )-EijR(8,3)+EijR(9,3)+EijR(14,3)-EijR(14,4)+3*EijR(16,3)-Eij
     -   R(16,4)-4*EijR(20,3)+EijR(21,4)+EijR(26,4)-EijR(27,4)-EijR(2
     -   8,4)+EijR(33,4)-2*(EijR(4,3)+EijR(10,2)-EijR(15,3)+EijR(15,4
     -   )-EijR(17,3)+EijR(19,3)-EijR(32,4)-EijR(34,4)+EijR(35,4)))
       FI(313) = -16*(EijI(4,2)+EijI(4,4)-EijI(5,2)+EijI(6,2)+EijI(8,2
     -   )-EijI(8,3)+EijI(9,3)+EijI(14,3)-EijI(14,4)+3*EijI(16,3)-Eij
     -   I(16,4)-4*EijI(20,3)+EijI(21,4)+EijI(26,4)-EijI(27,4)-EijI(2
     -   8,4)+EijI(33,4)-2*(EijI(4,3)+EijI(10,2)-EijI(15,3)+EijI(15,4
     -   )-EijI(17,3)+EijI(19,3)-EijI(32,4)-EijI(34,4)+EijI(35,4)))
       F(313)=DCMPLX(FR(313),FI(313))
       FR(314) = 16*(EijR(4,2)+EijR(4,3)-EijR(4,4)+EijR(10,2)+EijR(14,
     -   4)+EijR(15,4)-2*(EijR(7,2)+EijR(16,3))+EijR(16,4)-EijR(18,3)
     -   +EijR(19,3)+EijR(20,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR
     -   (35,4))
       FI(314) = 16*(EijI(4,2)+EijI(4,3)-EijI(4,4)+EijI(10,2)+EijI(14,
     -   4)+EijI(15,4)-2*(EijI(7,2)+EijI(16,3))+EijI(16,4)-EijI(18,3)
     -   +EijI(19,3)+EijI(20,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI
     -   (35,4))
       F(314)=DCMPLX(FR(314),FI(314))
       FR(315) = 16*(EijR(4,3)-EijR(4,4)-EijR(7,2)-EijR(7,3)+EijR(10,2
     -   )+EijR(16,4)-EijR(19,4)+EijR(25,4)+2*(EijR(14,4)-EijR(16,3)+
     -   EijR(19,3)-EijR(33,4)))
       FI(315) = 16*(EijI(4,3)-EijI(4,4)-EijI(7,2)-EijI(7,3)+EijI(10,2
     -   )+EijI(16,4)-EijI(19,4)+EijI(25,4)+2*(EijI(14,4)-EijI(16,3)+
     -   EijI(19,3)-EijI(33,4)))
       F(315)=DCMPLX(FR(315),FI(315))
       FR(316) = -16*(EijR(3,1)-EijR(4,1)+EijR(4,4)+EijR(5,2)+EijR(5,3
     -   )+EijR(6,2)-EijR(7,3)+EijR(8,2)+3*(EijR(4,2)-EijR(4,3)+EijR(
     -   15,3))-EijR(15,4)-EijR(16,4)+4*(EijR(14,3)-EijR(18,3))+EijR(
     -   19,4)+EijR(23,4)-EijR(24,4)-EijR(25,4)+EijR(34,4)-2*(EijR(7,
     -   2)+EijR(9,2)+EijR(10,2)+EijR(14,4)-EijR(16,3)-EijR(17,3)+Eij
     -   R(19,3)+EijR(20,3)-EijR(32,4)-EijR(33,4)+EijR(35,4)))
       FI(316) = -16*(EijI(3,1)-EijI(4,1)+EijI(4,4)+EijI(5,2)+EijI(5,3
     -   )+EijI(6,2)-EijI(7,3)+EijI(8,2)+3*(EijI(4,2)-EijI(4,3)+EijI(
     -   15,3))-EijI(15,4)-EijI(16,4)+4*(EijI(14,3)-EijI(18,3))+EijI(
     -   19,4)+EijI(23,4)-EijI(24,4)-EijI(25,4)+EijI(34,4)-2*(EijI(7,
     -   2)+EijI(9,2)+EijI(10,2)+EijI(14,4)-EijI(16,3)-EijI(17,3)+Eij
     -   I(19,3)+EijI(20,3)-EijI(32,4)-EijI(33,4)+EijI(35,4)))
       F(316)=DCMPLX(FR(316),FI(316))
       FR(317) = 16*(EijR(3,1)-EijR(4,1)-EijR(4,2)-EijR(4,4)+EijR(6,2)
     -   +EijR(7,2)-EijR(10,2)-EijR(14,3)+EijR(14,4)+2*(EijR(4,3)-Eij
     -   R(15,3))+EijR(15,4)-EijR(16,3)+EijR(16,4)+EijR(18,3)+EijR(20
     -   ,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4))
       FI(317) = 16*(EijI(3,1)-EijI(4,1)-EijI(4,2)-EijI(4,4)+EijI(6,2)
     -   +EijI(7,2)-EijI(10,2)-EijI(14,3)+EijI(14,4)+2*(EijI(4,3)-Eij
     -   I(15,3))+EijI(15,4)-EijI(16,3)+EijI(16,4)+EijI(18,3)+EijI(20
     -   ,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4))
       F(317)=DCMPLX(FR(317),FI(317))
       FR(318) = -16*(EijR(3,1)-EijR(4,1)+EijR(4,4)+EijR(9,3)-EijR(10,
     -   3)+EijR(14,3)-EijR(14,4)-EijR(16,4)+EijR(17,3)-EijR(18,3)-Ei
     -   jR(19,3)+3*(EijR(15,3)-EijR(20,3))+EijR(21,4)+EijR(26,4)-Eij
     -   R(27,4)-EijR(28,4)+EijR(33,4)+2*(EijR(4,2)-EijR(4,3)+EijR(8,
     -   2)-EijR(9,2)-EijR(10,2)-EijR(15,4)+EijR(16,3)+EijR(32,4)+Eij
     -   R(34,4)-EijR(35,4)))
       FI(318) = -16*(EijI(3,1)-EijI(4,1)+EijI(4,4)+EijI(9,3)-EijI(10,
     -   3)+EijI(14,3)-EijI(14,4)-EijI(16,4)+EijI(17,3)-EijI(18,3)-Ei
     -   jI(19,3)+3*(EijI(15,3)-EijI(20,3))+EijI(21,4)+EijI(26,4)-Eij
     -   I(27,4)-EijI(28,4)+EijI(33,4)+2*(EijI(4,2)-EijI(4,3)+EijI(8,
     -   2)-EijI(9,2)-EijI(10,2)-EijI(15,4)+EijI(16,3)+EijI(32,4)+Eij
     -   I(34,4)-EijI(35,4)))
       F(318)=DCMPLX(FR(318),FI(318))
       FR(319) = -16*(EijR(2,1)-EijR(4,1)-EijR(4,3)+EijR(4,4)+EijR(5,2
     -   )+EijR(7,2)-EijR(9,2)-EijR(10,2)-EijR(14,4)+EijR(15,3)-EijR(
     -   15,4)+EijR(16,3)-EijR(16,4)-EijR(20,3)+EijR(32,4)+EijR(33,4)
     -   +EijR(34,4)-EijR(35,4))
       FI(319) = -16*(EijI(2,1)-EijI(4,1)-EijI(4,3)+EijI(4,4)+EijI(5,2
     -   )+EijI(7,2)-EijI(9,2)-EijI(10,2)-EijI(14,4)+EijI(15,3)-EijI(
     -   15,4)+EijI(16,3)-EijI(16,4)-EijI(20,3)+EijI(32,4)+EijI(33,4)
     -   +EijI(34,4)-EijI(35,4))
       F(319)=DCMPLX(FR(319),FI(319))
       FR(320) = -16*(EijR(3,1)+EijR(3,2)-EijR(4,1)+3*(EijR(4,2)-EijR(
     -   4,3))+EijR(4,4)+EijR(6,3)-EijR(7,3)-4*(EijR(10,2)-EijR(14,3)
     -   )+5*EijR(16,3)+EijR(18,4)-6*EijR(19,3)+EijR(19,4)+EijR(22,4)
     -   +2*(EijR(6,2)-EijR(7,2)+EijR(11,3)-EijR(13,3)-EijR(14,4)-Eij
     -   R(16,4)-EijR(25,4)-EijR(30,4)+2*EijR(33,4)))
       FI(320) = -16*(EijI(3,1)+EijI(3,2)-EijI(4,1)+3*(EijI(4,2)-EijI(
     -   4,3))+EijI(4,4)+EijI(6,3)-EijI(7,3)-4*(EijI(10,2)-EijI(14,3)
     -   )+5*EijI(16,3)+EijI(18,4)-6*EijI(19,3)+EijI(19,4)+EijI(22,4)
     -   +2*(EijI(6,2)-EijI(7,2)+EijI(11,3)-EijI(13,3)-EijI(14,4)-Eij
     -   I(16,4)-EijI(25,4)-EijI(30,4)+2*EijI(33,4)))
       F(320)=DCMPLX(FR(320),FI(320))
       FR(321) = -16*(EijR(3,2)+EijR(4,2)+EijR(4,4)+4*(EijR(16,3)-EijR
     -   (20,3))+EijR(20,4)+EijR(21,4)+EijR(22,4)-2*(EijR(4,3)+EijR(1
     -   0,2)-EijR(12,3)+EijR(13,3)-EijR(15,3)+EijR(15,4)+EijR(16,4)+
     -   EijR(28,4)+EijR(31,4)-2*EijR(34,4)))
       FI(321) = -16*(EijI(3,2)+EijI(4,2)+EijI(4,4)+4*(EijI(16,3)-EijI
     -   (20,3))+EijI(20,4)+EijI(21,4)+EijI(22,4)-2*(EijI(4,3)+EijI(1
     -   0,2)-EijI(12,3)+EijI(13,3)-EijI(15,3)+EijI(15,4)+EijI(16,4)+
     -   EijI(28,4)+EijI(31,4)-2*EijI(34,4)))
       F(321)=DCMPLX(FR(321),FI(321))
       FR(322) = 16*(EijR(3,1)-EijR(4,1)-EijR(4,4)+EijR(6,3)-EijR(7,2)
     -   -EijR(10,2)-EijR(16,3)+EijR(16,4)-EijR(19,4)+EijR(25,4)+2*(E
     -   ijR(4,3)+EijR(6,2)-EijR(14,3)+EijR(14,4)-EijR(33,4)))
       FI(322) = 16*(EijI(3,1)-EijI(4,1)-EijI(4,4)+EijI(6,3)-EijI(7,2)
     -   -EijI(10,2)-EijI(16,3)+EijI(16,4)-EijI(19,4)+EijI(25,4)+2*(E
     -   ijI(4,3)+EijI(6,2)-EijI(14,3)+EijI(14,4)-EijI(33,4)))
       F(322)=DCMPLX(FR(322),FI(322))
       FR(323) = 16*(EijR(2,1)+EijR(3,1)+EijR(4,3)-EijR(4,4)+EijR(5,2)
     -   +EijR(6,2)+EijR(8,2)-2*(EijR(4,1)+EijR(9,2))-EijR(10,2)+EijR
     -   (14,4)+EijR(15,4)-EijR(16,3)+EijR(16,4)+EijR(17,3)-EijR(18,3
     -   )-EijR(32,4)-EijR(33,4)-EijR(34,4)+EijR(35,4))
       FI(323) = 16*(EijI(2,1)+EijI(3,1)+EijI(4,3)-EijI(4,4)+EijI(5,2)
     -   +EijI(6,2)+EijI(8,2)-2*(EijI(4,1)+EijI(9,2))-EijI(10,2)+EijI
     -   (14,4)+EijI(15,4)-EijI(16,3)+EijI(16,4)+EijI(17,3)-EijI(18,3
     -   )-EijI(32,4)-EijI(33,4)-EijI(34,4)+EijI(35,4))
       F(323)=DCMPLX(FR(323),FI(323))
       FR(324) = 16*(EijR(2,1)+EijR(2,2)-EijR(4,1)+EijR(4,2)+EijR(4,3)
     -   -EijR(4,4)+EijR(8,2)-3*EijR(9,2)+EijR(9,3)-EijR(10,3)-EijR(1
     -   6,3)+EijR(16,4)-EijR(21,4)+EijR(28,4)+2*(EijR(15,4)-EijR(34,
     -   4)))
       FI(324) = 16*(EijI(2,1)+EijI(2,2)-EijI(4,1)+EijI(4,2)+EijI(4,3)
     -   -EijI(4,4)+EijI(8,2)-3*EijI(9,2)+EijI(9,3)-EijI(10,3)-EijI(1
     -   6,3)+EijI(16,4)-EijI(21,4)+EijI(28,4)+2*(EijI(15,4)-EijI(34,
     -   4)))
       F(324)=DCMPLX(FR(324),FI(324))
       FR(325) = EijR(3,2)-EijR(4,3)+EijR(4,4)-EijR(6,2)+EijR(7,2)-Eij
     -   R(10,2)+EijR(11,3)+EijR(12,3)-EijR(14,4)-EijR(15,4)+3*EijR(1
     -   6,3)-EijR(17,3)+EijR(18,3)-EijR(19,3)-EijR(20,3)+EijR(22,4)+
     -   EijR(29,4)-EijR(30,4)-EijR(31,4)+EijR(32,4)-2*(EijR(13,3)+Ei
     -   jR(16,4)-EijR(33,4)-EijR(34,4)+EijR(35,4))
       FI(325) = EijI(3,2)-EijI(4,3)+EijI(4,4)-EijI(6,2)+EijI(7,2)-Eij
     -   I(10,2)+EijI(11,3)+EijI(12,3)-EijI(14,4)-EijI(15,4)+3*EijI(1
     -   6,3)-EijI(17,3)+EijI(18,3)-EijI(19,3)-EijI(20,3)+EijI(22,4)+
     -   EijI(29,4)-EijI(30,4)-EijI(31,4)+EijI(32,4)-2*(EijI(13,3)+Ei
     -   jI(16,4)-EijI(33,4)-EijI(34,4)+EijI(35,4))
       F(325)=DCMPLX(FR(325),FI(325))
       FR(326) = EijR(4,2)+EijR(4,4)+EijR(6,2)-EijR(7,2)-EijR(10,2)+Ei
     -   jR(11,3)+EijR(12,3)+EijR(14,3)-EijR(14,4)+EijR(15,3)-EijR(15
     -   ,4)+4*EijR(16,3)+EijR(22,4)+EijR(29,4)-EijR(30,4)-EijR(31,4)
     -   +EijR(32,4)-2*(EijR(4,3)+EijR(13,3)+EijR(16,4)+EijR(19,3)+Ei
     -   jR(20,3)-EijR(33,4)-EijR(34,4)+EijR(35,4))
       FI(326) = EijI(4,2)+EijI(4,4)+EijI(6,2)-EijI(7,2)-EijI(10,2)+Ei
     -   jI(11,3)+EijI(12,3)+EijI(14,3)-EijI(14,4)+EijI(15,3)-EijI(15
     -   ,4)+4*EijI(16,3)+EijI(22,4)+EijI(29,4)-EijI(30,4)-EijI(31,4)
     -   +EijI(32,4)-2*(EijI(4,3)+EijI(13,3)+EijI(16,4)+EijI(19,3)+Ei
     -   jI(20,3)-EijI(33,4)-EijI(34,4)+EijI(35,4))
       F(326)=DCMPLX(FR(326),FI(326))
       FR(327) = EijR(3,3)-EijR(4,3)+EijR(4,4)+EijR(11,4)-EijR(13,4)-E
     -   ijR(14,4)-3*(EijR(13,3)-EijR(16,3)+EijR(16,4)-EijR(22,4)+Eij
     -   R(30,4)-EijR(33,4))
       FI(327) = EijI(3,3)-EijI(4,3)+EijI(4,4)+EijI(11,4)-EijI(13,4)-E
     -   ijI(14,4)-3*(EijI(13,3)-EijI(16,3)+EijI(16,4)-EijI(22,4)+Eij
     -   I(30,4)-EijI(33,4))
       F(327)=DCMPLX(FR(327),FI(327))
       FR(328) = -EijR(4,3)+EijR(4,4)-EijR(13,3)-EijR(14,4)+EijR(22,4)
     -   -EijR(30,4)+2*(EijR(16,3)-EijR(16,4)+EijR(33,4))
       FI(328) = -EijI(4,3)+EijI(4,4)-EijI(13,3)-EijI(14,4)+EijI(22,4)
     -   -EijI(30,4)+2*(EijI(16,3)-EijI(16,4)+EijI(33,4))
       F(328)=DCMPLX(FR(328),FI(328))
       FR(329) = EijR(4,2)-EijR(4,3)+EijR(4,4)-EijR(7,2)-EijR(14,4)-Ei
     -   jR(15,4)+2*EijR(16,3)-EijR(16,4)+EijR(18,3)-EijR(19,3)-EijR(
     -   20,3)+EijR(32,4)+EijR(33,4)+EijR(34,4)-EijR(35,4)
       FI(329) = EijI(4,2)-EijI(4,3)+EijI(4,4)-EijI(7,2)-EijI(14,4)-Ei
     -   jI(15,4)+2*EijI(16,3)-EijI(16,4)+EijI(18,3)-EijI(19,3)-EijI(
     -   20,3)+EijI(32,4)+EijI(33,4)+EijI(34,4)-EijI(35,4)
       F(329)=DCMPLX(FR(329),FI(329))
       FR(330) = -EijR(4,4)+EijR(13,3)+EijR(14,3)+EijR(14,4)-EijR(16,3
     -   )-EijR(19,3)-EijR(22,4)+EijR(30,4)+2*(EijR(16,4)-EijR(33,4))
       FI(330) = -EijI(4,4)+EijI(13,3)+EijI(14,3)+EijI(14,4)-EijI(16,3
     -   )-EijI(19,3)-EijI(22,4)+EijI(30,4)+2*(EijI(16,4)-EijI(33,4))
       F(330)=DCMPLX(FR(330),FI(330))
       FR(331) = -EijR(4,4)+EijR(14,3)+EijR(14,4)-EijR(16,3)+EijR(16,4
     -   )-EijR(33,4)
       FI(331) = -EijI(4,4)+EijI(14,3)+EijI(14,4)-EijI(16,3)+EijI(16,4
     -   )-EijI(33,4)
       F(331)=DCMPLX(FR(331),FI(331))
       FR(332) = EijR(4,4)+EijR(5,2)-EijR(7,2)+EijR(8,2)+EijR(8,3)+Eij
     -   R(9,3)-EijR(10,2)-EijR(14,4)+5*EijR(15,3)-EijR(16,4)+EijR(17
     -   ,3)-EijR(19,3)-3*(EijR(4,3)+EijR(18,3)+EijR(20,3))+EijR(21,4
     -   )+EijR(26,4)-EijR(27,4)-EijR(28,4)+EijR(33,4)+2*(EijR(4,2)-E
     -   ijR(9,2)-EijR(10,3)+EijR(14,3)-EijR(15,4)+EijR(16,3)+EijR(32
     -   ,4)+EijR(34,4)-EijR(35,4))
       FI(332) = EijI(4,4)+EijI(5,2)-EijI(7,2)+EijI(8,2)+EijI(8,3)+Eij
     -   I(9,3)-EijI(10,2)-EijI(14,4)+5*EijI(15,3)-EijI(16,4)+EijI(17
     -   ,3)-EijI(19,3)-3*(EijI(4,3)+EijI(18,3)+EijI(20,3))+EijI(21,4
     -   )+EijI(26,4)-EijI(27,4)-EijI(28,4)+EijI(33,4)+2*(EijI(4,2)-E
     -   ijI(9,2)-EijI(10,3)+EijI(14,3)-EijI(15,4)+EijI(16,3)+EijI(32
     -   ,4)+EijI(34,4)-EijI(35,4))
       F(332)=DCMPLX(FR(332),FI(332))
       FR(333) = EijR(4,2)+EijR(4,4)+EijR(6,2)-EijR(7,2)-EijR(10,2)+Ei
     -   jR(12,3)-EijR(13,3)+EijR(14,3)-EijR(14,4)-EijR(15,4)+EijR(17
     -   ,3)-EijR(18,3)-EijR(19,3)+3*(EijR(16,3)-EijR(20,3))+EijR(22,
     -   4)+EijR(29,4)-EijR(30,4)-EijR(31,4)+EijR(32,4)-2*(EijR(4,3)-
     -   EijR(15,3)+EijR(16,4)-EijR(33,4)-EijR(34,4)+EijR(35,4))
       FI(333) = EijI(4,2)+EijI(4,4)+EijI(6,2)-EijI(7,2)-EijI(10,2)+Ei
     -   jI(12,3)-EijI(13,3)+EijI(14,3)-EijI(14,4)-EijI(15,4)+EijI(17
     -   ,3)-EijI(18,3)-EijI(19,3)+3*(EijI(16,3)-EijI(20,3))+EijI(22,
     -   4)+EijI(29,4)-EijI(30,4)-EijI(31,4)+EijI(32,4)-2*(EijI(4,3)-
     -   EijI(15,3)+EijI(16,4)-EijI(33,4)-EijI(34,4)+EijI(35,4))
       F(333)=DCMPLX(FR(333),FI(333))
       FR(334) = EijR(3,1)-EijR(4,1)+EijR(4,2)+EijR(4,3)-EijR(4,4)+Eij
     -   R(8,2)-EijR(9,2)-EijR(10,2)-EijR(15,3)-EijR(16,3)+EijR(16,4)
     -   +EijR(20,3)-EijR(21,4)+EijR(28,4)+2*(EijR(15,4)-EijR(34,4))
       FI(334) = EijI(3,1)-EijI(4,1)+EijI(4,2)+EijI(4,3)-EijI(4,4)+Eij
     -   I(8,2)-EijI(9,2)-EijI(10,2)-EijI(15,3)-EijI(16,3)+EijI(16,4)
     -   +EijI(20,3)-EijI(21,4)+EijI(28,4)+2*(EijI(15,4)-EijI(34,4))
       F(334)=DCMPLX(FR(334),FI(334))
       FR(335) = EijR(3,2)+EijR(4,2)+EijR(4,4)-EijR(6,2)+EijR(7,2)+Eij
     -   R(8,2)-EijR(9,2)+EijR(12,3)-EijR(13,3)+EijR(14,3)-EijR(14,4)
     -   -EijR(15,4)+EijR(17,3)-EijR(18,3)-EijR(19,3)+3*(EijR(16,3)-E
     -   ijR(20,3))+EijR(22,4)+EijR(29,4)-EijR(30,4)-EijR(31,4)+EijR(
     -   32,4)-2*(EijR(4,3)+EijR(10,2)-EijR(15,3)+EijR(16,4)-EijR(33,
     -   4)-EijR(34,4)+EijR(35,4))
       FI(335) = EijI(3,2)+EijI(4,2)+EijI(4,4)-EijI(6,2)+EijI(7,2)+Eij
     -   I(8,2)-EijI(9,2)+EijI(12,3)-EijI(13,3)+EijI(14,3)-EijI(14,4)
     -   -EijI(15,4)+EijI(17,3)-EijI(18,3)-EijI(19,3)+3*(EijI(16,3)-E
     -   ijI(20,3))+EijI(22,4)+EijI(29,4)-EijI(30,4)-EijI(31,4)+EijI(
     -   32,4)-2*(EijI(4,3)+EijI(10,2)-EijI(15,3)+EijI(16,4)-EijI(33,
     -   4)-EijI(34,4)+EijI(35,4))
       F(335)=DCMPLX(FR(335),FI(335))
       FR(336) = EijR(3,1)+EijR(3,2)-EijR(4,1)+EijR(4,3)-EijR(4,4)-Eij
     -   R(10,2)-EijR(15,3)+EijR(15,4)-EijR(16,3)+EijR(20,3)-EijR(22,
     -   4)+EijR(31,4)+2*(EijR(16,4)-EijR(34,4))
       FI(336) = EijI(3,1)+EijI(3,2)-EijI(4,1)+EijI(4,3)-EijI(4,4)-Eij
     -   I(10,2)-EijI(15,3)+EijI(15,4)-EijI(16,3)+EijI(20,3)-EijI(22,
     -   4)+EijI(31,4)+2*(EijI(16,4)-EijI(34,4))
       F(336)=DCMPLX(FR(336),FI(336))
       FR(337) = EijR(2,1)+EijR(2,2)-EijR(4,1)+EijR(4,2)-EijR(4,3)+Eij
     -   R(4,4)+EijR(15,3)+EijR(16,3)-EijR(16,4)-EijR(20,3)+EijR(21,4
     -   )-EijR(28,4)-2*(EijR(9,2)+EijR(15,4)-EijR(34,4))
       FI(337) = EijI(2,1)+EijI(2,2)-EijI(4,1)+EijI(4,2)-EijI(4,3)+Eij
     -   I(4,4)+EijI(15,3)+EijI(16,3)-EijI(16,4)-EijI(20,3)+EijI(21,4
     -   )-EijI(28,4)-2*(EijI(9,2)+EijI(15,4)-EijI(34,4))
       F(337)=DCMPLX(FR(337),FI(337))
       FR(338) = EijR(3,1)-EijR(4,1)+EijR(4,4)+EijR(8,2)-EijR(9,2)-Eij
     -   R(15,4)+EijR(22,4)-EijR(31,4)-2*(EijR(16,4)-EijR(34,4))
       FI(338) = EijI(3,1)-EijI(4,1)+EijI(4,4)+EijI(8,2)-EijI(9,2)-Eij
     -   I(15,4)+EijI(22,4)-EijI(31,4)-2*(EijI(16,4)-EijI(34,4))
       F(338)=DCMPLX(FR(338),FI(338))
       FR(339) = EijR(4,4)-EijR(9,2)+EijR(10,2)-EijR(15,4)-EijR(16,4)+
     -   EijR(34,4)
       FI(339) = EijI(4,4)-EijI(9,2)+EijI(10,2)-EijI(15,4)-EijI(16,4)+
     -   EijI(34,4)
       F(339)=DCMPLX(FR(339),FI(339))
       FR(340) = EijR(3,2)+EijR(4,4)+EijR(8,2)-EijR(9,2)+EijR(11,3)+Ei
     -   jR(12,3)-EijR(14,4)-EijR(15,4)+5*EijR(16,3)+EijR(17,3)-EijR(
     -   18,3)-3*(EijR(4,3)+EijR(10,2)+EijR(19,3)+EijR(20,3))+EijR(22
     -   ,4)+EijR(29,4)-EijR(30,4)-EijR(31,4)+EijR(32,4)+2*(EijR(4,2)
     -   -EijR(13,3)+EijR(14,3)+EijR(15,3)-EijR(16,4)+EijR(33,4)+EijR
     -   (34,4)-EijR(35,4))
       FI(340) = EijI(3,2)+EijI(4,4)+EijI(8,2)-EijI(9,2)+EijI(11,3)+Ei
     -   jI(12,3)-EijI(14,4)-EijI(15,4)+5*EijI(16,3)+EijI(17,3)-EijI(
     -   18,3)-3*(EijI(4,3)+EijI(10,2)+EijI(19,3)+EijI(20,3))+EijI(22
     -   ,4)+EijI(29,4)-EijI(30,4)-EijI(31,4)+EijI(32,4)+2*(EijI(4,2)
     -   -EijI(13,3)+EijI(14,3)+EijI(15,3)-EijI(16,4)+EijI(33,4)+EijI
     -   (34,4)-EijI(35,4))
       F(340)=DCMPLX(FR(340),FI(340))
       FR(341) = EijR(3,2)+EijR(3,3)+EijR(4,2)+EijR(4,4)+EijR(11,3)+Ei
     -   jR(11,4)-4*EijR(13,3)-EijR(13,4)+EijR(14,3)-EijR(14,4)+5*Eij
     -   R(16,3)-2*(EijR(4,3)+EijR(10,2)+EijR(19,3))-3*(EijR(16,4)-Ei
     -   jR(22,4)+EijR(30,4)-EijR(33,4))
       FI(341) = EijI(3,2)+EijI(3,3)+EijI(4,2)+EijI(4,4)+EijI(11,3)+Ei
     -   jI(11,4)-4*EijI(13,3)-EijI(13,4)+EijI(14,3)-EijI(14,4)+5*Eij
     -   I(16,3)-2*(EijI(4,3)+EijI(10,2)+EijI(19,3))-3*(EijI(16,4)-Ei
     -   jI(22,4)+EijI(30,4)-EijI(33,4))
       F(341)=DCMPLX(FR(341),FI(341))
       FR(342) = EijR(4,2)+EijR(4,4)-EijR(10,2)-EijR(13,3)+EijR(14,3)-
     -   EijR(14,4)+3*EijR(16,3)-EijR(19,3)+EijR(22,4)-EijR(30,4)-2*(
     -   EijR(4,3)+EijR(16,4)-EijR(33,4))
       FI(342) = EijI(4,2)+EijI(4,4)-EijI(10,2)-EijI(13,3)+EijI(14,3)-
     -   EijI(14,4)+3*EijI(16,3)-EijI(19,3)+EijI(22,4)-EijI(30,4)-2*(
     -   EijI(4,3)+EijI(16,4)-EijI(33,4))
       F(342)=DCMPLX(FR(342),FI(342))
       FR(343) = EijR(3,1)+EijR(3,2)-EijR(4,1)+EijR(4,4)+EijR(6,2)-Eij
     -   R(7,2)-3*EijR(10,2)+EijR(11,3)+EijR(12,3)+EijR(14,3)-EijR(14
     -   ,4)+EijR(15,3)-EijR(15,4)+4*EijR(16,3)+EijR(22,4)+EijR(29,4)
     -   -EijR(30,4)-EijR(31,4)+EijR(32,4)+2*(EijR(4,2)-EijR(4,3)-Eij
     -   R(13,3)-EijR(16,4)-EijR(19,3)-EijR(20,3)+EijR(33,4)+EijR(34,
     -   4)-EijR(35,4))
       FI(343) = EijI(3,1)+EijI(3,2)-EijI(4,1)+EijI(4,4)+EijI(6,2)-Eij
     -   I(7,2)-3*EijI(10,2)+EijI(11,3)+EijI(12,3)+EijI(14,3)-EijI(14
     -   ,4)+EijI(15,3)-EijI(15,4)+4*EijI(16,3)+EijI(22,4)+EijI(29,4)
     -   -EijI(30,4)-EijI(31,4)+EijI(32,4)+2*(EijI(4,2)-EijI(4,3)-Eij
     -   I(13,3)-EijI(16,4)-EijI(19,3)-EijI(20,3)+EijI(33,4)+EijI(34,
     -   4)-EijI(35,4))
       F(343)=DCMPLX(FR(343),FI(343))
       FR(344) = EijR(3,3)-EijR(4,3)+EijR(4,4)+EijR(12,4)-EijR(13,4)-E
     -   ijR(15,4)-3*(EijR(13,3)-EijR(16,3)+EijR(16,4)-EijR(22,4)+Eij
     -   R(31,4)-EijR(34,4))
       FI(344) = EijI(3,3)-EijI(4,3)+EijI(4,4)+EijI(12,4)-EijI(13,4)-E
     -   ijI(15,4)-3*(EijI(13,3)-EijI(16,3)+EijI(16,4)-EijI(22,4)+Eij
     -   I(31,4)-EijI(34,4))
       F(344)=DCMPLX(FR(344),FI(344))
       FR(345) = EijR(4,2)+EijR(4,3)-EijR(4,4)-EijR(10,2)+EijR(13,3)+E
     -   ijR(15,4)-EijR(22,4)+EijR(31,4)-2*(EijR(16,3)-EijR(16,4)+Eij
     -   R(34,4))
       FI(345) = EijI(4,2)+EijI(4,3)-EijI(4,4)-EijI(10,2)+EijI(13,3)+E
     -   ijI(15,4)-EijI(22,4)+EijI(31,4)-2*(EijI(16,3)-EijI(16,4)+Eij
     -   I(34,4))
       F(345)=DCMPLX(FR(345),FI(345))
       FR(346) = EijR(3,2)+EijR(3,3)+EijR(4,2)+EijR(4,4)+EijR(12,3)+Ei
     -   jR(12,4)-4*EijR(13,3)-EijR(13,4)+EijR(15,3)-EijR(15,4)+5*Eij
     -   R(16,3)-2*(EijR(4,3)+EijR(10,2)+EijR(20,3))-3*(EijR(16,4)-Ei
     -   jR(22,4)+EijR(31,4)-EijR(34,4))
       FI(346) = EijI(3,2)+EijI(3,3)+EijI(4,2)+EijI(4,4)+EijI(12,3)+Ei
     -   jI(12,4)-4*EijI(13,3)-EijI(13,4)+EijI(15,3)-EijI(15,4)+5*Eij
     -   I(16,3)-2*(EijI(4,3)+EijI(10,2)+EijI(20,3))-3*(EijI(16,4)-Ei
     -   jI(22,4)+EijI(31,4)-EijI(34,4))
       F(346)=DCMPLX(FR(346),FI(346))
       FR(347) = EijR(3,3)+EijR(3,4)-EijR(4,3)+EijR(4,4)-3*(EijR(13,3)
     -   -EijR(16,3))-4*(EijR(13,4)+EijR(16,4))+6*EijR(22,4)
       FI(347) = EijI(3,3)+EijI(3,4)-EijI(4,3)+EijI(4,4)-3*(EijI(13,3)
     -   -EijI(16,3))-4*(EijI(13,4)+EijI(16,4))+6*EijI(22,4)
       F(347)=DCMPLX(FR(347),FI(347))
       FR(348) = -EijR(4,3)+EijR(4,4)-EijR(13,3)-EijR(13,4)+2*EijR(16,
     -   3)-3*(EijR(16,4)-EijR(22,4))
       FI(348) = -EijI(4,3)+EijI(4,4)-EijI(13,3)-EijI(13,4)+2*EijI(16,
     -   3)-3*(EijI(16,4)-EijI(22,4))
       F(348)=DCMPLX(FR(348),FI(348))
       FR(349) = EijR(4,2)-EijR(4,3)+EijR(4,4)-EijR(10,2)-EijR(13,3)-E
     -   ijR(15,4)+EijR(22,4)-EijR(31,4)+2*(EijR(16,3)-EijR(16,4)+Eij
     -   R(34,4))
       FI(349) = EijI(4,2)-EijI(4,3)+EijI(4,4)-EijI(10,2)-EijI(13,3)-E
     -   ijI(15,4)+EijI(22,4)-EijI(31,4)+2*(EijI(16,3)-EijI(16,4)+Eij
     -   I(34,4))
       F(349)=DCMPLX(FR(349),FI(349))
       FR(350) = EijR(4,2)+2*EijR(4,3)-EijR(4,4)-EijR(7,2)+EijR(8,2)-E
     -   ijR(9,2)-EijR(14,3)+EijR(14,4)-EijR(15,3)+EijR(15,4)-EijR(16
     -   ,3)+EijR(16,4)+EijR(17,3)-EijR(32,4)-EijR(33,4)-EijR(34,4)+E
     -   ijR(35,4)
       FI(350) = EijI(4,2)+2*EijI(4,3)-EijI(4,4)-EijI(7,2)+EijI(8,2)-E
     -   ijI(9,2)-EijI(14,3)+EijI(14,4)-EijI(15,3)+EijI(15,4)-EijI(16
     -   ,3)+EijI(16,4)+EijI(17,3)-EijI(32,4)-EijI(33,4)-EijI(34,4)+E
     -   ijI(35,4)
       F(350)=DCMPLX(FR(350),FI(350))
       FR(351) = EijR(3,2)+EijR(4,2)+EijR(4,3)-EijR(4,4)+EijR(11,3)+Ei
     -   jR(14,4)-EijR(16,3)-EijR(19,3)-EijR(22,4)+EijR(30,4)-2*(EijR
     -   (10,2)-EijR(16,4)+EijR(33,4))
       FI(351) = EijI(3,2)+EijI(4,2)+EijI(4,3)-EijI(4,4)+EijI(11,3)+Ei
     -   jI(14,4)-EijI(16,3)-EijI(19,3)-EijI(22,4)+EijI(30,4)-2*(EijI
     -   (10,2)-EijI(16,4)+EijI(33,4))
       F(351)=DCMPLX(FR(351),FI(351))
       FR(352) = EijR(4,2)+EijR(4,3)-EijR(4,4)-EijR(10,2)+EijR(14,4)+E
     -   ijR(16,4)-EijR(19,3)-EijR(33,4)
       FI(352) = EijI(4,2)+EijI(4,3)-EijI(4,4)-EijI(10,2)+EijI(14,4)+E
     -   ijI(16,4)-EijI(19,3)-EijI(33,4)
       F(352)=DCMPLX(FR(352),FI(352))
       FR(353) = EijR(3,1)+EijR(3,2)-EijR(4,1)+EijR(4,2)-EijR(4,4)+Eij
     -   R(8,2)-EijR(9,2)+EijR(12,3)+EijR(15,3)+EijR(15,4)-EijR(22,4)
     -   +EijR(31,4)-2*(EijR(10,2)-EijR(16,4)+EijR(20,3)+EijR(34,4))
       FI(353) = EijI(3,1)+EijI(3,2)-EijI(4,1)+EijI(4,2)-EijI(4,4)+Eij
     -   I(8,2)-EijI(9,2)+EijI(12,3)+EijI(15,3)+EijI(15,4)-EijI(22,4)
     -   +EijI(31,4)-2*(EijI(10,2)-EijI(16,4)+EijI(20,3)+EijI(34,4))
       F(353)=DCMPLX(FR(353),FI(353))
       FR(354) = 2*EijR(4,2)-EijR(4,4)-EijR(9,2)-EijR(10,2)+EijR(15,3)
     -   +EijR(15,4)+EijR(16,4)-EijR(20,3)-EijR(34,4)
       FI(354) = 2*EijI(4,2)-EijI(4,4)-EijI(9,2)-EijI(10,2)+EijI(15,3)
     -   +EijI(15,4)+EijI(16,4)-EijI(20,3)-EijI(34,4)
       F(354)=DCMPLX(FR(354),FI(354))
       FR(355) = EijR(3,1)+EijR(3,2)-EijR(4,1)+EijR(4,2)+EijR(4,3)-Eij
     -   R(4,4)+EijR(6,2)-EijR(7,2)+EijR(11,3)+EijR(14,4)-EijR(16,3)-
     -   EijR(19,3)-EijR(22,4)+EijR(30,4)-2*(EijR(10,2)-EijR(16,4)+Ei
     -   jR(33,4))
       FI(355) = EijI(3,1)+EijI(3,2)-EijI(4,1)+EijI(4,2)+EijI(4,3)-Eij
     -   I(4,4)+EijI(6,2)-EijI(7,2)+EijI(11,3)+EijI(14,4)-EijI(16,3)-
     -   EijI(19,3)-EijI(22,4)+EijI(30,4)-2*(EijI(10,2)-EijI(16,4)+Ei
     -   jI(33,4))
       F(355)=DCMPLX(FR(355),FI(355))
       FR(356) = EijR(4,2)+EijR(4,3)-EijR(4,4)+EijR(8,2)-EijR(9,2)-Eij
     -   R(10,2)+EijR(12,3)+EijR(15,4)-EijR(16,3)-EijR(20,3)-EijR(22,
     -   4)+EijR(31,4)+2*(EijR(16,4)-EijR(34,4))
       FI(356) = EijI(4,2)+EijI(4,3)-EijI(4,4)+EijI(8,2)-EijI(9,2)-Eij
     -   I(10,2)+EijI(12,3)+EijI(15,4)-EijI(16,3)-EijI(20,3)-EijI(22,
     -   4)+EijI(31,4)+2*(EijI(16,4)-EijI(34,4))
       F(356)=DCMPLX(FR(356),FI(356))
       FR(357) = EijR(3,2)+EijR(3,3)+EijR(4,2)-EijR(4,4)-2*(EijR(10,2)
     -   +EijR(13,3))+EijR(13,4)+EijR(16,3)+3*(EijR(16,4)-EijR(22,4))
       FI(357) = EijI(3,2)+EijI(3,3)+EijI(4,2)-EijI(4,4)-2*(EijI(10,2)
     -   +EijI(13,3))+EijI(13,4)+EijI(16,3)+3*(EijI(16,4)-EijI(22,4))
       F(357)=DCMPLX(FR(357),FI(357))
       FR(358) = EijR(4,2)-EijR(4,4)-EijR(10,2)-EijR(13,3)+EijR(16,3)+
     -   2*EijR(16,4)-EijR(22,4)
       FI(358) = EijI(4,2)-EijI(4,4)-EijI(10,2)-EijI(13,3)+EijI(16,3)+
     -   2*EijI(16,4)-EijI(22,4)
       F(358)=DCMPLX(FR(358),FI(358))
       FR(359) = EijR(4,2)-EijR(4,4)-EijR(10,2)+EijR(14,3)+EijR(14,4)+
     -   EijR(16,4)-EijR(19,3)-EijR(33,4)
       FI(359) = EijI(4,2)-EijI(4,4)-EijI(10,2)+EijI(14,3)+EijI(14,4)+
     -   EijI(16,4)-EijI(19,3)-EijI(33,4)
       F(359)=DCMPLX(FR(359),FI(359))
       FR(360) = EijR(4,4)-EijR(15,3)-EijR(15,4)-EijR(16,4)+EijR(20,3)
     -   +EijR(34,4)
       FI(360) = EijI(4,4)-EijI(15,3)-EijI(15,4)-EijI(16,4)+EijI(20,3)
     -   +EijI(34,4)
       F(360)=DCMPLX(FR(360),FI(360))
       FR(361) = -EijR(4,3)-EijR(4,4)-EijR(13,3)+2*(EijR(16,3)+EijR(16
     -   ,4))-EijR(22,4)
       FI(361) = -EijI(4,3)-EijI(4,4)-EijI(13,3)+2*(EijI(16,3)+EijI(16
     -   ,4))-EijI(22,4)
       F(361)=DCMPLX(FR(361),FI(361))
       FR(362) = -EijR(4,3)-EijR(4,4)+EijR(16,3)+EijR(16,4)
       FI(362) = -EijI(4,3)-EijI(4,4)+EijI(16,3)+EijI(16,4)
       F(362)=DCMPLX(FR(362),FI(362))
       P(392) = p2sq-s12-s15+s23+s45-2*P(375)
       FR(363) = -8*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)-2*(p3sq*EijR(
     -   3,2)+4*EijR(11,2)-EijR(22,3)-EijR(23,3)+2*EijR(24,3)-EijR(39
     -   ,4)-EijR(43,4)+EijR(44,4)+EijR(45,4))-EijR(6,2)*P(14)+(EijR(
     -   2,1)-EijR(4,1))*P(55)-EijR(2,2)*P(166)-EijR(4,2)*P(213)+EijR
     -   (5,2)*P(269)-EijR(7,2)*P(270)+EijR(10,2)*P(379)+EijR(8,2)*P(
     -   385)+EijR(9,2)*P(392))
       FI(363) = -8*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)-2*(p3sq*EijI(
     -   3,2)+4*EijI(11,2)-EijI(22,3)-EijI(23,3)+2*EijI(24,3)-EijI(39
     -   ,4)-EijI(43,4)+EijI(44,4)+EijI(45,4))-EijI(6,2)*P(14)+(EijI(
     -   2,1)-EijI(4,1))*P(55)-EijI(2,2)*P(166)-EijI(4,2)*P(213)+EijI
     -   (5,2)*P(269)-EijI(7,2)*P(270)+EijI(10,2)*P(379)+EijI(8,2)*P(
     -   385)+EijI(9,2)*P(392))
       F(363)=DCMPLX(FR(363),FI(363))
       FR(364) = -8*(s12*(EijR(3,1)-EijR(4,1)-EijR(6,2)+EijR(7,2)+2*(E
     -   ijR(8,2)-EijR(9,2)))+s34*(-EijR(3,1)+EijR(4,1)-EijR(8,2)+Eij
     -   R(9,2))+2*(EijR(23,3)-EijR(24,3)+EijR(38,4)+EijR(39,4)-2*Eij
     -   R(45,4))+(EijR(6,2)-EijR(7,2)-EijR(8,2)+EijR(9,2))*P(147))
       FI(364) = -8*(s12*(EijI(3,1)-EijI(4,1)-EijI(6,2)+EijI(7,2)+2*(E
     -   ijI(8,2)-EijI(9,2)))+s34*(-EijI(3,1)+EijI(4,1)-EijI(8,2)+Eij
     -   I(9,2))+2*(EijI(23,3)-EijI(24,3)+EijI(38,4)+EijI(39,4)-2*Eij
     -   I(45,4))+(EijI(6,2)-EijI(7,2)-EijI(8,2)+EijI(9,2))*P(147))
       F(364)=DCMPLX(FR(364),FI(364))
       P(393) = p2sq-s15+s34
       FR(365) = -8*(D02345R+Dij2345R(2,1)+s15*EijR(7,2)-p2sq*(EijR(6,
     -   2)+EijR(7,2)-EijR(8,2))+s12*(EijR(2,1)-EijR(3,2)-EijR(4,1)+E
     -   ijR(6,2)+EijR(7,2)+EijR(8,2)-2*EijR(9,2))-2*(EijR(11,2)+EijR
     -   (24,3)-EijR(39,4)+EijR(45,4))+EijR(10,2)*P(92)-(EijR(3,1)+Ei
     -   jR(3,2)-EijR(4,1))*P(114)+EijR(9,2)*P(393))
       FI(365) = -8*(D02345I+Dij2345I(2,1)+s15*EijI(7,2)-p2sq*(EijI(6,
     -   2)+EijI(7,2)-EijI(8,2))+s12*(EijI(2,1)-EijI(3,2)-EijI(4,1)+E
     -   ijI(6,2)+EijI(7,2)+EijI(8,2)-2*EijI(9,2))-2*(EijI(11,2)+EijI
     -   (24,3)-EijI(39,4)+EijI(45,4))+EijI(10,2)*P(92)-(EijI(3,1)+Ei
     -   jI(3,2)-EijI(4,1))*P(114)+EijI(9,2)*P(393))
       F(365)=DCMPLX(FR(365),FI(365))
       P(394) = p2sq-p3sq-2*s12+s45
       P(395) = p2sq-p3sq-s12-s15+s34+s45
       FR(366) = 8*(D02345R+Dij2345R(2,1)+s15*(EE0R+2*EijR(1,1)+EijR(1
     -   ,2))+s12*EijR(3,1)-2*(3*EijR(11,2)+EijR(21,3)+EijR(23,3)-2*E
     -   ijR(24,3)+EijR(39,4)+EijR(41,4)-EijR(42,4)-EijR(45,4))+EijR(
     -   4,1)*P(55)-EijR(10,2)*P(90)-EijR(3,2)*P(128)+EijR(8,2)*P(152
     -   )-EijR(9,2)*P(283)-EijR(2,1)*P(289)-EijR(5,2)*P(292)-EijR(6,
     -   2)*P(394)+EijR(7,2)*P(395))
       FI(366) = 8*(D02345I+Dij2345I(2,1)+s15*(EE0I+2*EijI(1,1)+EijI(1
     -   ,2))+s12*EijI(3,1)-2*(3*EijI(11,2)+EijI(21,3)+EijI(23,3)-2*E
     -   ijI(24,3)+EijI(39,4)+EijI(41,4)-EijI(42,4)-EijI(45,4))+EijI(
     -   4,1)*P(55)-EijI(10,2)*P(90)-EijI(3,2)*P(128)+EijI(8,2)*P(152
     -   )-EijI(9,2)*P(283)-EijI(2,1)*P(289)-EijI(5,2)*P(292)-EijI(6,
     -   2)*P(394)+EijI(7,2)*P(395))
       F(366)=DCMPLX(FR(366),FI(366))
       FR(367) = -8*(2*(EijR(38,4)+EijR(39,4)-2*EijR(45,4))+(EijR(3,1)
     -   -EijR(4,1)+EijR(6,2)-EijR(7,2))*P(9)+(EijR(4,2)-EijR(10,2))*
     -   P(283))
       FI(367) = -8*(2*(EijI(38,4)+EijI(39,4)-2*EijI(45,4))+(EijI(3,1)
     -   -EijI(4,1)+EijI(6,2)-EijI(7,2))*P(9)+(EijI(4,2)-EijI(10,2))*
     -   P(283))
       F(367)=DCMPLX(FR(367),FI(367))
       P(396) = p2sq+s34
       P(397) = -p3sq-3*s12-s15+s45+2*P(396)
       FR(368) = 8*(D02345R+Dij2345R(2,1)-p2sq*EijR(2,1)+s12*EijR(3,1)
     -   -2*(EijR(11,2)-EijR(21,3)+EijR(23,3)+EijR(39,4)+EijR(41,4)-E
     -   ijR(42,4)-EijR(45,4))+EijR(4,1)*P(9)-EijR(10,2)*P(90)-EijR(3
     -   ,2)*P(128)+(-EijR(5,2)+EijR(8,2))*P(152)-EijR(4,2)*P(283)-Ei
     -   jR(6,2)*P(394)+EijR(7,2)*P(397))
       FI(368) = 8*(D02345I+Dij2345I(2,1)-p2sq*EijI(2,1)+s12*EijI(3,1)
     -   -2*(EijI(11,2)-EijI(21,3)+EijI(23,3)+EijI(39,4)+EijI(41,4)-E
     -   ijI(42,4)-EijI(45,4))+EijI(4,1)*P(9)-EijI(10,2)*P(90)-EijI(3
     -   ,2)*P(128)+(-EijI(5,2)+EijI(8,2))*P(152)-EijI(4,2)*P(283)-Ei
     -   jI(6,2)*P(394)+EijI(7,2)*P(397))
       F(368)=DCMPLX(FR(368),FI(368))
       P(398) = p2sq-p4sq-2*s12-s15+s34-s45
       P(399) = p2sq-s12+s23-s45
       P(400) = p2sq+2*p3sq-3*s12-s15-s23-s45
       FR(369) = -8*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)+s15*(EE0R+Eij
     -   R(1,1))-2*(p3sq*EijR(3,2)+5*EijR(11,2)-EijR(39,4)-EijR(43,4)
     -   +EijR(44,4)+EijR(45,4))+(EijR(2,1)-EijR(4,1)+2*EijR(5,2))*P(
     -   9)-EijR(6,2)*P(14)-EijR(2,2)*P(152)+EijR(10,2)*P(379)+EijR(8
     -   ,2)*P(385)+EijR(4,2)*P(398)-EijR(7,2)*P(399)-EijR(9,2)*P(400
     -   ))
       FI(369) = -8*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)+s15*(EE0I+Eij
     -   I(1,1))-2*(p3sq*EijI(3,2)+5*EijI(11,2)-EijI(39,4)-EijI(43,4)
     -   +EijI(44,4)+EijI(45,4))+(EijI(2,1)-EijI(4,1)+2*EijI(5,2))*P(
     -   9)-EijI(6,2)*P(14)-EijI(2,2)*P(152)+EijI(10,2)*P(379)+EijI(8
     -   ,2)*P(385)+EijI(4,2)*P(398)-EijI(7,2)*P(399)-EijI(9,2)*P(400
     -   ))
       F(369)=DCMPLX(FR(369),FI(369))
       FR(370) = -8*(D02345R+Dij2345R(2,1)+s12*EijR(2,1)-2*(EijR(11,2)
     -   +EijR(23,3)-EijR(24,3)-EijR(39,4)+EijR(45,4))-(EijR(6,2)+Eij
     -   R(7,2))*P(9)+EijR(10,2)*P(92)-EijR(3,1)*P(114)-EijR(3,2)*P(1
     -   28)+EijR(8,2)*P(152)-EijR(4,1)*P(195)+EijR(4,2)*P(283))
       FI(370) = -8*(D02345I+Dij2345I(2,1)+s12*EijI(2,1)-2*(EijI(11,2)
     -   +EijI(23,3)-EijI(24,3)-EijI(39,4)+EijI(45,4))-(EijI(6,2)+Eij
     -   I(7,2))*P(9)+EijI(10,2)*P(92)-EijI(3,1)*P(114)-EijI(3,2)*P(1
     -   28)+EijI(8,2)*P(152)-EijI(4,1)*P(195)+EijI(4,2)*P(283))
       F(370)=DCMPLX(FR(370),FI(370))
       FR(371) = -8*(s12*(-EijR(3,1)+EijR(4,1)+EijR(6,2)-EijR(7,2)-2*(
     -   EijR(3,2)-EijR(10,2)))+p2sq*(EijR(3,2)-EijR(6,2)+EijR(7,2)-E
     -   ijR(10,2))+2*(EijR(23,3)-EijR(24,3)+EijR(38,4)+EijR(39,4)-2*
     -   EijR(45,4))+(-EijR(3,1)-EijR(3,2)+EijR(4,1)+EijR(10,2))*P(36
     -   1))
       FI(371) = -8*(s12*(-EijI(3,1)+EijI(4,1)+EijI(6,2)-EijI(7,2)-2*(
     -   EijI(3,2)-EijI(10,2)))+p2sq*(EijI(3,2)-EijI(6,2)+EijI(7,2)-E
     -   ijI(10,2))+2*(EijI(23,3)-EijI(24,3)+EijI(38,4)+EijI(39,4)-2*
     -   EijI(45,4))+(-EijI(3,1)-EijI(3,2)+EijI(4,1)+EijI(10,2))*P(36
     -   1))
       F(371)=DCMPLX(FR(371),FI(371))
       P(401) = p3sq-s15+s34-s45
       P(402) = p2sq-p3sq-s12+s34+s45
       P(403) = p2sq-p3sq-s15+s45-2*P(55)
       FR(372) = -8*(D02345R+Dij2345R(2,1)-p2sq*EijR(2,1)-2*(EijR(11,2
     -   )-EijR(21,3)-EijR(23,3)+2*EijR(24,3)-EijR(39,4)-EijR(41,4)+E
     -   ijR(42,4)+EijR(45,4))+EijR(3,1)*P(119)-EijR(3,2)*P(128)+EijR
     -   (4,1)*P(141)+(-EijR(5,2)+EijR(8,2))*P(152)+EijR(6,2)*P(401)+
     -   EijR(7,2)*P(402)-EijR(10,2)*P(403))
       FI(372) = -8*(D02345I+Dij2345I(2,1)-p2sq*EijI(2,1)-2*(EijI(11,2
     -   )-EijI(21,3)-EijI(23,3)+2*EijI(24,3)-EijI(39,4)-EijI(41,4)+E
     -   ijI(42,4)+EijI(45,4))+EijI(3,1)*P(119)-EijI(3,2)*P(128)+EijI
     -   (4,1)*P(141)+(-EijI(5,2)+EijI(8,2))*P(152)+EijI(6,2)*P(401)+
     -   EijI(7,2)*P(402)-EijI(10,2)*P(403))
       F(372)=DCMPLX(FR(372),FI(372))
       P(404) = -2*p3sq-3*s12-s15+s23+s34+s45
       P(405) = p2sq+2*p3sq+p4sq-3*s12-s15+s45
       FR(373) = 8*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)+s15*(EE0R+EijR
     -   (1,1))-2*(p3sq*EijR(3,2)+5*EijR(11,2)+EijR(23,3)-EijR(24,3)+
     -   EijR(39,4)+EijR(43,4)-EijR(44,4)-EijR(45,4))+(EijR(2,1)+2*Ei
     -   jR(5,2))*P(9)-EijR(6,2)*P(14)-EijR(4,1)*P(141)-EijR(2,2)*P(1
     -   52)-EijR(4,2)*P(213)-EijR(3,1)*P(283)+EijR(9,2)*P(384)-EijR(
     -   7,2)*P(399)-EijR(8,2)*P(404)+EijR(10,2)*P(405))
       FI(373) = 8*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)+s15*(EE0I+EijI
     -   (1,1))-2*(p3sq*EijI(3,2)+5*EijI(11,2)+EijI(23,3)-EijI(24,3)+
     -   EijI(39,4)+EijI(43,4)-EijI(44,4)-EijI(45,4))+(EijI(2,1)+2*Ei
     -   jI(5,2))*P(9)-EijI(6,2)*P(14)-EijI(4,1)*P(141)-EijI(2,2)*P(1
     -   52)-EijI(4,2)*P(213)-EijI(3,1)*P(283)+EijI(9,2)*P(384)-EijI(
     -   7,2)*P(399)-EijI(8,2)*P(404)+EijI(10,2)*P(405))
       F(373)=DCMPLX(FR(373),FI(373))
       FR(374) = 8*(D02345R+Dij2345R(2,1)+s12*EijR(2,1)-2*(EijR(11,2)+
     -   EijR(39,4)-EijR(45,4))-(EijR(6,2)+EijR(7,2))*P(9)-EijR(3,1)*
     -   P(114)+EijR(10,2)*P(127)-EijR(3,2)*P(128)+EijR(8,2)*P(152)-E
     -   ijR(4,1)*P(195))
       FI(374) = 8*(D02345I+Dij2345I(2,1)+s12*EijI(2,1)-2*(EijI(11,2)+
     -   EijI(39,4)-EijI(45,4))-(EijI(6,2)+EijI(7,2))*P(9)-EijI(3,1)*
     -   P(114)+EijI(10,2)*P(127)-EijI(3,2)*P(128)+EijI(8,2)*P(152)-E
     -   ijI(4,1)*P(195))
       F(374)=DCMPLX(FR(374),FI(374))
       P(406) = 2*p3sq-p4sq+s15-4*s23+s34+3*P(45)
       P(407) = 3*p2sq+p3sq-s23-2*P(141)
       P(408) = p2sq-s23-2*P(86)
       P(409) = s23-s34
       P(410) = p3sq+p4sq-s15+s45+3*P(9)-2*P(409)
       P(411) = p3sq-s12-s15-s45
       P(412) = s12*P(200)+p2sq*P(411)
       P(413) = s15-s34-2*P(18)
       P(414) = -p3sq+2*s12+s15+s23-s34
       P(415) = -3*p2sq**2+p2sq*P(413)+s12*P(414)
       P(416) = s45**2+p2sq*P(8)-s12*P(202)+p3sq*P(203)-s45*P(204)
       P(417) = p4sq+5*s23-s34+3*P(7)
       P(418) = -5*p3sq-s12+s45
       P(419) = -3*p3sq**2+P(7)*P(207)+p3sq*P(417)+p2sq*P(418)
       P(420) = -3*s12-s15+s34+2*P(18)
       P(421) = p3sq-3*s12-s15-2*P(409)
       P(422) = s45**2+3*p2sq*P(213)+p4sq*P(420)+s45*P(421)
       P(423) = -3*p3sq+4*s12+5*s23-s34-2*s45
       P(424) = s15+s23
       P(425) = 3*p3sq-7*s12+s34-s45-2*P(424)
       P(426) = 3*p2sq**2+s12*P(423)+p2sq*P(425)
       P(427) = 3*s12+s15
       P(428) = -3*p3sq+p4sq+7*s12+5*s23-s34-6*s45
       P(429) = s12**2-s23*s45
       P(430) = s23**2+s45**2
       P(431) = -3*p2sq**2-p4sq*s12+s12*s15-5*s12*s23+s12*s34+6*s12*s4
     -   5-s15*s45+p3sq*P(427)+p2sq*P(428)-4*P(429)-2*P(430)
       P(432) = p2sq-s12+s34
       P(433) = s23*P(432)
       P(434) = 3*s34+4*P(43)
       P(435) = p2sq*p4sq-p4sq*s12+s23**2+s45**2
       P(436) = -(p3sq*s15)+s15*s23+s15*s34-3*P(433)+s45*P(434)+2*P(43
     -   5)
       P(437) = p2sq**2+p3sq**2
       P(438) = 7*p3sq-p4sq+s12-5*s23+s34
       P(439) = s12*s23-s23**2-s12*s34
       P(440) = 2*s12+s15+5*s23-s34+s45
       P(441) = -(p4sq*s12)-4*s12**2-s12*s15+s45*P(227)+3*P(437)+p2sq*
     -   P(438)-2*P(439)-p3sq*P(440)
       P(442) = p3sq**2+s12**2
       P(443) = s12+s15+s23-s34
       P(444) = 7*s12+2*s15+5*s23-s34+s45
       P(445) = 6*p3sq+p4sq-3*P(95)
       P(446) = s12*s15+5*s12*s23+s15*s23+2*s23**2+s12*s34+s15*s34-s34
     -   **2+3*P(442)+s45*P(443)-p3sq*P(444)+p2sq*P(445)
       P(447) = p3sq**2+s12**2+s23*s34
       P(448) = s15+s23+2*P(71)
       P(449) = -2*s23+s34+s45
       P(450) = p4sq-7*s12-s15-5*s23+s34+4*s45
       P(451) = 6*p3sq-3*P(71)+2*P(213)
       P(452) = s12*s15+2*s12*s23+s12*s34-s34**2-s45**2+3*P(447)-s45*P
     -   (448)+p4sq*P(449)+p3sq*P(450)+p2sq*P(451)
       P(453) = 14*p2sq+9*p3sq-p4sq-11*s12-s15-12*s23+3*s34+4*s45
       P(454) = p2sq-p3sq+2*s12+s23
       P(455) = p3sq-p4sq-s34-2*s45
       FR(375) = 2*(-C0234R+EE0R*s15*P(2)+4*(Dij2345R(7,2)-EijR(23,3)*
     -   P(7))+s15*EijR(1,2)*P(9)+s15*EijR(1,1)*P(199)+EijR(4,1)*P(30
     -   2)+D02345R*P(406)-Dij2345R(1,1)*P(407)+Dij2345R(2,1)*P(408)+
     -   Dij2345R(3,1)*P(410)+EijR(2,1)*P(412)+EijR(2,2)*P(415)+EijR(
     -   3,1)*P(416)+EijR(3,2)*P(419)-EijR(4,2)*P(422)+EijR(5,2)*P(42
     -   6)+EijR(6,2)*P(431)+EijR(7,2)*P(436)+EijR(8,2)*P(441)-EijR(9
     -   ,2)*P(446)+EijR(10,2)*P(452)-2*(EijR(21,3)*P(399)+EijR(11,2)
     -   *P(453)-EijR(22,3)*P(454)-EijR(24,3)*P(455)))
       FI(375) = 2*(-C0234I+EE0I*s15*P(2)+4*(Dij2345I(7,2)-EijI(23,3)*
     -   P(7))+s15*EijI(1,2)*P(9)+s15*EijI(1,1)*P(199)+EijI(4,1)*P(30
     -   2)+D02345I*P(406)-Dij2345I(1,1)*P(407)+Dij2345I(2,1)*P(408)+
     -   Dij2345I(3,1)*P(410)+EijI(2,1)*P(412)+EijI(2,2)*P(415)+EijI(
     -   3,1)*P(416)+EijI(3,2)*P(419)-EijI(4,2)*P(422)+EijI(5,2)*P(42
     -   6)+EijI(6,2)*P(431)+EijI(7,2)*P(436)+EijI(8,2)*P(441)-EijI(9
     -   ,2)*P(446)+EijI(10,2)*P(452)-2*(EijI(21,3)*P(399)+EijI(11,2)
     -   *P(453)-EijI(22,3)*P(454)-EijI(24,3)*P(455)))
       F(375)=DCMPLX(FR(375),FI(375))
       P(456) = -p2sq+p3sq+3*p4sq+s23+2*P(280)
       P(457) = 2*p2sq+p3sq-s23
       P(458) = s23-2*P(27)+3*P(86)
       P(459) = 2*p3sq-s23-s34-s45
       P(460) = 2*p3sq-s12+s45
       P(461) = p2sq+p3sq-s12-s34-s45
       P(462) = p2sq-p3sq+s12+s34+s45
       FR(376) = -4*(-Dij2345R(3,2)+Dij2345R(5,2)+p3sq*EijR(3,1)+3*p3s
     -   q*EijR(6,2)+4*(EijR(11,2)-EijR(21,3))+2*(Dij2345R(1,1)-Dij23
     -   45R(3,1)-EijR(22,3))+6*EijR(24,3)-EijR(2,1)*P(2)+(EijR(5,3)-
     -   EijR(7,3))*P(9)+EijR(4,1)*P(43)+EijR(4,3)*P(90)+EijR(18,3)*P
     -   (117)-(EijR(8,2)-EijR(16,3)-EijR(17,3)+EijR(19,3)+EijR(20,3)
     -   )*P(128)+(EijR(2,2)-EijR(8,3)+EijR(10,3))*P(152)+EijR(4,2)*P
     -   (456)-EijR(5,2)*P(457)-EijR(7,2)*P(458)+EijR(9,2)*P(459)-Eij
     -   R(10,2)*P(460)+EijR(14,3)*P(461)-EijR(15,3)*P(462))
       FI(376) = -4*(-Dij2345I(3,2)+Dij2345I(5,2)+p3sq*EijI(3,1)+3*p3s
     -   q*EijI(6,2)+4*(EijI(11,2)-EijI(21,3))+2*(Dij2345I(1,1)-Dij23
     -   45I(3,1)-EijI(22,3))+6*EijI(24,3)-EijI(2,1)*P(2)+(EijI(5,3)-
     -   EijI(7,3))*P(9)+EijI(4,1)*P(43)+EijI(4,3)*P(90)+EijI(18,3)*P
     -   (117)-(EijI(8,2)-EijI(16,3)-EijI(17,3)+EijI(19,3)+EijI(20,3)
     -   )*P(128)+(EijI(2,2)-EijI(8,3)+EijI(10,3))*P(152)+EijI(4,2)*P
     -   (456)-EijI(5,2)*P(457)-EijI(7,2)*P(458)+EijI(9,2)*P(459)-Eij
     -   I(10,2)*P(460)+EijI(14,3)*P(461)-EijI(15,3)*P(462))
       F(376)=DCMPLX(FR(376),FI(376))
       P(463) = p4sq+s12-s34
       P(464) = -p2sq+p3sq+s23-3*s45+2*P(463)
       P(465) = -p2sq+p3sq+s12+s23-3*s45+2*P(79)
       P(466) = -p3sq+s34+s45+2*P(152)
       P(467) = p2sq+s12-2*P(92)
       FR(377) = -4*(-Dij2345R(1,2)-Dij2345R(3,2)+2*Dij2345R(5,2)-s12*
     -   EijR(2,2)+p3sq*EijR(8,2)-p3sq*EijR(10,2)-6*(EijR(22,3)-EijR(
     -   24,3))+(EijR(8,3)+EijR(14,3)-2*EijR(18,3))*P(9)-(EijR(2,1)+E
     -   ijR(5,2))*P(14)-(EijR(4,1)+EijR(7,2))*P(60)+EijR(4,3)*P(90)+
     -   (EijR(9,3)+EijR(16,3)-2*EijR(20,3))*P(128)-EijR(2,3)*P(152)+
     -   EijR(4,2)*P(464)-EijR(9,2)*P(465)+EijR(10,3)*P(466)-EijR(15,
     -   3)*P(467))
       FI(377) = -4*(-Dij2345I(1,2)-Dij2345I(3,2)+2*Dij2345I(5,2)-s12*
     -   EijI(2,2)+p3sq*EijI(8,2)-p3sq*EijI(10,2)-6*(EijI(22,3)-EijI(
     -   24,3))+(EijI(8,3)+EijI(14,3)-2*EijI(18,3))*P(9)-(EijI(2,1)+E
     -   ijI(5,2))*P(14)-(EijI(4,1)+EijI(7,2))*P(60)+EijI(4,3)*P(90)+
     -   (EijI(9,3)+EijI(16,3)-2*EijI(20,3))*P(128)-EijI(2,3)*P(152)+
     -   EijI(4,2)*P(464)-EijI(9,2)*P(465)+EijI(10,3)*P(466)-EijI(15,
     -   3)*P(467))
       F(377)=DCMPLX(FR(377),FI(377))
       P(468) = p2sq-2*s12-s23+s34+4*s45
       P(469) = 2*P(9)-3*P(13)
       P(470) = p2sq+3*p3sq-s34-s45-2*P(48)
       FR(378) = -4*(-Dij2345R(1,1)+Dij2345R(3,1)-Dij2345R(3,2)+Dij234
     -   5R(5,2)-s12*EijR(2,1)+p3sq*EijR(3,1)-4*(EijR(11,2)-EijR(24,3
     -   ))-EijR(4,1)*P(3)+(EijR(5,2)+EijR(14,3)-EijR(18,3))*P(9)+Eij
     -   R(4,3)*P(90)+(EijR(8,2)+EijR(16,3)-EijR(20,3))*P(128)+(-EijR
     -   (2,2)+EijR(10,3))*P(152)+EijR(10,2)*P(460)-EijR(15,3)*P(462)
     -   -EijR(4,2)*P(468)+EijR(7,2)*P(469)-EijR(9,2)*P(470))
       FI(378) = -4*(-Dij2345I(1,1)+Dij2345I(3,1)-Dij2345I(3,2)+Dij234
     -   5I(5,2)-s12*EijI(2,1)+p3sq*EijI(3,1)-4*(EijI(11,2)-EijI(24,3
     -   ))-EijI(4,1)*P(3)+(EijI(5,2)+EijI(14,3)-EijI(18,3))*P(9)+Eij
     -   I(4,3)*P(90)+(EijI(8,2)+EijI(16,3)-EijI(20,3))*P(128)+(-EijI
     -   (2,2)+EijI(10,3))*P(152)+EijI(10,2)*P(460)-EijI(15,3)*P(462)
     -   -EijI(4,2)*P(468)+EijI(7,2)*P(469)-EijI(9,2)*P(470))
       F(378)=DCMPLX(FR(378),FI(378))
       P(471) = p2sq-p3sq-s23+5*s45-2*P(55)
       P(472) = p2sq-s12-3*P(13)
       P(473) = 4*p3sq+s12+2*P(286)
       P(474) = s23+s34+s45
       P(475) = -4*p3sq+3*s12+2*P(474)
       P(476) = p2sq+2*p3sq-4*s12-s23+5*s45
       P(477) = s12-s34+2*P(114)
       P(478) = p2sq+s34-2*P(114)
       FR(379) = -4*(2*(D02345R-Dij2345R(1,1)+Dij2345R(3,1))-Dij2345R(
     -   3,2)-Dij2345R(4,2)+Dij2345R(5,2)+Dij2345R(6,2)+s15*(EE0R+Eij
     -   R(1,1))-3*p3sq*EijR(3,2)-18*EijR(11,2)-4*(EijR(23,3)-EijR(24
     -   ,3))+(2*EijR(5,2)+EijR(14,3)+EijR(17,3)-EijR(18,3)-EijR(19,3
     -   ))*P(9)-(EijR(3,1)+3*EijR(6,2))*P(14)-EijR(4,1)*P(60)+EijR(4
     -   ,3)*P(90)+(EijR(12,3)-EijR(13,3))*P(128)+(-2*EijR(2,2)-EijR(
     -   9,3)+EijR(10,3))*P(152)-EijR(15,3)*P(462)-EijR(4,2)*P(471)+E
     -   ijR(7,2)*P(472)+EijR(8,2)*P(473)+EijR(9,2)*P(475)+EijR(10,2)
     -   *P(476)+EijR(16,3)*P(477)+EijR(20,3)*P(478))
       FI(379) = -4*(2*(D02345I-Dij2345I(1,1)+Dij2345I(3,1))-Dij2345I(
     -   3,2)-Dij2345I(4,2)+Dij2345I(5,2)+Dij2345I(6,2)+s15*(EE0I+Eij
     -   I(1,1))-3*p3sq*EijI(3,2)-18*EijI(11,2)-4*(EijI(23,3)-EijI(24
     -   ,3))+(2*EijI(5,2)+EijI(14,3)+EijI(17,3)-EijI(18,3)-EijI(19,3
     -   ))*P(9)-(EijI(3,1)+3*EijI(6,2))*P(14)-EijI(4,1)*P(60)+EijI(4
     -   ,3)*P(90)+(EijI(12,3)-EijI(13,3))*P(128)+(-2*EijI(2,2)-EijI(
     -   9,3)+EijI(10,3))*P(152)-EijI(15,3)*P(462)-EijI(4,2)*P(471)+E
     -   ijI(7,2)*P(472)+EijI(8,2)*P(473)+EijI(9,2)*P(475)+EijI(10,2)
     -   *P(476)+EijI(16,3)*P(477)+EijI(20,3)*P(478))
       F(379)=DCMPLX(FR(379),FI(379))
       P(479) = p2sq-3*p3sq-p4sq-s23+s34-2*P(7)
       P(480) = p4sq+2*s45
       P(481) = s15+3*P(14)
       P(482) = s15+2*P(9)-3*P(13)
       P(483) = 3*p3sq-s45+2*P(35)
       P(484) = p2sq-3*p3sq-s12-s23-2*P(79)
       FR(380) = 4*(D02345R-Dij2345R(1,1)-Dij2345R(2,1)+2*Dij2345R(3,1
     -   )-12*EijR(11,2)+4*(EijR(23,3)-EijR(24,3))+EijR(5,2)*P(9)+(-E
     -   ijR(3,1)+EijR(4,1))*P(130)-EijR(2,2)*P(152)-EijR(9,2)*P(470)
     -   +EijR(3,2)*P(479)-EijR(4,2)*P(480)-EijR(6,2)*P(481)+EijR(7,2
     -   )*P(482)+EijR(8,2)*P(483)-EijR(10,2)*P(484))
       FI(380) = 4*(D02345I-Dij2345I(1,1)-Dij2345I(2,1)+2*Dij2345I(3,1
     -   )-12*EijI(11,2)+4*(EijI(23,3)-EijI(24,3))+EijI(5,2)*P(9)+(-E
     -   ijI(3,1)+EijI(4,1))*P(130)-EijI(2,2)*P(152)-EijI(9,2)*P(470)
     -   +EijI(3,2)*P(479)-EijI(4,2)*P(480)-EijI(6,2)*P(481)+EijI(7,2
     -   )*P(482)+EijI(8,2)*P(483)-EijI(10,2)*P(484))
       F(380)=DCMPLX(FR(380),FI(380))
       P(485) = s15-s23+s45+2*P(9)
       P(486) = p2sq-p4sq-s23+s34-3*P(7)
       P(487) = p2sq+p4sq
       FR(381) = -4*(D02345R+s15*(EE0R+EijR(1,1))-p2sq*EijR(2,2)-p4sq*
     -   EijR(4,2)-6*EijR(11,2)-4*(EijR(22,3)-EijR(24,3))-EijR(4,1)*P
     -   (130)-EijR(3,1)*P(169)+EijR(2,1)*P(390)+(EijR(5,2)-EijR(7,2)
     -   )*P(485)+(-EijR(8,2)+EijR(10,2))*P(486)+EijR(9,2)*P(487))
       FI(381) = -4*(D02345I+s15*(EE0I+EijI(1,1))-p2sq*EijI(2,2)-p4sq*
     -   EijI(4,2)-6*EijI(11,2)-4*(EijI(22,3)-EijI(24,3))-EijI(4,1)*P
     -   (130)-EijI(3,1)*P(169)+EijI(2,1)*P(390)+(EijI(5,2)-EijI(7,2)
     -   )*P(485)+(-EijI(8,2)+EijI(10,2))*P(486)+EijI(9,2)*P(487))
       F(381)=DCMPLX(FR(381),FI(381))
       P(488) = 2*p4sq+s45
       P(489) = 3*p2sq+p3sq-s12-s23
       P(490) = p2sq+3*p3sq-p4sq-s23+s34-2*P(7)
       P(491) = p2sq-p3sq-p4sq-s12+s34
       P(492) = s15+s45+2*P(491)
       P(493) = p2sq-3*p3sq+s12+s23+2*P(58)
       P(494) = p2sq+p3sq-p4sq-s23+s34-4*P(7)
       FR(382) = -4*(Dij2345R(3,1)+s15*(EE0R+2*EijR(1,1)+EijR(1,2))-10
     -   *EijR(11,2)-4*(EijR(21,3)-EijR(24,3))+EijR(2,1)*P(2)-2*(Dij2
     -   345R(1,1)-EijR(8,2)*P(128)+EijR(2,2)*P(152))-EijR(3,1)*P(160
     -   )-EijR(4,1)*P(388)-EijR(4,2)*P(488)+EijR(5,2)*P(489)-EijR(6,
     -   2)*P(490)-EijR(7,2)*P(492)+EijR(9,2)*P(493)+EijR(10,2)*P(494
     -   ))
       FI(382) = -4*(Dij2345I(3,1)+s15*(EE0I+2*EijI(1,1)+EijI(1,2))-10
     -   *EijI(11,2)-4*(EijI(21,3)-EijI(24,3))+EijI(2,1)*P(2)-2*(Dij2
     -   345I(1,1)-EijI(8,2)*P(128)+EijI(2,2)*P(152))-EijI(3,1)*P(160
     -   )-EijI(4,1)*P(388)-EijI(4,2)*P(488)+EijI(5,2)*P(489)-EijI(6,
     -   2)*P(490)-EijI(7,2)*P(492)+EijI(9,2)*P(493)+EijI(10,2)*P(494
     -   ))
       F(382)=DCMPLX(FR(382),FI(382))
       P(495) = p2sq-2*p3sq-p4sq-s12-s23+s34+s45
       FR(383) = -4*(Dij2345R(1,1)-2*Dij2345R(3,1)+s12*EijR(2,1)-p3sq*
     -   EijR(3,1)+4*(EijR(11,2)+EijR(24,3))-EijR(5,2)*P(9)-EijR(8,2)
     -   *P(128)+EijR(2,2)*P(152)-EijR(4,1)*P(187)+EijR(9,2)*P(470)+E
     -   ijR(4,2)*P(480)-EijR(7,2)*P(482)+EijR(10,2)*P(495))
       FI(383) = -4*(Dij2345I(1,1)-2*Dij2345I(3,1)+s12*EijI(2,1)-p3sq*
     -   EijI(3,1)+4*(EijI(11,2)+EijI(24,3))-EijI(5,2)*P(9)-EijI(8,2)
     -   *P(128)+EijI(2,2)*P(152)-EijI(4,1)*P(187)+EijI(9,2)*P(470)+E
     -   ijI(4,2)*P(480)-EijI(7,2)*P(482)+EijI(10,2)*P(495))
       F(383)=DCMPLX(FR(383),FI(383))
       P(496) = 2*p2sq+p3sq+s12-s23
       P(497) = 2*p2sq+p3sq+p4sq+s12-s23-s34-s45
       P(498) = p2sq-p3sq-p4sq-s12+s34+s45
       P(499) = s23-2*P(498)
       P(500) = p2sq-p3sq-p4sq-s12
       P(501) = -s23+3*P(58)+2*P(500)
       P(502) = p4sq-s23+2*s45
       P(503) = p3sq-p4sq-s12-s23-s45
       P(504) = p3sq-p4sq-s12-2*s45
       FR(384) = -4*(-Dij2345R(3,2)+Dij2345R(5,2)+p3sq*EijR(3,1)+p4sq*
     -   EijR(4,2)+3*p3sq*EijR(6,2)+4*(EijR(11,2)-EijR(21,3))+2*(Dij2
     -   345R(1,1)-Dij2345R(3,1)-EijR(22,3))+6*EijR(24,3)+EijR(10,2)*
     -   P(7)+(EijR(5,3)-EijR(7,3))*P(13)+(EijR(8,3)-EijR(10,3))*P(26
     -   )+(-EijR(16,3)-EijR(17,3)+EijR(19,3)+EijR(20,3))*P(38)-EijR(
     -   8,2)*P(208)+EijR(4,3)*P(213)-EijR(2,1)*P(371)+EijR(4,1)*P(37
     -   2)+EijR(2,2)*P(496)-EijR(5,2)*P(497)-EijR(7,2)*P(499)-EijR(9
     -   ,2)*P(501)-EijR(14,3)*P(502)+EijR(15,3)*P(503)-EijR(18,3)*P(
     -   504))
       FI(384) = -4*(-Dij2345I(3,2)+Dij2345I(5,2)+p3sq*EijI(3,1)+p4sq*
     -   EijI(4,2)+3*p3sq*EijI(6,2)+4*(EijI(11,2)-EijI(21,3))+2*(Dij2
     -   345I(1,1)-Dij2345I(3,1)-EijI(22,3))+6*EijI(24,3)+EijI(10,2)*
     -   P(7)+(EijI(5,3)-EijI(7,3))*P(13)+(EijI(8,3)-EijI(10,3))*P(26
     -   )+(-EijI(16,3)-EijI(17,3)+EijI(19,3)+EijI(20,3))*P(38)-EijI(
     -   8,2)*P(208)+EijI(4,3)*P(213)-EijI(2,1)*P(371)+EijI(4,1)*P(37
     -   2)+EijI(2,2)*P(496)-EijI(5,2)*P(497)-EijI(7,2)*P(499)-EijI(9
     -   ,2)*P(501)-EijI(14,3)*P(502)+EijI(15,3)*P(503)-EijI(18,3)*P(
     -   504))
       F(384)=DCMPLX(FR(384),FI(384))
       P(505) = p3sq-p4sq-3*s12+s34+2*P(30)
       P(506) = -p4sq-s45+2*P(26)
       P(507) = p3sq-s12-s23-2*P(213)
       FR(385) = 4*(Dij2345R(1,2)+Dij2345R(3,2)+s45*EijR(4,2)+p3sq*Eij
     -   R(8,2)-p3sq*EijR(10,2)+6*(EijR(22,3)-EijR(24,3))-(EijR(8,3)+
     -   EijR(14,3)-2*EijR(18,3))*P(13)-2*(Dij2345R(5,2)-(EijR(5,2)-E
     -   ijR(7,2))*P(14))-EijR(2,3)*P(26)+(EijR(9,3)+EijR(16,3)-2*Eij
     -   R(20,3))*P(38)+EijR(9,2)*P(52)+EijR(2,1)*P(86)-EijR(4,3)*P(2
     -   13)+EijR(4,1)*P(259)-EijR(2,2)*P(505)+EijR(10,3)*P(506)-EijR
     -   (15,3)*P(507))
       FI(385) = 4*(Dij2345I(1,2)+Dij2345I(3,2)+s45*EijI(4,2)+p3sq*Eij
     -   I(8,2)-p3sq*EijI(10,2)+6*(EijI(22,3)-EijI(24,3))-(EijI(8,3)+
     -   EijI(14,3)-2*EijI(18,3))*P(13)-2*(Dij2345I(5,2)-(EijI(5,2)-E
     -   ijI(7,2))*P(14))-EijI(2,3)*P(26)+(EijI(9,3)+EijI(16,3)-2*Eij
     -   I(20,3))*P(38)+EijI(9,2)*P(52)+EijI(2,1)*P(86)-EijI(4,3)*P(2
     -   13)+EijI(4,1)*P(259)-EijI(2,2)*P(505)+EijI(10,3)*P(506)-EijI
     -   (15,3)*P(507))
       F(385)=DCMPLX(FR(385),FI(385))
       P(508) = -p4sq-2*s12+s34+s45+3*P(16)
       P(509) = 2*p2sq-p4sq-4*s12+s45+3*P(18)
       P(510) = 3*p3sq+p4sq-s12-s34+2*s45
       P(511) = p3sq+p4sq-s12+2*s45
       P(512) = p4sq+s23+2*s45
       FR(386) = -4*(D02345R-Dij2345R(1,1)+Dij2345R(3,1)-Dij2345R(3,2)
     -   -Dij2345R(4,2)+Dij2345R(5,2)+Dij2345R(6,2)-2*(6*EijR(11,2)+E
     -   ijR(22,3)+2*EijR(23,3)-3*EijR(24,3))+EijR(5,2)*P(9)+(EijR(14
     -   ,3)+EijR(17,3)-EijR(18,3)-EijR(19,3))*P(13)-3*(p3sq*EijR(3,2
     -   )+EijR(6,2)*P(14))+(EijR(9,3)-EijR(10,3))*P(26)+(-EijR(12,3)
     -   +EijR(13,3))*P(38)-EijR(3,1)*P(86)-EijR(2,2)*P(152)+EijR(4,3
     -   )*P(213)-EijR(4,1)*P(259)+EijR(7,2)*P(469)-EijR(4,2)*P(480)+
     -   EijR(15,3)*P(503)+EijR(8,2)*P(508)-EijR(9,2)*P(509)+EijR(10,
     -   2)*P(510)-EijR(16,3)*P(511)+EijR(20,3)*P(512))
       FI(386) = -4*(D02345I-Dij2345I(1,1)+Dij2345I(3,1)-Dij2345I(3,2)
     -   -Dij2345I(4,2)+Dij2345I(5,2)+Dij2345I(6,2)-2*(6*EijI(11,2)+E
     -   ijI(22,3)+2*EijI(23,3)-3*EijI(24,3))+EijI(5,2)*P(9)+(EijI(14
     -   ,3)+EijI(17,3)-EijI(18,3)-EijI(19,3))*P(13)-3*(p3sq*EijI(3,2
     -   )+EijI(6,2)*P(14))+(EijI(9,3)-EijI(10,3))*P(26)+(-EijI(12,3)
     -   +EijI(13,3))*P(38)-EijI(3,1)*P(86)-EijI(2,2)*P(152)+EijI(4,3
     -   )*P(213)-EijI(4,1)*P(259)+EijI(7,2)*P(469)-EijI(4,2)*P(480)+
     -   EijI(15,3)*P(503)+EijI(8,2)*P(508)-EijI(9,2)*P(509)+EijI(10,
     -   2)*P(510)-EijI(16,3)*P(511)+EijI(20,3)*P(512))
       F(386)=DCMPLX(FR(386),FI(386))
       FR(387) = -4*(-Dij2345R(1,1)+Dij2345R(3,1)-Dij2345R(3,2)+Dij234
     -   5R(5,2)-s12*EijR(2,1)+p3sq*EijR(3,1)-4*(EijR(11,2)-EijR(24,3
     -   ))-EijR(4,1)*P(3)+EijR(5,2)*P(9)+(EijR(14,3)-EijR(18,3))*P(1
     -   3)-EijR(10,3)*P(26)+(-EijR(16,3)+EijR(20,3))*P(38)+EijR(8,2)
     -   *P(128)-EijR(2,2)*P(152)+EijR(4,3)*P(213)+EijR(10,2)*P(460)+
     -   EijR(7,2)*P(469)-EijR(4,2)*P(480)+EijR(15,3)*P(503)-EijR(9,2
     -   )*P(509))
       FI(387) = -4*(-Dij2345I(1,1)+Dij2345I(3,1)-Dij2345I(3,2)+Dij234
     -   5I(5,2)-s12*EijI(2,1)+p3sq*EijI(3,1)-4*(EijI(11,2)-EijI(24,3
     -   ))-EijI(4,1)*P(3)+EijI(5,2)*P(9)+(EijI(14,3)-EijI(18,3))*P(1
     -   3)-EijI(10,3)*P(26)+(-EijI(16,3)+EijI(20,3))*P(38)+EijI(8,2)
     -   *P(128)-EijI(2,2)*P(152)+EijI(4,3)*P(213)+EijI(10,2)*P(460)+
     -   EijI(7,2)*P(469)-EijI(4,2)*P(480)+EijI(15,3)*P(503)-EijI(9,2
     -   )*P(509))
       F(387)=DCMPLX(FR(387),FI(387))
       FR(388) = 8*(EijR(3,2)-EijR(4,2)+EijR(4,3)-EijR(11,3)+EijR(13,3
     -   )-EijR(14,3)-2*(EijR(6,2)-EijR(7,2)+EijR(16,3)-EijR(19,3)))
       FI(388) = 8*(EijI(3,2)-EijI(4,2)+EijI(4,3)-EijI(11,3)+EijI(13,3
     -   )-EijI(14,3)-2*(EijI(6,2)-EijI(7,2)+EijI(16,3)-EijI(19,3)))
       F(388)=DCMPLX(FR(388),FI(388))
       FR(389) = -8*(2*EijR(4,2)-EijR(4,3)-3*EijR(7,2)+EijR(10,2)+EijR
     -   (14,3)+EijR(16,3)-EijR(19,3))
       FI(389) = -8*(2*EijI(4,2)-EijI(4,3)-3*EijI(7,2)+EijI(10,2)+EijI
     -   (14,3)+EijI(16,3)-EijI(19,3))
       F(389)=DCMPLX(FR(389),FI(389))
       FR(390) = -8*(-EijR(4,3)+2*(EijR(4,2)-EijR(9,2))+EijR(15,3)+Eij
     -   R(16,3)-EijR(20,3))
       FI(390) = -8*(-EijI(4,3)+2*(EijI(4,2)-EijI(9,2))+EijI(15,3)+Eij
     -   I(16,3)-EijI(20,3))
       F(390)=DCMPLX(FR(390),FI(390))
       FR(391) = 8*(EijR(4,3)+EijR(13,3)-2*(EijR(4,2)-EijR(10,2)+EijR(
     -   16,3)))
       FI(391) = 8*(EijI(4,3)+EijI(13,3)-2*(EijI(4,2)-EijI(10,2)+EijI(
     -   16,3)))
       F(391)=DCMPLX(FR(391),FI(391))
       FR(392) = -8*(EijR(3,1)-EijR(4,1)-EijR(4,2)-EijR(4,3)+EijR(6,2)
     -   +2*(EijR(7,2)-EijR(10,2))+EijR(14,3)+EijR(16,3)-EijR(19,3))
       FI(392) = -8*(EijI(3,1)-EijI(4,1)-EijI(4,2)-EijI(4,3)+EijI(6,2)
     -   +2*(EijI(7,2)-EijI(10,2))+EijI(14,3)+EijI(16,3)-EijI(19,3))
       F(392)=DCMPLX(FR(392),FI(392))
       FR(393) = -8*(EijR(3,1)-EijR(4,1)-EijR(4,3)+EijR(8,2)+EijR(9,2)
     -   -2*EijR(10,2)+EijR(15,3)+EijR(16,3)-EijR(20,3))
       FI(393) = -8*(EijI(3,1)-EijI(4,1)-EijI(4,3)+EijI(8,2)+EijI(9,2)
     -   -2*EijI(10,2)+EijI(15,3)+EijI(16,3)-EijI(20,3))
       F(393)=DCMPLX(FR(393),FI(393))
       FR(394) = EijR(3,1)-EijR(4,1)+EijR(4,3)-EijR(6,3)+EijR(7,3)-Eij
     -   R(16,3)-2*(EijR(14,3)-EijR(19,3))
       FI(394) = EijI(3,1)-EijI(4,1)+EijI(4,3)-EijI(6,3)+EijI(7,3)-Eij
     -   I(16,3)-2*(EijI(14,3)-EijI(19,3))
       F(394)=DCMPLX(FR(394),FI(394))
       FR(395) = EijR(2,1)-EijR(3,1)+EijR(4,2)-EijR(4,3)+EijR(5,2)+Eij
     -   R(6,2)-2*EijR(7,2)-EijR(8,2)+EijR(14,3)+EijR(15,3)+EijR(16,3
     -   )+EijR(17,3)-EijR(18,3)-EijR(19,3)-EijR(20,3)
       FI(395) = EijI(2,1)-EijI(3,1)+EijI(4,2)-EijI(4,3)+EijI(5,2)+Eij
     -   I(6,2)-2*EijI(7,2)-EijI(8,2)+EijI(14,3)+EijI(15,3)+EijI(16,3
     -   )+EijI(17,3)-EijI(18,3)-EijI(19,3)-EijI(20,3)
       F(395)=DCMPLX(FR(395),FI(395))
       FR(396) = EijR(2,1)-EijR(4,1)+EijR(4,3)+EijR(5,2)-EijR(7,2)-Eij
     -   R(8,2)+EijR(10,2)-EijR(14,3)-EijR(15,3)-EijR(16,3)-EijR(17,3
     -   )+EijR(18,3)+EijR(19,3)+EijR(20,3)
       FI(396) = EijI(2,1)-EijI(4,1)+EijI(4,3)+EijI(5,2)-EijI(7,2)-Eij
     -   I(8,2)+EijI(10,2)-EijI(14,3)-EijI(15,3)-EijI(16,3)-EijI(17,3
     -   )+EijI(18,3)+EijI(19,3)+EijI(20,3)
       F(396)=DCMPLX(FR(396),FI(396))
       FR(397) = EijR(4,2)-EijR(4,3)+EijR(8,2)-EijR(9,2)+EijR(9,3)-Eij
     -   R(10,2)-EijR(10,3)+EijR(16,3)+2*(EijR(15,3)-EijR(20,3))
       FI(397) = EijI(4,2)-EijI(4,3)+EijI(8,2)-EijI(9,2)+EijI(9,3)-Eij
     -   I(10,2)-EijI(10,3)+EijI(16,3)+2*(EijI(15,3)-EijI(20,3))
       F(397)=DCMPLX(FR(397),FI(397))
       FR(398) = EijR(4,2)-EijR(4,3)+EijR(8,2)-EijR(9,2)-EijR(10,2)+Ei
     -   jR(12,3)-EijR(13,3)+EijR(15,3)+2*(EijR(16,3)-EijR(20,3))
       FI(398) = EijI(4,2)-EijI(4,3)+EijI(8,2)-EijI(9,2)-EijI(10,2)+Ei
     -   jI(12,3)-EijI(13,3)+EijI(15,3)+2*(EijI(16,3)-EijI(20,3))
       F(398)=DCMPLX(FR(398),FI(398))
       FR(399) = EijR(3,1)-EijR(3,2)-EijR(4,1)+EijR(4,3)+EijR(6,2)-Eij
     -   R(7,2)+EijR(10,2)-EijR(11,3)+EijR(13,3)-EijR(14,3)-2*(EijR(1
     -   6,3)-EijR(19,3))
       FI(399) = EijI(3,1)-EijI(3,2)-EijI(4,1)+EijI(4,3)+EijI(6,2)-Eij
     -   I(7,2)+EijI(10,2)-EijI(11,3)+EijI(13,3)-EijI(14,3)-2*(EijI(1
     -   6,3)-EijI(19,3))
       F(399)=DCMPLX(FR(399),FI(399))
       FR(400) = EijR(3,2)+EijR(4,2)-EijR(4,3)+EijR(12,3)-EijR(13,3)+E
     -   ijR(15,3)-2*(EijR(10,2)-EijR(16,3)+EijR(20,3))
       FI(400) = EijI(3,2)+EijI(4,2)-EijI(4,3)+EijI(12,3)-EijI(13,3)+E
     -   ijI(15,3)-2*(EijI(10,2)-EijI(16,3)+EijI(20,3))
       F(400)=DCMPLX(FR(400),FI(400))
       FR(401) = EijR(3,2)+EijR(3,3)+EijR(4,2)-EijR(4,3)-2*EijR(10,2)-
     -   3*(EijR(13,3)-EijR(16,3))
       FI(401) = EijI(3,2)+EijI(3,3)+EijI(4,2)-EijI(4,3)-2*EijI(10,2)-
     -   3*(EijI(13,3)-EijI(16,3))
       F(401)=DCMPLX(FR(401),FI(401))
       FR(402) = EijR(3,1)+EijR(3,2)-EijR(4,1)-EijR(4,3)-EijR(10,2)-Ei
     -   jR(13,3)+2*EijR(16,3)
       FI(402) = EijI(3,1)+EijI(3,2)-EijI(4,1)-EijI(4,3)-EijI(10,2)-Ei
     -   jI(13,3)+2*EijI(16,3)
       F(402)=DCMPLX(FR(402),FI(402))
       FR(403) = EijR(4,2)-EijR(4,3)-EijR(10,2)+EijR(16,3)
       FI(403) = EijI(4,2)-EijI(4,3)-EijI(10,2)+EijI(16,3)
       F(403)=DCMPLX(FR(403),FI(403))
       FR(404) = -4*(s15*(EE0R+EijR(1,1))+2*(D02345R-Dij2345R(1,1)+Dij
     -   2345R(3,1)-2*p3sq*EijR(3,2)-9*EijR(11,2)+EijR(23,3)-EijR(24,
     -   3)+EijR(5,2)*P(9)-EijR(7,2)*P(13)-EijR(6,2)*P(14)-EijR(2,2)*
     -   P(152)-EijR(4,2)*P(213)+EijR(10,2)*P(379)+EijR(9,2)*P(384)+E
     -   ijR(8,2)*P(385)))
       FI(404) = -4*(s15*(EE0I+EijI(1,1))+2*(D02345I-Dij2345I(1,1)+Dij
     -   2345I(3,1)-2*p3sq*EijI(3,2)-9*EijI(11,2)+EijI(23,3)-EijI(24,
     -   3)+EijI(5,2)*P(9)-EijI(7,2)*P(13)-EijI(6,2)*P(14)-EijI(2,2)*
     -   P(152)-EijI(4,2)*P(213)+EijI(10,2)*P(379)+EijI(9,2)*P(384)+E
     -   ijI(8,2)*P(385)))
       F(404)=DCMPLX(FR(404),FI(404))
       FR(405) = -8*(EijR(2,2)-EijR(4,2)+EijR(4,3)-EijR(8,3)+EijR(10,3
     -   )-EijR(14,3)-2*(EijR(5,2)-EijR(7,2)+EijR(15,3)-EijR(18,3)))
       FI(405) = -8*(EijI(2,2)-EijI(4,2)+EijI(4,3)-EijI(8,3)+EijI(10,3
     -   )-EijI(14,3)-2*(EijI(5,2)-EijI(7,2)+EijI(15,3)-EijI(18,3)))
       F(405)=DCMPLX(FR(405),FI(405))
       FR(406) = 8*(EijR(4,2)-EijR(4,3)+2*(EijR(6,2)-EijR(7,2))-EijR(8
     -   ,2)+EijR(9,2)-EijR(10,2)+EijR(14,3)+EijR(15,3)+EijR(16,3)+Ei
     -   jR(17,3)-EijR(18,3)-EijR(19,3)-EijR(20,3))
       FI(406) = 8*(EijI(4,2)-EijI(4,3)+2*(EijI(6,2)-EijI(7,2))-EijI(8
     -   ,2)+EijI(9,2)-EijI(10,2)+EijI(14,3)+EijI(15,3)+EijI(16,3)+Ei
     -   jI(17,3)-EijI(18,3)-EijI(19,3)-EijI(20,3))
       F(406)=DCMPLX(FR(406),FI(406))
       FR(407) = 8*(2*EijR(4,2)-EijR(4,3)-3*EijR(7,2)+EijR(9,2)+EijR(1
     -   4,3)+EijR(15,3)-EijR(18,3))
       FI(407) = 8*(2*EijI(4,2)-EijI(4,3)-3*EijI(7,2)+EijI(9,2)+EijI(1
     -   4,3)+EijI(15,3)-EijI(18,3))
       F(407)=DCMPLX(FR(407),FI(407))
       FR(408) = 8*(-EijR(4,3)-EijR(10,3)+2*(EijR(4,2)-EijR(9,2)+EijR(
     -   15,3)))
       FI(408) = 8*(-EijI(4,3)-EijI(10,3)+2*(EijI(4,2)-EijI(9,2)+EijI(
     -   15,3)))
       F(408)=DCMPLX(FR(408),FI(408))
       FR(409) = 8*(-EijR(4,3)+2*(EijR(4,2)-EijR(10,2))+EijR(15,3)+Eij
     -   R(16,3)-EijR(20,3))
       FI(409) = 8*(-EijI(4,3)+2*(EijI(4,2)-EijI(10,2))+EijI(15,3)+Eij
     -   I(16,3)-EijI(20,3))
       F(409)=DCMPLX(FR(409),FI(409))
       FR(410) = 8*(EijR(2,1)-EijR(4,1)-EijR(4,2)-EijR(4,3)+EijR(5,2)+
     -   2*(EijR(7,2)-EijR(9,2))+EijR(14,3)+EijR(15,3)-EijR(18,3))
       FI(410) = 8*(EijI(2,1)-EijI(4,1)-EijI(4,2)-EijI(4,3)+EijI(5,2)+
     -   2*(EijI(7,2)-EijI(9,2))+EijI(14,3)+EijI(15,3)-EijI(18,3))
       F(410)=DCMPLX(FR(410),FI(410))
       FR(411) = EijR(2,1)-EijR(4,1)+EijR(4,3)-EijR(5,3)+EijR(7,3)-Eij
     -   R(15,3)-2*(EijR(14,3)-EijR(18,3))
       FI(411) = EijI(2,1)-EijI(4,1)+EijI(4,3)-EijI(5,3)+EijI(7,3)-Eij
     -   I(15,3)-2*(EijI(14,3)-EijI(18,3))
       F(411)=DCMPLX(FR(411),FI(411))
       FR(412) = EijR(2,1)-EijR(2,2)-EijR(4,1)+EijR(4,3)+EijR(5,2)-Eij
     -   R(7,2)-EijR(8,3)+EijR(9,2)+EijR(10,3)-EijR(14,3)-2*(EijR(15,
     -   3)-EijR(18,3))
       FI(412) = EijI(2,1)-EijI(2,2)-EijI(4,1)+EijI(4,3)+EijI(5,2)-Eij
     -   I(7,2)-EijI(8,3)+EijI(9,2)+EijI(10,3)-EijI(14,3)-2*(EijI(15,
     -   3)-EijI(18,3))
       F(412)=DCMPLX(FR(412),FI(412))
       FR(413) = EijR(2,2)+EijR(2,3)+EijR(4,2)-EijR(4,3)-2*EijR(9,2)-3
     -   *(EijR(10,3)-EijR(15,3))
       FI(413) = EijI(2,2)+EijI(2,3)+EijI(4,2)-EijI(4,3)-2*EijI(9,2)-3
     -   *(EijI(10,3)-EijI(15,3))
       F(413)=DCMPLX(FR(413),FI(413))
       FR(414) = EijR(3,1)-EijR(4,1)+EijR(4,3)+EijR(6,2)-EijR(7,2)-Eij
     -   R(8,2)+EijR(9,2)-EijR(14,3)-EijR(15,3)-EijR(16,3)-EijR(17,3)
     -   +EijR(18,3)+EijR(19,3)+EijR(20,3)
       FI(414) = EijI(3,1)-EijI(4,1)+EijI(4,3)+EijI(6,2)-EijI(7,2)-Eij
     -   I(8,2)+EijI(9,2)-EijI(14,3)-EijI(15,3)-EijI(16,3)-EijI(17,3)
     -   +EijI(18,3)+EijI(19,3)+EijI(20,3)
       F(414)=DCMPLX(FR(414),FI(414))
       FR(415) = EijR(2,1)+EijR(2,2)-EijR(4,1)-EijR(4,3)-EijR(9,2)-Eij
     -   R(10,3)+2*EijR(15,3)
       FI(415) = EijI(2,1)+EijI(2,2)-EijI(4,1)-EijI(4,3)-EijI(9,2)-Eij
     -   I(10,3)+2*EijI(15,3)
       F(415)=DCMPLX(FR(415),FI(415))
       FR(416) = EijR(3,1)-EijR(4,1)-EijR(4,3)+EijR(8,2)-EijR(9,2)+Eij
     -   R(15,3)+EijR(16,3)-EijR(20,3)
       FI(416) = EijI(3,1)-EijI(4,1)-EijI(4,3)+EijI(8,2)-EijI(9,2)+Eij
     -   I(15,3)+EijI(16,3)-EijI(20,3)
       F(416)=DCMPLX(FR(416),FI(416))
       FR(417) = EijR(4,2)-EijR(4,3)-EijR(9,2)+EijR(15,3)
       FI(417) = EijI(4,2)-EijI(4,3)-EijI(9,2)+EijI(15,3)
       F(417)=DCMPLX(FR(417),FI(417))
       FR(418) = 4*(s15*(EE0R+EijR(1,1))+2*(D02345R-Dij2345R(1,1)+Dij2
     -   345R(3,1)-2*p3sq*EijR(3,2)-9*EijR(11,2)+EijR(22,3)-EijR(24,3
     -   )+EijR(5,2)*P(9)-EijR(7,2)*P(13)+EijR(6,2)*P(60)-EijR(2,2)*P
     -   (152)-EijR(4,2)*P(213)+EijR(10,2)*P(379)+EijR(9,2)*P(384)+Ei
     -   jR(8,2)*P(385)))
       FI(418) = 4*(s15*(EE0I+EijI(1,1))+2*(D02345I-Dij2345I(1,1)+Dij2
     -   345I(3,1)-2*p3sq*EijI(3,2)-9*EijI(11,2)+EijI(22,3)-EijI(24,3
     -   )+EijI(5,2)*P(9)-EijI(7,2)*P(13)+EijI(6,2)*P(60)-EijI(2,2)*P
     -   (152)-EijI(4,2)*P(213)+EijI(10,2)*P(379)+EijI(9,2)*P(384)+Ei
     -   jI(8,2)*P(385)))
       F(418)=DCMPLX(FR(418),FI(418))
       FR(419) = 8*(EijR(2,2)-2*(EijR(5,2)-EijR(7,2))+EijR(8,2)-EijR(8
     -   ,3)-EijR(9,2)-EijR(10,2)+EijR(10,3)-EijR(15,3)+EijR(16,3)+Ei
     -   jR(17,3)+EijR(18,3)-EijR(19,3)-EijR(20,3))
       FI(419) = 8*(EijI(2,2)-2*(EijI(5,2)-EijI(7,2))+EijI(8,2)-EijI(8
     -   ,3)-EijI(9,2)-EijI(10,2)+EijI(10,3)-EijI(15,3)+EijI(16,3)+Ei
     -   jI(17,3)+EijI(18,3)-EijI(19,3)-EijI(20,3))
       F(419)=DCMPLX(FR(419),FI(419))
       FR(420) = 8*(EijR(3,2)-2*(EijR(6,2)-EijR(7,2))+EijR(8,2)-EijR(9
     -   ,2)-EijR(10,2)+EijR(11,3)-EijR(13,3)-EijR(15,3)+EijR(16,3)-E
     -   ijR(17,3)+EijR(18,3)-EijR(19,3)+EijR(20,3))
       FI(420) = 8*(EijI(3,2)-2*(EijI(6,2)-EijI(7,2))+EijI(8,2)-EijI(9
     -   ,2)-EijI(10,2)+EijI(11,3)-EijI(13,3)-EijI(15,3)+EijI(16,3)-E
     -   ijI(17,3)+EijI(18,3)-EijI(19,3)+EijI(20,3))
       F(420)=DCMPLX(FR(420),FI(420))
       FR(421) = 8*(EijR(3,1)-EijR(4,1)+EijR(6,2)-EijR(9,2)+2*(EijR(7,
     -   2)-EijR(10,2))-EijR(15,3)+EijR(16,3)+EijR(18,3)-EijR(19,3))
       FI(421) = 8*(EijI(3,1)-EijI(4,1)+EijI(6,2)-EijI(9,2)+2*(EijI(7,
     -   2)-EijI(10,2))-EijI(15,3)+EijI(16,3)+EijI(18,3)-EijI(19,3))
       F(421)=DCMPLX(FR(421),FI(421))
       FR(422) = 8*(EijR(3,1)-EijR(4,1)+EijR(8,2)+EijR(9,2)-2*EijR(10,
     -   2)+EijR(10,3)-EijR(15,3)+EijR(16,3)-EijR(20,3))
       FI(422) = 8*(EijI(3,1)-EijI(4,1)+EijI(8,2)+EijI(9,2)-2*EijI(10,
     -   2)+EijI(10,3)-EijI(15,3)+EijI(16,3)-EijI(20,3))
       F(422)=DCMPLX(FR(422),FI(422))
       FR(423) = -8*(EijR(3,2)-2*(EijR(6,2)-EijR(7,2))+EijR(8,2)-EijR(
     -   9,2)-EijR(10,2)-EijR(11,3)+EijR(13,3)+EijR(15,3)-EijR(16,3)+
     -   EijR(17,3)-EijR(18,3)+EijR(19,3)-EijR(20,3))
       FI(423) = -8*(EijI(3,2)-2*(EijI(6,2)-EijI(7,2))+EijI(8,2)-EijI(
     -   9,2)-EijI(10,2)-EijI(11,3)+EijI(13,3)+EijI(15,3)-EijI(16,3)+
     -   EijI(17,3)-EijI(18,3)+EijI(19,3)-EijI(20,3))
       F(423)=DCMPLX(FR(423),FI(423))
       FR(424) = -8*(EijR(2,1)-EijR(4,1)+EijR(5,2)+2*(EijR(7,2)-EijR(9
     -   ,2))-EijR(10,2)+EijR(15,3)-EijR(16,3)-EijR(18,3)+EijR(19,3))
       FI(424) = -8*(EijI(2,1)-EijI(4,1)+EijI(5,2)+2*(EijI(7,2)-EijI(9
     -   ,2))-EijI(10,2)+EijI(15,3)-EijI(16,3)-EijI(18,3)+EijI(19,3))
       F(424)=DCMPLX(FR(424),FI(424))
       FR(425) = EijR(2,1)-EijR(3,1)-EijR(5,3)+EijR(6,3)-EijR(15,3)+Ei
     -   jR(16,3)+2*(EijR(18,3)-EijR(19,3))
       FI(425) = EijI(2,1)-EijI(3,1)-EijI(5,3)+EijI(6,3)-EijI(15,3)+Ei
     -   jI(16,3)+2*(EijI(18,3)-EijI(19,3))
       F(425)=DCMPLX(FR(425),FI(425))
       FR(426) = EijR(2,1)-EijR(2,2)-EijR(3,1)+EijR(5,2)+EijR(6,2)-Eij
     -   R(8,2)-EijR(8,3)-2*(EijR(7,2)-EijR(9,2))+EijR(10,3)-EijR(15,
     -   3)+EijR(16,3)+EijR(17,3)+EijR(18,3)-EijR(19,3)-EijR(20,3)
       FI(426) = EijI(2,1)-EijI(2,2)-EijI(3,1)+EijI(5,2)+EijI(6,2)-Eij
     -   I(8,2)-EijI(8,3)-2*(EijI(7,2)-EijI(9,2))+EijI(10,3)-EijI(15,
     -   3)+EijI(16,3)+EijI(17,3)+EijI(18,3)-EijI(19,3)-EijI(20,3)
       F(426)=DCMPLX(FR(426),FI(426))
       FR(427) = EijR(2,2)+EijR(2,3)-EijR(8,2)-EijR(9,2)-EijR(9,3)+Eij
     -   R(10,2)+EijR(15,3)-EijR(16,3)-2*(EijR(10,3)-EijR(20,3))
       FI(427) = EijI(2,2)+EijI(2,3)-EijI(8,2)-EijI(9,2)-EijI(9,3)+Eij
     -   I(10,2)+EijI(15,3)-EijI(16,3)-2*(EijI(10,3)-EijI(20,3))
       F(427)=DCMPLX(FR(427),FI(427))
       FR(428) = EijR(3,2)-EijR(8,2)+EijR(9,2)-EijR(9,3)-EijR(10,2)+Ei
     -   jR(10,3)+EijR(12,3)-EijR(13,3)-EijR(15,3)+EijR(16,3)
       FI(428) = EijI(3,2)-EijI(8,2)+EijI(9,2)-EijI(9,3)-EijI(10,2)+Ei
     -   jI(10,3)+EijI(12,3)-EijI(13,3)-EijI(15,3)+EijI(16,3)
       F(428)=DCMPLX(FR(428),FI(428))
       FR(429) = EijR(9,3)-EijR(10,3)-EijR(12,3)+EijR(13,3)+EijR(15,3)
     -   -EijR(16,3)
       FI(429) = EijI(9,3)-EijI(10,3)-EijI(12,3)+EijI(13,3)+EijI(15,3)
     -   -EijI(16,3)
       F(429)=DCMPLX(FR(429),FI(429))
       FR(430) = -EijR(3,3)+EijR(12,3)+EijR(15,3)-EijR(16,3)+2*(EijR(1
     -   3,3)-EijR(20,3))
       FI(430) = -EijI(3,3)+EijI(12,3)+EijI(15,3)-EijI(16,3)+2*(EijI(1
     -   3,3)-EijI(20,3))
       F(430)=DCMPLX(FR(430),FI(430))
       FR(431) = EijR(3,1)+EijR(3,2)-EijR(4,1)-EijR(10,2)-EijR(13,3)-E
     -   ijR(15,3)+EijR(16,3)+EijR(20,3)
       FI(431) = EijI(3,1)+EijI(3,2)-EijI(4,1)-EijI(10,2)-EijI(13,3)-E
     -   ijI(15,3)+EijI(16,3)+EijI(20,3)
       F(431)=DCMPLX(FR(431),FI(431))
       FR(432) = EijR(2,1)+EijR(2,2)-EijR(4,1)-EijR(9,2)-EijR(10,3)+Ei
     -   jR(15,3)-EijR(16,3)+EijR(20,3)
       FI(432) = EijI(2,1)+EijI(2,2)-EijI(4,1)-EijI(9,2)-EijI(10,3)+Ei
     -   jI(15,3)-EijI(16,3)+EijI(20,3)
       F(432)=DCMPLX(FR(432),FI(432))
       FR(433) = EijR(3,1)-EijR(4,1)+EijR(8,2)-EijR(9,2)+EijR(13,3)+Ei
     -   jR(15,3)-EijR(16,3)-EijR(20,3)
       FI(433) = EijI(3,1)-EijI(4,1)+EijI(8,2)-EijI(9,2)+EijI(13,3)+Ei
     -   jI(15,3)-EijI(16,3)-EijI(20,3)
       F(433)=DCMPLX(FR(433),FI(433))
       FR(434) = -EijR(9,2)+EijR(10,2)+EijR(15,3)-EijR(16,3)
       FI(434) = -EijI(9,2)+EijI(10,2)+EijI(15,3)-EijI(16,3)
       F(434)=DCMPLX(FR(434),FI(434))
       FR(435) = -4*(s15*(EE0R+EijR(1,1))+2*(D02345R-Dij2345R(1,1)+Dij
     -   2345R(3,1)-2*p3sq*EijR(3,2)-9*EijR(11,2)+EijR(22,3)-EijR(23,
     -   3)+EijR(5,2)*P(9)-EijR(7,2)*P(13)-EijR(6,2)*P(14)-EijR(2,2)*
     -   P(152)-EijR(4,2)*P(213)+EijR(10,2)*P(379)+EijR(9,2)*P(384)+E
     -   ijR(8,2)*P(385)))
       FI(435) = -4*(s15*(EE0I+EijI(1,1))+2*(D02345I-Dij2345I(1,1)+Dij
     -   2345I(3,1)-2*p3sq*EijI(3,2)-9*EijI(11,2)+EijI(22,3)-EijI(23,
     -   3)+EijI(5,2)*P(9)-EijI(7,2)*P(13)-EijI(6,2)*P(14)-EijI(2,2)*
     -   P(152)-EijI(4,2)*P(213)+EijI(10,2)*P(379)+EijI(9,2)*P(384)+E
     -   ijI(8,2)*P(385)))
       F(435)=DCMPLX(FR(435),FI(435))
       P(513) = -p2sq+s15+2*P(144)
       P(514) = 2*s12-s34+3*P(114)
       P(515) = s12-2*s34+3*P(114)
       FR(436) = 4*(Dij2345R(2,1)+Dij2345R(2,2)-Dij2345R(3,1)+Dij2345R
     -   (3,2)+s12*EijR(8,2)-s12*EijR(9,2)-2*(Dij2345R(6,2)+EijR(23,3
     -   )-EijR(24,3))+(EijR(3,1)+EijR(6,2)-EijR(11,3)-EijR(14,3)+2*E
     -   ijR(19,3))*P(9)+(EijR(4,1)+EijR(7,2))*P(47)-EijR(4,3)*P(90)-
     -   EijR(3,2)*P(114)-EijR(3,3)*P(128)+(EijR(12,3)+EijR(15,3)-2*E
     -   ijR(20,3))*P(152)+EijR(4,2)*P(403)+EijR(10,2)*P(513)+EijR(13
     -   ,3)*P(514)-EijR(16,3)*P(515))
       FI(436) = 4*(Dij2345I(2,1)+Dij2345I(2,2)-Dij2345I(3,1)+Dij2345I
     -   (3,2)+s12*EijI(8,2)-s12*EijI(9,2)-2*(Dij2345I(6,2)+EijI(23,3
     -   )-EijI(24,3))+(EijI(3,1)+EijI(6,2)-EijI(11,3)-EijI(14,3)+2*E
     -   ijI(19,3))*P(9)+(EijI(4,1)+EijI(7,2))*P(47)-EijI(4,3)*P(90)-
     -   EijI(3,2)*P(114)-EijI(3,3)*P(128)+(EijI(12,3)+EijI(15,3)-2*E
     -   ijI(20,3))*P(152)+EijI(4,2)*P(403)+EijI(10,2)*P(513)+EijI(13
     -   ,3)*P(514)-EijI(16,3)*P(515))
       F(436)=DCMPLX(FR(436),FI(436))
       P(516) = p2sq+2*s12
       P(517) = p2sq-2*P(128)
       P(518) = p2sq-p3sq+s45
       P(519) = -s15-3*P(55)+2*P(518)
       P(520) = p2sq-s34+2*P(114)
       FR(437) = -4*(-Dij2345R(3,2)+Dij2345R(6,2)-p2sq*EijR(2,1)+s12*(
     -   EijR(3,1)+EijR(9,2))-s34*EijR(10,2)-6*EijR(11,2)+2*(D02345R+
     -   Dij2345R(2,1)+EijR(21,3)-EijR(23,3))+(EijR(4,1)+EijR(6,3)-Ei
     -   jR(7,3))*P(9)+EijR(4,3)*P(90)-(EijR(3,2)-EijR(11,3)+EijR(13,
     -   3))*P(128)+(EijR(8,2)-EijR(15,3)-EijR(17,3)+EijR(18,3)+EijR(
     -   20,3))*P(152)-EijR(4,2)*P(403)+EijR(14,3)*P(461)+EijR(16,3)*
     -   P(477)-EijR(5,2)*P(516)-EijR(6,2)*P(517)+EijR(7,2)*P(519)-Ei
     -   jR(19,3)*P(520))
       FI(437) = -4*(-Dij2345I(3,2)+Dij2345I(6,2)-p2sq*EijI(2,1)+s12*(
     -   EijI(3,1)+EijI(9,2))-s34*EijI(10,2)-6*EijI(11,2)+2*(D02345I+
     -   Dij2345I(2,1)+EijI(21,3)-EijI(23,3))+(EijI(4,1)+EijI(6,3)-Ei
     -   jI(7,3))*P(9)+EijI(4,3)*P(90)-(EijI(3,2)-EijI(11,3)+EijI(13,
     -   3))*P(128)+(EijI(8,2)-EijI(15,3)-EijI(17,3)+EijI(18,3)+EijI(
     -   20,3))*P(152)-EijI(4,2)*P(403)+EijI(14,3)*P(461)+EijI(16,3)*
     -   P(477)-EijI(5,2)*P(516)-EijI(6,2)*P(517)+EijI(7,2)*P(519)-Ei
     -   jI(19,3)*P(520))
       F(437)=DCMPLX(FR(437),FI(437))
       P(521) = p2sq+s12-s23
       P(522) = 3*p3sq-s45+2*P(521)
       P(523) = 2*p2sq+s12
       P(524) = -p2sq+p3sq+s15+s45+2*P(463)
       P(525) = p2sq-s15-2*s23-s45+3*P(3)
       P(526) = p4sq-s12-s34
       P(527) = 5*p3sq+s45+2*P(526)
       FR(438) = 4*(-Dij2345R(1,1)+Dij2345R(3,1)+Dij2345R(3,2)+Dij2345
     -   R(4,2)-Dij2345R(5,2)-Dij2345R(6,2)+s15*(EE0R+EijR(1,1))-4*p3
     -   sq*EijR(3,2)-18*EijR(11,2)+(EijR(2,1)+3*EijR(5,2)-EijR(14,3)
     -   -EijR(17,3)+EijR(18,3)+EijR(19,3))*P(9)+2*(D02345R-EijR(22,3
     -   )+EijR(23,3)-EijR(6,2)*P(14))+EijR(4,1)*P(47)-EijR(4,3)*P(90
     -   )+(-EijR(12,3)+EijR(13,3))*P(128)+(EijR(9,3)-EijR(10,3))*P(1
     -   52)-EijR(7,2)*P(181)+EijR(15,3)*P(462)-EijR(16,3)*P(477)-Eij
     -   R(20,3)*P(478)+EijR(8,2)*P(522)-EijR(2,2)*P(523)-EijR(4,2)*P
     -   (524)-EijR(9,2)*P(525)+EijR(10,2)*P(527))
       FI(438) = 4*(-Dij2345I(1,1)+Dij2345I(3,1)+Dij2345I(3,2)+Dij2345
     -   I(4,2)-Dij2345I(5,2)-Dij2345I(6,2)+s15*(EE0I+EijI(1,1))-4*p3
     -   sq*EijI(3,2)-18*EijI(11,2)+(EijI(2,1)+3*EijI(5,2)-EijI(14,3)
     -   -EijI(17,3)+EijI(18,3)+EijI(19,3))*P(9)+2*(D02345I-EijI(22,3
     -   )+EijI(23,3)-EijI(6,2)*P(14))+EijI(4,1)*P(47)-EijI(4,3)*P(90
     -   )+(-EijI(12,3)+EijI(13,3))*P(128)+(EijI(9,3)-EijI(10,3))*P(1
     -   52)-EijI(7,2)*P(181)+EijI(15,3)*P(462)-EijI(16,3)*P(477)-Eij
     -   I(20,3)*P(478)+EijI(8,2)*P(522)-EijI(2,2)*P(523)-EijI(4,2)*P
     -   (524)-EijI(9,2)*P(525)+EijI(10,2)*P(527))
       F(438)=DCMPLX(FR(438),FI(438))
       P(528) = -p3sq+s12+s34+s45
       FR(439) = 4*(D02345R+Dij2345R(2,1)+Dij2345R(3,2)-Dij2345R(6,2)+
     -   s12*EijR(2,1)+p2sq*EijR(9,2)-2*(EijR(11,2)-EijR(24,3))-(EijR
     -   (6,2)+2*EijR(7,2)+EijR(14,3)-EijR(19,3))*P(9)-EijR(4,3)*P(90
     -   )-EijR(3,1)*P(114)+(-EijR(3,2)+EijR(13,3))*P(128)+(EijR(8,2)
     -   +EijR(15,3)-EijR(20,3))*P(152)-EijR(4,1)*P(195)+EijR(4,2)*P(
     -   283)-EijR(16,3)*P(477)-EijR(10,2)*P(528))
       FI(439) = 4*(D02345I+Dij2345I(2,1)+Dij2345I(3,2)-Dij2345I(6,2)+
     -   s12*EijI(2,1)+p2sq*EijI(9,2)-2*(EijI(11,2)-EijI(24,3))-(EijI
     -   (6,2)+2*EijI(7,2)+EijI(14,3)-EijI(19,3))*P(9)-EijI(4,3)*P(90
     -   )-EijI(3,1)*P(114)+(-EijI(3,2)+EijI(13,3))*P(128)+(EijI(8,2)
     -   +EijI(15,3)-EijI(20,3))*P(152)-EijI(4,1)*P(195)+EijI(4,2)*P(
     -   283)-EijI(16,3)*P(477)-EijI(10,2)*P(528))
       F(439)=DCMPLX(FR(439),FI(439))
       P(529) = p2sq-3*s12-s15+s34
       FR(440) = 4*(p2sq*EijR(8,2)-p2sq*EijR(9,2)+2*(EijR(23,3)-EijR(2
     -   4,3)+EijR(7,2)*P(9)+EijR(6,2)*P(47))+(-EijR(3,1)+EijR(4,1))*
     -   P(141)+(EijR(3,2)-EijR(10,2))*P(529))
       FI(440) = 4*(p2sq*EijI(8,2)-p2sq*EijI(9,2)+2*(EijI(23,3)-EijI(2
     -   4,3)+EijI(7,2)*P(9)+EijI(6,2)*P(47))+(-EijI(3,1)+EijI(4,1))*
     -   P(141)+(EijI(3,2)-EijI(10,2))*P(529))
       F(440)=DCMPLX(FR(440),FI(440))
       P(530) = p2sq-s12+s15-s34-2*P(114)
       P(531) = p2sq-s15-3*P(55)-2*P(114)
       FR(441) = 4*(-(p2sq*(EijR(2,1)+EijR(9,2)))+EijR(3,1)*P(119)+Eij
     -   R(4,1)*P(141)+2*(D02345R+Dij2345R(2,1)-3*EijR(11,2)+EijR(21,
     -   3)-EijR(24,3)-EijR(3,2)*P(128)+EijR(8,2)*P(152)+EijR(7,2)*P(
     -   402))-EijR(5,2)*P(516)-EijR(6,2)*P(530)-EijR(10,2)*P(531))
       FI(441) = 4*(-(p2sq*(EijI(2,1)+EijI(9,2)))+EijI(3,1)*P(119)+Eij
     -   I(4,1)*P(141)+2*(D02345I+Dij2345I(2,1)-3*EijI(11,2)+EijI(21,
     -   3)-EijI(24,3)-EijI(3,2)*P(128)+EijI(8,2)*P(152)+EijI(7,2)*P(
     -   402))-EijI(5,2)*P(516)-EijI(6,2)*P(530)-EijI(10,2)*P(531))
       F(441)=DCMPLX(FR(441),FI(441))
       P(532) = 3*p3sq-s12+s45
       P(533) = s15-2*s23-s34-s45+3*P(106)
       P(534) = p2sq-3*p3sq+s34+s45+2*P(48)
       P(535) = p2sq+2*p4sq-5*s12-s15+3*P(6)
       FR(442) = -4*(-Dij2345R(1,1)-Dij2345R(2,1)+s15*(EE0R+EijR(1,1))
     -   -18*EijR(11,2)+EijR(2,1)*P(9)-3*EijR(5,2)*P(47)-EijR(4,1)*P(
     -   141)-EijR(6,2)*P(180)-EijR(3,1)*P(283)+2*(D02345R+Dij2345R(3
     -   ,1)-EijR(22,3)+EijR(24,3)-EijR(4,2)*P(213)-EijR(7,2)*P(399))
     -   -EijR(2,2)*P(523)-EijR(3,2)*P(532)+EijR(8,2)*P(533)+EijR(9,2
     -   )*P(534)+EijR(10,2)*P(535))
       FI(442) = -4*(-Dij2345I(1,1)-Dij2345I(2,1)+s15*(EE0I+EijI(1,1))
     -   -18*EijI(11,2)+EijI(2,1)*P(9)-3*EijI(5,2)*P(47)-EijI(4,1)*P(
     -   141)-EijI(6,2)*P(180)-EijI(3,1)*P(283)+2*(D02345I+Dij2345I(3
     -   ,1)-EijI(22,3)+EijI(24,3)-EijI(4,2)*P(213)-EijI(7,2)*P(399))
     -   -EijI(2,2)*P(523)-EijI(3,2)*P(532)+EijI(8,2)*P(533)+EijI(9,2
     -   )*P(534)+EijI(10,2)*P(535))
       F(442)=DCMPLX(FR(442),FI(442))
       FR(443) = -4*(D02345R+Dij2345R(2,1)+s12*EijR(2,1)+p2sq*EijR(9,2
     -   )-2*(EijR(11,2)-EijR(24,3))-(EijR(6,2)+2*EijR(7,2))*P(9)-Eij
     -   R(3,1)*P(114)-EijR(3,2)*P(128)+EijR(8,2)*P(152)+EijR(10,2)*P
     -   (174)-EijR(4,1)*P(195))
       FI(443) = -4*(D02345I+Dij2345I(2,1)+s12*EijI(2,1)+p2sq*EijI(9,2
     -   )-2*(EijI(11,2)-EijI(24,3))-(EijI(6,2)+2*EijI(7,2))*P(9)-Eij
     -   I(3,1)*P(114)-EijI(3,2)*P(128)+EijI(8,2)*P(152)+EijI(10,2)*P
     -   (174)-EijI(4,1)*P(195))
       F(443)=DCMPLX(FR(443),FI(443))
       P(536) = 3*p2sq-s12-s15+s34
       P(537) = -s15+3*P(9)
       P(538) = -s15+s23-s45+2*P(9)
       P(539) = -p2sq+s23+s45-2*P(106)
       P(540) = p2sq-p3sq+s34
       P(541) = -s12-s15+s23+s45+2*P(540)
       P(542) = p4sq-s34+s45+2*P(3)
       FR(444) = 4*(D02345R+Dij2345R(3,2)+Dij2345R(4,2)-Dij2345R(5,2)-
     -   Dij2345R(6,2)-12*EijR(11,2)+6*EijR(22,3)-2*(Dij2345R(1,1)-Di
     -   j2345R(3,1)+p3sq*EijR(3,2)-EijR(23,3))-8*EijR(24,3)+(-EijR(1
     -   4,3)-EijR(17,3)+EijR(18,3)+EijR(19,3))*P(13)-EijR(6,2)*P(14)
     -   +(-EijR(9,3)+EijR(10,3))*P(26)+(EijR(12,3)-EijR(13,3))*P(38)
     -   +EijR(2,1)*P(55)-(EijR(4,2)+EijR(4,3))*P(213)+EijR(4,1)*P(26
     -   7)-EijR(15,3)*P(503)+EijR(16,3)*P(511)-EijR(20,3)*P(512)-Eij
     -   R(2,2)*P(536)+EijR(5,2)*P(537)-EijR(7,2)*P(538)-EijR(8,2)*P(
     -   539)+EijR(9,2)*P(541)+EijR(10,2)*P(542))
       FI(444) = 4*(D02345I+Dij2345I(3,2)+Dij2345I(4,2)-Dij2345I(5,2)-
     -   Dij2345I(6,2)-12*EijI(11,2)+6*EijI(22,3)-2*(Dij2345I(1,1)-Di
     -   j2345I(3,1)+p3sq*EijI(3,2)-EijI(23,3))-8*EijI(24,3)+(-EijI(1
     -   4,3)-EijI(17,3)+EijI(18,3)+EijI(19,3))*P(13)-EijI(6,2)*P(14)
     -   +(-EijI(9,3)+EijI(10,3))*P(26)+(EijI(12,3)-EijI(13,3))*P(38)
     -   +EijI(2,1)*P(55)-(EijI(4,2)+EijI(4,3))*P(213)+EijI(4,1)*P(26
     -   7)-EijI(15,3)*P(503)+EijI(16,3)*P(511)-EijI(20,3)*P(512)-Eij
     -   I(2,2)*P(536)+EijI(5,2)*P(537)-EijI(7,2)*P(538)-EijI(8,2)*P(
     -   539)+EijI(9,2)*P(541)+EijI(10,2)*P(542))
       F(444)=DCMPLX(FR(444),FI(444))
       P(543) = 2*p2sq-P(218)
       P(544) = s15-2*P(9)
       P(545) = -s15+s34+2*P(9)
       P(546) = p4sq+3*s45+2*P(3)
       P(547) = p3sq+2*p4sq-s12+3*s45
       FR(445) = 4*(-Dij2345R(2,1)+Dij2345R(2,2)+Dij2345R(3,1)+Dij2345
     -   R(3,2)-2*Dij2345R(6,2)+s12*EijR(3,2)-s12*EijR(10,2)+8*(EijR(
     -   23,3)-EijR(24,3))-(EijR(11,3)+EijR(14,3)-2*EijR(19,3))*P(13)
     -   -(EijR(12,3)+EijR(15,3)-2*EijR(20,3))*P(26)+EijR(3,3)*P(38)+
     -   EijR(3,1)*P(55)-EijR(4,3)*P(213)+EijR(4,1)*P(267)+EijR(6,2)*
     -   P(543)+EijR(7,2)*P(544)+(-EijR(8,2)+EijR(9,2))*P(545)-EijR(1
     -   3,3)*P(546)+EijR(16,3)*P(547))
       FI(445) = 4*(-Dij2345I(2,1)+Dij2345I(2,2)+Dij2345I(3,1)+Dij2345
     -   I(3,2)-2*Dij2345I(6,2)+s12*EijI(3,2)-s12*EijI(10,2)+8*(EijI(
     -   23,3)-EijI(24,3))-(EijI(11,3)+EijI(14,3)-2*EijI(19,3))*P(13)
     -   -(EijI(12,3)+EijI(15,3)-2*EijI(20,3))*P(26)+EijI(3,3)*P(38)+
     -   EijI(3,1)*P(55)-EijI(4,3)*P(213)+EijI(4,1)*P(267)+EijI(6,2)*
     -   P(543)+EijI(7,2)*P(544)+(-EijI(8,2)+EijI(9,2))*P(545)-EijI(1
     -   3,3)*P(546)+EijI(16,3)*P(547))
       F(445)=DCMPLX(FR(445),FI(445))
       FR(446) = 4*(D02345R+Dij2345R(2,1)+Dij2345R(3,1)+Dij2345R(3,2)-
     -   Dij2345R(6,2)+s12*EijR(2,1)-2*EijR(11,2)-6*EijR(24,3)-EijR(6
     -   ,2)*P(9)+(-EijR(14,3)+EijR(19,3))*P(13)+(-EijR(15,3)+EijR(20
     -   ,3))*P(26)-EijR(13,3)*P(38)-EijR(3,1)*P(114)-EijR(3,2)*P(128
     -   )+EijR(8,2)*P(152)-EijR(4,3)*P(213)-EijR(7,2)*P(269)-EijR(4,
     -   1)*P(281)+EijR(16,3)*P(511)-EijR(10,2)*P(528)+EijR(9,2)*P(54
     -   5))
       FI(446) = 4*(D02345I+Dij2345I(2,1)+Dij2345I(3,1)+Dij2345I(3,2)-
     -   Dij2345I(6,2)+s12*EijI(2,1)-2*EijI(11,2)-6*EijI(24,3)-EijI(6
     -   ,2)*P(9)+(-EijI(14,3)+EijI(19,3))*P(13)+(-EijI(15,3)+EijI(20
     -   ,3))*P(26)-EijI(13,3)*P(38)-EijI(3,1)*P(114)-EijI(3,2)*P(128
     -   )+EijI(8,2)*P(152)-EijI(4,3)*P(213)-EijI(7,2)*P(269)-EijI(4,
     -   1)*P(281)+EijI(16,3)*P(511)-EijI(10,2)*P(528)+EijI(9,2)*P(54
     -   5))
       F(446)=DCMPLX(FR(446),FI(446))
       P(548) = 4*s12+s15-s34
       P(549) = -s23-3*P(7)+2*P(27)
       P(550) = p2sq-p3sq-s12+s34
       P(551) = s15+s23-3*s45-2*P(550)
       P(552) = 2*p2sq+p3sq-3*s12-s15-s23+s34
       P(553) = 3*p3sq+p4sq-2*s34
       P(554) = p3sq+p4sq-s12-s23+3*s45
       FR(447) = -4*(-Dij2345R(3,1)-Dij2345R(3,2)+Dij2345R(6,2)+s15*(E
     -   E0R+2*EijR(1,1)+EijR(1,2))+s12*EijR(3,1)+2*(D02345R+Dij2345R
     -   (2,1)-6*EijR(11,2)-3*EijR(21,3)-EijR(23,3)+4*EijR(24,3))+(Ei
     -   jR(6,3)-EijR(7,3))*P(13)+(EijR(15,3)+EijR(17,3)-EijR(18,3)-E
     -   ijR(20,3))*P(26)+(-EijR(11,3)+EijR(13,3))*P(38)+EijR(4,1)*P(
     -   55)-EijR(3,2)*P(208)+(-EijR(4,2)+EijR(4,3))*P(213)-EijR(2,1)
     -   *P(289)+EijR(8,2)*P(496)-EijR(14,3)*P(502)-EijR(16,3)*P(511)
     -   -EijR(5,2)*P(548)-EijR(6,2)*P(549)-EijR(7,2)*P(551)-EijR(9,2
     -   )*P(552)+EijR(10,2)*P(553)+EijR(19,3)*P(554))
       FI(447) = -4*(-Dij2345I(3,1)-Dij2345I(3,2)+Dij2345I(6,2)+s15*(E
     -   E0I+2*EijI(1,1)+EijI(1,2))+s12*EijI(3,1)+2*(D02345I+Dij2345I
     -   (2,1)-6*EijI(11,2)-3*EijI(21,3)-EijI(23,3)+4*EijI(24,3))+(Ei
     -   jI(6,3)-EijI(7,3))*P(13)+(EijI(15,3)+EijI(17,3)-EijI(18,3)-E
     -   ijI(20,3))*P(26)+(-EijI(11,3)+EijI(13,3))*P(38)+EijI(4,1)*P(
     -   55)-EijI(3,2)*P(208)+(-EijI(4,2)+EijI(4,3))*P(213)-EijI(2,1)
     -   *P(289)+EijI(8,2)*P(496)-EijI(14,3)*P(502)-EijI(16,3)*P(511)
     -   -EijI(5,2)*P(548)-EijI(6,2)*P(549)-EijI(7,2)*P(551)-EijI(9,2
     -   )*P(552)+EijI(10,2)*P(553)+EijI(19,3)*P(554))
       F(447)=DCMPLX(FR(447),FI(447))
       FR(448) = 2*(s15*(EE0R+EijR(1,1))+3*(D02345R-Dij2345R(1,1)+Dij2
     -   345R(3,1)-2*p3sq*EijR(3,2)-10*EijR(11,2)+EijR(5,2)*P(9)-EijR
     -   (7,2)*P(13)-EijR(6,2)*P(14)-EijR(2,2)*P(152)-EijR(4,2)*P(213
     -   )+EijR(10,2)*P(379)+EijR(9,2)*P(384)+EijR(8,2)*P(385)))
       FI(448) = 2*(s15*(EE0I+EijI(1,1))+3*(D02345I-Dij2345I(1,1)+Dij2
     -   345I(3,1)-2*p3sq*EijI(3,2)-10*EijI(11,2)+EijI(5,2)*P(9)-EijI
     -   (7,2)*P(13)-EijI(6,2)*P(14)-EijI(2,2)*P(152)-EijI(4,2)*P(213
     -   )+EijI(10,2)*P(379)+EijI(9,2)*P(384)+EijI(8,2)*P(385)))
       F(448)=DCMPLX(FR(448),FI(448))
       KI(1) = -Dij1345I(3,3)-Dij1345I(8,3)+2*Dij1345I(9,3)
       KR(1) = -Dij1345R(3,3)-Dij1345R(8,3)+2*Dij1345R(9,3)
       K(1)=DCMPLX(KR(1),KI(1))
       KI(2) = EijI(4,4)-2*EijI(15,4)+EijI(21,4)
       KR(2) = EijR(4,4)-2*EijR(15,4)+EijR(21,4)
       K(2)=DCMPLX(KR(2),KI(2))
       KI(3) = EijI(4,4)-EijI(13,4)-3*(EijI(16,4)-EijI(22,4))
       KR(3) = EijR(4,4)-EijR(13,4)-3*(EijR(16,4)-EijR(22,4))
       K(3)=DCMPLX(KR(3),KI(3))
       KI(4) = EijI(4,4)-2*EijI(16,4)+EijI(22,4)
       KR(4) = EijR(4,4)-2*EijR(16,4)+EijR(22,4)
       K(4)=DCMPLX(KR(4),KI(4))
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
       call FaFunctionABEmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2
     -   mup4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,
     -   p5mup3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
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
c$$$       Print*," SMB(1) ", SMB(1)
c$$$       Print*," SMB(2) ", SMB(2)
c$$$       Print*," SMB(3) ", SMB(3)
c$$$       Print*," SMB(4) ", SMB(4)
c$$$       Print*," SMB(5) ", SMB(5)
c$$$       Print*," SMB(6) ", SMB(6)
c$$$       Print*," SMB(7) ", SMB(7)
c$$$       Print*," SMB(8) ", SMB(8)
c$$$       Print*," SMB(9) ", SMB(9)
c$$$       Print*," SMB(10) ", SMB(10)
c$$$       Print*," SMB(11) ", SMB(11)
c$$$       Print*," SMB(12) ", SMB(12)
c$$$       Print*," SMB(13) ", SMB(13)
c$$$       Print*," SMB(14) ", SMB(14)
c$$$       Print*," SMB(15) ", SMB(15)
c$$$       Print*," SMB(16) ", SMB(16)
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
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       result(2) = Fa(1)*SMB(1)+Fa(2)*SMB(2)+Fa(3)*SMB(3)+Fa(4)*SMB(4)
     -   +Fa(5)*SMB(5)+F(3)*SMB(6)-2*p3mup4*F(2)*SMB(7)+Fa(6)*SMB(8)-
     -   2*mup3mup4*F(2)*SMB(10)+2*mup2mup4*F(2)*SMB(11)-2*mup2mup3*F
     -   (2)*SMB(12)+2*p2mup4*F(2)*SMB(13)-2*p2mup3*F(2)*SMB(14)-2*p1
     -   mup2*F(2)*SMB(15)+F(2)*SMB(16)
       result(4) = Fa(7)*SMB(1)+Fa(8)*SMB(2)+Fa(9)*SMB(3)+Fa(10)*SMB(4
     -   )+Fa(11)*SMB(5)+F(16)*SMB(6)+p3mup4*F(17)*SMB(7)+Fa(12)*SMB(
     -   8)+mup3mup4*F(17)*SMB(10)+mup2mup4*F(19)*SMB(11)+mup2mup3*F(
     -   17)*SMB(12)+p2mup4*F(19)*SMB(13)+Fa(13)*SMB(14)+p1mup2*F(17)
     -   *SMB(15)+F(22)*SMB(16)
       result(3) = Fa(14)*SMB(1)+Fa(15)*SMB(2)+Fa(16)*SMB(3)+Fa(17)*SM
     -   B(4)+Fa(18)*SMB(5)+F(31)*SMB(6)+p3mup4*F(32)*SMB(7)+Fa(19)*S
     -   MB(8)+Fa(20)*SMB(9)+mup3mup4*F(32)*SMB(10)+mup2mup4*F(33)*SM
     -   B(11)+mup2mup3*F(32)*SMB(12)+p2mup4*F(33)*SMB(13)+p2mup3*F(3
     -   2)*SMB(14)+Fa(21)*SMB(15)+F(35)*SMB(16)
       result(5) = Fa(22)*SMB(1)+Fa(23)*SMB(2)+Fa(24)*SMB(3)+Fa(25)*SM
     -   B(4)+Fa(26)*SMB(5)+F(44)*SMB(6)+Fa(27)*SMB(7)+Fa(28)*SMB(8)+
     -   mup3mup4*F(46)*SMB(10)+mup2mup4*F(47)*SMB(11)+mup2mup3*F(46)
     -   *SMB(12)+p2mup4*F(47)*SMB(13)+p2mup3*F(46)*SMB(14)+p1mup2*F(
     -   46)*SMB(15)+F(48)*SMB(16)
       result(6) = Fa(29)*SMB(1)+Fa(30)*SMB(2)+Fa(31)*SMB(3)+Fa(32)*SM
     -   B(4)+Fa(33)*SMB(5)+F(102)*SMB(6)+p3mup4*F(103)*SMB(7)+Fa(34)
     -   *SMB(8)+Fa(35)*SMB(9)+mup3mup4*F(103)*SMB(10)+mup2mup4*F(110
     -   )*SMB(11)+Fa(36)*SMB(12)+p2mup4*F(110)*SMB(13)+Fa(37)*SMB(14
     -   )+Fa(38)*SMB(15)+F(124)*SMB(16)
       result(7) = Fa(39)*SMB(1)+Fa(40)*SMB(2)+Fa(41)*SMB(3)+Fa(42)*SM
     -   B(4)+Fa(43)*SMB(5)+F(173)*SMB(6)+Fa(44)*SMB(7)+Fa(45)*SMB(8)
     -   +Fa(46)*SMB(10)+mup2mup4*F(182)*SMB(11)+mup2mup3*F(183)*SMB(
     -   12)+Fa(47)*SMB(13)+Fa(48)*SMB(14)+p1mup2*F(183)*SMB(15)+F(18
     -   6)*SMB(16)
       result(8) = Fa(49)*SMB(1)+Fa(50)*SMB(2)+Fa(51)*SMB(3)+Fa(52)*SM
     -   B(4)+Fa(53)*SMB(5)+F(375)*SMB(6)+Fa(54)*SMB(7)+Fa(55)*SMB(8)
     -   +Fa(56)*SMB(9)+Fa(57)*SMB(10)+Fa(58)*SMB(11)+Fa(59)*SMB(12)+
     -   Fa(60)*SMB(13)+Fa(61)*SMB(14)+Fa(62)*SMB(15)+F(448)*SMB(16)

c************************************************************************************
c** Selecting different colour factor****************
c************************************************************************************
       If (ngluon.eq.0) then
       result(1)=-(result(2)+result(3)+result(4)+result(5)+result(6)+r
     -   esult(7)+result(8))
       elseif (ngluon.eq.1) then
       If (posgluon.eq.2) then
       result(1)=-(result(2)+result(4)+result(5)+result(7))
       result(2)=-(result(3)+result(6)+result(8))
       elseIf (posgluon.eq.3) then
       result(1)=-(result(2)+result(3)+result(5))
       result(2)=-(result(4)+result(6)+result(7)+result(8))
       elseIf (posgluon.eq.4) then
       result(1)=-(result(2)+result(3)+result(4)+result(6))
       result(2)=-(result(5)+result(7)+result(8))
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     -   . Look to the heading for explanation" 
       endif
       elseif (ngluon.eq.2) then
c posgluon specifies the postion without a gluon!
c result(1): coefficient of t^b t^a color structure (i.e. the born color structure)
c result(2): coefficient of delta^ab color structure
       if(posgluon.eq.2) then
       result(1)=+(result(4)+result(5)+result(6)+result(7)+result(8)-8d0*(result(2)+result(3)))/6d0
       result(2)=-(result(7)+result(8))/4d0
       elseif(posgluon.eq.3) then
       result(1)=+(result(3)+result(5)+result(6)+result(7)+result(8)-8d0*(result(2)+result(4)))/6d0
       result(2)=-result(8)/4d0
       elseif(posgluon.eq.4) then
       result(1)=+(result(3)+result(4)+result(6)+result(7)+result(8)-8d0*(result(2)+result(5)))/6d0
       result(2)=-(result(6)+result(8))/4d0
       else
       write(*,*) "Error: Wrong position of gluon in PenlineAbe"
       endif
       endif
       resultb = Is45*(2*(mup3mup4*SMB(1)-mup2mup4*SMB(2)+mup2mup3*SMB
     -   (3))-SMB(6))+Is12s45*(4*((mup2mup4*(p1mup3+p2mup3)+mup2mup3*
     -   p3mup4)*SMB(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2mup3*p2m
     -   up4)*SMB(5))+(p2sq+p3sq-s23)*SMB(6)-2*((2*(p1mup3*p2mup4+p2m
     -   up3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMB(1)-(-2*p1
     -   mup2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMB(2)+(2*p1mup2*(p1mu
     -   p3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMB(3)+p3mup4*SMB(7)+(p
     -   1mup3+p2mup3)*SMB(8)+mup3mup4*SMB(10)-mup2mup4*SMB(11)+mup2m
     -   up3*SMB(12)-p2mup4*SMB(13)+p2mup3*SMB(14)+p1mup2*SMB(15))+SM
     -   B(16))
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
       call FaFunctionABEmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2
     -   mup4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,
     -   p5mup3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
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
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*SMG(4
     -   )+Fa(5)*SMG(5)+F(3)*SMG(6)-2*p3mup4*F(2)*SMG(7)+Fa(6)*SMG(8)
     -   -2*mup3mup4*F(2)*SMG(10)+2*mup2mup4*F(2)*SMG(11)-2*mup2mup3*
     -   F(2)*SMG(12)+2*p2mup4*F(2)*SMG(13)-2*p2mup3*F(2)*SMG(14)-2*p
     -   1mup2*F(2)*SMG(15)+F(2)*SMG(16)
       resultg(4) = Fa(7)*SMG(1)+Fa(8)*SMG(2)+Fa(9)*SMG(3)+Fa(10)*SMG(
     -   4)+Fa(11)*SMG(5)+F(16)*SMG(6)+p3mup4*F(17)*SMG(7)+Fa(12)*SMG
     -   (8)+mup3mup4*F(17)*SMG(10)+mup2mup4*F(19)*SMG(11)+mup2mup3*F
     -   (17)*SMG(12)+p2mup4*F(19)*SMG(13)+Fa(13)*SMG(14)+p1mup2*F(17
     -   )*SMG(15)+F(22)*SMG(16)
       resultg(3) = Fa(14)*SMG(1)+Fa(15)*SMG(2)+Fa(16)*SMG(3)+Fa(17)*S
     -   MG(4)+Fa(18)*SMG(5)+F(31)*SMG(6)+p3mup4*F(32)*SMG(7)+Fa(19)*
     -   SMG(8)+Fa(20)*SMG(9)+mup3mup4*F(32)*SMG(10)+mup2mup4*F(33)*S
     -   MG(11)+mup2mup3*F(32)*SMG(12)+p2mup4*F(33)*SMG(13)+p2mup3*F(
     -   32)*SMG(14)+Fa(21)*SMG(15)+F(35)*SMG(16)
       resultg(5) = Fa(22)*SMG(1)+Fa(23)*SMG(2)+Fa(24)*SMG(3)+Fa(25)*S
     -   MG(4)+Fa(26)*SMG(5)+F(44)*SMG(6)+Fa(27)*SMG(7)+Fa(28)*SMG(8)
     -   +mup3mup4*F(46)*SMG(10)+mup2mup4*F(47)*SMG(11)+mup2mup3*F(46
     -   )*SMG(12)+p2mup4*F(47)*SMG(13)+p2mup3*F(46)*SMG(14)+p1mup2*F
     -   (46)*SMG(15)+F(48)*SMG(16)
       resultg(6) = Fa(29)*SMG(1)+Fa(30)*SMG(2)+Fa(31)*SMG(3)+Fa(32)*S
     -   MG(4)+Fa(33)*SMG(5)+F(102)*SMG(6)+p3mup4*F(103)*SMG(7)+Fa(34
     -   )*SMG(8)+Fa(35)*SMG(9)+mup3mup4*F(103)*SMG(10)+mup2mup4*F(11
     -   0)*SMG(11)+Fa(36)*SMG(12)+p2mup4*F(110)*SMG(13)+Fa(37)*SMG(1
     -   4)+Fa(38)*SMG(15)+F(124)*SMG(16)
       resultg(7) = Fa(39)*SMG(1)+Fa(40)*SMG(2)+Fa(41)*SMG(3)+Fa(42)*S
     -   MG(4)+Fa(43)*SMG(5)+F(173)*SMG(6)+Fa(44)*SMG(7)+Fa(45)*SMG(8
     -   )+Fa(46)*SMG(10)+mup2mup4*F(182)*SMG(11)+mup2mup3*F(183)*SMG
     -   (12)+Fa(47)*SMG(13)+Fa(48)*SMG(14)+p1mup2*F(183)*SMG(15)+F(1
     -   86)*SMG(16)
       resultg(8) = Fa(49)*SMG(1)+Fa(50)*SMG(2)+Fa(51)*SMG(3)+Fa(52)*S
     -   MG(4)+Fa(53)*SMG(5)+F(375)*SMG(6)+Fa(54)*SMG(7)+Fa(55)*SMG(8
     -   )+Fa(56)*SMG(9)+Fa(57)*SMG(10)+Fa(58)*SMG(11)+Fa(59)*SMG(12)
     -   +Fa(60)*SMG(13)+Fa(61)*SMG(14)+Fa(62)*SMG(15)+F(448)*SMG(16)
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(1)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5)+re
     -   sultg(6)+resultg(7)+resultg(8))/Sqrt(abs(p2sq)*abs(s15)*abs(
     -   p4sq))
       resultgaugeb(1) = Is45*(2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2
     -   mup3*SMG(3))-SMG(6))+Is12s45*(4*((mup2mup4*(p1mup3+p2mup3)+m
     -   up2mup3*p3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2
     -   mup3*p2mup4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2
     -   mup4+p2mup3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1
     -   )-(-2*p1mup2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mu
     -   p2*(p1mup3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*S
     -   MG(7)+(p1mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(1
     -   1)+mup2mup3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG
     -   (15))+SMG(16))
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
       call FaFunctionABEmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2
     -   mup4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,
     -   p5mup3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
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
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*SMG(4
     -   )+Fa(5)*SMG(5)+F(3)*SMG(6)-2*p3mup4*F(2)*SMG(7)+Fa(6)*SMG(8)
     -   -2*mup3mup4*F(2)*SMG(10)+2*mup2mup4*F(2)*SMG(11)-2*mup2mup3*
     -   F(2)*SMG(12)+2*p2mup4*F(2)*SMG(13)-2*p2mup3*F(2)*SMG(14)-2*p
     -   1mup2*F(2)*SMG(15)+F(2)*SMG(16)
       resultg(4) = Fa(7)*SMG(1)+Fa(8)*SMG(2)+Fa(9)*SMG(3)+Fa(10)*SMG(
     -   4)+Fa(11)*SMG(5)+F(16)*SMG(6)+p3mup4*F(17)*SMG(7)+Fa(12)*SMG
     -   (8)+mup3mup4*F(17)*SMG(10)+mup2mup4*F(19)*SMG(11)+mup2mup3*F
     -   (17)*SMG(12)+p2mup4*F(19)*SMG(13)+Fa(13)*SMG(14)+p1mup2*F(17
     -   )*SMG(15)+F(22)*SMG(16)
       resultg(3) = Fa(14)*SMG(1)+Fa(15)*SMG(2)+Fa(16)*SMG(3)+Fa(17)*S
     -   MG(4)+Fa(18)*SMG(5)+F(31)*SMG(6)+p3mup4*F(32)*SMG(7)+Fa(19)*
     -   SMG(8)+Fa(20)*SMG(9)+mup3mup4*F(32)*SMG(10)+mup2mup4*F(33)*S
     -   MG(11)+mup2mup3*F(32)*SMG(12)+p2mup4*F(33)*SMG(13)+p2mup3*F(
     -   32)*SMG(14)+Fa(21)*SMG(15)+F(35)*SMG(16)
       resultg(5) = Fa(22)*SMG(1)+Fa(23)*SMG(2)+Fa(24)*SMG(3)+Fa(25)*S
     -   MG(4)+Fa(26)*SMG(5)+F(44)*SMG(6)+Fa(27)*SMG(7)+Fa(28)*SMG(8)
     -   +mup3mup4*F(46)*SMG(10)+mup2mup4*F(47)*SMG(11)+mup2mup3*F(46
     -   )*SMG(12)+p2mup4*F(47)*SMG(13)+p2mup3*F(46)*SMG(14)+p1mup2*F
     -   (46)*SMG(15)+F(48)*SMG(16)
       resultg(6) = Fa(29)*SMG(1)+Fa(30)*SMG(2)+Fa(31)*SMG(3)+Fa(32)*S
     -   MG(4)+Fa(33)*SMG(5)+F(102)*SMG(6)+p3mup4*F(103)*SMG(7)+Fa(34
     -   )*SMG(8)+Fa(35)*SMG(9)+mup3mup4*F(103)*SMG(10)+mup2mup4*F(11
     -   0)*SMG(11)+Fa(36)*SMG(12)+p2mup4*F(110)*SMG(13)+Fa(37)*SMG(1
     -   4)+Fa(38)*SMG(15)+F(124)*SMG(16)
       resultg(7) = Fa(39)*SMG(1)+Fa(40)*SMG(2)+Fa(41)*SMG(3)+Fa(42)*S
     -   MG(4)+Fa(43)*SMG(5)+F(173)*SMG(6)+Fa(44)*SMG(7)+Fa(45)*SMG(8
     -   )+Fa(46)*SMG(10)+mup2mup4*F(182)*SMG(11)+mup2mup3*F(183)*SMG
     -   (12)+Fa(47)*SMG(13)+Fa(48)*SMG(14)+p1mup2*F(183)*SMG(15)+F(1
     -   86)*SMG(16)
       resultg(8) = Fa(49)*SMG(1)+Fa(50)*SMG(2)+Fa(51)*SMG(3)+Fa(52)*S
     -   MG(4)+Fa(53)*SMG(5)+F(375)*SMG(6)+Fa(54)*SMG(7)+Fa(55)*SMG(8
     -   )+Fa(56)*SMG(9)+Fa(57)*SMG(10)+Fa(58)*SMG(11)+Fa(59)*SMG(12)
     -   +Fa(60)*SMG(13)+Fa(61)*SMG(14)+Fa(62)*SMG(15)+F(448)*SMG(16)
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(1)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5)+re
     -   sultg(6)+resultg(7)+resultg(8))
       resultgaugeb(1) = Is45*(2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2
     -   mup3*SMG(3))-SMG(6))+Is12s45*(4*((mup2mup4*(p1mup3+p2mup3)+m
     -   up2mup3*p3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2
     -   mup3*p2mup4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2
     -   mup4+p2mup3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1
     -   )-(-2*p1mup2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mu
     -   p2*(p1mup3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*S
     -   MG(7)+(p1mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(1
     -   1)+mup2mup3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG
     -   (15))+SMG(16))
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
       call FaFunctionABEmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2
     -   mup4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,
     -   p5mup3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
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
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*SMG(4
     -   )+Fa(5)*SMG(5)+F(3)*SMG(6)-2*p3mup4*F(2)*SMG(7)+Fa(6)*SMG(8)
     -   -2*mup3mup4*F(2)*SMG(10)+2*mup2mup4*F(2)*SMG(11)-2*mup2mup3*
     -   F(2)*SMG(12)+2*p2mup4*F(2)*SMG(13)-2*p2mup3*F(2)*SMG(14)-2*p
     -   1mup2*F(2)*SMG(15)+F(2)*SMG(16)
       resultg(4) = Fa(7)*SMG(1)+Fa(8)*SMG(2)+Fa(9)*SMG(3)+Fa(10)*SMG(
     -   4)+Fa(11)*SMG(5)+F(16)*SMG(6)+p3mup4*F(17)*SMG(7)+Fa(12)*SMG
     -   (8)+mup3mup4*F(17)*SMG(10)+mup2mup4*F(19)*SMG(11)+mup2mup3*F
     -   (17)*SMG(12)+p2mup4*F(19)*SMG(13)+Fa(13)*SMG(14)+p1mup2*F(17
     -   )*SMG(15)+F(22)*SMG(16)
       resultg(3) = Fa(14)*SMG(1)+Fa(15)*SMG(2)+Fa(16)*SMG(3)+Fa(17)*S
     -   MG(4)+Fa(18)*SMG(5)+F(31)*SMG(6)+p3mup4*F(32)*SMG(7)+Fa(19)*
     -   SMG(8)+Fa(20)*SMG(9)+mup3mup4*F(32)*SMG(10)+mup2mup4*F(33)*S
     -   MG(11)+mup2mup3*F(32)*SMG(12)+p2mup4*F(33)*SMG(13)+p2mup3*F(
     -   32)*SMG(14)+Fa(21)*SMG(15)+F(35)*SMG(16)
       resultg(5) = Fa(22)*SMG(1)+Fa(23)*SMG(2)+Fa(24)*SMG(3)+Fa(25)*S
     -   MG(4)+Fa(26)*SMG(5)+F(44)*SMG(6)+Fa(27)*SMG(7)+Fa(28)*SMG(8)
     -   +mup3mup4*F(46)*SMG(10)+mup2mup4*F(47)*SMG(11)+mup2mup3*F(46
     -   )*SMG(12)+p2mup4*F(47)*SMG(13)+p2mup3*F(46)*SMG(14)+p1mup2*F
     -   (46)*SMG(15)+F(48)*SMG(16)
       resultg(6) = Fa(29)*SMG(1)+Fa(30)*SMG(2)+Fa(31)*SMG(3)+Fa(32)*S
     -   MG(4)+Fa(33)*SMG(5)+F(102)*SMG(6)+p3mup4*F(103)*SMG(7)+Fa(34
     -   )*SMG(8)+Fa(35)*SMG(9)+mup3mup4*F(103)*SMG(10)+mup2mup4*F(11
     -   0)*SMG(11)+Fa(36)*SMG(12)+p2mup4*F(110)*SMG(13)+Fa(37)*SMG(1
     -   4)+Fa(38)*SMG(15)+F(124)*SMG(16)
       resultg(7) = Fa(39)*SMG(1)+Fa(40)*SMG(2)+Fa(41)*SMG(3)+Fa(42)*S
     -   MG(4)+Fa(43)*SMG(5)+F(173)*SMG(6)+Fa(44)*SMG(7)+Fa(45)*SMG(8
     -   )+Fa(46)*SMG(10)+mup2mup4*F(182)*SMG(11)+mup2mup3*F(183)*SMG
     -   (12)+Fa(47)*SMG(13)+Fa(48)*SMG(14)+p1mup2*F(183)*SMG(15)+F(1
     -   86)*SMG(16)
       resultg(8) = Fa(49)*SMG(1)+Fa(50)*SMG(2)+Fa(51)*SMG(3)+Fa(52)*S
     -   MG(4)+Fa(53)*SMG(5)+F(375)*SMG(6)+Fa(54)*SMG(7)+Fa(55)*SMG(8
     -   )+Fa(56)*SMG(9)+Fa(57)*SMG(10)+Fa(58)*SMG(11)+Fa(59)*SMG(12)
     -   +Fa(60)*SMG(13)+Fa(61)*SMG(14)+Fa(62)*SMG(15)+F(448)*SMG(16)
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(1)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5)+re
     -   sultg(6)+resultg(7)+resultg(8))
       resultgaugeb(1) = Is45*(2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2
     -   mup3*SMG(3))-SMG(6))+Is12s45*(4*((mup2mup4*(p1mup3+p2mup3)+m
     -   up2mup3*p3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2
     -   mup3*p2mup4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2
     -   mup4+p2mup3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1
     -   )-(-2*p1mup2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mu
     -   p2*(p1mup3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*S
     -   MG(7)+(p1mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(1
     -   1)+mup2mup3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG
     -   (15))+SMG(16))
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
       call FaFunctionABEmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2
     -   mup4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,
     -   p5mup3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
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
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*SMG(4
     -   )+Fa(5)*SMG(5)+F(3)*SMG(6)-2*p3mup4*F(2)*SMG(7)+Fa(6)*SMG(8)
     -   -2*mup3mup4*F(2)*SMG(10)+2*mup2mup4*F(2)*SMG(11)-2*mup2mup3*
     -   F(2)*SMG(12)+2*p2mup4*F(2)*SMG(13)-2*p2mup3*F(2)*SMG(14)-2*p
     -   1mup2*F(2)*SMG(15)+F(2)*SMG(16)
       resultg(4) = Fa(7)*SMG(1)+Fa(8)*SMG(2)+Fa(9)*SMG(3)+Fa(10)*SMG(
     -   4)+Fa(11)*SMG(5)+F(16)*SMG(6)+p3mup4*F(17)*SMG(7)+Fa(12)*SMG
     -   (8)+mup3mup4*F(17)*SMG(10)+mup2mup4*F(19)*SMG(11)+mup2mup3*F
     -   (17)*SMG(12)+p2mup4*F(19)*SMG(13)+Fa(13)*SMG(14)+p1mup2*F(17
     -   )*SMG(15)+F(22)*SMG(16)
       resultg(3) = Fa(14)*SMG(1)+Fa(15)*SMG(2)+Fa(16)*SMG(3)+Fa(17)*S
     -   MG(4)+Fa(18)*SMG(5)+F(31)*SMG(6)+p3mup4*F(32)*SMG(7)+Fa(19)*
     -   SMG(8)+Fa(20)*SMG(9)+mup3mup4*F(32)*SMG(10)+mup2mup4*F(33)*S
     -   MG(11)+mup2mup3*F(32)*SMG(12)+p2mup4*F(33)*SMG(13)+p2mup3*F(
     -   32)*SMG(14)+Fa(21)*SMG(15)+F(35)*SMG(16)
       resultg(5) = Fa(22)*SMG(1)+Fa(23)*SMG(2)+Fa(24)*SMG(3)+Fa(25)*S
     -   MG(4)+Fa(26)*SMG(5)+F(44)*SMG(6)+Fa(27)*SMG(7)+Fa(28)*SMG(8)
     -   +mup3mup4*F(46)*SMG(10)+mup2mup4*F(47)*SMG(11)+mup2mup3*F(46
     -   )*SMG(12)+p2mup4*F(47)*SMG(13)+p2mup3*F(46)*SMG(14)+p1mup2*F
     -   (46)*SMG(15)+F(48)*SMG(16)
       resultg(6) = Fa(29)*SMG(1)+Fa(30)*SMG(2)+Fa(31)*SMG(3)+Fa(32)*S
     -   MG(4)+Fa(33)*SMG(5)+F(102)*SMG(6)+p3mup4*F(103)*SMG(7)+Fa(34
     -   )*SMG(8)+Fa(35)*SMG(9)+mup3mup4*F(103)*SMG(10)+mup2mup4*F(11
     -   0)*SMG(11)+Fa(36)*SMG(12)+p2mup4*F(110)*SMG(13)+Fa(37)*SMG(1
     -   4)+Fa(38)*SMG(15)+F(124)*SMG(16)
       resultg(7) = Fa(39)*SMG(1)+Fa(40)*SMG(2)+Fa(41)*SMG(3)+Fa(42)*S
     -   MG(4)+Fa(43)*SMG(5)+F(173)*SMG(6)+Fa(44)*SMG(7)+Fa(45)*SMG(8
     -   )+Fa(46)*SMG(10)+mup2mup4*F(182)*SMG(11)+mup2mup3*F(183)*SMG
     -   (12)+Fa(47)*SMG(13)+Fa(48)*SMG(14)+p1mup2*F(183)*SMG(15)+F(1
     -   86)*SMG(16)
       resultg(8) = Fa(49)*SMG(1)+Fa(50)*SMG(2)+Fa(51)*SMG(3)+Fa(52)*S
     -   MG(4)+Fa(53)*SMG(5)+F(375)*SMG(6)+Fa(54)*SMG(7)+Fa(55)*SMG(8
     -   )+Fa(56)*SMG(9)+Fa(57)*SMG(10)+Fa(58)*SMG(11)+Fa(59)*SMG(12)
     -   +Fa(60)*SMG(13)+Fa(61)*SMG(14)+Fa(62)*SMG(15)+F(448)*SMG(16)
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(1)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5)+re
     -   sultg(6)+resultg(7)+resultg(8))
       resultgaugeb(1) = Is45*(2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2
     -   mup3*SMG(3))-SMG(6))+Is12s45*(4*((mup2mup4*(p1mup3+p2mup3)+m
     -   up2mup3*p3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2
     -   mup3*p2mup4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2
     -   mup4+p2mup3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1
     -   )-(-2*p1mup2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mu
     -   p2*(p1mup3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*S
     -   MG(7)+(p1mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(1
     -   1)+mup2mup3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG
     -   (15))+SMG(16))
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
       call FaFunctionABEmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2
     -   mup4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,
     -   p5mup3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
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
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*SMG(4
     -   )+Fa(5)*SMG(5)+F(3)*SMG(6)-2*p3mup4*F(2)*SMG(7)+Fa(6)*SMG(8)
     -   -2*mup3mup4*F(2)*SMG(10)+2*mup2mup4*F(2)*SMG(11)-2*mup2mup3*
     -   F(2)*SMG(12)+2*p2mup4*F(2)*SMG(13)-2*p2mup3*F(2)*SMG(14)-2*p
     -   1mup2*F(2)*SMG(15)+F(2)*SMG(16)
       resultg(4) = Fa(7)*SMG(1)+Fa(8)*SMG(2)+Fa(9)*SMG(3)+Fa(10)*SMG(
     -   4)+Fa(11)*SMG(5)+F(16)*SMG(6)+p3mup4*F(17)*SMG(7)+Fa(12)*SMG
     -   (8)+mup3mup4*F(17)*SMG(10)+mup2mup4*F(19)*SMG(11)+mup2mup3*F
     -   (17)*SMG(12)+p2mup4*F(19)*SMG(13)+Fa(13)*SMG(14)+p1mup2*F(17
     -   )*SMG(15)+F(22)*SMG(16)
       resultg(3) = Fa(14)*SMG(1)+Fa(15)*SMG(2)+Fa(16)*SMG(3)+Fa(17)*S
     -   MG(4)+Fa(18)*SMG(5)+F(31)*SMG(6)+p3mup4*F(32)*SMG(7)+Fa(19)*
     -   SMG(8)+Fa(20)*SMG(9)+mup3mup4*F(32)*SMG(10)+mup2mup4*F(33)*S
     -   MG(11)+mup2mup3*F(32)*SMG(12)+p2mup4*F(33)*SMG(13)+p2mup3*F(
     -   32)*SMG(14)+Fa(21)*SMG(15)+F(35)*SMG(16)
       resultg(5) = Fa(22)*SMG(1)+Fa(23)*SMG(2)+Fa(24)*SMG(3)+Fa(25)*S
     -   MG(4)+Fa(26)*SMG(5)+F(44)*SMG(6)+Fa(27)*SMG(7)+Fa(28)*SMG(8)
     -   +mup3mup4*F(46)*SMG(10)+mup2mup4*F(47)*SMG(11)+mup2mup3*F(46
     -   )*SMG(12)+p2mup4*F(47)*SMG(13)+p2mup3*F(46)*SMG(14)+p1mup2*F
     -   (46)*SMG(15)+F(48)*SMG(16)
       resultg(6) = Fa(29)*SMG(1)+Fa(30)*SMG(2)+Fa(31)*SMG(3)+Fa(32)*S
     -   MG(4)+Fa(33)*SMG(5)+F(102)*SMG(6)+p3mup4*F(103)*SMG(7)+Fa(34
     -   )*SMG(8)+Fa(35)*SMG(9)+mup3mup4*F(103)*SMG(10)+mup2mup4*F(11
     -   0)*SMG(11)+Fa(36)*SMG(12)+p2mup4*F(110)*SMG(13)+Fa(37)*SMG(1
     -   4)+Fa(38)*SMG(15)+F(124)*SMG(16)
       resultg(7) = Fa(39)*SMG(1)+Fa(40)*SMG(2)+Fa(41)*SMG(3)+Fa(42)*S
     -   MG(4)+Fa(43)*SMG(5)+F(173)*SMG(6)+Fa(44)*SMG(7)+Fa(45)*SMG(8
     -   )+Fa(46)*SMG(10)+mup2mup4*F(182)*SMG(11)+mup2mup3*F(183)*SMG
     -   (12)+Fa(47)*SMG(13)+Fa(48)*SMG(14)+p1mup2*F(183)*SMG(15)+F(1
     -   86)*SMG(16)
       resultg(8) = Fa(49)*SMG(1)+Fa(50)*SMG(2)+Fa(51)*SMG(3)+Fa(52)*S
     -   MG(4)+Fa(53)*SMG(5)+F(375)*SMG(6)+Fa(54)*SMG(7)+Fa(55)*SMG(8
     -   )+Fa(56)*SMG(9)+Fa(57)*SMG(10)+Fa(58)*SMG(11)+Fa(59)*SMG(12)
     -   +Fa(60)*SMG(13)+Fa(61)*SMG(14)+Fa(62)*SMG(15)+F(448)*SMG(16)
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(1)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5)+re
     -   sultg(6)+resultg(7)+resultg(8))
       resultgaugeb(1) = Is45*(2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2
     -   mup3*SMG(3))-SMG(6))+Is12s45*(4*((mup2mup4*(p1mup3+p2mup3)+m
     -   up2mup3*p3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2
     -   mup3*p2mup4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2
     -   mup4+p2mup3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1
     -   )-(-2*p1mup2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mu
     -   p2*(p1mup3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*S
     -   MG(7)+(p1mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(1
     -   1)+mup2mup3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG
     -   (15))+SMG(16))
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
       call FaFunctionABEmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2
     -   mup4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,
     -   p5mup3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
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
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*SMG(4
     -   )+Fa(5)*SMG(5)+F(3)*SMG(6)-2*p3mup4*F(2)*SMG(7)+Fa(6)*SMG(8)
     -   -2*mup3mup4*F(2)*SMG(10)+2*mup2mup4*F(2)*SMG(11)-2*mup2mup3*
     -   F(2)*SMG(12)+2*p2mup4*F(2)*SMG(13)-2*p2mup3*F(2)*SMG(14)-2*p
     -   1mup2*F(2)*SMG(15)+F(2)*SMG(16)
       resultg(4) = Fa(7)*SMG(1)+Fa(8)*SMG(2)+Fa(9)*SMG(3)+Fa(10)*SMG(
     -   4)+Fa(11)*SMG(5)+F(16)*SMG(6)+p3mup4*F(17)*SMG(7)+Fa(12)*SMG
     -   (8)+mup3mup4*F(17)*SMG(10)+mup2mup4*F(19)*SMG(11)+mup2mup3*F
     -   (17)*SMG(12)+p2mup4*F(19)*SMG(13)+Fa(13)*SMG(14)+p1mup2*F(17
     -   )*SMG(15)+F(22)*SMG(16)
       resultg(3) = Fa(14)*SMG(1)+Fa(15)*SMG(2)+Fa(16)*SMG(3)+Fa(17)*S
     -   MG(4)+Fa(18)*SMG(5)+F(31)*SMG(6)+p3mup4*F(32)*SMG(7)+Fa(19)*
     -   SMG(8)+Fa(20)*SMG(9)+mup3mup4*F(32)*SMG(10)+mup2mup4*F(33)*S
     -   MG(11)+mup2mup3*F(32)*SMG(12)+p2mup4*F(33)*SMG(13)+p2mup3*F(
     -   32)*SMG(14)+Fa(21)*SMG(15)+F(35)*SMG(16)
       resultg(5) = Fa(22)*SMG(1)+Fa(23)*SMG(2)+Fa(24)*SMG(3)+Fa(25)*S
     -   MG(4)+Fa(26)*SMG(5)+F(44)*SMG(6)+Fa(27)*SMG(7)+Fa(28)*SMG(8)
     -   +mup3mup4*F(46)*SMG(10)+mup2mup4*F(47)*SMG(11)+mup2mup3*F(46
     -   )*SMG(12)+p2mup4*F(47)*SMG(13)+p2mup3*F(46)*SMG(14)+p1mup2*F
     -   (46)*SMG(15)+F(48)*SMG(16)
       resultg(6) = Fa(29)*SMG(1)+Fa(30)*SMG(2)+Fa(31)*SMG(3)+Fa(32)*S
     -   MG(4)+Fa(33)*SMG(5)+F(102)*SMG(6)+p3mup4*F(103)*SMG(7)+Fa(34
     -   )*SMG(8)+Fa(35)*SMG(9)+mup3mup4*F(103)*SMG(10)+mup2mup4*F(11
     -   0)*SMG(11)+Fa(36)*SMG(12)+p2mup4*F(110)*SMG(13)+Fa(37)*SMG(1
     -   4)+Fa(38)*SMG(15)+F(124)*SMG(16)
       resultg(7) = Fa(39)*SMG(1)+Fa(40)*SMG(2)+Fa(41)*SMG(3)+Fa(42)*S
     -   MG(4)+Fa(43)*SMG(5)+F(173)*SMG(6)+Fa(44)*SMG(7)+Fa(45)*SMG(8
     -   )+Fa(46)*SMG(10)+mup2mup4*F(182)*SMG(11)+mup2mup3*F(183)*SMG
     -   (12)+Fa(47)*SMG(13)+Fa(48)*SMG(14)+p1mup2*F(183)*SMG(15)+F(1
     -   86)*SMG(16)
       resultg(8) = Fa(49)*SMG(1)+Fa(50)*SMG(2)+Fa(51)*SMG(3)+Fa(52)*S
     -   MG(4)+Fa(53)*SMG(5)+F(375)*SMG(6)+Fa(54)*SMG(7)+Fa(55)*SMG(8
     -   )+Fa(56)*SMG(9)+Fa(57)*SMG(10)+Fa(58)*SMG(11)+Fa(59)*SMG(12)
     -   +Fa(60)*SMG(13)+Fa(61)*SMG(14)+Fa(62)*SMG(15)+F(448)*SMG(16)
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(2)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5)+re
     -   sultg(6)+resultg(7)+resultg(8))
       resultgaugeb(2) = Is45*(2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2
     -   mup3*SMG(3))-SMG(6))+Is12s45*(4*((mup2mup4*(p1mup3+p2mup3)+m
     -   up2mup3*p3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2
     -   mup3*p2mup4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2
     -   mup4+p2mup3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1
     -   )-(-2*p1mup2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mu
     -   p2*(p1mup3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*S
     -   MG(7)+(p1mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(1
     -   1)+mup2mup3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG
     -   (15))+SMG(16))
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
       call FaFunctionABEmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2mup3,p2
     -   mup4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p5mup2,
     -   p5mup3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
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
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+Fa(4)*SMG(4
     -   )+Fa(5)*SMG(5)+F(3)*SMG(6)-2*p3mup4*F(2)*SMG(7)+Fa(6)*SMG(8)
     -   -2*mup3mup4*F(2)*SMG(10)+2*mup2mup4*F(2)*SMG(11)-2*mup2mup3*
     -   F(2)*SMG(12)+2*p2mup4*F(2)*SMG(13)-2*p2mup3*F(2)*SMG(14)-2*p
     -   1mup2*F(2)*SMG(15)+F(2)*SMG(16)
       resultg(4) = Fa(7)*SMG(1)+Fa(8)*SMG(2)+Fa(9)*SMG(3)+Fa(10)*SMG(
     -   4)+Fa(11)*SMG(5)+F(16)*SMG(6)+p3mup4*F(17)*SMG(7)+Fa(12)*SMG
     -   (8)+mup3mup4*F(17)*SMG(10)+mup2mup4*F(19)*SMG(11)+mup2mup3*F
     -   (17)*SMG(12)+p2mup4*F(19)*SMG(13)+Fa(13)*SMG(14)+p1mup2*F(17
     -   )*SMG(15)+F(22)*SMG(16)
       resultg(3) = Fa(14)*SMG(1)+Fa(15)*SMG(2)+Fa(16)*SMG(3)+Fa(17)*S
     -   MG(4)+Fa(18)*SMG(5)+F(31)*SMG(6)+p3mup4*F(32)*SMG(7)+Fa(19)*
     -   SMG(8)+Fa(20)*SMG(9)+mup3mup4*F(32)*SMG(10)+mup2mup4*F(33)*S
     -   MG(11)+mup2mup3*F(32)*SMG(12)+p2mup4*F(33)*SMG(13)+p2mup3*F(
     -   32)*SMG(14)+Fa(21)*SMG(15)+F(35)*SMG(16)
       resultg(5) = Fa(22)*SMG(1)+Fa(23)*SMG(2)+Fa(24)*SMG(3)+Fa(25)*S
     -   MG(4)+Fa(26)*SMG(5)+F(44)*SMG(6)+Fa(27)*SMG(7)+Fa(28)*SMG(8)
     -   +mup3mup4*F(46)*SMG(10)+mup2mup4*F(47)*SMG(11)+mup2mup3*F(46
     -   )*SMG(12)+p2mup4*F(47)*SMG(13)+p2mup3*F(46)*SMG(14)+p1mup2*F
     -   (46)*SMG(15)+F(48)*SMG(16)
       resultg(6) = Fa(29)*SMG(1)+Fa(30)*SMG(2)+Fa(31)*SMG(3)+Fa(32)*S
     -   MG(4)+Fa(33)*SMG(5)+F(102)*SMG(6)+p3mup4*F(103)*SMG(7)+Fa(34
     -   )*SMG(8)+Fa(35)*SMG(9)+mup3mup4*F(103)*SMG(10)+mup2mup4*F(11
     -   0)*SMG(11)+Fa(36)*SMG(12)+p2mup4*F(110)*SMG(13)+Fa(37)*SMG(1
     -   4)+Fa(38)*SMG(15)+F(124)*SMG(16)
       resultg(7) = Fa(39)*SMG(1)+Fa(40)*SMG(2)+Fa(41)*SMG(3)+Fa(42)*S
     -   MG(4)+Fa(43)*SMG(5)+F(173)*SMG(6)+Fa(44)*SMG(7)+Fa(45)*SMG(8
     -   )+Fa(46)*SMG(10)+mup2mup4*F(182)*SMG(11)+mup2mup3*F(183)*SMG
     -   (12)+Fa(47)*SMG(13)+Fa(48)*SMG(14)+p1mup2*F(183)*SMG(15)+F(1
     -   86)*SMG(16)
       resultg(8) = Fa(49)*SMG(1)+Fa(50)*SMG(2)+Fa(51)*SMG(3)+Fa(52)*S
     -   MG(4)+Fa(53)*SMG(5)+F(375)*SMG(6)+Fa(54)*SMG(7)+Fa(55)*SMG(8
     -   )+Fa(56)*SMG(9)+Fa(57)*SMG(10)+Fa(58)*SMG(11)+Fa(59)*SMG(12)
     -   +Fa(60)*SMG(13)+Fa(61)*SMG(14)+Fa(62)*SMG(15)+F(448)*SMG(16)
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(3)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5)+re
     -   sultg(6)+resultg(7)+resultg(8))
       resultgaugeb(3) = Is45*(2*(mup3mup4*SMG(1)-mup2mup4*SMG(2)+mup2
     -   mup3*SMG(3))-SMG(6))+Is12s45*(4*((mup2mup4*(p1mup3+p2mup3)+m
     -   up2mup3*p3mup4)*SMG(4)+(mup3mup4*p1mup2+mup2mup4*p2mup3-mup2
     -   mup3*p2mup4)*SMG(5))+(p2sq+p3sq-s23)*SMG(6)-2*((2*(p1mup3*p2
     -   mup4+p2mup3*(p2mup4+p3mup4))+mup3mup4*(p2sq+p3sq-s23))*SMG(1
     -   )-(-2*p1mup2*p3mup4+mup2mup4*(p2sq+p3sq-s23))*SMG(2)+(2*p1mu
     -   p2*(p1mup3+p2mup3)+mup2mup3*(p2sq+p3sq-s23))*SMG(3)+p3mup4*S
     -   MG(7)+(p1mup3+p2mup3)*SMG(8)+mup3mup4*SMG(10)-mup2mup4*SMG(1
     -   1)+mup2mup3*SMG(12)-p2mup4*SMG(13)+p2mup3*SMG(14)+p1mup2*SMG
     -   (15))+SMG(16))
           endif      
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
c************** Calling the Fa functions**********************************************************************
       subroutine FaFunctionABEmmm(p1mup2, p1mup3, p1mup4, p2mup2, p2m
     -   up3,p2mup4,p3mup2, p3mup3, p3mup4 , p4mup2, p4mup3, p4mup4,p
     -   5mup2,p5mup3, p5mup4, mup2mup4, mup2mup3, mup3mup4,Fa)
       IMPLICIT NONE
      Complex*16   p1mup2, p1mup3, p1mup4, p2mup2, p2mup3, p2mup4, 
     -          p3mup2, p3mup3, p3mup4, p4mup2, p4mup3, p4mup4, 
     -          p5mup2, p5mup3, p5mup4
       Complex*16 mup2mup3,mup2mup4,mup3mup4
        common/invarianetsmmm/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,Is12
      Real*8 p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,Is12
       common/Ffunctionsmmm/F
       common/Pfunctionsmmm/P
       common/Kfunctionsmmm/K
      Complex*16 Fa(62),F(448),K(4)
       Real*8 P(554) 
        Fa(1) = 2*mup3mup4*F(1)-4*(p1mup3*p2mup4*F(2)+p2mup3*(p2mup4*F(
     -   2)+p3mup4*F(2)))
       Fa(2) = -4*p1mup2*p3mup4*F(2)+2*mup2mup4*F(3)
       Fa(3) = 2*mup2mup3*F(1)-4*p1mup2*(p1mup3*F(2)+p2mup3*F(2))
       Fa(4) = 4*(mup2mup3*p3mup4*F(2)+mup2mup4*(p1mup3*F(2)+p2mup3*F(
     -   2)))
       Fa(5) = 4*(mup3mup4*p1mup2*F(2)+mup2mup4*p2mup3*F(2)-mup2mup3*p
     -   2mup4*F(2))
       Fa(6) = -2*(p1mup3*F(2)+p2mup3*F(2))
       Fa(7) = mup3mup4*F(4)+p1mup3*(p2mup4*F(5)+p3mup4*F(6))+p2mup3*(
     -   p2mup4*F(5)+p3mup4*F(7))-8*p3mup3*(p2mup4*F(8)+p3mup4*F(9))
       Fa(8) = p1mup2*p3mup4*F(10)+mup2mup4*F(11)
       Fa(9) = mup2mup3*F(4)+p1mup2*(p1mup3*F(5)+p2mup3*F(5))-8*p1mup2
     -   *p3mup3*F(8)
       Fa(10) = 8*mup2mup4*p3mup3*F(8)+mup2mup3*p3mup4*F(12)+mup2mup4*
     -   (p1mup3*F(13)+p2mup3*F(13))
       Fa(11) = 8*mup2mup4*p3mup3*F(9)+mup2mup3*p2mup4*F(10)+mup3mup4*
     -   p1mup2*F(12)+mup2mup4*(p1mup3*F(14)+p2mup3*F(15))
       Fa(12) = -4*p3mup3*F(8)+p1mup3*F(18)+p2mup3*F(18)
       Fa(13) = -4*p3mup3*F(9)+p1mup3*F(20)+p2mup3*F(21)
       Fa(14) = p1mup3*p2mup4*F(23)+p2mup3*(p2mup4*F(23)+p3mup4*F(23))
     -   +mup3mup4*F(24)
       Fa(15) = mup2mup4*F(25)+p1mup2*(p3mup4*F(26)-8*p2mup4*F(27))-8*
     -   p2mup2*(p2mup4*F(28)+p3mup4*F(28))
       Fa(16) = mup2mup3*F(24)+p1mup2*(p2mup3*F(23)+p1mup3*F(26))-8*p1
     -   mup3*p2mup2*F(28)
       Fa(17) = 8*(mup3mup4*p1mup2*F(27)+mup3mup4*p2mup2*F(28))+mup2mu
     -   p3*p3mup4*F(29)+mup2mup4*(p1mup3*F(29)+p2mup3*F(29))
       Fa(18) = mup2mup3*p2mup4*F(23)+8*mup3mup4*p2mup2*F(28)+mup2mup4
     -   *p2mup3*F(29)+mup3mup4*p1mup2*F(30)
       Fa(19) = p1mup3*F(32)+p2mup3*F(32)
       Fa(20) = -4*(p1mup2*F(27)+p2mup2*F(28))
       Fa(21) = -4*p2mup2*F(28)+p1mup2*F(34)
       Fa(22) = p1mup3*p2mup4*F(38)+p2mup3*(-8*(p1mup4*F(36)-p5mup4*F(
     -   37))+p2mup4*F(39)+p3mup4*F(39))+mup3mup4*F(40)
       Fa(23) = p1mup2*(-8*(p1mup4*F(36)+p2mup4*F(36)-p5mup4*F(37))+p3
     -   mup4*F(39))+mup2mup4*F(41)
       Fa(24) = p1mup2*(p1mup3*F(38)+p2mup3*F(38))+mup2mup3*F(40)
       Fa(25) = mup2mup4*(p1mup3*F(42)+p2mup3*F(42))+mup2mup3*(8*(p1mu
     -   p4*F(36)+p2mup4*F(36)-p5mup4*F(37))+p3mup4*F(43))
       Fa(26) = mup2mup3*p2mup4*F(38)+mup3mup4*p1mup2*F(42)+mup2mup4*p
     -   2mup3*F(42)
       Fa(27) = -4*(p1mup4*F(36)+p2mup4*F(36)-p5mup4*F(37))+p3mup4*F(4
     -   5)
       Fa(28) = p1mup3*F(46)+p2mup3*F(46)
       Fa(29) = p2mup3*(p2mup4*F(49)+p3mup4*F(50))+p3mup3*(p2mup4*F(52
     -   )+p3mup4*F(53))+p1mup3*(p2mup4*F(51)+p3mup4*F(54))+mup3mup4*
     -   F(55)
       Fa(30) = p2mup2*(p2mup4*F(56)+p3mup4*F(57))+p3mup2*(p3mup4*F(59
     -   )+p2mup4*F(60))+p1mup2*(p2mup4*F(58)+p3mup4*F(61))+mup2mup4*
     -   F(62)
       Fa(31) = mup2mup3*F(63)+p1mup2*(p1mup3*F(64)+p2mup3*F(65))+p2mu
     -   p2*(p1mup3*F(66)+p2mup3*F(67)+8*p3mup3*F(70))+p3mup2*(p1mup3
     -   *F(69)+8*p2mup3*F(71))+p3mup3*(p1mup2*F(68)+8*p3mup2*F(72))
       Fa(32) = mup2mup3*p3mup4*F(73)+mup2mup4*(p2mup3*F(74)+p1mup3*F(
     -   75))+mup3mup4*(p1mup2*F(78)+p3mup2*F(79))+p1mup2*(p1mup3*p3m
     -   up4*F(80)+p2mup3*p3mup4*F(81))+p2mup2*(mup3mup4*F(77)+p1mup3
     -   *p3mup4*F(82)+p2mup3*p3mup4*F(83)-16*p3mup3*p3mup4*F(86))+p3
     -   mup2*(p1mup3*p3mup4*F(85)-16*p2mup3*p3mup4*F(87))+p3mup3*(mu
     -   p2mup4*F(76)+p1mup2*p3mup4*F(84)-16*p3mup2*p3mup4*F(88))
       Fa(33) = mup2mup3*p2mup4*F(89)+mup2mup4*(p2mup3*F(90)+p1mup3*F(
     -   92))+mup3mup4*(p3mup2*F(94)+p1mup2*F(95))+p1mup2*(p1mup3*p2m
     -   up4*F(96)+p2mup3*p2mup4*F(97))+p2mup2*(16*p2mup4*p3mup3*F(86
     -   )+mup3mup4*F(93)+p1mup3*p2mup4*F(98)+p2mup3*p2mup4*F(99))+p3
     -   mup3*(16*p2mup4*p3mup2*F(88)+mup2mup4*F(91)+p1mup2*p2mup4*F(
     -   100))+p3mup2*(16*p2mup3*p2mup4*F(87)+p1mup3*p2mup4*F(101))
       Fa(34) = p2mup3*F(104)+p1mup3*F(105)+p3mup3*F(106)
       Fa(35) = p2mup2*F(107)+p1mup2*F(108)+p3mup2*F(109)
       Fa(36) = mup2mup3*F(111)+p1mup2*(p1mup3*F(112)+p2mup3*F(113))+p
     -   2mup2*(8*p3mup3*F(86)+p1mup3*F(114)+p2mup3*F(115))+p3mup3*(8
     -   *p3mup2*F(88)+p1mup2*F(116))+p3mup2*(8*p2mup3*F(87)+p1mup3*F
     -   (117))
       Fa(37) = p2mup3*F(118)+p3mup3*F(119)+p1mup3*F(120)
       Fa(38) = p2mup2*F(121)+p3mup2*F(122)+p1mup2*F(123)
       Fa(39) = p2mup4*p5mup3*F(132)+mup3mup4*F(133)+p1mup3*(p1mup4*F(
     -   125)+p5mup4*F(128)+p2mup4*F(130)+8*p3mup4*F(136))+p3mup3*(p2
     -   mup4*F(127)-8*(p1mup4*F(134)+p3mup4*F(137)-p5mup4*F(138)))-8
     -   *(p1mup4*p5mup3*F(135)-p3mup4*p5mup3*F(139)-p5mup3*p5mup4*F(
     -   140))+p2mup3*(p1mup4*F(126)+p2mup4*F(129)+p3mup4*F(131)+p5mu
     -   p4*F(141))
       Fa(40) = mup2mup4*F(144)+p1mup2*(p1mup4*F(142)+p2mup4*F(142)+p3
     -   mup4*F(143)+p5mup4*F(145))
       Fa(41) = p1mup2*p3mup3*F(146)+p1mup2*(p1mup3*F(147)+p2mup3*F(14
     -   7)+p5mup3*F(148))+mup2mup3*F(149)
       Fa(42) = mup2mup4*p3mup3*F(150)+mup2mup4*(p1mup3*F(151)+p2mup3*
     -   F(151)+p5mup3*F(152))+mup2mup3*(p1mup4*F(153)+p2mup4*F(153)+
     -   p3mup4*F(154)+p5mup4*F(155))
       Fa(43) = mup3mup4*p1mup2*F(164)+mup2mup3*(p2mup4*F(165)+8*(p1mu
     -   p4*F(166)+p3mup4*F(167)+p5mup4*F(168)))+p3mup3*(-16*p1mup2*(
     -   p1mup4*F(157)+p2mup4*F(157)+p3mup4*F(160)-p5mup4*F(162))+8*m
     -   up2mup4*F(171))+mup2mup4*(p2mup3*F(169)-8*(p1mup3*F(170)+p5m
     -   up3*F(172)))+p1mup2*(p1mup3*(p1mup4*F(156)+p2mup4*F(156)-16*
     -   (p3mup4*F(159)-p5mup4*F(161)))+p2mup3*(p1mup4*F(156)+p2mup4*
     -   F(156)-16*(p3mup4*F(159)-p5mup4*F(161)))+16*(p1mup4*p5mup3*F
     -   (158)+p2mup4*p5mup3*F(158)+p5mup3*p5mup4*F(163)-Is12*p3mup4*
     -   p5mup3*K(1)))
       Fa(44) = p1mup4*F(174)+p2mup4*F(174)+p3mup4*F(175)+p5mup4*F(176
     -   )
       Fa(45) = p3mup3*F(177)+p1mup3*F(178)+p2mup3*F(178)+p5mup3*F(179
     -   )
       Fa(46) = p1mup3*(8*(p3mup4*F(159)-p5mup4*F(161))+p1mup4*F(180)+
     -   p2mup4*F(180))+p2mup3*(8*(p3mup4*F(159)-p5mup4*F(161))+p1mup
     -   4*F(180)+p2mup4*F(180))+mup3mup4*F(181)+8*(-(p1mup4*p5mup3*F
     -   (158))-p2mup4*p5mup3*F(158)+p3mup3*(p1mup4*F(157)+p2mup4*F(1
     -   57)+p3mup4*F(160)-p5mup4*F(162))-p5mup3*p5mup4*F(163)+Is12*p
     -   3mup4*p5mup3*K(1))
       Fa(47) = -4*(p1mup4*F(166)+p3mup4*F(167)+p5mup4*F(168))+p2mup4*
     -   F(184)
       Fa(48) = 4*(p1mup3*F(170)-p3mup3*F(171)+p5mup3*F(172))+p2mup3*F
     -   (185)
       Fa(49) = p3mup3*(p2mup4*F(189)+p3mup4*F(190)+p5mup4*F(191)+p1mu
     -   p4*F(193))+p1mup4*p5mup3*F(197)+p1mup3*(p1mup4*F(187)+p3mup4
     -   *F(192)+p2mup4*F(194)+p5mup4*F(198))+p2mup4*p5mup3*F(199)+p2
     -   mup3*(p1mup4*F(188)+p2mup4*F(195)+p3mup4*F(196)+p5mup4*F(200
     -   ))+p3mup4*p5mup3*F(201)+mup3mup4*F(202)+p5mup3*p5mup4*F(203)
       Fa(50) = p2mup2*(p1mup4*F(205)+p2mup4*F(207)+p3mup4*F(210)+p5mu
     -   p4*F(212))+p3mup2*(p2mup4*F(208)+p3mup4*F(211)+p5mup4*F(213)
     -   +p1mup4*F(216))+p1mup2*(p1mup4*F(204)+p2mup4*F(206)+p3mup4*F
     -   (217)+p5mup4*F(218))+p5mup2*(p2mup4*F(209)+p5mup4*F(214)+p3m
     -   up4*F(215)+p1mup4*F(219))+mup2mup4*F(220)
       Fa(51) = p2mup2*(p2mup3*F(221)+p3mup3*F(225)+p5mup3*F(229)+p1mu
     -   p3*F(230))+p3mup3*(p1mup2*F(224)+p3mup2*F(226)+p5mup2*F(231)
     -   )+p3mup2*(p2mup3*F(222)+p5mup3*F(227)+p1mup3*F(232))+p5mup2*
     -   (p2mup3*F(223)+p5mup3*F(228)+p1mup3*F(234))+p1mup2*(p2mup3*F
     -   (233)+p1mup3*F(235)+p5mup3*F(236))+mup2mup3*F(237)
       Fa(52) = mup3mup4*(p1mup2*F(238)+p3mup2*F(240)+p5mup2*F(241))+p
     -   2mup2*(mup3mup4*F(239)+p2mup3*(p2mup4*F(250)+p3mup4*F(254)+p
     -   5mup4*F(258)-16*p1mup4*F(260))+p1mup3*(p1mup4*F(243)-16*(p2m
     -   up4*F(260)+p3mup4*F(278)-p5mup4*F(289)))-16*(p1mup4*p5mup3*F
     -   (266)+p2mup4*p5mup3*F(271)-p3mup4*p5mup3*F(286)-p5mup3*p5mup
     -   4*F(292)+p3mup3*(p1mup4*F(263)+p2mup4*F(273)+p3mup4*F(283)-p
     -   5mup4*F(294))))+p1mup2*(p1mup4*p5mup3*F(248)+p3mup4*p5mup3*F
     -   (255)+p1mup3*(p1mup4*F(242)+p2mup4*F(246)+p3mup4*F(251)+p5mu
     -   p4*F(256))+p2mup3*(p1mup4*F(246)+p2mup4*F(249)+p3mup4*F(253)
     -   +p5mup4*F(257))-16*(p2mup4*p5mup3*F(276)-p5mup3*p5mup4*F(297
     -   )))+p3mup2*(p2mup3*(p5mup4*F(259)+16*(p1mup4*F(261)-p2mup4*F
     -   (270)-p3mup4*F(280)))+p1mup3*(p1mup4*F(244)-16*(p2mup4*F(269
     -   )-p3mup4*F(279)-p5mup4*F(290)))+16*(p1mup4*p5mup3*F(267)+p2m
     -   up4*p5mup3*F(277)+p3mup4*p5mup3*F(287)+p5mup3*p5mup4*F(298))
     -   )+mup2mup3*(p1mup4*F(300)+p2mup4*F(301)+p5mup4*F(302)+p3mup4
     -   *F(303))+p3mup3*(p1mup2*(p1mup4*F(247)-16*(p2mup4*F(272)+p3m
     -   up4*F(282)-p5mup4*F(293)))-16*(p3mup2*(p1mup4*F(264)+p2mup4*
     -   F(274)+p3mup4*F(284)-p5mup4*F(295))+p5mup2*(p1mup4*F(265)+p2
     -   mup4*F(275)+p3mup4*F(285)+p5mup4*F(296)))+mup2mup4*F(304))+m
     -   up2mup4*(p1mup3*F(305)+p2mup3*F(306)+p5mup3*F(307))+p5mup2*(
     -   p1mup3*(p1mup4*F(245)+p3mup4*F(252)-16*(p2mup4*F(266)-p5mup4
     -   *F(291)))+16*(p1mup4*p5mup3*F(268)+p3mup4*p5mup3*F(288)+p2mu
     -   p3*(p1mup4*F(262)-p2mup4*F(271)-p3mup4*F(281)+p5mup4*F(292))
     -   +p5mup3*p5mup4*F(299)-p2mup4*p5mup3*K(2)))
       Fa(53) = p1mup2*(p1mup4*p5mup3*F(315)+p2mup4*p5mup3*F(319)+p1mu
     -   p3*(p1mup4*F(308)+p2mup4*F(316)+p3mup4*F(320)+p5mup4*F(322))
     -   +p2mup3*(p1mup4*F(312)+p2mup4*F(318)+p5mup4*F(323)-16*p3mup4
     -   *F(343))-16*(p3mup4*p5mup3*F(328)-p5mup3*p5mup4*F(359)))+p2m
     -   up2*(p2mup3*(p2mup4*F(254)+p1mup4*F(313)+p3mup4*F(321)+p5mup
     -   4*F(324))+p1mup3*(p1mup4*F(309)-16*(p2mup4*F(332)+p3mup4*F(3
     -   40)-p5mup4*F(350)))-16*(p1mup4*p5mup3*F(329)+p2mup4*p5mup3*F
     -   (337)+p3mup4*p5mup3*F(349)+p3mup3*(p2mup4*F(283)+p1mup4*F(32
     -   6)+p3mup4*F(346)-p5mup4*F(356))+p5mup3*p5mup4*F(360))+mup3mu
     -   p4*F(363))+mup3mup4*(p3mup2*F(364)+p5mup2*F(365)+p1mup2*F(36
     -   6))+mup2mup3*(p3mup4*F(367)+p1mup4*F(368)+p2mup4*F(369)+p5mu
     -   p4*F(370))+p3mup3*(p1mup2*(p1mup4*F(310)-16*(p2mup4*F(335)+p
     -   3mup4*F(341)-p5mup4*F(355)))-16*(p3mup2*(p2mup4*F(284)+p1mup
     -   4*F(327)+p3mup4*F(347)-p5mup4*F(357))+p5mup2*(p1mup4*F(328)-
     -   p2mup4*F(336)+p3mup4*F(348)-p5mup4*F(358)))+mup2mup4*F(371))
     -   +mup2mup4*(p1mup3*F(372)+p2mup3*F(373)+p5mup3*F(374))+p3mup2
     -   *(p1mup3*(p1mup4*F(310)-16*(p2mup4*F(333)+p3mup4*F(341)-p5mu
     -   p4*F(351)))-16*(-(p1mup4*p5mup3*F(330))+p2mup4*p5mup3*F(338)
     -   +p2mup3*(p2mup4*F(280)+p1mup4*F(325)+p3mup4*F(344)-p5mup4*F(
     -   353))-p5mup3*p5mup4*F(361)+p3mup4*p5mup3*K(3)))+p5mup2*(p1mu
     -   p3*(p1mup4*F(311)+p2mup4*F(317)-16*(p3mup4*F(342)-p5mup4*F(3
     -   52)))+p2mup3*(p1mup4*F(314)+16*(p2mup4*F(334)+p3mup4*F(345)+
     -   p5mup4*F(354)))+16*(p1mup4*p5mup3*F(331)-p2mup4*p5mup3*F(339
     -   )+p5mup3*p5mup4*F(362)-p3mup4*p5mup3*K(4)))
       Fa(54) = p1mup4*F(376)+p2mup4*F(377)+p5mup4*F(378)+p3mup4*F(379
     -   )
       Fa(55) = p3mup3*F(380)+p2mup3*F(381)+p1mup3*F(382)+p5mup3*F(383
     -   )
       Fa(56) = p1mup2*F(384)+p2mup2*F(385)+p3mup2*F(386)+p5mup2*F(387
     -   )
       Fa(57) = p1mup4*p5mup3*F(389)+p2mup4*p5mup3*F(390)+p3mup4*p5mup
     -   3*F(391)+p1mup3*(p5mup4*F(392)+8*(p1mup4*F(394)+p2mup4*F(396
     -   )+p3mup4*F(399)))+p2mup3*(p5mup4*F(393)-8*(p1mup4*F(395)+p2m
     -   up4*F(397)+p3mup4*F(400)))+p3mup3*(p1mup4*F(388)-8*(p2mup4*F
     -   (398)+p3mup4*F(401)+p5mup4*F(402)))-8*p5mup3*p5mup4*F(403)+m
     -   up3mup4*F(404)
       Fa(58) = p1mup2*(p5mup4*F(410)-8*(p1mup4*F(411)+p2mup4*F(412)+p
     -   3mup4*F(414)))+p2mup2*(p1mup4*F(405)+8*(p3mup4*F(397)+p2mup4
     -   *F(413)+p5mup4*F(415)))+p3mup2*(p1mup4*F(406)+8*(p2mup4*F(39
     -   7)+p3mup4*F(400)+p5mup4*F(416)))+p5mup2*(p1mup4*F(407)+p2mup
     -   4*F(408)+p3mup4*F(409)+8*p5mup4*F(417))+mup2mup4*F(418)
       Fa(59) = p1mup2*(p5mup3*F(424)+8*(p1mup3*F(425)+p2mup3*F(426)))
     -   +p3mup3*(p1mup2*F(423)-8*(p3mup2*F(430)-p5mup2*F(431)))+p2mu
     -   p2*(p1mup3*F(419)-8*(p2mup3*F(427)+p3mup3*F(429)+p5mup3*F(43
     -   2)))+p3mup2*(p1mup3*F(420)+8*(p2mup3*F(428)-p5mup3*F(433)))+
     -   p5mup2*(p1mup3*F(421)+p2mup3*F(422)-8*p5mup3*F(434))+mup2mup
     -   3*F(435)
       Fa(60) = p3mup4*F(436)+p1mup4*F(437)+p2mup4*F(438)+p5mup4*F(439
     -   )
       Fa(61) = p3mup3*F(440)+p1mup3*F(441)+p2mup3*F(442)+p5mup3*F(443
     -   )
       Fa(62) = p2mup2*F(444)+p3mup2*F(445)+p5mup2*F(446)+p1mup2*F(447
     -   )
       Return
       End
