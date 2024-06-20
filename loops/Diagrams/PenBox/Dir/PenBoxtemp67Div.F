c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:4/6/2012
c ************************************************************************************
c determine the  finite virtual corrections of 
c psi(p6) ---->--$$$$photon$$$$$$--->---   psi(p5)
c                   |              |              
c                   |              |$$$$$$$$$$$$  p4,mu_p4 
c                   |              |             
c                   |              |$$$$$$$$$$$$  p3, mu_p3 
c                   |              |             
c                   |              |              
c                   |              |             
c barpsi(p1)-->-$$$photon$$$$$$------>---   bar_psi(p2)
c                                               
c Note: To make it shorter in the promgram: mu_p3,...->mup3,... 
c Notation of External momenta: p1+p2+p3+p4+p5+p6=0 
c mu_p3, should be think as external current 
c alpha is the helicity of the initial spinor 
c musq is the renormalization scale energy  
c comp: integer value.The first time called with p1...p5, comp=1
c ATTENTION: ONLY!!!If you have to call the subroutine consecutively with the same arguments
c(p1,p2,p3,p4,p5). Then, comp=-1 (it safes 4000 lines of code) 
c This applies when you have for examples the same diagram for an off-shell photon
c and a Z boson. The differences are in the coupling and  the part that depends on the
c polarization vector that are calculated at the end of this program.
c************************************************************************************
c************************************************************************************
       subroutine PenBox67Div(p1,p2,p3,p4,p5,barpsi_p1,psi_p5,barpsi_p
     &   2,psi_p4,mup3,musq,comp,ngluon,posgluon,result,resultn,resul
     &   tborn,Div)
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
      Complex*16 SMB(39),F(86)
       Complex*16 SMB0(8),SMB11(0:3),SMB12(0:3)
      Complex*16 SMB1(13)

      Complex*16 barpsi_p1_P(2),psi_p5_P(2),barpsi_p2_P(2),psi_p4_P(2) 
       Complex*16 barpsi_p1_M(2),psi_p5_M(2),barpsi_p2_M(2),psi_p4_M(2)
       Complex*16 barpsi_p1(4),psi_p5(4),barpsi_p2(4),psi_p4(4)
       Real*8 P(48)
       Complex*16  SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr
       EXTERNAL    SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr 
       Real*8 delta
       External  delta
       Complex*16 ten2cc
       External  ten2cc
       Complex*16 v1(0:3) 
       Integer i,comp
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c    Declaration of variables 
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16   mup3(0:3)
       Complex*16   p1mup3, p2mup3, p3mup3, p4mup3, p5mup3
       Real*8 dotrr
       Complex*16 B0finDiv,C0finDiv,D0finDiv,E0fin
       EXTERNAL dotrr,B0finDiv,C0finDiv,D0finDiv,E0fin
        Real*8   p1sq, p1p2, p1p3, p1p4, p1p5 
       Real*8   p2sq, p2p3, p2p4, p2p5 
       Real*8   p3sq, p3p4, p3p5 
       Real*8   p4sq, p4p5 
       Real*8   p5sq 
       Real*8   s12, s13, s14, s15 
       Real*8   s23, s24, s25 
       Real*8   s34, s35 
       Real*8   s45 
       Complex*16  B012,B013,B014,B015 
       Complex*16  B023,B024,B025 
       Complex*16  B034,B035 
       Complex*16  B045 
       Real*8  B012R,B013R,B014R,B015R 
       Real*8  B023R,B024R,B025R 
       Real*8  B034R,B035R 
       Real*8  B045R 
       Real*8  B012I,B013I,B014I,B015I 
       Real*8  B023I,B024I,B025I 
       Real*8  B034I,B035I 
       Real*8  B045I 
       Complex*16 C0123,C0124,C0125 
       Complex*16 C0134,C0135 
       Complex*16 C0145 
       Complex*16 C0234,C0235 
       Complex*16 C0245 
       Complex*16 C0345 
       Real*8 C0123R,C0124R,C0125R 
       Real*8 C0134R,C0135R 
       Real*8 C0145R 
       Real*8 C0234R,C0235R 
       Real*8 C0245R 
       Real*8 C0345R 
       Real*8 C0123I,C0124I,C0125I 
       Real*8 C0134I,C0135I 
       Real*8 C0145I 
       Real*8 C0234I,C0235I 
       Real*8 C0245I 
       Real*8 C0345I 
       Real*8 Cij123R(4,2),Cij124R(4,2),Cij125R(4,2) 
       Real*8 Cij134R(4,2),Cij135R(4,2) 
       Real*8 Cij145R(4,2) 
       Real*8 Cij234R(4,2),Cij235R(4,2) 
       Real*8 Cij245R(4,2) 
       Real*8 Cij345R(4,2) 
       Real*8 Cij123I(4,2),Cij124I(4,2),Cij125I(4,2) 
       Real*8 Cij134I(4,2),Cij135I(4,2) 
       Real*8 Cij145I(4,2) 
       Real*8 Cij234I(4,2),Cij235I(4,2) 
       Real*8 Cij245I(4,2) 
       Real*8 Cij345I(4,2) 
       Complex*16  D01234,D01235 
       Complex*16 D01245 
       Complex*16 D01345 
       Complex*16 D02345 
       Real*8 D01234R,D01235R 
       Real*8 D01245R 
       Real*8 D01345R 
       Real*8 D02345R 
       Real*8 D01234I,D01235I 
       Real*8 D01245I 
       Real*8 D01345I 
       Real*8 D02345I 
       Real*8 Dij1234R(13,3),Dij1235R(13,3) 
       Real*8 Dij1245R(13,3) 
       Real*8 Dij1345R(13,3) 
       Real*8 Dij2345R(13,3) 
       Real*8 Dij1234I(13,3),Dij1235I(13,3) 
       Real*8 Dij1245I(13,3) 
       Real*8 Dij1345I(13,3) 
       Real*8 Dij2345I(13,3) 
       Complex*16 EE0 
       Real*8 EE0R 
       Real*8 EE0I 
       Real*8 EijR(46,4) 
       Real*8 EijI(46,4) 
       Real*8 Invs23MU,Invs34MU 
       Logical PrintB,PrintC,PrintD,PrintE 
      Complex*16 dotrc,dotcc,result(3),resultn,resultborn
       Real*8 musq
      EXTERNAL   dotrc,dotcc

       COMMON/PenBoxFhlFunctions/F
      Save/PenBoxFhlFunctions/
       COMMON/PenBoxPFunctions/P
      Save/PenBoxPFunctions/
       COMMON/PenBoxInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s4
     &   5,s15
       COMMON/PenBoxEVALFunctions/ EE0R,EijR,EE0I,EijI 
       COMMON/PenBoxDVALFunctions/ D01234R,Dij1234R,D01234I,Dij1234I, 
     &   D01235R,Dij1235R,D01235I,Dij1235I, D01245R,Dij1245R,D01245I,
     &   Dij1245I, D01345R,Dij1345R,D01345I,Dij1345I, D02345R,Dij2345
     &   R,D02345I,Dij2345I 
       COMMON/PenBoxCVALFunctions/ C0123R,Cij123R,C0123I,Cij123I, C012
     &   4R,Cij124R,C0124I,Cij124I, C0125R,Cij125R,C0125I,Cij125I, C0
     &   134R,Cij134R,C0134I,Cij134I, C0135R,Cij135R,C0135I,Cij135I, 
     &   C0145R,Cij145R,C0145I,Cij145I, C0234R,Cij234R,C0234I,Cij234I
     &   , C0235R,Cij235R,C0235I,Cij235I, C0245R,Cij245R,C0245I,Cij24
     &   5I, C0345R,Cij345R,C0345I,Cij345I 
       COMMON/PenBoxBVALFunctions/ B012R,B012I, B013R,B013I, B014R,B01
     &   4I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,B025I, B034
     &   R,B034I, B035R,B035I, B045R,B045I
       Integer ngluon, posgluon,Div
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the scalar products. Not inlcueded the contraction of the
c       moments with the external currents  
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
      Invs23MU=1d0/s23
        Invs34MU=1d0/s34
      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False.
c************************************************************************************
c************************************************************************************
        if (comp.gt.0) then    
c    Calling C_ij,D_ij Functions    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       B012=B0finDiv(p1sq,musq,Div)   
       B013=B0finDiv(s12,musq,Div)
       B014=B0finDiv(s45,musq,Div)   
       B015=B0finDiv(p5sq,musq,Div) 
       B023=B0finDiv(p2sq,musq,Div) 
       B024=B0finDiv(s23,musq,Div) 
       B025=B0finDiv(s15,musq,Div)   
       B034=B0finDiv(p3sq,musq,Div) 
       B035=B0finDiv(s34,musq,Div) 
       B045=B0finDiv(p4sq,musq,Div)
       B012R=Dble(B012)          
       B013R=Dble(B013) 
       B014R=Dble(B014)   
       B015R=Dble(B015)           
       B023R=Dble(B023 ) 
       B024R=Dble(B024)  
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
       C0123=C0finDiv(p1sq,p2sq,s12,musq,Div)
       C0124=C0finDiv(p1sq,s23,s45,musq,Div)
       C0125=C0finDiv(p1sq,s15,p5sq,musq,Div)
       C0134=C0finDiv(s12,p3sq,s45,musq,Div)
       C0135=C0finDiv(s12,s34,p5sq,musq,Div)
       C0145=C0finDiv(p4sq,p5sq,s45,musq,Div)
       C0234=C0finDiv(p2sq,p3sq,s23,musq,Div)
       C0235=C0finDiv(p2sq,s34,s15,musq,Div)
       C0245=C0finDiv(s23,p4sq,s15,musq,Div)
       C0345=C0finDiv(p3sq,p4sq,s34,musq,Div)
c************************************************************************************
c************************************************************************************
       call tens_red3_new_Re_Com_Div(p1sq,p2sq,s12, B023,B013,B012, C0
     &   123,C0123R,C0123I,Cij123R,Cij123I) 
       call tens_red3_new_Re_Com_Div(p1sq,s23,s45, B024,B014,B012, C01
     &   24,C0124R,C0124I,Cij124R,Cij124I) 
       call tens_red3_new_Re_Com_Div(p1sq,s15,p5sq, B025,B015,B012, C0
     &   125,C0125R,C0125I,Cij125R,Cij125I) 
       call tens_red3_new_Re_Com_Div(s12,p3sq,s45, B034,B014,B013, C01
     &   34,C0134R,C0134I,Cij134R,Cij134I) 
       call tens_red3_new_Re_Com_Div(s12,s34,p5sq, B035,B015,B013, C01
     &   35,C0135R,C0135I,Cij135R,Cij135I) 
       call tens_red3_new_Re_Com_Div(s45,p4sq,p5sq, B045,B015,B014, C0
     &   145,C0145R,C0145I,Cij145R,Cij145I) 
       call tens_red3_new_Re_Com_Div(p2sq,p3sq,s23, B034,B024,B023, C0
     &   234,C0234R,C0234I,Cij234R,Cij234I) 
       call tens_red3_new_Re_Com_Div(p2sq,s34,s15, B035,B025,B023, C02
     &   35,C0235R,C0235I,Cij235R,Cij235I) 
       call tens_red3_new_Re_Com_Div(s23,p4sq,s15, B045,B025,B024, C02
     &   45,C0245R,C0245I,Cij245R,Cij245I) 
       call tens_red3_new_Re_Com_Div(p3sq,p4sq,s34, B045,B035,B034, C0
     &   345,C0345R,C0345I,Cij345R,Cij345I) 
c************************************************************************************
c************************************************************************************
       D02345=D0finDiv(s23,s34,p2sq,p3sq,p4sq,s15,musq,Div)
       D01345=D0finDiv(s45,s34,s12,p3sq,p4sq,p5sq,musq,Div)
       D01245=D0finDiv(s45,s15,p1sq,s23,p4sq,p5sq,musq,Div)
       D01235=D0finDiv(s12,s15,p1sq,p2sq,s34,p5sq,musq,Div)
       D01234=D0finDiv(s12,s23,p1sq,p2sq,p3sq,s45,musq,Div)
c$$$       Print*," D02345",D02345
c$$$       Print*," D01345",D01345
c$$$       Print*," D01245",D01245
c$$$       Print*," D01235",D01235
c$$$       Print*," D01234",D01234
c************************************************************************************
c************************************************************************************
       call tens_red4_new_Re_Com_Div(p2sq,p3sq,p4sq,p2p3,p2p4,p3p4,C03
     &   45R,C0245R,C0235R,C0234R,Cij345R,Cij245R,Cij235R,Cij234R,C03
     &   45I,C0245I,C0235I,C0234I,Cij345I,Cij245I,Cij235I,Cij234I,D02
     &   345,D02345R,D02345I,Dij2345R,Dij2345I)
        call tens_red4_new_Re_Com_Div(s12,p3sq,p4sq,p1p3+p2p3,p1p4+p2p4
     &   ,p3p4,C0345R,C0145R,C0135R,C0134R,Cij345R,Cij145R,Cij135R,Ci
     &   j134R,C0345I,C0145I,C0135I,C0134I,Cij345I,Cij145I,Cij135I,Ci
     &   j134I,D01345,D01345R,D01345I,Dij1345R,Dij1345I)
        call tens_red4_new_Re_Com_Div(p1sq,s23,p4sq,p1p2+p1p3,p1p4,p2p4
     &   +p3p4,C0245R,C0145R,C0125R,C0124R,Cij245R,Cij145R,Cij125R,Ci
     &   j124R,C0245I,C0145I,C0125I,C0124I,Cij245I,Cij145I,Cij125I,Ci
     &   j124I,D01245,D01245R,D01245I,Dij1245R,Dij1245I)
        call tens_red4_new_Re_Com_Div(p1sq,p2sq,s34,p1p2,p1p3+p1p4,p2p3
     &   +p2p4,C0235R,C0135R,C0125R,C0123R,Cij235R,Cij135R,Cij125R,Ci
     &   j123R,C0235I,C0135I,C0125I,C0123I,Cij235I,Cij135I,Cij125I,Ci
     &   j123I,D01235,D01235R,D01235I,Dij1235R,Dij1235I)
        call tens_red4_new_Re_Com_Div(p1sq,p2sq,p3sq,p1p2,p1p3,p2p3,C02
     &   34R,C0134R,C0124R,C0123R,Cij234R,Cij134R,Cij124R,Cij123R,C02
     &   34I,C0134I,C0124I,C0123I,Cij234I,Cij134I,Cij124I,Cij123I,D01
     &   234,D01234R,D01234I,Dij1234R,Dij1234I)
c************************************************************************************
c************************************************************************************
       EE0=E0fin(p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15,D02345,D
     &   01345,D01245,D01235,D01234)
       EE0R=Dble(EE0) 
       EE0I=DIMAG(EE0) 
c************************************************************************************
c************************************************************************************
       call tens_red5_new_Re_Com(p1sq,p2sq,p3sq,p4sq,p1p2,p1p3,p1p4,p2
     &   p3,p2p4,p3p4,D02345R,D01345R,D01245R,D01235R,D01234R,Dij2345
     &   R,Dij1345R,Dij1245R,Dij1235R,Dij1234R,D02345I,D01345I,D01245
     &   I,D01235I,D01234I,Dij2345I,Dij1345I,Dij1245I,Dij1235I,Dij123
     &   4I,EijR,EijI)
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call PenBoxFFhl1(F(1))
       call PenBoxFFhl2(F(44))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup3 = dotrc(p1,mup3)
       p2mup3 = dotrc(p2,mup3)
       p3mup3 = dotrc(p3,mup3)
       p4mup3 = dotrc(p4,mup3)
       p5mup3 = dotrc(p5,mup3)
c************** Calling the Fa functions**********************************************************************
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
       psi_p5_M(1)=psi_p5(1)
       psi_p5_M(2)=psi_p5(2)
       psi_p5_P(1)=psi_p5(3)
       psi_p5_P(2)=psi_p5(4)
       psi_p4_M(1)=psi_p4(1)
       psi_p4_M(2)=psi_p4(2)
       psi_p4_P(1)=psi_p4(3)
       psi_p4_P(2)=psi_p4(4)
       barpsi_p1_P(1)=barpsi_p1(1)
       barpsi_p1_P(2)=barpsi_p1(2)
       barpsi_p1_M(1)=barpsi_p1(3)
       barpsi_p1_M(2)=barpsi_p1(4)
       barpsi_p2_P(1)=barpsi_p2(1)
       barpsi_p2_P(2)=barpsi_p2(2)
       barpsi_p2_M(1)=barpsi_p2(3)
       barpsi_p2_M(2)=barpsi_p2(4)
c************************************************************************************
c************************************************************************************
       SMB0(1) = SC1c(barpsi_p1_P,mup3,psi_p5_P,+1)
       SMB0(2) = SC1r(barpsi_p2_M,p1,psi_p4_M,-1)
       SMB0(3) = SC1r(barpsi_p2_M,p5,psi_p4_M,-1)
       SMB0(4) = SC1r(barpsi_p1_P,p2,psi_p5_P,+1)
       SMB0(5) = SC1c(barpsi_p2_M,mup3,psi_p4_M,-1)
       SMB0(6) = SC3crr(barpsi_p2_M,mup3,p5,p1,psi_p4_M,-1)
       SMB0(7) = SC1r(barpsi_p1_P,p4,psi_p5_P,+1)
       SMB0(8) = SC3crr(barpsi_p1_P,mup3,p4,p2,psi_p5_P,+1)
c************************************************************************************
c************************************************************************************
      do i=0,3
        v1(0)=delta(i,0)
        v1(1)=delta(i,1)
        v1(2)=delta(i,2)
        v1(3)=delta(i,3)
       SMB11(i) = SC1c(barpsi_p1_P,v1,psi_p5_P,+1)
       SMB12(i) = SC1c(barpsi_p2_M,v1,psi_p4_M,-1)
      enddo
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       SMB1(1) = SMB0(1)*SMB0(2)
       SMB1(2) = SMB0(1)*SMB0(3)
       SMB1(4) = SMB0(4)*SMB0(5)
       SMB1(5) = SMB0(2)*SMB0(4)
       SMB1(6) = SMB0(3)*SMB0(4)
       SMB1(7) = SMB0(4)*SMB0(6)
       SMB1(8) = SMB0(5)*SMB0(7)
       SMB1(9) = SMB0(2)*SMB0(7)
       SMB1(10) = SMB0(3)*SMB0(7)
       SMB1(11) = SMB0(6)*SMB0(7)
       SMB1(12) = SMB0(2)*SMB0(8)
       SMB1(13) = SMB0(3)*SMB0(8)
c************************************************************************************
       SMB1(3) = dotcc(SMB11,SMB12)
c************************************************************************************
c************************************************************************************
       SMB(1) = SMB1(4)
       SMB(2) = SMB1(8)
       SMB(3) = SMB1(3)
       SMB(4) = SMB1(1)
       SMB(5) = -2*SMB1(1)+2*p1mup3*SMB1(3)
       SMB(6) = SMB1(2)
       SMB(7) = 0
       SMB(8) = SMB1(6)
       SMB(9) = SMB1(5)
       SMB(10) = SMB1(7)
       SMB(11) = s15*SMB1(3)
       SMB(12) = SMB1(10)
       SMB(13) = SMB1(11)
       SMB(14) = (p3sq+s15-s23-s34)*SMB1(3)
       SMB(15) = 4*SMB1(3)
       SMB(16) = 0
       SMB(17) = 0
       SMB(18) = 0
       SMB(19) = 0
       SMB(20) = 0
       SMB(21) = 2*s15*(p2mup3*SMB1(3)-SMB1(4))
       SMB(22) = 0
       SMB(23) = 0
       SMB(24) = SMB1(9)
       SMB(25) = 0
       SMB(26) = 0
       SMB(27) = 0
       SMB(28) = 0
       SMB(29) = 0
       SMB(30) = SMB1(13)
       SMB(31) = 0
       SMB(32) = 0
       SMB(33) = -8*SMB1(1)+8*p1mup3*SMB1(3)
       SMB(34) = 0
       SMB(35) = 2*p2mup3*SMB1(3)-2*SMB1(4)
       SMB(36) = -2*p5mup3*SMB1(5)+2*p1mup3*SMB1(6)+SMB1(7)
       SMB(37) = -2*(p3sq+s15-s23-s34)*(SMB1(1)-p1mup3*SMB1(3))
       SMB(38) = SMB1(12)
       SMB(39) = 0
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
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Amplitude                         
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************


c The Finite virtual piece should be multiplied for (-1)  since 
c The born factor is -1. The virtual (I)
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*(I)=(-1)!!!
       result(1) = F(25)*SMB(2)+p2mup3*F(26)*SMB(3)+p4mup3*F(27)*SMB(3
     &   )+p1mup3*F(28)*SMB(3)+p5mup3*F(29)*SMB(3)+F(30)*SMB(5)+F(31)
     &   *SMB(6)+F(32)*SMB(7)+p2mup3*F(33)*SMB(8)+p5mup3*F(34)*SMB(8)
     &   +p1mup3*F(37)*SMB(11)+p5mup3*F(38)*SMB(11)+p2mup3*F(39)*SMB(
     &   11)+p4mup3*F(40)*SMB(11)+p4mup3*F(43)*SMB(14)+p2mup3*F(44)*S
     &   MB(14)+p5mup3*F(45)*SMB(14)+p1mup3*F(46)*SMB(14)+F(47)*SMB(1
     &   6)+F(49)*SMB(18)+p2mup3*F(52)*SMB(19)+p1mup3*F(53)*SMB(19)+p
     &   2mup3*F(54)*SMB(20)+p5mup3*F(55)*SMB(20)+F(58)*SMB(22)+F(59)
     &   *SMB(23)+p1mup3*F(60)*SMB(24)+8*(p1mup3*F(35)*SMB(8)-p4mup3*
     &   F(36)*SMB(8)+p1mup3*F(41)*SMB(12)-p2mup3*F(42)*SMB(12)-p4mup
     &   3*F(61)*SMB(24)-p5mup3*F(62)*SMB(24)-p2mup3*F(63)*SMB(24))+p
     &   1mup3*F(64)*SMB(25)+p5mup3*F(68)*SMB(26)+F(71)*(2*SMB(27)+SM
     &   B(28))+F(72)*SMB(29)+2*(p1mup3*F(73)*SMB(15)+p2mup3*F(74)*SM
     &   B(15)-p4mup3*F(75)*SMB(15)-p5mup3*F(76)*SMB(15)+p5mup3*F(50)
     &   *SMB(19)-p4mup3*F(51)*SMB(19)-p4mup3*F(56)*SMB(20)-p1mup3*F(
     &   57)*SMB(20)-F(77)*SMB(21)-p5mup3*F(65)*SMB(25)+p4mup3*F(66)*
     &   SMB(25)+p2mup3*F(67)*SMB(25)-p2mup3*F(36)*SMB(26)+p1mup3*F(6
     &   9)*SMB(26)-p4mup3*F(70)*SMB(26)+F(78)*SMB(32))-F(48)*(SMB(17
     &   )+SMB(33))+F(79)*SMB(37)-4*(F(80)*SMB(1)-F(81)*SMB(4)-F(82)*
     &   SMB(10)+F(83)*SMB(13)+F(42)*SMB(30)+F(84)*SMB(35)+F(85)*SMB(
     &   38)+F(86)*SMB(39))
       result(2) = Invs34MU*(F(1)*SMB(1)+p2mup3*F(2)*SMB(3)+p5mup3*F(3
     &   )*SMB(3)+p1mup3*F(4)*SMB(3)+F(5)*SMB(4)+F(6)*SMB(5)+F(7)*SMB
     &   (6)+F(8)*SMB(7)-16*p1mup3*F(9)*SMB(8)+p5mup3*F(10)*SMB(9)+F(
     &   11)*SMB(10)-F(12)*SMB(17)+F(13)*SMB(18)+2*(p1mup3*F(12)*SMB(
     &   15)+p5mup3*F(12)*SMB(15)+p5mup3*F(14)*SMB(19)-p1mup3*F(14)*S
     &   MB(20))+F(15)*SMB(22)-6*(F(12)*SMB(16)+F(12)*SMB(32))-F(12)*
     &   SMB(33)+8*(F(12)*SMB(31)+F(12)*SMB(34))+F(16)*SMB(35)+4*(p2m
     &   up3*F(14)*SMB(11)-F(14)*SMB(21)-F(14)*SMB(36)))
       result(3) = Invs23MU*(F(17)*SMB(2)+F(18)*SMB(5)+F(19)*SMB(7)+F(
     &   20)*SMB(13)+F(22)*SMB(23)+F(23)*SMB(28)-F(21)*(SMB(17)+SMB(3
     &   3))+F(24)*SMB(39))
       result(1) =-result(1)
       result(2) =-result(2)
       result(3) =-result(3)
       If (ngluon.eq.0) then
       resultn=(result(1)+result(2)+result(3))
       elseif (ngluon.eq.1) then
       If (posgluon.eq.3) then
       result(1)=(result(1)+result(2))
       result(2)=(result(3))
       elseIf (posgluon.eq.4) then
       result(1)=(result(1)+result(3))
       result(2)=(result(2))
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     &   . Look to the heading for explanation" 
       endif
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     &   . Look to the heading for explanation" 
       endif
       resultborn = -((2*SMB(2)+SMB(5)+SMB(7))/(s15*s23))-(-2*SMB(1)-2
     &   *(p1mup3+p5mup3)*SMB(3)+2*SMB(4)+SMB(5)+2*SMB(6)+SMB(7))/(s1
     &   5*s34)
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
