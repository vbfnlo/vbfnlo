c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 18/01/2010
c Modified:11/6/2012
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
       subroutine NoAbe4vertBoxLine_QUAD(p1t,p2t,p3t,p4t,p5t,barpsi_p3,psi_p1,mu
     &   p2,mup4,mup5,alpha,musqin,ngluon,posgluon,comp,result_graph,re
     &   sult_color)
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
      Complex*16 SMB(8),Fa(27),F(21)


      Complex*16 barpsi_p3(2),psi_p1(2)
       Real*8 P(3)
       Complex*16  SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr,SC5ccrrc
       EXTERNAL    SC1c,SC1r, SC3crr,SC3ccr,SC3crc,SC3ccc,SC3rrc,SC3rrr,SC5ccrrc
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
      Real*8   p1t(0:3),p2t(0:3),p3t(0:3),p4t(0:3),p5t(0:3)
       Complex*16   mup2(0:3),mup4(0:3),mup5(0:3)
       Complex*16   p1mup2, p1mup4, p1mup5, p2mup2, p2mup4, p2mup5, 
     -          p3mup2, p3mup4, p3mup5, p4mup2, p4mup4, p4mup5, 
     -          p5mup2, p5mup4, p5mup5
       Complex*16   mup2mup4, mup2mup5, mup4mup5
       Real*8 dotrr
       Complex*16 B0fin,C0fin,D0fin,E0fin
       EXTERNAL dotrr,B0fin,C0fin,D0fin,E0fin
        Real*8   p1sq 
       Real*8   p2sq 
       Real*8   p3sq 
       Real*8   p4sq 
       Real*8   p5sq 
       Real*8   s12, s15 
       Real*8   s23 
       Real*8   s34 
       Real*8   s45


       Real*8  B012R,B013R,B014R,B015R 
       Real*8  B023R,B024R,B025R 
       Real*8  B034R,B035R 
       Real*8  B045R 
       Real*8  B012I,B013I,B014I,B015I 
       Real*8  B023I,B024I,B025I 
       Real*8  B034I,B035I 
       Real*8  B045I


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

       Real*8 EE0R 
       Real*8 EE0I 
       Real*8 EijR(46,4) 
       Real*8 EijI(46,4) 
       Real*8 Invs23MU,Invs12MU 
       Logical PrintB,PrintC,PrintD,PrintE 
      Complex*16 dotrc,dotcc,result_graph(3,3),result_color(3)
       Real*8 musq
      EXTERNAL   dotrc,dotcc
       COMMON/NoAbe4vertBoxLineFaFunctions/Fa
       COMMON/NoAbe4vertBoxLineFhlFunctions/F
      Save/NoAbe4vertBoxLineFhlFunctions/
       COMMON/NoAbe4vertBoxLinePFunctions/P
      Save/NoAbe4vertBoxLinePFunctions/
       COMMON/NoAbe4vertBoxLineInvariants/p1sq,p2sq,p3sq,p4sq,p5sq,s12
     &   ,s23,s34,s45,s15
       COMMON/NoAbe4vertBoxLineEVALFunctions/ EE0R,EijR,EE0I,EijI 
       COMMON/NoAbe4vertBoxLineDVALFunctions/ D01234R,Dij1234R,D01234I
     &   ,Dij1234I, D01235R,Dij1235R,D01235I,Dij1235I, D01245R,Dij124
     &   5R,D01245I,Dij1245I, D01345R,Dij1345R,D01345I,Dij1345I, D023
     &   45R,Dij2345R,D02345I,Dij2345I 
       COMMON/NoAbe4vertBoxLineCVALFunctions/ C0123R,Cij123R,C0123I,Ci
     &   j123I, C0124R,Cij124R,C0124I,Cij124I, C0125R,Cij125R,C0125I,
     &   Cij125I, C0134R,Cij134R,C0134I,Cij134I, C0135R,Cij135R,C0135
     &   I,Cij135I, C0145R,Cij145R,C0145I,Cij145I, C0234R,Cij234R,C02
     &   34I,Cij234I, C0235R,Cij235R,C0235I,Cij235I, C0245R,Cij245R,C
     &   0245I,Cij245I, C0345R,Cij345R,C0345I,Cij345I 
       COMMON/NoAbe4vertBoxLineBVALFunctions/ B012R,B012I, B013R,B013I
     &   , B014R,B014I, B015R,B015I, B023R,B023I, B024R,B024I, B025R,
     &   B025I, B034R,B034I, B035R,B035I, B045R,B045I
       Integer alpha,comp,ngluon, posgluon
       real*8 musqin
       integer jtemp
ccccccccccccccccc
       musq=musqIn

c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the scalar products. Not inlcueded the contraction of the
c       moments with the external currents  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
 

      PrintB=.False. 
      PrintC=.False. 
      PrintD=.False. 
      PrintE=.False.
c************************************************************************************
c************************************************************************************

       do jtemp=0,3
         p1(jtemp)=p1t(jtemp)
         p2(jtemp)=p2t(jtemp)
         p3(jtemp)=p3t(jtemp)
         p4(jtemp)=p4t(jtemp)
         p5(jtemp)=p5t(jtemp)
      enddo

        if (comp.gt.0) then    
c    Calling C_ij,D_ij Functions    
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************

      call PrintIntPen4gMassLess_QUAD(p1,p2,p3,p4,p5,musq)
        Invs23MU=1d0/s23
        Invs12MU=1d0/s12

c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       call NoAbe4vertBoxLineFFhl1(F(1))
c************************************************************************************
c************************************************************************************
c************************************************************************************
       else
        Invs23MU=1d0/s23
        Invs12MU=1d0/s12
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup2 = dotrc(p1,mup2)
       p1mup4 = dotrc(p1,mup4)
       p1mup5 = dotrc(p1,mup5)
       p2mup2 = dotrc(p2,mup2)
       p2mup4 = dotrc(p2,mup4)
       p2mup5 = dotrc(p2,mup5)
       p3mup2 = dotrc(p3,mup2)
       p3mup4 = dotrc(p3,mup4)
       p3mup5 = dotrc(p3,mup5)
       p4mup2 = dotrc(p4,mup2)
       p4mup4 = dotrc(p4,mup4)
       p4mup5 = dotrc(p4,mup5)
       p5mup2 = dotrc(p5,mup2)
       p5mup4 = dotrc(p5,mup4)
       p5mup5 = dotrc(p5,mup5)
       mup2mup4 = dotcc(mup2,mup4)
       mup2mup5 = dotcc(mup2,mup5)
       mup4mup5 = dotcc(mup4,mup5)
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call NoAbe4vertBoxLineFa1(p1mup2,p1mup4,p1mup5,p2mup2,p2mup4,p2
     &   mup5,p3mup2,p3mup4,p3mup5,p4mup2,p4mup4,p4mup5,p5mup2,p5mup4
     &   ,p5mup5,mup2mup4,mup2mup5,mup4mup5,Fa(1))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c       Definition of the Matrix Element  
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       SMB(1) = SC1c(barpsi_p3,mup2,psi_p1,alpha)
       SMB(2) = SC1r(barpsi_p3,p2,psi_p1,alpha)
       SMB(3) = SC1c(barpsi_p3,mup4,psi_p1,alpha)
       SMB(4) = SC1c(barpsi_p3,mup5,psi_p1,alpha)
       SMB(5) = SC3ccc(barpsi_p3,mup5,mup4,mup2,psi_p1,alpha)
       SMB(6) = SC3crc(barpsi_p3,mup4,p2,mup2,psi_p1,alpha)
       SMB(7) = SC3crc(barpsi_p3,mup5,p2,mup2,psi_p1,alpha)
       SMB(8) = SC3ccr(barpsi_p3,mup5,mup4,p2,psi_p1,alpha)
c       Print*," SMB(1) ", SMB(1)
c       Print*," SMB(2) ", SMB(2)
c       Print*," SMB(3) ", SMB(3)
c       Print*," SMB(4) ", SMB(4)
c       Print*," SMB(5) ", SMB(5)
c       Print*," SMB(6) ", SMB(6)
c       Print*," SMB(7) ", SMB(7)
c       Print*," SMB(8) ", SMB(8)
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


c The Finite virtual piece should be multiplied for (-1)  since 
c The born factor is -1. The virtual (I)
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*(I)=(-1)!!!
       result_graph(1,1) = 4*p3mup5*Fa(1)*SMB(3)+Fa(5)*SMB(4)+F(18)*SM
     &   B(5)+2*(p3mup5*F(16)*(-2*(p2mup4*SMB(1)-mup2mup4*SMB(2))+SMB
     &   (6))+p1mup4*F(19)*SMB(7)-Fa(6)*SMB(8))
       result_graph(1,2) = mup4mup5*(F(3)*SMB(1)-4*Fa(7)*SMB(2))
       result_graph(1,3) = Fa(8)*SMB(1)+Fa(12)*SMB(3)-4*(Fa(13)*SMB(2)
     &   -p3mup4*Fa(1)*SMB(4))+F(21)*SMB(5)+2*(p1mup5*F(19)*SMB(6)+p3
     &   mup4*F(16)*SMB(7)+Fa(6)*SMB(8))
       result_graph(2,1) = (-4*(Fa(15)*SMB(1)-s23*Fa(16)*SMB(2))+F(1)*
     &   (s23*SMB(5)+2*Fa(14)*SMB(6))-2*(F(1)*Fa(18)*SMB(3)+s23*Fa(19
     &   )*SMB(4)-p1mup4*s23*F(7)*SMB(7)))/s23
       result_graph(2,2) = 2*mup4mup5*F(1)*SMB(1)
       result_graph(2,3) = (Fa(22)*SMB(1)+s23*(4*Fa(23)*SMB(2)-F(1)*SM
     &   B(5)-2*(Fa(24)*SMB(3)-p1mup5*F(7)*SMB(6)))-2*(F(1)*Fa(25)*SM
     &   B(4)-F(1)*Fa(20)*SMB(7)))/s23
       result_graph(3,1) = (2*p3mup5*F(8)*(2*p1mup2*SMB(3)+SMB(6))+s12
     &   *(F(10)*SMB(5)+2*F(9)*Fa(26)*(2*p1mup2*SMB(4)+SMB(7))))/s12
       result_graph(3,2) = -2*mup4mup5*F(2)*SMB(1)
       result_graph(3,3) = -((s12*F(2)*(2*mup4mup5*SMB(1)-SMB(5))-2*(s
     &   12*F(9)*Fa(27)*(2*p1mup2*SMB(3)+SMB(6))+p3mup4*F(8)*(2*p1mup
     &   2*SMB(4)+SMB(7))))/s12)
       result_graph(1,1) =-result_graph(1,1)
       result_graph(1,2) =-result_graph(1,2)
       result_graph(1,3) =-result_graph(1,3)
       result_graph(2,1) =-result_graph(2,1)
       result_graph(2,2) =-result_graph(2,2)
       result_graph(2,3) =-result_graph(2,3)
       result_graph(3,1) =-result_graph(3,1)
       result_graph(3,2) =-result_graph(3,2)
       result_graph(3,3) =-result_graph(3,3)
       If (ngluon.eq.2) then
       result_color(1)=(result_graph(1,1)+result_graph(2,1)+result_gra
     &   ph(3,1))
       result_color(2)=(result_graph(1,2)+result_graph(2,2)+result_gra
     &   ph(3,2))
       result_color(3)=(result_graph(1,3)+result_graph(2,3)+result_gra
     &   ph(3,3))
       elseif (ngluon.eq.3) then
       result_color(1)=(result_graph(1,1))
       result_color(2)=(result_graph(1,1))
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     &   . Look to the heading for explanation" 
       endif
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End
