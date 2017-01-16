       subroutine boxlineABEmm_QUAD(p1t,p2t,p3t,p4t,barpsi_p4,psi_p1,mup2,mup3,
     -   alpha,musqIn,ngluon,posgluon,gaugetest,comp,resultgauge,result
     -   ,resultgaugeb,resultb)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:25/7/2008
c ************************************************************************************
c determine the  finite virtual corrections along the quark line i.e 
c the sum of all virtual corrections to the Born diagram 
c psi(p1)    ---->------->------->---- bar_psi(p4)
c                   $       $                    
c                   $       $                     
c                   $       $                     
c                   V       V                     
c                   $       $                     
c                   $       $                     
c                  p2       p3                   
c                 mu_p2    mu_p3              
c Note: The vertices are just Gamma^(mu_p2)..terms. So the correct
c factor should be added by hand for each boson. i.e, Without adding
c anything it represents pp->3 photons up to the electromagnetic charge 
c (for uu->3 photons is missing (2/3)^3) from the coupling. For W and Z,
c we have to account for the Diracgamma_5 that are missing in the vertex
c Note: To make it shorter in the promgram: mu_p2,...->mup2,... 
c Notation of External momenta: p1+p2+p3+p4=0 
c mu_p2,mu_p3 should be think as external current 
c alpha is the helicity of the initial spinor 
c musq is the renormalization scale energy  
c gaugetest,integer value.Different gauge test: 
c gaugetest=0 should give zero 
c gaugetest=1 give the result after replacing mu_p2 with p2 
c gaugetest=2 give the result after replacing mu_p3 with p3 
c gaugetest=3 calculates the two different gaugue tests 
c gaugetest>5 no calculation of gauge test in case you
c have already done it before(it safes times)
c comp: integer value.The first time called with p1...p4, comp=1
c ATTENTION: ONLY!!!If you have to call the subroutine consecutively with the same arguments
c(p1,p2,p3,p4). Then, comp=-1 (it safes 400 lines of code) 
c This applies when you have for examples the same diagram for an off-shell photon
c and a Z boson. The differences are in the coupling and  the part that depends on the
c polarization vector that are calculated at the end of this program.
c resultgauge and resultgaugeb are arrays of dimension three. 
c In case you use gaugetest=(0,2). The result is given in the first argument:
c resultgauge(1) and resultgaugeb(1).The argument 2 is set to zero 
c In case you use gaugetest=3. resultgauge(1) is the result of gaugetest=1,
c resultgauge(2) of gaugetest=2 
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
      Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
      Real*8   p1t(0:3),p2t(0:3),p3t(0:3),p4t(0:3)
       Real*8   p1sq, p1p2, p1p3, p2sq, p2p3, p2p4, p3sq, 
     -          p3p4, p4sq
       Real*8   s12, s23
       Complex*16   p1mup2, p1mup3, p2mup2, p2mup3, p3mup2, p3mup3, 
     -          p4mup2, p4mup3
       Complex*16 mup2mup3!,mup2mup4,mup3mup4
       Real*8 dotrr
       Complex*16 B0fin,C0fin,D0fin
       EXTERNAL dotrr,B0fin,C0fin,D0fin
c       Complex*16 B012,B013,B014,B023,B024,B034
      Real*8 B012R,B013R,B014R,B023R,B024R,B034R
      Real*8 B012I,B013I,B014I,B023I,B024I,B034I

       Real*8 C0123R,C0124R,C0134R,C0234R
       Real*8 C0123I,C0124I,C0134I,C0234I
       Real*8 Cij123R(4,2),Cij124R(4,2),Cij134R(4,2),Cij234R(4,2)
       Real*8 Cij123I(4,2),Cij124I(4,2),Cij134I(4,2),Cij234I(4,2)

        Real*8 D01234R
        Real*8 D01234I
        Real*8 Dij1234R(13,3)
        Real*8 Dij1234I(13,3)
       Complex*16 SMB(4),SMG(4) ,Fa(5),F(27),K(1)
       Real*8 FI(27),FR(27),PaR(9),PaI(9),KI(1),KR(1)
       Complex*16 barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3)
       Complex*16 SC1c,SC3rcc,SC1r,dotrc,dotcc,result(5),resultg(5)
       Complex*16 resultgauge(2),resultgaugeb(2),resultb 
       Real*8 musq,Is12,P(18) 
       EXTERNAL   dotrc,dotcc,SC1c,SC1r,SC3rcc
       Integer alpha,comp,gaugetest,ngluon,posgluon
       common/invarianetsmm/p1sq,p2sq,p3sq,p4sq,s12,s23,Is12
       common/FfunctionsABEmm/F
       common/PfunctionsABEmm/P
       common/KfunctionsABEmm/K
       SAVE/FfunctionsABEmm/
       SAVE/PfunctionsABEmm/
       SAVE/KfunctionsABEmm/
       Real*8 musqIn
       Common/musqInv/musq
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      End QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       COMMON/DVALFunctionsD/ D01234R,DIJ1234R,D01234I,DIJ1234I

       COMMON/CVALFunctionsD/ C0123R,CIJ123R,C0123I,CIJ123I, 
     &    C0124R,CIJ124R,C0124I,CIJ124I,
     &    C0134R,CIJ134R,C0134I,CIJ134I,
     &    C0234R,CIJ234R,C0234I,CIJ234I
  
       COMMON/BVALFunctionsD/ B012R,B012I, B013R,B013I, B014R,B014I, 
     &      B023R,B023I,  B024R,B024I,
     &    B034R,B034I
       real*8 p1sqT,p2sqT,p3sqT,p4sqT,s12T,s23T
       COMMON/InvariantsD/p1sqT,p2sqT,p3sqT,p4sqT,s12T,s23T

       common/invariantsrest/p1p2
       integer jtemp
ccccccccccccccccc
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
       If(comp.gt.0) then

       do jtemp=0,3
         p1(jtemp)=p1t(jtemp)
         p2(jtemp)=p2t(jtemp)
         p3(jtemp)=p3t(jtemp)
         p4(jtemp)=p4t(jtemp)
      enddo
         call  PrintIntBoxMassless_QUAD(p1,p2,p3,p4,musq)
         p1p2=dotrr(p1,p2)
         p1p3=dotrr(p1,p3)
         p2p3=dotrr(p2,p3)
         p2p4=dotrr(p2,p4)
         p3p4=dotrr(p3,p4)
c         b013=DCMPLX(B013R,B013I)
          p1sq=p1sqT
          p2sq=p2sqT
          p3sq=p3sqT
          p4sq=p4sqT
          s12=s12T
          s23=s23T         
         Is12=1d0/s12


c************************************************************************************
c************************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c************************************************************************************
c************************************************************************************
c************************************************************************************
       P(1) = p2sq-s12
       PaR(1) = C0123R+Cij123R(1,1)
       PaI(1) = C0123I+Cij123I(1,1)
       PaR(2) = B023R+s12*Cij123R(2,1)-2*Cij123R(4,2)-P(1)*PaR(1)
       PaI(2) = B023I+s12*Cij123I(2,1)-2*Cij123I(4,2)-P(1)*PaI(1)
       FR(1) = -4*Is12*PaR(2)
       FI(1) = -4*Is12*PaI(2)
       F(1)=DCMPLX(FR(1),FI(1))
       PaR(3) = -B023R+s12*Cij123R(3,2)+2*Cij123R(4,2)+p2sq*PaR(1)
       PaI(3) = -B023I+s12*Cij123I(3,2)+2*Cij123I(4,2)+p2sq*PaI(1)
       FR(2) = 4*Is12*PaR(3)
       FI(2) = 4*Is12*PaI(3)
       F(2)=DCMPLX(FR(2),FI(2))
       PaR(4) = Cij123R(2,1)+Cij123R(2,2)
       PaI(4) = Cij123I(2,1)+Cij123I(2,2)
       FR(3) = 4*PaR(4)
       FI(3) = 4*PaI(4)
       F(3)=DCMPLX(FR(3),FI(3))
       FR(4) = 4*Is12*PaR(2)
       FI(4) = 4*Is12*PaI(2)
       F(4)=DCMPLX(FR(4),FI(4))
       FR(5) = -2*Is12*PaR(2)
       FI(5) = -2*Is12*PaI(2)
       F(5)=DCMPLX(FR(5),FI(5))
       PaR(5) = Cij134R(1,1)+Cij134R(1,2)-Cij134R(2,1)+Cij134R(2,2)-2*
     -   Cij134R(3,2)
       PaI(5) = Cij134I(1,1)+Cij134I(1,2)-Cij134I(2,1)+Cij134I(2,2)-2*
     -   Cij134I(3,2)
       FR(6) = 4*PaR(5)
       FI(6) = 4*PaI(5)
       F(6)=DCMPLX(FR(6),FI(6))
       PaR(6) = C0134R+Cij134R(1,1)-Cij134R(2,2)+Cij134R(3,2)
       PaI(6) = C0134I+Cij134I(1,1)-Cij134I(2,2)+Cij134I(3,2)
       FR(7) = -4*PaR(6)
       FI(7) = -4*PaI(6)
       F(7)=DCMPLX(FR(7),FI(7))
       P(2) = p3sq-s12
       PaR(7) = Cij134R(1,2)+Cij134R(2,2)-2*Cij134R(3,2)
       PaI(7) = Cij134I(1,2)+Cij134I(2,2)-2*Cij134I(3,2)
       PaR(8) = -B034R+2*Cij134R(4,2)+(C0134R+Cij134R(1,1))*P(2)+s12*P
     -   aR(7)
       PaI(8) = -B034I+2*Cij134I(4,2)+(C0134I+Cij134I(1,1))*P(2)+s12*P
     -   aI(7)
       FR(8) = 4*Is12*PaR(8)
       FI(8) = 4*Is12*PaI(8)
       F(8)=DCMPLX(FR(8),FI(8))
       P(3) = p3sq-2*s12
       PaR(9) = -B034R+s12*Cij134R(2,1)+2*Cij134R(4,2)+C0134R*P(2)+Cij
     -   134R(1,1)*P(3)
       PaI(9) = -B034I+s12*Cij134I(2,1)+2*Cij134I(4,2)+C0134I*P(2)+Cij
     -   134I(1,1)*P(3)
       FR(9) = 4*Is12*PaR(9)
       FI(9) = 4*Is12*PaI(9)
       F(9)=DCMPLX(FR(9),FI(9))
       FR(10) = -4*Is12*PaR(9)
       FI(10) = -4*Is12*PaI(9)
       F(10)=DCMPLX(FR(10),FI(10))
       FR(11) = 2*Is12*PaR(9)
       FI(11) = 2*Is12*PaI(9)
       F(11)=DCMPLX(FR(11),FI(11))
       P(4) = p2sq+p3sq-s12-s23
       FR(12) = 4*(p2sq*(Dij1234R(2,1)-Dij1234R(3,1)+Dij1234R(4,2)-Dij
     -   1234R(6,2))+2*(Dij1234R(7,2)+Dij1234R(11,3)-Dij1234R(13,3))+
     -   (Dij1234R(3,2)-Dij1234R(5,2))*P(4))
       FI(12) = 4*(p2sq*(Dij1234I(2,1)-Dij1234I(3,1)+Dij1234I(4,2)-Dij
     -   1234I(6,2))+2*(Dij1234I(7,2)+Dij1234I(11,3)-Dij1234I(13,3))+
     -   (Dij1234I(3,2)-Dij1234I(5,2))*P(4))
       F(12)=DCMPLX(FR(12),FI(12))
       P(5) = -p2sq+s12
       FR(13) = 4*(-Cij234R(1,1)+Cij234R(2,1)-s12*(Dij1234R(2,2)+2*Dij
     -   1234R(3,2)+Dij1234R(4,2)-Dij1234R(5,2)-3*Dij1234R(6,2))+p2sq
     -   *(Dij1234R(3,2)+Dij1234R(4,2)-Dij1234R(5,2)-Dij1234R(6,2))+s
     -   23*(D01234R+Dij1234R(1,1)-Dij1234R(3,2)+Dij1234R(6,2))-2*(Di
     -   j1234R(7,2)-Dij1234R(12,3)+Dij1234R(13,3))+Dij1234R(2,1)*P(1
     -   )+Dij1234R(3,1)*P(5))
       FI(13) = 4*(-Cij234I(1,1)+Cij234I(2,1)-s12*(Dij1234I(2,2)+2*Dij
     -   1234I(3,2)+Dij1234I(4,2)-Dij1234I(5,2)-3*Dij1234I(6,2))+p2sq
     -   *(Dij1234I(3,2)+Dij1234I(4,2)-Dij1234I(5,2)-Dij1234I(6,2))+s
     -   23*(D01234I+Dij1234I(1,1)-Dij1234I(3,2)+Dij1234I(6,2))-2*(Di
     -   j1234I(7,2)-Dij1234I(12,3)+Dij1234I(13,3))+Dij1234I(2,1)*P(1
     -   )+Dij1234I(3,1)*P(5))
       F(13)=DCMPLX(FR(13),FI(13))
       P(6) = s12+s23
       P(7) = p2sq-2*s12-s23
       FR(14) = 4*(C0234R+Cij234R(2,1)+s12*(Dij1234R(2,1)+Dij1234R(6,2
     -   ))-2*(Dij1234R(7,2)+Dij1234R(13,3))-Dij1234R(5,2)*P(1)-Dij12
     -   34R(3,1)*P(6)+Dij1234R(3,2)*P(7))
       FI(14) = 4*(C0234I+Cij234I(2,1)+s12*(Dij1234I(2,1)+Dij1234I(6,2
     -   ))-2*(Dij1234I(7,2)+Dij1234I(13,3))-Dij1234I(5,2)*P(1)-Dij12
     -   34I(3,1)*P(6)+Dij1234I(3,2)*P(7))
       F(14)=DCMPLX(FR(14),FI(14))
       P(8) = -p2sq+s12+s23
       P(9) = p2sq-s12-s23
       P(10) = p2sq+2*p3sq-s12-s23
       FR(15) = 4*(-C0234R+p3sq*Dij1234R(3,2)+4*Dij1234R(7,2)+2*(Dij12
     -   34R(12,3)-Dij1234R(13,3))+(Dij1234R(2,1)-Dij1234R(3,1))*P(2)
     -   +Dij1234R(2,2)*P(4)+Dij1234R(4,2)*P(8)+Dij1234R(5,2)*P(9)-Di
     -   j1234R(6,2)*P(10))
       FI(15) = 4*(-C0234I+p3sq*Dij1234I(3,2)+4*Dij1234I(7,2)+2*(Dij12
     -   34I(12,3)-Dij1234I(13,3))+(Dij1234I(2,1)-Dij1234I(3,1))*P(2)
     -   +Dij1234I(2,2)*P(4)+Dij1234I(4,2)*P(8)+Dij1234I(5,2)*P(9)-Di
     -   j1234I(6,2)*P(10))
       F(15)=DCMPLX(FR(15),FI(15))
       P(11) = p2sq+p3sq
       FR(16) = -4*(C0234R+s12*Dij1234R(2,1)-s12*Dij1234R(3,1)-p3sq*Di
     -   j1234R(3,2)-p2sq*Dij1234R(5,2)-2*(Dij1234R(7,2)-Dij1234R(13,
     -   3))+(Dij1234R(5,2)-Dij1234R(6,2))*P(6)+Dij1234R(6,2)*P(11))
       FI(16) = -4*(C0234I+s12*Dij1234I(2,1)-s12*Dij1234I(3,1)-p3sq*Di
     -   j1234I(3,2)-p2sq*Dij1234I(5,2)-2*(Dij1234I(7,2)-Dij1234I(13,
     -   3))+(Dij1234I(5,2)-Dij1234I(6,2))*P(6)+Dij1234I(6,2)*P(11))
       F(16)=DCMPLX(FR(16),FI(16))
       P(12) = p3sq-2*s12-s23
       P(13) = s12-s23
       P(14) = p3sq-s23
       P(15) = 2*s12-P(14)
       FR(17) = 4*(Cij234R(2,1)+s23*(D01234R+2*Dij1234R(1,1)+Dij1234R(
     -   1,2))-s12*Dij1234R(3,2)+2*(Dij1234R(11,3)-Dij1234R(13,3))-Di
     -   j1234R(3,1)*P(3)+(Dij1234R(2,1)+Dij1234R(4,2))*P(12)+Dij1234
     -   R(5,2)*P(13)+Dij1234R(6,2)*P(15))
       FI(17) = 4*(Cij234I(2,1)+s23*(D01234I+2*Dij1234I(1,1)+Dij1234I(
     -   1,2))-s12*Dij1234I(3,2)+2*(Dij1234I(11,3)-Dij1234I(13,3))-Di
     -   j1234I(3,1)*P(3)+(Dij1234I(2,1)+Dij1234I(4,2))*P(12)+Dij1234
     -   I(5,2)*P(13)+Dij1234I(6,2)*P(15))
       F(17)=DCMPLX(FR(17),FI(17))
       FR(18) = 8*(Dij1234R(2,1)-Dij1234R(3,1)-Dij1234R(3,3)+Dij1234R(
     -   4,3)-Dij1234R(5,3)+Dij1234R(9,3)+2*(Dij1234R(3,2)+Dij1234R(4
     -   ,2)-Dij1234R(5,2)-Dij1234R(6,2)+Dij1234R(7,3)-Dij1234R(10,3)
     -   ))
       FI(18) = 8*(Dij1234I(2,1)-Dij1234I(3,1)-Dij1234I(3,3)+Dij1234I(
     -   4,3)-Dij1234I(5,3)+Dij1234I(9,3)+2*(Dij1234I(3,2)+Dij1234I(4
     -   ,2)-Dij1234I(5,2)-Dij1234I(6,2)+Dij1234I(7,3)-Dij1234I(10,3)
     -   ))
       F(18)=DCMPLX(FR(18),FI(18))
       FR(19) = 8*(Dij1234R(2,2)+Dij1234R(3,2)-Dij1234R(3,3)+Dij1234R(
     -   6,3)+Dij1234R(7,3)-Dij1234R(8,3)-2*(Dij1234R(6,2)-Dij1234R(9
     -   ,3)+Dij1234R(10,3)))
       FI(19) = 8*(Dij1234I(2,2)+Dij1234I(3,2)-Dij1234I(3,3)+Dij1234I(
     -   6,3)+Dij1234I(7,3)-Dij1234I(8,3)-2*(Dij1234I(6,2)-Dij1234I(9
     -   ,3)+Dij1234I(10,3)))
       F(19)=DCMPLX(FR(19),FI(19))
       FR(20) = -8*(Dij1234R(3,3)-Dij1234R(5,2)+Dij1234R(6,2)-Dij1234R
     -   (7,3)-Dij1234R(9,3)+Dij1234R(10,3))
       FI(20) = -8*(Dij1234I(3,3)-Dij1234I(5,2)+Dij1234I(6,2)-Dij1234I
     -   (7,3)-Dij1234I(9,3)+Dij1234I(10,3))
       F(20)=DCMPLX(FR(20),FI(20))
       FR(21) = 8*(Dij1234R(2,1)+Dij1234R(2,2)-Dij1234R(3,1)-Dij1234R(
     -   3,3)+Dij1234R(4,2)-Dij1234R(5,2)-3*Dij1234R(6,2)+Dij1234R(6,
     -   3)+Dij1234R(7,3)-Dij1234R(8,3)+2*(Dij1234R(3,2)+Dij1234R(9,3
     -   )-Dij1234R(10,3)))
       FI(21) = 8*(Dij1234I(2,1)+Dij1234I(2,2)-Dij1234I(3,1)-Dij1234I(
     -   3,3)+Dij1234I(4,2)-Dij1234I(5,2)-3*Dij1234I(6,2)+Dij1234I(6,
     -   3)+Dij1234I(7,3)-Dij1234I(8,3)+2*(Dij1234I(3,2)+Dij1234I(9,3
     -   )-Dij1234I(10,3)))
       F(21)=DCMPLX(FR(21),FI(21))
       FR(22) = 8*(Dij1234R(2,2)+Dij1234R(2,3)+Dij1234R(3,2)-Dij1234R(
     -   3,3)-2*Dij1234R(6,2)-3*(Dij1234R(8,3)-Dij1234R(9,3)))
       FI(22) = 8*(Dij1234I(2,2)+Dij1234I(2,3)+Dij1234I(3,2)-Dij1234I(
     -   3,3)-2*Dij1234I(6,2)-3*(Dij1234I(8,3)-Dij1234I(9,3)))
       F(22)=DCMPLX(FR(22),FI(22))
       FR(23) = -8*(Dij1234R(2,1)-Dij1234R(3,1)-Dij1234R(3,2)+Dij1234R
     -   (3,3)+Dij1234R(4,2)-Dij1234R(7,3)-Dij1234R(9,3)+Dij1234R(10,
     -   3))
       FI(23) = -8*(Dij1234I(2,1)-Dij1234I(3,1)-Dij1234I(3,2)+Dij1234I
     -   (3,3)+Dij1234I(4,2)-Dij1234I(7,3)-Dij1234I(9,3)+Dij1234I(10,
     -   3))
       F(23)=DCMPLX(FR(23),FI(23))
       FR(24) = -8*(Dij1234R(2,1)+Dij1234R(2,2)-Dij1234R(3,1)+Dij1234R
     -   (3,3)-Dij1234R(6,2)+Dij1234R(8,3)-2*Dij1234R(9,3))
       FI(24) = -8*(Dij1234I(2,1)+Dij1234I(2,2)-Dij1234I(3,1)+Dij1234I
     -   (3,3)-Dij1234I(6,2)+Dij1234I(8,3)-2*Dij1234I(9,3))
       F(24)=DCMPLX(FR(24),FI(24))
       FR(25) = -8*(Dij1234R(3,2)+Dij1234R(3,3)-Dij1234R(6,2)-Dij1234R
     -   (9,3))
       FI(25) = -8*(Dij1234I(3,2)+Dij1234I(3,3)-Dij1234I(6,2)-Dij1234I
     -   (9,3))
       F(25)=DCMPLX(FR(25),FI(25))
       FR(26) = -4*(s23*(D01234R+Dij1234R(1,1))-2*(Dij1234R(7,2)+Dij12
     -   34R(12,3)-Dij1234R(13,3)))
       FI(26) = -4*(s23*(D01234I+Dij1234I(1,1))-2*(Dij1234I(7,2)+Dij12
     -   34I(12,3)-Dij1234I(13,3)))
       F(26)=DCMPLX(FR(26),FI(26))
       P(16) = p2sq+s12
       P(17) = p3sq+s12
       P(18) = p2sq+p3sq+2*s12
       FR(27) = 2*(C0234R-Cij234R(1,1)+Cij234R(2,1)+s23*(D01234R+Dij12
     -   34R(1,1))-6*Dij1234R(7,2)+(Dij1234R(4,2)-Dij1234R(5,2))*P(1)
     -   -Dij1234R(2,2)*P(16)-Dij1234R(3,2)*P(17)+Dij1234R(6,2)*P(18)
     -   )
       FI(27) = 2*(C0234I-Cij234I(1,1)+Cij234I(2,1)+s23*(D01234I+Dij12
     -   34I(1,1))-6*Dij1234I(7,2)+(Dij1234I(4,2)-Dij1234I(5,2))*P(1)
     -   -Dij1234I(2,2)*P(16)-Dij1234I(3,2)*P(17)+Dij1234I(6,2)*P(18)
     -   )
       F(27)=DCMPLX(FR(27),FI(27))
       KI(1) = -Dij1234I(3,3)-Dij1234I(8,3)+2*Dij1234I(9,3)
       KR(1) = -Dij1234R(3,3)-Dij1234R(8,3)+2*Dij1234R(9,3)
       K(1)=DCMPLX(KR(1),KI(1))
       else
         p1p2=dotrr(p1,p2)
         p1p3=dotrr(p1,p3)
         p2p3=dotrr(p2,p3)
         p2p4=dotrr(p2,p4)
         p3p4=dotrr(p3,p4)
         Is12=1d0/s12
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       endif  
c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c************************************************************************************
c************************************************************************************
c************************************************************************************
c Computation of Fa fucntions. Depends on the external currents, through
c the contraction of the moments with the currents
c************************************************************************************
c************************************************************************************
c************************************************************************************
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrc(p1,mup3)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrc(p2,mup3)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrc(p3,mup3)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrc(p4,mup3)
       mup2mup3=dotcc(mup2,mup3)
c      Print*," p1mup2 ",  p1mup2 
c      Print*," p1mup3 ",  p1mup3 
c      Print*," p2mup2 ",  p2mup2 
c      Print*," p2mup3 ",  p2mup3 
c      Print*," p3mup2 ",  p3mup2 
c      Print*," p3mup3 ",  p3mup3 
c      Print*," p4mup2 ",  p4mup2 
c      Print*," p4mup3 ",  p4mup3 
c       Print*," mup2mup3", mup2mup3  
c************** Calling the Fa functions**********************************************************************
c************************************************************************************
c************************************************************************************
       call FaFunctionABEmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3m
     -   up3, p4mup2, p4mup3, mup2mup3,Fa)
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
       SMB(3) = SC1r(barpsi_p4,p2,psi_p1,alpha)
       SMB(4) = SC3rcc(barpsi_p4,p2,mup2,mup3,psi_p1,alpha)
c       Print*," SMB(1) ", SMB(1)
c       Print*," SMB(2) ", SMB(2)
c       Print*," SMB(3) ", SMB(3)
c       Print*," SMB(4) ", SMB(4)
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


       result(2) = DCMPLX(B013R,B013I)*Is12*(2*(p2mup3*SMB(1)+p1mup2*SMB(2)-mup2mup3*
     -   SMB(3))+SMB(4))
       result(3) = p2mup3*F(1)*SMB(1)+Fa(1)*SMB(2)+mup2mup3*F(4)*SMB(3
     -   )+F(5)*SMB(4)
       result(4) = Fa(2)*SMB(1)+p1mup2*F(9)*SMB(2)+mup2mup3*F(10)*SMB(
     -   3)+F(11)*SMB(4)
       result(5) = Fa(3)*SMB(1)+Fa(4)*SMB(2)+Fa(5)*SMB(3)+F(27)*SMB(4)
c************************************************************************************
c** Selecting different colour factor****************
c************************************************************************************
       If (ngluon.eq.0) then
       result(1)=-(result(2)+result(3)+result(4)+result(5))
       elseif (ngluon.eq.1) then
       If (posgluon.eq.2) then
       result(1)=-(result(2)+result(4))
       result(2)=-(result(3)+result(5))
       elseIf (posgluon.eq.3) then
       result(1)=-(result(2)+result(3))
       result(2)=-(result(4)+result(5))
       else
       Write(*,*) "Error: The position of the gluon is badly indicated
     -   . Look to the heading for explanation" 
       endif
       elseif (ngluon.eq.2) then
c result(1): coefficient of t^b t^a color structure (i.e. the born color structure)
c result(2): coefficient of delta^ab color structure
       result(1)=+(result(2)+result(3)+result(4)-8d0*result(5))/6d0
       result(2)=-result(2)/4d0
       endif
       resultb = Is12*(2*(p2mup3*SMB(1)+p1mup2*SMB(2)-mup2mup3*SMB(3))
     -   +SMB(4))
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
c      mup2->p2+p3,mup3->p3. SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2 = dotrr(p1,p2)+dotrr(p1,p3)
       p1mup3 = dotrr(p1,p3)
       p2mup2 = dotrr(p2,p2)+dotrr(p2,p3)
       p2mup3 = dotrr(p2,p3)
       p3mup2 = dotrr(p3,p2)+dotrr(p3,p3)
       p3mup3 = dotrr(p3,p3)
       p4mup2 = dotrr(p4,p2)+dotrr(p4,p3)
       p4mup3 = dotrr(p4,p3)
       mup2mup3=dotrr(p2,p3)+dotrr(p3,p3)
c       Print*," p1mup2 ",  p1mup2 
c       Print*," p1mup3 ",  p1mup3 
c       Print*," p2mup2 ",  p2mup2 
c       Print*," p2mup3 ",  p2mup3 
c       Print*," p3mup2 ",  p3mup2 
c       Print*," p3mup3 ",  p3mup3 
c       Print*," p4mup2 ",  p4mup2 
c       Print*," p4mup3 ",  p4mup3 
c       Print*," mup2mup3", mup2mup3  
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FaFunctionABEmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3m
     -   up3, p4mup2, p4mup3,  mup2mup3,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2,4) and mu_p3=p2+p3+p4
c************************************************************************************
c************************************************************************************
       SMG(1) = 0
       SMG(2) = -SMB(3)
       SMG(3) = SMB(3)
       SMG(4) = (-p2sq+p3sq)*SMB(3)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c************************************************************************************
c************************************************************************************
c    This should be zero:
c************************************************************************************
c************************************************************************************
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = DCMPLX(B013R,B013I)*Is12*(2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3
     -   *SMG(3))+SMG(4))
       resultg(3) = p2mup3*F(1)*SMG(1)+Fa(1)*SMG(2)+mup2mup3*F(4)*SMG(
     -   3)+F(5)*SMG(4)
       resultg(4) = Fa(2)*SMG(1)+p1mup2*F(9)*SMG(2)+mup2mup3*F(10)*SMG
     -   (3)+F(11)*SMG(4)
       resultg(5) = Fa(3)*SMG(1)+Fa(4)*SMG(2)+Fa(5)*SMG(3)+F(27)*SMG(4
     -   )
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(1)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5))/S
     -   qrt(abs(p2sq)*abs(s23))
       resultgaugeb(1) = Is12*(2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3
     -   *SMG(3))+SMG(4))
       resultgaugeb(1) =resultgaugeb(1)/Sqrt(abs(p2sq)*abs(s23))
           resultgauge(2)=(0d0,0d0)      
           resultgaugeb(2)=(0d0,0d0)      
           return      
           endif      
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
       mup2mup3=dotrc(p2,mup3)
c       Print*," p1mup2 ",  p1mup2 
c       Print*," p1mup3 ",  p1mup3 
c       Print*," p2mup2 ",  p2mup2 
c       Print*," p2mup3 ",  p2mup3 
c       Print*," p3mup2 ",  p3mup2 
c       Print*," p3mup3 ",  p3mup3 
c       Print*," p4mup2 ",  p4mup2 
c       Print*," p4mup3 ",  p4mup3 
c       Print*," mup2mup3", mup2mup3  
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FaFunctionABEmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3m
     -   up3,  p4mup2, p4mup3,  mup2mup3,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(3)
       SMG(2) = SMB(2)
       SMG(3) = SMB(3)
       SMG(4) = p2sq*SMB(2)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c************************************************************************************
c************************************************************************************
c    This should be equal to the boxline with momenta p1,p2+p3,p4,p5:
c************************************************************************************
c************************************************************************************
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = DCMPLX(B013R,B013I)*Is12*(2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3
     -   *SMG(3))+SMG(4))
       resultg(3) = p2mup3*F(1)*SMG(1)+Fa(1)*SMG(2)+mup2mup3*F(4)*SMG(
     -   3)+F(5)*SMG(4)
       resultg(4) = Fa(2)*SMG(1)+p1mup2*F(9)*SMG(2)+mup2mup3*F(10)*SMG
     -   (3)+F(11)*SMG(4)
       resultg(5) = Fa(3)*SMG(1)+Fa(4)*SMG(2)+Fa(5)*SMG(3)+F(27)*SMG(4
     -   )
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(1)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5))
       resultgaugeb(1) = Is12*(2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3
     -   *SMG(3))+SMG(4))
           resultgauge(2)=(0d0,0d0)      
           resultgaugeb(2)=(0d0,0d0)      
           return      
           endif      
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
       mup2mup3=dotrc(p3,mup2)
c       Print*," p1mup2 ",  p1mup2 
c       Print*," p1mup3 ",  p1mup3 
c       Print*," p2mup2 ",  p2mup2 
c       Print*," p2mup3 ",  p2mup3 
c       Print*," p3mup2 ",  p3mup2 
c       Print*," p3mup3 ",  p3mup3 
c       Print*," p4mup2 ",  p4mup2 
c       Print*," p4mup3 ",  p4mup3 
c       Print*," mup2mup3", mup2mup3  
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FaFunctionABEmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3m
     -   up3, p4mup2, p4mup3,   mup2mup3,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:3)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = -SMB(3)
       SMG(3) = SMB(3)
       SMG(4) = (p2sq+p3sq-s12-s23)*SMB(1)-2*(p2mup2+p4mup2)*SMB(3)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c************************************************************************************
c************************************************************************************
c    This should be equal to the subsctraction of boxline with momenta (p1,p2,p3+p4,p5) -(p1,p2+p3,p4,p5):
c************************************************************************************
c************************************************************************************
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = DCMPLX(B013R,B013I)*Is12*(2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3
     -   *SMG(3))+SMG(4))
       resultg(3) = p2mup3*F(1)*SMG(1)+Fa(1)*SMG(2)+mup2mup3*F(4)*SMG(
     -   3)+F(5)*SMG(4)
       resultg(4) = Fa(2)*SMG(1)+p1mup2*F(9)*SMG(2)+mup2mup3*F(10)*SMG
     -   (3)+F(11)*SMG(4)
       resultg(5) = Fa(3)*SMG(1)+Fa(4)*SMG(2)+Fa(5)*SMG(3)+F(27)*SMG(4
     -   )
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(1)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5))
       resultgaugeb(1) = Is12*(2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3
     -   *SMG(3))+SMG(4))
           resultgauge(2)=(0d0,0d0)      
           resultgaugeb(2)=(0d0,0d0)      
          return      
           endif      
c************************************************************************************
c************************************************************************************
      If(gaugetest.eq.3) then 
cc      The mu_p2 is replaced for the incoming moment 
c      mup2->p2 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2=dotrr(p1,p2)
       p2mup2=dotrr(p2,p2)
       p3mup2=dotrr(p3,p2)
       p4mup2=dotrr(p4,p2)
       mup2mup3=dotrc(p2,mup3)
c       Print*," p1mup2 ",  p1mup2 
c       Print*," p1mup3 ",  p1mup3 
c       Print*," p2mup2 ",  p2mup2 
c       Print*," p2mup3 ",  p2mup3 
c       Print*," p3mup2 ",  p3mup2 
c       Print*," p3mup3 ",  p3mup3 
c       Print*," p4mup2 ",  p4mup2 
c       Print*," p4mup3 ",  p4mup3 
c       Print*," mup2mup3", mup2mup3  
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FaFunctionABEmm(p1mup2, p1mup3,p2mup2, p2mup3,p3mup2, p3mu
     -   p3, p4mup2, p4mup3,mup2mup3,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(3)
       SMG(2) = SMB(2)
       SMG(3) = SMB(3)
       SMG(4) = p2sq*SMB(2)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c************************************************************************************
c************************************************************************************
c    This should be equal to the boxline with momenta p1,p2+p3,p4,p5:
c************************************************************************************
c************************************************************************************
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = DCMPLX(B013R,B013I)*Is12*(2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3
     -   *SMG(3))+SMG(4))
       resultg(3) = p2mup3*F(1)*SMG(1)+Fa(1)*SMG(2)+mup2mup3*F(4)*SMG(
     -   3)+F(5)*SMG(4)
       resultg(4) = Fa(2)*SMG(1)+p1mup2*F(9)*SMG(2)+mup2mup3*F(10)*SMG
     -   (3)+F(11)*SMG(4)
       resultg(5) = Fa(3)*SMG(1)+Fa(4)*SMG(2)+Fa(5)*SMG(3)+F(27)*SMG(4
     -   )
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(1)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5))
       resultgaugeb(1) = Is12*(2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3
     -   *SMG(3))+SMG(4))
c************************************************************************************
c************************************************************************************
c************************************************************************************
c      The mu_p3 is replaced for the incoming moment 
c      mup3->p3 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = dotrr(p1,p3)
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = dotrr(p2,p3)
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = dotrr(p3,p3)
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = dotrr(p4,p3)
       mup2mup3=dotrc(p3,mup2)
c       Print*," p1mup2 ",  p1mup2 
c       Print*," p1mup3 ",  p1mup3 
c       Print*," p2mup2 ",  p2mup2 
c       Print*," p2mup3 ",  p2mup3 
c       Print*," p3mup2 ",  p3mup2 
c       Print*," p3mup3 ",  p3mup3 
c       Print*," p4mup2 ",  p4mup2 
c       Print*," p4mup3 ",  p4mup3 
c       Print*," mup2mup3", mup2mup3  
c************************************************************************************
c************************************************************************************
c************** Calling the Fa functions**********************************************************************
       call FaFunctionABEmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3m
     -   up3,p4mup2, p4mup3, mup2mup3,Fa)
c************************************************************************************
c************************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:3)
c************************************************************************************
c************************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = -SMB(3)
       SMG(3) = SMB(3)
       SMG(4) = (p2sq+p3sq-s12-s23)*SMB(1)-2*(p2mup2+p4mup2)*SMB(3)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c************************************************************************************
c************************************************************************************
c    This should be equal to the subsctraction of boxline with momenta (p1,p2,p3+p4,p5) -(p1,p2+p3,p4,p5):
c************************************************************************************
c************************************************************************************
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F's and k's without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!


       resultg(2) = DCMPLX(B013R,B013I)*Is12*(2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3
     -   *SMG(3))+SMG(4))
       resultg(3) = p2mup3*F(1)*SMG(1)+Fa(1)*SMG(2)+mup2mup3*F(4)*SMG(
     -   3)+F(5)*SMG(4)
       resultg(4) = Fa(2)*SMG(1)+p1mup2*F(9)*SMG(2)+mup2mup3*F(10)*SMG
     -   (3)+F(11)*SMG(4)
       resultg(5) = Fa(3)*SMG(1)+Fa(4)*SMG(2)+Fa(5)*SMG(3)+F(27)*SMG(4
     -   )
c************************************************************************************
c** ****************
c************************************************************************************
       resultgauge(2)=-(resultg(2)+resultg(3)+resultg(4)+resultg(5))
       resultgaugeb(2) = Is12*(2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3
     -   *SMG(3))+SMG(4))
        Endif
c************************************************************************************
c************************************************************************************
c************************************************************************************
c************************************************************************************
       Return
       End

