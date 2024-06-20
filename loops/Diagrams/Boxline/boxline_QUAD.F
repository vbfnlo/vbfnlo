       subroutine boxlinemm_QUAD(p1t,p2t,p3t,p4t,barpsi_p4,psi_p1,mup2,mup3,alpha
     -   ,musqIn,gaugetest,comp,resultgauge,result,resultgaugeb,resultb
     -   )
c ***********************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 13/2/2008
c Modified:20/6/2008
c ***********************************************************************
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
c anything it represents pp->2 photons up to the electromagnetic charge 
c (for uu->2 photons is missing (2/3)^2) from the coupling. For W and Z,
c we have to account for the Diracgamma_5 that are missing in the vertex
c Note: To make it shorter in the promgram: mu_p2,...->mup2,... 
c Notation of External momenta: p1+p2+p3+p4=0 
c mu_p2,mu_p3, should be think as external current 
c alpha is the helicity of the initial spinor 
c musq is the renormalization scale energy  
c gaugetest,integer value.Different gauge test: 
c gaugetest=0 should give zero: 
c resultgauge(1)=0.This is the result of the box-line 
c when making the replacement (mup2->p2+p3,mup3->p3)
c (For orientation, in a good point zero means 10^-12-10^-14
c ,but this condition is quite strong. 
c For example, only when this zero is of order 1*10^-3 
c in some of the boxes that goes in the penline gauge check
c we start to see problems, some points pass the test and others not
c with delta=0.1,delta=penli/box-1.
c Increasing the accuracy of this zero,  the acuracy of the gauge test
c increases too. Typically: 10^(-n)->delta=10^-(n+2)
c gaugetest=1 give the result after replacing mu_p2 with p2 
c gaugetest=2 give the result after replacing mu_p3 with p3 
c gaugetest=3 calculates the two(1,2) different gaugue tests 
c gaugetest>3 no calculation of gauge test in case you
c have already done it before(it safes times)
c comp: integer value.The first time called with p1...p4, comp=1
c ATTENTION: ONLY!!!If you have to call the subroutine consecutively 
c with the same arguments(p1,p2,p3,p4). Then, comp=-1 
c (it safes some lines of code) 
c OUTPUT:
c resultgauge and resultgaugeb(corresponding born amplitude contracted)
c are arrays of dimension two. 
c In case you use gaugetest=(0,2). The result is given in the first argument:
c resultgauge(1) and resultgaugeb(1).The argument 2 is set to zero 
c In case you use gaugetest=3. resultgauge(1) is the result of gaugetest=1,
c resultgauge(2) of gaugetest=2
c        result:result of the boxline 
c        resulb:result of the born amplitude 
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c    Declaration of variables 
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
c********************************************************************
       IMPLICIT NONE
       DOUBLE PRECISION   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       DOUBLE PRECISION   p1t(0:3),p2t(0:3),p3t(0:3),p4t(0:3)
       DOUBLE PRECISION p1sq, p1p2, p1p3, p2sq, p2p3, p2p4, p3sq, 
     -          p3p4, p4sq
       DOUBLE PRECISION   s12, s23
       DOUBLE COMPLEX p1mup2, p1mup3, p2mup2, p2mup3, p3mup2, p3mup3, 
     -          p4mup2, p4mup3
       DOUBLE COMPLEX mup2mup3
       DOUBLE PRECISION dotrr
       DOUBLE COMPLEX dotrc,dotcc,B0fin,C0fin,D02m_fin,D0fin
       EXTERNAL dotrr,dotrc,dotcc,B0fin,C0fin,D02m_fin,D0fin

       DOUBLE PRECISION B012R,B013R,B014R,B023R,B024R,B034R
       DOUBLE PRECISION B012I,B013I,B014I,B023I,B024I,B034I

       DOUBLE PRECISION C0123R,C0124R,C0134R,C0234R
       DOUBLE PRECISION C0123I,C0124I,C0134I,C0234I

       DOUBLE PRECISION Cij123R(4,2),Cij124R(4,2),Cij134R(4,2),Cij234R(4,2)
       DOUBLE PRECISION Cij123I(4,2),Cij124I(4,2),Cij134I(4,2),Cij234I(4,2)

       DOUBLE PRECISION D0R,D0I

       DOUBLE PRECISION DijR(13,3),DijI(13,3)
       DOUBLE COMPLEX SMB(4), SMG(4),F(16),Fa(3),K(1)
       DOUBLE PRECISION FR(16),FI(16),KR(1),KI(1)
       DOUBLE COMPLEX barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3) 
       DOUBLE COMPLEX SC1c,SC3ccc,SC3rcc,SC3rrc,SC1r,SC3rrr
       DOUBLE COMPLEX resultgauge(2),result,resultb,resultgaugeb(2) 
       DOUBLE PRECISION musq, P(18) ,musqIn
       EXTERNAL SC1c,SC1r,SC3ccc,SC3rcc,SC3rrc,SC3rrr
       INTEGER alpha,comp,gaugetest 
       COMMON/Ffunctionsmm/F,P
       COMMON/Kfunctionsmm/K
       SAVE /Ffunctionsmm/
       SAVE /Kfunctionsmm/
       Common/musqInv/musq

ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c      End QUAD precision
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

       COMMON/DVALFunctionsD/ D0R,DIJR,D0I,DIJI

       COMMON/CVALFunctionsD/ C0123R,CIJ123R,C0123I,CIJ123I, 
     &    C0124R,CIJ124R,C0124I,CIJ124I,
     &    C0134R,CIJ134R,C0134I,CIJ134I,
     &    C0234R,CIJ234R,C0234I,CIJ234I
  
       COMMON/BVALFunctionsD/ B012R,B012I, B013R,B013I, B014R,B014I, 
     &      B023R,B023I,  B024R,B024I,
     &    B034R,B034I
       COMMON/InvariantsD/p1sq,p2sq,p3sq,p4sq,s12,s23

       common/invariantsrest/p1p2
       integer jtemp
ccccccccccccccccc
       musq=musqIn
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c       Definition of the scalar products. Not inlcuded the contraction of the
c       moments with the external currents  
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
 
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

c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c       Definition of the F,P functions:Independent of the currents    
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
       P(1) = p2sq+p3sq-s12-s23
       FR(1) = 4*(Cij134R(1,1)+Cij134R(1,2)-Cij134R(2,1)+Cij134R(2,2)+
     -   p2sq*(DijR(2,1)-DijR(3,1)+DijR(4,2)-DijR(6,2))-2*(Cij134R(3,
     -   2)-DijR(7,2)-DijR(11,3)+DijR(13,3))+(DijR(3,2)-DijR(5,2))*P(
     -   1))
       FI(1) = 4*(Cij134I(1,1)+Cij134I(1,2)-Cij134I(2,1)+Cij134I(2,2)+
     -   p2sq*(DijI(2,1)-DijI(3,1)+DijI(4,2)-DijI(6,2))-2*(Cij134I(3,
     -   2)-DijI(7,2)-DijI(11,3)+DijI(13,3))+(DijI(3,2)-DijI(5,2))*P(
     -   1))
       F(1)=DCMPLX(FR(1),FI(1))
       P(2) = p2sq-s12
       P(3) = -p2sq+2*s12+s23
       P(4) = -p2sq+s12
       P(5) = -p2sq+3*s12+s23
       FR(2) = (2*(B013R+4*(Cij123R(4,2)+Cij134R(4,2))-2*(B023R+B034R-
     -   p2sq*Cij123R(1,1)-p3sq*(C0134R+Cij134R(1,1))-C0123R*P(2))))/
     -   s12+4*(-C0134R-Cij123R(1,1)-Cij123R(2,1)-Cij134R(1,1)+Cij134
     -   R(1,2)+Cij134R(2,2)-Cij234R(1,1)+Cij234R(2,1)+s23*(D0R+DijR(
     -   1,1))-s12*DijR(2,2)-2*(Cij134R(3,2)+DijR(7,2)-DijR(12,3)+Dij
     -   R(13,3))+DijR(2,1)*P(2)-DijR(3,1)*P(2)+DijR(4,2)*P(2)-DijR(3
     -   ,2)*P(3)+DijR(5,2)*P(4)+DijR(6,2)*P(5))
       FI(2) = (2*(B013I+4*(Cij123I(4,2)+Cij134I(4,2))-2*(B023I+B034I-
     -   p2sq*Cij123I(1,1)-p3sq*(C0134I+Cij134I(1,1))-C0123I*P(2))))/
     -   s12+4*(-C0134I-Cij123I(1,1)-Cij123I(2,1)-Cij134I(1,1)+Cij134
     -   I(1,2)+Cij134I(2,2)-Cij234I(1,1)+Cij234I(2,1)+s23*(D0I+DijI(
     -   1,1))-s12*DijI(2,2)-2*(Cij134I(3,2)+DijI(7,2)-DijI(12,3)+Dij
     -   I(13,3))+DijI(2,1)*P(2)-DijI(3,1)*P(2)+DijI(4,2)*P(2)-DijI(3
     -   ,2)*P(3)+DijI(5,2)*P(4)+DijI(6,2)*P(5))
       F(2)=DCMPLX(FR(2),FI(2))
       P(6) = s12+s23
       P(7) = p2sq-2*s12-s23
       FR(3) = 4*(-C0134R+C0234R-Cij134R(1,1)+Cij134R(2,2)-Cij134R(3,2
     -   )+Cij234R(2,1)+s12*(DijR(2,1)+DijR(6,2))-2*(DijR(7,2)+DijR(1
     -   3,3))-DijR(5,2)*P(2)-DijR(3,1)*P(6)+DijR(3,2)*P(7))
       FI(3) = 4*(-C0134I+C0234I-Cij134I(1,1)+Cij134I(2,2)-Cij134I(3,2
     -   )+Cij234I(2,1)+s12*(DijI(2,1)+DijI(6,2))-2*(DijI(7,2)+DijI(1
     -   3,3))-DijI(5,2)*P(2)-DijI(3,1)*P(6)+DijI(3,2)*P(7))
       F(3)=DCMPLX(FR(3),FI(3))
       P(8) = p2sq-s12-s23
       FR(4) = 4*(-C0234R+s12*(-DijR(2,1)+DijR(3,1))+p3sq*DijR(3,2)+2*
     -   (DijR(7,2)-DijR(13,3))-DijR(6,2)*P(1)+DijR(5,2)*P(8))
       FI(4) = 4*(-C0234I+s12*(-DijI(2,1)+DijI(3,1))+p3sq*DijI(3,2)+2*
     -   (DijI(7,2)-DijI(13,3))-DijI(6,2)*P(1)+DijI(5,2)*P(8))
       F(4)=DCMPLX(FR(4),FI(4))
       P(9) = p3sq-s12
       P(10) = P(9)/s12
       P(11) = p3sq-2*s12
       P(12) = P(11)/s12
       P(13) = p3sq-2*s12-s23
       P(14) = s12-s23
       FR(5) = (2*(B013R-2*(B023R+B034R-p2sq*(C0123R+Cij123R(1,1)))+4*
     -   (Cij123R(4,2)+Cij134R(4,2))))/s12+4*(Cij123R(3,2)+Cij134R(2,
     -   1)+Cij234R(2,1)+s23*(D0R+2*DijR(1,1)+DijR(1,2))-s12*DijR(3,2
     -   )+2*(DijR(11,3)-DijR(13,3))+C0134R*P(10)-DijR(3,1)*P(11)+Cij
     -   134R(1,1)*P(12)+(DijR(2,1)+DijR(4,2)-DijR(6,2))*P(13)+DijR(5
     -   ,2)*P(14))
       FI(5) = (2*(B013I-2*(B023I+B034I-p2sq*(C0123I+Cij123I(1,1)))+4*
     -   (Cij123I(4,2)+Cij134I(4,2))))/s12+4*(Cij123I(3,2)+Cij134I(2,
     -   1)+Cij234I(2,1)+s23*(D0I+2*DijI(1,1)+DijI(1,2))-s12*DijI(3,2
     -   )+2*(DijI(11,3)-DijI(13,3))+C0134I*P(10)-DijI(3,1)*P(11)+Cij
     -   134I(1,1)*P(12)+(DijI(2,1)+DijI(4,2)-DijI(6,2))*P(13)+DijI(5
     -   ,2)*P(14))
       F(5)=DCMPLX(FR(5),FI(5))
       P(15) = -p2sq+s12+s23
       P(16) = p2sq+2*p3sq-s12-s23
       FR(6) = 4*(-C0234R+Cij123R(2,1)+Cij123R(2,2)+p3sq*DijR(3,2)+4*D
     -   ijR(7,2)+2*(DijR(12,3)-DijR(13,3))+DijR(2,2)*P(1)+DijR(5,2)*
     -   P(8)+(DijR(2,1)-DijR(3,1))*P(9)+DijR(4,2)*P(15)-DijR(6,2)*P(
     -   16))
       FI(6) = 4*(-C0234I+Cij123I(2,1)+Cij123I(2,2)+p3sq*DijI(3,2)+4*D
     -   ijI(7,2)+2*(DijI(12,3)-DijI(13,3))+DijI(2,2)*P(1)+DijI(5,2)*
     -   P(8)+(DijI(2,1)-DijI(3,1))*P(9)+DijI(4,2)*P(15)-DijI(6,2)*P(
     -   16))
       F(6)=DCMPLX(FR(6),FI(6))
       FR(7) = 8*(DijR(2,1)-DijR(3,1)-DijR(3,3)+DijR(4,3)-DijR(5,3)+Di
     -   jR(9,3)+2*(DijR(3,2)+DijR(4,2)-DijR(5,2)-DijR(6,2)+DijR(7,3)
     -   -DijR(10,3)))
       FI(7) = 8*(DijI(2,1)-DijI(3,1)-DijI(3,3)+DijI(4,3)-DijI(5,3)+Di
     -   jI(9,3)+2*(DijI(3,2)+DijI(4,2)-DijI(5,2)-DijI(6,2)+DijI(7,3)
     -   -DijI(10,3)))
       F(7)=DCMPLX(FR(7),FI(7))
       FR(8) = 8*(DijR(2,2)+DijR(3,2)-DijR(3,3)+DijR(6,3)+DijR(7,3)-Di
     -   jR(8,3)-2*(DijR(6,2)-DijR(9,3)+DijR(10,3)))
       FI(8) = 8*(DijI(2,2)+DijI(3,2)-DijI(3,3)+DijI(6,3)+DijI(7,3)-Di
     -   jI(8,3)-2*(DijI(6,2)-DijI(9,3)+DijI(10,3)))
       F(8)=DCMPLX(FR(8),FI(8))
       FR(9) = 8*(-DijR(3,3)+DijR(5,2)-DijR(6,2)+DijR(7,3)+DijR(9,3)-D
     -   ijR(10,3))
       FI(9) = 8*(-DijI(3,3)+DijI(5,2)-DijI(6,2)+DijI(7,3)+DijI(9,3)-D
     -   ijI(10,3))
       F(9)=DCMPLX(FR(9),FI(9))
       FR(10) = 8*(DijR(2,1)+DijR(2,2)-DijR(3,1)-DijR(3,3)+DijR(4,2)-D
     -   ijR(5,2)-3*DijR(6,2)+DijR(6,3)+DijR(7,3)-DijR(8,3)+2*(DijR(3
     -   ,2)+DijR(9,3)-DijR(10,3)))
       FI(10) = 8*(DijI(2,1)+DijI(2,2)-DijI(3,1)-DijI(3,3)+DijI(4,2)-D
     -   ijI(5,2)-3*DijI(6,2)+DijI(6,3)+DijI(7,3)-DijI(8,3)+2*(DijI(3
     -   ,2)+DijI(9,3)-DijI(10,3)))
       F(10)=DCMPLX(FR(10),FI(10))
       FR(11) = 8*(DijR(2,2)+DijR(2,3)+DijR(3,2)-DijR(3,3)-2*DijR(6,2)
     -   -3*(DijR(8,3)-DijR(9,3)))
       FI(11) = 8*(DijI(2,2)+DijI(2,3)+DijI(3,2)-DijI(3,3)-2*DijI(6,2)
     -   -3*(DijI(8,3)-DijI(9,3)))
       F(11)=DCMPLX(FR(11),FI(11))
       FR(12) = -8*(DijR(2,1)-DijR(3,1)-DijR(3,2)+DijR(3,3)+DijR(4,2)-
     -   DijR(7,3)-DijR(9,3)+DijR(10,3))
       FI(12) = -8*(DijI(2,1)-DijI(3,1)-DijI(3,2)+DijI(3,3)+DijI(4,2)-
     -   DijI(7,3)-DijI(9,3)+DijI(10,3))
       F(12)=DCMPLX(FR(12),FI(12))
       FR(13) = -8*(DijR(2,1)+DijR(2,2)-DijR(3,1)+DijR(3,3)-DijR(6,2)+
     -   DijR(8,3)-2*DijR(9,3))
       FI(13) = -8*(DijI(2,1)+DijI(2,2)-DijI(3,1)+DijI(3,3)-DijI(6,2)+
     -   DijI(8,3)-2*DijI(9,3))
       F(13)=DCMPLX(FR(13),FI(13))
       FR(14) = 8*(-DijR(3,2)-DijR(3,3)+DijR(6,2)+DijR(9,3))
       FI(14) = 8*(-DijI(3,2)-DijI(3,3)+DijI(6,2)+DijI(9,3))
       F(14)=DCMPLX(FR(14),FI(14))
       FR(15) = -4*(-C0134R-Cij123R(1,1)-Cij123R(2,1)+Cij134R(2,1)+s23
     -   *(D0R+DijR(1,1))-2*(Cij134R(1,1)+DijR(7,2)+DijR(12,3)-DijR(1
     -   3,3)))-(2*(B013R+4*(Cij123R(4,2)+Cij134R(4,2))-2*(B023R+B034
     -   R-p2sq*Cij123R(1,1)-p3sq*(C0134R+Cij134R(1,1))-C0123R*P(2)))
     -   )/s12
       FI(15) = -4*(-C0134I-Cij123I(1,1)-Cij123I(2,1)+Cij134I(2,1)+s23
     -   *(D0I+DijI(1,1))-2*(Cij134I(1,1)+DijI(7,2)+DijI(12,3)-DijI(1
     -   3,3)))-(2*(B013I+4*(Cij123I(4,2)+Cij134I(4,2))-2*(B023I+B034
     -   I-p2sq*Cij123I(1,1)-p3sq*(C0134I+Cij134I(1,1))-C0123I*P(2)))
     -   )/s12
       F(15)=DCMPLX(FR(15),FI(15))
       P(17) = p2sq+s12
       P(18) = p3sq+s12
       FR(16) = (B013R+4*(Cij123R(4,2)+Cij134R(4,2))-2*(B023R+B034R-p2
     -   sq*Cij123R(1,1)-p3sq*(C0134R+Cij134R(1,1))-C0123R*P(2)))/s12
     -   +2*(-C0134R+C0234R-Cij123R(1,1)-Cij123R(2,1)-2*Cij134R(1,1)+
     -   Cij134R(2,1)-Cij234R(1,1)+Cij234R(2,1)+s23*(D0R+DijR(1,1))-p
     -   2sq*DijR(5,2)+(p2sq+p3sq)*DijR(6,2)+s12*(DijR(5,2)+2*DijR(6,
     -   2))-6*DijR(7,2)+DijR(4,2)*P(2)-DijR(2,2)*P(17)-DijR(3,2)*P(1
     -   8))
       FI(16) = (B013I+4*(Cij123I(4,2)+Cij134I(4,2))-2*(B023I+B034I-p2
     -   sq*Cij123I(1,1)-p3sq*(C0134I+Cij134I(1,1))-C0123I*P(2)))/s12
     -   +2*(-C0134I+C0234I-Cij123I(1,1)-Cij123I(2,1)-2*Cij134I(1,1)+
     -   Cij134I(2,1)-Cij234I(1,1)+Cij234I(2,1)+s23*(D0I+DijI(1,1))-p
     -   2sq*DijI(5,2)+(p2sq+p3sq)*DijI(6,2)+s12*(DijI(5,2)+2*DijI(6,
     -   2))-6*DijI(7,2)+DijI(4,2)*P(2)-DijI(2,2)*P(17)-DijI(3,2)*P(1
     -   8))
       F(16)=DCMPLX(FR(16),FI(16))
       KR(1) = -DijR(3,3)-DijR(8,3)+2*DijR(9,3)
       KI(1) = -DijI(3,3)-DijI(8,3)+2*DijI(9,3)
       K(1)=DCMPLX(KR(1),KI(1))
       else
         p1p2=dotrr(p1,p2)
         p1p3=dotrr(p1,p3)
         p2p3=dotrr(p2,p3)
         p2p4=dotrr(p2,p4)
         p3p4=dotrr(p3,p4)

       endif


c               PART THAT DEPENDS ON THE EXTERNAL CURRENT
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c Computation of Fa fucntions. Depends on the external currents, through
c the contraction of the moments with the currents
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
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
c************** Calling the Fa functions***************************************
c******************************************************************************
c******************************************************************************
       call FaFunctionmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3mup3
     -   , p4mup2, p4mup3,mup2mup3,Fa)
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c       Definition of the Matrix Element  
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
       SMB(1) = SC1c(barpsi_p4,mup2,psi_p1,alpha)
       SMB(2) = SC1c(barpsi_p4,mup3,psi_p1,alpha)
       SMB(3) = SC1r(barpsi_p4,p2,psi_p1,alpha)
       SMB(4) = SC3rcc(barpsi_p4,p2,mup2,mup3,psi_p1,alpha)
c$$$       Print*," SMB(1) ", SMB(1)
c$$$       Print*," SMB(2) ", SMB(2)
c$$$       Print*," SMB(3) ", SMB(3)
c$$$       Print*," SMB(4) ", SMB(4)
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c       Amplitude                         
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************


       result = Fa(1)*SMB(1)+Fa(2)*SMB(2)+Fa(3)*SMB(3)+F(16)*SMB(4)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F"s and k"s without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       result =-result
       resultb = (2*(p2mup3*SMB(1)+p1mup2*SMB(2)-mup2mup3*SMB(3))+SMB(
     -   4))/s12

c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c                         CHECK OF GAUGE INVARIANCE             
      If(gaugetest.eq.0) then 
cc      The mu_p2,mu_p3,mu_p4 must be replaced for the incoming moment 
c      mup2->p2+p3,mup3->p3. SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2 = p1p2+p1p3
       p1mup3 = p1p3
       p2mup2 = p2sq+p2p3
       p2mup3 = p2p3
       p3mup2 = p2p3+p3sq
       p3mup3 = p3sq
       p4mup2 = p2p4+p3p4
       p4mup3 = p3p4
       mup2mup3=p2p3+p3sq
c       Print*," p1mup2 ",  p1mup2 
c       Print*," p1mup3 ",  p1mup3 
c       Print*," p2mup2 ",  p2mup2 
c       Print*," p2mup3 ",  p2mup3 
c       Print*," p3mup2 ",  p3mup2 
c       Print*," p3mup3 ",  p3mup3 
c       Print*," p4mup2 ",  p4mup2 
c       Print*," p4mup3 ",  p4mup3 
c       Print*," mup2mup3", mup2mup3  
c******************************************************************************
c************** Calling the Fa functions***************************************
       call FaFunctionmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3mup3
     -   , p4mup2, p4mup3,mup2mup3,Fa)
c******************************************************************************
c******************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2,4)
c    and mu_p3=p2+p3+p4
c******************************************************************************
c******************************************************************************
       SMG(1) = 0
       SMG(2) = -SMB(3)
       SMG(3) = SMB(3)
       SMG(4) = (-p2sq+p3sq)*SMB(3)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c******************************************************************************
c******************************************************************************
c    This should be zero:
c******************************************************************************
c******************************************************************************
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+F(16)*S
     -   MG(4)
       resultgauge(1) =resultgauge(1)/(Sqrt(abs(p3sq)*abs(s23)))
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F"s and k"s without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps
c   Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
       resultgaugeb(1) = (2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3*SMG(
     -   3))+SMG(4))/s12
       resultgaugeb(1) =resultgaugeb(1)/(Sqrt(abs(p3sq)*abs(s23)))
           resultgauge(2)=(0d0,0d0)      
           resultgaugeb(2)=(0d0,0d0)      
           return      
           endif      
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
      If(gaugetest.eq.1) then
cc      The mu_p2 is replaced for the incoming moment 
c      mup2->p2 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2=p1p2
       p2mup2=p2sq
       p3mup2=p2p3
       p4mup2=p2p4
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
c******************************************************************************
c******************************************************************************
c***************calling the Fa functions***************************************
       call FaFunctionmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3mup3
     -   , p4mup2, p4mup3,mup2mup3,Fa)
c******************************************************************************
c******************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2)
c******************************************************************************
c******************************************************************************
       SMG(1) = SMB(3)
       SMG(2) = SMB(2)
       SMG(3) = SMB(3)
       SMG(4) = p2sq*SMB(2)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c******************************************************************************
c******************************************************************************
c    This should be equal to vertex with momenta p1,p2+p3,p4:
c******************************************************************************
c******************************************************************************
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+F(16)*S
     -   MG(4)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F"s and k"s without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps
c    Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
       resultgaugeb(1) = (2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3*SMG(
     -   3))+SMG(4))/s12
           resultgauge(2)=(0d0,0d0)      
           resultgaugeb(2)=(0d0,0d0)      
           return      
           endif      
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
      If(gaugetest.eq.2) then 
cc      The mu_p3 is replaced for the incoming moment 
c      mup3->p3 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup3=p1p3
       p2mup3=p2p3
       p3mup3=p3sq
       p4mup3=p3p4
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
c       Print*," mup2mup4",mup2mup4   
c       Print*,"mup3mup4", mup3mup4  
c******************************************************************************
c******************************************************************************
c************** Calling the Fa functions***************************************
       call FaFunctionmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3mup3
     -   , p4mup2, p4mup3,mup2mup3,Fa)
c******************************************************************************
c******************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:3)
c******************************************************************************
c******************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = -SMB(3)
       SMG(3) = SMB(3)
       SMG(4) = (p2sq+p3sq-s12-s23)*SMB(1)-2*(p2mup2+p4mup2)*SMB(3)
c       Print*," SMG(1) ", SMG(1)
c       Print*," SMG(2) ", SMG(2)
c       Print*," SMG(3) ", SMG(3)
c       Print*," SMG(4) ", SMG(4)
c******************************************************************************
c******************************************************************************
c This should be equal to the vertex line(graph) with momenta-(p1,p2+p3,p4):*
c******************************************************************************
c******************************************************************************
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+F(16)*S
     -   MG(4)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F"s and k"s without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps 
c  Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
       resultgaugeb(1) = (2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3*SMG(
     -   3))+SMG(4))/s12
           resultgauge(2)=(0d0,0d0)      
           resultgaugeb(2)=(0d0,0d0)      
          return      
           endif      
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
      If(gaugetest.eq.3) then 
cc      The mu_p2 is replaced for the incoming moment 
c      mup2->p2 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2=p1p2
       p2mup2=p2sq
       p3mup2=p2p3
       p4mup2=p2p4
       mup2mup3=dotrc(p2,mup3)
cfc       Print*," p1mup2 ",  p1mup2 
cfc       Print*," p1mup3 ",  p1mup3 
cfc       Print*," p2mup2 ",  p2mup2 
cfc       Print*," p2mup3 ",  p2mup3 
cfc       Print*," p3mup2 ",  p3mup2 
cfc       Print*," p3mup3 ",  p3mup3 
cfc       Print*," p4mup2 ",  p4mup2 
cfc       Print*," p4mup3 ",  p4mup3 
cfc       Print*," mup2mup3", mup2mup3  
c******************************************************************************
c******************************************************************************
c************** Calling the Fa functions***************************************
       call FaFunctionmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3mup3
     -   , p4mup2, p4mup3,mup2mup3,Fa)
c******************************************************************************
c******************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:2)
c******************************************************************************
c******************************************************************************
       SMG(1) = SMB(3)
       SMG(2) = SMB(2)
       SMG(3) = SMB(3)
       SMG(4) = p2sq*SMB(2)
cfc       Print*," SMG(1) ", SMG(1)
cfc       Print*," SMG(2) ", SMG(2)
cfc       Print*," SMG(3) ", SMG(3)
cfc       Print*," SMG(4) ", SMG(4)
c******************************************************************************
c******************************************************************************
c    This should be equal to vertex with momenta p1,p2+p3,p4:
c******************************************************************************
c******************************************************************************
       resultgauge(1) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+F(16)*S
     -   MG(4)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F"s and k"s without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps 
c   Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(1) =-resultgauge(1)
       resultgaugeb(1) = (2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3*SMG(
     -   3))+SMG(4))/s12

cfc       print*, "Fa",Fa(1:3)
cfc       print*, "F",F(16)
cfc       print*, "QUAD",resultgauge(1), Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+F(16)*SMG(4)
cfc       print*, "QUAD",resultgaugeb(1)
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c******************************************************************************
c      The mu_p3 is replaced for the incoming moment 
c      mup3->p3 SO:
c     I define the different contraction of the external
c      momenta with the currents making those replacements
       p1mup2 = dotrc(p1,mup2)
       p1mup3 = p1p3
       p2mup2 = dotrc(p2,mup2)
       p2mup3 = p2p3
       p3mup2 = dotrc(p3,mup2)
       p3mup3 = p3sq
       p4mup2 = dotrc(p4,mup2)
       p4mup3 = p3p4
       mup2mup3=dotrc(p3,mup2)
cfc       Print*," p1mup2 ",  p1mup2 
cfc       Print*," p1mup3 ",  p1mup3 
cfc       Print*," p2mup2 ",  p2mup2 
cfc       Print*," p2mup3 ",  p2mup3 
cfc       Print*," p3mup2 ",  p3mup2 
cfc       Print*," p3mup3 ",  p3mup3 
cfc       Print*," p4mup2 ",  p4mup2 
cfc       Print*," p4mup3 ",  p4mup3 
cfc       Print*," mup2mup3", mup2mup3  
c******************************************************************************
c******************************************************************************
c************** Calling the Fa functions***************************************
       call FaFunctionmm(p1mup2, p1mup3, p2mup2, p2mup3,p3mup2, p3mup3
     -   , p4mup2, p4mup3,mup2mup3,Fa)
c******************************************************************************
c******************************************************************************
c    Here, the matrix element with the replacement mu_p_j->pj,(j:3)
c******************************************************************************
c*****************************************************************************
       SMG(1) = SMB(1)
       SMG(2) = -SMB(3)
       SMG(3) = SMB(3)
       SMG(4) = (p2sq+p3sq-s12-s23)*SMB(1)-2*(p2mup2+p4mup2)*SMB(3)
cfc       Print*, "Gauge 2"
cfc       Print*," SMG(1) ", SMG(1)
cfc       Print*," SMG(2) ", SMG(2)
cfc       Print*," SMG(3) ", SMG(3)
cfc       Print*," SMG(4) ", SMG(4)
c******************************************************************************
c******************************************************************************
c    This should be equal to the vertex line  with momenta  -(p1,p2+p3,p4):
c******************************************************************************
c******************************************************************************
       resultgauge(2) = Fa(1)*SMG(1)+Fa(2)*SMG(2)+Fa(3)*SMG(3)+F(16)*S
     -   MG(4)
c The Finite virtual piece should be multiplied for (-1)  since 
c I have multiplied by (-I) to get the F"s and k"s without (I) factor
c . The factorization from the B_ij is Fact=(I/(4\[Pi])^2 (4 \[Pi])^Eps 
c   Gamma[1+Eps] (musq)^(-Eps))
c  c So, I*I=(-1)!!!
       resultgauge(2) =-resultgauge(2)
       resultgaugeb(2) = (2*(p2mup3*SMG(1)+p1mup2*SMG(2)-mup2mup3*SMG(
     -   3))+SMG(4))/s12
       
           Return      
           endif      
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
c*****************************************************************************
       Return
       End
c************** Calling the Fa functions**************************************



