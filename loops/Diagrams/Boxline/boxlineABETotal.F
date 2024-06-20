       subroutine boxlineABETotal(p1,p2,p3,p4,barpsi_p4,psi_p1,mup2,mup3,
     -   alpha,musqIn,ngluon,posgluon,gaugetest,comp,resultgauge,result
     -   ,resultgaugeb,resultb,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:25/7/2008

       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3)
       Complex*16 result(5)
       Complex*16 resultgauge(2),resultgaugeb(2),resultb 
       Real*8 musq,P(18) 
       Complex*16 F(27),K(1)
       Integer alpha,comp,gaugetest,ngluon,posgluon,Div
       common/FfunctionsABEmm/F
       common/PfunctionsABEmm/P
       common/KfunctionsABEmm/K
       SAVE/FfunctionsABEmm/
       SAVE/PfunctionsABEmm/
       SAVE/KfunctionsABEmm/
       Real*8 musqIn
       Common/musqInv/musq
       musq=musqIn      
       
       if(Div.eq.0)then
        call boxlineABEmm(p1,p2,p3,p4,barpsi_p4,psi_p1,mup2,mup3,
     -   alpha,musq,ngluon,posgluon,gaugetest,comp,resultgauge,result
     -   ,resultgaugeb,resultb) 
       else
        call boxlineABEmmDiv(p1,p2,p3,p4,barpsi_p4,psi_p1,mup2,mup3,
     -   alpha,musq,ngluon,posgluon,gaugetest,comp,resultgauge,result
     -   ,resultgaugeb,resultb,Div)
       endif


       End
