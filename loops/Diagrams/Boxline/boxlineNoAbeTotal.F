       subroutine boxlineNoAbeTotal(p1,p2,p3,p4,barpsi_p3,psi_p1,mup2,mup
     -   4,alpha,musqIn,gaugetest,comp,resultgauge,result,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 28/7/2008
c Modified:4/8/2008
c*****************************************************************************
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3)
       Complex*16 barpsi_p3(2),psi_p1(2),mup2(0:3),mup4(0:3)
       Complex*16 F(17)
       Complex*16 result
       Complex*16 resultgauge(2) 
       Real*8 musq,P(35) 
       Integer alpha,comp,gaugetest,Div
       common/FfunctionsNoAbemm/F
       common/PfunctionsNoAbemm/P
       SAVE/FfunctionsNoAbemm/
       SAVE/PfunctionsNoAbemm/
       Real*8 musqIn
       Common/musqInv/musq
       musq=musqIn
       
       if(Div.eq.0) then
       call boxlineNoAbemm(p1,p2,p3,p4,barpsi_p3,psi_p1,mup2,mup
     -   4,alpha,musq,gaugetest,comp,resultgauge,result) 
          
       else
       call boxlineNoAbemmDiv(p1,p2,p3,p4,barpsi_p3,psi_p1,mup2,mup
     -   4,alpha,musq,gaugetest,comp,resultgauge,result,Div)
       endif


       End
