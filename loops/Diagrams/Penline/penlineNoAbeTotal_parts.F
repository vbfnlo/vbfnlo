       subroutine penlineNoAbe_pT(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mu
     -   p2,mup3,mup5,alpha,musqIn,ngluon,posgluon,gaugetest,comp,resul
     -   tgauge,result,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 31/05/2012
c Modified:31/05/2012
c************************************************************************************
      IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3),mup5(0:3)
       integer alpha,ngluon,posgluon,gaugetest,comp,Div
       real*8 musqIn,musq
       complex*16 result(6),resultgauge(3)
       real*8 P(759) 
       complex*16 F(438),K(7)
       common/PartsFfunctionsNOABEmmm/F
       common/PartsPfunctionsNOABEmmm/P
       common/PartsKfunctionsNOABEmmm/K
       SAVE/PartsFfunctionsNOABEmmm/
       SAVE/PartsPfunctionsNOABEmmm/
       SAVE/PartsKfunctionsNOABEmmm/
       common/PartsKinematicsNOABEmmm/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
       real*8 p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
       Common/musqInv/musq
       musq=musqIn

       If(Div.eq.0) then
        call penlineNoAbe_pmmm(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mu
     -   p2,mup3,mup5,alpha,musq,ngluon,posgluon,gaugetest,comp,resul
     -   tgauge,result)
       else
        call penlineNoAbe_pmmmDiv(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mu
     -   p2,mup3,mup5,alpha,musq,ngluon,posgluon,gaugetest,comp,resul
     -   tgauge,result,Div)
       endif

       End


       subroutine penlineNoAbe_pTQ(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mu
     -   p2,mup3,mup5,alpha,musqIn,ngluon,posgluon,gaugetest,comp,resul
     -   tgauge,result,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 31/05/2012
c Modified:31/05/2012
c************************************************************************************
      IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3),mup5(0:3)
       integer alpha,ngluon,posgluon,gaugetest,comp,Div
       real*8 musqIn,musq
       complex*16 result(6),resultgauge(3)
       real*8 P(759) 
       complex*16 F(438),K(7)
       common/PartsFfunctionsNOABEmmm/F
       common/PartsPfunctionsNOABEmmm/P
       common/PartsKfunctionsNOABEmmm/K
       SAVE/PartsFfunctionsNOABEmmm/
       SAVE/PartsPfunctionsNOABEmmm/
       SAVE/PartsKfunctionsNOABEmmm/
       common/PartsKinematicsNOABEmmm/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
       real*8 p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
       Common/musqInv/musq
c cccccccccccccccccccccccccccccccccccccc 
       common/pengaugeValue/t
       real*8 t,tdble,tquad,td
       logical Ward,ward_q
       complex*16 result_Q(6),resultgauge_Q(3)
       integer jtemp
       real*8 accuracy
       logical gaugetrue
       logical wardpAbe(2)
       common/wardpenNoAbe/wardpAbe

        
       tdble=1d9
       tquad=1d9
       
      
       if(gaugetest.lt.5.and.gaugetest.gt.0) then
       gaugetrue=.true.
       else
       gaugetrue=.false.
       endif
 
          musq=musqIn
          accuracy=1d-2

        call penlineNoAbe_pT(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mu
     -   p2,mup3,mup5,alpha,musq,ngluon,posgluon,gaugetest,comp,resul
     -   tgauge,result,Div)

       if(gaugetrue) then

       wardpAbe(1)=.false.
       wardpAbe(2)=.false.

       call pentNoAbe_gauge(Accuracy,p1,p2,p3,p4,p5,barpsi_p4,psi_p1,
     -  mup2,mup3,mup5,alpha,musq,
     - resultgauge,
     - td,Ward)

       tdble=td

       if(tdble+1.eq.tdble) then
          tdble=1d9
       endif   

       if(tdble.ne.tdble) then
          tdble=1d9
       endif   

c       print*,"Inside new"
       if(tdble.gt.GaugeAccuracyPenLnA) then
c       print*, "gauge Accu Dble:P1g", tdble
c       print*, "WardP Dble:1g", Ward
       else
       wardpAbe(1)=.true.
       wardpAbe(2)=.true.
       endif

#ifdef WITH_QUAD
       if(tdble.gt.GaugeAccuracyPenLnA) then
 
        call penlineNoAbe_pT_QUAD(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mu
     -   p2,mup3,mup5,alpha,musq,ngluon,posgluon,gaugetest,comp,resul
     -   tgauge_Q,result_Q,Div)


       call pentNoAbe_gauge_QUAD(Accuracy,p1,p2,p3,p4,p5,barpsi_p4,psi_p1,
     -  mup2,mup3,mup5,alpha,musq,
     - resultgauge_Q,
     - td,Ward_Q)

       tquad=td
       if(tquad+1.eq.tquad) then
	tquad = 1d9
	print*, "NaN in line", __LINE__, __FILE__
       endif

       if(tquad.ne.tquad) then
	tquad = 1d9
	print*, "NaN in line", __LINE__, __FILE__
       endif

       if(tquad+1.eq.tquad) then
          tquad=1d9
          print*,"QUAD:PenNo1g_p"
       endif   

       if(tquad.ne.tquad) then
          tquad=1d9
          print*,"QUAD:PenNo1g_p"
       endif   


       if(tdble.gt.GaugeAccuracyPenLnA) then
c       print*, "gauge Accu QUAD:P1g", tquad
c       print*, "Ward QUAD:1g", Ward_q
c       do jtemp=1,6
c       print*, "Dble", result(jtemp)
c       print*, "QUAD", result_Q(jtemp)
c       print*, "ratio:P1g", abs(result(jtemp)/result_Q(jtemp) -1)
c       enddo
       endif ! printing

c to rescue Dble precision
       if (tquad.lt.tdble) then


          if(tquad.lt.GaugeAccuracyPenLnA) wardpAbe(2)=.true.
     
         do jtemp=1,6
           result(jtemp)=result_Q(jtemp)
           if(jtemp.le.3) then
           resultgauge(jtemp)=resultgauge_Q(jtemp)
  
           endif
        enddo
       endif

       endif  ! tdble

#endif
       endif ! gaugetrue

       End


#ifdef WITH_QUAD
       subroutine penlineNoAbe_pT_QUAD(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mu
     -   p2,mup3,mup5,alpha,musqIn,ngluon,posgluon,gaugetest,comp,resul
     -   tgauge,result,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 31/05/2012
c Modified:31/05/2012
c************************************************************************************
      IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3),mup5(0:3)
       integer alpha,ngluon,posgluon,gaugetest,comp,Div
       real*8 musqIn,musq
       complex*16 result(6),resultgauge(3)
       real*8 P(759) 
       complex*16 F(438),K(7)
       common/PartsFfunctionsNOABEmmm/F
       common/PartsPfunctionsNOABEmmm/P
       common/PartsKfunctionsNOABEmmm/K
       SAVE/PartsFfunctionsNOABEmmm/
       SAVE/PartsPfunctionsNOABEmmm/
       SAVE/PartsKfunctionsNOABEmmm/
       common/PartsKinematicsNOABEmmm/p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
       real*8 p1sq,p2sq,p3sq,p4sq,p5sq,s12,s23,s34,s45,s15
       Common/musqInv/musq
       musq=musqIn

       If(Div.eq.0) then
        call penlineNoAbe_pmmm_QUAD(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mu
     -   p2,mup3,mup5,alpha,musq,ngluon,posgluon,gaugetest,comp,resul
     -   tgauge,result)
       else
       print*,"NO QUAD DIV so far:Pen1g"
        call penlineNoAbe_pmmmDiv(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mu
     -   p2,mup3,mup5,alpha,musq,ngluon,posgluon,gaugetest,comp,resul
     -   tgauge,result,Div)
       endif

       End
#endif
