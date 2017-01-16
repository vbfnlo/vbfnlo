       subroutine penlineABETotal(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,m
     -   up3,mup4,alpha,musqIn,ngluon,posgluon,gaugetest,comp,resultgau
     -   ge,result,resultgaugeb,resultb,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:30/7/2008
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p5(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3) 
       Complex*16 result(8)
       Complex*16 resultgauge(3),resultgaugeb(3),resultb 
       Real*8 musq,musqIn,P(554)
       Integer alpha,comp,gaugetest,ngluon,posgluon, Div
       Complex*16 F(448),K(4)
       common/Ffunctionsmmm/F
       common/Pfunctionsmmm/P
       common/Kfunctionsmmm/K
       SAVE/Ffunctionsmmm/
       SAVE/Pfunctionsmmm/
       SAVE/Kfunctionsmmm/
       Common/musqInv/musq
       musq=musqIn
c************************************************************************************
       If(Div.eq.0) then
       call penlineABEmmm(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,m
     -   up3,mup4,alpha,musq,ngluon,posgluon,gaugetest,comp,resultgau
     -   ge,result,resultgaugeb,resultb) 
       else
       call penlineABEmmmDiv(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,m
     -   up3,mup4,alpha,musq,ngluon,posgluon,gaugetest,comp,resultgau
     -   ge,result,resultgaugeb,resultb,Div) 

       endif

       End


       subroutine penlineABETotalQ(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,m
     -   up3,mup4,alpha,musqIn,ngluon,posgluon,gaugetest,comp,resultgau
     -   ge,result,resultgaugeb,resultb,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:30/7/2008
       IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p5(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3) 
       Complex*16 result(8)
       Complex*16 resultgauge(3),resultgaugeb(3),resultb 
       Real*8 musq,musqIn,P(554)
       Integer alpha,comp,gaugetest,ngluon,posgluon, Div
       Complex*16 F(448),K(4)
       common/Ffunctionsmmm/F
       common/Pfunctionsmmm/P
       common/Kfunctionsmmm/K
       SAVE/Ffunctionsmmm/
       SAVE/Pfunctionsmmm/
       SAVE/Kfunctionsmmm/
       Common/musqInv/musq

      common/pengaugeValue/t
       real*8 t,tdble,tquad
       logical Ward,ward_q
       Complex*16 result_Q(8)
       Complex*16 resultgauge_Q(3),resultgaugeb_Q(3),resultb_Q
       integer jtemp
       logical gaugetrue
       logical wardpAbe(2)
       common/wardpenAbe/wardpAbe
      
       musq=musqIn
c************************************************************************************

       if(gaugetest.lt.5.and.gaugetest.gt.0) then
       gaugetrue=.true.
       else
       gaugetrue=.false.
       endif
       
       call  penlineABETotal(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,m
     -   up3,mup4,alpha,musq,ngluon,posgluon,gaugetest,comp,resultgau
     -   ge,result,resultgaugeb,resultb,Div)

       if(gaugetrue) then

       wardpAbe(1)=.false.
       wardpAbe(2)=.false.

       call pent_gauge(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,
     -  mup2,mup3,mup4,alpha,musq,resultgauge,resultgaugeb,Ward)

       tdble=t

       if(tdble.gt.GaugeAccuracyPenL) then
c        print*, "gauge Accu Dble:P", t
c       print*,  "WardP Dble", Ward
       else
       wardpAbe(1)=.true.
       wardpAbe(2)=.true.
       endif

#ifdef WITH_QUAD

       if(tdble.gt.GaugeAccuracyPenL) then
 
       call  penlineABETotal_QUAD(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,m
     -   up3,mup4,alpha,musq,ngluon,posgluon,gaugetest,comp,resultgau
     -   ge_Q,result_Q,resultgaugeb_Q,resultb_Q,Div)

       call pent_gauge_QUAD(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,
     -  mup2,mup3,mup4,alpha,musq,resultgauge_Q,resultgaugeb_Q,Ward_q)

       tquad=t

       if(tdble.gt.GaugeAccuracyPenL) then
c          print*, "COMP",COMP
c          print*, "gaugetest",gaugetest
c       print*, "gauge Accu QUAD:P", tquad
c       print*, "Ward QUAD", Ward_q
c       do jtemp=1,8
c       print*, "Dble", resultf(jtemp)
c       print*, "QUAD", resultf_Q(jtemp)
c       print*, "ratio:P", abs(result(jtemp)/result_Q(jtemp) -1)
c       enddo
       endif ! printing

c to rescue Dble precision
       if (tquad.lt.tdble) then
         if(tquad.lt.GaugeAccuracyPenL) wardpAbe(2)=.true.

          resultb= resultb_Q
          
        do jtemp=1,8
           result(jtemp)=result_Q(jtemp)
           if(jtemp.le.3) then
           resultgauge(jtemp)=resultgauge_Q(jtemp)
           resultgaugeb(jtemp)  = resultgaugeb_Q(jtemp)
           endif
        enddo
       endif


       endif  ! tdble
#endif


       endif ! gaugetrue


       End

#ifdef WITH_QUAD


       subroutine penlineABETotal_QUAD(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,m
     -   up3,mup4,alpha,musqIn,ngluon,posgluon,gaugetest,comp,resultgau
     -   ge,result,resultgaugeb,resultb,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:30/7/2008
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p5(2),psi_p1(2),mup2(0:3),mup3(0:3),mup4(0:3) 
       Complex*16 result(8)
       Complex*16 resultgauge(3),resultgaugeb(3),resultb 
       Real*8 musq,musqIn,P(554)
       Integer alpha,comp,gaugetest,ngluon,posgluon, Div
       Complex*16 F(448),K(4)
       common/Ffunctionsmmm/F
       common/Pfunctionsmmm/P
       common/Kfunctionsmmm/K
       SAVE/Ffunctionsmmm/
       SAVE/Pfunctionsmmm/
       SAVE/Kfunctionsmmm/
       Common/musqInv/musq
       musq=musqIn
c************************************************************************************
       If(Div.eq.0) then
       call penlineABEmmm_QUAD(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,m
     -   up3,mup4,alpha,musq,ngluon,posgluon,gaugetest,comp,resultgau
     -   ge,result,resultgaugeb,resultb) 
       else
       print*,"NO QUAD DIV so far:PenAbe"
       call penlineABEmmmDiv(p1,p2,p3,p4,p5,barpsi_p5,psi_p1,mup2,m
     -   up3,mup4,alpha,musq,ngluon,posgluon,gaugetest,comp,resultgau
     -   ge,result,resultgaugeb,resultb,Div) 

       endif

       End
#endif
