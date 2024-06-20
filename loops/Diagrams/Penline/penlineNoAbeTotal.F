       subroutine penlineNoAbeTotal(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mup2
     -   ,mup3,mup5,alpha,musqIn,gaugetest,comp,resultgauge,result,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 3/7/2008
c Modified:4/8/2008
c ******************************************************************************
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3),mup5(0:3) 
       Complex*16 result
       Complex*16 resultgauge(3),F(280)
       Real*8 musq
       Integer alpha,comp,gaugetest,Div
       Real*8 musqIn,P(685) 

       common/FfunctionsNOABEmmm/F,P
       SAVE/FfunctionsNOABEmmm/
       Common/musqInv/musq
       musq=musqIn
       If(Div.eq.0) then

       call penlineNoAbemmm(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mup2
     -   ,mup3,mup5,alpha,musq,gaugetest,comp,resultgauge,result)

       else
       call penlineNoAbemmmDiv(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mup2
     -   ,mup3,mup5,alpha,musq,gaugetest,comp,resultgauge,result,Div)

       endif

       End


       subroutine penlineNoAbeTotalQ(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mup2
     -   ,mup3,mup5,alpha,musqIn,gaugetest,comp,resultgauge,result,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 3/7/2008
c Modified:4/8/2008
c ******************************************************************************
       IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3),mup5(0:3) 
       Complex*16 result
       Complex*16 resultgauge(3),F(280)
       Real*8 musq
       Integer alpha,comp,gaugetest,Div
       Real*8 musqIn,P(685) 

       common/FfunctionsNOABEmmm/F,P
       SAVE/FfunctionsNOABEmmm/
       Common/musqInv/musq
c cccccccccccccccccccccccccccccccccccccc 
       common/pengaugeValue/t
       real*8 t,tdble,tquad,td
       logical Ward,ward_q
       complex*16 result_Q,resultgauge_Q(3)
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

       call penlineNoAbeTotal(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mup2
     -   ,mup3,mup5,alpha,musq,gaugetest,comp,resultgauge,result,Div)


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

       call penlineNoAbeTotal_QUAD(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mup2
     -   ,mup3,mup5,alpha,musq,gaugetest,comp,resultgauge_Q,result_Q,Div)

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
          print*,"QUAD:PenNo1g"
       endif   

       if(tquad.ne.tquad) then
          tquad=1d9
          print*,"QUAD:PenNo1g"
       endif   


       if(tdble.gt.GaugeAccuracyPenLnA) then
c       print*, "gauge Accu QUAD:P1g", tquad
c       print*, "Ward QUAD:1g", Ward_q

c       print*, "Dble", result
c       print*, "QUAD", result_Q
c       print*, "ratio:P1g", abs(result/result_Q -1)

       endif ! printing

c to rescue Dble precision
       if (tquad.lt.tdble) then


          if(tquad.lt.GaugeAccuracyPenLnA) wardpAbe(2)=.true.

           result=result_Q     
            do jtemp=1,3
             resultgauge(jtemp)=resultgauge_Q(jtemp)
           enddo
      endif

       endif  ! tdble

#endif
       endif ! gaugetrue


       End


#ifdef WITH_QUAD

       subroutine penlineNoAbeTotal_QUAD(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mup2
     -   ,mup3,mup5,alpha,musqIn,gaugetest,comp,resultgauge,result,Div)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 3/7/2008
c Modified:4/8/2008
c ******************************************************************************
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p4(2),psi_p1(2),mup2(0:3),mup3(0:3),mup5(0:3) 
       Complex*16 result
       Complex*16 resultgauge(3),F(280)
       Real*8 musq
       Integer alpha,comp,gaugetest,Div
       Real*8 musqIn,P(685) 

       common/FfunctionsNOABEmmm/F,P
       SAVE/FfunctionsNOABEmmm/
       Common/musqInv/musq
       musq=musqIn
       If(Div.eq.0) then
       call penlineNoAbemmm_QUAD(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mup2
     -   ,mup3,mup5,alpha,musq,gaugetest,comp,resultgauge,result)

       else
               print*,"NO QUAD DIV so far:Pen1g"
       call penlineNoAbemmmDiv(p1,p2,p3,p4,p5,barpsi_p4,psi_p1,mup2
     -   ,mup3,mup5,alpha,musq,gaugetest,comp,resultgauge,result,Div)

       endif

       End
#endif
