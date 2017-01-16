       subroutine NoAbe2gaPenLineTotal(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_graph,result_color,resultgauge,Div,gauge)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:30/7/2008
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p3(2),psi_p1(2),mup2(0:3),mup4(0:3),mup5(0:3) 
       Complex*16 result_graph(3),result_color(2)
       real*8 resultgauge(2)
       Real*8 musq,musqIn
       Integer alpha,comp,ngluon,posgluon, Div,gauge
       Common/musqInv/musq
       musq=musqIn
c************************************************************************************


       if(Div.eq.0) then
       call  NoAbe2gaPenLineT(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_graph,result_color,resultgauge,gauge)

       else

       call  NoAbe2gaPenLineDivT(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_graph,result_color,resultgauge,Div,gauge)

       Endif

       end

       subroutine NoAbe2gaPenLineT(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_graph,result_color,resultgauge,gauge)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:30/7/2008
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p3(2),psi_p1(2),mup2(0:3),mup4(0:3),mup5(0:3) 
       complex*16 mup2test(0:3)
       Complex*16 result_graph(3),result_color(2)
       Complex*16 result_color_g(2),resulth(3)
       real*8 resultgauge(2)
       Real*8 musq,musqIn
       Integer alpha,comp,ngluon,posgluon,gauge,i
       Common/musqInv/musq
       musq=musqIn
c************************************************************************************

      call  NoAbe2gaPenLine(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_graph,result_color)


        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
        enddo
    
        call  NoAbe2gaPenLine(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2test
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,-1,resulth,result_color_g)


        if(abs(resulth(3)).gt.abs(resulth(1))) then
           if(abs(resulth(3)).gt.abs(resulth(2))) then
              resultgauge(1)=(abs((resulth(1)+resulth(2))/(-resulth(3)) -1d0))
            ELSE
             resultgauge(1)=(abs((resulth(1)+resulth(3))/(-resulth(2)) -1d0)) 
            ENDIF
        ELSE
           if(abs(resulth(1)).gt.abs(resulth(2))) then
              resultgauge(1)=(abs((resulth(2)+resulth(3))/(-resulth(1)) -1d0))
            ELSE
             resultgauge(1)=(abs((resulth(1)+resulth(3))/(-resulth(2)) -1d0)) 
            ENDIF
        ENDIF

C        resultgauge(1)=(abs((resulth(1)+resulth(2))/(-resulth(3)) -1d0))
        
       endif


        end


       subroutine NoAbe2gaPenLineDivT(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_graph,result_color,resultgauge,Div,gauge)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:30/7/2008
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p3(2),psi_p1(2),mup2(0:3),mup4(0:3),mup5(0:3) 
       complex*16 mup2test(0:3)
       Complex*16 result_graph(3),result_color(2)
       real*8 resultgauge(2)
       Complex*16 result_color_g(2),resulth(3)
       Real*8 musq,musqIn
       Integer alpha,comp,ngluon,posgluon, Div,gauge,i
       Common/musqInv/musq
       musq=musqIn
c************************************************************************************

      call  NoAbe2gaPenLineDiv(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_graph,result_color,Div)


        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
        enddo
    
        call  NoAbe2gaPenLineDiv(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2test
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,-1,resulth,result_color_g,Div)


        if(abs(resulth(3)).gt.abs(resulth(1))) then
           if(abs(resulth(3)).gt.abs(resulth(2))) then
              resultgauge(1)=(abs((resulth(1)+resulth(2))/(-resulth(3)) -1d0))
            ELSE
             resultgauge(1)=(abs((resulth(1)+resulth(3))/(-resulth(2)) -1d0)) 
            ENDIF
        ELSE
           if(abs(resulth(1)).gt.abs(resulth(2))) then
              resultgauge(1)=(abs((resulth(2)+resulth(3))/(-resulth(1)) -1d0))
            ELSE
             resultgauge(1)=(abs((resulth(1)+resulth(3))/(-resulth(2)) -1d0)) 
            ENDIF
        ENDIF

C        resultgauge(1)=(abs((resulth(1)+resulth(2))/(-resulth(3)) -1d0))
        
       endif


        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c                             QUAD
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       subroutine NoAbe2gaPenLineTotalQ(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_graph,result_color,resultgauge,Div,gauge)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:30/7/2008
       IMPLICIT NONE
#include "VBFNLO/loops/Utilities/GaugeTest.inc"
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p3(2),psi_p1(2),mup2(0:3),mup4(0:3),mup5(0:3) 
       Complex*16 result_graph(3),result_color(2)
       real*8 resultgauge(2)
       Real*8 musq,musqIn
       Integer alpha,comp,ngluon,posgluon, Div,gauge
       Common/musqInv/musq

cccccccccccccccccccccccccccccccccccccc 
       real*8 tdble,tquad
       Complex*16 result_graph_Q(3),result_color_Q(2)
       real*8 resultgauge_Q(2)
       integer jtemp
       logical gaugetrue
       logical wardpAbe(2)
       common/wardpenNoAbe/wardpAbe

       if(gauge.eq.1) then
       gaugetrue=.true.
       else
       gaugetrue=.false.
       endif

       musq=musqIn
       resultgauge(1)=1d9
       resultgauge(2)=1d9
       tdble=1d9
       tquad=1d9
c************************************************************************************
       CALL  NoAbe2gaPenLineTotal(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_graph,result_color,resultgauge,Div,gauge)


       tdble=abs(resultgauge(1))

       if(tdble+1.eq.tdble) then
          tdble=1d9
       endif   

       if(tdble.ne.tdble) then
          tdble=1d9
       endif   

       if(gaugetrue) then
       wardpAbe(1)=.false.
       wardpAbe(2)=.false.

c       print*,"Inside new"
       if(tdble.gt.GaugeAccuracyPenLnA) then
c       print*, "gauge Accu Dble:P2g", tdble
c       print*, "WardP Dble:2g", Ward
        else
       wardpAbe(1)=.true.
       wardpAbe(2)=.true.
       endif
       if(tdble+1.eq.tdble) then
	tdble = 1d9
	print*, "NaN in line", __LINE__, __FILE__
       endif

       if(tdble.ne.tdble) then
	tdble = 1d9
	print*, "NaN in line", __LINE__, __FILE__
       endif


#ifdef WITH_QUAD
       if(tdble.gt.GaugeAccuracyPenLnA) then

        CALL  NoAbe2gaPenLineTotal_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_graph_Q,result_color_Q,resultgauge_Q,Div,gauge)

       tquad=abs(resultgauge_Q(1))
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
          print*,"QUAD:PenNo2g"
       endif   

       if(tquad.ne.tquad) then
          tquad=1d9
          print*,"QUAD:PenNo2g"
       endif   

       if(tdble.gt.GaugeAccuracyPenLnA) then
c          print*, "COMP", COMP
c          print*, "gauge",gauge
c       print*, "gauge Accu QUAD:P2g", tquad
c       print*, "Ward QUAD:2g", Ward_q
       do jtemp=1,3
c       print*, "Dble", result(jtemp)
c       print*, "QUAD", result_Q(jtemp)
c       print*, "ratio:P2g", abs(result_graph(jtemp)/result_graph_Q(jtemp) -1)
       enddo
       endif ! printing

c to rescue Dble precision
       if (tquad.lt.tdble) then

           if(tquad.lt.GaugeAccuracyPenLnA) wardpAbe(2)=.true.

        do jtemp=1,3
           result_graph(jtemp)=result_graph_Q(jtemp)
           if(jtemp.le.2) then
           result_color(jtemp)=result_color_Q(jtemp)
           resultgauge(jtemp)=resultgauge_Q(jtemp)
             endif
        enddo
        else
c           print*,"pen:2g"
c          stop
       endif
       endif  ! tdble


#endif
       endif ! gaugetrue
       end


#ifdef WITH_QUAD

       subroutine NoAbe2gaPenLineTotal_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_graph,result_color,resultgauge,Div,gauge)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:30/7/2008
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p3(2),psi_p1(2),mup2(0:3),mup4(0:3),mup5(0:3) 
       Complex*16 result_graph(3),result_color(2)
       real*8 resultgauge(2)
       Real*8 musq,musqIn
       Integer alpha,comp,ngluon,posgluon, Div,gauge
       Common/musqInv/musq
       musq=musqIn
c************************************************************************************


       if(Div.eq.0) then
       call  NoAbe2gaPenLineT_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_graph,result_color,resultgauge,gauge)

       else

       print*,"NO QUAD DIV so far:Pen2g"
       call  NoAbe2gaPenLineDivT(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_graph,result_color,resultgauge,Div,gauge)

       Endif

       end


       subroutine NoAbe2gaPenLineT_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musqIn,ngluon,posgluon,comp,result_graph,result_color,resultgauge,gauge)
c ************************************************************************************
c ************************************************************************************
c Author: Francisco Campanario
c E-mail: francam@particle.uni-karlsruhe.de
c Date: 1/4/2008
c Modified:30/7/2008
       IMPLICIT NONE
       Real*8   p1(0:3),p2(0:3),p3(0:3),p4(0:3),p5(0:3)
       Complex*16 barpsi_p3(2),psi_p1(2),mup2(0:3),mup4(0:3),mup5(0:3) 
       complex*16 mup2test(0:3)
       Complex*16 result_graph(3),result_color(2)
       Complex*16 result_color_g(2),resulth(3)
       real*8 resultgauge(2)
       Real*8 musq,musqIn
       Integer alpha,comp,ngluon,posgluon,gauge,i
       Common/musqInv/musq
       musq=musqIn
c************************************************************************************

      call  NoAbe2gaPenLine_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,comp,result_graph,result_color)


        if(gauge.eq.1)then

        do i=0,3
          mup2test(i)=p2(i)
        enddo
    
        call  NoAbe2gaPenLine_QUAD(p1,p2,p3,p4,p5,barpsi_p3,psi_p1,mup2test
     &   ,mup4,mup5,alpha,musq,ngluon,posgluon,-1,resulth,result_color_g)
c        print*, "comp",comp
c        print*,"result2g",resulth(1)
c        print*,"result2g",resulth(2)
c        print*,"result2g",resulth(3)

        if(abs(resulth(3)).gt.abs(resulth(1))) then
           if(abs(resulth(3)).gt.abs(resulth(2))) then
              resultgauge(1)=(abs((resulth(1)+resulth(2))/(-resulth(3)) -1d0))
            ELSE
             resultgauge(1)=(abs((resulth(1)+resulth(3))/(-resulth(2)) -1d0)) 
            ENDIF
        ELSE
           if(abs(resulth(1)).gt.abs(resulth(2))) then
              resultgauge(1)=(abs((resulth(2)+resulth(3))/(-resulth(1)) -1d0))
            ELSE
             resultgauge(1)=(abs((resulth(1)+resulth(3))/(-resulth(2)) -1d0)) 
            ENDIF
        ENDIF


       endif


        end

#endif
